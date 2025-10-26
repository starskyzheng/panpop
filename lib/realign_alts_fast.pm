package realign_alts_fast;

use strict;
use warnings;

use Exporter 'import';

our @EXPORT_OK = qw(
    sarray_is_same_miss_fast
    cal_sarray_is_compatible_fast
    clone_sarray_fast
    available
);

use Inline ( C => Config => CLEAN_AFTER_BUILD => 0 );
use Inline C => <<'END_C';
#include <string.h>
#include <stdlib.h>

typedef struct {
    size_t count;
    const char **items;
    STRLEN *lens;
    unsigned char *defined;
} SArray;

static void sarray_init(SArray *arr) {
    arr->count = 0;
    arr->items = NULL;
    arr->lens = NULL;
    arr->defined = NULL;
}

static void sarray_free(SArray *arr) {
    if (arr->items) {
        free((void *)arr->items);
    }
    if (arr->lens) {
        free((void *)arr->lens);
    }
    if (arr->defined) {
        free((void *)arr->defined);
    }
    sarray_init(arr);
}

static void sarray_from_av(SV *sarray_ref, SArray *arr, IV max_alts) {
    sarray_init(arr);
    if (!SvROK(sarray_ref) || SvTYPE(SvRV(sarray_ref)) != SVt_PVAV) {
        croak("expected array reference");
    }
    AV *av = (AV *)SvRV(sarray_ref);
    SSize_t highest = av_len(av);
    if (highest < 0) {
        arr->count = 0;
        return;
    }
    size_t total = (size_t)highest + 1;
    if (max_alts >= 0 && (size_t)(max_alts + 1) < total) {
        total = (size_t)max_alts + 1;
    }
    arr->count = total;
    arr->items = (const char **)malloc(total * sizeof(const char *));
    arr->lens = (STRLEN *)malloc(total * sizeof(STRLEN));
    arr->defined = (unsigned char *)malloc(total * sizeof(unsigned char));
    if (!arr->items || !arr->lens || !arr->defined) {
        sarray_free(arr);
        croak("memory allocation failed");
    }
    for (size_t i = 0; i < total; ++i) {
        SV **elem = av_fetch(av, (SSize_t)i, 0);
        if (!elem) {
            sarray_free(arr);
            croak("missing element when materialising sarray");
        }
        if (!SvOK(*elem)) {
            arr->defined[i] = 0;
            arr->items[i] = "";
            arr->lens[i] = 0;
        } else {
            arr->defined[i] = 1;
            arr->items[i] = SvPVbyte(*elem, arr->lens[i]);
        }
    }
}

static int is_missing_value(const SArray *arr, size_t idx) {
    if (!arr->defined[idx]) {
        return 1;
    }
    STRLEN len = arr->lens[idx];
    if (len == 0) {
        return 1;
    }
    if (len == 1 && arr->items[idx][0] == '-') {
        return 1;
    }
    return 0;
}

static int values_equal(const SArray *arr, size_t a, size_t b) {
    if (arr->defined[a] != arr->defined[b]) {
        return 0;
    }
    if (arr->lens[a] != arr->lens[b]) {
        return 0;
    }
    if (arr->lens[a] == 0) {
        return 1;
    }
    return memcmp(arr->items[a], arr->items[b], arr->lens[a]) == 0;
}

static void build_pattern(const SArray *arr, int *pattern) {
    int next_id = 0;
    for (size_t i = 0; i < arr->count; ++i) {
        int id = -1;
        for (size_t j = 0; j < i; ++j) {
            if (values_equal(arr, i, j)) {
                id = pattern[j];
                break;
            }
        }
        if (id == -1) {
            id = next_id++;
        }
        pattern[i] = id;
    }
}

void sarray_is_same_miss_fast(SV *sarray_ref, IV max_alts) {
    Inline_Stack_Void;
    SArray arr;
    sarray_from_av(sarray_ref, &arr, max_alts);

    int is_same = 1;
    int is_miss = 0;
    int is_ref_miss = 0;

    if (arr.count == 0) {
        is_same = 1;
        is_miss = 0;
        is_ref_miss = 0;
    } else {
        is_ref_miss = is_missing_value(&arr, 0);
        is_miss = is_ref_miss;
        for (size_t i = 1; i < arr.count; ++i) {
            if (is_missing_value(&arr, i)) {
                is_miss = 1;
            }
            if (is_same && !values_equal(&arr, 0, i)) {
                is_same = 0;
            }
        }
    }

    sarray_free(&arr);

    Inline_Stack_Push(sv_2mortal(newSViv(is_same)));
    Inline_Stack_Push(sv_2mortal(newSViv(is_miss)));
    Inline_Stack_Push(sv_2mortal(newSViv(is_ref_miss)));
    Inline_Stack_Done;
}

SV *clone_sarray_fast(SV *sarray_ref) {
    SArray arr;
    sarray_from_av(sarray_ref, &arr, -1);

    AV *copy = newAV();
    if (arr.count > 0) {
        av_extend(copy, (SSize_t)arr.count - 1);
        for (size_t i = 0; i < arr.count; ++i) {
            SV *value;
            if (!arr.defined[i]) {
                value = newSV(0);
            } else {
                value = newSVpvn(arr.items[i], arr.lens[i]);
            }
            av_store(copy, (SSize_t)i, value);
        }
    }

    sarray_free(&arr);
    return newRV_noinc((SV *)copy);
}

IV cal_sarray_is_compatible_fast(SV *sarray1_ref, SV *sarray2_ref) {
    SArray arr1;
    SArray arr2;
    sarray_from_av(sarray1_ref, &arr1, -1);
    sarray_from_av(sarray2_ref, &arr2, -1);

    IV result = 1;

    if (arr1.count != arr2.count) {
        result = 0;
    } else if (arr1.count > 0) {
        size_t n = arr1.count;
        int *pattern1 = (int *)malloc(n * sizeof(int));
        int *pattern2 = (int *)malloc(n * sizeof(int));
        if (!pattern1 || !pattern2) {
            if (pattern1) {
                free(pattern1);
            }
            if (pattern2) {
                free(pattern2);
            }
            sarray_free(&arr1);
            sarray_free(&arr2);
            croak("memory allocation failed");
        }
        build_pattern(&arr1, pattern1);
        build_pattern(&arr2, pattern2);
        for (size_t i = 0; i < n; ++i) {
            if (pattern1[i] != pattern2[i]) {
                result = 0;
                break;
            }
        }
        free(pattern1);
        free(pattern2);
    }

    sarray_free(&arr1);
    sarray_free(&arr2);
    return result;
}
END_C

sub available { 1 }

1;
