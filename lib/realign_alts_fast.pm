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

static int is_missing(const char* str, STRLEN len) {
    return (len == 0) || (len == 1 && str[0] == '-');
}

void sarray_is_same_miss_fast(SV* sarray_ref, IV max_alts) {
    Inline_Stack_Void;
    if (!SvROK(sarray_ref) || SvTYPE(SvRV(sarray_ref)) != SVt_PVAV) {
        croak("sarray_is_same_miss_fast expects an array reference");
    }
    AV* sarray = (AV*)SvRV(sarray_ref);
    SSize_t last_index = (max_alts >= 0) ? max_alts : av_len(sarray);
    if (last_index < 0) {
        Inline_Stack_Push(sv_2mortal(newSViv(1)));
        Inline_Stack_Push(sv_2mortal(newSViv(0)));
        Inline_Stack_Push(sv_2mortal(newSViv(0)));
        Inline_Stack_Done;
        return;
    }
    SV** elem = av_fetch(sarray, 0, 0);
    if (!elem) {
        croak("sarray_is_same_miss_fast: missing element 0");
    }
    STRLEN first_len = 0;
    const char* first_ptr = SvPVbyte(*elem, first_len);
    int is_ref_miss = is_missing(first_ptr, first_len);
    int is_miss = is_ref_miss;
    int is_same = 1;
    for (SSize_t i = 1; i <= last_index; i++) {
        elem = av_fetch(sarray, i, 0);
        if (!elem) {
            croak("sarray_is_same_miss_fast: missing element");
        }
        STRLEN len = 0;
        const char* ptr = SvPVbyte(*elem, len);
        if (is_missing(ptr, len)) {
            is_miss = 1;
        }
        if (is_same) {
            if (len != first_len) {
                is_same = 0;
            } else if (len > 0 && memcmp(ptr, first_ptr, len) != 0) {
                is_same = 0;
            }
        }
    }
    Inline_Stack_Push(sv_2mortal(newSViv(is_same)));
    Inline_Stack_Push(sv_2mortal(newSViv(is_miss)));
    Inline_Stack_Push(sv_2mortal(newSViv(is_ref_miss)));
    Inline_Stack_Done;
}

SV* clone_sarray_fast(SV* sarray_ref) {
    if (!SvROK(sarray_ref) || SvTYPE(SvRV(sarray_ref)) != SVt_PVAV) {
        croak("clone_sarray_fast expects an array reference");
    }
    AV* source = (AV*)SvRV(sarray_ref);
    AV* copy = newAV();
    SSize_t last = av_len(source);
    for (SSize_t i = 0; i <= last; i++) {
        SV** elem = av_fetch(source, i, 0);
        SV* value;
        if (!elem || !SvOK(*elem)) {
            value = newSV(0);
        } else {
            STRLEN len = 0;
            const char* ptr = SvPVbyte(*elem, len);
            value = newSVpvn(ptr, len);
        }
        av_push(copy, value);
    }
    return newRV_noinc((SV*)copy);
}

IV cal_sarray_is_compatible_fast(SV* sarray1_ref, SV* sarray2_ref) {
    if (!SvROK(sarray1_ref) || SvTYPE(SvRV(sarray1_ref)) != SVt_PVAV) {
        croak("cal_sarray_is_compatible_fast expects array references");
    }
    if (!SvROK(sarray2_ref) || SvTYPE(SvRV(sarray2_ref)) != SVt_PVAV) {
        croak("cal_sarray_is_compatible_fast expects array references");
    }
    AV* arr1 = (AV*)SvRV(sarray1_ref);
    AV* arr2 = (AV*)SvRV(sarray2_ref);
    SSize_t len1 = av_len(arr1);
    SSize_t len2 = av_len(arr2);
    if (len1 != len2) {
        return 0;
    }
    HV* map1 = newHV();
    HV* map2 = newHV();
    IV next1 = 0;
    IV next2 = 0;
    int result = 1;
    for (SSize_t i = 0; i <= len1; i++) {
        SV** elem1 = av_fetch(arr1, i, 0);
        SV** elem2 = av_fetch(arr2, i, 0);
        if (!elem1 || !elem2) {
            result = 0;
            break;
        }
        STRLEN len_a = 0;
        const char* str_a = SvPVbyte(*elem1, len_a);
        STRLEN len_b = 0;
        const char* str_b = SvPVbyte(*elem2, len_b);

        HE* he1 = hv_fetch(map1, str_a, (I32)len_a, 0);
        HE* he2 = hv_fetch(map2, str_b, (I32)len_b, 0);
        IV idx1;
        IV idx2;
        if (he1) {
            idx1 = SvIV(HeVAL(he1));
        } else {
            idx1 = next1++;
            hv_store(map1, str_a, (I32)len_a, newSViv(idx1), 0);
        }
        if (he2) {
            idx2 = SvIV(HeVAL(he2));
        } else {
            idx2 = next2++;
            hv_store(map2, str_b, (I32)len_b, newSViv(idx2), 0);
        }
        if (idx1 != idx2) {
            result = 0;
            break;
        }
    }
    SvREFCNT_dec((SV*)map1);
    SvREFCNT_dec((SV*)map2);
    return result;
}
END_C

sub available { 1 }

1;
