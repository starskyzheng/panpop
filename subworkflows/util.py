import re
import os
import sys
import shutil


def exe_check_ret_abs_path(exe):
    if exe is None:
        raise Exception("Error!!! config error: {}".format(exe))
    if os.path.exists(exe) and os.path.isfile(exe) and os.access(exe, os.X_OK):
        abspath = os.path.abspath(exe)
        return abspath
    elif shutil.which(exe) is not None:
        abspath = shutil.which(exe)
        return abspath
    else:
        raise Exception("Can't find executable: {}".format(exe))



def cal_fa_len(fa_file):
    """
    Calculate the length of fasta file
    """
    if type(fa_file) == type(set()):
        fa_file = list(fa_file)[0]
    #print(fa_file)
    fa_len = 0
    with open(fa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            fa_len += len(line.strip())
    return fa_len

def fa2chrs(fa_file):
    """
    Read chrs of fasta file
    """
    chrs = []
    with open(fa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chr = re.search('>(\S+)', line).group(1)
                # chr must be numberic only
                try:
                    int(chr)
                except:
                    print('Error: chr is not numeric in {}: {}'.format(fa_file, chr), file=sys.stderr)
                    # die
                    sys.exit(1)
                chrs.append(chr)
    if len(chrs) == 0:
        print('Error: No chr found in {}! which might be no `P line` in gfa file'.format(fa_file), file=sys.stderr)
        # die
        sys.exit(1)
    return chrs
    
def gfa2chrs(gfa_file):
    """
    Read chrs from gfa file
    """
    chrs = []
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('P'):
                chr = re.search('^P\t(\S+)\t', line).group(1)
                # chr must be numberic only
                try:
                    int(chr)
                except:
                    print('Error: chr is not numeric in {}: {}'.format(gfa_file, chr), file=sys.stderr)
                    # die
                    sys.exit(1)
                chrs.append(chr)   
    if len(chrs) == 0:
        print('Error: No chr found in {}! which might be no `P line` in gfa file'.format(gfa_file), file=sys.stderr)
        # die
        sys.exit(1)             
    return chrs

def rgfa_read_chr(rgfa_file):
    print('Now Read chrs from rgfa file: {}'.format(rgfa_file) )
    chrs = dict()
    with open(rgfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                re_results = re.search( r'\t(\S+)\t\S+\t(\S+)$', line.strip() )
                if not re_results:
                    print("????")
                    print(line)
                SR_raw = re_results.group(2)
                SR_now = re.search(r'^SR:i:(\d+)$', SR_raw).group(1)
                #print(SR_now)
                if SR_now != '0': continue
                chr_raw = re_results.group(1)
                chr_now = re.search(r'^SN:Z:(\S+)$', chr_raw).group(1)
                if chr_now not in chrs:
                    print(chr_now, end=' ')
                    chrs[chr_now] = 1
                # chrs.append(chr)            
    print('')    
    return chrs.keys()

def read_chrsfile(chrs_file):
    print('Now Read chrs from rgfa.chrs file: {}'.format(chrs_file) )
    chrs = []
    with open(chrs_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line != '':
                chrs.append(line)
    return(chrs)

def rgfa2chrs(rgfa_file):
    chr_file = "{}.chrs".format(rgfa_file)
    if os.path.exists(chr_file):
        chrs = read_chrsfile(chr_file)
        return(chrs)
    else:
        chrs = rgfa_read_chr(rgfa_file)
        with open(chr_file, 'w') as f:
            for chr_now in chrs:
                f.write(chr_now)
                f.write('\n')
        return(chrs)

def get_assb_fasta(wildcards):
    sid = wildcards.sample
    try:
        #fa = S2T.iloc[S2T.index==sid, :].values[0][1]
        fa = S2T.loc[sid, 'AssmbFasta']
    except:
        fa = '-'
    #print("get_assb_fasta::: {} : {}".format(sid, fa), file=sys.stderr)
    if fa == '-':
        return '-'
    else:
        fa_path = os.path.abspath(fa)
        return fa_path

def get_platform(wildcards):
    sid = wildcards.sample
    try:
        #fa = S2T.iloc[S2T.index==sid, :].values[0][1]
        platform = S2T.loc[sid, 'platform']
        platform = platform.lower()
    except:
        print("Error! no platform for sid: {}".format(sid))
        1
    return platform



rule sort_vcf:
    input:
        vcf = '{dir_now}/{filename}.unsorted.vcf.gz'
    output:
        vcf = '{dir_now}/{filename}.sorted.vcf.gz'
    wildcard_constraints:
        dir_now='^[01234678].+'
    threads: 1
    shell:
        """
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf} -O z {input.vcf}
        """

rule sort_vcf_tbi:
    input:
        vcf = '{dir}/{filename}.vcf.gz'
    output:
        vcf = '{dir}/{filename}.vcf.gz.tbi'
    threads: 1
    shell:
        """
        tabix {input.vcf}
        """