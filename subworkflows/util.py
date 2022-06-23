import re
import os

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
                chrs.append(chr)                
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
                chrs.append(chr)                
    return chrs

def rgfa_read_chr(rgfa_file):
    print('Now Read chrs from rgfa file: {}'.format(rgfa_file) )
    chrs = dict()
    with open(rgfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                re_results = re.search( r'\t(\S+)\t\S+\t(\S+)$', line.strip() )
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
