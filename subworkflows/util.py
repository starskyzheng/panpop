import re

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
    