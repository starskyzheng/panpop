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
