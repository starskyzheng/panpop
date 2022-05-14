# PanPop
  Application of pan-genome to population based on resequnces

# Usage
  A graph-based genome in GFA format and resequnce reads were needed. 
  Reference gfa should be renamed as Ref.gfa.
  Created a list file contains paired-end seq reads (like example/1.sample.reads.list).
  Modify `workdir` and `sample_reads_list_file` in config.yaml.
  You are ready to go!
## example:
    git clone https://github.com/StarSkyZheng/panpop.git
    cd panpop
    snakemake -j 3 --reason --printshellcmds
  results located in example/5.final_result
 
# Install:
## Dependencies:
  Needed: python3 & perl>=5.24   
    `pip install snakemake`  
    `cpanm Data::Dumper MCE::Flow MCE::Candy Getopt::Long List::Util Carp  File::Spec YAML Coro::Generator MCE::Channel Tie::CharArray  IPC::Open2 File::Temp`  
  
  
# Parameters:
```
MAP_MINQ: minimal mapping-quality of each reads for giraffe-mapping in VG. Default is 5
MAF: Remove alleles than lower than this value. This is necessary and will greatly reduce the complexity of SVs. Default is 0.01
max_missing_rate: Max missing rate allowed for each SNPs and SVs. 0 means no missing allowed and 1 means not filter by missing. Default is 0.3 (remove this site if more than 30% samples were missing)
dp_min_fold/dp_max_fold: Hard-filter based on depth. Min depth for each site is XX fold of mean depth of this sample. Default is 1/3 ~ 3
mad_min_fold: Hard-filter based on depth: Min MAD(Minimum site allele depth) for each site is XX fold of mean depth of this sample. Default is 1/3 * 1/3
realign_extend: During realign progress, the near by SNP/SVs will be merged together before realign. Values that are too low or too large may lose some accuracy. Default is 10 and 1.
realign_extend_bp_max: 10
realign_extend_bp_min: 1
```
