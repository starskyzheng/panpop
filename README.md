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
  Results located in `example/5.final_result`
 
# Install:
  Dependencies: python3 & perl>=5.24   
    `pip install snakemake`  
    `cpanm Data::Dumper MCE::Flow MCE::Candy Getopt::Long List::Util Carp  File::Spec YAML Coro::Generator MCE::Channel Tie::CharArray  IPC::Open2 File::Temp`  
  
  
# Parameters:
All parameters were defined in `config.yaml` file  
### I/O File paths
`workdir`: Dir contains sample_read.list and will store result files.

`sample_reads_list_file`: File name of sample_read.list. Must located in wirkdir

### Filter Parameters:
  `MAP_MINQ`: minimal mapping-quality of each reads for giraffe-mapping in VG. Default is 5  

  `MAF`: Remove alleles than lower than this value. This is necessary and will greatly reduce the complexity of SVs. Default is 0.01  

  `max_missing_rate`: Max missing rate allowed for each SNPs and SVs. 0 means no missing allowed and 1 means not filter by missing. Default is 0.3 (remove this site if more than 30% samples were missing)  

  `dp_min_fold` & `dp_max_fold`: Hard-filter based on depth. Min depth for each site is XX fold of mean depth of this sample. Default is 1/3 ~ 3  

  `mad_min_fold`: Hard-filter based on depth: Min MAD(Minimum site allele depth) for each site is XX fold of mean depth of this sample. Default is 1/3 * 1/3  

### Realign Parameters:
  `realign_extend_bp_max` & `realign_extend_bp_min`: During realign progress, the near by SNP/SVs will be merged together before realign. Values that are too low or too large may lose some accuracy. Default is 10 and 1.  

### More parameters
  `SV_min_length`: Minimal length for SVs. Variat more than one base and smaller than this value will treated as InDels. Default is 50.  
  
  `realign_max_try_times_per_method`: Max try-times of each align software. Default is 3.  

  `memory_tmp_dir`: Temprory directory in memory. Must be a very fast disk. Left space can be smaller than 100Mb. Default is /run/user/USERID
