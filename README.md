# <a name="intro"></a>PanPop

### Application of graph-based pan-genome to population based on resequnces.  
  Panpop aims to combination `genome pan-graph` and `population resequnces` together. By using graph-genome to get more accurate population SVs and even SNPs in standard [VCF format][vcf], which is convenient for subsequent analysis.  
  PanPop contains two mode. `genotype` mode will genotype SVs in the graph-genome which means no novel SVs nor SNPs were available. `augment` mode will also extract novel SNPs and SVs for each sample.  
  PanPop compatibility to single server, [ssh-based cluster][snakemake_ssh], [slurm-based cluster][snakemake_slurm], grid and cloud environments. Most part of PanPop was written in parallel which is also easy to install.


- [Usage](#Usage)
  - [Quick Start](#example)
  - [Installation](#install)
  - [Parameters](#parameters)
  - [Notice](#notice)
- [Citations](#cite)


## <a name="usage"></a>Usage
  A graph-based genome in GFA format and resequnce reads were needed. 
  Reference gfa should be renamed as Ref.gfa.
  Created a list file contains paired-end seq reads (like example/1.sample.reads.list).
  Modify `workdir` and `sample_reads_list_file` in config.yaml.
  You are ready to go!
### <a name="example"></a>Example:
```sh
git clone https://github.com/StarSkyZheng/panpop.git
cd panpop
snakemake -j 3 --reason --printshellcmds
```
  Results located in `example/5.final_result` or `example/9.aug_final_result` for augment mode.
 
## <a name="install"></a>Installation:
  Dependencies: python3 & perl>=5.24   
```sh
pip3 install --user snakemake
curl -L https://cpanmin.us | perl - App::cpanminus
cpanm Data::Dumper MCE::Flow MCE::Candy MCE::Channel MCE::Shared Getopt::Long List::Util Carp File::Spec YAML Tie::CharArray IPC::Open2 File::Temp  
```

## <a name="parameters"></a>Parameters:
All parameters were defined in `config.yaml` file  
#### Basic Parameters
  `workdir`: Dir contains sample_read.list and will store result files.

  `sample_reads_list_file`: File name of sample_read.list. Must located in wirkdir

  `split_chr`: (False or True). Weather split by chromosome for parallel running. This is useful for multi-node clster. This option can greatly reduce the memory usage in augment mode. If you only own one computer-node the recommended setting is False. Default is False.

  `mode`: (genotype or augment). `genotype` mode will genotype SVs in the graph-genome. `augment` mode will extract novel SNPs and SVs for each sample. Default is `genotype`.


#### Filter Parameters:
  `MAP_MINQ`: minimal mapping-quality of each reads for giraffe-mapping in VG. Default is 5  

  `MAF`: Remove alleles than lower than this value. This is necessary and will greatly reduce the complexity of SVs. Default is 0.01  

  `max_missing_rate`: Max missing rate allowed for each SNPs and SVs. 0 means no missing allowed and 1 means not filter by missing. Default is 0.3 (remove this site if more than 30% samples were missing)  

  `dp_min_fold` & `dp_max_fold`: Hard-filter based on depth. Min depth for each site is XX fold of mean depth of this sample. Default is 1/3 ~ 3  

  `mad_min_fold`: Hard-filter based on depth: Min MAD(Minimum site allele depth) for each site is XX fold of mean depth of this sample. Default is 1/3 * 1/3  

#### Realign Parameters:
  `realign_extend_bp_max` & `realign_extend_bp_min`: During realign progress, the near by SNP/SVs will be merged together before realign. Values that are too low or too large may lose some accuracy. Default is 10 and 1.  

#### More parameters
  `SV_min_length`: Minimal length for SVs. Variat more than one base and smaller than this value will treated as InDels. Default is 50.  

  `realign_max_try_times_per_method`: Max try-times of each align software. Default is 3.  

  `memory_tmp_dir`: Temprory directory in memory. Must be a very fast disk. Left space can be smaller than 100Mb. Default is /run/user/USERID

  `mapper`: Maping algorithm. Can be either 'map' or 'gaffe'

  `aug_nonmut_min_cov`: In augment mode, if the coverage of SV in reference allele is greater than this value will be treated as exists. Note this value is the first filter parameter, the further filter based on depth will be perfomed. Defalut is 0.8.

  `aug_nomut_min_dp`: In augment mode, if the average of depth of SV in reference allele is greater than this value will be treated as exists. Note this value is the first filter parameter, the further filter based on depth will be perfomed. Default is 3.

## <a name=notice></a>Notice
  Graph-genome should be [GFA format][gfa] and the name of chromosome must be specified.  
  Chromesome name cannot be `all`.  
  Backbone genome should be numeric only. Non-Backbone genome should NOT be numeric only.  
  `scripts/build_graph.pl` can be used to build suitable gfa format from genome fasta files.   

## <a name=cite></a>Citations
  This software used `HAlign`, `bcftools`, `vg`, `muscle` and `snakemake` software: 
- [1] Shixiang Wan and Quan Zou, HAlign-II: efficient ultra-large multiple sequence alignment and phylogenetic tree reconstruction with distributed and parallel computing, Algorithms for Molecular Biology, 2017, 12:25.
- [2] Danecek P, Bonfield JK, et al. Twelve years of SAMtools and BCFtools. Gigascience (2021) 10(2):giab008.
- [3] Hickey, G., Heller, D., Monlong, J. et al. Genotyping structural variants in pangenome graphs using the vg toolkit. Genome Biol 21, 35 (2020).
- [4] Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797.
- [5] Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.

[gfa]: http://gfa-spec.github.io/GFA-spec/GFA1.html
[rgfa]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
[snakemake_ssh]: https://github.com/StarSkyZheng/snakemake_ssh
[snakemake_slurm]: https://snakemake.readthedocs.io/en/stable/executing/cluster.html
[vcf]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
