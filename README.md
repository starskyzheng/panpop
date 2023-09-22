# <a name="intro"></a>PanPop

### Sequence-aware SV merging and processing pipeline designed for a single individual, as well as TGS-based and NGS-based populations.

  &emsp;Merging structural variations (SVs) at the population level presents a significant challenge, yet it is essential for conducting thorough phenotypic analyses, especially in the era of pangenomics. Here we introduce PanPop, a new tool that utilizes an advanced sequence-aware SV merging algorithm to efficiently merge SVs of various types. PanPop demonstrates exceptional expertise in merging and optimizing the majority of multiallelic SVs into informative biallelic variants, showcasing superior precision and reduced rates of missing data when compared to alternative software solutions.
  
  &emsp;PanPop compatibility to single server, [ssh-based cluster][snakemake_ssh], [slurm-based cluster][snakemake_slurm], grid and cloud environments. Most part of PanPop was written in parallel which is also easy to install.


- [Usage](#Usage)
  - [Quick Start](#example)
  - [Installation](#install)
  - [Parameters For TGS data](#parametersTGS)
  - [Parameters For NGS data](#parametersNGS)
  - [Detailed parameters](#parametersDetailed)
  - [Input GFA Reference file](#gfa_ref)
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
snakemake -j 3 --reason --printshellcmds -s Snakefile_TGS # For TGS dataset
snakemake -j 3 --reason --printshellcmds -s Snakefile_NGS # For NGS dataset
```
  Results located in `example/example_long_reads` and `example/example_short_reads`.

  This was a minimal example, which should be completed within 1 minute and produce vcfs with no SVs/SNPs.
 
## <a name="install"></a>Installation:
### Installation using conda (recommended)
```sh
conda env create -f conda.panpop.yaml
conda activate panpop
```

### Installation manually
  Install Software, and then add path to `config/software.yaml`. Some software were available in `misc/`
  
  `famsa, muscle3, bcftools, vg, tabix, bgzip, minimap2, ngmlr, sniffles, cutesv, svim, picard, pbsv, nucmer, Assemblytics, bedtools, samtools`

  Runtime dependencies: python3 & perl>=5.24   
```sh
pip3 install --user snakemake
curl -L https://cpanmin.us | perl - App::cpanminus
cpanm Data::Dumper MCE::Flow MCE::Candy MCE::Channel MCE::Shared Getopt::Long List::Util Carp File::Spec YAML Tie::CharArray IPC::Open2 File::Temp  
```

## <a name="parametersTGS"></a>Parameters For TGS data:
  This part were defined in `configs/config.TGS.yaml`
  
  `workdir`: Dir contains sample_read.list and will store result files.

  `sample_reads_list_file`: File name of sample_read.list. Must located in wirkdir
  
  `mode`: (simple or enhanced). `simple` mode will call SVs through 5 SV-caller and merge using PanPop `enhanced` recalculated depth information and then merge SVs using PanPop. Default is `enhanced`.

  `ref`: reference file in fasta format

  `MIN_SUPPORT_CALLER`: minimum support caller number for SVs. Default is 2

  `SV_CALLER_LIST`: SV caller list to be used. Default is all 5 SV callers. Options: [sniffles, cuteSV, svim, pbsv, Assemblytics]

  `mapper`: Long reads mapper, default is `minimap2`. Options: [minimap2, ngmlr]. Default is 'minimap2'.


## <a name="parametersNGS"></a>Parameters For NGS data:
  This part were defined in `configs/config.NGS.yaml`
  
  `workdir`: Dir contains sample_read.list and will store result files.

  `sample_reads_list_file`: File name of sample_read.list. Must located in wirkdir

  `split_chr`: (False or True). Weather split by chromosome for parallel running. This is useful for multi-node clster. This option can greatly reduce the memory usage in augment mode. If you only own one computer-node the recommended setting is False. Default is False.

  `mode`: (genotype or augment). `genotype` mode will genotype SVs in the graph-genome. `augment` mode will extract novel SNPs and SVs for each sample. Default is `genotype`.
  
  `graph`: Graph-based reference should be renamed as `Ref.gfa` and then placed in `workdir`. This value is not recommand for modifing! Default is `Ref`.

  `mapper`: Maping algorithm. Can be either 'map' or 'gaffe'. Default is 'gaffe'.

## <a name="parametersDetailed"></a>Detailed parameters
  This part were defined in `configs/base.yaml`. This config were used by both TGS and NGS pipeline.
  
### Filter Parameters:
  `MAP_MINQ`: minimal mapping-quality of each reads for giraffe-mapping in VG. Default is 5  

  `MAF`: Remove alleles than lower than this value. This is recommanded and will greatly reduce the complexity of SVs. Default is 0.01  

  `max_missing_rate`: Max missing rate allowed for each SNPs and SVs. 0 means no missing allowed and 1 means not filter by missing. Default is 0.3 (remove this site if more than 30% samples were missing)  

  `dp_min_fold` & `dp_max_fold`: Hard-filter based on depth. Min depth for each site is XX fold of mean depth of this sample. Default is 1/3 ~ 3  

  `mad_min_fold`: Hard-filter based on depth: Min MAD(Minimum site allele depth) for each site is XX fold of mean depth of this sample. Default is 1/3 * 1/3 

  `mad_min_abs`: Hard-filter based on depth: Min MAD(Minimum site allele depth) for each site. Default is 1

### Realign Parameters:
  `realign_extend_bp_max` & `realign_extend_bp_min`: During realign progress, the near by SNP/SVs will be merged together before realign. Values that are too low or too large may lose some accuracy. Default is 10 and 1.  

### More parameters
  `SV_min_length`: Minimal length for SVs. Variat more than one base and smaller than this value will treated as InDels. Default is 50.  

  `realign_max_try_times_per_method`: Max try-times of each align software. Default is 1.  

  `memory_tmp_dir`: Temprory directory in memory. Must be a very fast disk. Left space can be smaller than 100Mb. Default is /run/user/USERID

  `aug_nonmut_min_cov`: In augment mode, if the coverage of SV in reference allele is greater than this value will be treated as exists. Note this value is the first filter parameter, the further filter based on depth will be perfomed. Defalut is 0.8.

  `aug_nomut_min_dp`: In augment mode, if the average of depth of SV in reference allele is greater than this value will be treated as exists. Note this value is the first filter parameter, the further filter based on depth will be perfomed. Default is 3.

  `sv2pav_merge_identity_threshold`: During thin process, minimal sequence similarity of alleles to be merged. Default is 0.8

  `ssindel_maxlen` and `ssindel_maxallele`: Simple or small indels. default is 10 and 3

### Performance & Resource
  Cores and memorys for snakemake were defineded here

## <a name=gfa_ref></a>Input GFA Reference file
  Graph-genome should be [GFA format][gfa] and the name of the chromosome must be specified which should be numeric only.
  
  The GFA file to be run must contain path lines starting with `P`. Each `P` line represents one chromosome of the backbone genome.
  
### This following script can be used to generate or remove Non-Backbone paths:
  
  `scripts/build_graph.pl` can be used to build suitable gfa format from genome fasta files. Noted that Non-Backbone genome should NOT be numeric only and chromesome name cannot be `all`.
  
  `scripts/build_graph.pl` can be used to parse [rGFA][rgfa] format (generated by minigraph) using `--rgfa`.
  
  `scripts/build_graph.pl` can be used to remove Non-Backbone paths of GFA format by using `--gfa`.

## <a name=cite></a>Citations
  This software used `FAMSA`, `HAlign`, `bcftools`, `vg`, `muscle` and `snakemake` software: 
- [1].Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific Reports, 6, 33964
- [2] Shixiang Wan and Quan Zou, HAlign-II: efficient ultra-large multiple sequence alignment and phylogenetic tree reconstruction with distributed and parallel computing, Algorithms for Molecular Biology, 2017, 12:25.
- [3] Danecek P, Bonfield JK, et al. Twelve years of SAMtools and BCFtools. Gigascience (2021) 10(2):giab008.
- [4] Hickey, G., Heller, D., Monlong, J. et al. Genotyping structural variants in pangenome graphs using the vg toolkit. Genome Biol 21, 35 (2020).
- [5] Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797.
- [6] Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.


[gfa]: http://gfa-spec.github.io/GFA-spec/GFA1.html
[rgfa]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
[snakemake_ssh]: https://github.com/StarSkyZheng/snakemake_ssh
[snakemake_slurm]: https://snakemake.readthedocs.io/en/stable/executing/cluster.html
[vcf]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
