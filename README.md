# <a name="intro"></a>PanPop

[![DOI](https://zenodo.org/badge/485607258.svg)](https://zenodo.org/doi/10.5281/zenodo.10458868)

### Sequence-aware SV merging and processing pipeline designed for a single individual, as well as TGS-based and NGS-based populations.

  &emsp;Merging structural variations (SVs) at the population level presents a significant challenge, yet it is essential for conducting thorough phenotypic analyses, especially in the era of pangenomics. Here we introduce PanPop, a new tool that utilizes an advanced sequence-aware SV merging algorithm to efficiently merge SVs of various types. PanPop demonstrates exceptional expertise in merging and optimizing the majority of multiallelic SVs into informative biallelic variants, showcasing superior precision and reduced rates of missing data when compared to alternative software solutions.
  
  &emsp;PanPop compatibility to single server, [ssh-based cluster][snakemake_ssh], [slurm-based cluster][snakemake_slurm], grid and cloud environments. Most part of PanPop was written in parallel which is also easy to install.


- [Usage](#Usage)
  - [Quick Start: Pipeline](#pipe_example)
  - [Quick Start: Standalone](#s_example)
  - [Installation](#install)
  - [Parameters For TGS data](#parametersTGS)
  - [Parameters For NGS data](#parametersNGS)
  - [Detailed parameters](#parametersDetailed)
  - [Input GFA Reference file](#gfa_ref)
- [Citations](#cite)


## <a name="usage"></a>Usage
### <a name="pipe_example"></a>Pipeline Example:
#### For TGS data: 
  You will require TGS sequencing reads for each sample and a reference genome in FASTA format. Optionally, you can also provide genome assemblies for each sample to enable SV-caller based on assemblies.
  
  Create a list file containing paired-end TGS sequencing reads (refer to the example at `example/example_long_reads/1.sample.reads.tsv`).
  
#### For NGS data: 
  A graph-based genome in GFA format and NGS resequnce reads for each samples were needed. 
  
  Create a list file containing paired-end NGS sequencing reads (as shown in `example/example_short_reads/1.sample.reads.list`).
  
#### Then (For both TGS and NGS data):
  Modify `workdir` and `sample_reads_list_file` in config.yaml.
  
  You are ready to go!
  
```sh
git clone https://github.com/StarSkyZheng/panpop.git
cd panpop
snakemake -j 3 --reason --printshellcmds -s Snakefile_TGS # For TGS dataset, using configs/config.TGS.yaml
snakemake -j 3 --reason --printshellcmds -s Snakefile_NGS # For NGS dataset, using configs/config.NGS.yaml
```
  Results located in `example/example_long_reads` or `example/example_short_reads`.

  This was a minimal example, which should be completed within 1 minute and produce vcfs with no SVs/SNPs.

  This example declares the use of up to 3 threads. You can add threads via `-j`.

  Please note that the default temporary directory is configured as `/run/user/UID`. Should this directory not exist on your system, you will need to specify an alternative directory by modifying the `memory_tmp_dir` parameter in `configs/base.yaml`.

### <a name="s_example"></a>Standalone example:
  To facilitate the processing of VCFs originating from different pipelines, we also provide standalone access to the PART algorithm and the Fill-Depth-Information process.

  #### Standalone Entry for PART:
  There was an exapmle from CodeOcean: https://doi.org/10.24433/CO.1577027.v1
  
  PART (PAnpop Realign and Thin)'s standalone entrance: `panpop/bin/PART_run.pl` 
  
  Here was an example to run PART twice using 16 threads:
```sh
perl panpop/bin/PART_run.pl --in_vcf INPUT.vcf -o OUTDIR_RUN1 -r REFERENCE.FASTA  -t 16 --tmpdir TMPDIR
perl panpop/bin/PART_run.pl --in_vcf OUTDIR_RUN1/3.final.vcf.gz -o OUTDIR_RUN2 -r REFERENCE.FASTA -t 16 --tmpdir TMPDIR -not_first_merge
```
  Noted REF, ALT and GT (colunm 4, 5 and 9) need to be completed in the VCF file. (see https://github.com/starskyzheng/panpop/issues/18)
  #### Standalone Entry for Fill-depth-information:
  The standalone entry for fill-depth information is located at: `panpop/bin/Fill_DP.pl`.
  
  Please note that this process exclusively adds depth information to SV entries in the input VCF with './.'. Post-processing is required where the output VCF file must be filtered by depth to reduce the false-negative rate.
  
  This process needs a merged VCF file which contains all SVs of all samples, VCFs for each samples before merge and BAM file for each samples.
  
  Sample information must be compiled in a TSV format (tab-separated) with three columns: sample_name, BAM_path, VCF_path. (No header line)
```sh
perl panpop/bin/Fill_DP.pl --in_vcf INPUT.vcf --outdir OUTDIR --ref_fasta_file REFERENCE.FASTA --threads 16 --sid2files SAMPLE_INFORMATION.tsv
```

## <a name="install"></a>Installation:
### Installation using conda (recommended)
This will create a new conda env named `panpop`. This should be done less than 1 hour, but highly dependent on the network and hardware.
```sh
conda env create -f conda.panpop.yaml
conda activate panpop
conda install -c bioconda svim-asm
cpanm MCE::Channel Tie::CharArray
```

### Installation manually
  Install Software, and then add path to `config/software.yaml`. Some software were available in `misc/`
  
  `famsa, muscle3, bcftools, vg, tabix, bgzip, minimap2, ngmlr, sniffles, cutesv, svim, picard, pbsv, nucmer, svim-asm, Assemblytics, bedtools, samtools`

  Runtime dependencies: python3 & perl>=5.24   
```sh
pip3 install --user snakemake
curl -L https://cpanmin.us | perl - App::cpanminus
cpanm Data::Dumper MCE::Flow MCE::Candy MCE::Channel MCE::Shared Getopt::Long List::Util Carp File::Spec YAML Tie::CharArray IPC::Open2 File::Temp  
```

## <a name="parametersTGS"></a>Parameters For TGS data:
  This part were defined in `configs/config.TGS.yaml`
  
  `workdir`: Dir contains sample_read.list and will store result files.

  `sample_reads_list_file`: File name of sample_read.list. Must located in wirkdir. If assembly fasta is absent, just replace with `-`.
  
  `mode`: (simple or enhanced). `simple` mode will call SVs through 5 SV-caller and merge using PanPop `enhanced` recalculated depth information and then merge SVs using PanPop. Default is `enhanced`.

  `ref`: reference file in fasta format

  `MIN_SUPPORT_CALLER`: minimum support caller number for SVs. Default is 2

  `SV_CALLER_LIST`: SV caller list to be used. Default is all 5 SV callers. Options: [sniffles, cuteSV, svim, pbsv, svim-asm, Assemblytics]

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
  `memory_tmp_dir`: Temprory directory in memory. This directory must exists and must be a very fast disk. Left space can be smaller than 100Mb. Default is /run/user/USERID

  `SV_min_length`: Minimal length for SVs. Variat more than one base and smaller than this value will treated as InDels. Default is 50.  

  `realign_max_try_times_per_method`: Max try-times of each align software. Default is 1.  

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

### TGS SVs -> NGS SVs (Both performed by PanPop)
  `scripts/prase_vcf_rm_asterisk.pl` can be used to remove asterisk in the ALT colums of VCF file, and then you were able to use VG to consturct gfa file. 
  
  If you needs to process NGS data by using TGS SVs generated by PanPop (in VCF format), this process is necessary. Noted this process will increase multi-allelic. Usage: `perl scripts/prase_vcf_rm_asterisk.pl --infile IN.vcf.gz -r REF.fa -o OUT1.vcf.gz`. 
  
  Then using VG to consturct gfa file: `vg construct -r REF.fa -v OUT1.vcf.gz -t 4 > OUT2.vg; vg view --vg-in --gfa OUT2.vg > OUT3.gfa`.

## <a name=cite></a>Citations
  Zheng, Z. et al. A sequence-aware merger of genomic structural variations at population scale. Nat Commun (2024). https://doi.org/10.1038/s41467-024-45244-9


[gfa]: http://gfa-spec.github.io/GFA-spec/GFA1.html
[rgfa]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
[snakemake_ssh]: https://github.com/StarSkyZheng/snakemake_ssh
[snakemake_slurm]: https://snakemake.readthedocs.io/en/stable/executing/cluster.html
[vcf]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
