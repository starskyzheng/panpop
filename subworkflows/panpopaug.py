import os

def getvcf_aug(wildcards):
    GRAPH = config['graph']
    MAP = config['mapper']
    MINQ = str(config['MAP_MINQ'])
    return('2.callSV/{sample}/{sample}-' + GRAPH + '.' + MAP + '.aug.q' + MINQ + '.call.ext.vcf.gz')

def getpack_aug(wildcards):
    GRAPH = config['graph']
    MAP = config['mapper']
    MINQ = str(config['MAP_MINQ'])
    return('2.callSV/{sample}/{sample}-' + GRAPH + '.' + MAP + '.aug.q' + MINQ + '.pack')

def getpg_aug(wildcards):
    GRAPH = config['graph']
    MAP = config['mapper']
    MINQ = str(config['MAP_MINQ'])
    return('2.callSV/{sample}/{sample}-' + GRAPH + '.' + MAP + '.aug.pg')


# Realign
rule aug_realign0:
    input:
        vcf = '3.merge_rawvcf/4.filter_raw_vcf.{chrm}.vcf.gz',
        ref_fasta_file = GRAPH + '.gfa.fa'
    output:
        vcf = '4.realign/1.realign1.{chrm}.vcf.gz',
        vcf_sorted = '4.realign/1.realign1.{chrm}.sorted.vcf.gz',
    log:
        'logs/4.1.realign1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=config['mem_realign']
    params:
        realign_extend_bp_max = config['realign_extend_bp_max'],
        realign_extend_bp_min = config['realign_extend_bp_min'],
        tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir} --level 4 >> {log} 2>&1
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """


rule vcf2poss:
    input: '4.realign/1.realign1.{chrm}.sorted.vcf.gz',
    #input: '3.merge_rawvcf/3.merge_same_pos.{chrm}.vcf.gz',
    output: '6.aug_dp/1.{chrm}.poss',
    log:  'logs/6.1.vcf2poss.{chrm}.log',
    threads: 1
    shell:
        """
        perl {workflow.basedir}/scripts/vcf2pos_range.pl {input} {output} >> {log} 2>&1
        """

rule cal_aug_dp:
    input:
        poss = '6.aug_dp/1.{chrm}.poss',
        pack = getpack_aug,
        pg = getpg_aug,
    output: 
        list_packpg = '6.aug_dp/2.dps.{chrm}/{sample}.list_packpg',
        dpfile = '6.aug_dp/2.dps.{chrm}/{sample}.depth.txt.gz',
    log: 'logs/6.2.cal_aug_dp.{chrm}.{sample}.log',
    threads: config['cores_aug_filldp'],
    resources:
        mem_mb=config['mem_aug_filldp'],
    params:
        outdir = '6.aug_dp/2.dps.{chrm}',
        realthreads = 1,
    shell:
        """
        echo '{wildcards.sample} {input.pack} {input.pg}' > {output.list_packpg} && \
        perl {workflow.basedir}/scripts/cal_range_depth_aug.pl --vcfposs {input.poss} --list_packpg {output.list_packpg} --outdir {params.outdir} -t {params.realthreads} --chr {wildcards.chrm} >> {log} 2>&1
        """

rule merge_dp_vcf:
    input:
        dpfile = '6.aug_dp/2.dps.{chrm}/{sample}.depth.txt.gz',
        vcf = getvcf_aug,
        ref = GRAPH + '.gfa.fa',
    output:
        vcf = '6.aug_dp/3.vcf_with_dp.{chrm}/{sample}.vcf.gz',
        vcf_sorted = '6.aug_dp/3.vcf_with_dp.{chrm}/{sample}.sort.vcf.gz',
        vcf_sortedtbi = '6.aug_dp/3.vcf_with_dp.{chrm}/{sample}.sort.vcf.gz.tbi',
    params:
        min_dp = config['aug_nomut_min_dp'],
        min_cov = config['aug_nonmut_min_cov'],
    log: 'logs/6.3.merge_dp_vcf.{chrm}.{sample}.log',
    threads: config['cores_aug_filldp']
    resources:
        mem_mb=config['mem_aug_filldp']
    shell:
        """
        perl {workflow.basedir}/scripts/cal_range_depth_aug.fillvcf.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref {input.ref} --min_dp {params.min_dp} --min_cov {params.min_cov} --chr {wildcards.chrm} --dp_file {input.dpfile} >> {log} 2>&1 && \
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1 && \
        {TABIX} {output.vcf_sorted}
        """



rule aug_merge_dp2_infos:
    input:
        DPinfos = expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack.DPinfo', graph=GRAPH, map=MAPPER, minq=MINQ, sample=SAMPLES),
        ref_fasta_file = GRAPH + '.gfa.fa',
    output:
        merged='2.callSV.aug.DPinfos.txt',
    log: '2.callSV/logs/DPinfo.merged.log.txt',
    run:
        reflen = cal_fa_len(input.ref_fasta_file)
        #print(reflen)
        DPsums = dict()
        for DPinfo in input.DPinfos:
            with open(DPinfo, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    else:
                        line = line.strip().split('\t')
                        sample = line[0]
                        DPsum = line[2]
                        DPsums[sample] = DPsum
        OF = open(output.merged, 'w')
        OF.write('sample\tDP\n')
        for sample in DPsums:
            DPsum = int(DPsums[sample])
            DP = DPsum / reflen
            OF.write('{sample}\t{DP}\n'.format(sample=sample, DP=DP))
        #exit(-1)


rule aug_merge_rawvcfs_gen_list:
    input:
        vcfs = expand('6.aug_dp/3.vcf_with_dp.{{chrm}}/{sample}.sort.vcf.gz', sample=SAMPLES),
    output:
        outfile = '7.aug_merge_rawvcf/1.inputvcfs.{chrm}.list',
    run:
        with open(output.outfile, 'w') as f:
            for vcf in input.vcfs:
                f.write(vcf + '\n')

rule aug_merge_rawvcfs:
    input:
        vcfslist = '7.aug_merge_rawvcf/1.inputvcfs.{chrm}.list'
    output:
        vcf = '7.aug_merge_rawvcf/2.merge_rawvcf.{chrm}.vcf.gz'
    log:
        'logs/7.2.merge_rawvcf.{chrm}.log'
    threads: 44
    shell:
        # {BCFTOOLS} merge -m none --non_normalize_alleles -o {output.vcf} -O z --threads {threads} -l {input.vcfslist} >> {log} 2>&1
        """
        perl {workflow.basedir}/scripts/merge_vcf.pl --inlist {input.vcfslist} --out {output.vcf} --tmp_dir {ZTMPDIR} --vcfs_per_run 10000 --threads 11 --bcftools_threads 4 --m_none 1 >> {log} 2>&1
        """

rule aug_merge_same_pos:
    input:
        vcf = '7.aug_merge_rawvcf/2.merge_rawvcf.{chrm}.vcf.gz'
    output:
        vcf = '7.aug_merge_rawvcf/3.merge_same_pos.{chrm}.vcf.gz'
    log:
        'logs/7.3.merge_same_pos.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=2000
    shell:
        """
        perl {workflow.basedir}/scripts/merge_vcf_same_pos.pl --invcf {input.vcf} --outvcf {output.vcf} --threads {threads} >> {log} 2>&1
        """

rule aug_filter_raw_vcf:
    input:
        vcf = '7.aug_merge_rawvcf/3.merge_same_pos.{chrm}.vcf.gz',
        depthfile = '2.callSV.aug.DPinfos.txt',
    output:
        vcf = '7.aug_merge_rawvcf/4.filter_raw_vcf.{chrm}.vcf.gz',
    log:
        'logs/7.4.filter_raw_vcf.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=2000
    params:
        dp_min_fold = config['dp_min_fold'],
        dp_max_fold = config['dp_max_fold'],
        mad_min_fold = config['mad_min_fold'],
    shell:
        """
        perl {workflow.basedir}/scripts/flt_raw_vcfs.pl --invcf {input.vcf} --outvcf {output.vcf} --miss_threshold 1 --threads {threads} --dp_min_fold {params.dp_min_fold} --dp_max_fold {params.dp_max_fold} --mad_min_fold {params.mad_min_fold} --depth_file {input.depthfile} --remove_fail >> {log} 2>&1
        """

# Realign
rule aug_realign1:
    input:
        vcf = '7.aug_merge_rawvcf/4.filter_raw_vcf.{chrm}.vcf.gz',
        ref_fasta_file = GRAPH + '.gfa.fa'
    output:
        vcf = '8.aug_realign/1.realign1.{chrm}.vcf.gz',
        vcf_sorted = '8.aug_realign/1.realign1.{chrm}.sorted.vcf.gz',
    log:
        'logs/8.1.realign1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=config['mem_realign']
    params:
        realign_extend_bp_max = config['realign_extend_bp_max'],
        realign_extend_bp_min = config['realign_extend_bp_min'],
        tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir}  >> {log} 2>&1 && \
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """


rule aug_split_svp_sv_realign1:
    input:
        vcf = '8.aug_realign/1.realign1.{chrm}.sorted.vcf.gz',
    output:
        SNPs_vcf = '8.aug_realign/1.SV.split.{chrm}.snp.vcf.gz',
        SVs_vcf = '8.aug_realign/1.SV.split.{chrm}.sv.vcf.gz', # Augment output SV VCF
        INDELs_vcf = '8.aug_realign/1.SV.split.{chrm}.indel.vcf.gz',
    params:
        SV_min_length = config['SV_min_length'],
        outprefix = "8.aug_realign/1.SV.split.{chrm}",
    log:
        'logs/8.1x.split_snp_sv.{chrm}.vcf.gz.log',
    shell:
        """
	perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf}  {params.outprefix} {params.SV_min_length} >> {log} 2>&1
        """

rule aug_realign2:
    input:
        vcf = '8.aug_realign/1.realign1.{chrm}.sorted.vcf.gz',
        ref_fasta_file = GRAPH + '.gfa.fa'
    output:
        vcf = '8.aug_realign/2.realign2.{chrm}.vcf.gz',
        vcf_sorted = '8.aug_realign/2.realign2.{chrm}.sorted.vcf.gz',
    log:
        'logs/8.2.realign2.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=config['mem_realign']
    params:
        realign_extend_bp_max = config['realign_extend_bp_max'],
        realign_extend_bp_min = config['realign_extend_bp_min'],
        tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir} --level 6 >> {log} 2>&1 && \
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """


rule aug_split_svp_sv_realign2:
    input:
        vcf = '8.aug_realign/2.realign2.{chrm}.sorted.vcf.gz',
    output:
        SNPs_vcf = '8.aug_realign/2.SNP.split.{chrm}.snp.vcf.gz', # Augment output SNP/InDel VCF
        SVs_vcf = '8.aug_realign/2.SNP.split.{chrm}.sv.vcf.gz',
        INDELs_vcf = '8.aug_realign/2.SNP.split.{chrm}.indel.vcf.gz',
    params:
        SV_min_length = config['SV_min_length'],
        outprefix = "8.aug_realign/2.SNP.split.{chrm}",
    log:
        'logs/8.2x.split_snp_sv.{chrm}.vcf.gz.log',
    shell:
        """
	    perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf}  {params.outprefix} {params.SV_min_length} >> {log} 2>&1
        """

rule aug_thin_complexVCF:
    input:
        vcf = '8.aug_realign/1.realign1.{chrm}.sorted.vcf.gz'
    output:
        vcf1 = '8.aug_realign/3.thin1.{chrm}.vcf.gz',
        vcf2 = '8.aug_realign/3.thin2.{chrm}.vcf.gz',
        vcf2_sorted = '8.aug_realign/3.thin2.{chrm}.sorted.vcf.gz'
    log:
        'logs/8.3.aug_thin_complexVCF.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    params:
        sv2pav_merge_identity_threshold = config['sv2pav_merge_identity_threshold'],
        tmpdir = config['memory_tmp_dir'],
        max_len_tomerge = 5,
        sv_min_dp = config['SV_min_length'],
        sv2pav_merge_diff_threshold = 40,
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl --type 3 --invcf {input.vcf} --outvcf {output.vcf1} --threads {threads} --sv2pav_merge_identity_threshold {params.sv2pav_merge_identity_threshold} --tmpdir {params.tmpdir} --sv2pav_merge_diff_threshold {params.sv2pav_merge_diff_threshold} >> {log} 2>&1 && \
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {output.vcf1} --outvcf {output.vcf2} --threads {threads} --sv_min_dp {params.sv_min_dp} --max_len_tomerge {params.max_len_tomerge} >> {log} 2>&1 && \
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf2_sorted} -O z {output.vcf2} >> {log} 2>&1
        """

rule aug_thin1_realign2:
    input:
        vcf = '8.aug_realign/3.thin2.{chrm}.sorted.vcf.gz',
        ref_fasta_file = GRAPH + '.gfa.fa'
    output:
        vcf = '8.aug_realign/4.realign1.{chrm}.vcf.gz',
        vcf_sorted = '8.aug_realign/4.realign1.{chrm}.sorted.vcf.gz',
    log:
        'logs/8.4.realign1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=config['mem_realign']
    params:
        realign_extend_bp_max = config['realign_extend_bp_max'],
        realign_extend_bp_min = config['realign_extend_bp_min'],
        tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir} --level 1 >> {log} 2>&1
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """

rule aug_thin2_complexVCF:
    input:
        vcf = '8.aug_realign/4.realign1.{chrm}.sorted.vcf.gz'
    output:
        vcf1 = '8.aug_realign/5.thin1.{chrm}.vcf.gz',
        vcf2 = '8.aug_realign/5.thin2.{chrm}.vcf.gz',
        vcf2_sorted = '8.aug_realign/5.thin2.{chrm}.sorted.vcf.gz'
    log:
        'logs/8.5.aug_thin2_complexVCF.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    params:
        sv2pav_merge_identity_threshold = config['sv2pav_merge_identity_threshold'],
        tmpdir = config['memory_tmp_dir'],
        max_len_tomerge = 5,
        sv_min_dp = config['SV_min_length'], 
        sv2pav_merge_diff_threshold = 40,
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl --type 3 --invcf {input.vcf} --outvcf {output.vcf1} --threads {threads} --sv2pav_merge_identity_threshold {params.sv2pav_merge_identity_threshold} --tmpdir {params.tmpdir} --sv2pav_merge_diff_threshold {params.sv2pav_merge_diff_threshold} >> {log} 2>&1 && \
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {output.vcf1} --outvcf {output.vcf2} --threads {threads} --sv_min_dp {params.sv_min_dp} --max_len_tomerge {params.max_len_tomerge} >> {log} 2>&1 && \
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf2_sorted} -O z {output.vcf2} >> {log} 2>&1
        """


rule aug_filter_maf2:
    input:
        vcf = '8.aug_realign/5.thin2.{chrm}.sorted.vcf.gz'
    output:
        vcf = '8.aug_realign/6.filter_maf1.{chrm}.vcf.gz'
    log:
        'logs/8.6.filter_maf1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    params:
        min_maf = config['genotype_MAF'],
        max_miss_freq = config['max_missing_rate']
    shell:
        """
        perl {workflow.basedir}/scripts/flt_vcf_maf_by_allele.pl --in {input.vcf} --out {output.vcf} --min_maf {params.min_maf} --max_miss_freq {params.max_miss_freq} --threads {threads} >> {log} 2>&1
        """


rule aug_split_vcf_by_type_splitchr:
    input:
        vcf = '8.aug_realign/6.filter_maf1.{chrm}.vcf.gz'
    output:
        vcf_snp = '8.aug_realign/7.PAV.{chrm}.snp.vcf.gz',
        vcf_indel = '8.aug_realign/7.PAV.{chrm}.indel.vcf.gz',
        vcf_sv = '8.aug_realign/7.PAV.{chrm}.sv.vcf.gz',
    log:
        'logs/8.7.{chrm}.split_vcf.log'
    params:
        outprefix = '8.aug_realign/7.PAV.{chrm}',
        min_sv_len = config['SV_min_length']
    shell:
        """
        perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf} {params.outprefix} {params.min_sv_len} >> {log} 2>&1
        """


# Finally
rule aug_merge_vcf_nosplitchrs:
    input:
        vcf = expand('8.aug_realign/{{prefix}}.{chrm}.vcf.gz', chrm=CHRS),
    output:
        vcf = '9.aug_final_result/{prefix}.all.all.vcf.gz'
    threads:
        1
    run:
        invcf = os.path.abspath(input.vcfs)
        outvcf = os.path.abspath(output.vcf)
        os.system("ln -s {} {}".format(invcf, outvcf))


rule aug_merge_vcf_splitchrs:
    input:
        vcfs = expand('8.aug_realign/{{prefix}}.{chrm}.vcf.gz', chrm=CHRS),
    output:
        vcf = '9.aug_final_result/{prefix}.final_mergechr.all.vcf.gz'
    threads:
        6
    log:
        'logs/9.{prefix}.log'
    shell:
        """
        {BCFTOOLS} concat -o {output.vcf} -O z --threads {threads} {input.vcfs} >> {log} 2>&1
        """

rule aug_merge_vcf_splitchrs2:
    input:
        vcfs = expand('8.aug_realign/{{prefix}}.{chrm}.{{prefix1}}.vcf.gz', chrm=CHRS),
    output:
        vcf = '9.aug_final_result/{prefix}.final_mergechr.all.{prefix1}.vcf.gz'
    threads:
        6
    log:
        'logs/9.{prefix}.{prefix1}.log'
    shell:
        """
        {BCFTOOLS} concat -o {output.vcf} -O z --threads {threads} {input.vcfs} >> {log} 2>&1
        """

