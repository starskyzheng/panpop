
rule merge_samples:
    input:
        vcfs = expand("04_consensus_vcf/{sample}/09.thin2.sort.vcf.gz", sample=SAMPLE_INDEX),
        listfile = '05_merge_samples/01.merge_samples.vcfs.list',
    output:
        vcf = '05_merge_samples/01.merge_samples.vcf.gz'
    log: "logs/5.01.merge_samples.log"
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        """
        {BCFTOOLS} merge -m none -l {input.listfile} | {BCFTOOLS} sort | bgzip -c > {output.vcf} && \
        {TABIX} {output.vcf} >>{log} 2>&1
        """

# Realign
rule aug_realign0:
    input:
        vcf = '05_merge_samples/01.merge_samples.vcf.gz',
        ref_fasta_file = INDEX_REF
    output:
        vcf = '05_merge_samples/02.realign0.vcf.gz',
        vcf_sorted = '05_merge_samples/02.realign0.sorted.vcf.gz',
    log:
        'logs/5.02.realign1.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb = 4000
    params:
        realign_extend_bp_max = 500,
        realign_extend_bp_min = 50,
        #tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --skip_mut_at_same_pos 2 --level 1 --first_merge >> {log} 2>&1
        {BCFTOOLS} sort -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """

rule vcf2poss:
    input: '05_merge_samples/02.realign0.vcf.gz',
    output: '05_merge_samples/03.vcf2pos.tsv', # 1-based
    log:  'logs/5.03.vcf2poss.log',
    threads: 1
    shell:
        """
        perl {workflow.basedir}/scripts/vcf2pos_range.pl {input} {output} >> {log} 2>&1
        """

rule cal_aug_dp:
    input:
        poss = '05_merge_samples/03.vcf2pos.tsv',
	bam = lambda wildcards: "02_bam/{sample}.{platform}.ngmlr.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards)),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.ngmlr.sort.bam.bai".format(sample=wildcards.sample, platform=get_platform(wildcards))
    output:
        list_packpg = '05_merge_samples/04.dps/{sample}.list_packpg',
        dpfile = '05_merge_samples/04.dps/{sample}.depth.txt.gz',
        avgdp = '05_merge_samples/04.dps/{sample}.avgdp.txt',
    log: 'logs/5.04.cal_aug_dp.{sample}.log',
    threads: config['cores_aug_filldp'],
    resources:
        mem_mb=config['mem_aug_filldp'],
    params:
        outdir = '05_merge_samples/04.dps',
        realthreads = 1,
    shell:
        """
        echo '{wildcards.sample} {input.bam}' > {output.list_packpg} && \
        perl {workflow.basedir}/scripts/cal_range_depth_aug.pl --vcfposs {input.poss} --list_packpg {output.list_packpg} --outdir {params.outdir} --threads {params.realthreads} --input_format bam --output_avgdpinfo {output.avgdp} >> {log} 2>&1
        """

rule fill_aug_dp:
    input:
        vcf = '04_consensus_vcf/{sample}/09.thin2.sort.vcf.gz',
        dpfile = '05_merge_samples/04.dps/{sample}.depth.txt.gz',
        ref = INDEX_REF,
    output:
        vcf = '05_merge_samples/05.fill_aug_dp/{sample}.vcf.gz',
        vcf_sorted = '05_merge_samples/05.fill_aug_dp/{sample}.sorted.vcf.gz',
        vcf_sortedtbi = '05_merge_samples/05.fill_aug_dp/{sample}.sorted.vcf.gz.tbi',
    log: 'logs/5.05.fill_aug_dp.{sample}.log',
    threads: 1,
    resources:
        mem_mb=config['mem_aug_filldp'],
    params:
        min_dp = config['aug_nomut_min_dp'],
        min_cov = config['aug_nonmut_min_cov'],
    shell:
        """
        perl {workflow.basedir}/scripts/cal_range_depth_aug.fillvcf.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref {input.ref} --min_cov {params.min_cov} --min_dp {params.min_dp} --dp_file {input.dpfile} >> {log} 2>&1 && \
        {BCFTOOLS} sort -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1 && \
        {TABIX} {output.vcf_sorted} >> {log} 2>&1
        """


rule aug_merge_dp2_infos:
    input:
        DPinfos = expand('05_merge_samples/04.dps/{sample}.avgdp.txt', sample=SAMPLE_INDEX),
        ref_fasta_file = INDEX_REF,
    output:
        merged='05_merge_samples/04.DPinfos.txt',
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


rule merge_samples_genlist_filldp:
    input:
        vcfs = expand("05_merge_samples/05.fill_aug_dp/{sample}.sorted.vcf.gz", sample=SAMPLE_INDEX)
    output:
        listfile = '05_merge_samples/06.merge_samples.vcfs.list',
    run:
        with open(output.listfile, 'w') as f:
            for vcf in input.vcfs:
                f.write(vcf + '\n')

rule merge_samples_filldp:
    input:
        vcfs = expand("05_merge_samples/05.fill_aug_dp/{sample}.sorted.vcf.gz", sample=SAMPLE_INDEX),
        listfile = '05_merge_samples/06.merge_samples.vcfs.list',
    output:
        vcf = '05_merge_samples/06.merge_samples.vcf.gz'
    log: "logs/5.06.merge_samples.log"
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        """
        {BCFTOOLS} merge -m none -l {input.listfile} | {BCFTOOLS} sort | bgzip -c > {output.vcf} && \
        {TABIX} {output.vcf} >>{log} 2>&1
        """


rule filter_raw_vcf_by_depth:
    input:
        vcf = '05_merge_samples/06.merge_samples.vcf.gz',
        depthfile = '05_merge_samples/04.DPinfos.txt',
    output:
        vcf = '05_merge_samples/07.filter_raw_vcf_by_depth.vcf.gz',
    log: 'logs/5.07.filter_raw_vcf_by_depth.log'
    threads: config['cores_realign']
    resources:
        mem_mb = 4000
    params:
        dp_min_fold = config['dp_min_fold'],
        dp_max_fold = config['dp_max_fold'],
        mad_min_fold = config['mad_min_fold'],
    shell:
        """
        perl {workflow.basedir}/scripts/flt_raw_vcfs.pl --invcf {input.vcf} --outvcf {output.vcf} --miss_threshold 1 --threads {threads} --dp_min_fold {params.dp_min_fold} --dp_max_fold {params.dp_max_fold} --mad_min_fold {params.mad_min_fold} --depth_file {input.depthfile} --remove_fail --flt_append_only >> {log} 2>&1
        """


# Realign
rule aug_realign0_filldp:
    input:
        vcf = '05_merge_samples/07.filter_raw_vcf_by_depth.vcf.gz',
        ref_fasta_file = INDEX_REF
    output:
        vcf = '05_merge_samples/08.realign0.vcf.gz',
        vcf_sorted = '05_merge_samples/08.realign0.sorted.vcf.gz',
    log:
        'logs/5.08.realign1.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb = 4000
    params:
        realign_extend_bp_max = 500,
        realign_extend_bp_min = 50,
        #tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --skip_mut_at_same_pos 2 --level 1 --first_merge >> {log} 2>&1
        {BCFTOOLS} sort -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """



rule pop_thin11:
    input:
        vcf = "05_merge_samples/08.realign0.sorted.vcf.gz",
    output:
        vcf = "05_merge_samples/09.thin1.vcf.gz"
    threads: config['cores_realign']
    resources:
        mem_mb = 10000
    log: "logs/5.09.thin1.log"
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl --type 3 --invcf {input.vcf} --outvcf {output.vcf} --sv2pav_merge_diff_threshold 40 --sv2pav_merge_identity_threshold 0.6 --threads {threads} >>{log} 2>&1
        """

# perl {workflow.basedir}/scripts/sv2pav.pl --invcf 7.thin1.vcf.gz --outvcf 7.thin2.vcf.gz --max_len_tomerge 5 --sv_min_dp 50
rule pop_thin12:
    input:
        vcf = "05_merge_samples/09.thin1.vcf.gz",
    output:
        vcf = "05_merge_samples/10.thin2.vcf.gz"
    threads: config['cores_realign']
    log: "logs/5.10.thin2.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 40 --threads {threads} >>{log} 2>&1
        """


rule pop_realign2:
    input:
        vcf = "05_merge_samples/10.thin2.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "05_merge_samples/11.realign2.vcf.gz",
        vcfsorted = "05_merge_samples/11.realign2.sort.vcf.gz",
    threads: config['cores_realign']
    log: "logs/5.11.realign2.log"
    resources:
        mem_mb = 20000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 1 --skip_mut_at_same_pos 2 --ext_bp_min 100 --ext_bp_max 400 >>{log} 2>&1 && \
        {BCFTOOLS} sort -O z -o {output.vcfsorted} {output.vcf} >>{log} 2>&1
        """

rule pop_thin21:
    input:
        vcf = "05_merge_samples/11.realign2.sort.vcf.gz",
    output:
        vcf = "05_merge_samples/12.thin1.vcf.gz",
        vcfsorted = "05_merge_samples/12.thin1.sort.vcf.gz",
    threads: config['cores_realign']
    log: "logs/5.12.thin1.log"
    resources:
        mem_mb = 10000
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl --type 3 --invcf {input.vcf} --outvcf {output.vcf} --sv2pav_merge_diff_threshold 40 --sv2pav_merge_identity_threshold 0.5 --threads {threads} >>{log} 2>&1 && \
        {BCFTOOLS} sort -O z -o {output.vcfsorted} {output.vcf} >>{log} 2>&1
        """
    
rule pop_thin22:
    input:
        vcf = "05_merge_samples/12.thin1.sort.vcf.gz",
    output:
        vcf = "05_merge_samples/13.thin2.vcf.gz",
        vcfsorted = "05_merge_samples/13.thin2.sort.vcf.gz"
    threads: config['cores_realign']
    log: "logs/5.13.thin2.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 30 --threads {threads} >>{log} 2>&1 && \
        {BCFTOOLS} sort -O z -o {output.vcfsorted} {output.vcf} >>{log} 2>&1
        """

rule pop_realign3:
    input:
        vcf = "05_merge_samples/13.thin2.sort.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "05_merge_samples/14.realign3.vcf.gz",
        vcfsorted = "05_merge_samples/14.realign3.sort.vcf.gz",
    threads: config['cores_realign']
    log: "logs/5.14.realign3.log"
    resources:
        mem_mb = 20000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 6 --skip_mut_at_same_pos 2 --not_use_merge_alle_afterall 1 >>{log} 2>&1 && \
        {BCFTOOLS} sort -O z -o {output.vcfsorted} {output.vcf} >>{log} 2>&1
        """

rule pop_thin23:
    input:
        vcf = "05_merge_samples/14.realign3.sort.vcf.gz",
    output:
        vcf = "05_merge_samples/14.thin3.vcf.gz",
    threads: config['cores_realign']
    log: "logs/5.14.thin3.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 30 --threads {threads} --force_pav --mode by_type >>{log} 2>&1
        """

rule pop_split_by_type:
    input:
        vcf = "05_merge_samples/14.thin3.vcf.gz",
    output:
        vcf_snp = "05_merge_samples/15.thin3.snp.vcf.gz",
        vcf_indel = "05_merge_samples/15.thin3.indel.vcf.gz",
        vcf_sv = "05_merge_samples/15.thin3.sv.vcf.gz",
    threads: config['cores_realign']
    log: "logs/5.15.split_by_type.log"
    resources:
        mem_mb = 1000
    params:
        outprefix = '05_merge_samples/15.thin3',
        min_sv_len = config['SV_min_length']
    shell:
        """
        perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf} {params.outprefix} {params.min_sv_len} >> {log} 2>&1
        """

