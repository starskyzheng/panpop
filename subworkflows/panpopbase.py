def get_depth_file(wildcards):
    if config['mode'] == 'genotype':
        return('2.callSV.DPinfos.txt')
    elif config['mode'] == 'augment':
        return('2.callSV.aug.DPinfos.txt')
    else:
        print("Error mode!")
        exit(-1)


def get_rawvcf_files(wildcards):
    if config['mode'] == 'genotype':
        return expand('2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.ext.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ)
    elif config['mode'] == 'augment':
        return expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.ext.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ)
    else:
        print("Error mode!")
        exit(-1)

rule merge_dp2_infos:
    input:
        DPinfos = expand('2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack.DPinfo', graph=GRAPH, map=MAPPER, minq=MINQ, sample=SAMPLES),
        ref_fasta_file = GRAPH + '.gfa.fa',
    output:
        merged='2.callSV.DPinfos.txt',
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


# merge raw vcfs
rule merge_rawvcfs_gen_list:
    input:
        vcfs = get_rawvcf_files,
    output:
        outfile = '3.merge_rawvcf/1.inputvcfs.list'
    run:
        with open(output.outfile, 'w') as f:
            for vcf in input.vcfs:
                f.write(vcf + '\n')

if config['split_chr']==True:
    rule merge_rawvcfs:
        input:
            vcfslist = '3.merge_rawvcf/1.inputvcfs.list'
        output:
            vcf = '3.merge_rawvcf/2.merge_rawvcf.{chrm}.vcf.gz'
        log:
            'logs/3.2.merge_rawvcf.{chrm}.log'
        threads: 4
        shell:
            """
            {BCFTOOLS} merge -m none -o {output.vcf} -O z --threads {threads} -l {input.vcfslist} -r {wildcards.chrm} > {log} 2>&1
            """
else : # no chr split
    rule merge_rawvcfs:
        input:
            vcfslist = '3.merge_rawvcf/1.inputvcfs.list'
        output:
            vcf = '3.merge_rawvcf/2.merge_rawvcf.{chrm}.vcf.gz'
        log:
            'logs/3.2.merge_rawvcf.{chrm}.log'
        threads: 4
        shell:
            """
            {BCFTOOLS} merge -m none -o {output.vcf} -O z --threads {threads} -l {input.vcfslist} > {log} 2>&1
            """

rule merge_same_pos:
    input:
        vcf = '3.merge_rawvcf/2.merge_rawvcf.{chrm}.vcf.gz'
    output:
        vcf = '3.merge_rawvcf/3.merge_same_pos.{chrm}.vcf.gz'
    log:
        'logs/3.3.merge_same_pos.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=2000
    shell:
        """
        perl {workflow.basedir}/scripts/merge_vcf_same_pos.pl --invcf {input.vcf} --outvcf {output.vcf} --threads {threads} > {log} 2>&1
        """

rule filter_raw_vcf:
    input:
        vcf = '3.merge_rawvcf/3.merge_same_pos.{chrm}.vcf.gz',
        depthfile = get_depth_file,
    output:
        vcf = '3.merge_rawvcf/4.filter_raw_vcf.{chrm}.vcf.gz',
    log:
        'logs/3.4.filter_raw_vcf.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=2000
    params:
        dp_min_fold = config['dp_min_fold'],
        dp_max_fold = config['dp_max_fold'],
        mad_min_fold = config['mad_min_fold'],
    shell:
        """
        perl {workflow.basedir}/scripts/flt_raw_vcfs.pl --invcf {input.vcf} --outvcf {output.vcf} --miss_threshold 1 --threads {threads} --dp_min_fold {params.dp_min_fold} --dp_max_fold {params.dp_max_fold} --mad_min_fold {params.mad_min_fold} --depth_file {input.depthfile} --remove_fail > {log} 2>&1
        """

# Realign
rule realign1:
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
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir}  > {log} 2>&1
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} >> {log} 2>&1
        """

rule filter_maf1:
    input:
        vcf = '4.realign/1.realign1.{chrm}.sorted.vcf.gz'
    output:
        vcf = '4.realign/2.filter_maf1.{chrm}.vcf.gz'
    log:
        'logs/4.2.filter_maf1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    params:
        min_maf = config['genotype_MAF'],
        max_miss_freq = config['max_missing_rate']
    shell:
        """
        perl {workflow.basedir}/scripts/flt_vcf_maf_by_allele.pl --in {input.vcf} --out {output.vcf} --min_maf {params.min_maf} --max_miss_freq {params.max_miss_freq} --threads {threads} > {log} 2>&1
        """


rule realign2:
    input:
        vcf = '4.realign/2.filter_maf1.{chrm}.vcf.gz',
        ref_fasta_file = GRAPH + '.gfa.fa'
    output:
        vcf = '4.realign/3.realign2.{chrm}.vcf.gz',
        vcf_sorted = '4.realign/3.realign2.{chrm}.sorted.vcf.gz',
    log:
        'logs/4.3.realign1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    resources:
        mem_mb=config['mem_realign']
    params:
        realign_extend_bp_max = config['realign_extend_bp_max'],
        realign_extend_bp_min = config['realign_extend_bp_min'],
        tmpdir = config['memory_tmp_dir']
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref_fasta_file} --threads {threads} --ext_bp_max {params.realign_extend_bp_max} --ext_bp_min {params.realign_extend_bp_min} --tmpdir {params.tmpdir}  > {log} 2>&1
        {BCFTOOLS} sort --temp-dir {ZTMPDIR}/ -o {output.vcf_sorted} -O z {output.vcf} > {log} 2>&1
        """

rule filter_maf2:
    input:
        vcf = '4.realign/3.realign2.{chrm}.sorted.vcf.gz'
    output:
        vcf = '4.realign/4.filter_maf1.{chrm}.vcf.gz'
    log:
        'logs/4.4.filter_maf1.{chrm}.vcf.gz.log'
    threads: config['cores_realign']
    params:
        min_maf = config['genotype_MAF'],
        max_miss_freq = config['max_missing_rate']
    shell:
        """
        perl {workflow.basedir}/scripts/flt_vcf_maf_by_allele.pl --in {input.vcf} --out {output.vcf} --min_maf {params.min_maf} --max_miss_freq {params.max_miss_freq} --threads {threads} > {log} 2>&1
        """

# Finally
rule split_vcf_by_type:
    input:
        vcf = '4.realign/4.filter_maf1.{chrm}.vcf.gz'
    output:
        vcf_all = '5.final_result/1.final.{chrm}.all.vcf.gz',
        vcf_snp = '5.final_result/2.final.{chrm}.snp.vcf.gz',
        vcf_indel = '5.final_result/2.final.{chrm}.indel.vcf.gz',
        vcf_sv = '5.final_result/2.final.{chrm}.sv.vcf.gz'
    log:
        'logs/5.split_vcf.{chrm}.log'
    params:
        outprefix = '5.final_result/2.final.{chrm}',
        min_sv_len = config['SV_min_length']
    shell:
        """
        cp {input.vcf} {output.vcf_all}
        perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf} {params.outprefix} {params.min_sv_len} > {log} 2>&1
        """


rule merge_vcf_splitchrs:
    input:
        vcfs = expand('4.realign/4.filter_maf1.{chrm}.vcf.gz', chrm=CHRS),
    output:
        vcf = '5.final_result/1.final_mergechr.all.vcf.gz'
    threads:
        6
    log:
        'logs/5.1.merge.log'
    shell:
        """
        {BCFTOOLS} concat -o {output.vcf} -O z --threads {threads} {input.vcfs} > {log} 2>&1
        """

rule split_vcf_by_type_splitchr:
    input:
        vcf = '5.final_result/1.final_mergechr.all.vcf.gz'
    output:
        vcf_snp = '5.final_result/2.final_mergechr.snp.vcf.gz',
        vcf_indel = '5.final_result/2.final_mergechr.indel.vcf.gz',
        vcf_sv = '5.final_result/2.final_mergechr.sv.vcf.gz'
    log:
        'logs/5.2.split_vcf.log'
    params:
        outprefix = '5.final_result/2.final_mergechr',
        min_sv_len = config['SV_min_length']
    shell:
        """
        perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf} {params.outprefix} {params.min_sv_len} > {log} 2>&1
        """
