# include: "util.py"


rule merge_vcf_same_pos:
    input:
        vcf = "04_consensus_vcf/{sample}/01.merged.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "04_consensus_vcf/{sample}/02.merge_vcf_same_pos.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.02.{sample}.merge_vcf_same_pos.log"
    resources:
        mem_mb = 8000
    shell:
        """
        perl {workflow.basedir}/scripts/merge_vcf_same_pos.pl --invcf {input.vcf} --outvcf {output.vcf} --ref {input.ref} --skip_mut_at_same_pos 2 --threads {threads} --ignore_dp >>{log} 2>&1
        """


rule gen_merge_mask:
    input:
        vcf = "04_consensus_vcf/{sample}/02.merge_vcf_same_pos.sorted.vcf.gz",
    output:
        bed = "04_consensus_vcf/{sample}/02.merge_vcf_same_pos.sorted.vcf.gz.mask.bed",
    log: "logs/3.02.{sample}.gen_merge_mask.log"
    resources:
        mem_mb = 8000
    shell:
        """
        perl {workflow.basedir}/scripts/realign_gen_mask.pl -i {input.vcf} -o {output.bed} >>{log} 2>&1
        """

rule realign1:
    input:
        vcf = "04_consensus_vcf/{sample}/02.merge_vcf_same_pos.sorted.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "04_consensus_vcf/{sample}/03.realign1.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.03.{sample}.realign1.log"
    resources:
        mem_mb = 40000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 1 --skip_mut_at_same_pos 2 --ext_bp_min 100 --ext_bp_max 400 --first_merge >>{log} 2>&1
        """

rule realign12:
    input:
        vcf = "04_consensus_vcf/{sample}/03.realign1.sorted.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "04_consensus_vcf/{sample}/03.realign2.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.03.{sample}.realign12.log"
    resources:
        mem_mb = 40000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 6 --skip_mut_at_same_pos 2 >>{log} 2>&1
        """


rule thin11:
    input:
        vcf = "04_consensus_vcf/{sample}/03.realign2.sorted.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/04.thin1.vcf.gz"
    threads: config['cores_realign']
    resources:
        mem_mb = 10000
    log: "logs/3.04.{sample}.thin11.log"
    params:
        sv2pav_merge_diff_threshold = config['sv2pav_merge_diff_threshold'],
        sv2pav_merge_identity_threshold = config['sv2pav_merge_identity_threshold'],
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl  --invcf {input.vcf} --outvcf {output.vcf} --sv2pav_merge_diff_threshold {params.sv2pav_merge_diff_threshold} --sv2pav_merge_identity_threshold {params.sv2pav_merge_identity_threshold} --threads {threads} >>{log} 2>&1
        """

# perl {workflow.basedir}/scripts/sv2pav.pl --invcf 7.thin1.vcf.gz --outvcf 7.thin2.vcf.gz --max_len_tomerge 5 --sv_min_dp 50
rule thin12:
    input:
        vcf = "04_consensus_vcf/{sample}/04.thin1.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/04.thin2.vcf.gz"
    threads: config['cores_realign']
    log: "logs/3.04.{sample}.thin12.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 40 --threads {threads} >>{log} 2>&1
        """

rule realign2:
    input:
        vcf = "04_consensus_vcf/{sample}/04.thin2.vcf.gz",
        ref = INDEX_REF
    output:
        vcf = "04_consensus_vcf/{sample}/05.realign2.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.05.{sample}.realign2.log"
    resources:
        mem_mb = 20000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 4 --skip_mut_at_same_pos 2 --ext_bp_min 100 --ext_bp_max 400 >>{log} 2>&1
        """

rule thin21:
    input:
        vcf = "04_consensus_vcf/{sample}/05.realign2.sorted.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/06.thin1.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.06.{sample}.thin21.log"
    params:
        sv2pav_merge_diff_threshold = config['sv2pav_merge_diff_threshold'],
        sv2pav_merge_identity_threshold = config['sv2pav_merge_identity_threshold'],
    resources:
        mem_mb = 10000
    shell:
        """
        perl {workflow.basedir}/scripts/merge_similar_allele.pl --invcf {input.vcf} --outvcf {output.vcf} --sv2pav_merge_identity_threshold {params.sv2pav_merge_identity_threshold} --sv2pav_merge_diff_threshold {params.sv2pav_merge_diff_threshold} --threads {threads} --sv2pav_merge_identity_threshold {params.sv2pav_merge_identity_threshold} --sv2pav_merge_diff_threshold {params.sv2pav_merge_diff_threshold} >>{log} 2>&1
        """
    
rule thin22:
    input:
        vcf = "04_consensus_vcf/{sample}/06.thin1.sorted.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/06.thin2.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.06.{sample}.thin22.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 40 --threads {threads} >>{log} 2>&1
        """

def get_min_support_caller(wildcards):
    fa = get_assb_fasta(wildcards)
    if fa == '-':
        return MIN_SUPPORT_CALLER_WITHOUT_ASSB
    else:
        return MIN_SUPPORT_CALLER_WITH_ASSB

rule merge_callers:
    input:
        vcf = "04_consensus_vcf/{sample}/06.thin2.sorted.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/07.merge_callers.vcf.gz"
    threads: 1
    log: "logs/3.07.{sample}.merge_callers.log"
    resources:
        mem_mb = 4000
    params:
        ext_parm = " --force_merge_insertion ",
        min_support_callers = get_min_support_caller
    shell:
        """
        perl {workflow.basedir}/scripts/long_caller_merger.pl --in {input.vcf} --min_supporting {params.min_support_callers} -s {wildcards.sample} --out {output.vcf} {params.ext_parm} >>{log} 2>&1
        """

rule realign3:
    input:
        vcf = "04_consensus_vcf/{sample}/07.merge_callers.vcf.gz",
        ref = INDEX_REF,
        mask_bed = "04_consensus_vcf/{sample}/02.merge_vcf_same_pos.sorted.vcf.gz.mask.bed"
    output:
        vcf = "04_consensus_vcf/{sample}/08.realign_merge.vcf.gz"
    threads: config['cores_realign']
    log: "logs/3.08.{sample}.realign3.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/realign.pl --chr_tolerance --in_vcf {input.vcf} --out_vcf {output.vcf} --ref_fasta_file {input.ref} --threads {threads} --level 1 --skip_mut_at_same_pos 2 --ext_bp_min 1000 --ext_bp_max 1000 --mask_bed_file {input.mask_bed} --skip_snp 1 >>{log} 2>&1
        """

rule thin32:
    input:
        vcf = "04_consensus_vcf/{sample}/08.realign_merge.vcf.gz",
    output:
        vcf = "04_consensus_vcf/{sample}/09.thin2.unsorted.vcf.gz",
    threads: config['cores_realign']
    log: "logs/3.06.{sample}.thin22.log"
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/sv2pav.pl --invcf {input.vcf} --outvcf {output.vcf} --enable_norm_alle 1 --sv_min_dp 40 --threads {threads} --max_len_tomerge 0 >>{log} 2>&1
        """


rule concat_inv:
    input:
        vcf1 = "04_consensus_vcf/{sample}/09.thin2.sorted.vcf.gz",
        vcf2 = "04_consensus_vcf/{sample}/00.inv.vcf.gz",
        tbi1 = "04_consensus_vcf/{sample}/09.thin2.sorted.vcf.gz.tbi",
        tbi2 = "04_consensus_vcf/{sample}/00.inv.vcf.gz.tbi"
    log: "logs/3.10.{sample}.concat_inv.log"
    output:
        vcf = "04_consensus_vcf/{sample}/10.concat_inv.vcf.gz", 
    shell:
        """
        {BCFTOOLS} concat -a {input.vcf1} {input.vcf2} 2>>{log} | {BCFTOOLS} sort 2>>{log} | bgzip -c > {output.vcf}
        """

rule split_snp_indel_sv:
    input:
        vcf = "04_consensus_vcf/{sample}/10.concat_inv.vcf.gz"
    output:
        sv = "04_consensus_vcf/{sample}/10.concat_inv.split.sv.vcf.gz",
        indel = "04_consensus_vcf/{sample}/10.concat_inv.split.indel.vcf.gz",
        snp = "04_consensus_vcf/{sample}/10.concat_inv.split.snp.vcf.gz"
    log: "logs/3.10.{sample}.split_snp_indel_sv.log"
    threads: 1
    resources:
        mem_mb = 4000
    params:
        prefix = "04_consensus_vcf/{sample}/10.concat_inv.split"
    shell:
        """
        perl {workflow.basedir}/scripts/vcf_split_snp_indel_sv.pl {input.vcf} {params.prefix} 50 >>{log} 2>&1
        """

rule merge_samples_genlist:
    input:
        vcfs = expand("04_consensus_vcf/{sample}/09.thin2.sorted.vcf.gz", sample=SAMPLE_INDEX)
    output:
        listfile = '05_merge_samples/01.merge_samples.vcfs.list',
    run:
        with open(output.listfile, 'w') as f:
            for vcf in input.vcfs:
                f.write(vcf + '\n')
