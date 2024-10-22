import os
# include: "util.py"

# def fq_make_softlink(wildcards):
#     sid = wildcards.sample
#     fq = S2T.loc[S2T.index==sid, 'fq']
#     platform = get_platform(wildcards)
#     fq_path = os.path.abspath(fq)
#     fqdir = "01_raw_data"
#     fq_dst_basename = "{sample}.{platform}.fq.gz".format(sample=sid, platform=platform)
#     fq_dst_path = os.path.join(fqdir, fq_dst_basename)
#     if os.path.exists(fq_dst_path):
#         return fq_dst_path
#     else:
#         os.symlink(fq_path, fq_dst_path)
#         return fq_dst_path

rule fq_make_softlink:
    input:
        lambda wildcards: S2T.loc[S2T.index==wildcards.sample, 'fq']
    output:
        "01_raw_data/{sample}.{platform}.fq.gz"
    run:
        fq_path = os.path.abspath(input[0])
        fq_dst_path = os.path.abspath(output[0])
        # weather to make softlink
        if os.path.exists(fq_dst_path):
            return fq_dst_path
        else:
            os.symlink(fq_path, fq_dst_path)
            return fq_dst_path

def get_ngmlr_parms(wildcards):
    platform = wildcards.platform
    if platform == 'ont':
        return " -x ont "
    elif platform == 'pb':
        return " -x pacbio "
    elif platform == 'hifi':
        return " -x pacbio "
    else:
        print("Error! platform error: {}".format(platform))
        1

def get_minimap2_parms(wildcards):
    platform = wildcards.platform
    if platform == 'ont':
        return " -ax map-ont "
    elif platform == 'pb':
        return " -ax map-pb "
    elif platform == 'hifi':
        return " -ax map-hifi "
    else:
        print("Error! platform error: {}".format(platform))
        1

rule minimap2_alignment:
    input:
        reads = "01_raw_data/{sample}.{platform}.fq.gz",
        ref = INDEX_REF
    output:
        "02_bam/{sample}.{platform}.minimap2.bam"
    threads: config['cores_ngmlr_map']
    resources:
        mem_mb = config['mem_ngmlr_map'],
    log: "logs/1.{sample}.minimap2.{platform}.log"
    params:
        extra_parm = lambda wildcards: get_minimap2_parms(wildcards)
    shell:
        "{MINIMAP2} {params.extra_parm} --MD -Y -t {threads} {input.ref} {input.reads} -R '@RG\\tID:4\\tSM:{output}\\tLB:library1\\tPL:illumina\\tPU:unit1' 2>>{log}  | {SAMTOOLS} view -h -o {output} --output-fmt BAM 2>>{log}"

rule ngmlr_alignment:
    input:
        reads = "01_raw_data/{sample}.{platform}.fq.gz",
        ref = INDEX_REF
    output:
        "02_bam/{sample}.{platform}.ngmlr.bam"
        #"02_bam/{sample}.ngmlr.sort.bam"
    threads: config['cores_ngmlr_map']
    resources:
        mem_mb = config['mem_ngmlr_map'],
    log: "logs/1.{sample}.ngmlr.{platform}.log"
    params:
        extra_parm = lambda wildcards: get_ngmlr_parms(wildcards)
    shell:
        "{NGMLR} --rg-id 4 --rg-sm '{output}' --rg-lb library1 --rg-pl illumina --rg-pu unit1 -t {threads} -r {input.ref} -q {input.reads} -o /dev/stdout {params.extra_parm} 2>>{log} | perl -alne 'print and next if /^\@/; next if $F[4]<0; next if length($F[9])!=length($F[10]); print' 2>>{log} | {SAMTOOLS} view -h -o {output} --output-fmt BAM 2>>{log}" # | {SAMTOOLS} sort --output-fmt BAM --threads 12 -o {output}"
        #"{NGMLR} -t {threads} -r {input.ref} -q {input.reads} -o /dev/stdout -x ont | perl -lne 'print and next if /^\@/; /^\S+\t\S+\t\S+\t\S+\t(\S+)/ or die; print if $1>=0' | {SAMTOOLS} view -h -o {output} --output-fmt BAM" # | {SAMTOOLS} sort --output-fmt BAM --threads 12 -o {output}"

rule sort_bam:
    input:
        lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        "02_bam/{sample}.{platform}.{mapper}.sort.bam"
    threads: 10
    resources:
        mem_mb = 10000
    shell:
        "{SAMTOOLS} sort -@ {threads} -O BAM -o {output} {input} "

rule index:
    input:
        #"02_bam/{sample}.{platform}.{mapper}.sort.bam"
        lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai"
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        "{SAMTOOLS} index {input} "

rule sniffles_call:
    input:
        #"02_bam/{sample}.{platform}.{mapper}.sort.bam",
        #"02_bam/{sample}.{platform}.{mapper}.sort.bam.bai"
        bam = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        vcf = "03_vcf/01_sniffles/{sample}.sniffles.vcf"
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.sniffles.log"
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        "{SNIFFLES} -t {threads} --input {input.bam} --vcf {output.vcf} 2>>{log}"


def get_cuteSV_call_parms(wildcards):
    platform = wildcards.platform
    if platform == 'ont':
        return " --min_read_len 500 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 "
    elif platform == 'pb':
        return " --min_read_len 500 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 "
    elif platform == 'hifi':
        return " --min_read_len 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 "
    else:
        print("Error! platform error: {}".format(platform))
        1


rule cuteSV_call:
    input:
        bam = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai".format( sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        "03_vcf/02_cuteSV/{sample}.{platform}.cuteSV.vcf"
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.cuteSV.{platform}.log"
    params:
        workdir = "03_vcf/02_cuteSV/{sample}.tmp",
        extra_parm = lambda wildcards: get_cuteSV_call_parms(wildcards)
    #shadow: "shallow"
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        """
        mkdir {params.workdir}
        {CUTESV} {input.bam} {INDEX_REF} {output} {params.workdir} --threads {threads} {params.extra_parm} --genotype 2>>{log}
        """



rule svim_call:
    input:
        #"02_bam/{sample}.{platform}.{mapper}.sort.bam",
        #"02_bam/{sample}.{platform}.{mapper}.sort.bam.bai"
        bam = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=S2T.loc[wildcards.sample, 'platform'], mapper=MAPPER),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai".format(sample=wildcards.sample, platform=S2T.loc[wildcards.sample, 'platform'], mapper=MAPPER)
    output:
        vcf1 = "03_vcf/03_svim/{sample}/variants.vcf",
        vcf2 = "03_vcf/03_svim/{sample}/pre-variants.vcf"
    log: "logs/2.{sample}.svim.log"
    params:
        DIR="03_vcf/03_svim/{sample}"
    threads: 1
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        """
        {SVIM} alignment {params.DIR} {input.bam} {INDEX_REF} --insertion_sequences 2>>{log} && \
        cp {output.vcf1} {output.vcf2}
        """
        #cat {output.vcf1} | grep -v '#' | awk '$6>=10' > {output.vcf2}

rule pbsv_align:
    input:
        bam = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        "02_bam/{sample}.pbsv.sort.bam",
    threads: 1
    log: "logs/2.{sample}.pbsv.1.picard.log"
    resources:
        mem_mb = 10000
    shell:
        "{PICARD} AddOrReplaceReadGroups I={input.bam} O={output} RGID=4 SORT_ORDER=coordinate RGSM={input.bam} RGLB=library1 RGPL=illumina RGPU=unit1 2>>{log}"

rule pbsv_discovery:
    input:
        #bam = "02_bam/{sample}.pbsv.sort.bam"
        bam = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER),
        bai = lambda wildcards: "02_bam/{sample}.{platform}.{mapper}.sort.bam.bai".format(sample=wildcards.sample, platform=get_platform(wildcards), mapper=MAPPER)
    output:
        outfile="03_vcf/04_pbsv/{sample}/pbsv.svsig.gz",
    threads: 1
    params:
        outdir="03_vcf/04_pbsv/{sample}"
    log: "logs/2.{sample}.pbsv.2.pbsv_discover.log"
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        """
        mkdir -p {params.outdir}
        {PBSV} discover {input.bam} {output.outfile} 2>>{log}
        """


rule pbsv_call:
    input:
        "03_vcf/04_pbsv/{sample}/pbsv.svsig.gz",
    output:
        "03_vcf/04_pbsv/{sample}/pbsv.vcf"
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.pbsv.3.pbsv_call.log"
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        "{PBSV} call --num-threads {threads} {INDEX_REF} {input} {output} 2>>{log}"



# nucmer + Assemblytics
rule Assemblytics1nucmer:
    input:
        ref = INDEX_REF,
        que = get_assb_fasta
    params:
        prefix = "03_vcf/05_Assemblytics/{sample}/1.nucmer",
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.Assemblytics.1.nucmer.log"
    output:
        delta = "03_vcf/05_Assemblytics/{sample}/1.nucmer.delta",
        dir_now = directory("03_vcf/05_Assemblytics/{sample}"),
        #dir_now = "03_vcf/05_Assemblytics/{sample}",
    resources:
        mem_mb = 20000
    shell:
        """
        mkdir -p {output.dir_now}
        {NUCMER} --maxmatch -c 500 -b 500 -l 100 -t {threads} {input.ref} {input.que} --prefix {params.prefix} >{log} 2>&1
        """

rule Assemblytics2filter:
    input:
        delta = "03_vcf/05_Assemblytics/{sample}/1.nucmer.delta",
    output:
        filterd = "03_vcf/05_Assemblytics/{sample}/2.nucmer.filter",
    log: "logs/2.{sample}.Assemblytics.2.filter.log"
    threads: 1
    shell:
        """
        {DELTAFILTER} -m -i 90 -l 100 {input.delta} > {output.filterd} 2>>{log}
        """
        
rule Assemblytics3showcoords:
    input:
        filterd = "03_vcf/05_Assemblytics/{sample}/2.nucmer.filter",
    output:
        coords = "03_vcf/05_Assemblytics/{sample}/3.nucmer.coords",
    log: "logs/2.{sample}.Assemblytics.3.showcoords.log"
    threads: 1
    shell:
        """
        {SHOWCOORDS} -THrd {input.filterd} > {output.coords} 2>>{log}
        """

rule Assemblytics4Assemblytics:
    input:
        filterd = "03_vcf/05_Assemblytics/{sample}/2.nucmer.filter",
        coords = "03_vcf/05_Assemblytics/{sample}/3.nucmer.coords",
    params:
        prefix = "03_vcf/05_Assemblytics/{sample}",
    threads: 1
    output:
        "03_vcf/05_Assemblytics/{sample}/4.Assemblytics.Assemblytics_structural_variants.bed"
    shell:
        """
        cd {params.prefix}; pwd; {ASSEMBLYTICS} ./2.nucmer.filter ./4.Assemblytics 1000 50 10000000; cd ..
        """

rule Assemblytics5parse:
    input:
        bed = "03_vcf/05_Assemblytics/{sample}/4.Assemblytics.Assemblytics_structural_variants.bed",
        ref = INDEX_REF,
        que = get_assb_fasta,
    output:
        vcf = "03_vcf/05_Assemblytics/{sample}/5.parse.vcf.gz",
    log: "logs/2.{sample}.Assemblytics.5.parse.log"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        perl {workflow.basedir}/scripts/long_assemblytics_parser.pl --in {input.bed} --out - --ref {input.ref} --query {input.que} 2>>{log} | {BCFTOOLS} sort -O z -o {output.vcf} >>{log} 2>&1 && \
        {TABIX} {output.vcf} 2>>{log}
        """

rule SVIM_asm_call_1_map:
    input:
        ref = INDEX_REF,
        que = get_assb_fasta
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.SVIM_asm.1.minimap2.log"
    output:
        bam = "03_vcf/06_SVIM_asm/{sample}/1.minimap.bam",
    resources:
        mem_mb = 20000
    shell:
        """
        {MINIMAP2} -a -x asm5 --cs -r2k -t {threads} {input.ref} {input.que} 2>>{log} | {SAMTOOLS} view -b -o {output.bam} - 2>>{log}
        """

rule SVIM_asm_call_2_sort:
    input:
        bam = "03_vcf/06_SVIM_asm/{sample}/1.minimap.bam",
    output:
        sort_bam = "03_vcf/06_SVIM_asm/{sample}/2.minimap.sort.bam",
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        """
        {SAMTOOLS} sort -@{threads} -m4G -O BAM -o {output.sort_bam} {input.bam}
        """

rule SVIM_asm_call_3_index:
    input:
        sort_bam = "03_vcf/06_SVIM_asm/{sample}/2.minimap.sort.bam",
    output:
        sort_bam_bai = "03_vcf/06_SVIM_asm/{sample}/2.minimap.sort.bam.bai",
    threads: 1
    log: "logs/2.{sample}.SVIM_asm.3.index.log"
    resources:
        mem_mb = 1000
    shell:
        """
        {SAMTOOLS} index {input.sort_bam} 2>>{log}
        """

rule SVIM_asm_call_4_haploidrun:
    input:
        sort_bam = "03_vcf/06_SVIM_asm/{sample}/2.minimap.sort.bam",
        sort_bam_bai = "03_vcf/06_SVIM_asm/{sample}/2.minimap.sort.bam.bai",
        ref = INDEX_REF,
    output:
        vcf = "03_vcf/06_SVIM_asm/{sample}/3.haploid/variants.vcf",
    threads: config['cores_sv3_call']
    log: "logs/2.{sample}.SVIM_asm.4.haploid.log"
    params:
        outdir = "03_vcf/06_SVIM_asm/{sample}/3.haploid",
    resources:
        mem_mb = config['mem_sv3_call']
    shell:
        """
        {SVIM_ASM} haploid {params.outdir} {input.sort_bam} {input.ref}  2>>{log}
        """

rule filter_parse:
    input:
        vcf = "03_vcf/{software}/{filename}.vcf",
        ref = INDEX_REF
    output:
        vcf = "03_vcf/{software}/{filename}.parse.vcf.gz",
        tbi = "03_vcf/{software}/{filename}.parse.vcf.gz.tbi",
        inv_bed = "03_vcf/{software}/{filename}.parse.inv.bed",
    threads: 1
    resources:
        mem_mb = 4000
    log: "logs/2.{software}.{filename}.parse.log"
    shell:
        """
        perl {workflow.basedir}/scripts/long_caller_parser.pl --in {input.vcf} --out - --out_inv {output.inv_bed} --ref {input.ref} 2>>{log} | {BCFTOOLS} sort -O z -o {output.vcf} >>{log} 2>&1 && \
        {TABIX} {output.vcf} 2>>{log}
        """


# rule consensus_SURVIVOR:
#     input:
#         "03_vcf/01_sniffles/{sample}.sniffles.vcf",
#         #"03_vcf/02_cuteSV/{sample}.cuteSV.vcf",
#         lambda wildcards: "03_vcf/02_cuteSV/{sample}.{platform}.cuteSV.vcf".format(sample=wildcards.sample, platform=S2T.loc[wildcards.sample, 'platform']),
#         "03_vcf/03_svim/{sample}/variants.vcf",
#         "03_vcf/04_pbsv/{sample}/pbsv.vcf"
#     output:
#         "04_consensus_vcf/{sample}.consensus_SURVIVOR.vcf"
#     params:
#         "03_vcf/05_sample_list/list_{sample}"
#     threads: 1
#     resources:
#         mem_mb = 4000
#     shell:
#         """
#         mkdir -p 03_vcf/05_sample_list/
#         ls {input[0]} {input[1]} {input[2]} {input[3]}> {params}
#         {SURVIVOR} merge {params} 1000 2 1 1 0 30 {output}
#         """




rule cal_genome_txt:
    # chr   length
    input:
        ref = INDEX_REF
    output:
        genome = INDEX_REF + ".genome"
    threads: 1
    run:
        chrlen = dict()
        chr_now = ''
        with open(input.ref) as f:
            for line in f:
                if line.startswith('>'):
                    chr_now = line.strip().split()[0][1:]
                    chrlen[chr_now] = 0
                else:
                    chrlen[chr_now] += len(line.strip())
        with open(output.genome, 'w') as f:
            for chr_now in chrlen:
                print(chr_now, chrlen[chr_now], sep='\t', file=f)
        

rule consensus_inv:
    input:
        #genome = INDEX_REF + ".genome",
        ref = INDEX_REF,
        bed1 = "03_vcf/01_sniffles/{sample}.sniffles.parse.inv.bed",
        bed2 = lambda wildcards: "03_vcf/02_cuteSV/{sample}.{platform}.cuteSV.parse.inv.bed".format(sample=wildcards.sample, platform=get_platform(wildcards)),
        bed3 = "03_vcf/03_svim/{sample}/variants.parse.inv.bed",
        bed4 = "03_vcf/04_pbsv/{sample}/pbsv.parse.inv.bed",
    output:
        vcf = "04_consensus_vcf/{sample}/00.inv.vcf.gz",
    threads: 1
    log: "logs/3.00.{sample}.consensus_inv.log"
    params:
        sample = "{sample}",
    resources:
        mem_mb = 1000,
    shell:
        """
        perl {workflow.basedir}/scripts/long_inv_bed2vcf.pl --ref {input.ref} --out {output.vcf} --sample {params.sample} {input.bed1} {input.bed2} {input.bed3} {input.bed4} >>{log} 2>&1
        """

def get_consensus_vcfs_to_merge(wildcards):
    sid = wildcards.sample
    fa = get_assb_fasta(wildcards)
    platform = get_platform(wildcards)
    vcfs = []
    # is in SV_CALLER_LIST (in lower case)
    SV_CALLER_LIST_lower = [x.lower() for x in config['SV_CALLER_LIST']]
    if 'sniffles' in SV_CALLER_LIST_lower:
        vcfs.append("03_vcf/01_sniffles/{}.sniffles.parse.vcf.gz".format(sid))
    if 'cutesv' in SV_CALLER_LIST_lower:
        vcfs.append("03_vcf/02_cuteSV/{}.{}.cuteSV.parse.vcf.gz".format(sid, platform))
    if 'svim' in SV_CALLER_LIST_lower:
        vcfs.append("03_vcf/03_svim/{}/variants.parse.vcf.gz".format(sid))
    if 'pbsv' in SV_CALLER_LIST_lower: 
        vcfs.append("03_vcf/04_pbsv/{}/pbsv.parse.vcf.gz".format(sid))
    if 'assemblytics' in SV_CALLER_LIST_lower and fa != '-':
        vcf = "03_vcf/05_Assemblytics/{}/5.parse.vcf.gz".format(sid)
        vcfs.append(vcf)
    if 'svim_asm' in SV_CALLER_LIST_lower and fa != '-':
        vcf = "03_vcf/06_SVIM_asm/{}/3.haploid/variants.parse.vcf.gz".format(sid)
        vcfs.append(vcf)
    return vcfs

rule consensus:
    input:
        get_consensus_vcfs_to_merge
    output:
        vcf = "04_consensus_vcf/{sample}/01.merged.vcf.gz",
    threads: 4
    resources:
        mem_mb = 4000,
    log: "logs/3.01.{sample}.consensus.log"
    params:
        dir_now = directory("04_consensus_vcf/{sample}"),
    shell:
        """
        mkdir -p {params.dir_now}
        {BCFTOOLS} merge -m none -O v -o {output.vcf} --threads {threads} {input} >>{log} 2>&1
        """
