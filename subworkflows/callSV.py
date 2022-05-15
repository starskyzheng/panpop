# Needs:: tabix bcftools vg:v1.35.0
import os

include: "gfa_construct_vg.py"



def list_prepaire_files(file):
    sids = []
    #svdir = zworkdir + '/2.callSV/'
    svdir = '2.callSV/'
    # create dir if not exists
    if not os.path.exists(svdir):
        os.makedirs(svdir)
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'): continue
            sid = line.split()[0]
            r1_ori = line.split()[1]
            r1_ori = os.path.abspath(r1_ori)
            r2_ori = line.split()[2]
            r2_ori = os.path.abspath(r2_ori)
            dirnow = svdir + sid + '/'
            if not os.path.exists(dirnow):
                os.makedirs(dirnow)
            # soft link
            r1 = dirnow + sid + '_1.fastq.gz'
            r2 = dirnow + sid + '_2.fastq.gz'
            if not os.path.exists(r1):
                os.symlink(r1_ori, r1)
            if not os.path.exists(r2):
                os.symlink(r2_ori, r2)
            sids.append(sid)
    return sids



##
## Config
##

# parse config values
GRAPH=config['graph']
MAPPER=config['mapper']

MINQ=config['MAP_MINQ']
if MAPPER=='GraphAligner': MINQ=0

# BASEDIR is '2.callSV'


##
## Main rules
##

# indexes used by the default mapper and variant caller
rule construct_all:
    input: expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'gcsa'])

# indexes used by the giraffe mapper and variant caller
rule construct_all_gaffe:
    input:
        expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'dist', 'min', 'dist','giraffe.gbz']),




#
# Map reads from a sample and call variants
#
read1_in = '2.callSV/{sample}/{sample}_1.fastq.gz'
read2_in = '2.callSV/{sample}/{sample}_2.fastq.gz'
map_lab = '{sample}'
map_out = '2.callSV/{sample}/{sample}-{genome}'



# map reads to the graph
rule map:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}.xg',
        gcsa='{genome}.gcsa',
        gcsalcp='{genome}.gcsa.lcp'
    output: map_out + '.map.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    log: '2.callSV/logs/' + map_lab + '-{genome}-map.log.txt'
    run:
        shell("{VG} map -t {threads} -x {input.xg} -g {input.gcsa} -f {input.r1} -f {input.r2} > {output} 2>> {log}")

# map reads to the graph using mpmap in single-path mode
rule map_mpmap:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}.xg',
        gcsa='{genome}.gcsa',
        gcsalcp='{genome}.gcsa.lcp'
    output: map_out + '.mpmap.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    log: '2.callSV/logs/' + map_lab + '-{genome}-mpmap.log.txt'
    run:
        shell("{VG} mpmap -S -t {threads} -x {input.xg} -g {input.gcsa} -f {input.r1} -f {input.r2} > {output} 2>> {log}")

# map reads to the graph using giraffe
rule map_gaffe:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}.xg',
        min='{genome}.min',
        dist='{genome}.dist',
        gbz='{genome}.giraffe.gbz'
    output: map_out + '.gaffe.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    log: '2.callSV/logs/' + map_lab + '-{genome}-gaffe.log.txt'
    run:
        shell("{VG} giraffe --sample {wildcards.sample} -Z {input.gbz} -m {input.min} -d {input.dist} -f {input.r1} -f {input.r2} -t {threads} > {output} 2>> {log}")

rule map_GraphAligner:
    input:
        r1=read1_in,
        #vg='{genome}.vg'
        gfa='{genome}.xg.gfa'
    output: map_out + '.GraphAligner.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    log: '2.callSV/logs/' + map_lab + '-{genome}-GraphAligner.log.txt'
    run:
        shell("GraphAligner -g {input.gfa} -t {threads} -f {input.r1} -a {output} -x vg  2>> {log}")


# compute packed coverage from aligned reads
rule pack:
    input:
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.gam',
        xg='{graph}.xg'
    output: '2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack'
    threads: config['cores_pack']
    resources:
        mem_mb=config['mem_pack']
    log: '2.callSV/logs/{sample}-{graph}-{map}-q{minq}-pack.log.txt'
    run:
        shell("{VG} pack -x {input.xg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2>> {log}")


# cal avg depth of gam
rule cal_avg_depth:
    input:
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.gam',
        xg='{graph}.xg'
    output:
        avgdp='2.callSV/{sample}/{sample}-{graph}.{map}.gam.avgdp'
    log: '2.callSV/logs/{sample}-{graph}-{map}-avgdp.log.txt'
    threads: config['cores_avgdp']
    resources:
        mem_mb=config['mem_avgdp']
    run:
        shell("{VG} depth --threads 4 -g {input.gam} {input.xg} > {output.avgdp} 2>> {log}")

rule cal_avg_depth2:
    input:
        pack='2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack',
        xg='{graph}.xg'
    output:
        DPinfo='2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack.DPinfo'
    log: '2.callSV/logs/{sample}-{graph}-{map}.q{minq}.pack.DPinfo.log.txt'
    threads: config['cores_avgdp']
    run:
        shell("perl {workflow.basedir}/scripts/cal_dp_from_pack.pl {input.pack} {input.xg} {output.DPinfo} {log} {wildcards.sample} 1 2>> {log}")


# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack',
        xg='{graph}.xg',
        snarls='{graph}.snarls'
    output:
        vcf_ext="2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.ext.vcf.gz",
        vcf='2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.vcf.gz',
        idx='2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.vcf.gz.tbi'
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    log: '2.callSV/logs/{sample}-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("{VG} call -k {input.pack} -t {threads} -s {wildcards.sample} --genotype-snarls --snarls {input.snarls} {input.xg} 2>> {log} | {BCFTOOLS} sort 2>> {log} | bgzip --threads {threads} > {output.vcf_ext}")
        shell("{TABIX} -f -p vcf {output.vcf_ext}")
        shell("{BCFTOOLS} view -e 'GT=\"0/0\" || GT=\"./.\"' {output.vcf_ext} | {BCFTOOLS} filter -i 'FILTER==\"PASS\"' | bgzip --threads {threads} > {output.vcf}")
        shell("{TABIX} -f -p vcf {output.vcf}")


#
# Augment graph and call variants
#

# convert a XG index to a PG index
rule pgconvert:
    input: '{graph}.xg'
    output: '{graph}.pg'
    threads: 1
    resources:
        mem_mb=config['mem_pgconvert']
    log: '2.callSV/logs/{graph}-pgconvert.log.txt'
    run:
        shell('{VG} convert {input} -p > {output} 2>> {log}')

# augment a graph with aligned reads
rule augment:
    input:
        pg='{graph}.pg',
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.gam'
    output:
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.aug.gam',
        pg='2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg'
    threads: config['cores_augment']
    resources:
        mem_mb=config['mem_augment']
    log: '2.callSV/logs/{sample}-{graph}-{map}-augment.log.txt'             
    run:
        shell('{VG} augment {input.pg} {input.gam} -t {threads} -m 4 -q 5 -Q 5 -A {output.gam} > {output.pg} 2>> {log}')

# cal avg depth of gam
rule cal_avg_depth_aug:
    input:
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.aug.gam',
        xg='2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg'
    output:
        avgdp='2.callSV/{sample}/{sample}-{graph}.{map}.aug.gam.avgdp'
    log: '2.callSV/logs/{sample}-{graph}-{map}-aug-avgdp.log.txt'
    threads: config['cores_avgdp']
    run:
        shell("{VG} depth -g {input.gam} {input.xg} > {output.avgdp} 2>> {log}")

rule cal_avg_depth_aug2:
    input:
        pack='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack',
        xg='2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg'
    output:
        DPinfo='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack.DPinfo'
    log: '2.callSV/logs/{sample}-{graph}-{map}-aug.q{minq}.pack.DPinfo.log.txt'
    threads: config['cores_avgdp']
    run:
        shell("perl {workflow.basedir}/scripts/cal_dp_from_pack.pl {input.pack} {input.xg} {output.DPinfo} {log} {wildcards.sample} 2>> {log}")



# prepare the snarls index for the augmented graph
rule index_snarls_aug:
    input: '2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '2.callSV/{sample}/{sample}-{graph}.{map}.aug.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    log: '2.callSV/logs/{sample}-{graph}-{map}-aug-snarls.log.txt'
    run:
        shell('{VG} snarls -t {threads} {input} > {output} 2>> {log}')

# compute packed coverage on augmented graph
rule pack_aug:
    input:
        gam='2.callSV/{sample}/{sample}-{graph}.{map}.aug.gam',
        pg='2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack'
    threads: config['cores_pack']
    resources:
        mem_mb=config['mem_pack']   
    log: '2.callSV/logs/{sample}-{graph}-{map}-aug-q{minq}-pack.log.txt'
    run:
        shell("{VG} pack -x {input.pg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2>> {log}")


# call variants from the packed read coverage
rule call_aug:
    input:
        pack='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack',
        pg='2.callSV/{sample}/{sample}-{graph}.{map}.aug.pg',
        snarls='2.callSV/{sample}/{sample}-{graph}.{map}.aug.snarls'
    output:
        vcf_ext='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.ext.vcf.gz',
        vcf='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz',
        idx='2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz.tbi'
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    log: '2.callSV/logs/{sample}-{graph}-{map}-aug-q{minq}-call.log.txt'
    run:
        shell("{VG} call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.pg}  2>> {log} | {BCFTOOLS} sort 2>> {log} | bgzip --threads {threads} > {output.vcf_ext}")
        shell("{TABIX} -f -p vcf {output.vcf_ext}")
        shell("{BCFTOOLS} view -e 'GT=\"0/0\" || GT=\"./.\"' {output.vcf_ext} | {BCFTOOLS} filter -i 'FILTER==\"PASS\"' | bgzip --threads {threads} > {output.vcf}")
        shell("{TABIX} -f -p vcf {output.vcf}")

