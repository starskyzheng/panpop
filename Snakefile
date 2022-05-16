import os
configfile: "config.yaml"
include: "subworkflows/util.py"

# parse config values
GRAPH=config['graph']
MAPPER=config['mapper']
BCFTOOLS = os.path.abspath(config['bcftools'])
VG = os.path.abspath(config["vg"])
TABIX = os.path.abspath(config["tabix"])
ZWORKDIR = config['workdir']
ZTMPDIR = 'tmp'

workdir: ZWORKDIR

include: "subworkflows/callSV.py"


SAMPLES = list_prepaire_files(config['sample_reads_list_file'])


if config['split_chr']==True:
    CHRS = gfa2chrs( GRAPH + '.gfa')
    if config['mode'] == 'genotype':
        rule all:
            input:
                expand('4.realign/4.filter_maf1.{chrm}.vcf.gz', chrm=CHRS),
                '5.final_result/1.final_mergechr.all.vcf.gz',
                '5.final_result/2.final_mergechr.sv.vcf.gz',
    elif config['mode'] == 'augment':
        rule all:
            input:
                #expand('4.realign/4.filter_maf1.{chrm}.vcf.gz', chrm=CHRS),
                '9.aug_final_result/1.final_mergechr.all.vcf.gz',
                '9.aug_final_result/2.final_mergechr.sv.vcf.gz',
elif config['split_chr']==False:
    CHRS=['all']
    if config['mode'] == 'genotype':
        rule all:
            input:
                '5.final_result/1.final.all.all.vcf.gz',
                '5.final_result/2.final.all.sv.vcf.gz'
    elif config['mode'] == 'augment':
        rule all:
            input:
                '9.aug_final_result/1.final.all.all.vcf.gz',
                '9.aug_final_result/2.final.all.sv.vcf.gz'

include: "subworkflows/panpopbase.py"
include: "subworkflows/panpopaug.py"



for dir in [ZTMPDIR]:
    if not os.path.exists(dir):
        os.makedirs(dir)


rule cleanall:
    params:
        rmfiles=expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'gcsa', 'gcsa.lcp', 'ids.mapping', 'dist', 'trivial.snarls', 'gfa.fa', 'min', 'snarls', 'giraffe.gbz', 'pg']) +
        [ZTMPDIR] +
        ["2.callSV", "3.merge_rawvcf", "4.realign", "5.final_result" ] +
        ["6.aug_dp", "7.aug_merge_rawvcf", "8.aug_realign", "9.aug_final_result"] +
        ["2.callSV.DPinfos.txt"] +
        ["logs"]
    shell:
        "rm -rf {params.rmfiles}"

# call variants in the graph
rule callSV:
    input: 
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.ext.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ),
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.call.ext.vcf.gz.tbi', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ),
        #expand('2.callSV/{sample}/{sample}-{graph}.{map}.gam.avgdp', sample=SAMPLES, graph=GRAPH, map=MAPPER)
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.q{minq}.pack.DPinfo', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ)

# call variants in the graph
rule callSVaug:
    input: 
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.ext.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ),
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.ext.vcf.gz.tbi', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ),
        #expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.gam.avgdp', sample=SAMPLES, graph=GRAPH, map=MAPPER)
        expand('2.callSV/{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack.DPinfo', sample=SAMPLES, graph=GRAPH, map=MAPPER, minq=MINQ)

