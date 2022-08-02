
#
# Graph construction
#
rule index_init:
    input:
        gfa='{genome}.gfa',
        #gfa='{genome}.rgfa',
    output: '{genome}.xg', '{genome}.gcsa', '{genome}.gcsa.lcp'
    threads: config['cores_xg']
    params:
        prefix='{genome}'
    resources:
        mem_mb=config['mem_xg']
    log: '2.callSV/logs/{genome}-autoindex-init.log.txt'
    shell:
        '{VG} autoindex --tmp-dir {ZTMPDIR} --workflow map --gfa {input.gfa} --prefix {params.prefix} -t {threads} >> {log} 2>&1'


rule index_giraffe:
    input:
        gfa='{genome}.gfa',
        #gfa='{genome}.rgfa',
    output: '{genome}.dist', '{genome}.giraffe.gbz', '{genome}.min'
    threads: config['cores_xg']
    params:
        prefix='{genome}'
    resources:
        mem_mb=config['mem_xg']
    log: '2.callSV/logs/{genome}-autoindex-giraffe.log.txt'
    shell:
        '{VG} autoindex --tmp-dir {ZTMPDIR} --workflow giraffe --gfa {input.gfa} --prefix {params.prefix} -t {threads} >> {log} 2>&1'


rule vg2gfa:
    input:
        '{graph}.xg'
    output:
        '{graph}.xg.gfa'
    shell:
        '{VG} view --gfa {input} > {output}'


rule gfa2fa:
    input: '{genome}.gfa'
    output: '{genome}.gfa.fa'
    threads: 1
    log: '2.callSV/logs/{genome}-gfa2fa.log.txt'
    shell:
        'perl {workflow.basedir}/scripts/gfa2fa.pl {input} {output} 2>> {log}'

# prepare the snarls index
rule index_snarls:
    input: '{genome}.xg'
    output: '{genome}.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    log: '2.callSV/logs/{genome}-snarls.log.txt'
    shell:
        '{VG} snarls -t {threads} {input} > {output} 2>> {log}'



rule index_trivial_snarls:
    input: '{genome}.xg'
    output: '{genome}.trivial.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    log: '2.callSV/logs/{genome}-trivialsnarls.log.txt'
    shell:
        '{VG} snarls -t {threads} --include-trivial {input} > {output} 2>> {log}'
