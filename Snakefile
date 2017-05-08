
import numpy as np
import pandas as pd
import os
import sh

configfile: "config.yaml"

def determine_regions(size, step, overlap):
    
    x = np.arange(1, size, step - overlap)
    w = np.column_stack([x, x + step])[:-1]
    w[-1][-1] = size
    return w

# based on chromosome size work out pir filenames
def determine_shapeit_jobs(chromosome, step, overlap):
    l = []
    fstem = "phased/{chrom}/haplotypes_ag1000g_{start}_{stop}.gz"
    
    size = config["contiglen"][chromosome]
    win = determine_regions(size, step, overlap)
    for _s, _e in win:
        l.append(fstem.format(chrom=chromosome, start=_s, stop=_e))
    return l

# load bam locations
bam_locations = pd.read_csv(config["bam_loc_table"], sep="\t", index_col=0)

rule all:
    input:
        expand("shapeit/{chrom}_haplotypes_ag1000g.h5", chrom=config["chromosomes"])

rule convert_hdf5:
    input:
        haplotypes="shapeit/{chrom}_haplotypes_ag1000g.gz",
        samples="shapeit/{chrom}_samples_ag1000g.gz"
    output:
        hdf5="shapeit/{chrom}_haplotypes_ag1000g.h5"
    conda:
        "pir_env.yaml"
    params:
        req="h_vmem=4G,h=india.well.ox.ac.uk",
        chunk_size=100000,
        max_sites=40000000
    script:
        "../autosomal_phasing/shapeit_2_hdf5.py"

rule ligate_haplotypes:
    input: 
        phasedchunks=lambda y: determine_shapeit_jobs(
            y.chrom, config["chunksize"], config["chunkoverlap"]),
        vcf=config["vcfstem"]
    output:
        haplotypes="shapeit/{chrom}_haplotypes_ag1000g.gz",
        samples="shapeit/{chrom}_samples_ag1000g.gz",
        vcf=temp("shapeit/ag1000g.{chrom}.all.vcf.gz")
    params:
        req="h_vmem=16G"
    shell:
        "gunzip -c {input.vcf} | gzip -c > {output.vcf}; "
        "/home/njh/exec/bin/ligateHAPLOTYPES "
        "--vcf {output.vcf} "
        "--chunks {input.phasedchunks} "
        "--output {output.haplotypes} {output.samples}"
        
rule shape_it:
    input:
        pir="PIR/{chrom}/{start}_{stop}_split.pir.gz",
        vcf="_build/{chrom}/{start}_{stop}_split.vcf.gz",
        tbi="_build/{chrom}/{start}_{stop}_split.vcf.gz.tbi",
        map=config["map_path"]
    output:
        haplotypes="phased/{chrom}/haplotypes_ag1000g_{start}_{stop}.gz",
        samples="phased/{chrom}/samples_ag1000g_{start}_{stop}.gz",
        log="phased/{chrom}/logging_shapeit_{start}_{stop}.log"
    params:
        thread=1,
        window=0.5,
        ne=1000000,
        req="h_vmem=8G"
    shell:
        "shapeit -assemble "
        "--input-pir {input.pir} "
        "--output-log {output.log} "
        "--thread {params.thread} "
        "--noped --window {params.window} "
        "--effective-size {params.ne} "
        "--aligned "
        "--input-vcf {input.vcf} "
        "--input-map {input.map} "
        "--output-max {output.haplotypes} {output.samples}"

rule check_samples:
    input:
        config["vcfstem"]
    output:
        "_build/allsamples_{chrom}.txt"
    params:
        req="h_vmem=2G"
    shell:
        "vcfsamplenames {input} > {output}"
        
rule verify_bams:
    input:
        "_build/allsamples_{chrom}.txt"
    output:
        "_build/bams_{chrom}.ok"
    params:
        req="h_vmem=2G"
    run:
        with open(input[0], "r") as fh:
            samples = [l.rstrip() for l in fh.readlines()]
        
        for sid in samples:
            fn = bam_locations.loc[sid].path
            assert os.path.isfile(fn), "Not found. {0} : {1}".format(sid, fn)
        sh.touch(output[0])
        
rule make_bamlist:
    output:
        "_build/bamlist_{chrom}.txt"
    params:
        req="h_vmem=2G"
    run: 
        _df = bam_locations.copy()
        _df["chrom"] = wildcards.chrom
        _df.to_csv("_build/bamlist_{chrom}.txt".format(chrom=wildcards.chrom),
                   sep=" ", header=False, index=True)
        
rule extract_pirs:
    input:
        vcf="_build/{chrom}/{start}_{stop}_split.vcf.gz",
        tbi="_build/{chrom}/{start}_{stop}_split.vcf.gz.tbi",
        baml="_build/bamlist_{chrom}.txt",
        tabl="_build/allsamples_{chrom}.txt",
        okbm="_build/bams_{chrom}.ok"
    output:
        "PIR/{chrom}/{start}_{stop}_split.pir.gz"
    conda:
        "pir_env.yaml"
    params:
        prefix="PIR/{chrom}/{start}_{stop}_split.pir",
        req="h_vmem=2G"
    shell:
        "extractPIRs --vcf {input.vcf} "
        "--bam {input.baml} "
        "--out {params.prefix}; "
        "gzip {params.prefix}"
        
rule split_vcf:
    input:
        config["vcfstem"]
    output:
        vcf=temp("_build/{chrom}/{start}_{stop}_split.vcf.gz"),
        tbi=temp("_build/{chrom}/{start}_{stop}_split.vcf.gz.tbi")
    conda:
        "pir_env.yaml"
    params:
        req="h_vmem=2G"
    shell:
        "tabix -h {input} {wildcards.chrom}:{wildcards.start}-{wildcards.stop} | "
        "bgzip -c > {output.vcf}; "
        "tabix -p vcf {output.vcf}"
