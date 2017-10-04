#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES = config["samples"]

rule all:
    input:
        expand("tss-v-tfiib-{norm}.png", norm=["libsizenorm","spikenorm"])

rule build_tss_regions:
    input:
        transcripts = config["genome"]["transcripts"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        "tss-regions.bed"
    params:
        windowsize = config["tss-windowsize"]
    shell: """
        bash scripts/makeStrandedBed.sh {input.transcripts} | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}}$6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools slop -b {params.windowsize} -g {input.chrsizes} -s -i stdin | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

#build BED file of regions upstream of TSS peaks to map TFIIB signal to
rule build_tfiib_regions:
    input:
        transcripts = config["genome"]["transcripts"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        "tfiib-regions.bed"
    params:
        upstr = config["tfiib-upstream"]
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}}$6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' {input.transcripts} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools slop -l {params.upstr} -r 0 -s -g {input.chrsizes} -i stdin | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule map_tss:
    input:
        bed = lambda wildcards: "tss-regions.bed" if wildcards.assay=="tss" else "tfiib-regions.bed",
        bg = lambda wildcards: config["tss-dir"][wildcards.norm] + SAMPLES[wildcards.sample]["tss"][wildcards.norm] if wildcards.assay=="tss" else config["tfiib-dir"][wildcards.norm] + SAMPLES[wildcards.sample]["tfiib"][wildcards.norm]
    output:
        temp("data/{sample}-{norm}-{assay}.tsv")
    shell: """
        bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | cut -f4,7 | sort -k1,1 | cat <(echo -e "gene\t{wildcards.sample}") - > {output}
        """

# rule map_tfiib:
#     input:
#         bed = "tfiib-regions.bed",
#         bg = lambda wildcards: config["tfiib-dir"][wildcards.norm] + SAMPLES[wildcards.sample]["tfiib"][wildcards.norm]
#     output:
#         temp("data/{sample}-{norm}-tfiib.tsv")
#     shell: """
#         bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | cut -f4,7 | sort -k1,1 | cat <(echo -e "gene\t{wildcards.sample}") - > {output}
#         """

rule paste_data:
    input:
        expand("data/{sample}-{{norm}}-{{assay}}.tsv", sample=SAMPLES)
    output:
        "data/allsamples-{norm}-{assay}.tsv"
    params:
        n = 2*len(SAMPLES)
    shell: """
        paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) > {output}
        """

rule plot:
    input:
        tss = "data/allsamples-{norm}-tss.tsv",
        tfiib = "data/allsamples-{norm}-tfiib.tsv"
    output:
        plot = "tss-v-tfiib-{norm}.png"
    script:
        "scripts/tss-v-tfiib.R"
