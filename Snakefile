#!/usr/bin/env python

configfile: "config.yaml"

CONDITIONS = config["conditions"] + config["controls"]
# COMPARISONS = config["comparisons"]
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

localrules: all,
    make_peak_tables, cat_peaks, plot_peaks


rule all:
    input:
        expand("tables/{condition}-tfiib-upstr-tss-{category}.tsv", condition=CONDITIONS, category=CATEGORIES+["all"]),
        "tables/all-tfiib-upstr-tss.tsv",
        "figures/peaks-tss-v-tfiib-overlap-mosaic.svg",

rule make_peak_tables:
    input:
        tss_peaks = lambda wildcards: config["tss"]["peaks-path"] + wildcards.condition + "-exp-idrpeaks.narrowPeak" if wildcards.category=="all" else config["tss"]["peaks-path"] + wildcards.category + "/" + wildcards.condition + "-exp-idrpeaks-" + wildcards.category + ".tsv",
        tss_coverage = lambda wildcards: config["tss"]["coverage"] + config["tss"]["norm"] + "/" + wildcards.condition + "-1-tss-" + config["tss"]["norm"] + "-SENSE.bedgraph",
        tfiib_peaks = lambda wildcards: config["nexus"]["peaks-path"] + wildcards.condition + "-" + config["genome"]["prefix"] + "_peaks.narrowPeak",
        tfiib_coverage = lambda wildcards: config["nexus"]["coverage"] + config["nexus"]["norm"] + "/" + wildcards.condition + "-1-tfiib-chipnexus-" + config["nexus"]["norm"] + "-combined.bedgraph",
        chrsizes = config["genome"]["chrsizes"]
    params: upstr = config["tss-upstr-dist"]
    output:
        "tables/{condition}-tfiib-upstr-tss-{category}.tsv"
    shell: """
        cut -f1-10 {input.tss_peaks} | awk 'BEGIN{{FS=OFS="\t"}}$6=="+"{{$1=$1"-plus"; print $0}}$6=="-"{{$1=$1"-minus"; print $0}}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.tss_coverage} -c 4 -o sum | sed "s/-minus//g;s/-plus//g" | sort -k1,1 -k2,2n | paste <(sort -k1,1 -k2,2n {input.tss_peaks} | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.tfiib_coverage} -c 4 -o sum | bedtools intersect -a stdin -b {input.tfiib_peaks} -wao | cut --complement -f1-3,7,8,10,11,16,19,20,21,23 | awk -v condition={wildcards.condition} -v category={wildcards.category} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, category}}' > {output}
        """

rule cat_peaks:
    input:
        expand("tables/{condition}-tfiib-upstr-tss-{category}.tsv", condition=CONDITIONS, category=CATEGORIES+["all"])
    output:
        "tables/all-tfiib-upstr-tss.tsv"
    shell: """
        cat <(echo -e "chrom\ttss_start\ttss_end\tstrand\ttss_idr\ttss_summit\ttss_expression\ttfiib_levels\ttfiib_start\ttfiib_end\ttfiib_enrichment\ttfiib_qval\ttfiib_summit\toverlap\tcondition\tcategory") {input} > {output}
        """

rule plot_peaks:
    input:
        "tables/all-tfiib-upstr-tss.tsv"
    params:
        dist = config["tss-upstr-dist"]
    output:
        mosaic = "figures/peaks-tss-v-tfiib-overlap-mosaic.svg",
        tss_expr = "figures/peaks-tss-v-tfiib-tss-expression.svg",
        tss_size = "figures/peaks-tss-v-tfiib-tss-sizes.svg",
        tss_idr = "figures/peaks-tss-v-tfiib-tss-idr.svg",
        tfiib_lvl = "figures/peaks-tss-v-tfiib-tfiib-levels.svg",
        contour = "figures/peaks-tss-v-tfiib-expression-contourplot.svg",
        scatter = "figures/peaks-tss-v-tfiib-expression-scatter.svg"
    script: "scripts/tss-v-tfiib.R"

