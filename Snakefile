#!/usr/bin/env python

configfile: "config.yaml"

CONDITIONS = config["conditions"] + config["controls"]
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

localrules: all,
    make_peak_tables, cat_peaks, plot_peaks

rule all:
    input:
        expand("tables/{condition}-tfiib-upstr-tss-{category}.tsv", condition=CONDITIONS, category=CATEGORIES+["all"]),
        "tables/all-tfiib-upstr-tss.tsv",
        "figures/peaks-tss-v-tfiib-overlap-mosaic.svg",
        expand(expand("diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-{{category}}.tsv", zip, condition=config["conditions"], control=config["controls"]), category=CATEGORIES+["all"]),
        "diffexp_tables/all-diffexp-tfiib-upstr-tss.tsv",
        "figures/tss-v-tfiib-diffexp-scatter.svg",
        expand("tables/bed/{condition}_{category}-tss-{zz}-tfiib-upstr.bed", condition=CONDITIONS, category=CATEGORIES+["all"], zz = ["with", "without"]),
        expand("figures/heatmaps/{condition}-both.tsv.gz", condition=CONDITIONS),
        expand("figures/heatmaps/{condition}-heatmap.svg", condition=CONDITIONS)

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

rule make_diffexp_tables_all:
    input:
        tss = lambda wildcards: config["tss"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-results-" + config["tss"]["norm"] + "-all.tsv",
        nexus = lambda wildcards: config["nexus"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-tfiib-chipnexus-results-" + config["nexus"]["norm"] + "-all.tsv",
        chrsizes = config["genome"]["chrsizes"]
    params: upstr = config["tss-upstr-dist"]
    output:
        "diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-all.tsv"
    shell: """
        tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $11, $12, $13, "na", "na", "na"}}' | sort -k1,1 -k2,2n | paste <(tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3}}' | sort -k1,1 -k2,2n | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | bedtools intersect -a stdin -b <(tail -n +2 {input.nexus} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n) -wao | cut --complement -f1-3,8,18,22,23,26-28 | awk -v condition={wildcards.condition} -v control={wildcards.control} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, control, "all"}}' > {output}
        """

rule make_diffexp_tables_category:
    input:
        tss = lambda wildcards: config["tss"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.category + "/" + wildcards.condition + "-v-" + wildcards.control + "-results-" + config["nexus"]["norm"] + "-all-" + wildcards.category + ".tsv",
        nexus = lambda wildcards: config["nexus"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-tfiib-chipnexus-results-" + config["nexus"]["norm"] + "-all.tsv",
        chrsizes = config["genome"]["chrsizes"]
    params:
        upstr = config["tss-upstr-dist"]
    output:
        "diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-{category}.tsv"
    shell: """
        tail -n +2 {input.tss} | awk -v ttype={wildcards.category} 'BEGIN{{FS=OFS="\t"}} ttype=="intergenic"{{$13="na"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n | paste <(tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3}}' | sort -k1,1 -k2,2n | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | bedtools intersect -a stdin -b <(tail -n +2 {input.nexus} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n) -wao | cut --complement -f1-3,8,18,22,23,26-28 | awk -v category={wildcards.category} -v condition={wildcards.condition} -v control={wildcards.control} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, control, category}}' > {output}
        """

rule cat_diffexp:
    input:
        expand(expand("diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-{{category}}.tsv", zip, condition=config["conditions"], control=config["controls"]), category=CATEGORIES+["all"])
    output:
        "diffexp_tables/all-diffexp-tfiib-upstr-tss.tsv"
    shell: """
        cat <(echo -e "chrom\ttss_start\ttss_end\ttss_name\tstrand\ttss_meanExpr\ttss_lfc\ttss_logpadj\ttss_cond\ttss_ctrl\tfeat_start\tfeat_end\tfeat_name\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_meanExpr\ttfiib_lfc\ttfiib_logpadj\ttfiib_cond\ttfiib_ctrl\toverlap\tcondition\tcontrol\tcategory") {input} > {output}
        """

rule plot_diffexp:
    input:
        "diffexp_tables/all-diffexp-tfiib-upstr-tss.tsv"
    output:
        scatter = "figures/tss-v-tfiib-diffexp-scatter.svg"
    script: "scripts/tss-v-tfiib-diffexp.R"

rule get_bed:
    input:
        "tables/{condition}-tfiib-upstr-tss-{category}.tsv"
    output:
        yes = "tables/bed/{condition}_{category}-tss-with-tfiib-upstr.bed",
        no = "tables/bed/{condition}_{category}-tss-without-tfiib-upstr.bed",
    script: "scripts/getbed.R"

rule compute_matrix:
    input:
        annotation = "tables/bed/{condition}-tss-{class}-tfiib-upstr.bed",
        bw = lambda wildcards: config["nexus"]["coverage"] + config["nexus"]["norm"] + "/bw/" + wildcards.condition + "-1-tfiib-chipnexus-" + config["nexus"]["norm"] + "-qfrags.bw",
    output:
        dtfile = temp("figures/heatmaps/{condition}-{class}.mat"),
        matrix = temp("figures/heatmaps/{condition}-{class}.tsv"),
        matrix_gz = temp("figures/heatmaps/{condition}-{class}.tsv.gz"),
    params:
        refpoint = "TSS",
        upstream = config["heatmaps"]["upstream"],
        dnstream = config["heatmaps"]["downstream"],
        binsize = config["heatmaps"]["binsize"],
        sort = "keep",
        binstat = config["heatmaps"]["binstat"]
    threads : config["threads"]
    log: "logs/compute_matrix-{condition}-{class}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}
        """

rule melt_matrix:
    input:
        matrix = "figures/heatmaps/{condition}-{class}.tsv.gz"
    output:
        temp("figures/heatmaps/{condition}-{class}-melted.tsv.gz")
    params:
        binsize = config["heatmaps"]["binsize"],
        upstream = config["heatmaps"]["upstream"],
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("figures/heatmaps/{{condition}}-{zz}-melted.tsv.gz", zz=["with","without"])
    output:
        "figures/heatmaps/{condition}-both.tsv.gz"
    shell: """
        cat {input} > {output}
        """

rule plot_heatmaps:
    input:
        "figures/heatmaps/{condition}-both.tsv.gz"
    output:
        "figures/heatmaps/{condition}-heatmap.svg"
    params:
        cutoff=config["heatmaps"]["pct_cutoff"],
        upstream=config["heatmaps"]["upstream"],
        dnstream=config["heatmaps"]["downstream"],
        cmap=config["heatmaps"]["colormap"]
    script: "scripts/plotHeatmaps.R"
