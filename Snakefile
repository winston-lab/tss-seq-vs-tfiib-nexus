#!/usr/bin/env python
from math import log2

configfile: "config.yaml"

subworkflow tss_pipe:
    workdir: config["tss_workflow"]

subworkflow nexus_pipe:
    workdir: config["nexus_workflow"]

subworkflow tfiib_pipe:
    workdir: config["tfiib_workflow"]

subworkflow build_annotations:
    workdir: config["annotation_workflow"]

configfile: build_annotations("config.yaml")

CONTROLS = [v for k,v in config["comparisons"].items()]
CONDITIONS = [k for k,v in config["comparisons"].items()]
TSS_NORMS = ["libsizenorm"] + (["spikenorm"] if config["spike_in_normalization"]["tss"] else [])
TFIIB_NORMS = ["libsizenorm"] + (["spikenorm"] if config["spike_in_normalization"]["tfiib"] else [])
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

localrules:
    target,
    match_tfiib_to_tss,
    match_tss_to_tfiib

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule target:
    input:
        "config.yaml",
        expand(expand("matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{{category}}/{condition}-v-{control}_{{category}}-tss-seq-{{tss_norm}}-tfiib-{{tfiib_norm}}-matched-peaks-distance{{distance}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), category=CATEGORIES, tss_norm=TSS_NORMS, tfiib_norm=TFIIB_NORMS, distance=config["search_distance"]),
        expand(expand("matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{{category}}/{condition}-v-{control}_{{category}}-tfiib-{{tfiib_norm}}-tss-seq-{{tss_norm}}-matched-peaks-distance{{distance}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), category=["genic", "intragenic", "intergenic"], tss_norm=TSS_NORMS, tfiib_norm=TFIIB_NORMS, distance=config["search_distance"])

rule match_tfiib_to_tss:
    input:
        tss_results = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{category}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-{category}-all.tsv"),
        tfiib_results = nexus_pipe("diff_binding/peaks/{condition}-v-{control}/{tfiib_norm}/{condition}-v-{control}_tfiib-chipnexus-{tfiib_norm}-peaks-diffbind-results-all.tsv"),
        genome_fasta = config["genome"]["fasta"]
    output:
        "matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{category}/{condition}-v-{control}_{category}-tss-seq-{tss_norm}-tfiib-{tfiib_norm}-matched-peaks-distance{distance}.tsv"
    params:
        search_distance = config["search_distance"],
        feature_header = lambda wc: "" + \
                ("feature_chrom\tfeature_start\tfeature_end\tfeature_name\tfeature_strand\t" if wc.category != "intergenic" else "") + \
                ("feature_distance\t" if wc.category not in ["genic", "intergenic"] else ""),
        drop_columns = lambda wc: {"genic": "26,36,37,39",
                                   "intergenic": "31,32,34"}.get(wc.category, "26,37,38,40")
    shell: """
        tail -n +2 {input.tss_results} | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2+$15, $2+$15+1, $4, $5, $6, $0}}' | \
        bedtools slop -l {params.search_distance} -r 0 -s -i stdin -g <(faidx -i chromsizes {input.genome_fasta}) | \
        bedtools intersect -loj -a stdin -b <(tail -n +2 {input.tfiib_results}) | \
        cut --complement -f1-6,15,16,18,{params.drop_columns} | \
        cat <(echo -e "tss_chrom\ttss_start\ttss_end\ttss_name\ttss_score\ttss_strand\ttss_lfc\ttss_lfc_SE\ttss_FDR\ttss_expr_condition\ttss_expr_control\ttss_summit\t{params.feature_header}tfiib_chrom\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_score\ttfiib_strand\ttfiib_lfc\ttfiib_lfc_SE\ttfiib_FDR\ttfiib_abundance_condition\ttfiib_abundance_control") - > {output}
        """

rule match_tss_to_tfiib:
    input:
        tfiib_results = tfiib_pipe("diff_binding/peaks/{condition}-v-{control}/{tfiib_norm}/{category}/{condition}-v-{control}_tfiib-chipnexus-{tfiib_norm}-peaks-diffbind-results-{category}-all.tsv"),
        tss_results = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-all.tsv"),
        genome_fasta = config["genome"]["fasta"]
    output:
        "matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{category}/{condition}-v-{control}_{category}-tfiib-{tfiib_norm}-tss-seq-{tss_norm}-matched-peaks-distance{distance}.tsv"
    params:
        search_distance = config["search_distance"],
        feature_header = lambda wc: "" + \
                ("feature_chrom\tfeature_start\tfeature_end\tfeature_name\tfeature_strand\t" if wc.category != "intergenic" else "") + \
                ("feature_distance\t" if wc.category not in ["genic", "intergenic"] else ""),
        drop_columns = lambda wc: {"genic": "26,36,37,39",
                                   "intergenic": "31,32,34"}.get(wc.category, "26,37,38,40")
    shell: """
        tail -n +2 {input.tfiib_results} | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2+$15, $2+$15+1, $4, $5, $6, $0}}' | \
        bedtools slop -l {params.search_distance} -r {params.search_distance} -i stdin -g <(faidx -i chromsizes {input.genome_fasta}) | \
        bedtools intersect -loj -a stdin -b <(tail -n +2 {input.tss_results}) | \
        cut --complement -f1-6,15,16,18,{params.drop_columns} | \
        cat <(echo -e "tfiib_chrom\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_score\ttfiib_strand\ttfiib_lfc\ttfiib_lfc_SE\ttfiib_FDR\ttfiib_abundance_condition\ttfiib_abundance_control\ttfiib_summit\t{params.feature_header}tss_chrom\ttss_start\ttss_end\ttss_name\ttss_score\ttss_strand\ttss_lfc\ttss_lfc_SE\ttss_FDR\ttss_expr_condition\ttss_expr_control") - > {output}
        """

