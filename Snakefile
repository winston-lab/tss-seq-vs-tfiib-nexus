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
configfile: nexus_pipe("config.yaml")

CONTROLS = [v for k,v in config["tss_vs_tfiib_comparisons"].items()]
CONDITIONS = [k for k,v in config["tss_vs_tfiib_comparisons"].items()]
TSS_NORMS = ["libsizenorm"] + (["spikenorm"] if config["spike_in_normalization"]["tss"] else [])
TFIIB_NORMS = ["libsizenorm"] + (["spikenorm"] if config["spike_in_normalization"]["tfiib"] else [])
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

norm_sample_dict = {"libsizenorm" : PASSING,
                    "spikenorm" : SIPASSING}

def get_samples(norm, groups):
    return([k for k,v in norm_sample_dict[norm].items() if v["group"] in groups])

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys())),

localrules:
    target,
    match_tfiib_to_tss,
    match_tss_to_tfiib,
    create_tfiib_windows,
    map_counts_to_windows,
    combine_window_counts,
    differential_binding,
    reclassify_windows,
    datavis

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule target:
    input:
        "config.yaml",
        expand(expand("matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{{category}}/{condition}-v-{control}_{{category}}-tss-seq-{{tss_norm}}-tfiib-{{tfiib_norm}}-matched-peaks-distance{{distance}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), category=CATEGORIES, tss_norm=TSS_NORMS, tfiib_norm=TFIIB_NORMS, distance=int(config["search_distance"])),
        expand(expand("matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{{category}}/{condition}-v-{control}_{{category}}-tfiib-{{tfiib_norm}}-tss-seq-{{tss_norm}}-matched-peaks-distance{{distance}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), category=["genic", "intragenic", "intergenic"], tss_norm=TSS_NORMS, tfiib_norm=TFIIB_NORMS, distance=int(config["search_distance"])),
        expand(expand("window_diff_binding/{condition}-v-{control}/tss-{{tss_norm}}/{condition}-v-{control}_tss-{{tss_norm}}-window-{{window_size}}-tfiib-chipnexus-{{norm}}-diffbind-results.tsv", zip, condition=CONDITIONS, control=CONTROLS), tss_norm=TSS_NORMS, norm=TFIIB_NORMS, window_size=int(config["window_size"])),
        expand(expand("window_diff_binding/{condition}-v-{control}/tss-{{tss_norm}}/{condition}-v-{control}_tss-{{tss_norm}}-window-{{window_size}}-tfiib-chipnexus-{{norm}}-diffbind-results-{{category}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), tss_norm=TSS_NORMS, norm=TFIIB_NORMS, window_size=int(config["window_size"]), category=CATEGORIES),
        expand(expand("datavis/{condition}-v-{control}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{condition}-v-{control}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tfiib-peaks-matched-to-tss-peaks-distance{{distance}}-scatter.svg", zip, condition=CONDITIONS, control=CONTROLS), tss_norm=TSS_NORMS, tfiib_norm=TFIIB_NORMS, distance=int(config["search_distance"]))

rule match_tfiib_to_tss:
    input:
        tss_results = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{category}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-{category}-all.tsv"),
        tfiib_results = nexus_pipe("diff_binding/peaks/{condition}-v-{control}/{tfiib_norm}/{condition}-v-{control}_tfiib-chipnexus-{tfiib_norm}-peaks-diffbind-results-all.tsv"),
        genome_fasta = config["genome"]["fasta"]
    output:
        all = "matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{category}/{condition}-v-{control}_{category}-tss-seq-{tss_norm}-tfiib-{tfiib_norm}-matched-peaks-distance{distance}.tsv",
        matched= "matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{category}/{condition}-v-{control}_{category}-tss-seq-{tss_norm}-tfiib-{tfiib_norm}-matched-peaks-distance{distance}-matched.tsv",
        unmatched = "matched_peaks/{condition}-v-{control}/tfiib_matched_to_tss/{category}/{condition}-v-{control}_{category}-tss-seq-{tss_norm}-tfiib-{tfiib_norm}-matched-peaks-distance{distance}-unmatched.tsv",
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
        cat <(echo -e "tss_chrom\ttss_start\ttss_end\ttss_name\ttss_score\ttss_strand\ttss_lfc\ttss_lfc_SE\ttss_FDR\ttss_expr_condition\ttss_expr_control\ttss_summit\t{params.feature_header}tfiib_chrom\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_score\ttfiib_strand\ttfiib_lfc\ttfiib_lfc_SE\ttfiib_FDR\ttfiib_abundance_condition\ttfiib_abundance_control") - | \
        tee {output.all} | \
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{for(i=13; i<=NF-9; i++) {{if($i=="tfiib_start") tfiib_start=i}}}} \
            (NR>1 && $tfiib_start>=0) {{print $0 > "{output.matched}"}} \
            (NR>1 && $tfiib_start<0)  {{print $0 > "{output.unmatched}"}}'
        """

rule match_tss_to_tfiib:
    input:
        tfiib_results = tfiib_pipe("diff_binding/peaks/{condition}-v-{control}/{tfiib_norm}/{category}/{condition}-v-{control}_tfiib-chipnexus-{tfiib_norm}-peaks-diffbind-results-{category}-all.tsv"),
        tss_results = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-all.tsv"),
        genome_fasta = config["genome"]["fasta"]
    output:
        all = "matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{category}/{condition}-v-{control}_{category}-tfiib-{tfiib_norm}-tss-seq-{tss_norm}-matched-peaks-distance{distance}.tsv",
        matched = "matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{category}/{condition}-v-{control}_{category}-tfiib-{tfiib_norm}-tss-seq-{tss_norm}-matched-peaks-distance{distance}-matched.tsv",
        unmatched = "matched_peaks/{condition}-v-{control}/tss_matched_to_tfiib/{category}/{condition}-v-{control}_{category}-tfiib-{tfiib_norm}-tss-seq-{tss_norm}-matched-peaks-distance{distance}-unmatched.tsv",
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
        cat <(echo -e "tfiib_chrom\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_score\ttfiib_strand\ttfiib_lfc\ttfiib_lfc_SE\ttfiib_FDR\ttfiib_abundance_condition\ttfiib_abundance_control\ttfiib_summit\t{params.feature_header}tss_chrom\ttss_start\ttss_end\ttss_name\ttss_score\ttss_strand\ttss_lfc\ttss_lfc_SE\ttss_FDR\ttss_expr_condition\ttss_expr_control") - | \
        tee {output.all} | \
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{for(i=13; i<=NF-9; i++) {{if($i=="tss_start") tss_start=i}}}} \
            (NR>1 && $tss_start>=0) {{print $0 > "{output.matched}"}} \
            (NR>1 && $tss_start<0)  {{print $0 > "{output.unmatched}"}}'
        """

rule create_tfiib_windows:
    input:
        tss_tsv = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-all.tsv"),
        tss_narrowpeak = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-all.narrowpeak"),
        genome_fasta = config["genome"]["fasta"]
    output:
        "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-upstream-window-{window_size}.tsv"
    shell: """
        cut -f5,9,10,12 --complement {input.tss_tsv} | \
        tail -n +2 | \
        join --nocheck-order -t $'\t' -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2 - <(cut -f4,10 {input.tss_narrowpeak}) | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2+$11, $2+$11+1, $4, 0, $5, $0}}' | \
        bedtools slop -l {wildcards.window_size} -r 0 -s -i stdin -g <(faidx -i chromsizes {input.genome_fasta}) | \
        sort -k1,1 -k2,2n > {output}
        """

rule map_counts_to_windows:
    input:
        bed = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-seq-{tss_norm}-upstream-window-{window_size}.tsv",
        bg = lambda wc: nexus_pipe(f"coverage/counts/{wc.sample}_tfiib-chipnexus-counts-midpoints.bedgraph")
    output:
        temp("window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{sample}_experimental-tfiib-chipnexus-counts-tss-window-{window_size}.tsv")
    shell: """
        bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}
        """

rule combine_window_counts:
    input:
        lambda wc: ["window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/".format(**wc) + sample + "_experimental-tfiib-chipnexus-counts-tss-window-{window_size}.tsv".format(**wc) for sample in get_samples("libsizenorm", [wc.control, wc.condition])]
    output:
        "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_allsamples-experimental-tfiib-chipnexus-counts-tss-window-{window_size}.tsv"
    params:
        n = lambda wc: 18*len(get_samples("libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("libsizenorm", [wc.control, wc.condition]))
    shell: """
        paste {input} | \
         cut -f$(paste -d, <(echo "1-4,6,8,9,12-17") <(seq -s, 18 18 {params.n})) | \
         cat <(echo -e "chrom\twindow_start\twindow_end\ttss_name\tstrand\ttss_start\ttss_end\ttss_lfc\ttss_lfc_SE\ttss_fdr\ttss_expr_condition\ttss_expr_control\ttss_summit\t{params.names}") - > {output}
        """

rule differential_binding:
    input:
        exp_counts = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_allsamples-experimental-tfiib-chipnexus-counts-tss-window-{window_size}.tsv",
        spike_counts = lambda wc: nexus_pipe("diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{factor}-chipnexus-counts-peaks.tsv.gz") if wc.norm=="spikenorm" else []
    output:
        results_all = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-{tss_norm}-window-{window_size}-tfiib-chipnexus-{norm}-diffbind-results.tsv",
        counts_norm = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-{tss_norm}-window-{window_size}-tfiib-chipnexus-{norm}-counts-sizefactornorm.tsv",
        counts_rlog = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-{tss_norm}-window-{window_size}-tfiib-chipnexus-{norm}-counts-rlogtransform.tsv"
    params:
        samples = lambda wc: get_samples(wc.norm, [wc.control, wc.condition]),
        groups = lambda wc: [PASSING[x]["group"] for x in get_samples(wc.norm, [wc.control, wc.condition])],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold_change_threshold"]),
    conda:
        "envs/diff_bind_minimal.yaml"
    script:
        "scripts/window_differential_binding.R"

rule reclassify_windows:
    input:
        window_results = "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-{tss_norm}-window-{window_size}-tfiib-chipnexus-{norm}-diffbind-results.tsv",
        classified_peaks = tss_pipe("diff_exp/peaks/{condition}-v-{control}/{tss_norm}/{category}/{condition}-v-{control}_tss-seq-{tss_norm}-peaks-diffexp-results-{category}-all.tsv"),
    output:
        "window_diff_binding/{condition}-v-{control}/tss-{tss_norm}/{condition}-v-{control}_tss-{tss_norm}-window-{window_size}-tfiib-chipnexus-{norm}-diffbind-results-{category}.tsv",
    shell: """
        tail -n +2 {input.window_results} | \
        sort -k4,4 | \
        join -t $'\t' -1 4 -2 1 - <(tail -n +2 {input.classified_peaks} | \
                                    cut -f1-3,5-15 --complement | \
                                    sort -k1,1) | \
        sort -k8,8nr | \
        cat <(paste <(head -n 1 {input.window_results} | \
                      cut -f4) \
                    <(head -n 1 {input.window_results} | \
                      cut -f4 --complement) \
                    <(head -n 1 {input.classified_peaks} | \
                      cut -f1-15 --complement)) - > {output}
        """

rule datavis:
    input:
        matched_tfiib_paths = expand("matched_peaks/{{condition}}-v-{{control}}/tfiib_matched_to_tss/{category}/{{condition}}-v-{{control}}_{category}-tss-seq-{{tss_norm}}-tfiib-{{tfiib_norm}}-matched-peaks-distance{distance}.tsv", category=CATEGORIES, distance=config["search_distance"]),
        matched_tss_paths = expand("matched_peaks/{{condition}}-v-{{control}}/tss_matched_to_tfiib/{category}/{{condition}}-v-{{control}}_{category}-tfiib-{{tfiib_norm}}-tss-seq-{{tss_norm}}-matched-peaks-distance{distance}.tsv", category=["genic", "intragenic", "intergenic"], distance=config["search_distance"]),
        window_paths = expand("window_diff_binding/{{condition}}-v-{{control}}/tss-{{tss_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-window-{window_size}-tfiib-chipnexus-{{tfiib_norm}}-diffbind-results-{category}.tsv", category=CATEGORIES, window_size=config["window_size"])
    output:
        tfiib_matched_to_tss = "datavis/{{condition}}-v-{{control}}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tfiib-peaks-matched-to-tss-peaks-distance{distance}-scatter.svg".format(distance=config["search_distance"]),
        tss_matched_to_tfiib = "datavis/{{condition}}-v-{{control}}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tss-peaks-matched-to-tfiib-peaks-distance{distance}-scatter.svg".format(distance=config["search_distance"]),
        tss_vs_tfiib_windows = "datavis/{{condition}}-v-{{control}}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tss-peaks-vs-tfiib-window{window_size}-scatter.svg".format(window_size=config["window_size"]),
        matched_tfiib_mosaic = "datavis/{{condition}}-v-{{control}}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tfiib-peaks-matched-to-tss-peaks-distance{distance}-mosaic.svg".format(distance=config["search_distance"]),
        matched_tss_mosaic = "datavis/{{condition}}-v-{{control}}/tss-{{tss_norm}}-vs-tfiib-{{tfiib_norm}}/{{condition}}-v-{{control}}_tss-{{tss_norm}}-v-tfiib-{{tfiib_norm}}-tss-peaks-matched-to-tfiib-peaks-distance{distance}-mosaic.svg".format(distance=config["search_distance"]),
    params:
        fdr_cutoff_tss = config["fdr_cutoff_tss"],
        fdr_cutoff_tfiib = config["fdr_cutoff_tfiib"],
        window_size = config["window_size"],
        search_distance = config["search_distance"]
    conda:
        "envs/diff_bind_minimal.yaml"
    script:
        "scripts/tss_vs_tfiib.R"

