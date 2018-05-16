#!/usr/bin/env python
from math import log2

configfile: "config.yaml"

CONTROLS = [v for k,v in config["comparisons"].items()]
CONDITIONS = [k for k,v in config["comparisons"].items()]
GROUPS = CONDITIONS + CONTROLS
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]
COMPARISONS = [k + "-v-" + v for k,v in config["comparisons"].items()]

localrules: all,
    make_peak_tables, cat_peaks, plot_peaks

rule all:
    input:
        #mean coverage for each assay
        expand("coverage/{group}_{assay}.bedgraph", group=GROUPS, assay=["tss-seq", "tfiib-nexus"]),
        #peaks called independently for each assay, with coverage mapped to them
        expand("peaks_mapped/{group}_tss-seq-peaks-mapped-{category}.tsv", group=GROUPS, category=CATEGORIES+["all"]),
        expand("peaks_mapped/{group}_tfiib-nexus-peaks-mapped-{category}.tsv", group=GROUPS, category=["all", "genic", "intragenic", "intergenic"]),
        #stats on peaks called independently for each assay
        expand("peak_distances/{group}-tfiib-relative-to-tss-allcategories.tsv", group=GROUPS),
        expand("peak_distances/{group}-tss-relative-to-tfiib-allcategories.tsv", group=GROUPS),
        expand("peak_distances/{group}-peak-distances.svg", group=GROUPS),
        #distributions of coverage in bins
        expand("bins_mapped/{group}-tfiib-nexus-bins.svg", group=GROUPS),
        #regions above certain percentile of background signal
        expand("signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.bed", group=GROUPS, binsize=config["bin_size"], step_size=config["step_size"], signal_cutoff=config["signal_cutoff"]),
        expand("signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}_sensitivity-specificity.tsv", group=CONTROLS, binsize=config["bin_size"], step_size=config["step_size"], signal_cutoff=config["roc_cutoffs"]),
        expand("signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_sensitivity-specificity-all.tsv", group=CONTROLS, binsize=config["bin_size"], step_size=config["step_size"]),
        expand("signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-ROC.svg", group=CONTROLS, binsize=config["bin_size"], step_size=config["step_size"]),
        expand("signal_region_distances/{group}-peak-distances-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg", group=GROUPS, binsize=config["bin_size"], stepsize=config["step_size"], signal_cutoff=config["signal_cutoff"]),
        #differential binding of window upstream of TSS
        expand(expand("diff_binding/{condition}-v-{control}-tfiib-window-results-{{category}}.tsv", zip, condition=CONDITIONS, control=CONTROLS), category=CATEGORIES),
        expand("diff_binding/{condition}-v-{control}-tss-v-tfiib.svg", zip, condition=CONDITIONS, control=CONTROLS),

rule get_mean_coverage:
    input:
        lambda wc: config["groups"][wc.group]["tss_coverage"] if wc.assay=="tss-seq" else config["groups"][wc.group]["tfiib_coverage"]
    output:
        "coverage/{group}_{assay}.bedgraph"
    shell: """
        bedtools unionbedg -i {input} | awk 'BEGIN{{FS=OFS="\t"}}{{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1, $2, $3, sum/(NF-3)}}' > {output}
        """

# rule map_tfiib_coverage_to_intragenic_tfiib_peaks:
#     input:
#         peaks =  lambda wc: config["groups"][wc.condition]["tfiib_peaks"]["intragenic"],
#         condition_coverage = "coverage/{condition}_tfiib-nexus.bedgraph",
#         control_coverage  = "coverage/{control}_tfiib-nexus.bedgraph",
#         genic_anno = config["genome"]["genic_regions"],
#         orfs_anno = config["genome"]["orfs"],
#     output:
#         peaks = "bins_mapped/{condition}-v-{control}-intragenic-bins-inside-tfiib-peaks.tsv",
#         complement = "bins_mapped/{condition}-v-{control}-intragenic-bins-outside-tfiib-peaks.tsv",
#     params:
#         bin_size = config["bin_size"],
#         step_size = config["step_size"]
#     shell: """
#         bedtools makewindows -w {params.bin_size} -s {params.step_size} -b {input.peaks} | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.condition_coverage} -c 4 -o sum | bedtools map -a stdin -b {input.control_coverage} -c 4 -o sum > {output.peaks}
#         bedtools subtract -a {input.orfs_anno} -b {input.genic_anno} -s | bedtools subtract -a stdin -b {input.peaks} | bedtools makewindows -w {params.bin_size} -s {params.step_size} -b stdin | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.condition_coverage} -c 4 -o sum | bedtools map -a stdin -b {input.control_coverage} -c 4 -o sum > {output.complement}
#         """

# rule map_tfiib_coverage_to_upstream_tss_peaks:
#     input:
#         peaks =  lambda wc: config["groups"][wc.condition]["tss_peaks"]["intragenic"],
#         condition_coverage = "coverage/{condition}_tfiib-nexus.bedgraph",
#         control_coverage  = "coverage/{control}_tfiib-nexus.bedgraph",
#         orfs_anno = config["genome"]["orfs"],
#         chrom_sizes = config["genome"]["chrsizes"]
#     output:
#         upstream = "bins_mapped/{condition}-v-{control}-bins-upstream-intragenic-tss-peaks.tsv",
#         complement = "bins_mapped/{condition}-v-{control}-intragenic-bins-not-upstream-intragenic-tss-peaks.tsv",
#     params:
#         bin_size = config["bin_size"],
#         step_size = config["step_size"],
#         search_dist = config["search-distance"]
#     shell: """
#         awk 'BEGIN{{FS=OFS="\t"}}{{summit=$2+$10; print $1, summit, summit+1, $4, $5, $6}}' {input.peaks} | bedtools slop -i stdin -l {params.search_dist} -r 0 -s -i stdin -g {input.chrom_sizes} | bedtools makewindows -w {params.bin_size} -s {params.step_size} -b stdin | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.condition_coverage} -c 4 -o sum | bedtools map -a stdin -b {input.control_coverage} -c 4 -o sum > {output.upstream}
#         awk 'BEGIN{{FS=OFS="\t"}}{{summit=$2+$10; print $1, summit, summit+1, $4, $5, $6}}' {input.peaks} | bedtools slop -i stdin -l {params.search_dist} -r 0 -s -i stdin -g {input.chrom_sizes} | bedtools subtract -a {input.orfs_anno} -b stdin -s | bedtools makewindows -w {params.bin_size} -s {params.step_size} -b stdin | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.condition_coverage} -c 4 -o sum | bedtools map -a stdin -b {input.control_coverage} -c 4 -o sum > {output.complement}
#         """

#in this section, look at metrics for TSS and TFIIB peaks called by conventional (narrow peak)
# peak-callers: distance to closest peak of the other class, etc.
rule map_coverage_to_tss_peaks:
    input:
        peaks = lambda wc: config["groups"][wc.group]["tss_peaks"][wc.category],
        coverage = "coverage/{group}_tss-seq.bedgraph"
    output:
        "peaks_mapped/{group}_tss-seq-peaks-mapped-{category}.tsv"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{$1=$1"-plus"}} $6=="-"{{$1=$1"-minus"}}{{print $0}}' {input.peaks} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.coverage} -c 4 -o sum | sed 's/-minus//g;s/-plus//g' | sort -k1,1 -k2,2n > {output}
        """

rule map_coverage_to_tfiib_peaks:
    input:
        peaks = lambda wc: config["groups"][wc.group]["tfiib_peaks"][wc.category],
        coverage = "coverage/{group}_tfiib-nexus.bedgraph"
    output:
        "peaks_mapped/{group}_tfiib-nexus-peaks-mapped-{category}.tsv"
    shell: """
        cut -f1-10 {input.peaks} | uniq | bedtools map -a stdin -b {input.coverage} -c 4 -o sum > {output}
        """

rule get_closest_upstream_tfiib_peak:
    input:
        tss_peaks = "peaks_mapped/{group}_tss-seq-peaks-mapped-{category}.tsv",
        tfiib_peaks = "peaks_mapped/{group}_tfiib-nexus-peaks-mapped-all.tsv"
    output:
        temp("peak_distances/.{group}-tfiib-relative-to-tss-{category}.tsv")
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tss_peaks} | sort -k1,1 -k2,2n | bedtools closest -D a -id -t first -k 2 -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tfiib_peaks} | sort -k1,1 -k2,2n) | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.category}"}}' > {output}
        """

rule cat_closest_upstream_tfiib_peak:
    input:
        expand("peak_distances/.{{group}}-tfiib-relative-to-tss-{category}.tsv", category=["all"]+CATEGORIES)
    output:
        "peak_distances/{group}-tfiib-relative-to-tss-allcategories.tsv"
    shell: """
        cat {input} > {output}
        """

rule get_closest_tss_peaks:
    input:
        tfiib_peaks = "peaks_mapped/{group}_tfiib-nexus-peaks-mapped-{category}.tsv",
        tss_peaks = "peaks_mapped/{group}_tss-seq-peaks-mapped-all.tsv"
    output:
        temp("peak_distances/.{group}-tss-relative-to-tfiib-{category}.tsv")
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tfiib_peaks} | sort -k1,1 -k2,2n | bedtools closest -D ref -t first -k 4 -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tss_peaks} | sort -k1,1 -k2,2n) | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.category}"}}' > {output}
        """

rule cat_closest_tss_peaks:
    input:
        expand("peak_distances/.{{group}}-tss-relative-to-tfiib-{category}.tsv", category=["all", "genic", "intragenic", "intergenic"])
    output:
        "peak_distances/{group}-tss-relative-to-tfiib-allcategories.tsv"
    shell: """
        cat {input} > {output}
        """

rule plot_peak_information:
    input:
        tfiib_rel_tss = "peak_distances/{group}-tfiib-relative-to-tss-allcategories.tsv",
        tss_rel_tfiib = "peak_distances/{group}-tss-relative-to-tfiib-allcategories.tsv"
    params:
        max_dist = config["search-distance"]
    output:
        peak_distances = "peak_distances/{group}-peak-distances.svg",
        tss_signal_violin = "peak_distances/{group}-TSS-signal-violin.svg",
        tss_signal_density = "peak_distances/{group}-TSS-signal-density.svg",
        tfiib_signal_violin = "peak_distances/{group}-tfiib-signal-violin.svg",
        tfiib_signal_density = "peak_distances/{group}-tfiib-signal-density.svg",
        tfiib_rel_tss_mosaic = "peak_distances/{group}-tfiib-rel-tss-mosaic.svg",
        tss_rel_tfiib_mosaic = "peak_distances/{group}-tss-rel-tfiib-mosaic.svg"
    script:
        "scripts/plot_peak_information.R"

#in this section, look at differentially expressed peaks called by conventional
# (narrow) peak callers

#in this section, look at TFIIB bins: those mapping within conventionally called peaks,
#compared to background (not in peaks)
rule map_coverage_to_tfiib_peak_bins:
    input:
        peaks = lambda wc: config["groups"][wc.group]["tfiib_peaks"][wc.category],
        coverage = "coverage/{group}_tfiib-nexus.bedgraph"
    output:
        temp("bins_mapped/{group}_tfiib-nexus-bins-in-{category}-tfiib-peaks.tsv")
    params:
        bin_size = config["bin_size"],
        step_size = config["step_size"],
    shell: """
        cut -f1-6 {input.peaks} | bedtools makewindows -w {params.bin_size} -s {params.step_size} -b stdin | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.coverage} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.category}"}}' > {output}
        """

rule map_coverage_to_background_bins:
    input:
        peaks = lambda wc: config["groups"][wc.group]["tfiib_peaks"]["all"],
        chrsizes = config["genome"]["chrsizes"],
        coverage = "coverage/{group}_tfiib-nexus.bedgraph"
    output:
        temp("bins_mapped/{group}_tfiib-nexus-background-bins.tsv")
    params:
        bin_size = config["bin_size"],
        step_size = config["step_size"],
    shell: """
        bedtools complement -i {input.peaks} -g <(sort -k1,1 {input.chrsizes}) | bedtools makewindows -w {params.bin_size} -s {params.step_size} -b stdin | bedtools map -a stdin -b {input.coverage} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "background"}}' > {output}
        """

rule cat_bin_coverage:
    input:
        "bins_mapped/{group}_tfiib-nexus-background-bins.tsv",
        expand("bins_mapped/{{group}}_tfiib-nexus-bins-in-{category}-tfiib-peaks.tsv", category=["all", "genic", "intragenic", "intergenic"])
    output:
        "bins_mapped/{group}-tfiib-nexus-bins-all.tsv"
    shell: """
        cat {input} > {output}
        """

rule plot_bin_coverage:
    input:
        "bins_mapped/{group}-tfiib-nexus-bins-all.tsv"
    output:
        "bins_mapped/{group}-tfiib-nexus-bins.svg"
    params:
        bin_size = config["bin_size"],
        quantile = config["background_quantile"]
    script:
        "scripts/tfiib_bin_signal_distribution.R"

rule get_signal_regions:
    input:
        coverage = "coverage/{group}_tfiib-nexus.bedgraph",
        chrsizes = config["genome"]["chrsizes"],
    output:
        "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.bed"
    params:
        bin_size = config["bin_size"],
        step_size = config["step_size"],
        # signal_cutoff = config["signal_cutoff"]
    shell: """
        bedtools makewindows -w {params.bin_size} -s {params.step_size} -g <(sort -k1,1 {input.chrsizes}) | bedtools map -a stdin -b {input.coverage} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}} ($4/($3-$2))>{wildcards.signal_cutoff}' | bedtools merge -d {params.bin_size} -c 4 -o mean > {output}
        """

rule calculate_signal_region_sensitivity_specificity:
    input:
        tfiib_regions = "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.bed",
        tfiib_peaks = lambda wc: config["groups"][wc.group]["tfiib_peaks"]["all"],
        chrsizes = config["genome"]["chrsizes"],
    output:
        "signal_regions/{{group}}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}_sensitivity-specificity.tsv"
    params:
        bin_size = config["bin_size"],
        step_size = config["step_size"],
    shell: """
        bedtools makewindows -w {params.bin_size} -s {params.step_size} -g <(sort -k1,1 {input.chrsizes}) | bedtools intersect -a stdin -b <(cut -f1-3 {input.tfiib_regions}) -loj | bedtools intersect -a stdin -b <(cut -f1-3 {input.tfiib_peaks} | sort -k1,1 -k2,2n) -loj | awk 'BEGIN{{FS=OFS="\t"; TN=FN=FP=TP=0}} {{ if($5==-1) {{$8==-1 ? TN++ : FN++}} else {{$8==-1? FP++ : TP++}} }} END {{print "{wildcards.signal_cutoff}", TN, FN, FP, TP}}' > {output}
        """

rule cat_sens_spec:
    input:
        expand("signal_regions/{{group}}_tfiib-signal-regions-binsize_{{binsize}}-stepsize_{{step_size}}-cutoff_{signal_cutoff}_sensitivity-specificity.tsv", signal_cutoff=config['roc_cutoffs'])
    output:
        "signal_regions/{{group}}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_sensitivity-specificity-all.tsv"
    shell: """
        cat <(echo -e "cutoff\tTN\tFN\tFP\tTP")  {input} > {output}
        """

rule plot_roc:
    input:
        "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_sensitivity-specificity-all.tsv"
    output:
        roc = "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-ROC.svg",
        roc_zoom = "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-ROCzoom.svg",
    script:
        "scripts/signal_region_roc.R"

rule get_closest_upstream_tfiib_signal_region:
    input:
        tss_peaks = "peaks_mapped/{group}_tss-seq-peaks-mapped-{category}.tsv",
        tfiib_regions = "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.bed",
        coverage = "coverage/{group}_tfiib-nexus.bedgraph",
    output:
        temp("signal_region_distances/.{group}-tfiib-relative-to-tss-{category}-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.tsv")
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tss_peaks} | sort -k1,1 -k2,2n | bedtools closest -D a -id -t first -k 2 -a stdin -b <(bedtools map -a {input.tfiib_regions} -b {input.coverage} -c 4 -o sum) | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.category}"}}' > {output}
        """

rule get_closest_tss_peaks_to_tfiib_signal_region:
    input:
        tfiib_regions = "signal_regions/{group}_tfiib-signal-regions-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.bed",
        tss_peaks = "peaks_mapped/{group}_tss-seq-peaks-mapped-all.tsv",
        coverage = "coverage/{group}_tfiib-nexus.bedgraph",
    output:
        "signal_region_distances/{group}-tss-relative-to-tfiib-signal-region-binsize_{binsize}-stepsize_{step_size}-cutoff_{signal_cutoff}.tsv"
    shell: """
        bedtools map -a {input.tfiib_regions} -b {input.coverage} -c 4 -o sum | bedtools closest -D ref -t first -k 4 -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tss_peaks} | sort -k1,1 -k2,2n)  > {output}
        """

rule cat_closest_upstream_tfiib_region:
    input:
        expand("signal_region_distances/.{{group}}-tfiib-relative-to-tss-{category}-binsize_{{binsize}}-stepsize_{{stepsize}}-cutoff_{{signal_cutoff}}.tsv", category=["all"]+CATEGORIES)
    output:
        "signal_region_distances/{group}-tfiib-signal-region-relative-to-tss-allcategories-binsize_{binsize}-stepsize_{stepsize}-cutoff_{signal_cutoff}.tsv"
    shell: """
        cat {input} > {output}
        """

rule plot_region_information:
    input:
        tfiib_rel_tss = "signal_region_distances/{group}-tfiib-signal-region-relative-to-tss-allcategories-binsize_{binsize}-stepsize_{stepsize}-cutoff_{signal_cutoff}.tsv",
        tss_rel_tfiib = "signal_region_distances/{group}-tss-relative-to-tfiib-signal-region-binsize_{binsize}-stepsize_{stepsize}-cutoff_{signal_cutoff}.tsv"
    params:
        max_dist = config["search-distance"]
    output:
        signal_region_distances = "signal_region_distances/{group}-peak-distances-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg",
        tss_signal_violin = "signal_region_distances/{group}-TSS-signal-violin-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg",
        tss_signal_density = "signal_region_distances/{group}-TSS-signal-density-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg",
        tfiib_signal_density = "signal_region_distances/{group}-tfiib-signal-density-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg",
        tfiib_rel_tss_mosaic = "signal_region_distances/{group}-tfiib-rel-tss-mosaic-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg",
        tss_rel_tfiib_mosaic = "signal_region_distances/{group}-tss-rel-tfiib-mosaic-tfiib-regions-binsize_{binsize}-stepsize-{stepsize}-cutoff-{signal_cutoff}.svg"
    script:
        "scripts/plot_region_information.R"

### do a differential binding analysis using the TFIIB signal upstream of TSSs
rule map_tfiib_counts_upstream_of_tss:
    input:
        tss_results = lambda wc: config["diff_exp_and_binding"][wc.condition][wc.control],
        tfiib_counts = lambda wc: config["tfiib_counts"][wc.sample]["path"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        "diff_binding/{condition}-v-{control}_{sample}-tfiib-counts.tsv"
    params:
        window_size = config["search-distance"]
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.tss_results} | bedtools slop -i stdin -g {input.chrsizes} -l {params.window_size} -r 0 -s | cut -f1-6 | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.tfiib_counts} -c 4 -o sum | sort -k4,4 | paste <(sort -k4,4 {input.tss_results}) - > {output}
        """

def get_samples(groups: list):
    return [k for k,v in config["tfiib_counts"].items() if v["group"] in groups]

rule window_differential_binding:
    input:
        controls = lambda wc: ["diff_binding/" + wc.condition + "-v-" + wc.control + "_" + sample + "-tfiib-counts.tsv" for sample in get_samples([wc.control])],
        conditions = lambda wc: ["diff_binding/" + wc.condition + "-v-" + wc.control + "_" + sample + "-tfiib-counts.tsv" for sample in get_samples([wc.condition])],
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv"
    params:
        alpha=0.1,
        lfc=log2(1.5)
    script:
        "scripts/call_de_peaks.R"

rule classify_diffbind_genic:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        genic_anno = config["genome"]["genic_regions"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-genic.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.genic_anno} -u -s | cat <(head -n 1 {input.results}) - > {output}
        """

rule classify_diffbind_intragenic:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        genic_anno = config["genome"]["genic_regions"],
        orf_anno = config["genome"]["orfs"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-intragenic.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orf_anno} -u -s | cat <(head -n 1 {input.results}) - > {output}
        """

rule classify_diffbind_antisense:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        transcript_anno = config["genome"]["transcripts"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-antisense.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.transcript_anno} -u -S | cat <(head -n 1 {input.results}) - > {output}
        """

rule classify_diffbind_convergent:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        genic_anno = config["genome"]["genic_regions"],
        conv_anno = config["genome"]["conv_regions"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-convergent.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -u -s | cat <(head -n 1 {input.results}) - > {output}
        """

rule classify_diffbind_divergent:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        genic_anno = config["genome"]["genic_regions"],
        div_anno = config["genome"]["div_regions"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-divergent.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -u -s | cat <(head -n 1 {input.results}) - > {output}
        """

rule classify_diffbind_intergenic:
    input:
        results = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        intergenic_anno = config["genome"]["intergenic_regions"],
        transcript_anno = config["genome"]["transcripts"],
        genic_anno = config["genome"]["genic_regions"],
        orf_anno = config["genome"]["orfs"]
    output:
        "diff_binding/{condition}-v-{control}-tfiib-window-results-intergenic.tsv",
    shell: """
        bedtools intersect -a {input.results} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | cat <(head -n 1 {input.results}) - > {output}
        """

rule plot_diffbind:
    input:
        all = "diff_binding/{condition}-v-{control}-tfiib-window-results.tsv",
        genic = "diff_binding/{condition}-v-{control}-tfiib-window-results-genic.tsv",
        intragenic = "diff_binding/{condition}-v-{control}-tfiib-window-results-intragenic.tsv",
        antisense = "diff_binding/{condition}-v-{control}-tfiib-window-results-antisense.tsv",
        convergent = "diff_binding/{condition}-v-{control}-tfiib-window-results-convergent.tsv",
        divergent = "diff_binding/{condition}-v-{control}-tfiib-window-results-divergent.tsv",
        intergenic = "diff_binding/{condition}-v-{control}-tfiib-window-results-intergenic.tsv",
    output:
        "diff_binding/{condition}-v-{control}-tss-v-tfiib.svg"
    script:
        "scripts/plot_tss_v_tfiib.R"

# rule make_peak_tables:
#     input:
#         tss_peaks = lambda wildcards: config["tss"]["peaks-path"] + wildcards.condition + "-exp-idrpeaks.narrowPeak" if wildcards.category=="all" else config["tss"]["peaks-path"] + wildcards.category + "/" + wildcards.condition + "-exp-idrpeaks-" + wildcards.category + ".tsv",
#         tss_coverage = lambda wildcards: config["tss"]["coverage"] + config["tss"]["norm"] + "/" + wildcards.condition + "-1-tss-" + config["tss"]["norm"] + "-SENSE.bedgraph",
#         tfiib_peaks = lambda wildcards: config["nexus"]["peaks-path"] + wildcards.condition + "-" + config["genome"]["prefix"] + "_peaks.narrowPeak",
#         tfiib_coverage = lambda wildcards: config["nexus"]["coverage"] + config["nexus"]["norm"] + "/" + wildcards.condition + "-1-tfiib-chipnexus-" + config["nexus"]["norm"] + "-combined.bedgraph",
#         chrsizes = config["genome"]["chrsizes"]
#     params: upstr = config["tss-upstr-dist"]
#     output:
#         "tables/{condition}-tfiib-upstr-tss-{category}.tsv"
#     shell: """
#         cut -f1-10 {input.tss_peaks} | awk 'BEGIN{{FS=OFS="\t"}}$6=="+"{{$1=$1"-plus"; print $0}}$6=="-"{{$1=$1"-minus"; print $0}}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.tss_coverage} -c 4 -o sum | sed "s/-minus//g;s/-plus//g" | sort -k1,1 -k2,2n | paste <(sort -k1,1 -k2,2n {input.tss_peaks} | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.tfiib_coverage} -c 4 -o sum | bedtools intersect -a stdin -b {input.tfiib_peaks} -wao | cut --complement -f1-3,7,8,10,11,16,19,20,21,23 | awk -v condition={wildcards.condition} -v category={wildcards.category} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, category}}' > {output}
#         """

# rule cat_peaks:
#     input:
#         expand("tables/{condition}-tfiib-upstr-tss-{category}.tsv", condition=CONDITIONS, category=CATEGORIES+["all"])
#     output:
#         "tables/all-tfiib-upstr-tss.tsv"
#     shell: """
#         cat <(echo -e "chrom\ttss_start\ttss_end\tstrand\ttss_idr\ttss_summit\ttss_expression\ttfiib_levels\ttfiib_start\ttfiib_end\ttfiib_enrichment\ttfiib_qval\ttfiib_summit\toverlap\tcondition\tcategory") {input} > {output}
#         """

# rule plot_peaks:
#     input:
#         "tables/all-tfiib-upstr-tss.tsv"
#     params:
#         dist = config["tss-upstr-dist"]
#     output:
#         mosaic = "figures/peaks-tss-v-tfiib-overlap-mosaic.svg",
#         tss_expr = "figures/peaks-tss-v-tfiib-tss-expression.svg",
#         tss_size = "figures/peaks-tss-v-tfiib-tss-sizes.svg",
#         tss_idr = "figures/peaks-tss-v-tfiib-tss-idr.svg",
#         tfiib_lvl = "figures/peaks-tss-v-tfiib-tfiib-levels.svg",
#         contour = "figures/peaks-tss-v-tfiib-expression-contourplot.svg",
#         scatter = "figures/peaks-tss-v-tfiib-expression-scatter.svg"
#     script: "scripts/tss-v-tfiib.R"

# rule make_diffexp_tables_all:
#     input:
#         tss = lambda wildcards: config["tss"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-results-" + config["tss"]["norm"] + "-all.tsv",
#         nexus = lambda wildcards: config["nexus"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-tfiib-chipnexus-results-" + config["nexus"]["norm"] + "-all.tsv",
#         chrsizes = config["genome"]["chrsizes"]
#     params: upstr = config["tss-upstr-dist"]
#     output:
#         "diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-all.tsv"
#     shell: """
#         tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $11, $12, $13, "na", "na", "na"}}' | sort -k1,1 -k2,2n | paste <(tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3}}' | sort -k1,1 -k2,2n | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | bedtools intersect -a stdin -b <(tail -n +2 {input.nexus} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n) -wao | cut --complement -f1-3,8,18,22,23,26-28 | awk -v condition={wildcards.condition} -v control={wildcards.control} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, control, "all"}}' > {output}
#         """

# rule make_diffexp_tables_category:
#     input:
#         tss = lambda wildcards: config["tss"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.category + "/" + wildcards.condition + "-v-" + wildcards.control + "-results-" + config["nexus"]["norm"] + "-all-" + wildcards.category + ".tsv",
#         nexus = lambda wildcards: config["nexus"]["diffexp-path"] + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-tfiib-chipnexus-results-" + config["nexus"]["norm"] + "-all.tsv",
#         chrsizes = config["genome"]["chrsizes"]
#     params:
#         upstr = config["tss-upstr-dist"]
#     output:
#         "diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-{category}.tsv"
#     shell: """
#         tail -n +2 {input.tss} | awk -v ttype={wildcards.category} 'BEGIN{{FS=OFS="\t"}} ttype=="intergenic"{{$13="na"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n | paste <(tail -n +2 {input.tss} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3}}' | sort -k1,1 -k2,2n | bedtools flank -i stdin -g {input.chrsizes} -l {params.upstr} -r 0 -s | cut -f1-3) - | bedtools intersect -a stdin -b <(tail -n +2 {input.nexus} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, 0, $3, $6, $7, $8, $9, $10, $11, $12, $13}}' | sort -k1,1 -k2,2n) -wao | cut --complement -f1-3,8,18,22,23,26-28 | awk -v category={wildcards.category} -v condition={wildcards.condition} -v control={wildcards.control} 'BEGIN{{FS=OFS="\t"}}{{print $0, condition, control, category}}' > {output}
#         """

# rule cat_diffexp:
#     input:
#         expand(expand("diffexp_tables/{condition}-v-{control}-tfiib-upstr-tss-{{category}}.tsv", zip, condition=config["conditions"], control=config["controls"]), category=CATEGORIES+["all"])
#     output:
#         "diffexp_tables/all-diffexp-tfiib-upstr-tss.tsv"
#     shell: """
#         cat <(echo -e "chrom\ttss_start\ttss_end\ttss_name\tstrand\ttss_meanExpr\ttss_lfc\ttss_logpadj\ttss_cond\ttss_ctrl\tfeat_start\tfeat_end\tfeat_name\ttfiib_start\ttfiib_end\ttfiib_name\ttfiib_meanExpr\ttfiib_lfc\ttfiib_logpadj\ttfiib_cond\ttfiib_ctrl\toverlap\tcondition\tcontrol\tcategory") {input} > {output}
#         """

# rule plot_diffexp:
#     input:
#         "diffexp_tables/all-diffexp-tfiib-upstr-tss.tsv"
#     output:
#         scatter = "figures/tss-v-tfiib-diffexp-scatter.svg"
#     script: "scripts/tss-v-tfiib-diffexp.R"

# rule get_bed:
#     input:
#         "tables/{condition}-tfiib-upstr-tss-{category}.tsv"
#     output:
#         yes = "tables/bed/{condition}_{category}-tss-with-tfiib-upstr.bed",
#         no = "tables/bed/{condition}_{category}-tss-without-tfiib-upstr.bed",
#     script: "scripts/getbed.R"

# rule compute_matrix:
#     input:
#         annotation = "tables/bed/{condition}-tss-{class}-tfiib-upstr.bed",
#         bw = lambda wildcards: config["nexus"]["coverage"] + config["nexus"]["norm"] + "/bw/" + wildcards.condition + "-1-tfiib-chipnexus-" + config["nexus"]["norm"] + "-qfrags.bw",
#     output:
#         dtfile = temp("figures/heatmaps/{condition}-{class}.mat"),
#         matrix = temp("figures/heatmaps/{condition}-{class}.tsv"),
#         matrix_gz = temp("figures/heatmaps/{condition}-{class}.tsv.gz"),
#     params:
#         refpoint = "TSS",
#         upstream = config["heatmaps"]["upstream"],
#         dnstream = config["heatmaps"]["downstream"],
#         binsize = config["heatmaps"]["binsize"],
#         sort = "keep",
#         binstat = config["heatmaps"]["binstat"]
#     threads : config["threads"]
#     log: "logs/compute_matrix-{condition}-{class}.log"
#     shell: """
#         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}
#         """

# rule melt_matrix:
#     input:
#         matrix = "figures/heatmaps/{condition}-{class}.tsv.gz"
#     output:
#         temp("figures/heatmaps/{condition}-{class}-melted.tsv.gz")
#     params:
#         binsize = config["heatmaps"]["binsize"],
#         upstream = config["heatmaps"]["upstream"],
#     script:
#         "scripts/melt_matrix.R"

# rule cat_matrices:
#     input:
#         expand("figures/heatmaps/{{condition}}-{zz}-melted.tsv.gz", zz=["with","without"])
#     output:
#         "figures/heatmaps/{condition}-both.tsv.gz"
#     shell: """
#         cat {input} > {output}
#         """

# rule plot_heatmaps:
#     input:
#         "figures/heatmaps/{condition}-both.tsv.gz"
#     output:
#         "figures/heatmaps/{condition}-heatmap.svg"
#     params:
#         cutoff=config["heatmaps"]["pct_cutoff"],
#         upstream=config["heatmaps"]["upstream"],
#         dnstream=config["heatmaps"]["downstream"],
#         cmap=config["heatmaps"]["colormap"]
#     script: "scripts/plotHeatmaps.R"
