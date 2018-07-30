#!/usr/bin/env python
from math import log2

configfile: "config.yaml"

CONTROLS = [v for k,v in config["comparisons"].items()]
CONDITIONS = [k for k,v in config["comparisons"].items()]
GROUPS = set(CONDITIONS + CONTROLS)
CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]
COMPARISONS = [f"{k}-v-{v}" for k,v in config["comparisons"].items()]

localrules: all,
    make_peak_tables, cat_peaks, plot_peaks

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

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

