
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('diff_binding/spt6-1004-37C-v-WT-37C_WT-37C-1-tfiib-counts.tsv', 'diff_binding/spt6-1004-37C-v-WT-37C_WT-37C-2-tfiib-counts.tsv', 'diff_binding/spt6-1004-37C-v-WT-37C_spt6-1004-37C-1-tfiib-counts.tsv', 'diff_binding/spt6-1004-37C-v-WT-37C_spt6-1004-37C-2-tfiib-counts.tsv', "controls" = c('diff_binding/spt6-1004-37C-v-WT-37C_WT-37C-1-tfiib-counts.tsv', 'diff_binding/spt6-1004-37C-v-WT-37C_WT-37C-2-tfiib-counts.tsv'), "conditions" = c('diff_binding/spt6-1004-37C-v-WT-37C_spt6-1004-37C-1-tfiib-counts.tsv', 'diff_binding/spt6-1004-37C-v-WT-37C_spt6-1004-37C-2-tfiib-counts.tsv')),
    output = list('diff_binding/spt6-1004-37C-v-WT-37C-tfiib-window-results.tsv'),
    params = list(0.1, 0.5849625007211562, "alpha" = 0.1, "lfc" = 0.5849625007211562),
    wildcards = list('spt6-1004-37C', 'WT-37C', "condition" = 'spt6-1004-37C', "control" = 'WT-37C'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("genome" = list("fasta" = '../genomefiles_cerevisiae/S_cerevisiae.R64-2-1.fa', "annotation_workflow" = '../build-annotations-cerevisiae/', "transcripts" = '../genomefiles_cerevisiae/annotations/Scer_polIItranscripts-adjustedTSS.bed', "orfs" = '../genomefiles_cerevisiae/annotations/Scer_nondubious_ORFs_and_blocked_reading_frames.bed', "prefix" = 'Scer_'), "search-distance" = 200, "bin_size" = 20, "step_size" = 10, "background_quantile" = 0.9, "signal_cutoff" = 0.081, "roc_cutoffs" = c(0, 0.01, 0.02, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.22, 0.3, 0.5, 1.0, 3.0), "groups" = list("WT-37C" = list("tss_peaks" = list("all" = '../tss-seq-publication/peakcalling/WT-37C/WT-37C_experimental-idrpeaks.narrowPeak', "genic" = '../tss-seq-publication/peakcalling/WT-37C/genic/WT-37C-experimental-idrpeaks-genic.narrowpeak', "intragenic" = '../tss-seq-publication/peakcalling/WT-37C/intragenic/WT-37C-experimental-idrpeaks-intragenic.narrowpeak', "antisense" = '../tss-seq-publication/peakcalling/WT-37C/antisense/WT-37C-experimental-idrpeaks-antisense.narrowpeak', "convergent" = '../tss-seq-publication/peakcalling/WT-37C/convergent/WT-37C-experimental-idrpeaks-convergent.narrowpeak', "divergent" = '../tss-seq-publication/peakcalling/WT-37C/divergent/WT-37C-experimental-idrpeaks-divergent.narrowpeak', "intergenic" = '../tss-seq-publication/peakcalling/WT-37C/intergenic/WT-37C-experimental-idrpeaks-intergenic.narrowpeak'), "tss_coverage" = c('../tss-seq-publication/coverage/spikenorm/WT-37C-1_tss-seq-spikenorm-SENSE.bedgraph', '../tss-seq-publication/coverage/spikenorm/WT-37C-2_tss-seq-spikenorm-SENSE.bedgraph'), "tfiib_peaks" = list("all" = '../chipnexus-tfiib-publication/peakcalling/macs/WT-37C/WT-37C_experimental-tfiib-chipnexus_peaks.narrowPeak', "genic" = '../chipnexus-tfiib-publication/peakcalling/macs/WT-37C/WT-37C_experimental-tfiib-chipnexus_peaks-genic.narrowpeak', "intragenic" = '../chipnexus-tfiib-publication/peakcalling/macs/WT-37C/WT-37C_experimental-tfiib-chipnexus_peaks-intragenic.narrowpeak', "intergenic" = '../chipnexus-tfiib-publication/peakcalling/macs/WT-37C/WT-37C_experimental-tfiib-chipnexus_peaks-intergenic.narrowpeak'), "tfiib_coverage" = c('../chipnexus-tfiib-publication/coverage/libsizenorm/WT-37C-1_tfiib-chipnexus-libsizenorm-midpoints.bedgraph', '../chipnexus-tfiib-publication/coverage/libsizenorm/WT-37C-2_tfiib-chipnexus-libsizenorm-midpoints.bedgraph')), "spt6-1004-37C" = list("tss_peaks" = list("all" = '../tss-seq-publication/peakcalling/spt6-1004-37C/spt6-1004-37C_experimental-idrpeaks.narrowPeak', "genic" = '../tss-seq-publication/peakcalling/spt6-1004-37C/genic/spt6-1004-37C-experimental-idrpeaks-genic.narrowpeak', "intragenic" = '../tss-seq-publication/peakcalling/spt6-1004-37C/intragenic/spt6-1004-37C-experimental-idrpeaks-intragenic.narrowpeak', "antisense" = '../tss-seq-publication/peakcalling/spt6-1004-37C/antisense/spt6-1004-37C-experimental-idrpeaks-antisense.narrowpeak', "convergent" = '../tss-seq-publication/peakcalling/spt6-1004-37C/convergent/spt6-1004-37C-experimental-idrpeaks-convergent.narrowpeak', "divergent" = '../tss-seq-publication/peakcalling/spt6-1004-37C/divergent/spt6-1004-37C-experimental-idrpeaks-divergent.narrowpeak', "intergenic" = '../tss-seq-publication/peakcalling/spt6-1004-37C/intergenic/spt6-1004-37C-experimental-idrpeaks-intergenic.narrowpeak'), "tss_coverage" = c('../tss-seq-publication/coverage/spikenorm/spt6-1004-37C-1_tss-seq-spikenorm-SENSE.bedgraph', '../tss-seq-publication/coverage/spikenorm/spt6-1004-37C-2_tss-seq-spikenorm-SENSE.bedgraph'), "tfiib_peaks" = list("all" = '../chipnexus-tfiib-publication/peakcalling/macs/spt6-1004-37C/spt6-1004-37C_experimental-tfiib-chipnexus_peaks.narrowPeak', "genic" = '../chipnexus-tfiib-publication/peakcalling/macs/spt6-1004-37C/spt6-1004-37C_experimental-tfiib-chipnexus_peaks-genic.narrowpeak', "intragenic" = '../chipnexus-tfiib-publication/peakcalling/macs/spt6-1004-37C/spt6-1004-37C_experimental-tfiib-chipnexus_peaks-intragenic.narrowpeak', "intergenic" = '../chipnexus-tfiib-publication/peakcalling/macs/spt6-1004-37C/spt6-1004-37C_experimental-tfiib-chipnexus_peaks-intergenic.narrowpeak'), "tfiib_coverage" = c('../chipnexus-tfiib-publication/coverage/libsizenorm/spt6-1004-37C-1_tfiib-chipnexus-libsizenorm-midpoints.bedgraph', '../chipnexus-tfiib-publication/coverage/libsizenorm/spt6-1004-37C-2_tfiib-chipnexus-libsizenorm-midpoints.bedgraph'))), "comparisons" = list("spt6-1004-37C" = 'WT-37C'), "diff_exp_and_binding" = list("spt6-1004-37C" = list("WT-37C" = '../tss-seq-publication/diff_exp/spt6-1004-37C-v-WT-37C/spikenorm/spt6-1004-37C-v-WT-37C_tss-seq-spikenorm-diffexp-results-all.narrowpeak')), "tfiib_counts" = list("spt6-1004-37C-1" = list("path" = '../chipnexus-tfiib-publication/coverage/counts/spt6-1004-37C-1_tfiib-chipnexus-counts-midpoints.bedgraph', "group" = 'spt6-1004-37C'), "spt6-1004-37C-2" = list("path" = '../chipnexus-tfiib-publication/coverage/counts/spt6-1004-37C-2_tfiib-chipnexus-counts-midpoints.bedgraph', "group" = 'spt6-1004-37C'), "WT-37C-1" = list("path" = '../chipnexus-tfiib-publication/coverage/counts/WT-37C-1_tfiib-chipnexus-counts-midpoints.bedgraph', "group" = 'WT-37C'), "WT-37C-2" = list("path" = '../chipnexus-tfiib-publication/coverage/counts/WT-37C-2_tfiib-chipnexus-counts-midpoints.bedgraph', "group" = 'WT-37C'))),
    rule = 'window_differential_binding'
)
######## Original script #########
library(tidyverse)
library(forcats)
library(DESeq2)

import_data = function(df, paths, group){
    for (i in seq_along(paths)){
        df = read_tsv(paths[i],
                      col_names = c("chrom", "tss_start", "tss_end", "tss_id", "tss_score", "tss_strand",
                                    "tss_lfc", "tss_logp", "tss_logq", "tss_summit",
                                    "a", "window_start", "window_end", "b", "c", "d", "tfiib_counts")) %>%
            select(-c(a,b,c,d)) %>%
            mutate(group=group,
                   replicate = i) %>%
            bind_rows(df, .)
    }
    return(df)
}

main = function(controls, conditions, control_id, condition_id,
                results_all, alpha, lfc){
    df = tibble()
    df = df %>% import_data(paths=controls, group=control_id) %>%
        import_data(paths=conditions, group=condition_id)

    countdata = df %>% select(tss_id, tfiib_counts, group, replicate) %>%
        unite(col="sample", group, replicate) %>%
        mutate(sample = fct_inorder(sample, ordered=TRUE)) %>%
        spread(key=sample, value=tfiib_counts) %>%
        remove_rownames() %>%
        column_to_rownames(var="tss_id") %>%
        as.data.frame()

    coldata = data.frame(condition=factor(c(rep(control_id, length(controls)),
                                            rep(condition_id, length(conditions))),
                                          levels = c(control_id, condition_id)),
                         row.names=names(countdata))

    dds = DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ condition)
    dds = dds %>%
        estimateSizeFactors() %>%
        estimateDispersions() %>%
        nbinomWaldTest()

    ncountsavg = dds %>%
        counts(normalized=TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var='name') %>%
        as_tibble() %>%
        gather(sample, value, -name) %>%
        separate(sample, into=c("group", "replicate"), sep="_") %>%
        mutate(group = if_else(group==condition_id, "condition_expr", "control_expr")) %>%
        group_by(name, group) %>% summarise(mean = mean(value)) %>% spread(group, mean) %>%
        ungroup()

    #extract DESeq2 results and write to file
    resdf = results(dds, alpha=alpha, lfcThreshold=lfc, altHypothesis="greaterAbs") %>%
        as_data_frame() %>%
        rownames_to_column(var='name') %>%
        inner_join(ncountsavg, by='name') %>%
        mutate_at(c('pvalue','padj'), funs(-log10(.))) %>%
        mutate_if(is.numeric, round, 3) %>%
        dplyr::rename(tfiib_lfc=log2FoldChange, tfiib_stat = stat, tfiib_lfcSE=lfcSE,
                      tfiib_logp=pvalue, tfiib_logq=padj, tfiib_meanLevels=baseMean,
                      tfiib_cond_levels=condition_expr, tfiib_ctrl_levels=control_expr) %>%
        left_join(df %>% select(-c(tfiib_counts, group, replicate)) %>% distinct(),
                  ., by=c("tss_id"="name")) %>%
        arrange(desc(tss_logq)) %>%
        write_tsv(path=results_all, col_names=TRUE)
}

main(controls = snakemake@input[["controls"]],
     conditions=snakemake@input[["conditions"]],
     condition_id= snakemake@wildcards[["condition"]],
     control_id= snakemake@wildcards[["control"]],
     results_all= snakemake@output[[1]],
     alpha= snakemake@params[["alpha"]],
     lfc= snakemake@params[["lfc"]])

