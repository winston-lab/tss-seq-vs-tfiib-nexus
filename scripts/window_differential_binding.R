library(tidyverse)
library(magrittr)
library(DESeq2)
# library(viridis)
# library(scales)
# library(gridExtra)

get_countdata = function(path, samples){
    df = read_tsv(path) %>%
        select(samples) %>%
        rowid_to_column(var="index") %>%
        column_to_rownames(var="index") %>%
        as.data.frame()
    df = df[rowSums(df)>1,]
    return(df)
}

initialize_dds = function(data_path,
                          samples,
                          groups,
                          condition_id,
                          control_id){
    DESeqDataSetFromMatrix(countData = get_countdata(data_path, samples),
                           colData = data.frame(condition = factor(groups,
                                                                   levels = c(control_id,
                                                                              condition_id)),
                                                row.names = samples),
                           design = ~ condition) %>%
        return()
}

extract_normalized_counts = function(dds){
    dds %>%
        counts(normalized=TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var="index") %>%
        as_tibble() %>%
        return()
}

extract_rlog_counts = function(dds){
    dds %>%
        rlog(blind=FALSE) %>%
        assay() %>%
        as.data.frame() %>%
        rownames_to_column(var="index") %>%
        as_tibble() %>%
        return()
}

# build_mean_sd_df_pre = function(dds){
#      dds %>%
#         normTransform() %>%
#         assay() %>%
#         as_tibble() %>%
#         rowid_to_column(var="index") %>%
#         gather(sample, signal, -index) %>%
#         group_by(index) %>%
#         summarise(mean = mean(signal),
#                   sd = sd(signal)) %>%
#         mutate(rank = min_rank(desc(mean))) %>%
#         return()
# }
#
# build_mean_sd_df_post = function(counts){
#     counts %>%
#         gather(sample, signal, -index) %>%
#         group_by(index) %>%
#         summarise(mean = mean(signal),
#                   sd = sd(signal)) %>%
#         mutate(rank = min_rank(desc(mean))) %>%
#         return()
# }
#
# reverselog_trans <- function(base = exp(1)) {
#     trans <- function(x) -log(x, base)
#     inv <- function(x) base^(-x)
#     trans_new(paste0("reverselog-", format(base)), trans, inv,
#               log_breaks(base = base),
#               domain = c(1e-100, Inf))
# }
#
# mean_sd_plot = function(df, ymax, title){
#     ggplot(data = df, aes(x=rank, y=sd)) +
#         geom_hex(aes(fill=log10(..count..), color=log10(..count..)), bins=100, size=0) +
#         geom_smooth(color="#4292c6") +
#         scale_fill_viridis(option="inferno", name=expression(log[10](count)), guide=FALSE) +
#         scale_color_viridis(option="inferno", guide=FALSE) +
#         scale_x_continuous(trans = reverselog_trans(10),
#                            name="rank(mean expression)",
#                            expand = c(0,0)) +
#         scale_y_continuous(limits = c(NA, ymax),
#                            name = "SD") +
#         theme_light() +
#         ggtitle(title) +
#         theme(text = element_text(size=8))
# }

get_mean_counts = function(counts_table,
                           samples,
                           groups,
                           condition_id){
    counts_table %>%
        # select(-c(1:6)) %>%
        # rownames_to_column(var = "index") %>%
        gather(sample, value, -index) %>%
        mutate(group = if_else(sample %in% samples[groups==condition_id],
                               "condition_expr",
                               "control_expr")) %>%
        group_by(index, group) %>%
        summarise(mean = mean(value)) %>%
        spread(group, mean) %>%
        ungroup() %>%
        return()
}

extract_deseq_results = function(dds,
                                 annotations,
                                 mean_counts_table,
                                 alpha,
                                 lfc){
    results(dds,
            alpha=alpha,
            lfcThreshold=lfc,
            altHypothesis="greaterAbs",
            tidy=TRUE) %>%
        as_tibble() %>%
        dplyr::rename(index = row) %>%
        left_join(annotations, ., by="index") %>%
        left_join(mean_counts_table, ., by="index") %>%
        arrange(padj) %>%
        mutate(name = if_else(tss_name==".",
                              paste0("peak_", row_number()),
                              tss_name),
               score = as.integer(pmin(-125*log2(padj), 1000))) %>%
        mutate_at(vars(pvalue, padj), funs(-log10(.))) %>%
        mutate_if(is.double, round, 3) %>%
        select(index, chrom, tss_start, tss_end, tss_name, strand,
               tss_lfc, tss_lfc_SE, tss_fdr, tss_expr_condition, tss_expr_control, tss_summit,
               window_start, window_end,
               tfiib_lfc=log2FoldChange, tfiib_lfc_SE=lfcSE,
               tfiib_fdr=padj,
               tfiib_abundance_condition=condition_expr,
               tfiib_abundance_control=control_expr) %>%
        return()
}

write_counts_table = function(results_df,
                              annotations,
                              counts_df,
                              output_path){
    results_df %>%
        select(index, chrom, window_start, window_end, tss_name, strand) %>%
        right_join(annotations %>% select(c(index, chrom, window_start, window_end, strand)),
                   by = c("index", "chrom", "window_start", "window_end", "strand")) %>%
        left_join(counts_df, by="index") %>%
        select(-index) %>%
        write_tsv(output_path) %>%
        return()
}

# plot_ma = function(df_sig = results_df_filtered_significant,
#                    df_nonsig = results_df_filtered_nonsignificant,
#                    xvar = mean_expr,
#                    yvar = log2_foldchange,
#                    lfc,
#                    condition,
#                    control){
#     xvar = enquo(xvar)
#     yvar = enquo(yvar)
#     ggplot() +
#         geom_hline(yintercept = 0, color="black", linetype="dashed") +
#         geom_hline(yintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
#         geom_point(data = df_nonsig,
#                    aes(x=!!xvar, y=!!yvar),
#                    color="black", alpha=0.3, stroke=0, size=0.7) +
#         geom_point(data = df_sig,
#                    aes(x=!!xvar, y=!!yvar),
#                    color="red", alpha=0.3, stroke=0, size=0.7) +
#         scale_x_log10(name="mean of normalized counts") +
#         ylab(bquote(log[2]~frac(.(condition),.(control)))) +
#         theme_light() +
#         theme(text = element_text(size=8, color="black"),
#               axis.text = element_text(color = "black"),
#               axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
# }
#
# plot_volcano = function(df = results_df_filtered,
#                         xvar = log2_foldchange,
#                         yvar = log10_padj,
#                         lfc,
#                         alpha,
#                         condition,
#                         control){
#     xvar = enquo(xvar)
#     yvar = enquo(yvar)
#     ggplot() +
#         geom_vline(xintercept = 0, color="black", linetype="dashed") +
#         geom_vline(xintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
#         geom_point(data = df,
#                    aes(x = !!xvar, y = !!yvar),
#                    alpha=0.3, stroke=0, size=0.7) +
#         geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
#         xlab(bquote(log[2] ~ frac(.(condition), .(control)))) +
#         ylab(expression(-log[10] ~ FDR)) +
#         theme_light() +
#         theme(text = element_text(size=8),
#               axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
# }

main = function(exp_table,
                spike_table,
                samples,
                groups,
                norm,
                condition,
                control,
                alpha,
                lfc,
                counts_norm_out,
                counts_rlog_out,
                results_all_out){

    annotations = read_tsv(exp_table) %>%
        select(1:13) %>%
        rownames_to_column(var="index") %>%
        mutate(chrom = str_replace(chrom, "-minus$|-plus$", ""))

    # initialize DESeq datasets
    dds = initialize_dds(data_path = exp_table,
                         samples = samples,
                         groups = groups,
                         condition_id = condition,
                         control_id = control)

    if (norm=="spikenorm"){
        dds_spike = initialize_dds(data_path = spike_table,
                                   samples = samples,
                                   groups = groups,
                                   condition_id = condition,
                                   control_id = control) %>%
            estimateSizeFactors()

        sizeFactors(dds) = sizeFactors(dds_spike)
    } else {
        dds %<>% estimateSizeFactors()
    }
    dds %<>% estimateDispersions() %>% nbinomWaldTest()

    #extract normalized counts and write to file
    counts_norm = extract_normalized_counts(dds = dds)
    counts_rlog = extract_rlog_counts(dds = dds)

    mean_counts_norm = get_mean_counts(counts_table = counts_norm,
                                       samples = samples,
                                       groups = groups,
                                       condition_id = condition)

    results_df = extract_deseq_results(dds = dds,
                                       annotations = annotations,
                                       mean_counts_table = mean_counts_norm,
                                       alpha = alpha,
                                       lfc = lfc) %>%
        mutate(chrom = str_replace(chrom, "-minus$|-plus$", ""))

    write_counts_table(results_df = results_df,
                       annotations = annotations,
                       counts_df = counts_norm,
                       output_path = counts_norm_out)
    write_counts_table(results_df = results_df,
                       annotations = annotations,
                       counts_df = counts_rlog,
                       output_path = counts_rlog_out)

    results_df %<>%
        select(-index) %>%
        write_tsv(results_all_out)
}

main(exp_table = snakemake@input[["exp_counts"]],
     spike_table = snakemake@input[["spike_counts"]],
     samples = snakemake@params[["samples"]],
     groups = snakemake@params[["groups"]],
     norm = snakemake@wildcards[["norm"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     alpha = snakemake@params[["alpha"]],
     lfc = snakemake@params[["lfc"]],
     counts_norm_out = snakemake@output[["counts_norm"]],
     counts_rlog_out = snakemake@output[["counts_rlog"]],
     results_all_out = snakemake@output[["results_all"]])

