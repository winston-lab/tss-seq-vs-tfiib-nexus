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

