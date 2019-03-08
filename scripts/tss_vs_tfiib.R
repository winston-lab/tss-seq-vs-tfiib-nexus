library(tidyverse)
library(viridis)

import_table = function(path, category_id){
    read_tsv(path,
             col_types = cols(tss_start = col_integer(),
                              tss_end = col_integer(),
                              tss_score = col_integer(),
                              tss_lfc = col_double(),
                              tss_lfc_SE = col_double(),
                              tss_FDR = col_double(),
                              tss_expr_condition = col_double(),
                              tss_expr_control = col_double(),
                              tss_summit = col_integer(),
                              tfiib_summit = col_integer(),
                              tfiib_start = col_integer(),
                              tfiib_end = col_integer(),
                              tfiib_score = col_integer(),
                              tfiib_lfc = col_double(),
                              tfiib_lfc_SE = col_double(),
                              tfiib_FDR = col_double(),
                              tfiib_abundance_condition = col_double(),
                              tfiib_abundance_control = col_double(),
                              feature_distance = col_integer())) %>%
        mutate(category = category_id) %>%
        return()
}

import_windows = function(path, category_id){
    read_tsv(path,
             col_types = cols(tss_name = col_character(),
                              chrom = col_character(),
                              tss_start = col_integer(),
                              tss_end = col_integer(),
                              strand = col_character(),
                              tss_lfc = col_double(),
                              tss_lfc_SE = col_double(),
                              tss_fdr = col_double(),
                              tss_expr_condition = col_double(),
                              tss_expr_control = col_double(),
                              tss_summit = col_integer(),
                              window_start = col_double(),
                              window_end = col_double(),
                              tfiib_lfc = col_double(),
                              tfiib_lfc_SE = col_double(),
                              tfiib_fdr = col_double(),
                              tfiib_abundance_condition = col_double(),
                              tfiib_abundance_control = col_double())) %>%
        select(1:18) %>%
        mutate(category = category_id) %>%
        return()
}

theme = theme_light() +
    theme(text = element_text(color="black") ,
          axis.text = element_text(color="black"),
          axis.title = element_text(size=10),
          axis.title.y = element_text(angle=0, vjust=0.5, hjust=1),
          plot.title = element_text(size=8),
          strip.text = element_text(color="black", hjust=0, size=10),
          strip.background = element_blank(),
          legend.position = "none")

main = function(matched_tfiib_paths,
                matched_tss_paths,
                window_paths,
                condition_id,
                control_id,
                search_distance,
                window_size,
                fdr_cutoff_tss,
                fdr_cutoff_tfiib,
                tfiib_peaks_matched_to_tss_out,
                tss_peaks_matched_to_tfiib_out,
                tss_peaks_vs_tfiib_windows_out,
                matched_tfiib_mosaic_out,
                matched_tss_mosaic_out){
    df_matched_tfiib = import_table(matched_tfiib_paths[1], "genic") %>%
        mutate(feature_distance = NA) %>%
        bind_rows(import_table(matched_tfiib_paths[2], "intragenic")) %>%
        bind_rows(import_table(matched_tfiib_paths[3], "antisense")) %>%
        bind_rows(import_table(matched_tfiib_paths[4], "convergent")) %>%
        bind_rows(import_table(matched_tfiib_paths[5], "divergent")) %>%
        bind_rows(import_table(matched_tfiib_paths[6], "intergenic")) %>%
        mutate(category = fct_inorder(category, ordered=TRUE)) %>%
        group_by(category, tss_chrom, tss_start, tss_end, tss_name) %>%
        arrange(abs((tss_start + tss_summit)-((tfiib_start + tfiib_end)/2)),
                .by_group=TRUE) %>%
        slice(1) %>%
        ungroup() %>%
        distinct(tss_name, tfiib_name, category, .keep_all = TRUE) %>%
        filter(! is.na(tss_FDR))

    df_matched_tss = import_table(matched_tss_paths[1], "genic") %>%
        mutate(feature_distance=NA) %>%
        bind_rows(import_table(matched_tss_paths[2], "intragenic")) %>%
        bind_rows(import_table(matched_tss_paths[3], "intergenic")) %>%
        mutate(category = fct_inorder(category, ordered=TRUE)) %>%
        group_by(category, tfiib_chrom, tfiib_start, tfiib_end, tfiib_name) %>%
        arrange(abs((tfiib_start + tfiib_summit)-((tss_start + tss_end)/2)),
                .by_group=TRUE) %>%
        slice(1) %>%
        ungroup() %>%
        distinct(tfiib_name, tss_name, category, .keep_all = TRUE) %>%
        filter(! is.na(tfiib_FDR))

    df_matched_tfiib_count = df_matched_tfiib %>%
        filter(tfiib_start >= 0) %>%
        count(category)

    df_matched_tss_count = df_matched_tss %>%
        filter(tss_start >= 0) %>%
        count(category)

    df_windows = import_windows(window_paths[1], "genic") %>%
        bind_rows(import_windows(window_paths[2], "intragenic")) %>%
        bind_rows(import_windows(window_paths[3], "antisense")) %>%
        bind_rows(import_windows(window_paths[4], "convergent")) %>%
        bind_rows(import_windows(window_paths[5], "divergent")) %>%
        bind_rows(import_windows(window_paths[6], "intergenic")) %>%
        filter(! is.na(tss_fdr)) %>%
        mutate(category = fct_inorder(category, ordered=TRUE))

    df_windows_count = df_windows %>%
        count(category)

    tfiib_peaks_matched_to_tss = ggplot(data = df_matched_tfiib) +
        geom_vline(xintercept=0,
                   color="grey65") +
        geom_hline(yintercept = 0,
                   color="grey65") +
        geom_segment(aes(x=tss_lfc - tss_lfc_SE,
                         xend=tss_lfc + tss_lfc_SE,
                         y=tfiib_lfc,
                         yend=tfiib_lfc),
                     alpha=0.2,
                     size=0.1) +
        geom_segment(aes(x=tss_lfc,
                         xend=tss_lfc,
                         y=tfiib_lfc - tfiib_lfc_SE,
                         yend=tfiib_lfc + tfiib_lfc_SE),
                     alpha=0.2,
                     size=0.1) +
        stat_bin_hex(geom="point",
                     binwidth = c(0.1, 0.1),
                     aes(x=tss_lfc,
                         y=tfiib_lfc,
                         color=..count..),
                     shape=16,
                     alpha=0.6,
                     size=0.3) +
        geom_label(data = df_matched_tfiib_count,
                   aes(label = paste0("n=", n)),
                   x = min(df_matched_tfiib[["tss_lfc"]] - df_matched_tfiib[["tss_lfc_SE"]], na.rm=TRUE),
                   y = max(df_matched_tfiib[["tfiib_lfc"]] + df_matched_tfiib[["tfiib_lfc_SE"]], na.rm=TRUE),
                   vjust=1,
                   hjust=0,
                   size=8/72*25.4,
                   label.r = unit(0, "lines"),
                   label.padding = unit(0.1, "lines"),
                   label.size= NA) +
        facet_wrap(~category) +
        scale_color_viridis(option="inferno")  +
        scale_x_continuous(name = bquote("TSS-seq" ~ log[2] ~ textstyle(frac(.(condition_id), .(control_id))))) +
        scale_y_continuous(name = bquote(atop("TFIIB ChIP-nexus", log[2] ~ textstyle(frac(.(condition_id), .(control_id)))))) +
        ggtitle(paste("TSS-seq peaks with TFIIB ChIP-nexus peaks within", search_distance, "nt upstream")) +
        theme

    tss_peaks_matched_to_tfiib = ggplot(data = df_matched_tss) +
        geom_vline(xintercept=0,
                   color="grey65") +
        geom_hline(yintercept = 0,
                   color="grey65") +
        geom_segment(aes(x=tfiib_lfc - tfiib_lfc_SE,
                         xend=tfiib_lfc + tfiib_lfc_SE,
                         y=tss_lfc,
                         yend=tss_lfc),
                     alpha=0.2,
                     size=0.1) +
        geom_segment(aes(x=tfiib_lfc,
                         xend=tfiib_lfc,
                         y=tss_lfc - tss_lfc_SE,
                         yend=tss_lfc + tss_lfc_SE),
                     alpha=0.2,
                     size=0.1) +
        stat_bin_hex(geom="point",
                     binwidth = c(0.1, 0.1),
                     aes(x=tfiib_lfc,
                         y=tss_lfc,
                         color=..count..),
                     shape=16,
                     alpha=0.6,
                     size=0.3) +
        geom_label(data = df_matched_tss_count,
                   aes(label = paste0("n=", n)),
                   x = min(df_matched_tfiib[["tfiib_lfc"]] - df_matched_tfiib[["tfiib_lfc_SE"]], na.rm=TRUE),
                   y = max(df_matched_tss[["tss_lfc"]] + df_matched_tss[["tss_lfc_SE"]], na.rm=TRUE),
                   vjust=1,
                   hjust=0,
                   size=8/72*25.4,
                   label.r = unit(0, "lines"),
                   label.padding = unit(0.1, "lines"),
                   label.size= NA) +
        facet_wrap(~category, ncol=2, nrow=2) +
        scale_color_viridis(option="inferno")  +
        scale_x_continuous(name = bquote("TFIIB ChIP-nexus" ~ log[2] ~ textstyle(frac(.(condition_id), .(control_id))))) +
        scale_y_continuous(name = bquote(atop("TSS-seq", log[2] ~ textstyle(frac(.(condition_id), .(control_id)))))) +
        ggtitle(bquote("TFIIB ChIP-nexus peaks with TSS-seq peaks within" %+-% .(search_distance) ~ nt)) +
        theme

    tss_peaks_vs_tfiib_windows =  ggplot(data = df_windows) +
        geom_vline(xintercept=0,
                   color="grey65") +
        geom_hline(yintercept = 0,
                   color="grey65") +
        geom_segment(aes(x=tss_lfc - tss_lfc_SE,
                         xend=tss_lfc + tss_lfc_SE,
                         y=tfiib_lfc,
                         yend=tfiib_lfc),
                     alpha=0.2,
                     size=0.1) +
        geom_segment(aes(x=tss_lfc,
                         xend=tss_lfc,
                         y=tfiib_lfc - tfiib_lfc_SE,
                         yend=tfiib_lfc + tfiib_lfc_SE),
                     alpha=0.2,
                     size=0.1) +
        stat_bin_hex(geom="point",
                     binwidth = c(0.1, 0.1),
                     aes(x=tss_lfc,
                         y=tfiib_lfc,
                         color=..count..),
                     shape=16,
                     alpha=0.6,
                     size=0.3) +
        geom_label(data = df_windows_count,
                   aes(label = paste0("n=", n)),
                   x = min(df_windows[["tss_lfc"]] - df_windows[["tss_lfc_SE"]], na.rm=TRUE),
                   y = max(df_windows[["tfiib_lfc"]] + df_windows[["tfiib_lfc_SE"]], na.rm=TRUE),
                   vjust=1,
                   hjust=0,
                   size=8/72*25.4,
                   label.r = unit(0, "lines"),
                   label.padding = unit(0.1, "lines"),
                   label.size= NA) +
        facet_wrap(~category) +
        scale_color_viridis(option="inferno")  +
        scale_x_continuous(name = bquote("TSS-seq" ~ log[2] ~ textstyle(frac(.(condition_id), .(control_id))))) +
        scale_y_continuous(name = bquote(atop("TFIIB ChIP-nexus", log[2] ~ textstyle(frac(.(condition_id), .(control_id)))))) +
        ggtitle(paste("TSS-seq in peaks vs. TFIIB ChIP-nexus in window extending", window_size, "bp upstream of TSS peak")) +
        theme

    ggsave(tfiib_peaks_matched_to_tss_out, plot=tfiib_peaks_matched_to_tss,
           width=16*1.25, height=9*1.25, units="cm", dpi=326)
    ggsave(tss_peaks_matched_to_tfiib_out, plot=tss_peaks_matched_to_tfiib,
           width=16*1.25, height=9*1.25, units="cm", dpi=326)
    ggsave(tss_peaks_vs_tfiib_windows_out, plot=tss_peaks_vs_tfiib_windows,
           width=16*1.25, height=9*1.25, units="cm", dpi=326)

    df_matched_tfiib_mosaic = df_matched_tfiib %>%
        mutate(tss_differential = if_else((tss_FDR <= -log10(fdr_cutoff_tss) | is.na(tss_FDR)),
                                          "no significant change",
                                          if_else(tss_lfc >=0, "upregulated", "downregulated")),
               tss_differential = ordered(tss_differential,
                                          levels = c("downregulated",
                                                     "no significant change",
                                                     "upregulated"))) %>%
        group_by(tss_differential, category) %>%
        count(tfiib_match = (tfiib_start >= 0)) %>%
        ungroup() %>%
        arrange(category, tss_differential, tfiib_match) %>%
        mutate(cumsum_category = cumsum(n),
               cumsum_category_lag = lag(cumsum_category, default=0)) %>%
        group_by(category) %>%
        mutate(category_max = max(cumsum_category),
               category_min = min(cumsum_category_lag),
               cumsum_differential = cumsum(n),
               cumsum_differential_lag = lag(cumsum_differential, default=0)) %>%
        select(-c(cumsum_category, cumsum_category_lag)) %>%
        group_by(category, tss_differential) %>%
        mutate(differential_max = max(cumsum_differential) / (category_max-category_min),
               differential_min = min(cumsum_differential_lag) / (category_max-category_min)) %>%
        select(-c(cumsum_differential, cumsum_differential_lag)) %>%
        arrange(desc(tfiib_match), .by_group=TRUE) %>%
        mutate(cumsum_match = cumsum(n),
               cumsum_match_lag = lag(cumsum_match, default=0),
               match_max = cumsum_match / max(cumsum_match) * (category_max - category_min) + category_min,
               match_min = cumsum_match_lag / max(cumsum_match) * (category_max - category_min) + category_min) %>%
        select(-c(cumsum_match, cumsum_match_lag))

    df_matched_tss_mosaic = df_matched_tss %>%
        mutate(tfiib_differential = if_else((tfiib_FDR <= -log10(fdr_cutoff_tfiib) | is.na(tfiib_FDR)),
                                            "no significant change",
                                            if_else(tfiib_lfc >=0, "upregulated", "downregulated")),
               tfiib_differential = ordered(tfiib_differential,
                                            levels = c("downregulated",
                                                       "no significant change",
                                                       "upregulated"))) %>%
        group_by(tfiib_differential, category) %>%
        count(tss_match = (tss_start >= 0)) %>%
        ungroup() %>%
        arrange(category, tfiib_differential, tss_match) %>%
        mutate(cumsum_category = cumsum(n),
               cumsum_category_lag = lag(cumsum_category, default=0)) %>%
        group_by(category) %>%
        mutate(category_max = max(cumsum_category),
               category_min = min(cumsum_category_lag),
               cumsum_differential = cumsum(n),
               cumsum_differential_lag = lag(cumsum_differential, default=0)) %>%
        select(-c(cumsum_category, cumsum_category_lag)) %>%
        group_by(category, tfiib_differential) %>%
        mutate(differential_max = max(cumsum_differential) / (category_max-category_min),
               differential_min = min(cumsum_differential_lag) / (category_max-category_min)) %>%
        select(-c(cumsum_differential, cumsum_differential_lag)) %>%
        arrange(desc(tss_match), .by_group=TRUE) %>%
        mutate(cumsum_match = cumsum(n),
               cumsum_match_lag = lag(cumsum_match, default=0),
               match_max = cumsum_match / max(cumsum_match) * (category_max - category_min) + category_min,
               match_min = cumsum_match_lag / max(cumsum_match) * (category_max - category_min) + category_min) %>%
        select(-c(cumsum_match, cumsum_match_lag))

    df_matched_tfiib_mosaic_labels = df_matched_tfiib_mosaic %>%
        ungroup() %>%
        distinct(category, category_min, category_max) %>%
        transmute(category=category,
                  x = (category_min + category_max)/2)

    df_matched_tss_mosaic_labels = df_matched_tss_mosaic %>%
        ungroup() %>%
        distinct(category, category_min, category_max) %>%
        transmute(category=category,
                  x = (category_min + category_max)/2)

    matched_tfiib_mosaic = ggplot(data = df_matched_tfiib_mosaic) +
        geom_rect(aes(xmin = match_min, xmax = match_max,
                      ymin = differential_min, ymax = differential_max,
                      fill=tss_differential, alpha=tfiib_match)) +
        geom_vline(aes(xintercept = category_max),
                   color="white", size=0.5) +
        geom_text(aes(x = (match_min + match_max)/2,
                      y = (differential_min + differential_max)/2,
                      label = n),
                  size=4/72*25.4,
                  alpha=0.5) +
        scale_x_continuous(breaks = df_matched_tfiib_mosaic_labels[["x"]],
                           labels = df_matched_tfiib_mosaic_labels[["category"]]) +
        scale_fill_brewer(palette = "Set1",
                          guide=guide_legend(title=paste0("TSS peak status,\n",
                                                        condition_id,
                                                        " vs. ",
                                                        control_id),
                                             reverse=TRUE,
                                             override.aes = list(alpha=0.8))) +
        scale_alpha_manual(values = c(0.4, 0.8),
                           labels = c("False", "True"),
                           name = paste0("TFIIB ChIP-nexus peak within\n", search_distance, " nt upstream of TSS"),
                           guide=guide_legend(reverse=TRUE)) +
        theme_void() +
        theme(axis.text.x = element_text(color="black", size=8, angle=30, vjust=1, hjust=1,
                                         margin = margin(t=-7, unit="pt")),
              legend.title = element_text(size=8))

    matched_tss_mosaic = ggplot(data = df_matched_tss_mosaic) +
        geom_rect(aes(xmin = match_min, xmax = match_max,
                      ymin = differential_min, ymax = differential_max,
                      fill=tfiib_differential, alpha=tss_match)) +
        geom_vline(aes(xintercept = category_max),
                   color="white", size=0.5) +
        geom_text(aes(x = (match_min + match_max)/2,
                      y = (differential_min + differential_max)/2,
                      label = n),
                  size=4/72*25.4,
                  alpha=0.5) +
        scale_x_continuous(breaks = df_matched_tss_mosaic_labels[["x"]],
                           labels = df_matched_tss_mosaic_labels[["category"]]) +
        scale_fill_brewer(palette = "Set1",
                          guide=guide_legend(title=paste0("TFIIB ChIP-nexus peak status,\n",
                                                        condition_id,
                                                        " vs. ",
                                                        control_id),
                                             reverse=TRUE,
                                             override.aes = list(alpha=0.8))) +
        scale_alpha_manual(values = c(0.4, 0.8),
                           labels = c("False", "True"),
                           name = paste0("TSS-seq peak within\n", search_distance, " bp of TFIIB peak"),
                           guide=guide_legend(reverse=TRUE)) +
        theme_void() +
        theme(axis.text.x = element_text(color="black", size=8, angle=30, vjust=1, hjust=1,
                                         margin = margin(t=-7, unit="pt")),
              legend.title = element_text(size=8))

    ggsave(matched_tfiib_mosaic_out, plot = matched_tfiib_mosaic, width=16, height=9, units="cm")
    ggsave(matched_tss_mosaic_out, plot = matched_tss_mosaic, width=16, height=9, units="cm")
}

main(matched_tfiib_paths = snakemake@input[["matched_tfiib_paths"]],
     matched_tss_paths = snakemake@input[["matched_tss_paths"]] ,
     window_paths = snakemake@input[["window_paths"]],
     condition_id = snakemake@wildcards[["condition"]],
     control_id = snakemake@wildcards[["control"]],
     search_distance = snakemake@params[["search_distance"]],
     window_size = snakemake@params[["window_size"]],
     fdr_cutoff_tss = snakemake@params[["fdr_cutoff_tss"]],
     fdr_cutoff_tfiib = snakemake@params[["fdr_cutoff_tfiib"]],
     tfiib_peaks_matched_to_tss_out = snakemake@output[["tfiib_matched_to_tss"]],
     tss_peaks_matched_to_tfiib_out = snakemake@output[["tss_matched_to_tfiib"]],
     tss_peaks_vs_tfiib_windows_out = snakemake@output[["tss_vs_tfiib_windows"]],
     matched_tfiib_mosaic_out = snakemake@output[["matched_tfiib_mosaic"]],
     matched_tss_mosaic_out = snakemake@output[["matched_tss_mosaic"]])
