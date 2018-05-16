library(tidyverse)
library(forcats)
library(gridExtra)
library(ggthemes)

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black", face="bold"),
          axis.text = element_text(size=12, color="black"),
          axis.text.y = element_text(face="plain"),
          strip.background = element_blank(),
          strip.text = element_text(size=12, color="black"),
          strip.text.y = element_text(angle=0, hjust=0),
          plot.title = element_text(size=12, face="bold"))

main = function(tfiib_rel_tss_path, tss_rel_tfiib_path, group, max_dist,
                rel_distances_out, tss_signal_violin_out, tss_signal_density_out, tfiib_signal_violin_out,
                tfiib_signal_density_out, tfiib_rel_tss_mosaic_out, tss_rel_tfiib_mosaic_out){
    tfiib_rel_tss = read_tsv(tfiib_rel_tss_path,
                  col_names = c("chrom", "tss_start", "tss_end", "a", "tss_score", "tss_strand",
                                "b", "tss_pval", "tss_qval", "c", "tss_coverage",
                                "d", "tfiib_start", "tfiib_end", "e", "tfiib_score", "f",
                                "tfiib_enrichment", "tfiib_pval", "tfiib_qval", "g", "tfiib_coverage",
                                "distance", "category")) %>%
        select(-c(a,b,c,d,e,f,g)) %>%
        mutate(category = fct_inorder(category, ordered=TRUE)) %>%
        group_by(chrom, tss_start, tss_end, tss_score, tss_strand, category) %>%
        arrange(desc(distance), .by_group=TRUE) %>%
        mutate(order = row_number())

    tss_rel_tfiib = read_tsv(tss_rel_tfiib_path,
                  col_names = c("chrom", "tfiib_start", "tfiib_end", "a", "tfiib_score", "b",
                                "tfiib_enrichment", "tfiib_pval", "tfiib_qval", "c", "tfiib_coverage",
                                "d", "tss_start", "tss_end", "e", "tss_score", "tss_strand",
                                "f", "tss_pval", "tss_qval", "g", "tss_coverage",
                                "distance", "category")) %>%
        select(-c(a,b,c,d,e,f,g)) %>%
        mutate(distance = abs(distance),
               category = fct_inorder(category, ordered=TRUE)) %>%
        group_by(chrom, tfiib_start, tfiib_end, tfiib_score, category) %>%
        arrange(distance, .by_group=TRUE) %>%
        mutate(order = row_number()) %>%
        ungroup()

    tfiib_rel_tss_distances = ggplot(data = tfiib_rel_tss, aes(x=-distance+0.1)) +
        geom_density(bw=5, fill=ptol_pal()(1)) +
        facet_grid(order~category, scales="free_y") +
        ggtitle(group) +
        scale_x_continuous(limits = c(NA, 300),
                           name = "(-) displacement from TSS summit to TFIIB summit (nt)") +
        scale_y_continuous(breaks = scales::pretty_breaks(2)) +
        theme_default

    tss_rel_tfiib_distances = ggplot(data = tss_rel_tfiib, aes(x=distance+0.1)) +
        geom_density(bw=5, fill=ptol_pal()(1)) +
        facet_grid(order~category, scales="free_y") +
        scale_x_continuous(limits = c(NA, 300),
                           name = "distance from TFIIB summit to TSS summit (bp)") +
        scale_y_continuous(breaks = scales::pretty_breaks(2)) +
        theme_default

    rel_distances_plot = arrangeGrob(tfiib_rel_tss_distances, tss_rel_tfiib_distances,
                                ncol=1, heights = c(1.5,2))

    tfiib_rel_tss_first = tfiib_rel_tss %>%
        filter(order==1) %>%
        select(-order) %>%
        mutate(match = if_else(distance>=-max_dist, TRUE, FALSE))

    tss_signal_group_match_violin = ggplot(data = tfiib_rel_tss_first,
                                           aes(x=fct_rev(category), y=tss_coverage, fill=match)) +
        geom_violin(bw=.05, draw_quantiles = c(0.5)) +
        scale_y_log10(name="normalized TSS-seq counts") +
        coord_flip() +
        scale_fill_ptol(guide=guide_legend(reverse=TRUE)) +
        ggtitle(bquote(bold(.(group) ~ TSS ~ signal %+-% matching ~ TFIIB ~ peak))) +
        theme_default +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(face="bold"))

    tss_signal_group_match_density = ggplot(data = tfiib_rel_tss_first,
                                            aes(x=tss_coverage, color=match)) +
        geom_density(bw=.05) +
        scale_color_ptol(guide=guide_legend(reverse=TRUE)) +
        facet_wrap(~category, nrow=2) +
        scale_x_log10("normalized TSS-seq counts") +
        ggtitle(bquote(bold(.(group) ~ TSS ~ signal %+-% matching ~ TFIIB ~ peak))) +
        theme_default

    tss_rel_tfiib_first = tss_rel_tfiib %>%
        filter(order==1) %>%
        select(-order) %>%
        mutate(match = if_else(distance<=max_dist, TRUE, FALSE))

    tfiib_signal_group_match_violin = ggplot(data = tss_rel_tfiib_first,
                                             aes(x=fct_rev(category), y=tfiib_coverage, fill=match)) +
        geom_violin(bw=.05, draw_quantiles = c(0.5)) +
        scale_y_log10(name="normalized TFIIB ChIP-nexus counts") +
        coord_flip() +
        scale_fill_ptol(guide=guide_legend(reverse=TRUE)) +
        ggtitle(bquote(bold(.(group) ~ TFIIB ~ signal %+-% matching ~ "TSS-seq" ~ peak))) +
        theme_default +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(face="bold"))

    tfiib_signal_group_match_density = ggplot(data = tss_rel_tfiib_first,
                                            aes(x=tfiib_coverage, color=match)) +
        geom_density(bw=.05) +
        scale_color_ptol(guide=guide_legend(reverse=TRUE)) +
        facet_wrap(~category, nrow=1) +
        scale_x_log10("normalized TFIIB ChIP-nexus counts") +
        ggtitle(bquote(bold(.(group) ~ TFIIB ~ signal %+-% matching ~ "TSS-seq" ~ peak))) +
        theme_default

    get_mosaic_df = function(df){
        df %>%
            filter(category != "all") %>%
            ungroup() %>%
            count(category, match) %>%
            mutate(xmax=cumsum(n), xmin=cumsum(n)-n) %>%
            group_by(category) %>%
            arrange(desc(match), .by_group=TRUE) %>%
            mutate(ymax=cumsum(n), ymin=(cumsum(n)-n),
                   xmax=max(xmax), xmin=min(xmin)) %>%
            mutate_at(vars(ymin, ymax), funs(./max(ymax))) %>%
            ungroup() %>%
            mutate_at(vars(xmin, xmax), funs(./max(xmax))) %>%
            mutate(x=(xmin+xmax)/2,
                   y=(ymin+ymax)/2) %>%
            return()
    }

    tfiib_rel_tss_counts = get_mosaic_df(tfiib_rel_tss_first)
    tss_rel_tfiib_counts = get_mosaic_df(tss_rel_tfiib_first)

    mosaic = function(df, title){
        ggplot() +
            geom_rect(data = df,
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=match),
                      size=2, color="white") +
            geom_text(data = df,
                      aes(x=x, y=y, label=n),
                      size=12/72*25.4) +
            geom_text(data = df %>% group_by(category) %>% summarise(x=first(x)),
                          aes(x=x, y=0.02, angle=30, hjust=1, vjust=1, label=category),
                      size=12/72*25.4, fontface="bold") +
            scale_fill_ptol() +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(limits = c(-0.2, 1), expand=c(0,0)) +
            ggtitle(title) +
            theme_void() +
            theme(text = element_text(size=12, color="black", face="bold"),
                  plot.title = element_text(size=12, face="bold"))
    }

    tfiib_rel_tss_mosaic = mosaic(tfiib_rel_tss_counts, title=paste0(group, ": TFIIB peak status upstream of TSS-seq peaks"))
    tss_rel_tfiib_mosaic = mosaic(tss_rel_tfiib_counts, title=paste0(group, ": TSS-seq peak status around TFIIB peaks"))

    ggsave(rel_distances_out, rel_distances_plot, width=24, height=20, units="cm")
    ggsave(tss_signal_violin_out, tss_signal_group_match_violin, width=14, height=12, units="cm")
    ggsave(tss_signal_density_out, tss_signal_group_match_density, width=18, height=10, units="cm")
    ggsave(tfiib_signal_violin_out, tfiib_signal_group_match_violin, width=14, height=8, units="cm")
    ggsave(tfiib_signal_density_out, tfiib_signal_group_match_density, width=20, height=6, units="cm")
    ggsave(tfiib_rel_tss_mosaic_out, tfiib_rel_tss_mosaic, width=18, height=9, units="cm")
    ggsave(tss_rel_tfiib_mosaic_out, tss_rel_tfiib_mosaic, width=14, height=9, units="cm")
}

main(tfiib_rel_tss_path = snakemake@input[["tfiib_rel_tss"]],
     tss_rel_tfiib_path = snakemake@input[["tss_rel_tfiib"]],
     group = snakemake@wildcards[["group"]],
     max_dist = snakemake@params[["max_dist"]],
     rel_distances_out= snakemake@output[["peak_distances"]],
     tss_signal_violin_out= snakemake@output[["tss_signal_violin"]],
     tss_signal_density_out= snakemake@output[["tss_signal_density"]],
     tfiib_signal_violin_out= snakemake@output[["tfiib_signal_violin"]],
     tfiib_signal_density_out= snakemake@output[["tfiib_signal_density"]],
     tfiib_rel_tss_mosaic_out= snakemake@output[["tfiib_rel_tss_mosaic"]],
     tss_rel_tfiib_mosaic_out= snakemake@output[["tss_rel_tfiib_mosaic"]])
