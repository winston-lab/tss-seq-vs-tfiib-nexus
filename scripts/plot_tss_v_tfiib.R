library(tidyverse)
library(forcats)
library(viridis)
library(grid)
library(gridExtra)

import = function(df, path, category){
    df = read_tsv(path) %>%
        mutate(category=category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(all, genic, intragenic, antisense, convergent,
                divergent, intergenic, condition, control, out_path){
    df = tibble() %>%
        import(all, 'all') %>%
        import(genic, 'genic') %>%
        import(intragenic, 'intragenic') %>%
        import(antisense, 'antisense') %>%
        import(convergent, 'convergent') %>%
        import(divergent, 'divergent') %>%
        import(intergenic, 'intergenic') %>%
        mutate(category=fct_inorder(category, ordered=TRUE))

    count_df = df %>% count(category)

    theme_default = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black"),
              strip.background = element_blank(),
              strip.text = element_text(size=12, color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))

    plot_all = ggplot(data = df %>% filter(category == "all")) +
        geom_hline(yintercept = 0, color="grey65") +
        geom_vline(xintercept = 0, color="grey65") +
        geom_abline(intercept = 0, slope=1, color="grey65") +
        stat_bin_hex(geom="point",
                     aes(x=tss_lfc, y=tfiib_lfc, color=(..count..)),
                     binwidth=c(0.08, 0.08), alpha=0.6, size=0.8, shape=16) +
        geom_text(data = count_df %>% filter(category == "all"),
                  aes(label=paste0("n=",n)),
                  x=-5, y=4, hjust=0, size=12/75*25.4) +
        scale_color_viridis(guide=FALSE) +
        scale_fill_viridis(guide=FALSE) +
        scale_y_continuous(limits = c(-4.5, 6),
                           name=bquote(bold(atop("TFIIB ChIP-nexus", log[2] ~ frac(.(condition), .(control)))))) +
        scale_x_continuous(limits = c(-6, 9),
                           name=bquote(bold("TSS-seq" ~ log[2] ~ textstyle(frac(.(condition), .(control)))))) +
        theme_default

    facet_category = ggplot(data = df %>% filter(category != "all")) +
        geom_hline(yintercept = 0, color="grey65") +
        geom_vline(xintercept = 0, color="grey65") +
        geom_abline(intercept = 0, slope=1, color="grey65") +
        stat_bin_hex(geom="point",
                     aes(x=tss_lfc, y=tfiib_lfc, color=(..count..)),
                     binwidth=c(0.105, 0.105), alpha=0.5, size=0.6, shape=16) +
        geom_text(data = count_df %>% filter(category != "all"),
                  aes(label=paste0("n=",n)),
                  x=-4.5, y=4, hjust=0, size=12/75*25.4) +
        facet_wrap(~category, nrow=2) +
        scale_color_viridis(guide=FALSE) +
        scale_fill_viridis(guide=FALSE) +
        scale_y_continuous(limits = c(-4.5, 6),
                           name=bquote(bold(atop("TFIIB ChIP-nexus", log[2] ~ frac(.(condition), .(control)))))) +
        scale_x_continuous(limits = c(-6, 9),
                           name=bquote(bold("TSS-seq" ~ log[2] ~ textstyle(frac(.(condition), .(control)))))) +
        theme_default

    plot= arrangeGrob(arrangeGrob(nullGrob(), plot_all, nullGrob(), nrow=1, widths=c(0.15, 0.7, 0.15)),
                      facet_category, ncol=1, heights=c(0.7, 1))

    ggsave(out_path, plot=plot, width=28, height=28, units="cm")
}

main(all = snakemake@input[["all"]],
     genic = snakemake@input[["genic"]],
     intragenic = snakemake@input[["intragenic"]],
     antisense = snakemake@input[["antisense"]],
     convergent = snakemake@input[["convergent"]],
     divergent = snakemake@input[["divergent"]],
     intergenic = snakemake@input[["intergenic"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     out_path = snakemake@output[[1]])
