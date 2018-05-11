library(tidyverse)
library(forcats)
library(gridExtra)

df = read_tsv(snakemake@input[["tfiib_rel_tss"]],
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

p1 = ggplot(data = df, aes(x=-distance+0.1)) +
    geom_density(bw=5) +
    facet_grid(order~category, scales="free_y") +
    ggtitle(snakemake@wildcards[["group"]]) +
    scale_x_continuous(limits = c(NA, 300),
                       name = "(-) displacement from TSS summit to TFIIB summit (nt)") +
    scale_y_continuous(breaks = scales::pretty_breaks(2))

df2 = read_tsv(snakemake@input[["tss_rel_tfiib"]],
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

p2 = ggplot(data = df2, aes(x=distance+0.1)) +
    geom_density(bw=5) +
    facet_grid(order~category, scales="free_y") +
    scale_x_continuous(limits = c(NA, 300),
                       name = "distance from TFIIB summit to TSS summit (bp)") +
    scale_y_continuous(breaks = scales::pretty_breaks(2))

p = arrangeGrob(p1, p2, ncol=1, heights = c(1.2,2))

ggsave(snakemake@output[[1]], plot=p, width=24, height=18, units="cm")
