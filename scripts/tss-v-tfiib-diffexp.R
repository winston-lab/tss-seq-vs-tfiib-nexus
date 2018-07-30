library(tidyverse)
library(forcats)
library(viridis)

df = read_tsv(snakemake@input[[1]]) %>%
        mutate_at(vars(tfiib_meanExpr, tfiib_lfc, tfiib_logpadj,
                       tfiib_cond, tfiib_ctrl),
                  as.numeric) %>%
        mutate(group = paste(condition, "vs.", control))

df$group = fct_inorder(df$group, ordered=TRUE)
df$category = fct_inorder(df$category, ordered=TRUE) %>%
                fct_relevel("all")

#for TSSs with more than one TFIIB peak, take the TFIIB peak with most significant change
df = df %>% group_by(chrom, tss_start, tss_end, group, category) %>%
        arrange(desc(tfiib_logpadj)) %>% slice(1) %>% ungroup()

xmin = min(df$tss_lfc, na.rm=TRUE)
ymax = max(df$tfiib_lfc, na.rm=TRUE)

countdf = df %>% group_by(group, category) %>% count()

scatter = ggplot() +
            geom_vline(xintercept=0, color="grey85") +
            geom_hline(yintercept=0, color="grey85") +
            stat_bin_hex(data = df, aes(x=tss_lfc, y=tfiib_lfc, color=log10(..count..)),
                         geom="point", binwidth=c(0.1,0.1),
                         shape=16, size=0.5, alpha=0.8, stroke=0) +
            geom_text(data=countdf, aes(label=paste0("n=",n)), x=.98*xmin,y=.9*ymax, hjust=0, size=4) +
            scale_color_viridis(option="inferno", guide=FALSE) +
            xlab(expression(bold(paste("TSS-seq ", log[2]~frac(condition, control))))) +
            ylab(expression(bold(paste("TFIIB ChIP-nexus ", log[2]~frac(condition, control))))) +
            facet_grid(group~category) +
            theme_bw() +
            theme(legend.position="none",
                  strip.background = element_blank(),
                  strip.text = element_text(size=12, face="bold", color="black"),
                  axis.text = element_text(size=10),
                  axis.title = element_text(size=12, face="bold"))

ggsave(snakemake@output[["scatter"]], scatter, width=28, height=12, units="cm")
