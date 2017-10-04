library(tidyverse)
library(gridExtra)

tss = read_tsv(snakemake@input[["tss"]])
tfiib = read_tsv(snakemake@input[["tfiib"]])

xmax = tss %>% select(-1) %>% max()
ymax = tfiib %>% select(-1) %>% max()

nsamples = ncol(tss)-1
plotlist = list()

for (i in 1:nsamples){
    sname = names(tss)[i+1]
    df = tss %>% select(c(1,i+1)) %>%
            inner_join(tfiib %>% select(c(1,i+1)), by="gene")
    names(df) = c("gene","tss","tfiib")
    df = df %>% filter(tss>0 & tfiib>0) %>%
            mutate_at(c("tss","tfiib"), funs(.+.1))
    
    r = cor(df$tss, df$tfiib)
    
    plot = ggplot(data = df, aes(x=tss, y = tfiib)) +
            geom_point(alpha=0.3, stroke=0, size=.8) +
            #geom_hex(aes(fill=log10(..count..)),bins=50) +
            #scale_fill_viridis(guide="none") +
            geom_smooth(method="lm", color="red", size=1) +
            annotate("text", x=1e3, y=10, label=paste("R=", sprintf("%.2f", r)))+
            scale_x_log10(limits = c(NA, xmax), name="TSS-seq signal") +
            scale_y_log10(limits = c(NA, ymax), name="TFIIB ChIP-nexus signal") +
            ggtitle(sname) +
            theme_light() +
            theme(text = element_text(size=8))
    
    plotlist[[i]] = plot
}

plots = grid.arrange(grobs = plotlist, nrow=2)
ggsave(snakemake@output[["plot"]], plot=plots, width=28, height=14, unit="cm")
