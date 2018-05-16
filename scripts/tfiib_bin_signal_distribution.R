library(tidyverse)

main = function(in_path, group_id, binsize, q, out_path){
    df = read_tsv(in_path,
                  col_names=c('chrom', 'start', 'end', 'signal', 'bin_type')) %>%
        mutate(signal=signal/(end-start))

    background_quantile = df %>%
        filter(bin_type=="background") %>%
        pull(signal) %>%
        quantile(q)

    label = paste(q, "~ quantile ==", round(background_quantile, 3))

    plot = ggplot(data = df %>% filter(bin_type %in% c("background", "genic")),
           aes(x=signal+0.01 , color=bin_type)) +
        geom_vline(xintercept = background_quantile,
                   linetype="dashed") +
        geom_density(bw=0.02, aes(y=..density..)) +
        annotate(geom="text",
                 x=1.1*background_quantile,
                 y=1.5,
                 hjust=0,
                 size=12/75*25.4,
                 label= label, parse=TRUE) +
        scale_x_log10(name="TFIIB signal") +
        ggtitle(paste(group_id, "TFIIB ChIP-nexus signal in", binsize, "bp bins")) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.y = element_text(size=10, face="plain"),
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              legend.position = c(0.99,0.99),
              legend.justification = c(1,1),
              plot.title = element_text(size=12))

    ggsave(out_path, plot=plot, width=16, height=10, units="cm")
}

main(in_path = snakemake@input[[1]],
     group_id = snakemake@wildcards[["group"]],
     binsize = snakemake@params[["bin_size"]],
     q = snakemake@params[["quantile"]],
     out_path = snakemake@output[[1]])
