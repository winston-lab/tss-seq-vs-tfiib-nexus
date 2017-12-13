library(tidyverse)
library(viridis)

axislabel = function(up, dn, xlab){
    if (up < (1/3)*dn){
        return( c('', xlab, (dn/1000)))
    }
    else if (dn < (1/3)*up){
        return( c((-up/1000), xlab, ''))
    }
    else {
        return( c((-up/1000), xlab, (dn/1000)))
    }
}

main = function(intable, tfiib_cutoff, condition, upstream, dnstream, refptlab, cmap, outpath){
    tfiib = read_tsv(intable, col_names=c("type","condition","index","position","signal"))
    
    nanno = tfiib %>% group_by(type) %>% summarise(n=max(index))
    
    #percentile cutoff for heatmap visualization
    tfiib_cutoff = quantile(tfiib$signal, probs=tfiib_cutoff, na.rm=TRUE)
    
    heatmap = ggplot() +
                geom_raster(data = tfiib %>% mutate_at(vars(signal), funs(pmin(tfiib_cutoff, .))),
                            aes(x=position, y=index, fill=signal)) +
                facet_grid(type~., scales="free_y", space="free_y", switch="y",
                           labeller = labeller(type=c(with=paste(nanno$n[1], "TSSs\nwith TFIIB peaks"),
                                                      without=paste(nanno$n[2], "TSSs\nwithout TFIIB peaks")))) +
                scale_y_reverse(expand=c(0.01, 0)) +
                scale_x_continuous(breaks = c(-upstream/1000, 0, dnstream/1000),
                                   labels= axislabel(up=upstream, dn=dnstream, xlab=refptlab), 
                                   minor_breaks = scales::pretty_breaks(n=10),
                                   name=paste("distance from", refptlab, "(kb)")) +
                scale_fill_viridis(option = cmap, na.value="FFFFFF00",
                                   name=('TFIIB ChIP-nexus signal'),
                                   guide=guide_colorbar(title.position="top",
                                                        barwidth=13, barheight=0.8,
                                                        title.hjust=0.5)) +
                ggtitle(condition) +
                theme_minimal() +
                theme(text = element_text(size=12, face="bold", color="black"),
                        legend.position = "top",
                        legend.title = element_text(size=12, face="bold", color="black"),
                        legend.text = element_text(size=8, face="plain"),
                        strip.text = element_text(size=12, face="bold", color="black"),
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        strip.text.y = element_text(angle=-180, hjust=1), 
                        axis.text.x = element_text(size=12, face="bold", color="black",
                                                   margin = unit(c(0,0,0,0),"cm")),
                        panel.grid.major.x = element_line(color="black"),
                        panel.grid.minor.x = element_line(color="grey80"),
                        panel.grid.major.y = element_line(color="grey80"),
                        panel.grid.minor.y = element_blank(),
                        panel.spacing.x = unit(.5, "cm"))
    
    hmap.width = max(12, .0008*(upstream+dnstream)+3.4)
    hmap.height = sum(nanno$n)*.0009+11.5
    
    ggsave(outpath, heatmap, height=hmap.height, width=hmap.width, units="cm")
}

main(intable=snakemake@input[[1]],
     tfiib_cutoff=snakemake@params[["cutoff"]],
     condition=snakemake@wildcards[["condition"]],
     upstream=snakemake@params[["upstream"]],
     dnstream=snakemake@params[["dnstream"]],
     refptlab="TSS",
     cmap=snakemake@params[["cmap"]],
     outpath=snakemake@output[[1]])
