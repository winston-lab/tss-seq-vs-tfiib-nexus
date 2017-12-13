library(tidyverse)
library(forcats)
library(viridis)

main = function(intable, dist, out.mosaic, out.tssexpr, out.tsssize, out.tssidr,
                out.tfiiblvl, out.contour, out.scatter){
    #import and set factor levels
    df = read_tsv(intable, col_types=cols(tfiib_enrichment=col_double(), tfiib_qval=col_double()))
    df$condition = factor(df$condition, levels=c("YPD", "diamide", "WT-37", "spt6-37"), ordered=TRUE)
    df$category = fct_inorder(df$category, ordered=TRUE)
   
    #mark if there is a TFIIB peak for a TSS. If more than one, keep the most significant 
    df = df %>% group_by(chrom, tss_start, tss_end, strand, condition, category) %>% arrange(desc(tfiib_qval)) %>%
            slice(1) %>% ungroup() %>% mutate(match=if_else(tfiib_start==-1, "no", "yes"))
    df$match = factor(df$match, levels=c("yes", "no"), ordered=TRUE)
   
    #build df for mosaic plot 
    countdf = df %>% filter(category!="all") %>%
                group_by(condition, category, match) %>% count() %>%
                group_by(condition, category) %>% 
                mutate(ymax=cumsum(n), ymin=(cumsum(n))-n) %>%
                mutate_at(vars(ymin, ymax), funs(./max(ymax)))
    csx = countdf %>% group_by(condition, category) %>%
            summarise(classn=sum(n)) %>% 
            mutate(xmax=cumsum(classn), xmin=cumsum(classn)-classn) %>%
            mutate_at(vars(xmin,xmax), funs(./max(xmax))) %>%
            select(-classn)
    countdf = countdf %>% left_join(csx, by=c('condition', 'category'))
    
    mplot = ggplot() +
                geom_rect(data=countdf, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=match),
                          color="white", alpha=0.9, size=1) +
                geom_text(data = countdf, aes(x=(xmin+xmax)/2, y=(ymin+ymax)/2, label=n),
                          size=4, color="black", fontface="bold") +
                geom_text(data = (countdf %>% summarise(x=(max(xmax)+min(xmin))/2)),
                          aes(x=x, label=category), y=1.1, angle=30, hjust=0.1, size=3) +
                scale_fill_brewer(palette = "Set1", direction=-1, name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_continuous(limits = c(NA, max(countdf$xmax)*1.05), expand=c(0,0)) +
                scale_y_continuous(limits = c(NA, 1.30)) +
                ggtitle("Overlap of TSS-seq peaks and TFIIB ChIP-nexus peaks") +
                guides(fill=guide_legend(reverse=TRUE)) +
                facet_grid(condition~., switch="y") +
                theme_void() +
                theme(legend.text = element_text(size=10, face="bold", color="black"),
                      legend.title = element_text(size=12, face="bold"),
                      legend.position = c(1, .5),
                      legend.justification = "left",
                      legend.key.size = unit(1, "cm"),
                      strip.text = element_text(size=12, face="bold", color="black"),
                      strip.text.y = element_text(angle=180, hjust=1),
                      plot.margin = unit(c(.5, 4, 0, 0.25), "cm"),
                      plot.title = element_text(size=12, face="bold", color="black"),
                      plot.subtitle = element_text(size=10))
    
    ggsave(out.mosaic, plot=mplot, height=18, width=18, units="cm")
    
    #compare TSS expression level, facetted by match status
    df$category = fct_relevel(df$category, "all")
    tssexp = ggplot(data = df, aes(x=tss_expression, fill=match)) +
                geom_density(bw=.05, alpha=0.5) +
                scale_fill_brewer(palette="Set1", direction=-1,
                                  name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_log10(name="TSS expression level (normalized counts)") +
                facet_grid(condition~category, switch="both", scales="free_y") +
                ggtitle("TSS expression level distributions") +
                theme_bw() +
                theme(strip.placement="outside",
                      text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      axis.title.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1))
    
    df = df %>% mutate(size=tss_end-tss_start)
    cutoff = as.integer(quantile(df$size, 0.98))
    
    tsssize= ggplot(data = df %>% mutate_at(vars(size), funs(if_else(.>cutoff, cutoff, .))),
                    aes(x=size, fill=match)) +
                geom_density(bw=5, alpha=0.5) +
                scale_fill_brewer(palette="Set1", direction=-1,
                                  name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_continuous(name="TSS peak size (nt)") +
                facet_grid(condition~category, switch="both", scales="free_y") +
                ggtitle("TSS peak size distributions") +
                theme_bw() +
                theme(strip.placement="outside",
                      text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      axis.title.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1))
    
    tssidr = ggplot(data = df, aes(x=tss_idr, fill=match)) +
                geom_density(bw=0.05, alpha=0.5) +
                scale_fill_brewer(palette="Set1", direction=-1,
                                  name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_continuous(name="TSS peak IDR") +
                facet_grid(condition~category, switch="both", scales="free_y") +
                ggtitle("TSS peak IDR distributions") +
                theme_bw() +
                theme(strip.placement="outside",
                      text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      axis.title.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1))
    
    tfiiblvl = ggplot(data = df, aes(x=tfiib_levels, fill=match)) +
                geom_density(bw=.05, alpha=0.5) +
                scale_fill_brewer(palette="Set1", direction=-1,
                                  name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_log10(name="TFIIB levels (normalized counts)") +
                facet_grid(condition~category, switch="both", scales="free_y") +
                ggtitle("TFIIB level distributions") +
                theme_bw() +
                theme(strip.placement="outside",
                      text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      axis.title.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1))
    ggsave(out.tssexpr, plot=tssexp, width=24, height=12, units="cm")
    ggsave(out.tsssize, plot=tsssize, width=24, height=12, units="cm")
    ggsave(out.tssidr, plot=tssidr, width=24, height=12, units="cm")
    ggsave(out.tfiiblvl, plot=tfiiblvl, width=24, height=12, units="cm")
    
    df = df %>% mutate(yfacet=paste(condition, match, sep="-"))
    df$yfacet = factor(df$yfacet, levels = c("spt6-37-yes","spt6-37-no","WT-37-yes","WT-37-no",
                                             "diamide-yes","diamide-no","YPD-yes", "YPD-no"),
                       ordered=TRUE)
    
    contour = ggplot(data = df, aes(x=tss_expression, y=tfiib_levels, color=match)) +
                geom_smooth(method="lm", size=0.3) +
                geom_density_2d(size=0.15, alpha=1) +
                scale_color_brewer(palette="Set1", direction=-1,
                                   name=paste0("TFIIB peak\nwithin ", dist, "nt\n5′ of TSS?")) +
                scale_x_log10(name="TSS expression level (normalized counts)") +
                scale_y_log10(name="TFIIB level (normalized counts)") +
                facet_grid(condition~category, switch="both") +
                ggtitle("TFIIB levels vs. TSS expression levels") +
                theme_bw() +
                theme(text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      strip.placement = "outside",
                      strip.background=element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1))
    ggsave(out.contour, plot=contour, width=24, height=12, units="cm")
    
    scatter = ggplot(data = df, aes(x=tss_expression, y=tfiib_levels)) +
                geom_smooth(method="lm", size=0.5) +
                #geom_point(shape=16, alpha=0.6) +
                stat_bin_hex(geom="point", aes(color=log10(..count..)),
                             binwidth=c(0.05, 0.05), size=.05, stroke=0) +
                scale_color_viridis(option="inferno") +
                scale_x_log10(name="TSS expression level (normalized counts)") +
                scale_y_log10(name="TFIIB level (normalized counts)") +
                facet_grid(yfacet~category, switch="both") +
                ggtitle("TFIIB levels vs. TSS expression levels") +
                theme_bw() +
                theme(text = element_text(size=12, face="bold", color="black"),
                      axis.text = element_text(size=10),
                      strip.placement = "outside",
                      strip.background=element_blank(),
                      strip.text.y = element_text(angle=-180, hjust=1),
                      legend.position ="none")
    ggsave(out.scatter, plot=scatter, width=24, height=24, units="cm")
}

main(intable = snakemake@input[[1]],
     dist = snakemake@params[["dist"]],
     out.mosaic = snakemake@output[["mosaic"]],
     out.tssexpr = snakemake@output[["tss_expr"]],
     out.tsssize = snakemake@output[["tss_size"]],
     out.tssidr= snakemake@output[["tss_idr"]],
     out.tfiiblvl = snakemake@output[["tfiib_lvl"]],
     out.contour = snakemake@output[["contour"]],
     out.scatter = snakemake@output[["scatter"]])
