library(tidyverse)
library(viridis)
library(ggrepel)

main = function(in_path, group, roc_out, roc_zoom_out){
    df = read_tsv(in_path) %>%
        mutate(FPR=FP/(FP+TN),
               TPR=TP/(TP+FN),
               dist_to_corner = sqrt((FPR^2+(1-TPR)^2)))

    theme_default = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black"),
              legend.position = c(0.75,0.3),
              legend.justification = c(0.5,0.5),
              legend.background = element_blank(),
              legend.text = element_text(size=10),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=12, face="plain"))

    roc = ggplot(data = df, aes(x=FPR, y=TPR)) +
        geom_hline(yintercept = 1, color="grey65") +
        geom_vline(xintercept = 0, color="grey65") +
        geom_abline(slope=1, intercept = 0, linetype="dashed", color="grey65") +
        geom_line(size=1, lineend="round", color="black") +
        geom_point(aes(color=log10(cutoff+0.01)), size=2) +
        scale_color_viridis(guide=guide_colorbar(title=expression(bold(log[10] ~ cutoff)),
                                                 barheight = 7)) +
        theme_light() +
        scale_x_continuous(limits = c(0,1),
                           name = "false positive rate") +
        scale_y_continuous(limits = c(0,1),
                           name = "true positive rate") +
        ggtitle(paste(group, "ROC curve:"),
                subtitle = "TFIIB \'signal regions\' vs. MACS2 peaks as ground truth") +
        theme_default

    ggsave(roc_out, plot=roc, width=16, height=12, units="cm")

    roc_zoom = ggplot(data = df %>% filter(FPR<=0.2, TPR>=0.85), aes(x=FPR, y=TPR)) +
        geom_hline(yintercept = 1, color="grey65") +
        geom_vline(xintercept = 0, color="grey65") +
        geom_abline(slope=1, intercept = 0, linetype="dashed") +
        geom_line(size=1, lineend="round", color="black") +
        geom_text_repel(aes(label=cutoff),
                        box.padding = unit(0.5, "pt"),
                        nudge_x=0.01, nudge_y=-0.008,
                        size=8/75*25.4) +
        geom_point(aes(color=cutoff), size=2) +
        scale_color_viridis(guide=guide_colorbar(title="cutoff",
                                                 barheight = 8)) +
        theme_light() +
        scale_x_continuous(limits = c(0,0.2),
                           name = "false positive rate") +
        scale_y_continuous(limits = c(0.85,1),
                           name = "true positive rate") +
        ggtitle(paste(group, "ROC curve:"),
                subtitle = "TFIIB \'signal regions\' vs. MACS2 peaks as ground truth") +
        theme_default

    ggsave(roc_zoom_out, plot=roc_zoom, width=16, height=12, units="cm")
}

main(in_path = snakemake@input[[1]],
     group = snakemake@wildcards[["group"]],
     roc_out = snakemake@output[["roc"]],
     roc_zoom_out = snakemake@output[["roc_zoom"]])
