

extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")

p <- ggplot(extra.rna.dtw.long.filt, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Expr') + xlab('Passage') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="white", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text = element_text(color ="white")) +
  #facet_grid(.~trending) + 
  facet_wrap( trending ~., scales = "fixed", ncol = 2)+
  theme(panel.spacing = unit(2.5, "lines"))+
  theme(
    strip.text.x = element_text(
      size = 14,  face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 14, face = "bold.italic"
    )
  ) +
  ggtitle("rna-seq") +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'black', hjust = 0.5),
    axis.title.x = element_blank(),
    #axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(axis.text.x =  element_blank(),
        axis.ticks.x =  element_blank())+
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))

plot(p)

ggsave(filename="../Output/compScBdTgPb/grant/lab_adapt_rna_up_down_trendingV2.png", 
       plot=p,
       width = 6, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p <- ggplot(extra.atac.dtw.long.filt, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Expr') + xlab('Passage') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="white", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text = element_text(color ="white")) +
  facet_grid(trending~.) + 
  theme(panel.spacing = unit(2, "lines"))+
  theme(
    strip.text.x = element_text(
      size = 14,  face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 14, face = "bold.italic"
    )
  ) +
  ggtitle("atac-seq") +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'black', hjust = 0.5),
    axis.title.x = element_blank(),
    #axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(axis.text.x =  element_blank(),
        axis.ticks.x =  element_blank() )+
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))

plot(p)



ggsave(filename="../Output/compScBdTgPb/grant/lab_adapt_atac_up_down_trending.png", 
       plot=p,
       width = 3, height = 5, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

#####

plot.atac.trend <- function(gene.id, gene.name,tc.logCPM){
  
  tmp <- tc.logCPM %>% dplyr::filter(GeneID == gene.id) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y,size = 2), color = 'cyan') + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + 
    theme_bw(base_size = 14) +
    ylab('Expr') + xlab('Passage') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #ggtitle(paste("atac", paste(gene.name, gene.id, sep = ':'), sep = ":")) +
    #ggtitle(paste("atac-seq", gene.name, sep = ": "))+
    ggtitle(gene.name)+
    theme(
      plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
      axis.title.y = element_blank(),
      #axis.title.y = element_text(size=14, face="bold"),
      legend.position = 'none'
    )+
    theme(axis.text.x = element_text(color ="black"), 
          axis.text.y = element_text(color = "black"),
          axis.ticks.x = element_blank())
  
  
  return(p)
  
}


plot.rna.trend <- function(gene.id, gene.name ,tc.logCPM){
  
  tmp <- tc.logCPM %>% dplyr::filter(GeneID == gene.id) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y,color = "red",size = 1.7)) + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + 
    ylab('Expr') + xlab('Passage') +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #ggtitle(paste("rna", paste(gene.name, gene.id, sep = ':'), sep = ":")) +
    #ggtitle(paste("rna-seq", gene.name, sep = ": "))+
    ggtitle(gene.name)+
    theme(
      plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
      axis.title.y = element_blank(),
      #axis.title.y = element_text(size=14, face="bold"),
      legend.position = 'none'
    ) +
    theme(axis.text.x = element_text(color ="black"), 
          axis.text.y = element_text(color = "black"),
          axis.ticks.x = element_blank())
  
  
  
  return(p)
  
}


extra.tc.logCPM.rna <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')
extra.tc.logCPM.atac <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')

# gene1 <- "TGGT1_299020" ; name <- "AP2III-4"
# gene <- "TGGT1_320680" ; name <- " AP2IV-2"
gene.x <- "TGGT1_268860" ; name.x <- "ENO1"
gene.y <- "TGGT1_268850" ;name.y <- "ENO2"

p1 <- plot.rna.trend(gene.id = gene.x, gene.name = name.x, tc.logCPM = extra.tc.logCPM.rna)
p2 <- plot.rna.trend(gene.id = gene.y, gene.name = name.y, tc.logCPM = extra.tc.logCPM.rna)

#p2 <- plot.atac.trend(gene.id = gene, gene.name = name,tc.logCPM = extra.tc.logCPM.atac) 

p1
p2
p <- do.call(grid.arrange, c(list(p1,p2), nrow=1))
ggsave(filename="../Output/compScBdTgPb/grant/ENO1_ENO2_rna_trend.png",
       plot=p,
       width = 6, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


# ggsave(filename="../Output/compScBdTgPb/grant/ENO1_ENO_atac_trend.png", 
#        plot=p,
#        width = 3, height = 3.5, 
#        units = "in", # other options are "in", "cm", "mm" 
#        dpi = 300
# )
