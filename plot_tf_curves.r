
library(gridExtra)
library(ggplot2)
library(gridExtra)


extra.tc.logCPM.rna <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')
source.tf.ids  <- unique(rna.tf.targ.info$source.tf.id)


ps <- lapply(1:length(genes), function(i){
  
  tmp <- extra.tc.logCPM.rna %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, color = 'red')) + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + 
    ylab('Expr') + xlab('Passage') +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 11, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 11, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste("rna",genes[i], sep = ":")) +
    theme(
      plot.title = element_text(size=11, face="bold"),
      axis.title.x = element_text(size=11, face="bold", hjust = 0.5),
      axis.title.y = element_text(size=11, face="bold"),
      legend.position = 'none'
    ) 
  #plot(p)
  
})

pdf("../Output/compScBdTgPb/fig_tfs/up_tfs_ecpresseion_curve_rna.pdf", onefile = TRUE, width = 9, height = 6)

do.call("grid.arrange", c(ps, ncol=4))

dev.off()



plot.rna.trend <- function(gene.id, gene.name ,tc.logCPM){
  
  tmp <- tc.logCPM %>% dplyr::filter(GeneID == gene.id) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, color = 'red')) + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + 
    ylab('Expr') + xlab('Passage') +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste("rna", paste(gene.name, gene.id, sep = ':'), sep = ":")) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
      axis.title.y = element_text(size=14, face="bold"),
      legend.position = 'none'
    ) 
  
  
  
  return(p)
  
}


gene <- "TGGT1_299020" ; name <- "AP2III-4"
gene <- "TGGT1_320680" ; name <- "AP2IV-2"
gene <- "TGGT1_268860" ; name <- "ENO1"
gene <- "TGGT1_268850" ; name <- "ENO2"
gene <- "TGGT1_200385" ; name <- "BFD1"
gene <- "TGGT1_311100" ; name <- "BFD2"
gene <- "TGGT1_293670" ; name <- "TFIIS"
gene <- "TGGT1_321450" ; name <- "Myb3"
gene <- "TGGT1_209030" ; name <- "actin ACT1"
gene <- "TGGT1_290180" ; name <- "AP2 domain transcription factor AP2IX-6"

p2 <- plot.rna.trend(gene.id = gene,gene.name = name, tc.logCPM = extra.tc.logCPM.rna)
