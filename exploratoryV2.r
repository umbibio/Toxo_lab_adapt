library(dtwclust)
library(bigmemory)
library(doParallel)
library(tidyverse)

library(tidytext)
library(gridExtra)
library(gt)

library(grid)
#library(fda)
#library(sme)
library(openxlsx)
#library(tvReg)

source("./util_funcs.R")



## match RNA with ATAC

#extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.rna.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_rna_dtw_trending_2_clusters_less_stringent.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_atac_dtw_trending_2_clusters.rds")

p <- ggplot(extra.rna.dtw.long.filt, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Expr') + xlab('Passage') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_grid(trending~.) + 
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
    axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))

plot(p)

ggsave(filename="../Output/Toxo_lab_adapt/figs/lab_adapt_rna_up_down_trending.png", 
       plot=p,
       width = 3, height = 6, 
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
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_grid(trending~.) + 
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
    axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))

plot(p)

ggsave(filename="../Output/Toxo_lab_adapt/figs/lab_adapt_atac_up_down_trending.pdf", 
       plot=p,
       width = 3, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


# agreement between rna and atac (not much)
rna.atac <- inner_join(extra.rna.dtw.long.filt, extra.atac.dtw.long.filt , by = "GeneID") %>%
  distinct(GeneID, .keep_all =T) %>% transmute(GeneID = GeneID, trending.rna  = trending.x , trending.atac = trending.y)
rna.atac <- tibble(rna.atac)
rna.atac.agreement <- table(rna.atac$trending.rna, rna.atac$trending.atac)
rna.atac.agreement


# match genes in rna and atac

dim(extra.rna.dtw.long.filt)
dim(extra.atac.dtw.long.filt)

extra.rna.summ <- extra.rna.dtw.long.filt %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% 
  summarise(genes = list(unique(GeneID)), total = n(), .groups ='drop')

extra.atac.summ <- extra.atac.dtw.long.filt %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% 
  summarise(genes = list(unique(GeneID)), total = n())

matched.rna.atac <- inner_join(extra.rna.dtw.long.filt, extra.atac.dtw.long.filt, by = "GeneID") %>% 
  dplyr::select(GeneID, everything())

saveRDS( matched.rna.atac,"../Input/Toxo_lab_adapt/RDS/matched_extra_rna_atac_dtw_trending_2_clusters.rds")


# ATAC expr for corresponding genes in RNA clusters
p <- ggplot(matched.rna.atac, aes(x = x.y, y = y.y, group = GeneID)) + 
  geom_line(aes(x = x.y, y = y.y, color = factor(cluster.x)), alpha = 0.4) + 
  #geom_smooth(aes(x = x.y, y = y.y, group = NA), method='lm', color = 'black') + 
  theme_bw() +
  facet_wrap(.~cluster.x)+
  ggtitle("matched genes ATAC-seq") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))

plot(p)

ggsave(filename="../Output/Toxo_lab_adapt/figs/lab_adapt_matched_rna_atac_expr_2cluster.pdf", 
       plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## plot a few curves
## Plot trends

#clust4.genes <- matched.rna.atac %>% filter(cluster.x == 4) %>% distinct(GeneID,  x.y, y.y)
#sig.genes <- unique(clust4.genes$GeneID)

# "TGGT1_234440"
# sig.genes <- c( "TGGT1_200270")
# par(mfrow = c(1,1))
# for(i in 1:length(sig.genes)){
#   tmp <- extra.atac.spline.fits %>% dplyr::filter(GeneID == sig.genes[i])
#   plot(tmp$x, tmp$y, type = 'l', lwd = 2, col = 'red')
#   Sys.sleep(0.8)
# }


plot.atac.trend <- function(gene.id, gene.name,tc.logCPM){
  
  tmp <- tc.logCPM %>% dplyr::filter(GeneID == gene.id) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, color = 'red')) + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + 
    theme_bw(base_size = 14) +
    ylab('Expr') + xlab('Passage') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste("atac", paste(gene.name, gene.id, sep = ':'), sep = ":")) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
      axis.title.y = element_text(size=14, face="bold"),
      legend.position = 'none'
    )
  
  
  return(p)
  
}

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

extra.tc.logCPM.rna <- readRDS('../Input/Toxo_lab_adapt/RDS/extra_tc_rna_logCPM.rds')
extra.tc.logCPM.atac <- readRDS('../Input/Toxo_lab_adapt/RDS/extra_tc_atac_logCPM.rds')

matched.trending.rna.ata <- matched.rna.atac %>% filter(trending.x == trending.y) %>% 
  distinct(GeneID, .keep_all = T) %>%
  group_by(trending.x, cluster.x)  %>%
  summarise(genes = list(GeneID), total = n())



gene <- "TGGT1_299020" ; name <- "AP2III-4"
gene <- "TGGT1_320680" ; name <- " AP2IV-2"
gene <- "TGGT1_268860" ; name <- "ENO1"
gene <- "TGGT1_268850" ;name <- "ENO2"

# down regulated upon iKD of AP2X-5
gene <- "TGGT1_318610" ; name <- "AP2IV-3"
gene <- "TGGT1_208020" ; name <- "AP2Ib-1"
gene <- "TGGT1_306620" ; name <- "AP2IX-9"
gene <- "TGGT1_320680" ; name <- "AP2IV-2"

p1 <- plot.atac.trend(gene.id = gene, gene.name = name,tc.logCPM = extra.tc.logCPM.atac)
p2 <- plot.rna.trend(gene.id = gene,gene.name = name, tc.logCPM = extra.tc.logCPM.rna)


p <- do.call(grid.arrange, c(list(p1,p2), nrow=1))

ggsave(filename="../Output/Toxo_lab_adapt/figs/ENO2_rna_atac_trend.png", 
       plot=p,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

genes <- c("TGGT1_268860","TGGT1_268850",  "TGGT1_299020", "TGGT1_320680",
           "TGGT1_208020","TGGT1_306620", "TGGT1_318610" )

names <- c("ENO1",  "ENO2", "AP2III-4", " AP2IV-2",
           "AP2Ib-1", "AP2IX-9", "AP2IV-3")


p <- lapply(1:length(genes), function(i){
  
  pp <- plot.rna.trend(genes[i], names[i], extra.tc.logCPM.rna)

})


pdf("../Output/Toxo_lab_adapt/figs/AP2sENO1_ENO2_trend.pdf", height = 10, width = 8 ,onefile = TRUE)

do.call("grid.arrange", c(p, ncol = 2))  

dev.off()





pdf(file = "../Output/Toxo_lab_adapt/figs/ENO2_AP2s_trend.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

#specify to save plots in 2x2 grid
par(mfrow = c(4,2))

#save plots to PDF
for (i in 1:8) {   
  
  plot.rna.trend(gene.id = genes[i], gene.name = names[i], tc.logCPM = extra.tc.logCPM.rna)
  
}

dev.off()

mismatched.trending.rna.ata <- matched.rna.atac %>% filter(trending.x != trending.y) %>% 
  distinct(GeneID, .keep_all = T) %>%
  group_by(trending.x, cluster.x)  %>%
  summarise(genes = list(GeneID), total = n())

# gene <- "TGGT1_201250"
# p1 <- plot.atac.trend(gene.id = gene, tc.logCPM = extra.tc.logCPM.atac)
# p2 <- plot.rna.trend(gene.id = gene, tc.logCPM = extra.tc.logCPM.rna)
# 
# p2 <- do.call(grid.arrange, c(list(p1,p2), nrow=2))
# ggsave(filename="../Output/compScBdTgPb/figs/mismatch_trending_rna_atac.pdf",
#        plot=p2,
#        width = 8, height = 8,
#        units = "in", # other options are "in", "cm", "mm"
#        dpi = 300
# )

## 

#############
#############
#############

ENO2_targ <- read.xlsx("../Input/Toxo_lab_adapt/genes/ENO2_chip.xlsx")
# ENO2_targ$ID <- gsub(" ", "", ENO2_targ$ID, )
orth  <- read.xlsx("../Input/Toxo_lab_adapt/genes/convertIDs.xlsx")
ENO2_targ <- inner_join(ENO2_targ, orth, by = c("ID" = "TGME49ID") )

extra.rna.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_rna_dtw_trending_2_clusters_less_stringent.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_atac_dtw_trending_2_clusters.rds")

extra.rna.dtw.long.filt <- left_join(extra.rna.dtw.long.filt, gene.desc, by = "GeneID")
extra.atac.dtw.long.filt <- left_join(extra.atac.dtw.long.filt, gene.desc, by = "GeneID")

ENO2.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% ENO2_targ$TGGT1ID, ]
ENO2.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% ENO2_targ$TGGT1ID, ]

write.xlsx(ENO2.rna, "../Output/Toxo_lab_adapt/table/ENO2_targets_expr_rna.xlsx")
write.xlsx(ENO2.atac, "../Output/Toxo_lab_adapt/table/ENO2_targets_acc_atac.xlsx")

ENO2.rna.stat <- ENO2.rna %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
ENO2.atac.stat <- ENO2.atac %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
#
colnames(ENO2.rna.stat) <- c("Trend", "#genes")
# ENO2.rna.stat.gt <- gt(ENO2.rna.stat)
# gtsave(data = ENO2.rna.stat.gt, filename = "ENO2_rna_stat.png", path = "../Output/compScBdTgPb/figs/")

p <- ggplot(ENO2.rna, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Expr') + xlab('Passage') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_grid(trending~.) + 
  theme(
    strip.text.x = element_text(
      size = 14,  face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 14, face = "bold.italic"
    )
  ) +
  ggtitle("rna: ENO2 Target Genes") +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'black', hjust = 0.5),
    axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/ENO2_target_up_down_rna.png", 
       plot=p,
       width = 3, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

# TF

TF.tab <- read.xlsx("../Input/Toxo_lab_adapt/genes/TF_Info_Updated.xlsx")
extra.tc.logCPM.rna.tf <- left_join(extra.rna.dtw.long.filt, TF.tab, by = c("GeneID" = "GeneName"))

# ENO2 tarhets
ENO2.rna.info <-  inner_join(ENO2.rna, ENO2_targ, by = c("GeneID"  ="TGGT1ID"))
ENO2.rna.info <- ENO2.rna.info %>% distinct(GeneID, .keep_all = T) %>% transmute(GeneID = GeneID, trending, Description) %>% arrange(trending)

# info on targets of ENO2 that are up or down trending in our rna seq data 
write.xlsx(ENO2.rna.info, "../Output/compScBdTgPb/table/ENO2_rna_info.xlsx") 

ENO2.atac.info <-  inner_join(ENO2.atac, ENO2_targ, by = c("GeneID"  ="TGGT1ID"))
ENO2.atac.info <- ENO2.atac.info %>% distinct(GeneID, .keep_all = T) %>% transmute(GeneID = GeneID, trending, Description) %>% arrange(trending)

#info on targets of ENO2 that are up or down trending in our atac seq data 
write.xlsx(ENO2.atac.info, "../Output/compScBdTgPb/table/ENO2_atac_info.xlsx")
#tmp.gt <- gt(tmp)
#gtsave(data = tmp.gt, filename = "ENO2_targets_Up_Down_rna.png", path = "../Output/compScBdTgPb/figs/")

#################
#################
#################

extra.rna.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_rna_dtw_trending_2_clusters_less_stringent.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_atac_dtw_trending_2_clusters.rds")
gene.desc <- read.xlsx("../Input/Toxo_lab_adapt/genes/ProductDescription_GT1.xlsx")

hyperlopid <- read.xlsx("../Input/Toxo_lab_adapt/genes/hyperlopit_YR.xlsx")
hyperlopid <- hyperlopid %>%  dplyr::select(Accession ,tagm.map.allocation.pred)
toxo.orth <- read.xlsx("../Input/Toxo_lab_adapt/genes/convertIDs.xlsx")
toxo.orth <- inner_join(toxo.orth, gene.desc, by = c("TGGT1ID" = "GeneID"))
hyperlopid <- left_join(hyperlopid, toxo.orth, by = c("Accession" = "TGME49ID")) %>% na.omit()
hyperlopid <- hyperlopid %>% transmute(GeneID = TGGT1ID, HL.tagm.map.allocation.pred =  tagm.map.allocation.pred)

rna.nuclear.localiz <- read.xlsx("../Input/compScBdTgPb/stress/arabbidopsis_thailiana_orth_up_down_trending_rna_nuclus_localization.xlsx")
rna.nuclear.localiz <- rna.nuclear.localiz %>% dplyr::select(GeneID, AT_orth_GO_cellular_component)

localization.df <- left_join(hyperlopid, rna.nuclear.localiz, by = "GeneID")


