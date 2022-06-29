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



extra.tc.logCPM.atac <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')
genes <- unique(extra.tc.logCPM.atac$GeneID)


extra.atac.spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- extra.tc.logCPM.atac %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  extra.atac.sp <- smooth.spline(tmp$x, tmp$y)
  extra.atac.sp <- predict(extra.atac.sp, seq(11, 210,length.out = 10)) 
  mu <- data.frame(x = extra.atac.sp$x, y = extra.atac.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  fit <- lm(mu$y~mu$x)
  res <- summary(fit)
  mu$adj.r.squared <- res$adj.r.squared
  mu$trend.pval <- tmp$trend.pval[1]
  mu$trend.fdr <- tmp$trend.fdr[1]
  
  return(mu)
}, mc.cores = num.cores)


extra.atac.spline.fits <- bind_rows(extra.atac.spline.fits)

extra.atac.spline.fits.filt <- extra.atac.spline.fits %>% dplyr::filter(adj.r.squared > 0.5 & trend.fdr < 0.35) ## 1791

extra.atac.dtw.wide <- extra.atac.spline.fits.filt %>% 
  pivot_wider(-c(adj.r.squared, trend.fdr, trend.pval), names_from = 'GeneID', values_from = 'y') %>%
  mutate_at(vars(matches('TGGT1')), scale) %>%
  as.data.frame()

## Generate the clusters
num.clust <- 2L

#extra.atac.hc_dtws <- dtwClustCurves(extra.atac.dtw.wide[2:ncol(extra.atac.dtw.wide)], nclust = num.clust)

#saveRDS(extra.atac.hc_dtws, '../Input/compScBdTgPb/LabAdaptationRNA/extra.atac.hc_dtws_all.rds')

## Run from here
extra.atac.hc_dtws <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra.atac.hc_dtws_all.rds')

plot(extra.atac.hc_dtws, type = 'sc')
plot(extra.atac.hc_dtws, type = "series", clus = 2L)
plot(extra.atac.hc_dtws, type = "centroids", clus = 2L)


## plot a few curves
## Plot trends
# 
# sig.genes <- c( "TGGT1_200270", "TGGT1_234440")
# 
# par(mfrow = c(1,2))
# for(i in 1:length(sig.genes)){
#   tmp <- extra.atac.spline.fits %>% dplyr::filter(GeneID == sig.genes[i])
#   plot(tmp$x, tmp$y, type = 'l', lwd = 2, col = 'red')
#   Sys.sleep(0.8)
# }



atac.clust.info <- data.frame(GeneID = colnames(extra.atac.dtw.wide)[2:ncol(extra.atac.dtw.wide)], 
                             cluster = cutree(extra.atac.hc_dtws, k = 2))

gene.groups <- atac.clust.info %>% group_by(cluster) %>% summarise(genes = list(GeneID))


extra.atac.dtw.long <- extra.atac.dtw.wide %>% pivot_longer(-x, names_to = 'GeneID', values_to = 'y')
extra.atac.dtw.long <- left_join(extra.atac.dtw.long, atac.clust.info, by = 'GeneID')

fits <- lapply(1:nrow(gene.groups), function(i){
  #tmp <- extra.atac.spline.fits[extra.atac.spline.fits$GeneID %in% unlist(gene.groups$genes[i]), ]
  tmp <- extra.atac.dtw.long[extra.atac.dtw.long$GeneID %in% unlist(gene.groups$genes[i]), ]
  fit <- lm(tmp$y~tmp$x)
  summary(fit)
})

atac.r2 <- data.frame(R2 = unlist(lapply(fits, function(x) x$adj.r.squared)), 
                     slope = unlist(lapply(fits, function(x) x$coefficients[2])),
                     cluster = 1:nrow(gene.groups))

extra.atac.dtw.long <- left_join(extra.atac.dtw.long, atac.r2, by = 'cluster')
extra.atac.dtw.long$is.trending <- ifelse(extra.atac.dtw.long$R2 > 0.4, 'yes', 'no')
extra.atac.dtw.long$trending <- ifelse(extra.atac.dtw.long$slope > 0, 'up', 'down')
extra.atac.dtw.long.filt <- extra.atac.dtw.long %>% dplyr::filter(is.trending == 'yes')

saveRDS(extra.atac.dtw.long.filt, "../Input/compScBdTgPb/RData/extra.atac.dtw.long.filt.RData")

p <- ggplot(extra.atac.dtw.long.filt, aes(x = x, y = y, group = GeneID)) +
  geom_line(aes(x = x, y = y, color = factor(cluster)), alpha = 0.4) +
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + theme_bw() +
  facet_wrap(.~cluster)+
  ggtitle("ATAc-Seq Trending Genes") +
  theme(plot.title = element_text(face = "italic", hjust = 0.5))

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_trending_atac_2_clusters.pdf", 
       plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


saveRDS(extra.atac.dtw.long.filt, "../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")



########### STOP here ##########

## after this point should be deleted when uploading on github

 
## match RNA with ATAC

extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")

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

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_rna_up_down_trending.png", 
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

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_atac_up_down_trending.pdf", 
       plot=p,
       width = 3, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)





rna.atac <- inner_join(extra.rna.dtw.long.filt, extra.atac.dtw.long.filt , by = "GeneID") %>%
  distinct(GeneID, .keep_all =T) %>% transmute(GeneID = GeneID, trending.rna  = trending.x , trending.atac = trending.y)
rna.atac <- tibble(rna.atac)
rna.atac.agreement <- table(rna.atac$trending.rna, rna.atac$trending.atac)


dim(extra.rna.dtw.long.filt)
dim(extra.atac.dtw.long.filt)

extra.rna.summ <- extra.rna.dtw.long.filt %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% 
  summarise(genes = list(unique(GeneID)), total = n(), .groups ='drop')

extra.atac.summ <- extra.atac.dtw.long.filt %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% 
  summarise(genes = list(unique(GeneID)), total = n())

matched.rna.atac <- inner_join(extra.rna.dtw.long.filt, extra.atac.dtw.long.filt, by = "GeneID") %>% 
  dplyr::select(GeneID, everything())

saveRDS( matched.rna.atac,"../Input/compScBdTgPb/RData/matched_extra_rna_atac_dtw_trending_2_clusters.rds")

# ATAC expr for corresponding genes in RNA clusters
p <- ggplot(matched.rna.atac, aes(x = x.y, y = y.y, group = GeneID)) + 
  geom_line(aes(x = x.y, y = y.y, color = factor(cluster.x)), alpha = 0.4) + 
  #geom_smooth(aes(x = x.y, y = y.y, group = NA), method='lm', color = 'black') + 
  theme_bw() +
  facet_wrap(.~cluster.x)+
  ggtitle("matched genes ATAC-seq") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_matched_rna_atac_expr_2cluster.pdf", 
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


extra.tc.logCPM.rna <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')
extra.tc.logCPM.atac <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')

matched.trending.rna.ata <- matched.rna.atac %>% filter(trending.x == trending.y) %>% 
  distinct(GeneID, .keep_all = T) %>%
  group_by(trending.x, cluster.x)  %>%
  summarise(genes = list(GeneID), total = n())


gene <- "TGGT1_299020" ; name <- "AP2III-4"
gene <- "TGGT1_320680" ; name <- " AP2IV-2"
gene <- "TGGT1_268860" ; name <- "ENO1"
gene <- "TGGT1_268850" ;name <- "ENO2"

p1 <- plot.atac.trend(gene.id = gene, gene.name = name,tc.logCPM = extra.tc.logCPM.atac)
p2 <- plot.rna.trend(gene.id = gene,gene.name = name, tc.logCPM = extra.tc.logCPM.rna)

p <- do.call(grid.arrange, c(list(p1,p2), nrow=2))

ggsave(filename="../Output/compScBdTgPb/figs/ENO2_rna_atac_trend.png", 
       plot=p,
       width = 4, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

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


###########################
ENO2_targ <- read.xlsx("../Input/compScBdTgPb/genes/ENO2_chip.xlsx")
ENO2_targ$ID <- gsub(" ", "", ENO2_targ$ID, )
orth  <- read.xlsx("../Input/compScBdTgPb/genes/convertIDs.xlsx")
ENO2_targ <- inner_join(ENO2_targ, orth, by = c("ID" = "TGME49ID") )

extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")


ENO2.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% ENO2_targ$TGGT1ID, ]
ENO2.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% ENO2_targ$TGGT1ID, ]

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

ENO2.rna.info <-  inner_join(ENO2.rna, ENO2_targ, by = c("GeneID"  ="TGGT1ID"))
ENO2.rna.info <- ENO2.rna.info %>% distinct(GeneID, .keep_all = T) %>% transmute(GeneID = GeneID, trending, Description) %>% arrange(trending)

write.xlsx(ENO2.rna.info, "../Output/compScBdTgPb/table/ENO2_rna_info.xlsx")

ENO2.atac.info <-  inner_join(ENO2.atac, ENO2_targ, by = c("GeneID"  ="TGGT1ID"))
ENO2.atac.info <- ENO2.atac.info %>% distinct(GeneID, .keep_all = T) %>% transmute(GeneID = GeneID, trending, Description) %>% arrange(trending)

write.xlsx(ENO2.atac.info, "../Output/compScBdTgPb/table/ENO2_atac_info.xlsx")
tmp.gt <- gt(tmp)
gtsave(data = tmp.gt, filename = "ENO2_targets_Up_Down_rna.png", path = "../Output/compScBdTgPb/figs/")


#################################################

## arap peptide orthologs in GT1
arap.orth <- readRDS( "../Input/compScBdTgPb/RData/rec_GT1_vs_Arab_pep.RData")
str <- strsplit(arap.orth$query_id, "-")
arap.orth$query_id <- str.ids <- unlist(lapply(str, "[", 1))


extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")


arap.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% arap.orth$query_id, ]
arap.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% arap.orth$query_id, ]

arap.rna.stat <- arap.rna %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
arap.atac.stat <- arap.atac %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
#
colnames(arap.rna.stat) <- c("Trend", "#genes")

p <- ggplot(arap.atac, aes(x = x, y = y, group = GeneID)) + 
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
  #ggtitle("rna: ENO2 Target Genes") +
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


arab.genes.strs.tf.cis <- readRDS( "../Input/compScBdTgPb/RData/arabidopsis_stress_genes_STIFDB.RData")
