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

#############
#############
#############

ENO2_targ <- read.xlsx("../Input/Toxo_lab_adapt/genes/ENO2_chip.xlsx")
# ENO2_targ$ID <- gsub(" ", "", ENO2_targ$ID, )
orth  <- read.xlsx("../Input/Toxo_lab_adapt/genes/convertIDs.xlsx")
ENO2_targ <- inner_join(ENO2_targ, orth, by = c("ID" = "TGME49ID") )

extra.rna.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_rna_dtw_trending_2_clusters_less_stringent.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/Toxo_lab_adapt/RDS/extra_atac_dtw_trending_2_clusters.rds")


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

### ran up to here 

## arab peptide orthologs in GT1
arab.pep.orth <- readRDS( "../Input/Toxo_lab_adapt/RData/rec_GT1_vs_Arab_pep.RData")
str <- strsplit(arab.pep.orth$query_id, "-")
arab.pep.orth$query_id <- unlist(lapply(str, "[", 1))

arab.pep.orth.desc <- left_join(arab.pep.orth, gene.desc, by = c("query_id" = "GeneID"))
arab.pep.orth.desc <- arab.pep.orth.desc %>% dplyr::select(query_id, subject_id, ProductDescription)
write.xlsx(arab.pep.orth.desc, "../Output/compScBdTgPb/table/arab_toxo_pep_orth_desc_215.xlsx")

arab.pep.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% arab.pep.orth.desc$query_id, ]

arab.pep.rna.orth <- inner_join(extra.rna.dtw.long.filt, arab.pep.orth.desc, by = c("GeneID" = "query_id")) %>% 
  distinct(GeneID, .keep_all = T) %>%
  transmute(GeneID = GeneID, arab_orth = subject_id,trending = trending, ProductDescription = ProductDescription) %>%
  mutate(data = "rna") 


arab.pep.rna.orth <- left_join(arab.pep.rna.orth, localization.df, by  = "GeneID")
write.xlsx(arab.pep.rna.orth, "../Output/compScBdTgPb/table/arabbidopsis_thailiana_orth_up_down_trending_rna.xlsx")

# arab.pep.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% arab.pep.orth.desc$query_id, ]
# 
# arab.pep.atac.orth <- inner_join(extra.atac.dtw.long.filt, arab.pep.orth.desc, by= c("GeneID" = "query_id")) %>%
#   distinct(GeneID, .keep_all = T) %>% 
#   transmute(GeneID = GeneID, arab_orth = subject_id,trending = trending, ProductDescription = ProductDescription) %>%
#   mutate(orth_to_arabidopsis = "yes", data = "atac") 
# 
# arab.pep.atac.orth <- left_join(arab.pep.atac.orth, localization.df, by = "GeneID" )
# write.xlsx(arab.pep.atac.orth, "../Output/compScBdTgPb/table/arabbidopsis_thailiana_orth_up_down_trending_atac.xlsx")
# 

arab.pep.rna.stat <- arab.pep.rna %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
#arab.pep.atac.stat <- arab.pep.atac %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
#
colnames(arab.pep.rna.stat) <- c("Trend.rna", "#genes")
colnames(arab.pep.atac.stat) <- c("Trend.atac", "#genes")

p <- ggplot(arab.pep.rna, aes(x = x, y = y, group = GeneID)) + 
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
  ggtitle("rna: arabidopsis thaliana stress genes (peptides orthologs)") +
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

ggsave(filename="../Output/compScBdTgPb/figs/stress_related_arab_pep_orth_rna_up_down_trend_genes.png", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

# # matched rna-atac trending (very a few)
# rna.atac <- inner_join(arab.pep.rna, arab.pep.atac, by = "GeneID") %>% 
#   distinct(GeneID, .keep_all =T) %>% transmute(GeneID = GeneID, trending.rna  = trending.x , trending.atac = trending.y)
# 
# rna.atac <- tibble(rna.atac)
# rna.atac.agreement <- table(rna.atac$trending.rna, rna.atac$trending.atac)
# rna.atac.agreement


 ############ Add STIFDB info 

cis.element.tf <- read.xlsx("../Input/compScBdTgPb/stress/STIFDB/arabidopsis_tf_cis_element_STIFDB.xlsx")
cis.elements <- unique(cis.element.tf$cis.element)
cis.elements <- gsub("\\/|\\,", "|", gsub("\\(", "\\[",gsub("\\)", "]",cis.elements)))
cis.elements

# orthologs
arab.pep.orth <- readRDS( "../Input/compScBdTgPb/RData/rec_GT1_vs_Arab_pep.RData")
str <- strsplit(arab.pep.orth$query_id, "-")
arab.pep.orth$query_id <- unlist(lapply(str, "[", 1))


# 717 stress related genes in arabidopsis (all types of stress)
arab.stress.genes <- read.xlsx("../Input/compScBdTgPb/stress/stress_related_genes_arabidopsis.xlsx", colNames  = F)
colnames(arab.stress.genes) <- "GeneID"

# 3110 stress related TFs from STIFDB
arab.genes.strs.tf.cis <- readRDS("../Input/compScBdTgPb/RData/arabidopsis_stress_genes_STIFDB.RData")
arab.genes.strs.tf.cis <- arab.genes.strs.tf.cis %>% na.omit()

# stress related  TFs
tf.stress.genes <- inner_join(arab.stress.genes, arab.genes.strs.tf.cis, by = "GeneID" )
length(unique(tf.stress.genes$GeneID)) # stress related  TF
arab.pep.rna.orth <- arab.pep.rna.orth %>% dplyr::mutate(arab_gene_id = gsub("\\..*", "", arab.pep.rna.orth$arab_orth))

# 12 TFs in rna 
tf.rna.df <- inner_join(arab.pep.rna.orth, tf.stress.genes , by = c("arab_gene_id" = "GeneID" ))
write.xlsx(tf.rna.df, "../Output/compScBdTgPb/table/arabbidopsis_thailiana_orth_up_down_trending_rna_tf_DTIFDB.xlsx")


cis.freq.up.down <- readRDS(  "../Input/compScBdTgPb/RData/cis_freq_under_atac_peaks_up_down_rna_genes.RData")

rna.df <- inner_join(tf.rna.df, cis.freq.up.down, by = "cis.element")  
colnames(rna.df) <- gsub("GeneID", "source.tf.id" , gsub("gene_name", "target.gene.id", colnames(rna.df) ))



rna.tf.targ <- rna.df %>% filter(cis.freq > 0) %>% 
  group_by(source.tf.id, stress, cis.element) %>% mutate(total_num_tf_targets = length(unique(target.gene.id)))


tf.info <- read.xlsx("../Input/compScBdTgPb/genes/TF_Info_Updated.xlsx")
tf.info <- tf.info %>% dplyr::select(GeneName, Ap2Name)

shared.brady.extra <- read.xlsx("../Input/compScBdTgPb/genes/shared_brady_extra_markers.xlsx")

shared.brady.ord <- shared.brady.extra %>% 
  filter(avg_log2FC.Brady >= 1 & p_val_adj.Brady < 0.01) %>% dplyr::select(GeneID) %>%
  mutate(is.brady.sig = "yes") %>% arrange("avg_log2FC.Brady") %>% slice(1:10)

shared.extra.ord <- shared.brady.extra %>% 
  filter(avg_log2FC.Extra >= 1 & p_val_adj.Extra < 0.01) %>% dplyr::select(GeneID) %>%
  mutate(is.extra.sig = "yes") %>% arrange("avg_log2FC.Extra") %>% slice(1:10) 

top.shared.markers <- unique(c(shared.brady.ord$GeneID, shared.extra.ord$GeneID))
intersect(tf.info$GeneName, top.shared.markers)

rna.tf.targ <- rna.tf.targ %>% mutate(is.tf = ifelse(target.gene.id %in% tf.info$GeneName , "yes", "no"))
rna.tf.targ <- rna.tf.targ %>% mutate(is.brady.marker = ifelse(target.gene.id %in% shared.brady.ord$GeneID, "yes", "no"))
rna.tf.targ <- rna.tf.targ %>% mutate(is.extra.marker = ifelse(target.gene.id %in% shared.extra.ord$GeneID, "yes", "no"))

colnames(rna.tf.targ) <- gsub("_", ".", colnames(rna.tf.targ))
rna.tf.targ <- rna.tf.targ %>% dplyr::select(-tf.family.y)

colnames(rna.tf.targ) <- gsub('trending.x', 'source.tf.trend', 
                              gsub('ProductDescription.x', 'source.ProductDescription', 
                                   gsub('HL.tagm.map.allocation.pred' , 'HyperLopit.tagm.map.allocation.pred', 
                                        gsub('AT.orth.GO.cellular.component', 'arab.orth.GO.cellular.component', 
                                                  gsub('^Description', 'arab.Description', 
                                                            gsub('tf.family.y', 'tf.family', 
                                                                 gsub('ProductDescription.y', 'target.ProductDescription',
                                                                      gsub('trending.y', 'target.trend', 
                                                                           gsub("tf.family.x", "tf.family", 
                                                                                gsub('peak.location',  'atac.peak.region',colnames(rna.tf.targ)))))))))))


write.xlsx(rna.tf.targ, "../Output/compScBdTgPb/table/tf.source.target.xlsx")


NLS <- read.xlsx("~/Downloads/TFs_NLS_predictedresults.xlsx")
NLS <- NLS %>% na.omit() 
NLS <- NLS[,c(1,3)]

rna.tf.targ.NLS <- left_join(rna.tf.targ, NLS, by = c("source.tf.id"= "gene_ID"))  


tfs.up <- rna.tf.targ.NLS %>% ungroup() %>%
  dplyr::select(source.tf.id, arab.orth,source.ProductDescription ,Predicted_results,source.tf.trend,cis.element, stress ,total.num.tf.targets) %>%
  filter(source.tf.trend == "up") %>% distinct()



write.xlsx(tfs.up, "../Output/compScBdTgPb/table/tfs_up_info_slide.xlsx")
tfs.up.gt <- gt(tfs.up)
gtsave(data = tfs.up.gt, filename = "tfs.up.png", path = "../Output/compScBdTgPb/figs/")

genes <- unique(tfs.up$source.tf.id)


tmp <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "TGTCTC") %>% 
  distinct()

tmp1 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "GCCGCC") %>% 
  distinct()

tmp2 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "CGCGTG") %>% 
  distinct()

tmp3 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "CATGTG") %>% 
  distinct()

tmp4 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "CACGTG") %>% 
  distinct()

tmp5 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "TAACTG") %>% 
  distinct()

tmp6 <- rna.tf.targ %>% ungroup() %>%
  dplyr::select(source.tf.id, source.ProductDescription, tf.family, stress, cis.element) %>% filter(cis.element == "CCACGTGG") %>% 
  distinct()




genes <- unique(tmp$source.tf.id) 

#shared.df <- full_join(shared.brady.ord,  shared.extra.ord, by = "GeneID") %>% distinct(GeneID, .keep_all  = T)


extra.ovlp.target  <- shared.extra.ord$GeneID[which(unique(shared.extra.ord$GeneID) %in%  unique(rna.tf.targ$target.gene.id))]
brady.ovrp.target <- shared.brady.ord$GeneID[which(unique(shared.brady.ord$GeneID) %in%  unique(rna.tf.targ$target.gene.id))]
tf.ovlp.target <- tf.info$GeneName[which(unique(tf.info$GeneName) %in%  unique(rna.tf.targ$target.gene.id))]

Reduce(intersect, list(extra.ovlp.target, brady.ovrp.target, tf.ovlp.target))


rna.tf.targ.info <- left_join(rna.tf.targ,  tf.info, by = c("source.tf.id" = "GeneName"))
rna.tf.targ.info <- left_join(rna.tf.targ.info, shared.df, by = c("source.tf.id" = "GeneID"))


unique(rna.tf.targ$total_num_tf_targets)

unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 1220])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 879])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 843])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 759])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 661])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 422])
unique(rna.tf.targ.info$cis.element[rna.tf.targ.info$total_num_tf_targets == 58])



write.xlsx(rna.tf.targ, "../Output/compScBdTgPb/table/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_rna.xlsx")
write.xlsx(atac.tf.targ, "../Output/compScBdTgPb/table/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_atac.xlsx")

saveRDS(rna.tf.targ, "../Input/compScBdTgPb/RData/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_rna.RData")
saveRDS(atac.tf.targ, "../Input/compScBdTgPb/RData/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_atac.RData")


rna.atac.tf.targ <- rbind(rna.tf.targ, atac.tf.targ)
write.xlsx(rna.atac.tf.targ, "../Output/compScBdTgPb/table/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_rna_atac.xlsx")
saveRDS(rna.atac.tf.targ, "../Input/compScBdTgPb/RData/stress_related_tfs_in_toxo_orth_to_arabidopsis_targets_cis_element_freq_rna_atac.RData")


gene.id <- "TGGT1_268850"
tmp <- rna.atac.tf.targ[which(rna.atac.tf.targ$gene_name == gene.id),]

rna.tf.genes <- unique(rna.tf.targ$GeneID)
atac.tf.genes <- unique(atac.tf.targ$GeneID)

extra.tc.logCPM.rna <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')
extra.tc.logCPM.atac <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')

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


BFD1 <- "TGGT1_200385"
BFD2 <- "TGGT1_311100"
ENO1 <- "TGGT1_268860"
ENO2 <- "TGGT1_268850"

plot.rna.trend("TGGT1_268850", "BFD2", extra.tc.logCPM.rna)


# ## summarized version 
# rna.cis.element.scanned <- inner_join(arab.pep.rna.orth, cis.freq.up.down , by = c("GeneID" = "gene_name" ))
# write.xlsx(rna.cis.element.scanned, "../Output/compScBdTgPb/table/rna_up_down_cis_element_scanned.xlsx")
# 
# atac.cis.element.scanned <- inner_join(arab.pep.atac.orth, cis.freq.up.down , by = c("GeneID" = "gene_name" ))
# write.xlsx(atac.cis.element.scanned, "../Output/compScBdTgPb/table/atac_up_down_cis_element_scanned.xlsx")

### yeast 

yeast.orth <- read.xlsx("../Input/compScBdTgPb/stress/orthologs/GT1_Yeast_pep_orth.xlsx")
gene.desc <- read.xlsx("../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx")
yeast.orth.desc <- inner_join(gene.desc, yeast.orth,by = c("GeneID" = "GT1"))


yeast <- read.xlsx("../Input/compScBdTgPb/stress/stress_related_genes_saccharomyces.xlsx", colNames  =T)
colnames(yeast)[1] <- "yeast"

yeast.stress.df <- inner_join(yeast.orth.desc, yeast, by = "yeast")


extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")

yeast.pep.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% yeast.stress.df$GeneID, ]

yeast.pep.rna.orth <- inner_join(extra.rna.dtw.long.filt, yeast.stress.df, by = "GeneID") %>% 
  distinct(GeneID, .keep_all = T) %>% dplyr::select(!c(x, y, cluster, slope, R2))

write.xlsx(yeast.pep.rna.orth, "../Output/compScBdTgPb/table/yeast_orth_up_down_trending_rna.xlsx")

yeast.pep.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% yeast.stress.df$GeneID, ]

yeast.pep.atac.orth <- inner_join(extra.atac.dtw.long.filt, yeast.stress.df, by= "GeneID" ) %>%
  distinct(GeneID, .keep_all = T) %>% dplyr::select(!c(x, y, cluster, slope, R2))

write.xlsx(yeast.pep.atac.orth, "../Output/compScBdTgPb/table/yeast_orth_up_down_trending_atac.xlsx")


## Ecoli 
Ecoli <- read.xlsx("../Input/compScBdTgPb/stress/stress_related_genes_Ecoli.xlsx")
Ecoli <- Ecoli %>% dplyr::select(-EB.number)

Ecoli.orth <- read.xlsx("../Input/compScBdTgPb/stress/orthologs/GT1_Ecoli_pep_orth.xlsx")
str <- strsplit(Ecoli.orth$Ecoli, "\\|")
Ecoli.orth$Ecoli <- unlist(lapply(str, "[", 3))
gene.desc <- read.xlsx("../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx")
Ecoli.orth.desc <- inner_join(gene.desc, Ecoli.orth, by = c("GeneID" = "GT1"))
Ecoli.stress.df <- inner_join(Ecoli.orth.desc, Ecoli, by= c("Ecoli" = "Blattner.number"))



extra.rna.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters.rds")
extra.atac.dtw.long.filt <- readRDS("../Input/compScBdTgPb/RData/extra_atac_dtw_trending_2_clusters.rds")

Ecoli.pep.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% Ecoli.stress.df$GeneID, ]

Ecoli.pep.rna.orth <- inner_join(extra.rna.dtw.long.filt, Ecoli.stress.df, by = "GeneID") %>% 
  distinct(GeneID, .keep_all = T) %>% dplyr::select(!c(x, y, cluster, slope, R2))

write.xlsx(Ecoli.pep.rna.orth, "../Output/compScBdTgPb/table/Ecoli_orth_up_down_trending_rna.xlsx")

Ecoli.pep.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% Ecoli.stress.df$GeneID, ]

Ecoli.pep.atac.orth <- inner_join(extra.atac.dtw.long.filt, Ecoli.stress.df, by= "GeneID" ) %>%
  distinct(GeneID, .keep_all = T) %>%  dplyr::select(!c(x, y, cluster, slope, R2))

write.xlsx(yeast.pep.atac.orth, "../Output/compScBdTgPb/table/Ecoli_orth_up_down_trending_atac.xlsx")


# CDS 
# at.cds.orth <- readRDS("../Input/compScBdTgPb/RData/rec_GT1_vs_Arab_cds.RData")
# str <- strsplit(at.cds.orth$query_id, "-")
# at.cds.orth$query_id <- unlist(lapply(str, "[", 1))
# 
# 
# at.cds.rna <- extra.rna.dtw.long.filt[extra.rna.dtw.long.filt$GeneID %in% at.cds.orth$query_id, ]
# at.cds.atac <- extra.atac.dtw.long.filt[extra.atac.dtw.long.filt$GeneID %in% at.cds.orth$query_id, ]
# 
# at.cds.rna.stat <- at.cds.rna %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
# at.cds.atac.stat <- at.cds.atac %>% distinct(GeneID, .keep_all = T) %>% group_by(trending) %>% summarise(total = n())
# #
# colnames(at.cds.rna.stat) <- c("Trend", "#genes")
# colnames(at.cds.atac.stat) <- c("Trend", "#genes")
# 
# p <- ggplot(at.cds.rna, aes(x = x, y = y, group = GeneID)) + 
#   geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
#   geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
#   theme_bw(base_size = 14) +
#   ylab('Expr') + xlab('Passage') +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
#   theme(strip.background = element_rect(colour="black", fill="white",
#                                         size=0.5, linetype="solid")) +
#   facet_grid(trending~.) + 
#   theme(
#     strip.text.x = element_text(
#       size = 14,  face = "bold.italic"
#     ),
#     strip.text.y = element_text(
#       size = 14, face = "bold.italic"
#     )
#   ) +
#   ggtitle("rna: arabidopsis thaliana stress genes (CDS orthologs)") +
#   theme(
#     plot.title = element_text(size=14, face = "bold.italic", color = 'black', hjust = 0.5),
#     axis.title.x = element_text(size=14, face="bold", hjust = 0.5),
#     axis.title.y = element_text(size=14, face="bold")
#   ) + 
#   theme(legend.position = "None",
#         #legend.position = c(0.88, 0.17),
#         legend.title = element_text(colour="black", size=10, 
#                                     face="bold"),
#         legend.text = element_text(colour="black", size=10, 
#                                    face="bold")) + 
#   guides(colour = guide_legend(override.aes = list(size=2)))
# 
# plot(p)
# 
# ggsave(filename="../Output/compScBdTgPb/figs/AT_cds_rna_up_down_trend.png", 
#        plot=p,
#        width = 6, height = 6, 
#        units = "in", # other options are "in", "cm", "mm" 
#        dpi = 300
# )
# 
# 
# 
# 
