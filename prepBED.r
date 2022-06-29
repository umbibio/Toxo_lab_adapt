


extra.rna.dtw.long <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_2_clusters_less_stringent.rds")
count.tab <- read.xlsx("../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")
Peaks_Genes_assigned <- count.tab %>% dplyr::select(gene_name, peak_location)

p <- ggplot(extra.rna.dtw.long, aes(x = x, y = y, group = GeneID)) +
  
  geom_line(aes(x = x, y = y, color = factor(cluster)), alpha = 0.4) +
  
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + theme_bw() +
  
  facet_wrap(.~cluster)


plot(p)


extra.rna.dtw.genes <- extra.rna.dtw.long %>% dplyr::select(GeneID, cluster) %>% distinct()
extra.rna.dtw.genes.peaks <- inner_join(extra.rna.dtw.genes, Peaks_Genes_assigned, by = c("GeneID" = "gene_name")) 
peak_regions_str <- strsplit(as.character(extra.rna.dtw.genes.peaks$peak_location), split="[:|-]+")

extra.rna.dtw.genes.peaks$V1 <- sapply(peak_regions_str, '[', 1)
extra.rna.dtw.genes.peaks$V2 <- sapply(peak_regions_str, '[', 2)
extra.rna.dtw.genes.peaks$V3 <- sapply(peak_regions_str, '[', 3)

cluster.list <- extra.rna.dtw.genes.peaks %>% split(f = extra.rna.dtw.genes.peaks$cluster)
names(cluster.list) <-  paste("cluster", names(cluster.list), sep = "_")
# 
# ENO1 manual modification of the atac peak
cluster.list$cluster_1$V2[cluster.list$cluster_1$GeneID == "TGGT1_268860"] <- "6246801" # end of intron
cluster.list$cluster_1$V3[cluster.list$cluster_1$GeneID == "TGGT1_268860"]  <- 6247335 + 1000 # 1000 upstream of gene
#

cluster.list$cluster_2$V2[cluster.list$cluster_2$GeneID == "TGGT1_268850"] 
cluster.list$cluster_2$V3[cluster.list$cluster_2$GeneID == "TGGT1_268850"]  <- 6250901 + 2000 # 1000 upstream of gene
out.dir <- "../Input/compScBdTgPb/BulkATACToxoPlasma/motif/bed_up_down_2_cluster_V6/"

cluster.bed.list <- lapply(1:length(cluster.list), function(i){
  
  name <- paste(names(cluster.list)[i], "bed", sep = ".")
  tmp <- cluster.list[[i]]
  cluster.bed <- tmp %>% dplyr::select(V1:V3, GeneID)
  write.table(cluster.bed, paste(out.dir, name, sep = ""), sep = "\t", quote = F, row.names = F, col.names = F)
  
  #return(tmp)
  return(cluster.bed)
})

names(cluster.bed.list) <- names(cluster.list)



##############################################################


matched.rna.atac <- readRDS("../Input/compScBdTgPb/RData/extra_rna_atac_stw_trending_2_clusters.rds")
count.tab <- read.xlsx("../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")
Peaks_Genes_assigned <- count.tab %>% dplyr::select(gene_name, peak_location)



matched.rna.atac.long <- matched.rna.atac %>% filter(trending.x == trending.y) %>% 
  distinct(GeneID, .keep_all = T) %>% transmute(GeneID = GeneID, trend = trending.x)

matched.rna.atac.gene.peak <- inner_join(matched.rna.atac.long, Peaks_Genes_assigned, by = c("GeneID" = "gene_name")) 
peak_regions_str <- strsplit(as.character(matched.rna.atac.gene.peak$peak_location), split="[:|-]+")

matched.rna.atac.gene.peak$V1 <- sapply(peak_regions_str, '[', 1)
matched.rna.atac.gene.peak$V2 <- sapply(peak_regions_str, '[', 2)
matched.rna.atac.gene.peak$V3 <- sapply(peak_regions_str, '[', 3)

trend.list <- matched.rna.atac.gene.peak %>% split(f = matched.rna.atac.gene.peak$trend)

out.dir <- "../Input/compScBdTgPb/BulkATACToxoPlasma/motif/bed_matched_Up_dow_rna_atac_2cluster/"
cluster.bed.list <- lapply(1:length(trend.list), function(i){
  
  name <- paste(names(trend.list)[i], "bed", sep = ".")
  tmp <- trend.list[[i]]
  cluster.bed <- tmp %>% dplyr::select(V1:V3)
  
  write.table(cluster.bed, paste(out.dir, name, sep = ""), sep = "\t", quote = F, row.names = F, col.names = F)
  
  return(cluster.bed)
})

names(cluster.bed.list) <- names(trend.list)



