


UV <- read.table("../Input/compScBdTgPb/stress/STIFDB/UVB.txt", sep = "\t") 
UV.df <- UV %>% dplyr::transmute(GeneID = V2) %>% mutate(UV = "UV")

HE <- read.table("../Input/compScBdTgPb/stress/STIFDB/Heat.txt", sep = "\t") 
HE.df <- HE %>% dplyr::transmute(GeneID = V2) %>% mutate(HE = "HE")

OSM <- read.table("../Input/compScBdTgPb/stress/STIFDB/Osmotic.txt", sep = "\t") 
OSM.df <- OSM %>% dplyr::transmute(GeneID = V2) %>% mutate(OSM = "OSM")

IRON <- read.table("../Input/compScBdTgPb/stress/STIFDB/IRON.txt", sep = "\t") 
IR.df <- IRON %>% dplyr::transmute(GeneID = V2) %>% mutate(IR = "IR")

ALUM <- read.table("../Input/compScBdTgPb/stress/STIFDB/AL.txt", sep = "\t") 
AL.df <- ALUM %>% dplyr::transmute(GeneID = V2) %>% mutate(AL = "AL")


WOUND <- read.table("../Input/compScBdTgPb/stress/STIFDB/WOUND.txt", sep = "\t") 
WO.df <- WOUND %>% dplyr::transmute(GeneID = V2) %>% mutate(WO = "WO")


ABA <- read.table("../Input/compScBdTgPb/stress/STIFDB/ABA.txt", sep = "\t") 
ABA.df <- ABA %>% dplyr::transmute(GeneID = V2) %>% mutate(AB = "AB")

COLD <- read.table("../Input/compScBdTgPb/stress/STIFDB/COLD.txt", sep = "\t") 
CO.df <- COLD %>% dplyr::transmute(GeneID = V2) %>% mutate(CO = "CO")


DROUGHT <- read.table("../Input/compScBdTgPb/stress/STIFDB/DROUGHT.txt", sep = "\t") 
DR.df <- DROUGHT %>% dplyr::transmute(GeneID = V2) %>% mutate(DR = "DR")

LIGHT <- read.table("../Input/compScBdTgPb/stress/STIFDB/LIGHT.txt", sep = "\t") 
LI.df <- LIGHT %>% dplyr::transmute(GeneID = V2) %>% mutate(LI = "LI")

NaCI <- read.table("../Input/compScBdTgPb/stress/STIFDB/NaCI.txt", sep = "\t") 
NaCI.df <- NaCI %>% dplyr::transmute(GeneID = V2) %>% mutate(NaCI = "NaCI")


OXIDATIV <- read.table("../Input/compScBdTgPb/stress/STIFDB/OXIDATIVE.txt", sep = "\t") 
OX.df <- OXIDATIV %>% dplyr::transmute(GeneID = V2) %>% mutate(OX = "OX")

DEHYD <- read.table("../Input/compScBdTgPb/stress/STIFDB/DEHYDRATION.txt", sep = "\t") 
DH.df <- DEHYD %>% dplyr::transmute(GeneID = V2) %>% mutate(DH = "DH")

CDS <- read.table("../Input/compScBdTgPb/stress/STIFDB/cold_drought_salt.txt", sep = "\t") 
CDS.df <- CDS %>% dplyr::transmute(GeneID = V2) %>% mutate(CDS = "CDS")



plant.stress.genes <- list(UV.df, HE.df, OSM.df, IR.df, AL.df, WO.df, ABA.df, 
                           CO.df, DR.df, LI.df, NaCI.df, OX.df, DH.df, CDS.df)


tmp.list <- lapply(1:length(plant.stress.genes), function(i){
  
  df <- plant.stress.genes[[i]]
  colnames(df)[2] <- "stress"
  return(df)
})

tmp.df <- do.call("rbind", tmp.list)

#plant.stress.genes.df <- purrr::reduce(plant.stress.genes, dplyr::left_join, by = 'GeneID' )


write.xlsx(tmp.df, "../Input/compScBdTgPb/stress//arabidopsis_stress_genes_STIFDB_YR.xlsx")



#test.df <- read.table("~/Desktop/tmp.txt", sep = "\t")


library(tidyverse)
library(rvest)


plant.urls <- read.table("../Input/compScBdTgPb/stress/STIFDB/URL_Stress_responsive_genes_predicted_TF_arabidopsis.txt", header = F)
str <- strsplit(plant.urls$V1, "=")
TF.str <- unlist(lapply(str, "[" , 3))


TF.Genes.list <- lapply(1:length(TF.str), function(i){
  
  url <- plant.urls[i, ]
  tmp.all <- url %>% read_html() %>% html_nodes("table") %>% html_table(fill = T)
  tf.genes.df <- tmp.all[[5]]
  
  names(tf.genes.df) <- tf.genes.df[1,]
  tf.genes.df <- tf.genes.df[-1,]
  
  return(tf.genes.df)
  
})

names(TF.Genes.list) <- TF.str

saveRDS(TF.Genes.list, "../Input/compScBdTgPb/RData/arabiidopsis_tf_genes_list.RData")
TF.Genes.list <- readRDS("../Input/compScBdTgPb/RData/arabiidopsis_tf_genes_list.RData")

Cis.urls <- read.table("../Input/compScBdTgPb/stress/STIFDB/URL_TF_family_Cis_element_arabidopsis.txt", header = F)
str <- strsplit(Cis.urls$V1, "=")
Cis.str <- unlist(lapply(str, "[" , 3))


Cis.list <- lapply(1:length(Cis.str), function(i){
  
  url <- Cis.urls[i, ]
  tmp.all <- url %>% read_html() %>% html_nodes("table") %>% html_table(fill = T)
  cis.df <- tmp.all[[5]]
  cis.df <- cis.df[6,-1]
   
  names(cis.df) <- "cis.element"
  
  
  return(cis.df)
  
})

names(Cis.list) <- Cis.str

cis.element.tf <- do.call('rbind', Cis.list) %>% rownames_to_column(var = "tf.family")
write.xlsx(cis.element.tf, "../Input/compScBdTgPb/stress/arabidopsis_tf_cis_element_STIFDB.xlsx")


arab.tf.genes <- readRDS("../Input/compScBdTgPb/RData/arabiidopsis_tf_genes_list.RData")
arab.tf.genes.all <- do.call("rbind", arab.tf.genes) %>% rownames_to_column(var = "tf.family")
arab.tf.genes.all$tf.family <- gsub("\\..*", "",arab.tf.genes.all$tf.family)
colnames(arab.tf.genes.all) <- gsub("Gene ID", "GeneID", colnames(arab.tf.genes.all))

cis.element.tf <- read.xlsx("../Input/compScBdTgPb/stress/STIFDB/arabidopsis_tf_cis_element_STIFDB.xlsx")
arab.stress.genes <- read.xlsx("../Input/compScBdTgPb/stress/STIFDB/arabidopsis_stress_genes_STIFDB_YR.xlsx")


arab.genes.strs.tf <- left_join(arab.stress.genes, arab.tf.genes.all, by = "GeneID") 
arab.genes.strs.tf.cis <- left_join(arab.genes.strs.tf, cis.element.tf, by = "tf.family")

saveRDS(arab.genes.strs.tf.cis, "../Input/compScBdTgPb/RData/arabidopsis_stress_genes_STIFDB.RData")


# arap <- read.xlsx("../Input/compScBdTgPb/stress/stress_related_genes_arabidopsis.xlsx", colNames  =F)
# write.table(arap, "../Input/compScBdTgPb/stress/PlantGenome/stress_related_genes_arabidopsis_IDs.txt", row.names = F, 
#             quote = FALSE)


yeast <- read.xlsx("../Input/compScBdTgPb/stress/stress_related_genes_saccharomyces.xlsx", colNames  =T)
yeast <- data.frame(yeast$Gene.Systematic.Name)
write.table(yeast, "../Input/compScBdTgPb/stress/YeastGenome/stress_related_genes_saccharomyces_IDs.txt", row.names = F, 
            col.names = F, 
            quote = FALSE)

#grep -A 1 -wFf list.txt sequences.fas > newfile2.fas


STIFDB <- read.table("~/Desktop/stress_related_STIFDBV2.txt", sep = "\t") 
NaCI.df <- NaCI %>% dplyr::transmute(GeneID = V2) %>% mutate(NaCI = "NaCI")
