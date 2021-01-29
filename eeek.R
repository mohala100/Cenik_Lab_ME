# Requirements for the ribo R to run.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ribor")

install.packages("devtools")
library("devtools")
install_github("ribosomeprofiling/ribor")
library(ribor)
library(tidyr)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(ggplot2)
library("ggpubr")
install.packages("lmodel2")
library(lmodel2)

# Creating the ribo R object
second.ribo <- Ribo("/Users/mohanadelchouemi/Desktop/Dr Cenik/RibosomeTranslation/all.ribo", rename = rename_default )
second.ribo 


# Getting the mRNA data

rna <- get_rnaseq(second.ribo,
           tidy = TRUE,
           region = c( "CDS"),
           experiment = "GSM1691202",
           compact = FALSE,
           alias = TRUE) %>% mutate(GENE_NAME3 = str_replace(transcript , "-.*",replacement = "")) %>% select(-c(transcript)) %>% select(count,GENE_NAME3) %>% rename(rna_count = count)

# Getting the Ribosomal profiling footprint data

rc <- get_region_counts(second.ribo,tidy = TRUE,
                        compact      = FALSE,
                        transcript  = FALSE,
                        alias       = TRUE,
                        region      = c("CDS"),
                        experiment ="GSM1691202") %>% mutate(GENE_NAME3 = str_replace(transcript , "-.*",replacement = "")) %>% select(-c(transcript)) %>% select(count,GENE_NAME3) %>%
  rename(rc_count = count)



# Downloading the proteomics data.


# https://www.biostars.org/p/199073/

# Information abour RMA http://strata.uga.edu/8370/lecturenotes/regression.html

data <- read_xlsx(path = "/Users/mohanadelchouemi/Desktop/Dr Cenik/RibosomeTranslation/mmc4.xlsx", sheet = 2) %>% mutate_at(vars(2:3),~as.numeric(.x ))
data[is.na(data)] <- 0
data$prot_count <- (data$`iBAQ 1 protein copies per cell` + data$`iBAQ 2 protein copies per cell`) / 2
prot <- data %>% select(`gene symbol`,prot_count) %>% rename(GENE_NAME3 = `gene symbol`)

# RNA-Seq vs Proteomics
rna_prot <- left_join(x = rna,y = prot, by = "GENE_NAME3") %>% na.omit() %>% .[.$rna_count > 2 & .$prot_count >2,]
plot(log2(rna_prot$rna_count), log2(rna_prot$prot_count), xlab = "RNA-Seq", ylab = "Proteomics", main = "log2 CDS Read Counts", pch = 19, cex = 0.5)
abline(lm(log2(rna_prot$prot_count) ~ log2(rna_prot$rna_count), data = rc_prot), col = "blue")
cor(log2(rna_prot$rna_count), log2(rna_prot$prot_count), method = c("pearson"))
# Ranged Major Axis Regression for RNA-SEQ and Proteomics
lmodel2(log2(prot_count) ~ log2(rna_count),  data = rna_prot,range.y = "relative",range.x = "relative",nperm = 1000)
#the graph of this!
abline(-14.424822 ,3.5503772, col = "black")


#RPF vs Proteomics
rc_prot <- left_join(x = rc,y = prot, by = "GENE_NAME3") %>% na.omit() %>% .[.$rc_count > 2 & .$prot_count >2,]
plot(log2(rc_prot$rc_count), log2(rc_prot$prot_count), xlab = "Ribosomal Profiling FootPrint", ylab = "Proteomics", main = "log2 CDS Read Counts", pch = 19, cex = 0.5)
abline(lm(log2(rc_prot$prot_count) ~ log2(rc_prot$rc_count), data = rc_prot), col = "blue")
cor(log2(rc_prot$rc_count), log2(rc_prot$prot_count), method = c("pearson"))
# Ranged Major Axis Regression for RPF and Proteomics
lmodel2(log2(prot_count) ~ log2(rc_count),  data = rc_prot,range.y = "relative",range.x = "relative" )
#the graph of this!
abline(-5.471448,2.3189631 , col = "black")





