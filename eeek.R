# Requirements for the ribo R to run.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ribor")
install.packages("lmodel2")
install.packages("devtools")
install_github("ribosomeprofiling/ribor")
library("devtools")
library(ribor)
library(tidyr)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(ggplot2)
library("ggpubr")
library(lmodel2)

# Creating the ribo R object
second.ribo <- Ribo("/Users/mohanadelchouemi/Desktop/Dr Cenik/RibosomeTranslation/Cenik_Lab_ME/all.ribo", rename = rename_default )
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



# Downloading the proteomics data.v


# https://www.biostars.org/p/199073/

# Information abour RMA http://strata.uga.edu/8370/lecturenotes/regression.html

data <- read_xlsx(path = "/Users/mohanadelchouemi/Desktop/Dr Cenik/RibosomeTranslation/Cenik_Lab_ME/mmc4.xlsx", sheet = 2) %>% mutate_at(vars(2:3),~as.numeric(.x ))
data[is.na(data)] <- 0
data$prot_count <- (data$`iBAQ 1 protein copies per cell` + data$`iBAQ 2 protein copies per cell`) / 2
prot <- data %>% select(`gene symbol`,prot_count) %>% rename(GENE_NAME3 = `gene symbol`)

# RNA-Seq vs Proteomics
rna_prot <- left_join(x = rna,y = prot, by = "GENE_NAME3") %>% na.omit() %>% .[.$rna_count > 2 & .$prot_count >2,]
plot(log2(rna_prot$rna_count), log2(rna_prot$prot_count), xlab = "RNA-Seq", ylab = "Proteomics", main = "log2 CDS Read Counts", pch = 19, cex = 0.5)
abline(lm(log2(rna_prot$prot_count) ~ log2(rna_prot$rna_count), data = rna_prot), col = "blue")
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

#RPF_ALL_REGIONS vs Proteomics



?get_coverage
# aaaaabbbbbbbbbcccddfdgdg
# A dataframe of the ribosomal profiling for second.ribo

rc_raw <- get_region_counts(second.ribo,tidy = TRUE,
                        compact      = FALSE,
                        transcript  = FALSE,
                        alias       = TRUE,
                        experiment ="GSM1691202") %>% group_by(transcript) %>% summarise(total_count = sum(count))  #19736 unique transcripts


rc_raw[rc_raw$transcript == "SDF4-202",]

experiment.info <- get_info(ribo.object = second.ribo)[['experiment.info']]

# A data frame for the coverage of genes experiment GSM1691202
cov <- get_coverage(ribo.object = second.ribo,
                    name = "SDF4-202",
                    length      = TRUE,
                    alias       = TRUE,
                    tidy        = TRUE,
                    experiment = "GSM1691202",
                    compact = FALSE) 
# Big code to makes a get_gene_no_outliers
cov_new <- tibble(
  transcript = NA,
  count = NA
  
)

transcript_names <- unique(rc_raw$transcript)



for (i in 1:length(unique(rc_raw$transcript))) {
  naming_variable = transcript_names[i] 
  cov_temp <- get_coverage(ribo.object = second.ribo,
                           name = naming_variable,
                           length      = TRUE,
                           alias       = TRUE,
                           tidy        = TRUE,
                           experiment = "GSM1691202",
                           compact = FALSE) 
  nom <- sum(cov_temp$count[cov_temp$count > (mean(cov_temp$count) - 5*sd(cov_temp$count)) &
                       cov_temp$count < (mean(cov_temp$count) + 5*sd(cov_temp$count))] )
  cov_new <- add_row(.data = cov_new,transcript = transcript_names[i], count = nom)
  
}
# Verify that the code is running correctly, make a training data set! Always!
# USE ONLY PROTEOMICS TRANSCRIPTS
# REMOVE 5' AND 3' AND ANY OTHER NONE CDS REGION.
# Remember that 



# Do not run this again unless you want it to take 2 hours. This is the data with the gene coverage checked for and outlier removed.
cov_new2 <-cov_new[-1,]
write.csv(cov_new2, file = "gene_coverage_sum_sd_5.csv")

# Rerun the RPF_edited with the Proteomics data.

rpf_edited <- read.csv("gene_coverage_sum_sd_5.csv")  %>% mutate(GENE_NAME3 = str_replace(transcript , "-.*",replacement = "")) %>%
  select(-c(transcript)) %>% select(count,GENE_NAME3) %>% rename(rpf_edited_count = count)

rpf_edited_prot <- left_join(x = rpf_edited,y = prot, by = "GENE_NAME3") %>% na.omit() %>% .[.$rpf_edited_count > 2 & .$prot_count >2,]
plot(log2(rpf_edited_prot$rpf_edited_count), log2(rpf_edited_prot$prot_count), xlab = "Ribosomal Profiling FootPrint edited", ylab = "Proteomics", main = "log2 Read Counts", pch = 19, cex = 0.5)
abline(lm(log2(rpf_edited_prot$prot_count) ~ log2(rpf_edited_prot$rpf_edited_count), data = rpf_edited_prot), col = "blue")
cor(log2(rpf_edited_prot$rpf_edited_count), log2(rpf_edited_prot$prot_count), method = c("pearson"))
# Ranged Major Axis Regression for RPF and Proteomics
lmodel2(log2(prot_count) ~ log2(rpf_edited_count),  data = rpf_edited_prot,range.y = "relative",range.x = "relative" )
#the graph of this!
abline(-6.643835,2.438341 , col = "black")

