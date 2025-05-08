# SIG AND SUG FILES -------------------------------------------------------
library(dplyr)
library(data.table)
library(readxl)
library(here)
library(ggpubr)
library(plotly)
library(reshape2)

#collecting all the files together
raw_gwas_files <- list.files()
raw_gwas_files_100kb_bw <- raw_gwas_files[raw_gwas_files%like%'100Kb' & raw_gwas_files %like% '_bw' & raw_gwas_files %like% 'assoc']

raw_gwas_files_100kb_nbw <- raw_gwas_files[raw_gwas_files%like%'100Kb' & raw_gwas_files %like% '_nbw' & raw_gwas_files %like% 'assoc']

raw_gwas_files_assoc <- raw_gwas_files[raw_gwas_files %like% 'assoc']
raw_gwas_files_bw <- raw_gwas_files_assoc[raw_gwas_files_assoc%like%'_bw'&!(raw_gwas_files_assoc %in% raw_gwas_files_100kb_bw)]
raw_gwas_files_nbw <- raw_gwas_files_assoc[raw_gwas_files_assoc%like%'_nbw'&!(raw_gwas_files_assoc %in% raw_gwas_files_100kb_nbw)]

# now we have data files with the names of the different analyses types (bw and no bw, 100kb interval and just sweep interval) need to wrote a loop to go through each
#use suggestive level of 5.6
sug_level <- 5.6

for (i in 1:24){
  file_name <- raw_gwas_files_100kb_bw[i]
  raw <- fread(file_name, sep = "\t", stringsAsFactors = FALSE)
  reduce <- as.data.frame(cbind(raw$chr, raw$ps, raw$rs, raw$p_lrt))
  reduce$V4 <- as.numeric(reduce$V4)
  ready <-reduce %>%
    filter(-log10(V4) > 1)
  ready$V3 <- paste(ready$V1, ready$V2, sep="_")
  ready$V1 <- ready$V1 %>%
    gsub("NC_006088.5", "1", .) %>% gsub("NC_006089.5", "2", .) %>%
    gsub("NC_006090.5", "3", .) %>% gsub("NC_006091.5", "4", .) %>%
    gsub("NC_006092.5", "5", .) %>% gsub("NC_006093.5", "6", .) %>%
    gsub("NC_006094.5", "7", .) %>% gsub("NC_006095.5", "8", .) %>%
    gsub("NC_006096.5", "9", .) %>% gsub("NC_006097.5", "10", .) %>%
    gsub("NC_006098.5", "11", .) %>% gsub("NC_006099.5", "12", .) %>%
    gsub("NC_006100.5", "13", .) %>% gsub("NC_006101.5", "14", .) %>%
    gsub("NC_006102.5", "15", .) %>% gsub("NC_006103.5", "16", .) %>%
    gsub("NC_006104.5", "17", .) %>% gsub("NC_006105.5", "18", .) %>%
    gsub("NC_006106.5", "19", .) %>% gsub("NC_006107.5", "20", .) %>%
    gsub("NC_006108.5", "21", .) %>% gsub("NC_006109.5", "22", .) %>%
    gsub("NC_006110.5", "23", .) %>% gsub("NC_006111.5", "24", .) %>%
    gsub("NC_006112.4", "25", .) %>% gsub("NC_006113.5", "26", .) %>%
    gsub("NC_006114.5", "27", .) %>% gsub("NC_006115.5", "28", .) %>%
    gsub("NC_028739.2", "30", .) %>% gsub("NC_028740.2", "31", .) %>%
    gsub("NC_006119.4", "32", .) %>% gsub("NC_008465.4", "33", .) %>%
    gsub("NC_006126.5", "34", .) %>% gsub("NC_006127.5", "35", .)
  #W=34 Z=35
  
  colnames(ready) <- c("CHR", "BP", "SNP", "P")
  
  ready$CHR <- as.numeric(ready$CHR)
  ready$BP <- as.numeric(ready$BP)
  
  ready$log <- -log10(ready$P)
  ready_suggestive <- subset(ready, log >= sug_level)
  #ready_suggestive$file_name <- file_name
  

    if (nrow(ready_suggestive) > 0){
      ready_suggestive$file_name <- file_name
      if (exists("summary_suggestive")){
        summary_suggestive <- rbind(summary_suggestive, ready_suggestive)
      } else { summary_suggestive <- ready_suggestive
      }
    }
}
  
summary_suggestive_100kb_bw <- summary_suggestive
rm(summary_suggestive)
# repeat above for 100kb but no body weight covariate

sug_level <- 5.6

for (i in 1:25){
  file_name <- raw_gwas_files_100kb_nbw[i]
  raw <- fread(file_name, sep = "\t", stringsAsFactors = FALSE)
  reduce <- as.data.frame(cbind(raw$chr, raw$ps, raw$rs, raw$p_lrt))
  reduce$V4 <- as.numeric(reduce$V4)
  ready <-reduce %>%
    filter(-log10(V4) > 1)
  ready$V3 <- paste(ready$V1, ready$V2, sep="_")
  ready$V1 <- ready$V1 %>%
    gsub("NC_006088.5", "1", .) %>% gsub("NC_006089.5", "2", .) %>%
    gsub("NC_006090.5", "3", .) %>% gsub("NC_006091.5", "4", .) %>%
    gsub("NC_006092.5", "5", .) %>% gsub("NC_006093.5", "6", .) %>%
    gsub("NC_006094.5", "7", .) %>% gsub("NC_006095.5", "8", .) %>%
    gsub("NC_006096.5", "9", .) %>% gsub("NC_006097.5", "10", .) %>%
    gsub("NC_006098.5", "11", .) %>% gsub("NC_006099.5", "12", .) %>%
    gsub("NC_006100.5", "13", .) %>% gsub("NC_006101.5", "14", .) %>%
    gsub("NC_006102.5", "15", .) %>% gsub("NC_006103.5", "16", .) %>%
    gsub("NC_006104.5", "17", .) %>% gsub("NC_006105.5", "18", .) %>%
    gsub("NC_006106.5", "19", .) %>% gsub("NC_006107.5", "20", .) %>%
    gsub("NC_006108.5", "21", .) %>% gsub("NC_006109.5", "22", .) %>%
    gsub("NC_006110.5", "23", .) %>% gsub("NC_006111.5", "24", .) %>%
    gsub("NC_006112.4", "25", .) %>% gsub("NC_006113.5", "26", .) %>%
    gsub("NC_006114.5", "27", .) %>% gsub("NC_006115.5", "28", .) %>%
    gsub("NC_028739.2", "30", .) %>% gsub("NC_028740.2", "31", .) %>%
    gsub("NC_006119.4", "32", .) %>% gsub("NC_008465.4", "33", .) %>%
    gsub("NC_006126.5", "34", .) %>% gsub("NC_006127.5", "35", .)
  #W=34 Z=35
  
  colnames(ready) <- c("CHR", "BP", "SNP", "P")
  
  ready$CHR <- as.numeric(ready$CHR)
  ready$BP <- as.numeric(ready$BP)
  
  ready$log <- -log10(ready$P)
  ready_suggestive <- subset(ready, log >= sug_level)
  # ready_suggestive$file_name <- file_name
  
  
  if (nrow(ready_suggestive) > 0){
    ready_suggestive$file_name <- file_name
    if (exists("summary_suggestive")){
      summary_suggestive <- rbind(summary_suggestive, ready_suggestive)
    } else { summary_suggestive <- ready_suggestive
    }
  }
}

summary_suggestive_100kb_nbw <- summary_suggestive
rm(summary_suggestive)

# same again but for sweep regions alone (no extra 100kb intervals) with body weight

sug_level <- 5.0

for (i in 1:24){
  file_name <- raw_gwas_files_bw[i]
  raw <- fread(file_name, sep = "\t", stringsAsFactors = FALSE)
  reduce <- as.data.frame(cbind(raw$chr, raw$ps, raw$rs, raw$p_lrt))
  reduce$V4 <- as.numeric(reduce$V4)
  ready <-reduce %>%
    filter(-log10(V4) > 1)
  ready$V3 <- paste(ready$V1, ready$V2, sep="_")
  ready$V1 <- ready$V1 %>%
    gsub("NC_006088.5", "1", .) %>% gsub("NC_006089.5", "2", .) %>%
    gsub("NC_006090.5", "3", .) %>% gsub("NC_006091.5", "4", .) %>%
    gsub("NC_006092.5", "5", .) %>% gsub("NC_006093.5", "6", .) %>%
    gsub("NC_006094.5", "7", .) %>% gsub("NC_006095.5", "8", .) %>%
    gsub("NC_006096.5", "9", .) %>% gsub("NC_006097.5", "10", .) %>%
    gsub("NC_006098.5", "11", .) %>% gsub("NC_006099.5", "12", .) %>%
    gsub("NC_006100.5", "13", .) %>% gsub("NC_006101.5", "14", .) %>%
    gsub("NC_006102.5", "15", .) %>% gsub("NC_006103.5", "16", .) %>%
    gsub("NC_006104.5", "17", .) %>% gsub("NC_006105.5", "18", .) %>%
    gsub("NC_006106.5", "19", .) %>% gsub("NC_006107.5", "20", .) %>%
    gsub("NC_006108.5", "21", .) %>% gsub("NC_006109.5", "22", .) %>%
    gsub("NC_006110.5", "23", .) %>% gsub("NC_006111.5", "24", .) %>%
    gsub("NC_006112.4", "25", .) %>% gsub("NC_006113.5", "26", .) %>%
    gsub("NC_006114.5", "27", .) %>% gsub("NC_006115.5", "28", .) %>%
    gsub("NC_028739.2", "30", .) %>% gsub("NC_028740.2", "31", .) %>%
    gsub("NC_006119.4", "32", .) %>% gsub("NC_008465.4", "33", .) %>%
    gsub("NC_006126.5", "34", .) %>% gsub("NC_006127.5", "35", .)
  #W=34 Z=35
  
  colnames(ready) <- c("CHR", "BP", "SNP", "P")
  
  ready$CHR <- as.numeric(ready$CHR)
  ready$BP <- as.numeric(ready$BP)
  
  ready$log <- -log10(ready$P)
  ready_suggestive <- subset(ready, log >= sug_level)
  #ready_suggestive$file_name <- file_name
  
  
  if (nrow(ready_suggestive) > 0){
    ready_suggestive$file_name <- file_name
    if (exists("summary_suggestive")){
      summary_suggestive <- rbind(summary_suggestive, ready_suggestive)
    } else { summary_suggestive <- ready_suggestive
    }
  }
}

summary_suggestive_bw <- summary_suggestive
rm(summary_suggestive)

# sweep regions alone, no bw covariate

sug_level <- 5.0

for (i in 1:25){
  file_name <- raw_gwas_files_nbw[i]
  raw <- fread(file_name, sep = "\t", stringsAsFactors = FALSE)
  reduce <- as.data.frame(cbind(raw$chr, raw$ps, raw$rs, raw$p_lrt))
  reduce$V4 <- as.numeric(reduce$V4)
  ready <-reduce %>%
    filter(-log10(V4) > 1)
  ready$V3 <- paste(ready$V1, ready$V2, sep="_")
  ready$V1 <- ready$V1 %>%
    gsub("NC_006088.5", "1", .) %>% gsub("NC_006089.5", "2", .) %>%
    gsub("NC_006090.5", "3", .) %>% gsub("NC_006091.5", "4", .) %>%
    gsub("NC_006092.5", "5", .) %>% gsub("NC_006093.5", "6", .) %>%
    gsub("NC_006094.5", "7", .) %>% gsub("NC_006095.5", "8", .) %>%
    gsub("NC_006096.5", "9", .) %>% gsub("NC_006097.5", "10", .) %>%
    gsub("NC_006098.5", "11", .) %>% gsub("NC_006099.5", "12", .) %>%
    gsub("NC_006100.5", "13", .) %>% gsub("NC_006101.5", "14", .) %>%
    gsub("NC_006102.5", "15", .) %>% gsub("NC_006103.5", "16", .) %>%
    gsub("NC_006104.5", "17", .) %>% gsub("NC_006105.5", "18", .) %>%
    gsub("NC_006106.5", "19", .) %>% gsub("NC_006107.5", "20", .) %>%
    gsub("NC_006108.5", "21", .) %>% gsub("NC_006109.5", "22", .) %>%
    gsub("NC_006110.5", "23", .) %>% gsub("NC_006111.5", "24", .) %>%
    gsub("NC_006112.4", "25", .) %>% gsub("NC_006113.5", "26", .) %>%
    gsub("NC_006114.5", "27", .) %>% gsub("NC_006115.5", "28", .) %>%
    gsub("NC_028739.2", "30", .) %>% gsub("NC_028740.2", "31", .) %>%
    gsub("NC_006119.4", "32", .) %>% gsub("NC_008465.4", "33", .) %>%
    gsub("NC_006126.5", "34", .) %>% gsub("NC_006127.5", "35", .)
  #W=34 Z=35
  
  colnames(ready) <- c("CHR", "BP", "SNP", "P")
  
  ready$CHR <- as.numeric(ready$CHR)
  ready$BP <- as.numeric(ready$BP)
  
  ready$log <- -log10(ready$P)
  ready_suggestive <- subset(ready, log >= sug_level)
  #ready_suggestive$file_name <- file_name
  
  if (nrow(ready_suggestive) > 0){
    ready_suggestive$file_name <- file_name
    if (exists("summary_suggestive")){
      summary_suggestive <- rbind(summary_suggestive, ready_suggestive)
    } else { summary_suggestive <- ready_suggestive
    }
  }
}

summary_suggestive_nbw <- summary_suggestive
rm(summary_suggestive)
## Significance threshold
## Reference: Guo, Y., Huang, Y., Hou, L., Ma, J., Chen, C., Ai, H., ... & Ren, J. (2017). Genome-wide detection of genetic markers associated with growth and fatness in four pig populations using four approaches. Genetics Selection Evolution, 49(1), 21.

n <- 430000 # number of snps used in hypothalamus GWAS
sig_level <- -log10(0.05/n)
sug_level <- -log10(1/n)

## ANALYSIS PERFORMED USING P_SCORES FROM GEMMA STBD1
raw <- fread("sweeps_adults_bw_bill_length_100kb.assoc.txt", sep = "\t", stringsAsFactors = FALSE)

#raw %>% slice(which.min(`p_score`))
#raw %>% filter(chr == "NC_006091.5", ps == 49883749)

reduce <- as.data.frame(cbind(raw$chr, raw$ps, raw$rs, raw$p_lrt))
reduce$V4 <- as.numeric(reduce$V4)
ready <-reduce %>%
  filter(-log10(V4) > 1)
ready$V3 <- paste(ready$V1, ready$V2, sep="_")
ready$V1 <- ready$V1 %>%
  gsub("NC_006088.5", "1", .) %>% gsub("NC_006089.5", "2", .) %>%
  gsub("NC_006090.5", "3", .) %>% gsub("NC_006091.5", "4", .) %>%
  gsub("NC_006092.5", "5", .) %>% gsub("NC_006093.5", "6", .) %>%
  gsub("NC_006094.5", "7", .) %>% gsub("NC_006095.5", "8", .) %>%
  gsub("NC_006096.5", "9", .) %>% gsub("NC_006097.5", "10", .) %>%
  gsub("NC_006098.5", "11", .) %>% gsub("NC_006099.5", "12", .) %>%
  gsub("NC_006100.5", "13", .) %>% gsub("NC_006101.5", "14", .) %>%
  gsub("NC_006102.5", "15", .) %>% gsub("NC_006103.5", "16", .) %>%
  gsub("NC_006104.5", "17", .) %>% gsub("NC_006105.5", "18", .) %>%
  gsub("NC_006106.5", "19", .) %>% gsub("NC_006107.5", "20", .) %>%
  gsub("NC_006108.5", "21", .) %>% gsub("NC_006109.5", "22", .) %>%
  gsub("NC_006110.5", "23", .) %>% gsub("NC_006111.5", "24", .) %>%
  gsub("NC_006112.4", "25", .) %>% gsub("NC_006113.5", "26", .) %>%
  gsub("NC_006114.5", "27", .) %>% gsub("NC_006115.5", "28", .) %>%
  gsub("NC_028739.2", "30", .) %>% gsub("NC_028740.2", "31", .) %>%
  gsub("NC_006119.4", "32", .) %>% gsub("NC_008465.4", "33", .) %>%
  gsub("NC_006126.5", "34", .) %>% gsub("NC_006127.5", "35", .)
#W=34 Z=35

colnames(ready) <- c("CHR", "BP", "SNP", "P")

ready$CHR <- as.numeric(ready$CHR)
ready$BP <- as.numeric(ready$BP)

# Make the Manhattan plot
manhattan(na.omit(ready),
          main = paste("STBD1", "manhattan", sep="_"),
          col = c("blue4", "orange3"),
          suggestiveline = -log10(1/11991659), genomewideline = -log10(0.05/11991659))

# Make the Manhattan plot
qq(na.omit(ready$P))


# Table with significant and suggestive snps
ready$log <- -log10(ready$P)
ready_suggestive <- subset(ready, log >= sug_level)

write.table(ready_suggestive, "hypothalamus/STBD1_sug.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

ready_sig <-subset(ready, log >= sig_level)
write.table(ready_sig, "hypothalamus/STBD1_sig.txt", sep = "\t",row.names = TRUE, col.names = NA)

stbd1_gwas_sig <- ready_sig
stbd1_gwas_sug <- ready_suggestive
stbd1_eqtl <- stbd1_snps





# Merging collected results with sweep information
sweeps <- read.table(file="sweeps_to_extract.txt", header=FALSE)

# This merges the sweep table with the GWAS results, meaning that every sweep for a chromosome is saved to every combination of GWAS QTL, then we need to remove the ones that fall outside of the specific sweep
sweeps_GWAS <- merge(sweeps, summary_suggestive_nbw, by.x="V1", by.y="CHR")
unique(sweeps_GWAS$SNP)

# This then pipes through the sweeps_GWAS file to then retain only those with a BP>= V2 and <=V3 (i.e. those that fall within a sweep)
sweeps_GWAS_nbw <- sweeps_GWAS %>% 
  filter(BP >= V2, BP <= V3) %>%
  mutate(sweep = paste0(V2,"-",V3))

# for all traits
sweeps_GWAS <- merge(sweeps, summary_suggestive_bw, by.x="V1", by.y="CHR")
sweeps_GWAS_bw <- sweeps_GWAS %>% 
  filter(BP >= V2, BP <= V3) %>%
  mutate(sweep = paste0(V2,"-",V3))


# for 100kb sweeps
sweeps <- sweeps %>%
  mutate(start_100 = V2-100000, stop_100 =V3+100000)

sweeps_GWAS <- merge(sweeps, summary_suggestive_100kb_nbw, by.x="V1", by.y="CHR")
sweeps_GWAS_100kb_nbw <- sweeps_GWAS %>% 
  filter(BP >= start_100, BP <= stop_100) %>%
  mutate(sweep = paste0(V2,"-",V3))


sweeps_GWAS <- merge(sweeps, summary_suggestive_100kb_bw, by.x="V1", by.y="CHR")
sweeps_GWAS_100kb_bw <- sweeps_GWAS %>% 
  filter(BP >= start_100, BP <= stop_100) %>%
  mutate(sweep = paste0(V2,"-",V3))


# retain just the top SNP per trait per sweep
sweeps_GWAS_100kb_bw_top <- sweeps_GWAS_100kb_bw %>%
  group_by(V1,sweep,file_name) %>%
  slice_min(order_by = P) %>%
  slice_min(order_by= BP)

sweeps_GWAS_100kb_nbw_top <- sweeps_GWAS_100kb_nbw %>%
  group_by(V1,sweep,file_name) %>%
  slice_min(order_by = P) %>%
  slice_min(order_by= BP)

sweeps_GWAS_bw_top <- sweeps_GWAS_bw %>%
  group_by(V1,sweep,file_name) %>%
  slice_min(order_by = P) %>%
  slice_min(order_by= BP)

sweeps_GWAS_nbw_top <- sweeps_GWAS_nbw %>%
  group_by(V1,sweep,file_name) %>%
  slice_min(order_by = P) %>%
  slice_min(order_by= BP)



