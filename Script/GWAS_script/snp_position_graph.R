# libraries ---------------------------------------------------------------
library(tidyverse)
library(data.table)

# input -------------------------------------------------------------------
# snp file
df.snp <- fread("files/snp_morphology_gwas_012")
df.snp <- df.snp %>%
  mutate(snp_name = paste(CHROM, POS, sep = "_"))

# phenotype file
df.phenotype <- fread("files/02_adults_variant_sites_gwas_maf_miss_sweeps.fam")

first_colnames <- as.vector(paste0("V",2:6))
last_colnames <- as.vector(paste0("V",31:39))
colnames <- fread("files/colnames_morphology", header = FALSE)
colnames <- c("id", first_colnames, colnames$V1, last_colnames)
colnames(df.phenotype) <- colnames

# significant GWAS results
summary_table <- fread("files/summary_suggestive_bw.txt")
summary_table <- summary_table %>%
  rowwise() %>%
  mutate(snp_name = paste0("chr", CHR, "_", BP),
        test_name = str_split(file_name,".assoc")[[1]][1]) %>%
  filter(!(snp_name %like% "chr34") & !(snp_name %like% "chr35"))


# box plot ----------------------------------------------------------------
# prep box plot colours
level_colour <- c("0" = "#987AA1", "1" = "#89BEBD", "2" = "#FDF08B")

for (i in 1:nrow(summary_table)){
  #if (i%%10==0) {print(i)}
  print(i)
  str_snp_name <- summary_table[i,]$snp_name
  str_test_name <- summary_table[i,]$test_name
  #make the file name to save the plots into
  file_name <- paste0("plot_output_morphology_sweeps/", str_snp_name, "_", str_test_name, ".jpg")
  #make a new data frame comprising of just the id and the phentoype being used, all_of command just makes sure you are taking the string and not just the name (i.e. that it is a whole column)
  phenotype <- df.phenotype %>%
    select(id, all_of(str_test_name))
  #the SNP is a matrix but it needs to be transposed t= transpose, filter out just the snp we want, then removes the column chrom, pos and snp_name
  snp <- as.data.frame(t(df.snp %>% filter(snp_name == str_snp_name) %>% select(-CHROM, -POS, -snp_name))) %>%
    rownames_to_column(., var = "id")
  
  #left_join is just a way of merging, joining phenotype and genotype together
  selection <- phenotype %>%
    left_join(snp, by = "id") %>%
    setNames(c("id", "phenotype", "genotype"))
  
  selection$genotype <- as.factor(selection$genotype)
  
  selection <- na.omit(selection)
  
  #following is so you can see the number in each genotype class
  xlabs <- paste(levels(selection$genotype),"\n(N=",table(selection$genotype),")",sep="")
  
  p1 <- selection %>%
    ggplot(aes(x = genotype, y = phenotype, fill = genotype)) +
    geom_boxplot() +
    scale_fill_manual(values = level_colour) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(str_snp_name, str_test_name, sep = "_")) +
    scale_x_discrete(labels=xlabs)
  
  p1
  
  ggsave(file_name,
         width = 18, height = 13, units = 'cm', bg='#ffffff')
}
