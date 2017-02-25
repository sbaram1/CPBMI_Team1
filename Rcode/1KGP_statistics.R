
library(dplyr)
setwd("~/Dropbox/BIBook/CPBMI_Team1")

###################################################################################
## Count Major and Minor allele count #############################################

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching <- read.csv("data/ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching")
Sample_ID_Asian_European <- read.csv("data/Sample_ID_Asian_European.csv")
Num_Subject <- dim(Sample_ID_Asian_European)[1]

allele_count <- apply(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching[,10:719], 1, function(x){
  allele_1st <- as.integer(x / 10)
  allele_2nd <- x %% 10
  allele <- c(allele_1st, allele_2nd) %>%
    factor(levels = c(0, 1, 2, 3, 4),
           labels = c("Major", "minor1", "minor2", "minor3", "minor4"))   %>%
    factor() 
  Major_allele_count <- sum(allele == "Major")
  Minor_allele_count <- sum(allele != "Major")
  minor1_allele_count <- sum(allele == "minor1")
  minor2_allele_count <- sum(allele == "minor2")
  minor3_allele_count <- sum(allele == "minor3")
  minor4_allele_count <- sum(allele == "minor4")
  
  c(Major_allele_count, Minor_allele_count, minor1_allele_count, 
    minor2_allele_count, minor3_allele_count, minor4_allele_count)
})

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$Major_allele_count <- 
  allele_count[1,]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$Minor_allele_count <- 
  allele_count[2,]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor1_allele_count <- 
  allele_count[3,]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor2_allele_count <- 
  allele_count[4,]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor3_allele_count <- 
  allele_count[5,]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor4_allele_count <- 
  allele_count[6,]

###################################################################################
## Remove rows with  ##############################################################

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching %>% nrow() #nrow = 12322

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching <- 
  ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching %>% 
  filter(Major_allele_count != 0, Minor_allele_count != 0) #nrow = 12302
  
sum(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor2_allele_count != 0) #95
sum(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor3_allele_count != 0) #10
sum(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$minor4_allele_count != 0) #1

###################################################################################
## Fisher Exact Test  #############################################################

## 2 alleles only
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele <- 
  ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching %>%
  filter(minor1_allele_count != 0, minor2_allele_count == 0) #12205
  
# x <- ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele[1,10:719]
pvals <- 
  apply(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele[,10:719], 1, function(x){
    allele_1st <- as.integer(x / 10)
    allele_2nd <- x %% 10
    allele <- c(allele_1st, allele_2nd) %>%
      factor(levels = c(0, 1, 2, 3, 4),
             labels = c("Major", "minor1", "minor2", "minor3", "minor4")) %>%
      factor() 
    race <- c(Sample_ID_Asian_European$Race, Sample_ID_Asian_European$Race)
    
    fisher.test(table(allele, race))$p.value
})

adj_pvals <- p.adjust(pvals, method = "fdr")
hist(adj_pvals)
sum(adj_pvals < 0.05) #3519
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele$adj_pvals <- adj_pvals

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele %>%
  filter(adj_pvals < 0.05) %>%
  select(POS, ID, HGNC.symbol) %>% head()

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching_2allele %>%
  filter(adj_pvals < 0.05) %>%
  select(HGNC.symbol) %>% unique()

###################################################################################
## Comparison between European and Asian ##########################################



