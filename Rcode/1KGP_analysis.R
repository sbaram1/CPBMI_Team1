
library(dplyr)

###############################################################################
## Read ch19_Minor_Allele_present.csv #########################################

setwd("~/Dropbox/BIBook/CPBMI_Team1/Rcode")

# ch19_Minor_Allele_present <- read.csv("../data/ch19_Minor_Allele_present.csv",
#                                       stringsAsFactors = F)[-1]
# 
# Sample_ID_Asian_European <- read.csv("../data/Sample_ID_Asian_European.csv",
#                                      stringsAsFactors = F)[-1]
# 
# dim(Sample_ID_Asian_European)
# dim(ch19_Minor_Allele_present)[2]
# 
# ch19_Minor_Allele_present_rm_pipe <-
#   ch19_Minor_Allele_present %>%
#   apply(2, function(x) gsub("|", "", x, fixed = T))
# names(ch19_Minor_Allele_present_rm_pipe) <- names(ch19_Minor_Allele_present)
# 
# # dim(ch19_Minor_Allele_present)
# # dim(ch19_Minor_Allele_present_rm_pipe)
# 
# ## ch19_Minor_Allele_present_rm_pipe <- cbind(ch19_Minor_Allele_present[1:9], ch19_Minor_Allele_present_rm_pipe)
# 
# setwd("~/Dropbox/BIBook/CPBMI_Team1/data")
# write.csv(ch19_Minor_Allele_present_rm_pipe, file = "ch19_Minor_Allele_present_rm_pipe.csv", row.names = F)
# 
# rm(ch19_Minor_Allele_present)

###############################################################################
## select gene in ch19 ########################################################

ch19_Minor_Allele_present_rm_pipe <- read.csv("../data/ch19_Minor_Allele_present_rm_pipe.csv",
                                      stringsAsFactors = F)
a <- ch19_Minor_Allele_present_rm_pipe[1:10, ] 

VitD_gene_ch19 <- read.csv("../data/PLOSone_gene_sets.csv", stringsAsFactors = F) 
VitD_gene_ch19 <- VitD_gene_ch19 %>%
  filter(Chromosome.Name == 19)

# str(VitD_gene_ch19)
str(ch19_Minor_Allele_present_rm_pipe)

VitD_gene_ch19$GRCh37.p13_start <- NA
VitD_gene_ch19$GRCh37.p13_end <- NA
tmp <- NA

for (i in 1:nrow(VitD_gene_ch19)) {
  tmp <- strsplit(VitD_gene_ch19$GRCh37.p13[i], "[: -]") %>% 
    unlist %>% gsub(",", "", .) %>% as.numeric()
  VitD_gene_ch19$GRCh37.p13_start[i] <- tmp[4]
  VitD_gene_ch19$GRCh37.p13_end[i] <- tmp[5]
}

# vitD_gene_index <- sapply(VitD_gene_ch19$GRCh37.p13_start, function(x) {
#   any(x >= VitD_gene_ch19$GRCh37.p13_start &&
#           x <= VitD_gene_ch19$GRCh37.p13_end)
# })
# 
# sum(vitD_gene_index)

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
library(IRanges)

# rng <- IRanges(start=c(4,7, 3), end=c(13, 11, 7))
# point <- IRanges(start = 6, end = 6)
# findOverlaps(rng, point) %>% summary %>% str

gene_range <- IRanges(start = VitD_gene_ch19$GRCh37.p13_start, 
                      end = VitD_gene_ch19$GRCh37.p13_end)
position <- IRanges(start = ch19_Minor_Allele_present_rm_pipe$POS,
                   end = ch19_Minor_Allele_present_rm_pipe$POS)
hits <- findOverlaps(gene_range, position)

ch19_Minor_Allele_present_rm_pipe$POS[subjectHits(hits)] %>% length()

VitD_gene_ch19$HGNC.symbol[queryHits(hits)]

ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching <- 
  ch19_Minor_Allele_present_rm_pipe[subjectHits(hits),]
ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching$HGNC.symbol <- VitD_gene_ch19$HGNC.symbol[queryHits(hits)]

write.csv(ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching, 
          file = "ch19_Minor_Allele_present_rm_pipe_VitD_gene_matching", row.names = F)

