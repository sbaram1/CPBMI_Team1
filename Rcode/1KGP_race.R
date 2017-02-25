#####################################################################################
## Extract Asian and European from chra ch19 ########################################

library(dplyr)

setwd("~/Dropbox/BIBook/CPBMI_Team1/1KGP")

igsr_race <- read.csv("igsr_samples.tsv", sep = "\t")
head(igsr_race)
dim(igsr_race)
# names(igsr_race)
# unique(igsr_race$Superpopulation.code)

Asian <- igsr_race %>%
  filter(Population.code == "CHB" | Population.code == "JPT") %>%
  select(Sample.name) %>% unlist %>% as.character()


European <- igsr_race %>%
  filter(Superpopulation.code == "EUR") %>%
  select(Sample.name) %>% unlist %>% as.character()


heading <- read.csv("heading.txt", sep = "\t") %>% 
  colnames()
head(heading)
str(heading)

Asian_no <- which(heading %in% Asian)
# length(Asian_no) 207
# 1784-1886,1916-1985,2004-2037
# matrix(c(Asian_no, c(1, 1, 1)), nrow = 10)
names(Asian_no) <- rep(c("Asian"), length(Asian_no))
Asian_no

European_no <- which(heading %in% European)
head(European_no)
tail(European_no)
# length(European_no) 503
# 10-194,442,511-534,554-626,847-860,1666-1764,2304-2410
# matrix(c(European_no, rep(1, 7)), nrow = 10)
names(European_no) <- rep(c("European"), length(European_no))
European_no

Sample_ID_Asian_European <- c(Asian_no, European_no) %>% sort()
Sample_ID_Asian_European <- data.frame(ColID = Sample_ID_Asian_European, 
                                       Race = names(Sample_ID_Asian_European))
#####################################################################################
## Read ch19 ########################################################################

ch19 <- read.csv(file = "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.AsianEuropean.vcf",
                 sep = "\t", stringsAsFactors = F)

str(ch19)
dim(ch19)
names(ch19)[9]
names(ch19)[10]

ch19_mini <- ch19[1:5000,]

#####################################################################################
## Remove all 0|0 ###################################################################

# Major_Allele_only_index <- apply(ch19_mini[10:719], 1, function(x) all(x == "0|0"))
# ch19_mini_Minor_Allele_present <- ch19_mini[!Major_Allele_only_index,]

Major_Allele_only_index <- apply(ch19[10:719], 1, function(x) all(x == "0|0"))
ch19_Minor_Allele_present <- ch19[!Major_Allele_only_index,]

setwd("~/Dropbox/BIBook/CPBMI_Team1/data")

write.csv(ch19_Minor_Allele_present, file = "ch19_Minor_Allele_present.csv")
write.csv(Sample_ID_Asian_European, file = "Sample_ID_Asian_European.csv")

# # install.packages("vcfR")
# library(vcfR)
# 
# ch19 <- read.vcfR(file = "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.AsianEuropean.vcf.gz", verbose = F)
