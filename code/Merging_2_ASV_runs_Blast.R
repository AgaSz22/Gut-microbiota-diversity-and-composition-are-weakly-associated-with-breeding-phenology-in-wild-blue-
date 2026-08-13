#Date: 28/07/2025
#Authors: Katy Maher & Agnes Sz.
#Purpose:To combine ASV tables from both Illumina runs (Blast)
#and introduce suitable column names matching metadata

### Combine two runs from DADA2

rm(list = ls())
library(dplyr)

#### Combine ASV count files ####
# Import counts matrix
counts1<-read.table("../data/ps_countmat_Katerun_a.tsv", row.names=NULL)
counts2<-read.table("../data/ps_countmat_Katerun_b.tsv", row.names=NULL)


### added

colnames(counts2) <- gsub("^X1_(\\d+\\..+)", "X\\1", colnames(counts2))
colnames(counts2) <- gsub("^X(\\d+)_", "X\\1.", colnames(counts2))

# names_col1 <- colnames(counts1)
# names_col2 <- colnames(counts2)
# 
# setdiff(names_col1, names_col2)

# combine runs and sum counts
counts.total<-aggregate(. ~ row.names, rbind(counts1, counts2), sum)

# write output
write.csv(counts.total, "../data/06_ASV_counts_both_runs_Blast.csv", row.names = F)

#### Combine taxonomy files ####

# Import taxonomy datasets
tax1<-read.table("../data/ps_taxamat_Katerun_a.tsv", row.names=1)
tax2<-read.table("../data/ps_taxamat_Katerun_b.tsv", row.names=1)

# Combine two datasets
tax.all <- rbind(tax1, tax2)

# remove duplicate rows
tax.total<-unique(tax.all)

# write output

#put lines of code for the separacte corrected ASVs

write.csv(tax.total, "../data/tax.total_Blast.csv", row.names = TRUE)
ASV_table <- read.csv("../data/06_ASV_counts_both_runs_Blast.csv")

###Checking whether I have all the samples###

# colnames(ASV_table)
# array <- colnames(ASV_table)
# filtered <- array[grep("Field", array)]
# filtered2 <- array[grep("PCR", array)]
# filtered3 <- array[grep("Extraction", array)]
# all_filtered <- c(filtered3,filtered2, filtered)
# 
# array2 <- setdiff(array, all_filtered)
# array2 <- gsub(".*\\.", "", array2)
# array2 <- array2[-1]
# array2 <- as.numeric(array2)

#Assigning the new names

ASV_table_cleaned <- ASV_table

new_colnames <- ifelse(
  grepl("PCR", colnames(ASV_table)),
  colnames(ASV_table),  # Keep "PCR" columns unchanged
  ifelse(
    colnames(ASV_table) == "X304.70", "S70.2",
    ifelse(colnames(ASV_table) == "X312.71", "S71.2",
           gsub("^X[^.]*\\.", "", colnames(ASV_table)) # Remove unwanted prefixes for others
    )
  )
)

# Making sure that column names are accepted by R (column names cannot start with numbers)

new_colnames <- trimws(new_colnames)
new_colnames <- make.names(new_colnames, unique = TRUE) #it introduces Xs
new_colnames <- gsub("^X", "S", new_colnames) #But we don't want Xs. Who needs an eX?
 
colnames(ASV_table_cleaned) <- new_colnames

#Saving the cleaned data frame

write.csv(
  ASV_table_cleaned,
  "../data/ASV_table_Blast.csv",
  row.names = FALSE,  #changed this to true and phyloseq stopped working xd
  quote = FALSE)
