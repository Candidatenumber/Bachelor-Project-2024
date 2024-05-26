
####################
###load libraries###
####################

library(tidyverse)
library(data.table)


#################################################
### CCLE preprocessing ###
#################################################

########################
#read expression data###
########################
expr <- as.data.frame(fread("OmicsExpressionProteinCodingGenesTPMLogp1.csv", header=TRUE))
colnames(expr)[1] <- "id"

##########################
# filter expr CCLE file 
######################
expr <- column_to_rownames(expr, "id")

new_colnames <- gsub("\\ .*","",colnames(expr))
colnames(expr) <- new_colnames

####transpose expr 
expr <- t(expr)

####to get a dataframe instead of matrix
expr <- as.data.frame(expr)

# genes - rows to col 
expr <- rownames_to_column(expr, "gene")

#save as a file to use further 
fwrite(expr, "expression_data_ccle_edited.csv")
?fwrite
test_1 <- fread("expression_data_ccle_edited.csv")




#################################################
### TCGA preprocessing ###
#################################################

##Read expression file
expr <- as.data.frame(fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"))
expr <- column_to_rownames(expr, "gene_id")
x <- as.data.frame(names(expr))

#log2 transform the data 
expr <- log2(expr+1)

#check for duplicates 
test <- substr(names(expr),start = 1,stop = 15)

# gets 11060 that is unique, meaning 9 out of the 11069 is duplicates
length(unique(test))

#look at the duplicated ones
test <- test[duplicated(test)]
test #9 duplicates 

# identify the duplicates and then remove the second occurrance of the duplicate
cols_remove <- c("TCGA-21-1076-01A-02R-0692-07", "TCGA-DD-AACA-02B-11R-A41C-07", 
                 "TCGA-06-0156-01A-03R-1849-01", "TCGA-06-0211-01B-01R-1849-01", 
                 "TCGA-DU-6404-02B-11R-A36H-07", "TCGA-DU-6407-02B-11R-A36H-07", 
                 "TCGA-FG-5965-02B-11R-A29R-07","TCGA-TQ-A7RK-02B-11R-A40A-07", 
                 "TCGA-23-1023-01R-01R-1564-13")

# remove the duplicates
expr <- expr[, !(names(expr) %in% cols_remove), drop = FALSE]
# has the 11060 samples now 


# look at the content of expr 
# -01A = 9,570
# -01B = 128
# -01C = 4
# -02A = 40 
# -02B = 1
# -02C = 0
# -03A = 158
# -03B = 15 
# -03C = 0 
# -04A = 0
# -04B = 0
# -04C= 0 
# -05A = 11
# -05B = 0
# -05C = 0 
# -06A = 393
# -06B = 2
# -06C = 0
# -07A = 1
# -07B = 0
# -07C = 0
# -08A = 0
# -08B = 0 
# -08C = 0 
# -09A = 0
# -09B = 0
# -09C = 0
# -10A = 0
# -10B = 0
# -10C = 0
# -11A = 719 
# -11B = 17 
# -11C = 1
# total = 11060 

# use gsub or sub here? 
# The sub() function applies for the first match. 
# The gsub() function applies for all matches
# have to remove everything behind -01A... etc to be able to match the expr with the clinical
colnames(expr) <- sub("-01A.*", "", colnames(expr))
colnames(expr) <- sub("-01B.*", "", colnames(expr))
colnames(expr) <- sub("-01C.*", "", colnames(expr))
colnames(expr) <- sub("-02A.*", "", colnames(expr))
colnames(expr) <- sub("-02B.*", "", colnames(expr))
colnames(expr) <- sub("-03A.*", "", colnames(expr))
colnames(expr) <- sub("-03B.*", "", colnames(expr))
colnames(expr) <- sub("-05A.*", "", colnames(expr))
colnames(expr) <- sub("-06A.*", "", colnames(expr))
colnames(expr) <- sub("-06B.*", "", colnames(expr))
colnames(expr) <- sub("-07A.*", "", colnames(expr))

# remove samples with 11, because they are Solid Tissue Normal --> 737
# before removing: 11060
# should be 10323 after removing
expr <- expr[, !grepl("-11A", colnames(expr))] # now: 10341
expr <- expr[, !grepl("-11B", colnames(expr))] # now 10324
expr <- expr[, !grepl("-11C", colnames(expr))] # now 10323


#can see that we have some non-unique values when trying to remove |. 
# has to look at ? and SLC35E2
# has to remove the rows with ?
expr <- rownames_to_column(expr, "gene")

#rows containing ? is = 29, before removing: 20531, after:20502
expr <- subset(expr, !grepl("\\?.*\\|", expr$gene))

# two of SLC35E2, edits the names based on the entrez id
expr$gene[expr$gene=="SLC35E2|728661"]<-"SLC35E2B"
expr$gene[expr$gene=="SLC35E2|9906"]<-"SLC35E2A"

# remove everything after | in the row names 
expr$gene <- sub("\\|.*", "", expr$gene)

#save as a file to use further 
fwrite(expr, "expression_data_tcga_log2.csv")
?fwrite
test_1 <- fread("expression_data_tcga_log2.csv")

?gsub
?column_to_rownames
