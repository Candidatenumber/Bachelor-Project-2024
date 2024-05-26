#### CCLE tileplot ####

####################
###load libraries###
####################

library(tidyverse)
library(data.table)
library(ggpubr)

########################
#read expression data###
########################
expr <- as.data.frame(fread("expression_data_ccle_edited.csv", header=TRUE))


#####################
# read clinical data ccle (edited)
#####################
clin_ccle <- as.data.frame(fread("clinical_ccle_data.csv"))

######################
# read in TRIM file 
#####################
trim <- as.data.frame(fread("trim_overview_edit_CCLE.csv"))

##### merge files (trim and expr)
trim_expr <- merge(expr, trim, by = "gene")
#### 72 of 84 trim before editing, but after its 82 of 84 

#### look at which ones that was not included in the trim_expr
# nonMatch_Uniquedf1 <- trim %>% filter(!trim$gene %in% trim_expr$gene)

#edit exel file and then load the file again ^^^^

### pasients in TRIM_expr with clin 
##############################
#subset based on clinical data 
##############################
#column to row and row to column 
trim_expr <- column_to_rownames(trim_expr, "gene")
trim_expr <- as.data.frame(t(trim_expr))
trim_expr <- rownames_to_column(trim_expr, "id")

merged <- merge(trim_expr, clin_ccle, by.x = "id", by.y = "ModelID")

sort(unique(merged$OncotreePrimaryDisease))

#remove non-cancer
noncancer <- subset(merged, merged$OncotreePrimaryDisease == "Non-Cancerous")
#no non-cancer samples present 



#######################
## Clean merged file ##
#######################

table(merged$`tcga code`) #29 unique cancer codes 


#Remove the NAs
#Those cell lines that do not have any corresponding TCGA code 
merged <- merged[complete.cases(merged$`tcga code`), ]


#Test for duplicates 
dup <- as.data.frame(duplicated(merged$id))
## no duplicates present 
test <- distinct(merged, merged$id,.keep_all = TRUE)
rm(dup, test)
#No duplictaes are present

###########################
## Calculate mean/median ##
###########################

# data has wide table format 
# want a long format rather than wide 
# format wide df to long 
#Format wide df til long df 
merged_long <- pivot_longer(merged, cols = c(2:83))
# each patient for that gene gets a value ^


#Calculate the median expression for each trim protein in each of the TCGA cancer code 
median_data <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(MedianExpression = median(value))


#to get the same order in this dotplot as in tcga
level_order <- c("CMYA5", "MEFV", "MID1", "MID2","PML", "RNF152", "RNF207", "SPRYD5", "TRIML1", "TRIM2", "TRIML2", "TRIM3", "TRIM4", 
                 "TRIM5", "TRIM6", "TRIM7", "TRIM8", "TRIM9", "TRIM10", "TRIM11", 
                 "TRIM13", "TRIM14", "TRIM15", "TRIM16", "TRIM16L", "TRIM17", "TRIM21", "TRIM22", 
                 "TRIM23", "TRIM24", "TRIM25", "TRIM26", "TRIM27", "TRIM28", "TRIM29",
                 "TRIM31", "TRIM32", "TRIM33", "TRIM34", "TRIM35", "TRIM36", "TRIM37", "TRIM38",
                 "TRIM39", "TRIM40", "TRIM41", "TRIM42", "TRIM43", "TRIM43B", "TRIM44", "TRIM45", "TRIM46",
                 "TRIM47", "TRIM48", "TRIM49", "TRIM49B", "TRIM49C", "TRIM49D2", "TRIM50", "TRIM51", "TRIM51G", 
                 "TRIM52", "TRIM54", "TRIM55", "TRIM56", "TRIM58", "TRIM59", "TRIM60", "TRIM61", "TRIM62", 
                 "TRIM63", "TRIM64", "TRIM64B", "TRIM64C", "TRIM65","TRIM66", "TRIM67", "TRIM68", "TRIM69", 
                 "TRIM71", "TRIM72", "TRIM73", "TRIM74",
                 "TRIM75", "TRIM77")


#####dotplot#####
dotplot <- ggplot(median_data, aes(x = `tcga code`, y = factor(name, level = level_order), fill = MedianExpression, size = MedianExpression)) +
  geom_point(shape=21) + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  stat_compare_means() +
  labs(x = "Cancer Type", y = " ", fill = "Median Expression") + theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),  
    strip.background = element_blank())


pdf("dotplot_ccle.pdf", width = 15, height = 20, onefile = F)
print(dotplot)
dev.off()

###figure with tiles istead 
tileplot <- ggplot(median_data, aes(x = `tcga code`, y = factor(name, level = level_order), fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  stat_compare_means() +
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 11),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


pdf("tileplot_ccle.pdf", width = 15, height = 20, onefile = F)
print(tileplot)
dev.off()

#######################################
### sort based on median expression ###
#######################################
median_data$name <- reorder(median_data$name, median_data$MedianExpression)

tileplot_median <- ggplot(median_data, aes(x = `tcga code`, y = name, fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


pdf("tileplot_median_ccle.pdf", width = 15, height = 20, onefile = F)
print(tileplot_median)
dev.off()
#####################################################################
#### remove low median expression, based on looking at the figure ###
#####################################################################
median_data <- subset(median_data, !(name %in% c("TRIM42", "TRIM64C", "TRIM77", "TRIML1", 
                                                 "TRIMG4", "TRIM64B", "TRIM49D2", "TRIM49B", 
                                                 "TRIM40", "TRIM49C", "TRIM43B", "TRIM50", 
                                                 "TRIM49", "TRIM48", "TRIM43", "TRIM72", 
                                                 "TRIM67", "MEFV", "TRIM60", "TRIM10", 
                                                 "TRIM54", "TRIM51", "TRIM74", "TRIM61", "TRIM64")))

tileplot_high_median <- ggplot(median_data, aes(x = `tcga code`, y = name, fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  theme_bw()+
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())

pdf("tileplot_high_median_ccle.pdf", width = 15, height = 20, onefile = F)
print(tileplot_high_median)
dev.off()