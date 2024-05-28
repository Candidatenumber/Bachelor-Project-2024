#### CCLE tileplot ####

####################
###load libraries###
####################
library(tidyverse)
library(data.table)
library(ggpubr)


### read expression data
expr <- as.data.frame(fread("expression_data_ccle_edited.csv", header=TRUE))

### read clinical data ccle (edited)
clin_ccle <- as.data.frame(fread("clinical_ccle_data.csv"))

### read in TRIM file 
trim <- as.data.frame(fread("trim_overview_edit_CCLE.csv"))

############################
#### filter based on trim ##
############################

#### merge files (trim and expr)
trim_expr <- merge(expr, trim, by = "gene")
### 72 of 84 trim before editing, but after its 82 of 84 

#### look at which ones that was not included in the trim_expr
# nonMatch_Uniquedf1 <- trim %>% filter(!trim$gene %in% trim_expr$gene)
#edit exel file and then load the file again ^^^^
 

### subset based on clinical data 
### column to row and row to column 
trim_expr <- column_to_rownames(trim_expr, "gene")
trim_expr <- as.data.frame(t(trim_expr))
trim_expr <- rownames_to_column(trim_expr, "id")

merged <- merge(trim_expr, clin_ccle, by.x = "id", by.y = "ModelID")

sort(unique(merged$OncotreePrimaryDisease))

### remove non-cancer
noncancer <- subset(merged, merged$OncotreePrimaryDisease == "Non-Cancerous")
### no non-cancer samples present 

### Clean merged file 
table(merged$`tcga code`) #29 unique cancer codes 

### Remove the NAs
### Those cell lines that do not have any corresponding TCGA code 
merged <- merged[complete.cases(merged$`tcga code`), ]


### Test for duplicates 
dup <- as.data.frame(duplicated(merged$id))
### no duplicates present 
test <- distinct(merged, merged$id,.keep_all = TRUE)
rm(dup, test)
### No duplictaes are present

#######################
## Calculate median ##
######################

### data has wide table format 
### want a long format rather than wide 
### format wide df to long  
merged_long <- pivot_longer(merged, cols = c(2:83))
### each patient for that gene gets a value ^


### Calculate the median expression for each trim protein in each of the TCGA cancer code 
median_data <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(MedianExpression = median(value))

####################
##### tileplot #####
####################
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

#####################################
#### remove low median expression ###
#####################################
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

###################################################################
### look at the ones with the same expression in all cell types ###
###################################################################
### look at variance  
variance_data <- merged_long %>%
  group_by(name) %>%
  summarize(VarianceExpression = var(value)) 

### remove low variance values ##
variance_data <- merged_long %>% group_by(name) %>% 
  summarize(VarianceExpression = var(value)) %>% filter(VarianceExpression > 1)

### make plot with high variance 
merged_high <- merge(median_data, variance_data, by = "name")

tileplot_high_variance <- ggplot(merged_high, aes(x = `tcga code`, y = name, fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


pdf("tileplot_high_variance_ccle.pdf", width = 15, height = 20, onefile = F)
print(tileplot_high_variance)
dev.off()

##### figure with the trims that has high variance in ccle and tcga 
variance_both <- c("MID1", "MID2", "RNF207", "TRIM15", "TRIM2", "TRIM22", 
                   "TRIM29", "TRIM55", "TRIM58", "TRIM6", "TRIM7", "TRIM9")

filtered_variance <- merged_high[merged_high$name %in% variance_both, ]

tileplot_filtered_variance <- ggplot(filtered_variance, aes(x = `tcga code`,  y = factor(name, level = variance_both), fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), 
        strip.background = element_blank())

pdf("tileplot_filtered_variance_ccle.pdf", width = 10, height = 10, onefile = F)
print(tileplot_filtered_variance)
dev.off()


