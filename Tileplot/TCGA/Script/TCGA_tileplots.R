##### TCGA tileplots ####

######################
### load libraries ###
######################
library(tidyverse)
library(data.table)
library(ggpubr)
library(mclust)
library()


### load clinical tcga data 
clin_tcga <- as.data.frame(fread("clinical_tcga_data.csv"))

### load tcga data 
tcga <- as.data.frame(fread("expression_data_tcga_log2.csv"))

### load trim data 
trim <- as.data.frame(fread("trim_overview_final.csv"))

############################
#### filter based on trim ##
############################

### filter tcga 
trim_tcga <- trim %>% filter(trim$gene %in% tcga$gene)
### 68 of 84 matches

#### look at which ones that was not included in the trim_expr
nonMatch <- trim %>% filter(!trim$gene %in% trim_tcga$gene)

### edit exel file with trim proteins and load file again with the new names for some trims
trim_edited <- as.data.frame(fread("trim_overview_edit_TCGA.csv"))

### merge trim and tcga 
trim_tcga <- merge(tcga, trim_edited, by = "gene")

### to get id on column 1 and trims on rows
trim_tcga <- column_to_rownames(trim_tcga, "gene")
trim_tcga <- as.data.frame(t(trim_tcga))
trim_tcga <- rownames_to_column(trim_tcga, "id")


### merge with clinical data 
colnames(clin_tcga)[2] <- "id" #to get the same name for id as in the other data

merged <- merge(trim_tcga, clin_tcga, by = "id")

table(merged$`tcga code`) #to check if all codes are included --> 33 codes, all codes are included 

####Test for duplicates 
dup <- as.data.frame(duplicated(merged$id))

dup_rows <- merged[dup$`duplicated(merged$id)`, ]
###gets 89 duplicates

dup_rows$id
###looked at the id in tcga dataset and saw that there are two of the same id with same values

### the duplicated ones has the same values, can delete the second occurrance
dup_results <- duplicated(merged$id)
merged <- merged[!dup_results, ] #before removing:10312, after: 10223 --> 89 removed

### no matches for duplicates 
dup_test <- as.data.frame(duplicated(merged$id))


#######################
## Calculate median ##
######################
### data has wide table format 
### format wide df to long because it is better to use in analysis 
merged_long <- pivot_longer(merged, cols = c(2:76)) 
### selects the columns that has trims 




### Calculate the median expression for each trim protein in each of the TCGA cancer code 
median_data <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(MedianExpression = median(value))


################
### tileplot ###
################
tileplot <- ggplot(median_data, aes(x = `tcga code`, y = factor(name, level = level_order), fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  stat_compare_means() + 
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())

pdf("tileplot_tcga.pdf", width = 15, height = 20, onefile = F)
print(tileplot)
dev.off()

#######################################
### sort based on median expression ###
#######################################
median_data$name <- reorder(median_data$name, median_data$MedianExpression)

tileplot_median <- ggplot(median_data, aes(x = `tcga code`, y = name, fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw() + 
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


pdf("tileplot_median_tcga.pdf", width = 15, height = 20, onefile = F)
print(tileplot_median)
dev.off()

####################################
### remove low median expression ###
####################################
# removed by looking at the figure:
median_data <- subset(median_data, !(name %in% c("TRIML1", "TRIM77", "TRIM64", 
                                                 "TRIM60", "TRIM49", "TRIM48", 
                                                 "TRIM43", "TRIM42", "TRIM40", 
                                                 "SPRYD5", "TRIM72")))

tileplot_high_median <- ggplot(median_data, aes(x = `tcga code`, y = name, fill = MedianExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw()+ 
  labs(y = " ", x = "Cancer Type", fill = "Median Expression") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


pdf("tileplot_high_median_tcga.pdf", width = 15, height = 20, onefile = F)
print(tileplot_high_median)
dev.off()



###################################################################
### look at the ones with the same expression in all cell types ###
###################################################################
### look at variance 
variance_data <- merged_long %>%
  group_by(name) %>%
  summarize(VarianceExpression = var(value)) 

### remove values under 2 ##
variance_data <- merged_long %>% group_by(name) %>% 
  summarize(VarianceExpression = var(value)) %>% filter(VarianceExpression > 2)

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

pdf("tileplot_high_variance_tcga.pdf", width = 15, height = 20, onefile = F)
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

pdf("tileplot_filtered_variance_tcga.pdf", width = 10, height = 10, onefile = F)
print(tileplot_filtered_variance)
dev.off()

# save variance for subsequent correlation analysis 
fwrite(filtered_variance, "filtered_variance_tcga.csv", append=TRUE)
