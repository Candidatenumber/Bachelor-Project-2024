## boxplots ###

library(tidyverse)
library(data.table)
library(ggpubr)
library(mclust)



#load clinical tcga data 
clin_tcga <- as.data.frame(fread("clinical_tcga_data.csv"))

#load tcga data 
tcga <- as.data.frame(fread("expression_data_tcga_log2.csv"))

#load trim data 
trim <- as.data.frame(fread("trim_overview_final.csv"))

############################
#### filter based on trim ##
############################

#filter tcga 
trim_tcga <- trim %>% filter(trim$gene %in% tcga$gene)
#68 of 84 matches

##### look at which ones that was not included in the trim_expr
nonMatch <- trim %>% filter(!trim$gene %in% trim_tcga$gene)

### edit exel file with trim proteins and load file again with the new names for some trims
trim_edited <- as.data.frame(fread("trim_overview_edit_TCGA.csv"))

####################
#### merge files ####
#####################
#merge trim and tcga 
trim_tcga <- merge(tcga, trim_edited, by = "gene")

#to get id on column 1 and trims on rows
trim_tcga <- column_to_rownames(trim_tcga, "gene")
trim_tcga <- as.data.frame(t(trim_tcga))
trim_tcga <- rownames_to_column(trim_tcga, "id")


##############################
## subset for clinical data ##
#############################

#merge with clinical data 
colnames(clin_tcga)[2] <- "id" #to get the same name for id as in the other data

merged <- merge(trim_tcga, clin_tcga, by = "id")

#######################
## Clean merged file ##
#######################
table(merged$`tcga code`) #to check if all codes are included --> 33 codes, all codes are included 

#####Test for duplicates 
dup <- as.data.frame(duplicated(merged$id))

dup_rows <- merged[dup$`duplicated(merged$id)`, ]
#gets 89 duplicates

dup_rows$id
#looked at the id in tcga dataset and saw that there are two of the same id with same values

#the duplicated ones has the same values, can delete the second occurrance
dup_results <- duplicated(merged$id)
merged <- merged[!dup_results, ] #before removing:10312, after: 10223 --> 89 removed

#no matches for duplicates 
dup_test <- as.data.frame(duplicated(merged$id))


###########################
## Calculate mean/median ##
###########################
# data has wide table format --> meaning it has values that do not repeat in the first colum
# want a long format rather than wide --> values that do repeat in the first column.
# format wide df to long because it is better to use in analysis 
merged_long <- pivot_longer(merged, cols = c(2:76)) 
# selects the columns that has trims 



###############
#Calculate the median expression for each trim protein in each of the TCGA cancer code 
median_data <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(MedianExpression = median(value))



###################################################################
### look at the ones with the same expression in all cell types ###
###################################################################
# look at variance 
# gets variance for each trim 
variance_data <- merged_long %>%
  group_by(name) %>%
  summarize(VarianceExpression = var(value)) 

### remove values under 2 ##
variance_data <- merged_long %>% group_by(name) %>% 
  summarize(VarianceExpression = var(value)) %>% filter(VarianceExpression > 2)

# trim with high variance
merged_high <- merge(median_data, variance_data, by = "name")

#filter to get TRIM proteins with high variance in TCGA and CCLE
variance_both <- c("MID1", "MID2", "RNF207", "TRIM15", "TRIM2", "TRIM22", 
                   "TRIM29", "TRIM55", "TRIM58", "TRIM6", "TRIM7", "TRIM9")

filtered_variance <- merged_high[merged_high$name %in% variance_both, ]

#######################################################
### calculate std for each trim in each cancer type ###
#######################################################
std_trim <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(STDExpression = sd(value))

# calculate the std for the ones with high variance 
merged_std <- merge(std_trim, variance_data, by = "name")


# gets the 25 with high variance 
table(merged_std$name)


##########################################################
### STD for the 12 with high variance in CCLE and TCGA ###
##########################################################
high_var_std <- merge(filtered_variance, merged_std, by =c("name", "tcga code", "VarianceExpression"))

## look at which cancer type that has highest std in each of the 12 trim 

#TRIM55 --> UCS: 4.003219
#TRIM15 --> ESCA: 3.7849355
#TRIM58 --> SKCM: 2.7960201
#TRIM9 --> LIHC: 2.4005718
#TRIM7 --> COAD: 2.1105425
#TRIM6 --> SARC: 1.7701076
#TRIM29 --> SKCM: 3.7635371
#RNF207 --> GBM: 1.5523677
#MID2 --> SKCM: 2.2671384
#MID1 --> TGCT: 2.4981020
#TRIM22--> UVM: 1.7881138
#TRIM2 --> CESC: 1.7443290


######### want to look at expression value with std for the 12 trims of interest 
filtered_merged_long <- merged_long[merged_long$name %in% variance_both, ]

table(filtered_merged_long$name)
# gets the merged file for the 12 trims 

# to get std, value and variance together for the 12 trims
std_expr <- merge(filtered_merged_long, merged_std, by =c("name", "tcga code"))

table(std_expr$name) #check if all 12 trim was included (and only them)

std_expr <- std_expr[,c(1:5, 40:42)]


#make boxplots for the cancer type with the highest STD for each of the 12 trims 
#based on the values from merged long

#boxplot for each trim and the highest STD cancer type 
# histogram besides the boxplots to look at normal distribution

#######################
########trim15#########
#######################
t15 <- ggplot(std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")), aes(x = `tcga code`, y = value, fill = "TRIM15")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM15 Expression in ESCA") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#9B111E"))

histogram_t15 <- ggplot(std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")), aes(x = value)) +
  geom_density(fill = "#9B111E") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))




comb_t15 <- cowplot::plot_grid(t15 + theme(legend.position = "none"), 
                               histogram_t15 + coord_flip(), 
                               ncol = 2, align = 'h', rel_widths = c(2, 1))

?plot_grid
pdf("TRIM15_ESCA_STD.pdf", width = 4, height = 4, onefile = F)
print(t15)
dev.off()

pdf("TRIM15_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t15)
dev.off()

# to add a red line for cutoff value used in survival 
# GMM as cutoff: 
data <- std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")) %>%
  pull(value) # gets the data we want to look at

fit <- Mclust(data, G = 2) # modelling of data by using two gaussian mix models --> G is number of clusters 

means <- fit$parameters$mean # extracts the means of the two clusters from the fitted model 
#mean of each of the two models (in this case its two)

v <- mean(means) # mean of the cluster mean
# the value here is 4.928209 (does not directly get the minimum)
# it performs a clustering analysis and calculates the mean of the cluster means

hist_t15_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#40E0D0", size = 1) +
  geom_vline(xintercept = v, color = "#40E0D0", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) 


pdf("TRIM15_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t15_GMM)
dev.off()


#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 5.701184


hist_t15_median <- ggplot(std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#40E0D0", size = 1) +
  geom_vline(xintercept = median_value, color = "#40E0D0", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM15_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t15_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 5.068089


hist_t15_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM15") & `tcga code` %in% c("ESCA")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#40E0D0", size = 1 ) +
  geom_vline(xintercept = mean_value, color = "#40E0D0", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM15_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t15_mean)
dev.off()

############################
############trim55##########
############################
t55 <- ggplot(std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")), aes(x = `tcga code`, y = value, fill = "TRIM55")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM55 Expression in UCS") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#333ED4"))


histogram_t55 <- ggplot(std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")), aes(x = value)) +
  geom_density(fill = "#333ED4") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) 

comb_t55 <- cowplot::plot_grid(t55 + theme(legend.position = "none"), 
                               histogram_t55 + coord_flip(), 
                               ncol = 2, align = 'h', rel_widths = c(2, 1))

pdf("TRIM55_UCS_STD.pdf", width = 4, height = 4, onefile = F)
print(t55)
dev.off()

pdf("TRIM55_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t55)
dev.off()

# to add a red line for cutoff value used in survival 
# it is right skewed so has to: 
data <- std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")) %>%
  pull(value) # gets the data we want to look at

fit <- Mclust(data, G = 2) # modelling of data by using two gaussian mix models 

means <- fit$parameters$mean # mean of each of the two models 

v <- mean(means) # mean of the mean = 4.720237

hist_t55_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#333ED4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff9900", size = 1) +
  geom_vline(xintercept = v, color = "#ff9900", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) 


pdf("TRIM55_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t55_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 2.778082


hist_t55_median <- ggplot(std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#333ED4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff9900", size = 1) +
  geom_vline(xintercept = median_value, color = "#ff9900", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM55_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t55_median)
dev.off()


# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 4.520819


hist_t55_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM55") & `tcga code` %in% c("UCS")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#333ED4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff9900", size = 1) +
  geom_vline(xintercept = mean_value, color = "#ff9900", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM55_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t55_mean)
dev.off()

#########################
######## trim58 #########
#########################
t58 <- ggplot(std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")), aes(x = `tcga code`, y = value, fill = "TRIM58")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM58 Expression in SKCM") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#66ff00"))

histogram_t58 <- ggplot(std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_density(fill = "#66ff00") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t58 <- cowplot::plot_grid(t58 + theme(legend.position = "none"), 
                               histogram_t58 + coord_flip(), 
                               ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("TRIM58_SKCM_STD.pdf", width = 4, height = 4, onefile = F)
print(t58)
dev.off()

pdf("TRIM58_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t58)
dev.off()

# to add a red line for cutoff value used in survival 
data <- std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")) %>% 
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) # gets the value 5.285059


hist_t58_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#66ff00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff0099", size = 1) +
  geom_vline(xintercept = v, color = "#ff0099", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM58_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t58_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 4.025676


hist_t58_median <- ggplot(std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#66ff00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff0099", size = 1) +
  geom_vline(xintercept = median_value, color = "#ff0099", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM58_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t58_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 4.70558


hist_t58_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM58") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#66ff00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff0099", size = 1) +
  geom_vline(xintercept = mean_value, color = "#ff0099", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM58_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t58_mean)
dev.off()



####################
####### trim9 ######
####################
t9 <- ggplot(std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")), aes(x = `tcga code`, y = value, fill = "TRIM9")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM9 Expression in LIHC") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#cc00cc"))

histogram_t9 <- ggplot(std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")), aes(x = value)) +
  geom_density(fill = "#cc00cc") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t9 <- cowplot::plot_grid(t9 + theme(legend.position = "none"), 
                              histogram_t9 + coord_flip(), 
                              ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("TRIM9_LIHC_STD.pdf", width = 4, height = 4, onefile = F)
print(t9)
dev.off()

pdf("TRIM9_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t9)
dev.off()

# to add a red line for cutoff value used in survival 
data <- std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 3.729666

hist_t9_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#cc00cc", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#fff000", size = 1) +
  geom_vline(xintercept = v, color = "#fff000", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


pdf("TRIM9_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t9_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 3.375415


hist_t9_median <- ggplot(std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#cc00cc", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#fff000", size = 1) +
  geom_vline(xintercept = median_value, color = "#fff000", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM9_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t9_median)
dev.off()


# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 3.781376


hist_t9_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM9") & `tcga code` %in% c("LIHC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#cc00cc", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#fff000", size = 1) +
  geom_vline(xintercept = mean_value, color = "#fff000", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM9_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t9_mean)
dev.off()


#####################
###### trim7 #########
######################
t7 <- ggplot(std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")), aes(x = `tcga code`, y = value, fill = "TRIM7")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM7 Expression in COAD") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#ff9900"))

histogram_t7 <- ggplot(std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")), aes(x = value)) +
  geom_density(fill = "#ff9900") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 6.4),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t7 <- cowplot::plot_grid(t7 + theme(legend.position = "none"), 
                              histogram_t7 + coord_flip(), 
                              ncol = 2, align = 'h', rel_widths = c(2, 1))

pdf("TRIM7_COAD_STD.pdf", width = 4, height = 4, onefile = F)
print(t7)
dev.off()

pdf("TRIM7_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t7)
dev.off()

# to add a red line for cutoff value used in survival 
data <- std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 6.236155


hist_t7_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#ff9900", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "blue", size = 1) +
  geom_vline(xintercept = v, color = "blue", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


pdf("TRIM7_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t7_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 5.167258


hist_t7_median <- ggplot(std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#ff9900", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "blue", size = 1) +
  geom_vline(xintercept = median_value, color = "blue", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM7_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t7_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 5.673641


hist_t7_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM7") & `tcga code` %in% c("COAD")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#ff9900", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "blue", size = 1) +
  geom_vline(xintercept = mean_value, color = "blue", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM7_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t7_mean)
dev.off()



##################
####### trim6 ####
##################
t6 <- ggplot(std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")), aes(x = `tcga code`, y = value, fill = "TRIM6")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM6 Expression in SARC") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("red"))

histogram_t6 <- ggplot(std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")), aes(x = value)) +
  geom_density(fill = "red") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 6.4),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t6 <- cowplot::plot_grid(t6 + theme(legend.position = "none"), 
                              histogram_t6 + coord_flip(), 
                              ncol = 2, align = 'h', rel_widths = c(2, 1))

pdf("TRIM6_SARC_STD.pdf", width = 4, height = 4, onefile = F)
print(t6)
dev.off()

pdf("TRIM6_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t6)
dev.off()

# add a red line where the cutoff for survival is 
data <- std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 5.666906


hist_t6_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "red", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "green", size = 1) +
  geom_vline(xintercept = v, color = "green", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM6_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t6_GMM)
dev.off()


#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 6.480402


hist_t6_median <- ggplot(std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "red", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "green", size = 1) +
  geom_vline(xintercept = median_value, color = "green", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM6_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t6_median)
dev.off()


# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 6.215571


hist_t6_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM6") & `tcga code` %in% c("SARC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "red", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "green", size = 1) +
  geom_vline(xintercept = mean_value, color = "green", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM6_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t6_mean)
dev.off()


#################
#### TRIM29 #####
#################
t29 <- ggplot(std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")), aes(x = `tcga code`, y = value, fill = "TRIM29")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM29 Expression in SKCM") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#fff000"))

histogram_t29 <- ggplot(std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_density(fill = "#fff000") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t29 <- cowplot::plot_grid(t29 + theme(legend.position = "none"), 
                               histogram_t29 + coord_flip(), 
                               ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("TRIM29_SKCM_STD.pdf", width = 4, height = 4, onefile = F)
print(t29)
dev.off()

pdf("TRIM29_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t29)
dev.off()

# add a red line where the cutoff for survival is 
data <- std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 4.692118 


hist_t29_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#fff000", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff33ff", size = 1) +
  geom_vline(xintercept = v, color = "#ff33ff", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM29_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t29_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 2.151404


hist_t29_median <- ggplot(std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#fff000", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff33ff", size = 1) +
  geom_vline(xintercept = median_value, color = "#ff33ff", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM29_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t29_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 3.491103


hist_t29_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM29") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#fff000", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#ff33ff", size = 1) +
  geom_vline(xintercept = mean_value, color = "#ff33ff", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM29_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t29_mean)
dev.off()



##################
###### RNF207 ####
##################
rnf207 <- ggplot(std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")), aes(x = `tcga code`, y = value, fill = "RNF207")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "RNF207 Expression in GBM") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("cyan"))

histogram_rnf207 <- ggplot(std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")), aes(x = value)) +
  geom_density(fill = "cyan") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 7.4),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_rnf207 <- cowplot::plot_grid(rnf207 + theme(legend.position = "none"), 
                                  histogram_rnf207 + coord_flip(), 
                                  ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("RNF207_GBM_STD.pdf", width = 4, height = 4, onefile = F)
print(rnf207)
dev.off()

pdf("RNF207_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_rnf207)
dev.off()

# test for normality distribution
data <- std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")) %>%
  pull(value)

shapiro.test(data) # p value is 0,1573 so its greater than 0,05 and we can assume the normality.

#GMM method 
data <- std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 4.530978

hist_rnf207_GMM <- ggplot(std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "cyan", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "magenta", size = 1) +
  geom_vline(xintercept = v, color = "magenta", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("RNF207_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_rnf207_GMM)
dev.off()


#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 4.354782


hist_rnf207_median <- ggplot(std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "cyan", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "magenta", size = 1) +
  geom_vline(xintercept = median_value, color = "magenta", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("RNF207_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_rnf207_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 4.458414

hist_rnf207_mean <- ggplot(std_expr %>% filter(name %in% c("RNF207") & `tcga code` %in% c("GBM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "cyan", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "magenta", size = 1) +
  geom_vline(xintercept = mean_value, color = "magenta", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("RNF207_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_rnf207_mean)
dev.off()


####################
####### MID2 #######
####################
mid2 <- ggplot(std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")), aes(x = `tcga code`, y = value, fill = "MID2")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "MID2 Expression in SKCM") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#008080"))

histogram_mid2 <- ggplot(std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_density(fill = "#008080") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_mid2 <- cowplot::plot_grid(mid2 + theme(legend.position = "none"), 
                                histogram_mid2 + coord_flip(), 
                                ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("MID2_SKCM_STD.pdf", width = 4, height = 4, onefile = F)
print(mid2)
dev.off()

pdf("MID2_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_mid2)
dev.off()

#add a red line for cutoff value for survival 
# is bimodal distributed so has to calculate: 
data <- std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 6.460281


hist_mid2_GMM <- ggplot(std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#008080", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF7F50", size = 1) +
  geom_vline(xintercept = v, color = "#FF7F50", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("MID2_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_mid2_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 7.0506


hist_mid2_median <- ggplot(std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#008080", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF7F50", size = 1) +
  geom_vline(xintercept = median_value, color = "#FF7F50", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("MID2_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_mid2_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 6.688921


hist_mid2_mean <- ggplot(std_expr %>% filter(name %in% c("MID2") & `tcga code` %in% c("SKCM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#008080", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF7F50", size = 1) +
  geom_vline(xintercept = mean_value, color = "#FF7F50", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("MID2_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_mid2_mean)
dev.off()


###################
####### MID1 ######
###################
mid1 <- ggplot(std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")), aes(x = `tcga code`, y = value, fill = "MID1")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "MID1 Expression in TGCT") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#CCFF00"))

histogram_mid1 <- ggplot(std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")), aes(x = value)) +
  geom_density(fill = "#CCFF00") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_mid1 <- cowplot::plot_grid(mid1 + theme(legend.position = "none"), 
                                histogram_mid1 + coord_flip(), 
                                ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("MID1_TGCT_STD.pdf", width = 4, height = 4, onefile = F)
print(mid1)
dev.off()

pdf("MID1_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_mid1)
dev.off()

# add a red line to use for cutoff value for survival 
# has a bimodal distribution so: 
data <- std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) # gets the value 7.955581

hist_mid1_GMM <- ggplot(std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#CCFF00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF00FF", size = 1) +
  geom_vline(xintercept = v, color = "#FF00FF", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


pdf("MID1_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_mid1_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 7.708293


hist_mid1_median <- ggplot(std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#CCFF00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF00FF", size = 1) +
  geom_vline(xintercept = median_value, color = "#FF00FF", linetype = "dashed", size = 1.5) +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("MID1_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_mid1_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 7.842584


hist_mid1_mean <- ggplot(std_expr %>% filter(name %in% c("MID1") & `tcga code` %in% c("TGCT")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#CCFF00", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#FF00FF", size = 1) +
  geom_vline(xintercept = mean_value, color = "#FF00FF", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("MID1_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_mid1_mean)
dev.off()


#################
##### TRIM22 ####
#################
t22 <- ggplot(std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")), aes(x = `tcga code`, y = value, fill = "TRIM22")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM22 Expression in UVM") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#FF6F61"))

histogram_t22 <- ggplot(std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")), aes(x = value)) +
  geom_density(fill = "#FF6F61") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t22 <- cowplot::plot_grid(t22 + theme(legend.position = "none"), 
                               histogram_t22 + coord_flip(), 
                               ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("TRIM22_UVM_STD.pdf", width = 4, height = 4, onefile = F)
print(t22)
dev.off()

pdf("TRIM22_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t22)
dev.off()

# test for normality distribution
data <- (std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")) %>%
           pull(value))

shapiro.test(data) # p-value is 0,8876 and is therefore over 0,05, and we can assume the normality.

## GMM method 
data <- (std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")) %>%
           pull(value))

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) #gets the value 7.698415

hist_t22_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF6F61", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#6F00FF", size = 1) +
  geom_vline(xintercept = v, color = "#6F00FF", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


pdf("TRIM22_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t22_GMM)
dev.off()


#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 7.718867


hist_t22_median <- ggplot(std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF6F61", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#6F00FF", size = 1) +
  geom_vline(xintercept = median_value, color = "#6F00FF", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM22_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t22_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 7.674049


hist_t22_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM22") & `tcga code` %in% c("UVM")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF6F61", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#6F00FF", size = 1) +
  geom_vline(xintercept = mean_value, color = "#6F00FF", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM22_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t22_mean)
dev.off()


######################
####### TRIM2 ########
######################
t2 <- ggplot(std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")), aes(x = `tcga code`, y = value, fill = "TRIM2")) +
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA)+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.3, height = 0.5), shape = 21, size = 3) +
  labs(y = "Expression (log2+1)", title = "TRIM2 Expression in CESC") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text.x = element_text(size = 13), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#FF69B4"))

histogram_t2 <- ggplot(std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")), aes(x = value)) +
  geom_density(fill = "#FF69B4") +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))


comb_t2 <- cowplot::plot_grid(t2 + theme(legend.position = "none"), 
                              histogram_t2 + coord_flip(), 
                              ncol = 2, align = 'h', rel_widths = c(2, 1))


pdf("TRIM2_CESC_STD.pdf", width = 4, height = 4, onefile = F)
print(t2)
dev.off()

pdf("TRIM2_boxplot_hist.pdf", width = 5, height = 5, onefile = F)
print(comb_t2)
dev.off()

# to add a red line for cutoff value for survival 
data <- std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")) %>%
  pull(value)

fit <- Mclust(data, G = 2)

means <- fit$parameters$mean

v <- mean(means) # gets the value 8.602529

hist_t2_GMM <- ggplot(std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF69B4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#9F00FF", size = 1) +
  geom_vline(xintercept = v, color = "#9F00FF", linetype = "dashed", size = 1.5) +
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM2_histogram_GMM.pdf", width = 4, height = 4, onefile = F)
print(hist_t2_GMM)
dev.off()

#by using median as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")) %>% pull(value)
median_value <- median(data) # calculate the median, gives us the value 9.370452


hist_t2_median <- ggplot(std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF69B4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#9F00FF", size = 1) +
  geom_vline(xintercept = median_value, color = "#9F00FF", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM2_histogram_median.pdf", width = 4, height = 4, onefile = F)
print(hist_t2_median)
dev.off()

# using mean as the cutoff 
# to get the data we are interested in 
data <- std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")) %>% pull(value)
mean_value <- mean(data) # calculate the mean, gives us the value 9.184499


hist_t2_mean <- ggplot(std_expr %>% filter(name %in% c("TRIM2") & `tcga code` %in% c("CESC")), aes(x = value)) +
  geom_histogram(aes(y = ..density..), fill = "#FF69B4", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#9F00FF", size = 1) +
  geom_vline(xintercept = mean_value, color = "#9F00FF", linetype = "dashed", size = 1.5) +  # Add vertical line for v
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

pdf("TRIM2_histogram_mean.pdf", width = 4, height = 4, onefile = F)
print(hist_t2_mean)
dev.off()


