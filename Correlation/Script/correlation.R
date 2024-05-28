#### correlation plot ####

####################
###load libraries###
####################
library(tidyverse)
library(data.table)
library(ggpubr)



### load data 
expr_ccle <- as.data.frame(fread("filtered_variance_ccle.csv"))
expr_tcga <- as.data.frame(fread("filtered_variance_tcga.csv"))

### to onnly get median values in dataframes 
expr_ccle <- expr_ccle[,c(1:3)]
expr_tcga <- expr_tcga[,c(1:3)]

### rename column
expr_ccle <- expr_ccle %>% rename(cancer_type = "tcga code")
expr_tcga <- expr_tcga %>% rename(cancer_type = "tcga code")
expr_ccle <- expr_ccle %>% rename(gene = "name")
expr_tcga <- expr_tcga %>% rename(gene = "name")



### want to compare median in trim from CCLE with trim from TCGA for the specific cancer types 
### make a matrix containing cancer type, gene, one column for ccle median value and one for tcga median value 

### make a column containing the database
expr_ccle$database <- "CCLE"
expr_tcga$database <- "TCGA"

### combine data 
combined <- rbind(expr_ccle, expr_tcga)

### convert to wide format 
wide_df <- pivot_wider(data = combined, names_from = database, values_from = MedianExpression)

### remove so we only get data for the 29 cancer types in ccle 
nonMatch <- expr_tcga %>% filter(!expr_tcga$cancer_type %in% expr_ccle$cancer_type)

### look at which ones that are not included
cancer_types_to_remove <- c("READ", "PCPG", "STAD", "THYM")
filtered_df <- wide_df[!wide_df$cancer_type %in% cancer_types_to_remove, ]

### select for gene in question
cor_data <- filtered_df %>% filter(gene %in% "MID1") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)

############
### MID1 ###
############
mid_1 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "MID1",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 10),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13)) 



pdf("correlation_MID1.pdf", width = 4, height = 4, onefile = F)
print(mid_1)
dev.off()

############
### MID2 ###
############
# select the gene in question 
cor_data <- filtered_df %>% filter(gene %in% "MID2") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
mid_2 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "MID2",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 10),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_MID2.pdf", width = 4, height = 4, onefile = F)
print(mid_2)
dev.off()

##############
### RNF207 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "RNF207") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
rnf_207 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "RNF207",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 8.5),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_RNF207.pdf", width = 4, height = 4, onefile = F)
print(rnf_207)
dev.off()

##############
### trim15 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "TRIM15") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t15 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM15",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 13.5),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM15.pdf", width = 4, height = 4, onefile = F)
print(t15)
dev.off()

#############
### TRIM2 ###
#############
cor_data <- filtered_df %>% filter(gene %in% "TRIM2") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t2 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
                 title = "TRIM2",
                 color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE,
                 cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 12.5),
                 xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM2.pdf", width = 4, height = 4, onefile = F)
print(t2)
dev.off()

##############
### TRIM22 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "TRIM22") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t22 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM22",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 11.5),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM22.pdf", width = 4, height = 4, onefile = F)
print(t22)
dev.off()

##############
### TRIM29 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "TRIM29") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t29 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM29",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 17.5),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM29.pdf", width = 4, height = 4, onefile = F)
print(t29)
dev.off()

##############
### TRIM55 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "TRIM55") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t55 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM55",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 9),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM55.pdf", width = 4, height = 4, onefile = F)
print(t55)
dev.off()

##############
### TRIM58 ###
##############
cor_data <- filtered_df %>% filter(gene %in% "TRIM58") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t58 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM58",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 10.5),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM58.pdf", width = 4, height = 4, onefile = F)
print(t58)
dev.off()

#############
### TRIM6 ###
#############
cor_data <- filtered_df %>% filter(gene %in% "TRIM6") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t6 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM6",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 9),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM6.pdf", width = 4, height = 4, onefile = F)
print(t6)
dev.off()

#############
### TRIM7 ###
#############
cor_data <- filtered_df %>% filter(gene %in% "TRIM7") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t7 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM7",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 8),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM7.pdf", width = 4, height = 4, onefile = F)
print(t7)
dev.off()

#############
### TRIM9 ###
#############
cor_data <- filtered_df %>% filter(gene %in% "TRIM9") %>% select(cancer_type, gene, TCGA, CCLE)

### high correlation = expression of trim is from cancer 
### low correlation = trim is expressed in other cell types (eks immune cells)
t9 <- ggscatter(cor_data, x = "CCLE", y = "TCGA",
          title = "TRIM9",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x = 0.2, label.y= 11),
          xlab = "", ylab = "")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic", size = 14), #Italic if it is a gene. 
        axis.text.x = element_text(size=13), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 13))

pdf("correlation_TRIM9.pdf", width = 4, height = 4, onefile = F)
print(t9)
dev.off()




