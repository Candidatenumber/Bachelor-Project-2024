########################
### look at survival ###
########################

#load new libraries 
library(tidyverse)
library(data.table)
library(ggpubr)
library(survival)
library(survminer)


# read expression file 
expr <- as.data.frame(fread("expression_data_tcga_log2.csv"))

#filter based on our 12 trim 
trim_edited_hv <-  as.data.frame(fread("12-trim-high-variance.csv"))

#subset sample for 12 trim 
expr <- merge(expr, trim_edited_hv, by ="gene")

expr <- column_to_rownames(expr, "gene")
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "sample")

#check for duplicates 
dup <- as.data.frame(duplicated(expr$sample))


#subset the gene of interest 
genes <- expr[,c("sample","TRIM9")]

## read clnical data 
clin_tcga <- as.data.frame(fread("clinical_tcga_data.csv"))

clin_tcga <- clin_tcga[,c(2:38)]
colnames(clin_tcga)[1] <- "sample"


#Merge the genes file with clinical data 
gene_clin <- merge(genes, clin_tcga, by = "sample")


#Subset cancer type of interest (the cancer type that had highest variance)
gene_clin_LIHC <- subset(gene_clin, gene_clin$`cancer type abbreviation` == "LIHC") #373 patients

#Subset info file to obtain DSS and PFI 
info_surv_DSS <- gene_clin_LIHC[,c(1,2,3,28,29)]
info_surv_PFI <- gene_clin_LIHC[,c(1,2,3,32,33)]
info_surv_DSS$cancer_type <- as.factor(info_surv_DSS$`cancer type abbreviation`)
info_surv_PFI$cancer_type <- as.factor(info_surv_PFI$`cancer type abbreviation`)

#Define TRIM9 expression based on high/low
info_surv_DSS$TRIM9_expression <- ifelse(info_surv_DSS$TRIM9 >= 3.729666, 'High', "Low")
info_surv_PFI$TRIM9_expression <- ifelse(info_surv_PFI$TRIM9 >= 3.729666, 'High', "Low")

table(info_surv_DSS$TRIM9_expression) #high vs low 
table(info_surv_PFI$TRIM9_expression) #high vs low 

#### for DSS #### 
#Remove timepoints that start at 0 
info_surv_DSS <- subset(info_surv_DSS, info_surv_DSS$DSS.time > 0) #before: 373, after 367

#Convert DSS.time into years 
info_surv_DSS$years <- info_surv_DSS$DSS.time/365

#remove NAs
#before removing: 367, after 359
info_surv_DSS <- info_surv_DSS[complete.cases(info_surv_DSS),]

#Define survival
survival = Surv(time= info_surv_DSS$years, event = info_surv_DSS$DSS)

survival_fit<- survfit(formula = survival ~ info_surv_DSS$TRIM9_expression, data = info_surv_DSS)

table(info_surv_DSS$cancer_type) #1, the one we are interested in 
table(info_surv_DSS$TRIM9_expression) #high vs low 


dss_t9 <- ggsurvplot(fit= survival_fit, 
           pval = TRUE, 
           surv.median.line = "hv", legend = c(0.1,0.2),
           xlab = "DSS (Years)", 
           ylab = "DSS Probability",
           title = "TRIM9 DSS for LIHC",
           ylim=c(0.0,1), 
           xlim=c(0,10.3),
           palette = c("#E69F00","#0072B2"),
           pval.coord=c(-0.06,0.05), 
           break.x.by= 1,         
           conf.int = T,
           risk.table = T, risk.table.title="",
           risk.table.height = 0.15,
           ncensor.plot = TRUE,
           ncensor.plot.height = 0.15,
           #legend = c(0.5, 0.95),
           legend.labs=c("High","Low"),
           legend.title= "TRIM9", 
           tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                axis.title.x = element_blank(), axis.title.y = element_blank(),
                                axis.text.y = element_text(size = 14)),font.x=14, font.y=14, font.tickslab=14)

pdf("dss_TRIM9.pdf", width = 8, height = 8, onefile = F)
print(dss_t9)
dev.off()

### TRIM9 PFI ###
#Remove timepoints that start at 0 for PFI
info_surv_PFI <- subset(info_surv_PFI, info_surv_PFI$PFI.time > 0) #before: 373, after 368

#Convert PFI.time into years 
info_surv_PFI$years <- info_surv_PFI$PFI.time/365

#remove NAs
#before removing: 368, after: 368
info_surv_PFI <- info_surv_PFI[complete.cases(info_surv_PFI),]

#Define survival (switch to PFI)
survival = Surv(time= info_surv_PFI$years, event = info_surv_PFI$PFI)

survival_fit<- survfit(formula = survival ~ info_surv_PFI$TRIM9_expression, data = info_surv_PFI)

pfi_t9 <- ggsurvplot(fit= survival_fit, 
           pval = TRUE, 
           surv.median.line = "hv", legend = c(0.1,0.2),
           xlab = "PFI (Years)", 
           ylab = "PFI Probability",
           title = "TRIM9 PFI for LIHC",
           ylim=c(0.0,1), 
           xlim=c(0,10.3),
           palette = c("#E69F00","#0072B2"),
           pval.coord=c(-0.06,0.05), 
           break.x.by= 1,         
           conf.int = T,
           risk.table = T, risk.table.title="",
           risk.table.height = 0.15,
           ncensor.plot = TRUE,
           ncensor.plot.height = 0.15,
           #legend = c(0.5, 0.95),
           legend.labs=c("High","Low"),
           legend.title= "TRIM9", 
           tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                axis.title.x = element_blank(), axis.title.y = element_blank(),
                                axis.text.y = element_text(size = 14)),font.x=14, font.y=14, font.tickslab=14)

pdf("pfi_trim9.pdf", width = 8, height = 8, onefile = F)
print(pfi_t9)
dev.off()
