## ML workflow: model generation and prediction

## ---- Setup----------------------------------------------------------
setwd("~/Documents/Forskerlinje/Bioinformatics_mammatumor/Survival_paper")
library(tidyverse)
library(DESeq2)
library(splitstackshape)
library(survival)
library(glmnet)
library(c060)
library(caret)
source("helper.R") 
library(VIM)
library(faux)
library(corrplot)
library(riskRegression)
library(ggpubr)
load("data/tumor_samples_with_subtype.Rdata")
load("data/counts_tumors.Rdata")

## ------Formatting annotation file--------------------------------------
# Formatting: making sure variables that we will use later are the correct data type, imputing missing values with KNN and dropping dogs with no Survival days information. 
clin_vars = c("Age", "ER_status", "Lymphatic_invasion")

clinical_tbl <- tumor_samples_with_subtype %>% 
  mutate(subtype_MPAM50 = as.factor(subtype_MPAM50)) %>% 
  kNN(variable = clin_vars, k=5) %>% 
  drop_na("Survival_days") %>% 
  mutate(Age = as.vector(scale(Age))) %>%
  select(c(Sample_ID, all_of(clin_vars), Status_survival, Survival_days)) %>%  
  as_tibble()

# Onehot encoding
clinical_onehot <- clinical_tbl %>%
  mutate(ER_positive = ifelse(ER_status == "N", 0, 1)) %>%
  mutate(Lymphatic_invasion_present = ifelse(Lymphatic_invasion == "Absent", 0, 1)) %>%
  select(-ER_status, -Lymphatic_invasion)


## -------Formatting RNAseq count file ---------------------------------
# Making sure counts_tumors (Integer matrix with number of counts per gene in each sample) columns are the same as clinical_tbl Sample ID column - this is essential since DESeq doesn't make any guesses###
counts_tumors <- counts_tumors %>% 
  select(clinical_tbl$Sample_ID)
all(colnames(counts_tumors) == clinical_tbl$Sample_ID) #true --> same order of dogIDs

#Preprocessing: normalizing, variance filtering, removing highly correlated genes, log2 transforming, scaling and univariate cox filtering. 
preprocess_res <- preprocess_genes(x = counts_tumors, y=clinical_tbl, variance_cutoff = 1000, corr_cutoff = 0.95, fdr_cutoff = 0.1)

results_univariate_genes <- as.data.frame(preprocess_res$univariate_result) %>% 
  rownames_to_column(var="Gene names") %>% 
  rename(Cox.fdr = "preprocess_res$univariate_result") %>% 
  mutate(Cox.sign = ifelse(Cox.fdr<0.1, "Yes", "No"))
  
writexl::write_xlsx(results_univariate_genes, "plots/results_univariate_genes.xlsx")

gene_data <- preprocess_res$univariate_significant_genes

###### Checking correlations ######## 
clin_vars_corr <- c("Age" ,"ER_positive", "Lymphatic_invasion_present")
variables_both <- clinical_onehot %>% 
  column_to_rownames("Sample_ID") %>% 
  select(all_of(clin_vars_corr)) %>% 
  merge(gene_data, by=0) %>%
  column_to_rownames("Row.names")

corr_matrix <- matrix(nrow = ncol(variables_both), ncol = ncol(variables_both), dimnames = list(colnames(variables_both), colnames(variables_both)))
for (i in colnames(variables_both)) {
  for (j in colnames(variables_both)) {
  corr <- stats::cor(variables_both[,i], variables_both[,j])
  corr_matrix[j,i] <- corr
}
}

colnames(corr_matrix) <- colnames(corr_matrix) %>% str_replace_all(., "_", " ")
rownames(corr_matrix) <- rownames(corr_matrix) %>% str_replace_all(., "_", " ")

# Initialize file path
file_path= "plots/corr_plot_genclin.png"
png(height=1000, width=1000, file=file_path)

# Correlation plot
corr_plot <- corrplot(corr_matrix, 
                      order = "hclust", 
                      tl.col = "#E41A1CFF", 
                      tl.cex = 0.75, 
                      tl.srt = 45, 
                      type = "lower", 
                      diag = F)
colnames(corr_matrix)[1:3] <- ""
rownames(corr_matrix)[1:3] <- ""
corr_plot <- corrplot(corr_matrix, 
                      order = "hclust", 
                      tl.col = "black", 
                      tl.cex = 0.75, 
                      tl.srt = 45, 
                      type = "lower", 
                      diag = F, 
                      add = T)

# Then
dev.off()

## -------Model fitting ---------------------------------
### dataset splitting, alpha tuning (both and gene model), model fitting, 
# model prediction and C- index calculation (with a self made function model_evaluation() 
# (see helper.R for code)
Results_model_eval <- model_evaluation(x = clinical_onehot, 
                            y = gene_data, 
                            nrep = 100,                   # number of iterations
                            split.percentage = 0.8, 
                            method =  '1se', 
                            nfolds = 3,                 # number of folds
                            type.measure = "C", 
                            lambda = "lambda.min", 
                            pf = c(rep(0, length(clin_vars)), # Clinical variables are not penalized
                                rep(1, ncol(gene_data)))) 

save(Results_model_eval, file = "results_modelevaluation_2025_02_09.Rda") # saving results table 

### plotting obtained C-indexes for the three models.###
my_comparisons <- list( c("Cox.both", "Cox.clinical"), c("Cox.both", "Cox.genes"), c("Cox.clinical", "Cox.genes") )

# Overall C-indexes:
C_index_plot <- Results_model_eval %>% 
  filter(Metric_id == "Uno's C-index") %>% 
  mutate(across(Model_id, 
                ~factor(., levels=c("Cox.clinical", "Cox.genes", "Cox.both")))) %>% 
  ggplot(., aes(x = Model_id, y = value))+ 
  geom_boxplot(fill=c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF"))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1))+
  ylab("Uno's C-index")+
  xlab("")+
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  theme_pubr(base_size = 14)+
  scale_x_discrete(labels = c("Clinical","GEX" , "Clinical + GEX"))+
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(0.92, 1.0, 1.05), symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

#calculate medians for Uno's C-indexes
Results_model_eval %>% filter(Model_id=="Cox.clinical", Metric_id == "Uno's C-index") %>% pull(value) %>% median()
Results_model_eval %>% filter(Model_id=="Cox.genes", Metric_id == "Uno's C-index") %>% pull(value) %>% median()
Results_model_eval %>% filter(Model_id=="Cox.both", Metric_id == "Uno's C-index") %>% pull(value) %>% median()

## Time dependant C-index = ROC-AUC(t) ##
timepoints <- c("Time_dependant_AUC_6m" = "t = 6 months",
                "Time_dependant_AUC_1y" = "t = 1 year",
                "Time_dependant_AUC_2y" = "t = 2 years")

ROC_AUC_plot <- Results_model_eval %>% 
  filter(Metric_id %in% c("Time_dependant_AUC_6m", "Time_dependant_AUC_1y", "Time_dependant_AUC_2y")) %>% 
  mutate(across(Metric_id, 
                ~factor(., levels=c("Time_dependant_AUC_6m", "Time_dependant_AUC_1y", "Time_dependant_AUC_2y")))) %>%
  mutate(across(Model_id, 
                ~factor(., levels=c("Cox.clinical", "Cox.genes", "Cox.both")))) %>% 
  ggplot(., aes(x = Model_id, y = value, fill = Model_id))+ 
  geom_boxplot()+
  facet_wrap(~Metric_id, labeller = as_labeller(timepoints))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.25))+
  ylab("ROC-AUC(t)")+
  xlab("")+
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  theme_pubr(base_size = 14, legend = "bottom")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.background =element_rect(fill="white"))+
  scale_fill_manual(values = c("Cox.clinical" = "#E41A1CFF", "Cox.genes" = "#377EB8FF", "Cox.both" = "#4DAF4AFF"),     
                    labels = c("Cox.clinical" = "Clinical", "Cox.genes" = "GEX", "Cox.both" = "Clinical + GEX"), 
                    name = "Model") +  
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(0.99, 1.1, 1.15), label = "p.sign", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))


Results_model_eval %>% filter(Model_id=="Cox.clinical", Metric_id == "Time_dependant_AUC_6m") %>% pull(value) %>% median()
Results_model_eval %>% filter(Model_id=="Cox.genes", Metric_id == "Time_dependant_AUC_6m") %>% pull(value) %>% median()
Results_model_eval %>% filter(Model_id=="Cox.both", Metric_id == "Time_dependant_AUC_6m") %>% pull(value) %>% median()

ggsave(C_index_plot,file="plots/C_index_plot_2025_02_09.jpeg", height = 5, width = 5)
ggsave(ROC_AUC_plot,file="plots/ROC_AUC_plot_2025_02_09.jpeg", height = 4.5, width = 9.5)




