Interpretations
================
Mingzhou Fu
2024-08-19

``` r
rm(list = ls())
pacman::p_load(tidyverse, h2o, tidymodels, survival)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
'%!in%' = function(x,y)!('%in%'(x,y))
```

# Part 1. Shared features

## 1. Candidate SNPs

``` r
load(file = paste0(raw_data_path, "FUMA_outputs/overlap_fuma.rda"))
# load in csv file (there are some replication in the seed)
shared_features = read_csv(file = paste0(raw_data_path, "modeling/shared_features_elnet.csv"))
shared_features = shared_features %>% unique()
length(unique(shared_features$seed)) # 999
# read in a txt file with rsID
rsID = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/snps.txt"), 
                  header = T, sep = "\t", fill = T) %>% 
  select(uniqID, rsID) %>% unique()
```

``` r
shared_features %>% group_by(seed) %>% summarise(n = n()) %>% summary()
# median = 22, mean = 21.84
common_feature = shared_features %>% group_by(SNP) %>%
  summarise(n = n_distinct(seed)) %>% 
  left_join(overlap_fuma, by = c("SNP" = "chr_pos_name")) %>% 
  arrange(desc(n)) %>% filter(n >= 0.05*999) %>%
  select(SNP, n, rsID, phenotype, chr_num, pos, CADD, nearestGene, func, 
         posMapFilt, eqtlMapFilt, ciMapFilt) %>% unique()
dim(common_feature) # dim = (83,12)
# get coefficients
common_feature_coef = shared_features %>% 
  filter(SNP %in% common_feature$SNP) %>% group_by(SNP) %>% 
  summarise(mean_AD_coef = mean(AD_coef), 
            mean_LOE_coef = mean(LOE_coef)) %>% unique() %>% 
  mutate(signs = ifelse(mean_AD_coef * mean_LOE_coef > 0, "Same", "Different")) 
save(common_feature_coef, file = paste0(raw_data_path, "modeling/common_feature_coef.rda"))
# combine results
common_feature_final = common_feature %>% 
  left_join(common_feature_coef, by = "SNP")
dim(common_feature_final) # dim = (83,15)
# outputs
write.table(common_feature_final,
            file = paste0(output_path, "rsID_shared.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

## 2. Map to genes

``` r
AD_genes = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/genes.txt"), 
                      header = T, sep = "\t", fill = T)
pos_AD = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/snps.txt"), 
                    header = T, sep = "\t", fill = T) %>% 
  filter(posMapFilt == 1) %>% select(uniqID, rsID, nearestGene) %>% unique()
eqtl_AD = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/eqtl.txt"), 
                     header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique() %>% 
  mutate(chr_pos = str_extract(uniqID, "^\\d+:\\d+"))
AD_CI = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/ci.txt"), 
                   header = T, sep = "\t", fill = T) %>% select(SNPs, genes) %>% unique()

LOE_genes = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_GGE/genes.txt"), 
                       header = T, sep = "\t", fill = T)
pos_LOE = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_GGE/snps.txt"),
                     header = T, sep = "\t", fill = T) %>% 
  filter(posMapFilt == 1) %>% select(uniqID, rsID, nearestGene) %>% unique()
eqtl_LOE = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_GGE/eqtl.txt"),
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique() %>% 
  mutate(chr_pos = str_extract(uniqID, "^\\d+:\\d+"))
LOE_CI = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_GGE/ci.txt"), 
                    header = T, sep = "\t", fill = T) %>% select(SNPs, genes) %>% unique()
```

``` r
AD_SNPs_df = common_feature_final %>% filter(phenotype == "AD_EUR") %>%
  select(rsID, signs) %>% unique() %>% left_join(overlap_fuma) %>% 
  select(rsID, IndSigSNP, signs, posMapFilt, eqtlMapFilt, ciMapFilt) %>% 
  mutate(direct = if_else(rsID == IndSigSNP, 1, 0)) %>% unique()
dim(AD_SNPs_df) # dim = (69, 7)
LOE_SNPs_df = common_feature_final %>% filter(phenotype == "Epilepsy") %>%
  select(rsID, signs) %>% unique() %>% left_join(overlap_fuma) %>% 
  select(rsID, IndSigSNP, signs, posMapFilt, eqtlMapFilt, ciMapFilt) %>% 
  mutate(direct = if_else(rsID == IndSigSNP, 1, 0)) %>% unique()
dim(LOE_SNPs_df) # dim = (14, 7)
```

### 1) Same direction set

``` r
AD_SNPs_direct = AD_SNPs_df %>% filter(signs == "Same") %>%
  filter(direct == 1) %>% pull(rsID) %>% unique() 
length(AD_SNPs_direct) # 8
AD_genes_direct_mapped = AD_genes %>% 
  filter(str_detect(IndSigSNPs, str_c(AD_SNPs_direct, collapse = "|"))) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(AD_genes_direct_mapped) # dim = (30,11)

AD_SNP_indirect = AD_SNPs_df %>% filter(signs == "Same") %>% filter(direct == 0) 
# position only, no eqtl or ci
AD_indirect_genes = pos_AD %>% filter(rsID %in% AD_SNP_indirect$rsID) %>% 
  pull(nearestGene) %>% unique()
AD_genes_indirect_mapped = AD_genes %>% 
  filter(symbol %in% AD_indirect_genes) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(AD_genes_indirect_mapped) # dim = (7,11)
# Combine all
AD_genes_mapped = rbind(AD_genes_direct_mapped, AD_genes_indirect_mapped) %>% 
  unique() %>% mutate(phenotype = "AD")
dim(AD_genes_mapped) # dim = (34,12)
```

``` r
LOE_SNPs_direct = LOE_SNPs_df %>% filter(signs == "Same") %>%
  filter(direct == 1) %>% pull(rsID) %>% unique() 
length(LOE_SNPs_direct) # 0
LOE_SNP_indirect = LOE_SNPs_df %>% filter(signs == "Same") %>% filter(direct == 0) 
# no eqtl or ci
pos_SNP = LOE_SNP_indirect %>% filter(posMapFilt == 1)
LOE_indirect_genes = pos_LOE %>% filter(rsID %in% pos_SNP$rsID) %>% 
  pull(nearestGene) %>% unique()
LOE_genes_indirect_mapped = LOE_genes %>% 
  filter(symbol %in% LOE_indirect_genes) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(LOE_genes_indirect_mapped) # dim = (2,11)
# Combine all
LOE_genes_mapped = LOE_genes_indirect_mapped %>% mutate(phenotype = "LOE")
dim(LOE_genes_mapped) # dim = (2,12)
```

``` r
all_genes_mapped_same = rbind(AD_genes_mapped, LOE_genes_mapped) %>% unique()
dim(all_genes_mapped_same) # dim = (36,12)
# output
write.table(all_genes_mapped_same,
            file = paste0(output_path, "all_genes_mapped_same.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

### 2) Opposite direction set

``` r
AD_SNPs_direct = AD_SNPs_df %>% filter(signs == "Different") %>%
  filter(direct == 1) %>% pull(rsID) %>% unique() 
length(AD_SNPs_direct) # 9
AD_genes_direct_mapped = AD_genes %>% 
  filter(str_detect(IndSigSNPs, str_c(AD_SNPs_direct, collapse = "|"))) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(AD_genes_direct_mapped) # dim = (68,11)

AD_SNP_indirect = AD_SNPs_df %>% filter(signs == "Different") %>% filter(direct == 0) 
# positional
AD_indirect_pos = AD_SNP_indirect %>% filter(posMapFilt == 1) %>%
  select(rsID) %>% unique() %>% left_join(pos_AD) %>%
  pull(nearestGene) %>% unique()
# eQTL
eqtl_SNP = AD_SNP_indirect %>% filter(eqtlMapFilt == 1) %>% pull(rsID) %>% unique()
AD_indirect_eqtl = pos_AD %>% filter(rsID %in% eqtl_SNP) %>% 
  select(uniqID, rsID) %>% unique() %>% left_join(eqtl_AD) %>% 
  pull(symbol) %>% unique()
# Chrom interaction
ci_SNP = AD_SNP_indirect %>% filter(ciMapFilt == 1) %>% pull(rsID) %>% unique()
AD_indirect_ci = AD_CI %>% 
  filter(str_detect(SNPs, str_c(ci_SNP, collapse = "|"))) %>% drop_na() %>%
  unique() %>% left_join(AD_genes, by = c("genes" = "ensg")) %>% 
  filter(ciMap == "Yes") %>% pull(symbol) %>% unique()
AD_genes_indirect_mapped = AD_genes %>% 
  filter(symbol %in% unique(c(AD_indirect_pos, AD_indirect_eqtl, AD_indirect_ci))) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(AD_genes_indirect_mapped) # dim = (46,11)

# Combine all
AD_genes_mapped = rbind(AD_genes_direct_mapped, AD_genes_indirect_mapped) %>% 
  unique() %>% mutate(phenotype = "AD")
dim(AD_genes_mapped) # dim = (101,12)
```

``` r
LOE_SNPs_direct = LOE_SNPs_df %>% filter(signs == "Different") %>%
  filter(direct == 1) %>% pull(rsID) %>% unique() 
length(LOE_SNPs_direct) # 0
LOE_SNP_indirect = LOE_SNPs_df %>% filter(signs == "Different") %>% filter(direct == 0) 
# positional
LOE_indirect_pos = LOE_SNP_indirect %>% filter(posMapFilt == 1) %>%
  select(rsID) %>% unique() %>% left_join(pos_LOE) %>%
  pull(nearestGene) %>% unique()
# Chrom interaction
ci_SNP = LOE_SNP_indirect %>% filter(ciMapFilt == 1) %>% pull(rsID) %>% unique()
LOE_indirect_ci = LOE_CI %>% 
  filter(str_detect(SNPs, str_c(ci_SNP, collapse = "|"))) %>% drop_na() %>%
  unique() %>% left_join(LOE_genes, by = c("genes" = "ensg")) %>% 
  filter(ciMap == "Yes") %>% pull(symbol) %>% unique()
LOE_genes_indirect_mapped = LOE_genes %>% 
  filter(symbol %in% unique(c(LOE_indirect_pos, LOE_indirect_ci))) %>% 
  select(ensg, symbol, chr, start, end, type, posMapSNPs, posMapMaxCADD,
         eqtlMapSNPs, ciMap, minGwasP) %>% unique() 
dim(LOE_genes_indirect_mapped) # dim = (3,11)
# no mapped genes found for LOE
```

``` r
all_genes_mapped_diff = rbind(AD_genes_mapped, LOE_genes_mapped) %>% unique()
dim(all_genes_mapped_diff) # dim = (103,12)
# output
write.table(all_genes_mapped_diff,
            file = paste0(output_path, "all_genes_mapped_different.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

# Part 2. Genetic association tests

``` r
phe_table_extract_date = "2024-08-13"
# Load in datasets
load(file = paste0(raw_data_path, "atlas/pheno/final/sample_demo_eligible_", phe_table_extract_date, ".rda"))
load(file = paste0(raw_data_path, "atlas/pheno/final/sample_enc_eligible_", phe_table_extract_date, ".rda"))
```

## 1. Shared genetic risk score

``` r
load(file = paste0(raw_data_path, "modeling/sample_data_full.rda"))
load(file = paste0(raw_data_path, "modeling/genotype_df_filter.rda"))
load(file = paste0(raw_data_path, "modeling/common_feature_coef.rda"))
```

### 1) Retrain weights for shared SNPs

See details `4_MultitaskElasticNet.ipynb`.

``` r
# read in csv
shared_SNP_coef = read_csv(paste0(raw_data_path, "modeling/shared_coef_retrain.csv"))
names(shared_SNP_coef) = c("SNP", "new_coef")
shared_genotype = genotype_df_filter %>% 
  select(IID, any_of(shared_SNP_coef$SNP)) %>% unique()
dim(shared_genotype) # dim = (34463, 28)
```

### 2) Calculate shared-PRS

``` r
# Filter genotype info
genotype.targ = shared_genotype %>% select(all_of(shared_SNP_coef$SNP))
genotype.targ.matrix = as.matrix(genotype.targ[])
# Calculate PRS
pred_prs = genotype.targ.matrix %*% shared_SNP_coef$new_coef
pred_prs_final = cbind(pred_prs, shared_genotype[, 'IID']) %>% as.data.frame()
names(pred_prs_final) = c('AD_LOE_shared_PRS', 'UniqueSampleId')
mean_X = mean(pred_prs_final$AD_LOE_shared_PRS) # -0.01196118
sd_X = sd(pred_prs_final$AD_LOE_shared_PRS) # 0.08250386
# Normalize PRS
pred_prs_final = pred_prs_final %>% 
  mutate(AD_LOE_shared_PRS_normalized = (AD_LOE_shared_PRS - mean_X) / sd_X)
dim(pred_prs_final) # dim = (34463, 3)
save(pred_prs_final, file = paste0(raw_data_path, "modeling/AD_LOE_shared_PRS.rda"))
```

``` r
# Join sample demo
patient_genetic = sample_demo_eligible %>% 
  left_join(pred_prs_final) %>% select(-AD_LOE_shared_PRS)
dim(patient_genetic) # dim = (2179, 28)
```

## 2. Test for associations

### 1) AD \| LOE ~ AD_LOE_shared_PRS

``` r
# build models
glm_AD = linear_reg() %>% 
  set_engine("glm") %>% 
  set_mode("regression") %>% 
  fit(AD ~ AD_LOE_shared_PRS_normalized + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
      data = patient_genetic)
glm_LOE = linear_reg() %>% 
  set_engine("glm") %>% 
  set_mode("regression") %>% 
  fit(LOE ~ AD_LOE_shared_PRS_normalized + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
      data = patient_genetic)
# extract results
sum_res = rbind(tidy(glm_AD, conf.int = TRUE, exponentiate = TRUE)[2,],
                tidy(glm_LOE, conf.int = TRUE, exponentiate = TRUE)[2,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(OR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(OR, p_value)
sum_res
```

``` r
# Fit the full logistic regression model
full_model = glm(AD ~ AD_LOE_shared_PRS_normalized + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Fit the reduced model excluding the variable of interest (e.g., x2)
reduced_model = glm(AD ~ female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Calculate the deviance for both models
deviance_full = deviance(full_model)
deviance_reduced = deviance(reduced_model)
# Calculate the percentage of variation explained by the variable x2
percent_variation_explained = 100 * (deviance_reduced - deviance_full) / deviance_reduced
percent_variation_explained
```

``` r
# Fit the full logistic regression model
full_model = glm(LOE ~ AD_LOE_shared_PRS_normalized + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Fit the reduced model excluding the variable of interest (e.g., x2)
reduced_model = glm(LOE ~ female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Calculate the deviance for both models
deviance_full = deviance(full_model)
deviance_reduced = deviance(reduced_model)
# Calculate the percentage of variation explained by the variable x2
percent_variation_explained = 100 * (deviance_reduced - deviance_full) / deviance_reduced
percent_variation_explained
```

### 2) AD \| LOE ~ e4count (comparison)

``` r
glm_AD = linear_reg() %>% 
  set_engine("glm") %>% 
  set_mode("regression") %>% 
  fit(AD ~ e4count + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
      data = patient_genetic)
glm_LOE = linear_reg() %>% 
  set_engine("glm") %>% 
  set_mode("regression") %>% 
  fit(LOE ~ e4count + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
      data = patient_genetic)
# extract results
sum_res = rbind(tidy(glm_AD, conf.int = TRUE, exponentiate = TRUE)[2,],
                tidy(glm_LOE, conf.int = TRUE, exponentiate = TRUE)[2,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(OR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(OR, p_value)
sum_res
```

``` r
# Fit the full logistic regression model
full_model = glm(AD ~ e4count + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Fit the reduced model excluding the variable of interest (e.g., x2)
reduced_model = glm(AD ~ female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Calculate the deviance for both models
deviance_full = deviance(full_model)
deviance_reduced = deviance(reduced_model)
# Calculate the percentage of variation explained by the variable x2
percent_variation_explained = 100 * (deviance_reduced - deviance_full) / deviance_reduced
percent_variation_explained
```

``` r
# Fit the full logistic regression model
full_model = glm(LOE ~ e4count + female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Fit the reduced model excluding the variable of interest (e.g., x2)
reduced_model = glm(LOE ~ female + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        data = patient_genetic, family = binomial)
# Calculate the deviance for both models
deviance_full = deviance(full_model)
deviance_reduced = deviance(reduced_model)
# Calculate the percentage of variation explained by the variable x2
percent_variation_explained = 100 * (deviance_reduced - deviance_full) / deviance_reduced
percent_variation_explained
```
