Supplementary and senstivity analysis
================
Mingzhou Fu
2024-08-31

``` r
rm(list = ls())
pacman::p_load(tidyverse, gtsummary, PheWAS)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
```

# Part 1. Model evaluation

``` r
model_evaluation = read_csv(file = paste0(raw_data_path, "modeling/model_eval_elnet.csv"))
model_evaluation = model_evaluation %>% unique()
length(unique(model_evaluation$seed)) # 999
```

## 1. AD

``` r
aucpr_AD = model_evaluation %>%  
  select(model, seed, auprc_ad) %>% drop_na() %>% unique() %>%
  pivot_wider(names_from = model, values_from = auprc_ad)

wilcox_AD_aucpr = wilcox.test(aucpr_AD$`Multi-Task Elastic Net`, 
                              aucpr_AD$`Separate Elastic Net (AD)`, 
                              paired = TRUE, alternative = "greater")
wilcox_AD_aucpr
```

## 2. LOE

``` r
aucpr_LOE = model_evaluation %>%  
  select(model, seed, auprc_loe) %>% drop_na() %>% unique() %>%
  pivot_wider(names_from = model, values_from = auprc_loe)
wilcox_LOE_aucpr = wilcox.test(aucpr_LOE$`Multi-Task Elastic Net`, 
                              aucpr_LOE$`Separate Elastic Net (LOE)`, 
                              paired = TRUE, alternative = "greater")
wilcox_LOE_aucpr
```

``` r
aucpr_AD_sig = aucpr_AD %>% 
  mutate(diff_AD = `Multi-Task Elastic Net` - `Separate Elastic Net (AD)`) %>% 
  filter(diff_AD > 0) %>% select(seed, diff_AD)
aucpr_LOE_sig = aucpr_LOE %>%
  mutate(diff_LOE = `Multi-Task Elastic Net` - `Separate Elastic Net (LOE)`) %>% 
  filter(diff_LOE > 0) %>% select(seed, diff_LOE)
aucpr_sig_final = aucpr_AD_sig %>% inner_join(aucpr_LOE_sig) %>% 
  arrange(desc(diff_LOE)) %>% select(seed, diff_AD, diff_LOE)
```

# Part 2. LDSC

``` bash
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
conda activate ldsc
```

``` bash
python munge_sumstats.py \
    --sumstats /Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/sumstats/ldsc/AD_Wightman2021_EUR_hg19.txt \
    --merge-alleles 1000G_Phase3_baselineLDscore/w_hm3.snplist \
    --out formatted_AD_sumstats \
    --a1 A1 --a2 A2 --snp SNP --N-col N --p P 
    
# Writing summary statistics for 1217311 SNPs (1040320 with nonmissing beta) to formatted_AD_sumstats.sumstats.gz.

# Metadata:
# Mean chi^2 = 1.18
# Lambda GC = 1.098
# Max chi^2 = 1226.904
# 405 Genome-wide significant SNPs (some may have been removed by filtering).
```

``` bash
python munge_sumstats.py \
    --sumstats /Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/sumstats/ldsc/GGE_ILAEC2023_EUR_hg19.txt \
    --merge-alleles 1000G_Phase3_baselineLDscore/w_hm3.snplist \
    --out formatted_GGE_sumstats \
    --a1 A1 --a2 A2 --snp SNP --N-col N --p P 
    
# Writing summary statistics for 1217311 SNPs (903482 with nonmissing beta) to formatted_GGE_sumstats.sumstats.gz.

# Metadata:
# Mean chi^2 = 1.35
# Lambda GC = 1.293
# Max chi^2 = 78.289
# 299 Genome-wide significant SNPs (some may have been removed by filtering).
```

``` bash
python ldsc.py \
    --rg formatted_AD_sumstats.sumstats.gz,formatted_GGE_sumstats.sumstats.gz \
    --ref-ld-chr 1000G_Phase3_baselineLDscore/ \
    --w-ld-chr 1000G_Phase3_baselineLDscore/ \
    --out gencor_AD_GGE
```

# Part 3. Supplementary tables

## 1. Phecodes & ICD-10

``` r
epi_phecode = c("345", "345.1", "345.11", "345.12", "345.3")
icd_short = icd.data::icd10cm2016 %>% as.data.frame() %>% 
  mutate(code = as.character(code)) %>% 
  select(code, long_desc) %>% dplyr::rename("code_nodigit" = "code")
epi_tbl = PheWAS::phecode_map %>% 
  filter(vocabulary_id == "ICD10CM" & phecode %in% epi_phecode) %>% 
  mutate(code_nodigit = str_replace(code, "\\.", "")) %>% 
  left_join(icd_short, by = "code_nodigit") %>% 
  select(phecode, code, long_desc) %>% unique() %>% drop_na() %>% 
  mutate(long_desc = gsub("<0xa0>", " ", long_desc))
write_delim(epi_tbl, delim = "\t",
            file = paste0(output_path, "epi_phecode_tbl.txt"), quote = "none")
```

``` r
ad_tbl = PheWAS::phecode_map %>% 
  filter(vocabulary_id == "ICD10CM" & phecode == "290.11") %>% 
  mutate(code_nodigit = str_replace(code, "\\.", "")) %>% 
  left_join(icd_short, by = "code_nodigit") %>% 
  select(phecode, code, long_desc) %>% unique() %>% drop_na() %>% 
  mutate(long_desc = gsub("<0xa0>", " ", long_desc))
write_delim(ad_tbl, delim = "\t",
            file = paste0(output_path, "ad_phecode_tbl.txt"), quote = "none")
```

## 2. SNP coefficients figure

``` r
load(file = paste0(raw_data_path, "modeling/common_feature_coef.rda"))
shared_SNP_coef = common_feature_coef %>% filter(signs == "Same") %>%
  mutate(mean_coef = (mean_AD_coef + mean_LOE_coef) / 2) %>%
  select(SNP, mean_coef) %>% unique()
dim(shared_SNP_coef) # dim = (27,2)
# output to table
write_delim(shared_SNP_coef, delim = "\t",
            file = paste0(output_path, "shared_SNP_coef.txt"), quote = "none")
```

``` r
# plot figure
shared_SNP_coef_plt = common_feature_coef %>% filter(signs == "Same") %>% 
  mutate(chr = as.numeric(substr(str_extract(SNP, "chr[0-9]+"), 4, 1000000L)),
         pos = as.numeric(str_extract(SNP, "(?<=:)[0-9]+(?=:)"))) %>%
  arrange(desc(chr), desc(pos)) %>% dplyr::rename(AD = mean_AD_coef, LOE = mean_LOE_coef) %>%
  select(SNP, AD, LOE) %>% 
  pivot_longer(cols = c(AD, LOE), names_to = "Phenotype", 
               values_to = "Coefficient")
shared_SNP_coef_plt_add = shared_SNP_coef_plt %>% 
  mutate(SNP = factor(SNP, levels = unique(SNP)))

pdf(paste0(output_path, "shared_coef_plot.pdf"), width = 8, height = 6)
ggplot(shared_SNP_coef_plt_add, aes(x = SNP, y = Coefficient, fill = Coefficient > 0)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#ED706B", "FALSE" = "#4D4DF5")) +
  facet_wrap(~Phenotype, nrow = 1) + coord_cartesian(ylim = c(-0.02, 0.02)) + 
  coord_flip() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()
```
