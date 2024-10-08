Genotype data cleaning
================
Mingzhou Fu
2024-08-13

``` r
rm(list = ls())
pacman::p_load(tidyverse, GenomicRanges, rtracklayer, bigsnpr, bigreadr, 
               data.table, ieugwasr, R.utils)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
'%!in%' = function(x,y)!('%in%'(x,y))
```

# Part 0. GWAS preprocessing

Preprocessing of GWAS summary statistics, ran once

## 1. AD EUR

``` r
sumstats_file = paste0(raw_data_path, 
                       "sumstats/raw/Wightman2021_AD_EUR.txt.gz")
sumstats = fread(sumstats_file)
names(sumstats) = c("chr", "pos", "a1", "a0", "z", "p", "N")
sumstats_clean = sumstats %>% filter(!is.na(chr) & !is.na(pos)) %>% 
  mutate(p = as.numeric(p))
dim(sumstats_clean) # dim = (12688339,7)
write.table(sumstats_clean, 
            file = paste0(raw_data_path, 
                          'sumstats/fuma_upload/AD_Wightman2021_EUR_hg19.txt'),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

``` r
AD_sumstats_clean = sumstats_clean
# read in reference file
ref_file = paste0(raw_data_path, 'sumstats/raw/hm3_151_match.txt')
hm3_ref = fread(ref_file)
names(hm3_ref) = c("CHR", "BP", "SNP")
hm3_ref = hm3_ref %>% 
  mutate(chr = as.numeric(substring(CHR, 4, nchar(CHR)))) %>% 
  dplyr::rename("pos" = "BP") %>% select(chr, pos, SNP) %>% unique()
dim(hm3_ref) # dim = (1473450,3)
# prepare for ldsc
sumstats_ldsc = AD_sumstats_clean %>% 
  left_join(hm3_ref, by = c("chr" = "chr", "pos" = "pos")) %>%
  filter(!is.na(SNP)) %>% select(SNP, N, z, a1, a0, p) %>% unique()
names(sumstats_ldsc) = c("SNP", "N", "Z", "A1", "A2", "P")
dim(sumstats_ldsc) # dim = (1236234,6)
write.table(sumstats_ldsc, 
            file = paste0(raw_data_path, 
                          'sumstats/ldsc/AD_Wightman2021_EUR_hg19.txt'),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

## 2. ILAEC EUR (European descent (92%))

``` r
sumstats_file = paste0(raw_data_path, 
                       "sumstats/raw/ILAEC2023_GGE_EUR.tsv.gz")
sumstats = fread(sumstats_file)
names(sumstats) = c("chr", "pos", "a1", "a0", "beta", "se", "eaf", "p", "rsid")
sumstats_clean = sumstats %>% filter(!is.na(chr) & !is.na(pos)) %>% 
  mutate(p = as.numeric(p), a1 = toupper(a1), a0 = toupper(a0))
dim(sumstats_clean) # dim = (4860774,9)
write.table(sumstats_clean, 
            file = paste0(raw_data_path, 
                          'sumstats/fuma_upload/GGE_ILAEC2023_EUR_hg19.txt'),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

``` r
LOE_sumstats_clean = sumstats_clean
sumstats_ldsc = LOE_sumstats_clean %>% 
  mutate(N = 49843) %>% 
  mutate(z = abs(qnorm(1 - p / 2)) * sign(beta)) %>% 
  select(rsid, N, z, a1, a0, p) %>% unique()
names(sumstats_ldsc) = c("SNP", "N", "Z", "A1", "A2", "P")
dim(sumstats_ldsc) # dim = (4860774,6)
write.table(sumstats_ldsc, 
            file = paste0(raw_data_path, 
                          'sumstats/ldsc/GGE_ILAEC2023_EUR_hg19.txt'),
            sep = "\t", quote = F, row.names = F, col.names = T)
```

# Part 1. FUMA GWAS results

## 1. Clean FUMA results

### 1) AD Kunkle EUR

``` r
# Full list of independent significant SNPs + LD SNPs (candidate SNPs, catch as many SNPs as we can from ATLAS)
full_AD_EUR = read.table(file = paste0(raw_data_path, 
                                       "FUMA_outputs/FUMA_AD_EUR/snps.txt"), 
                         header = T, sep = "\t", fill = T)
lead_AD_EUR = read.table(file = paste0(raw_data_path, 
                                       "FUMA_outputs/FUMA_AD_EUR/leadsnps.txt"), 
                         header = T, sep = "\t", fill = T)
candi_AD_EUR = full_AD_EUR %>% filter(!is.na(gwasP)) %>% 
  mutate(lead = if_else(rsID %in% lead_AD_EUR$rsID, 1, 0)) %>% 
  mutate(phenotype = "AD_EUR") %>% 
  select(chr, pos, rsID, non_effect_allele, effect_allele, gwasP, IndSigSNP, 
         nearestGene, func, CADD, posMapFilt, eqtlMapFilt, ciMapFilt, lead, 
         phenotype) %>% unique()
dim(candi_AD_EUR) # dim = (4154,15)
```

### 2) ILAEC EUR (European descent (92%))

``` r
# Full list of independent significant SNPs + LD SNPs (candidate SNPs, catch as many SNPs as we can from ATLAS)
full_Epilepsy_EUR = read.table(file = paste0(raw_data_path, 
                                             "FUMA_outputs/FUMA_GGE/snps.txt"), 
                         header = T, sep = "\t", fill = T)
lead_Epilepsy_EUR = read.table(file = paste0(raw_data_path, 
                                             "FUMA_outputs/FUMA_GGE/leadsnps.txt"), 
                         header = T, sep = "\t", fill = T)
candi_Epilepsy_EUR = full_Epilepsy_EUR %>% filter(!is.na(gwasP)) %>% 
  mutate(lead = if_else(rsID %in% lead_Epilepsy_EUR$rsID, 1, 0)) %>% 
  mutate(phenotype = "Epilepsy") %>% 
  select(chr, pos, rsID, non_effect_allele, effect_allele, gwasP, IndSigSNP, 
         nearestGene, func, CADD, posMapFilt, eqtlMapFilt, ciMapFilt, lead, 
         phenotype) %>% unique()
dim(candi_Epilepsy_EUR) # dim = (2104,15)
```

## 2. Liftover + annotation

``` r
candi_full = rbind(candi_AD_EUR, candi_Epilepsy_EUR) %>% as.data.frame() %>% 
  mutate(index = seq(1, nrow(.)))
# liftover
chain_file = import.chain("hg19ToHg38.over.chain")
# Prepare Genomic Ranges
hg19.gr = GRanges(
  seqname = paste("chr", candi_full$chr, sep = ""),
  ranges = IRanges(start = as.numeric(candi_full$pos), 
                   end = as.numeric(candi_full$pos)),
  snp.name = candi_full$rsID)
candi_hg38_liftover = as.data.frame(liftOver(hg19.gr, chain_file)) %>% 
  mutate(chr = as.integer(substring(seqnames, 4, last = 1000000L))) %>% 
  select(group, chr, start) %>% unique() %>% drop_na()
liftover_full = candi_full %>% 
  left_join(candi_hg38_liftover, by = c("index" = "group", "chr" = "chr")) %>% 
  select(chr, start, rsID, non_effect_allele, effect_allele, gwasP, IndSigSNP, 
         nearestGene, func, CADD, posMapFilt, eqtlMapFilt, ciMapFilt, lead, phenotype) %>% 
  dplyr::rename(pos = start, a0 = non_effect_allele, a1 = effect_allele) %>% 
  filter(!is.na(chr) & !is.na(pos)) %>% unique()
dim(liftover_full) # dim = 6240,15
```

## 3. Extract and clean ATLAS data – Server step

``` r
# Check overlap with FUMA
liftover_1 = liftover_full %>% dplyr::rename("chr_num" = "chr") %>% 
  mutate(chr_pos_name1 = paste0("chr", chr_num, ":", pos, ":", a0, ":", a1)) %>% 
  dplyr::rename("chr_pos_name" = "chr_pos_name1")
liftover_2 = liftover_full %>% dplyr::rename("chr_num" = "chr") %>% 
  mutate(chr_pos_name2 = paste0("chr", chr_num, ":", pos, ":", a1, ":", a0)) %>% 
  dplyr::rename("chr_pos_name" = "chr_pos_name2")
overlap_fuma = rbind(liftover_1, liftover_2) %>% unique()
dim(overlap_fuma) # dim = 12480,16
# Save extract_fuma
extract_fuma = overlap_fuma %>% pull(chr_pos_name) %>% unique()
length(extract_fuma) # length = 12480
save(overlap_fuma, file = paste0(raw_data_path, "FUMA_outputs/overlap_fuma.rda"))
# Write output for snp extraction on server
write.table(extract_fuma,
            file = paste0(raw_data_path, "FUMA_outputs/extract_fuma.txt"),
            sep = "\t", quote = F, row.names = F, col.names = F)
```

``` r
# Read in the QCed genotype and sumstats files + preprocess the bed file (only need to do once for each data set)
geno_file = paste0(raw_data_path, 'atlas/geno/fuma.raw')
snp_readBed(paste0(geno_file, '.bed'))
obj.bigSNP = snp_attach(paste0(geno_file, '.rds'))
# extract the SNP information from the genotype
map = obj.bigSNP$map[-3]
names(map) = c("chr", "rsid", "pos", "a1", "a0")
# Rename rsid to keep consistent with bim file
map = map %>% mutate(RSID = paste0("chr", chr, ":", pos, ":", a0, ":", a1)) %>% 
  select(chr, RSID, pos, a1, a0) %>% dplyr::rename("rsid" = "RSID")
# Get fam order
fam.order = as.data.table(obj.bigSNP$fam)
# get genotype info
genotype = obj.bigSNP$genotypes
genotype_df = as.matrix(genotype[])[,] %>% as.data.frame()
colnames(genotype_df) = map$rsid
```

``` r
# The above SNPs have a missing, we decide to exclude those SNPs (N = 1)
missing_info = genotype_df %>% 
  summarise_all(function(x) sum(is.na(x))) %>% gather(variable, missing_count) 
columns_with_missing = missing_info %>% 
  filter(missing_count > 0) %>% pull(variable)
length(columns_with_missing) # 194
genotype_df_noNA = genotype_df %>% select(-all_of(columns_with_missing)) 
dim(genotype_df_noNA) # 34463,5639
genotype_df_fam = cbind(fam.order[,1:2], genotype_df_noNA) %>% 
  as.data.frame() %>% drop_na()
colnames(genotype_df_fam)[1:2] = c("FID", "IID")
dim(genotype_df_fam) # 34463,5641
# also need to update map file
map_noNA = map %>% filter(rsid %!in% columns_with_missing)
dim(map_noNA) # 5639,5
save(map_noNA, file = paste0(raw_data_path, "modeling/map_noNA_full_SNP.rda"))
save(genotype_df_fam, file = paste0(raw_data_path, 
                                    "modeling/freeze60k_geno_freq_full_SNP.rda"))
```

# Part 2. Prioritization of SNPs

## 1. Select mapped SNPs

``` r
candidate_snp_info = overlap_fuma %>% 
  select(chr_num, pos, rsID, gwasP, IndSigSNP, func, CADD, lead, phenotype, 
         nearestGene, posMapFilt, eqtlMapFilt, ciMapFilt) %>%
  mutate(chr_num = as.numeric(chr_num)) %>%
  inner_join(map_noNA, by = c("chr_num" = "chr", "pos" = "pos")) %>%
  mutate(mapping_num = posMapFilt + eqtlMapFilt + ciMapFilt) %>% 
  select(rsid, rsID, chr_num, pos, gwasP, IndSigSNP, func, CADD, mapping_num, 
         nearestGene, phenotype) %>% filter(mapping_num > 0 | CADD >= 12.37) %>% 
  unique()
dim(candidate_snp_info) # dim = 4454,11
save(candidate_snp_info, file = paste0(raw_data_path, "modeling/candidate_snp_info.rda"))
genotype_df_filter = genotype_df_fam %>% select(IID, all_of(candidate_snp_info$rsid))
dim(genotype_df_filter) # dim = 34463,4455
save(genotype_df_filter, file = paste0(raw_data_path, "modeling/genotype_df_filter.rda"))
```

## 2. Finalize modeling dataframe

``` r
phe_table_extract_date = "2024-08-13"
load(file = paste0(raw_data_path, "atlas/pheno/final/sample_demo_eligible_", 
                   phe_table_extract_date, ".rda"))
load(file = paste0(raw_data_path, "modeling/genotype_df_filter.rda"))
```

``` r
# make a full dataset
sample_data_full = sample_demo_eligible %>% 
  select(UniqueSampleId, female, age_last_visit, e4count, 
         AD_PRS_full_norm, AD_PRS_rmapoe_norm, AD, LOE, PC1:PC10) %>%
  inner_join(genotype_df_filter, by = c("UniqueSampleId" = "IID")) %>%
  select(-UniqueSampleId) %>% as.data.frame()
dim(sample_data_full) # dim = 2179,4471
save(sample_data_full, file = paste0(raw_data_path, "modeling/sample_data_full.rda"))
# # output to csv
write.table(sample_data_full,
            file = paste0(raw_data_path, "modeling/sample_data_full.csv"),
            sep = ",", quote = F, row.names = F)
```
