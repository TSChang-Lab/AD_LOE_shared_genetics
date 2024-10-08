Phenotype data cleaning
================
Mingzhou Fu
2024-08-13

# Part 1. Raw EHR data preprocessing (ATLAS patients)

``` r
rm(list = ls())
pacman::p_load(tidyverse, PheWAS)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
'%!in%' = function(x,y)!('%in%'(x,y))
```

## 1. Clean encounters info

``` r
phe_table_extract_date = "2024-08-13"
# Read in ATLAS patient data
atlas_full_raw = read_csv(file = paste0(raw_data_path, 
                                        "atlas/pheno/raw/atlas_full_raw_0813.csv"),
                          col_names = FALSE, na = c("NULL"))
names(atlas_full_raw) = c("PatientID", "Sex", "Ethnicity", "FirstRace", 
                          "BirthDate", "DeathDate", "PatientLivingStatus",
                          "DateFirst", "DateLast", "BankSampleID", 
                          "DiagnosisCode", "DiagnosisDate")
dim(atlas_full_raw) # dim = (18341456,12)
# Clean encounter info
atlas_full = atlas_full_raw %>% 
  mutate(BirthDate = as.Date(BirthDate), DeathDate = as.Date(DeathDate),
         DateFirst = as.Date(DateFirst), DateLast = as.Date(DateLast),
         DiagnosisDate = as.Date(DiagnosisDate)) %>% 
  mutate(valid_record = case_when(
    is.na(DeathDate) & DiagnosisDate >= BirthDate ~ 1,
    !is.na(DeathDate) & DiagnosisDate >= BirthDate & DiagnosisDate <= DeathDate ~ 1
  )) %>% filter(valid_record == 1) %>% select(-valid_record) %>% 
  filter(Sex != "X" & Sex != "Unknown")
dim(atlas_full) # dim = (18340024,12)
```

## 2. Clean diagnosis info

``` r
atlas_enc_phecode = atlas_full %>% 
  group_by(PatientID, DiagnosisCode) %>% 
  arrange(DiagnosisDate) %>% slice_head(n = 1) %>% ungroup() %>% 
  mutate(vocabulary_id = "ICD10CM") %>% 
  left_join(PheWAS::phecode_map, by = c("vocabulary_id" = "vocabulary_id",
                                        "DiagnosisCode" = "code")) %>% unique() %>% 
  group_by(PatientID, phecode) %>% 
  arrange(DiagnosisDate) %>% slice_head(n = 1) %>% ungroup() %>% 
  select(-DiagnosisCode) %>% unique()
dim(atlas_enc_phecode) # dim = (2676141,13)
# Map to unique sample ID
atlas_ID_raw = read_csv(file = paste0(raw_data_path, 
                                      "atlas/pheno/raw/atlas_ID_0813.csv"),
                        col_names = FALSE, na = c("NULL"))
names(atlas_ID_raw) = c("PatientID", "UniqueSampleId", "BankSampleID", "CollectDate")
atlas_ID = atlas_ID_raw %>% select(BankSampleID, UniqueSampleId) %>% 
  filter(!is.na(UniqueSampleId)) %>% unique()
# Merge patient data with unique sample ID
enc_demo = atlas_enc_phecode %>% 
  select(PatientID, BankSampleID, Sex, BirthDate, DeathDate, DateFirst, DateLast) %>% 
  inner_join(atlas_ID, by = "BankSampleID") %>% unique() %>%
  select(PatientID, UniqueSampleId, Sex, BirthDate, DeathDate, DateFirst, DateLast) %>%
  filter(!is.na(UniqueSampleId))
dim(enc_demo) # dim = (53123,7)
```

## 3. Finalize patients records and demo info

``` r
atlas_data_path = "/Users/joyfu/Desktop/Projects/Dementia-prediction/data/"
load(file = paste0(atlas_data_path, "atlas/geno/gen_ancestry.rda"))
enc_demo_clean = atlas_enc_phecode %>% group_by(PatientID) %>% 
  mutate(age_first_visit = as.numeric(DateFirst - BirthDate)/365.25) %>%
  mutate(age_last_visit = as.numeric(DateLast - BirthDate)/365.25) %>% 
  mutate(record_length = as.numeric(DateLast - DateFirst)/365.25) %>% 
  mutate(n_diagnosis = length(unique(phecode))) %>% 
  mutate(n_encounter = length(unique(DiagnosisDate))) %>% ungroup() %>% 
  select(PatientID, age_first_visit, age_last_visit, record_length, n_diagnosis, 
         n_encounter) %>% unique() %>% inner_join(enc_demo) %>%
  select(PatientID, UniqueSampleId, Sex, BirthDate, age_first_visit, 
         age_last_visit,record_length, n_diagnosis, n_encounter) %>% 
  unique() %>% inner_join(gen_ancestry) %>% drop_na()
dim(enc_demo_clean) # dim = (26484,10)
save(enc_demo_clean, file = paste0(raw_data_path, "atlas/pheno/mod/enc_demo_clean_",
                                   phe_table_extract_date, ".rda"))
enc_phecode_clean = atlas_enc_phecode %>% 
  filter(PatientID %in% enc_demo_clean$PatientID) %>%
  select(PatientID, DiagnosisDate, phecode) %>% unique() %>% 
  inner_join(enc_demo) %>% 
  select(PatientID, UniqueSampleId, DiagnosisDate, phecode) %>% 
  arrange(PatientID, DiagnosisDate) %>% filter(!is.na(phecode))
dim(enc_phecode_clean) # dim = (1425079,4)
save(enc_phecode_clean, file = paste0(raw_data_path, "atlas/pheno/mod/enc_phecode_clean_",
                                      phe_table_extract_date, ".rda"))
```

# Part 2. Case/control definition (ATLAS patients)

``` r
rm(list = ls())
pacman::p_load(tidyverse, PheWAS)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
phe_table_extract_date = "2024-08-13"
'%!in%' = function(x,y)!('%in%'(x,y))
load(file = paste0(raw_data_path, "atlas/pheno/mod/enc_demo_clean_", 
                   phe_table_extract_date, ".rda"))
load(file = paste0(raw_data_path, "atlas/pheno/mod/enc_phecode_clean_", 
                   phe_table_extract_date, ".rda"))
```

## 1. Define epilepsy

``` r
# Select epilepsy case and controls
epi_phecode = c("345", "345.1", "345.11", "345.12", "345.3")
all_epi_case = enc_phecode_clean %>% 
  filter(phecode %in% epi_phecode) %>% pull(PatientID) %>% unique()
length(all_epi_case) # N = 832
LOE_case_info = enc_phecode_clean %>% 
  filter(phecode %in% epi_phecode) %>% left_join(enc_demo_clean) %>%
  mutate(age_diagnosis = as.numeric(DiagnosisDate - BirthDate)/365.25) %>%
  select(PatientID, UniqueSampleId, DiagnosisDate, phecode, age_diagnosis, 
         age_first_visit) %>% group_by(PatientID, UniqueSampleId) %>% 
  mutate(epi_age_diagnosis = min(age_diagnosis), 
         EpilepsyDate = min(DiagnosisDate)) %>% 
  filter(epi_age_diagnosis >= 60 & age_first_visit < 60) %>% 
  select(PatientID, UniqueSampleId, EpilepsyDate, epi_age_diagnosis) %>% 
  unique() %>% left_join(enc_demo_clean) %>% mutate(epilepsy = 1)
dim(LOE_case_info) # dim = (334,13)
```

## 2. Define AD

``` r
# Select AD case and controls
ad_phecode = "290.11"
all_AD_case = enc_phecode_clean %>% 
  filter(phecode %in% ad_phecode) %>% pull(PatientID) %>% unique()
length(all_AD_case) # N = 453
AD_case_info = enc_phecode_clean %>% 
  filter(phecode %in% ad_phecode) %>% left_join(enc_demo_clean) %>%
  mutate(age_diagnosis = as.numeric(DiagnosisDate - BirthDate)/365.25) %>%
  select(PatientID, UniqueSampleId, DiagnosisDate, phecode, age_diagnosis) %>%
  group_by(PatientID, UniqueSampleId) %>% 
  mutate(AD_age_diagnosis = min(age_diagnosis), ADDate = min(DiagnosisDate)) %>% 
  select(PatientID, UniqueSampleId, ADDate, AD_age_diagnosis) %>% 
  unique() %>% left_join(enc_demo_clean) %>% mutate(AD = 1)
dim(AD_case_info) # dim = (453,12)
```

## 3. Define controls

``` r
exclude_epi_phecode = PheWAS::phecode_exclude %>% 
  filter(code %in% epi_phecode) %>% pull(exclusion_criteria) %>% unique()
exclude_ad_phecode = PheWAS::phecode_exclude %>%
  filter(code %in% ad_phecode) %>% pull(exclusion_criteria) %>% unique()
exclude_control = enc_phecode_clean %>%
  filter(phecode %in% c(epi_phecode, exclude_epi_phecode, ad_phecode, exclude_ad_phecode)) %>% 
  pull(PatientID) %>% unique()
length(exclude_control) # N = 8589
all_control_eligible = enc_phecode_clean %>%
  filter(PatientID %!in% exclude_control) %>% 
  left_join(enc_demo_clean) %>% 
  mutate(age_diagnosis = as.numeric(DiagnosisDate - BirthDate)/365.25, 
         diag_after60 = if_else(age_diagnosis >= 60, 1, 0)) %>%
  group_by(PatientID, UniqueSampleId) %>%
  mutate(n_diag_after60 = sum(diag_after60)) %>% ungroup() %>%
  filter(age_last_visit >= 75 & age_last_visit < 90 & n_diag_after60 >= 2 &
           record_length >= 5 & n_encounter/record_length >= 2) %>%
  pull(PatientID) %>% unique()
length(all_control_eligible) # N = 2031
```

``` r
# for a larger sample size
all_control_eligible = enc_phecode_clean %>%
  filter(PatientID %!in% exclude_control) %>% 
  left_join(enc_demo_clean) %>% 
  mutate(age_diagnosis = as.numeric(DiagnosisDate - BirthDate)/365.25, 
         diag_after60 = if_else(age_diagnosis >= 60, 1, 0)) %>%
  group_by(PatientID, UniqueSampleId) %>%
  mutate(n_diag_after60 = sum(diag_after60)) %>% ungroup() %>%
  filter(age_last_visit >= 60 & age_last_visit < 90 & n_diag_after60 >= 2) %>%
  pull(PatientID) %>% unique()
length(all_control_eligible) # N = 16344
```

## 4. Finalize samples

``` r
# Add PCs and APOE
geno_path = "/Users/joyfu/Desktop/Projects/Dementia-prediction/data/atlas/geno/"
# Add genetic info
load(paste0(geno_path, "apoe.rda"))
load(paste0(geno_path, "pcs/anc_pcs.rda"))
load(file = paste0(raw_data_path, "prs/mod/AD_PRS_final.rda"))
apoe$UniqueSampleId = as.numeric(apoe$UniqueSampleId)
anc_pcs$UniqueSampleId = as.numeric(anc_pcs$UniqueSampleId)
sample_demo_eligible = enc_demo_clean %>% 
  filter(PatientID %in% c(LOE_case_info$PatientID, all_AD_case, all_control_eligible)) %>% 
  left_join(apoe) %>% left_join(anc_pcs) %>% left_join(AD_PRS_final) %>%
  filter(gen_ancestry == "EUR") %>% 
  select(PatientID, UniqueSampleId, Sex, age_last_visit, record_length, 
         n_diagnosis, n_encounter, paste0("PC", c(1:10)), e4count, 
         AD_PRS_full_norm, AD_PRS_rmapoe_norm) %>% unique() %>% drop_na() %>% 
  mutate(AD = if_else(PatientID %in% all_AD_case, 1, 0),
         LOE = if_else(PatientID %in% LOE_case_info$PatientID, 1, 0),
         record_density = n_encounter/record_length) %>% 
  mutate(female = if_else(Sex == "Female", 1, 0)) %>% select(-Sex) %>% 
  left_join(LOE_case_info %>% select(PatientID, EpilepsyDate, epi_age_diagnosis)) %>%
  left_join(AD_case_info %>% select(PatientID, ADDate, AD_age_diagnosis)) %>% 
  select(PatientID, UniqueSampleId, female, age_last_visit, record_length, 
         n_diagnosis, n_encounter, record_density, PC1:PC10, e4count, 
         AD_PRS_full_norm, AD_PRS_rmapoe_norm, 
         AD, ADDate, AD_age_diagnosis, LOE, EpilepsyDate, epi_age_diagnosis) %>%
  unique() %>% filter(record_length > 0)
dim(sample_demo_eligible) # dim = (2179, 27) (12862, 27)
save(sample_demo_eligible, file = paste0(raw_data_path, "atlas/pheno/final/sample_demo_eligible_",
                                      phe_table_extract_date, ".rda"))
save(sample_demo_eligible, file = paste0(raw_data_path, "atlas/pheno/final/sample_demo_eligible_big.rda"))
```

``` r
sample_enc_eligible = enc_phecode_clean %>% 
  filter(PatientID %in% sample_demo_eligible$PatientID) %>% 
  select(PatientID, UniqueSampleId, DiagnosisDate, phecode) %>% unique() %>% drop_na()
dim(sample_enc_eligible) # dim = (166070,4) (592190,4)
save(sample_enc_eligible, file = paste0(raw_data_path, "atlas/pheno/final/sample_enc_eligible_",
                                      phe_table_extract_date, ".rda"))
save(sample_enc_eligible, file = paste0(raw_data_path, "atlas/pheno/final/sample_enc_eligible_big.rda"))
```

### New added - Death date

``` r
phe_table_extract_date = "2024-08-13"
# Read in ATLAS patient data
atlas_full_raw = read_csv(file = paste0(raw_data_path, 
                                        "atlas/pheno/raw/atlas_full_raw_0813.csv"),
                          col_names = FALSE, na = c("NULL"))
names(atlas_full_raw) = c("PatientID", "Sex", "Ethnicity", "FirstRace", 
                          "BirthDate", "DeathDate", "PatientLivingStatus",
                          "DateFirst", "DateLast", "BankSampleID", 
                          "DiagnosisCode", "DiagnosisDate")
dim(atlas_full_raw) # dim = (18341456,12)
```

``` r
atlas_death_info = atlas_full_raw %>% 
  mutate(BirthDate = as.Date(BirthDate), DeathDate = as.Date(DeathDate)) %>%
  mutate(death_age = as.numeric(DeathDate - BirthDate)/365.25) %>%
  select(PatientID, DeathDate, death_age) %>% unique()
dim(atlas_death_info) # dim = (56384,3)
save(atlas_death_info, file = paste0(raw_data_path, "atlas/pheno/mod/atlas_death_info_",
                                      phe_table_extract_date, ".rda"))
```
