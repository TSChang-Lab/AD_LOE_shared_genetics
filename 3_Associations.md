Association tests between AD and LOE
================
Mingzhou Fu
2024-08-29

``` r
rm(list = ls())
pacman::p_load(tidyverse, tidymodels, gtsummary, cmprsk)
raw_data_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/data/"
output_path = "/Users/joyfu/Desktop/Projects/AD-LOE-genetics/output/"
'%!in%' = function(x,y)!('%in%'(x,y))
phe_table_extract_date = "2024-08-13"
load(file = paste0(raw_data_path, "atlas/pheno/final/sample_enc_eligible_big.rda"))
load(file = paste0(raw_data_path, "atlas/pheno/final/sample_demo_eligible_big.rda"))
load(file = paste0(raw_data_path, "atlas/pheno/mod/atlas_death_info_", 
                   phe_table_extract_date, ".rda"))
```

# Part 1. Sample preparation

## 1. Add covariates

``` r
# HTN
htn_case_info = sample_enc_eligible %>% 
  filter(phecode %in% c("401.1")) %>% group_by(UniqueSampleId) %>% 
  mutate(HTNDate = min(DiagnosisDate), HTN = 1) %>% ungroup() %>%
  select(PatientID, UniqueSampleId, HTN, HTNDate) %>% unique()
dim(htn_case_info) # dim = (7238,4)
# Diabetes
diab_case_info = sample_enc_eligible %>% 
  filter(phecode %in% c("250.2")) %>% group_by(UniqueSampleId) %>% 
  mutate(DiabetesDate = min(DiagnosisDate), diabetes = 1) %>% ungroup() %>%
  select(PatientID, UniqueSampleId, diabetes, DiabetesDate) %>% unique()
dim(diab_case_info) # dim = (2147,4)
# Stroke
stroke_case_info = sample_enc_eligible %>% 
  filter(phecode %in% c("433")) %>% group_by(UniqueSampleId) %>% 
  mutate(StrokeDate = min(DiagnosisDate), stroke = 1) %>% ungroup() %>%
  select(PatientID, UniqueSampleId, stroke, StrokeDate) %>% unique()
dim(stroke_case_info) # dim = (111,4)
# Hyperlipid
hyperlipid_case_info = sample_enc_eligible %>% 
  filter(phecode %in% c("272.1")) %>% group_by(UniqueSampleId) %>% 
  mutate(HyperlipidDate = min(DiagnosisDate), hyperlipid = 1) %>% ungroup() %>%
  select(PatientID, UniqueSampleId, hyperlipid, HyperlipidDate) %>% unique()
dim(hyperlipid_case_info) # dim = (7170,4)
```

## 2. Finalize association sample

``` r
patient_final = sample_demo_eligible %>% 
  left_join(htn_case_info) %>% left_join(diab_case_info) %>% 
  left_join(stroke_case_info) %>% left_join(hyperlipid_case_info) %>% 
  mutate(HTN = if_else(!is.na(HTN), 1, 0),
         diabetes = if_else(!is.na(diabetes), 1, 0),
         stroke = if_else(!is.na(stroke), 1, 0),
         hyperlipid = if_else(!is.na(hyperlipid), 1, 0)) %>% unique() %>% 
  left_join(atlas_death_info) %>% unique() %>% 
  mutate(deceased = if_else(!is.na(DeathDate), 1, 0)) 
dim(patient_final) # dim = (12862, 38)
save(patient_final, file = paste0(raw_data_path, "modeling/patient_final.rda"))
```

# Part 2. Descriptive statistics

``` r
# Desc by epilepsy status
patient_final %>% mutate(LOE = as.factor(LOE)) %>% 
  select(LOE, AD, age_last_visit, female, record_length, record_density, 
         deceased, HTN, diabetes, stroke, hyperlipid, e4count) %>% 
  tbl_summary(by = LOE) %>% add_p()
```

``` r
# Desc by AD status
patient_final %>% mutate(AD = as.factor(AD)) %>% 
  select(AD, LOE, age_last_visit, female, record_length, record_density, 
         deceased, HTN, diabetes, stroke, hyperlipid, e4count) %>% 
  tbl_summary(by = AD) %>% add_p()
```

# Part 3. Survival modeling

## 1. AD ~ LOE

### Cox proportional-Hazards

``` r
patient_AD_longitudinal = patient_final %>% 
  mutate(LOE_status = case_when(
    LOE == 1 & (EpilepsyDate <= ADDate | is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(AD_status = case_when(
    AD == 1 & (ADDate >= EpilepsyDate | is.na(EpilepsyDate)) ~ 1,
    AD == 1 & ADDate < EpilepsyDate ~ 99,
    TRUE ~ 0
  )) %>% filter(AD_status != 99) %>%
  mutate(baseline_age = case_when(
    LOE_status == 1 ~ epi_age_diagnosis,
    TRUE ~ age_last_visit - record_length
  )) %>%
  mutate(follow_up_time = case_when(
    LOE_status == 1 & AD_status == 1 ~ AD_age_diagnosis - epi_age_diagnosis,
    LOE_status == 1 & AD_status != 1 ~ age_last_visit - epi_age_diagnosis,
    LOE_status == 0 & AD_status == 1 ~ AD_age_diagnosis - baseline_age,
    LOE_status == 0 & AD_status != 1 ~ age_last_visit - baseline_age
  )) %>%
  mutate(baseline_HTN = case_when(
    HTN == 1 & (HTNDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(baseline_diabetes = case_when(
    diabetes == 1 & (DiabetesDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_stroke = case_when(
    stroke == 1 & (StrokeDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_hyperlipid = case_when(
    hyperlipid == 1 & (HyperlipidDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  select(PatientID, UniqueSampleId, AD_status, LOE_status, baseline_age, female,
         follow_up_time, record_length, record_density, baseline_HTN, 
         baseline_diabetes, baseline_stroke, baseline_hyperlipid, e4count, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
         unique() %>% drop_na() %>% filter(follow_up_time > 0)
dim(patient_AD_longitudinal) # dim = (12854,24)
```

``` r
cox_demo = coxph(Surv(follow_up_time, AD_status) ~ LOE_status + baseline_age + 
                    female + record_length + record_density + PC1 + PC2 + PC3 + 
                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                    data = patient_AD_longitudinal)
cox_health = coxph(Surv(follow_up_time, AD_status) ~ LOE_status + baseline_age + 
                     female + record_length + record_density + PC1 + PC2 + PC3 + 
                     PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + baseline_HTN + 
                     baseline_diabetes + baseline_stroke, 
                   data = patient_AD_longitudinal)

sum_res = rbind(tidy(cox_demo, conf.int = TRUE, exponentiate = TRUE)[1,],
                tidy(cox_health, conf.int = TRUE, exponentiate = TRUE)[1,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(RR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(RR, p_value)
sum_res
```

### FG modeling

``` r
patient_AD_longitudinal = patient_final %>% 
  mutate(LOE_status = case_when(
    LOE == 1 & (EpilepsyDate <= ADDate | is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(AD_status = case_when(
    AD == 1 & (ADDate >= EpilepsyDate | is.na(EpilepsyDate)) ~ 1,
    AD == 1 & ADDate < EpilepsyDate ~ 99,
    AD == 0 & (!is.na(DeathDate)) ~ 2, 
    TRUE ~ 0
  )) %>% filter(AD_status != 99) %>%
  mutate(baseline_age = case_when(
    LOE_status == 1 ~ epi_age_diagnosis,
    TRUE ~ age_last_visit - record_length
  )) %>%
  mutate(follow_up_time = case_when(
    LOE_status == 1 & AD_status == 1 ~ AD_age_diagnosis - epi_age_diagnosis,
    LOE_status == 1 & AD_status != 1 ~ age_last_visit - epi_age_diagnosis,
    LOE_status == 0 & AD_status == 1 ~ AD_age_diagnosis - baseline_age,
    LOE_status == 0 & AD_status != 1 ~ age_last_visit - baseline_age
  )) %>%
  mutate(baseline_HTN = case_when(
    HTN == 1 & (HTNDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(baseline_diabetes = case_when(
    diabetes == 1 & (DiabetesDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_stroke = case_when(
    stroke == 1 & (StrokeDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_hyperlipid = case_when(
    hyperlipid == 1 & (HyperlipidDate <= EpilepsyDate | !is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  select(PatientID, UniqueSampleId, AD_status, LOE_status, baseline_age, female,
         follow_up_time, record_length, record_density, baseline_HTN, 
         baseline_diabetes, baseline_stroke, baseline_hyperlipid, e4count, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
         unique() %>% drop_na() %>% filter(follow_up_time > 0)
dim(patient_AD_longitudinal) # dim = (12854,24)
```

``` r
cov_demo = model.matrix(~ LOE_status + baseline_age + female + record_length + 
                          record_density + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                          PC7 + PC8 + PC9 + PC10, 
                        data = patient_AD_longitudinal)[, -1]
crr_demo = crr(cov1 = cov_demo, ftime = patient_AD_longitudinal$follow_up_time, 
               fstatus = patient_AD_longitudinal$AD_status, failcode = 1,
               cencode = 0)

cov_health = model.matrix(~ LOE_status + baseline_age + female + record_length + 
                            record_density + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                            PC7 + PC8 + PC9 + PC10 + baseline_HTN + 
                            baseline_diabetes + baseline_stroke, 
                          data = patient_AD_longitudinal)[, -1]
crr_health = crr(cov1 = cov_health, ftime = patient_AD_longitudinal$follow_up_time, 
                fstatus = patient_AD_longitudinal$AD_status, failcode = 1,
                cencode = 0)

sum_res = rbind(tidy(crr_demo, conf.int = TRUE, exponentiate = TRUE)[1,],
                tidy(crr_health, conf.int = TRUE, exponentiate = TRUE)[1,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(RR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(RR, p_value)
sum_res
```

## 2. LOE ~ AD

### Cox proportional-Hazards

``` r
patient_LOE_longitudinal = patient_final %>% 
  mutate(AD_status = case_when(
    AD == 1 & (ADDate <= EpilepsyDate | is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(LOE_status = case_when(
    LOE == 1 & (EpilepsyDate >= ADDate | is.na(ADDate)) ~ 1,
    LOE == 1 & EpilepsyDate < ADDate ~ 99,
    TRUE ~ 0
  )) %>% filter(LOE_status != 99) %>%
  mutate(baseline_age = case_when(
    AD_status == 1 ~ AD_age_diagnosis,
    TRUE ~ age_last_visit - record_length
  )) %>%
  mutate(follow_up_time = case_when(
    AD_status == 1 & LOE_status == 1 ~ epi_age_diagnosis - AD_age_diagnosis,
    AD_status == 1 & LOE_status != 1 ~ age_last_visit - AD_age_diagnosis,
    AD_status == 0 & LOE_status == 1 ~ epi_age_diagnosis - baseline_age,
    AD_status == 0 & LOE_status != 1 ~ age_last_visit - baseline_age
  )) %>%
  mutate(baseline_HTN = case_when(
    HTN == 1 & (HTNDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(baseline_diabetes = case_when(
    diabetes == 1 & (DiabetesDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_stroke = case_when(
    stroke == 1 & (StrokeDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_hyperlipid = case_when(
    hyperlipid == 1 & (HyperlipidDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  select(PatientID, UniqueSampleId, AD_status, LOE_status, baseline_age, female,
         follow_up_time, record_length, record_density, baseline_HTN, 
         baseline_diabetes, baseline_stroke, baseline_hyperlipid, e4count, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, 
         PC9, PC10) %>% unique() %>% drop_na() %>% filter(follow_up_time > 0)
dim(patient_LOE_longitudinal) # dim = (12852,25)
```

``` r
cox_demo = coxph(Surv(follow_up_time, LOE_status) ~ AD_status + baseline_age + 
                    female + record_length + record_density + PC1 + PC2 + PC3 + 
                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                    data = patient_LOE_longitudinal)
cox_health = coxph(Surv(follow_up_time, LOE_status) ~ AD_status + baseline_age + 
                     female + record_length + record_density + PC1 + PC2 + PC3 + 
                     PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + baseline_HTN + 
                     baseline_diabetes + baseline_stroke, 
                   data = patient_LOE_longitudinal)

sum_res = rbind(tidy(cox_demo, conf.int = TRUE, exponentiate = TRUE)[1,],
                tidy(cox_health, conf.int = TRUE, exponentiate = TRUE)[1,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(RR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(RR, p_value)
sum_res
```

### FG modeling

``` r
patient_LOE_longitudinal = patient_final %>% 
  mutate(AD_status = case_when(
    AD == 1 & (ADDate <= EpilepsyDate | is.na(EpilepsyDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(LOE_status = case_when(
    LOE == 1 & (EpilepsyDate >= ADDate | is.na(ADDate)) ~ 1,
    LOE == 1 & EpilepsyDate < ADDate ~ 99,
    LOE == 0 & (!is.na(DeathDate)) ~ 2, 
    TRUE ~ 0
  )) %>% filter(LOE_status != 99) %>%
  mutate(baseline_age = case_when(
    AD_status == 1 ~ AD_age_diagnosis,
    TRUE ~ age_last_visit - record_length
  )) %>%
  mutate(follow_up_time = case_when(
    AD_status == 1 & LOE_status == 1 ~ epi_age_diagnosis - AD_age_diagnosis,
    AD_status == 1 & LOE_status != 1 ~ age_last_visit - AD_age_diagnosis,
    AD_status == 0 & LOE_status == 1 ~ epi_age_diagnosis - baseline_age,
    AD_status == 0 & LOE_status != 1 ~ age_last_visit - baseline_age
  )) %>%
  mutate(baseline_HTN = case_when(
    HTN == 1 & (HTNDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(baseline_diabetes = case_when(
    diabetes == 1 & (DiabetesDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_stroke = case_when(
    stroke == 1 & (StrokeDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(baseline_hyperlipid = case_when(
    hyperlipid == 1 & (HyperlipidDate <= ADDate | !is.na(ADDate)) ~ 1,
    TRUE ~ 0
  )) %>%
  select(PatientID, UniqueSampleId, AD_status, LOE_status, baseline_age, female,
         follow_up_time, record_length, record_density, baseline_HTN, 
         baseline_diabetes, baseline_stroke, baseline_hyperlipid, e4count, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
         unique() %>% drop_na() %>% filter(follow_up_time > 0)
dim(patient_LOE_longitudinal) # dim = (12852,25)
```

``` r
cov_demo = model.matrix(~ AD_status + baseline_age + female + record_length + 
                           record_density + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                           PC7 + PC8 + PC9 + PC10, 
                         data = patient_LOE_longitudinal)[, -1]
crr_demo = crr(cov1 = cov_demo, ftime = patient_LOE_longitudinal$follow_up_time, 
                fstatus = patient_LOE_longitudinal$LOE_status, failcode = 1,
                cencode = 0)

cov_health = model.matrix(~ AD_status + baseline_age + female + record_length + 
                            record_density + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                            PC7 + PC8 + PC9 + PC10 + baseline_HTN + 
                            baseline_diabetes + baseline_stroke, 
                          data = patient_LOE_longitudinal)[, -1]
crr_health = crr(cov1 = cov_health, ftime = patient_LOE_longitudinal$follow_up_time, 
                fstatus = patient_LOE_longitudinal$LOE_status, failcode = 1,
                cencode = 0)

sum_res = rbind(tidy(crr_demo, conf.int = TRUE, exponentiate = TRUE)[1,],
                tidy(crr_health, conf.int = TRUE, exponentiate = TRUE)[1,]) %>% 
  select(estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 2))) %>% 
  mutate(RR = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  mutate(p_value = round(p.value, 3)) %>% 
  select(RR, p_value)
sum_res
```
