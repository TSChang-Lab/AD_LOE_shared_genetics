# This script includes all functions used in the dementia genetic prediction project
# Listed below
# 1. '%!in%': not in
# 2. zeroVar: find variables of zero variance in a dataset
# 3. extract_cv_results: extract cv results from h2o model object

#=================================================
'%!in%' = function(x,y)!('%in%'(x,y))

#=================================================
zeroVar = function(data, useNA = 'ifany') {
  out = apply(data, 2, function(x) {length(table(x, useNA = useNA))})
  which(out == 1)
}

#=================================================
extract_cv_results = function(model_name) {
  result_lst = list()
  cv_metrics = model_name@model$cross_validation_metrics@metrics
  
  MSE = cv_metrics$MSE
  RMSE = cv_metrics$RMSE
  LogLoss = cv_metrics$logloss
  mean_per_class_error = cv_metrics$mean_per_class_error
  AUC = cv_metrics$AUC
  AUCPR = cv_metrics$pr_auc
  Gini = cv_metrics$Gini
  R_2 = cv_metrics$r2
  max_criteria_and_metric_scores = cv_metrics$max_criteria_and_metric_scores
  threshold_f1_index = max_criteria_and_metric_scores %>% filter(metric == "max absolute_mcc") %>% pull(idx)
  max_f1 = max_criteria_and_metric_scores %>% filter(metric == "max absolute_mcc") %>% pull(value)
  # Prediction with threshold F1
  threshold_results = cv_metrics$thresholds_and_metric_scores[(threshold_f1_index+1),]
  
  model_vec = c(MSE, RMSE, LogLoss, mean_per_class_error, AUC, AUCPR, Gini, R_2)
  names(model_vec) = c("MSE", "RMSE", "LogLoss", "mean_per_class_error", "AUC", 
                       "AUCPR", "Gini", "R_2")
  result_lst = c(result_lst, list(model_vec))
  
  perform_vec = list(threshold_results)
  result_lst = c(result_lst, perform_vec)
  return(result_lst)
}

#=================================================
build_coef_table = function(var_names, model) {
  summary_results = lapply(model, summary)
  # Make to a full table
  coef_table = array(NA, dim = c(length(var_names), 6))
  for (i in 1:length(var_names)) {
    coeff = summary_results[[i]]$coefficients[2]
    confinterval = confint(model[[i]], level = 0.95)
    lower_CI = confinterval[2]
    upper_CI = confinterval[8]
    CI_95 = paste0(sprintf('%.2f', coeff), ' (', sprintf('%.2f',lower_CI), ', ', 
                   sprintf('%.2f',upper_CI), ')')
    coef_table[i,1] = var_names[i]
    coef_table[i,2] = length(summary_results[[i]]$residuals)
    coef_table[i,3] = coeff
    coef_table[i,4] = lower_CI
    coef_table[i,5] = upper_CI
    coef_table[i,6] = CI_95
    i = i + 1
  }
  coef_table_final = as.data.frame(coef_table)
  colnames(coef_table_final) = c('var_name', 'N', 'beta', 
                                 'lower_CI', 'upper_CI', 'text')
  return(coef_table_final)
}

#=================================================
loe_set = c("ABCA7", "AC011558.5", "CCDC89", "CNN2", "CREBZF","EED", "HMHA1", 
            "PICALM", "PTK2B", "RNU6-2", "RNU6-560P", "SETP17", "snoU13")
ad_set = c("ABCA7", "AC005779.1", "AC005779.2", "AC006126.3", "AC006126.4",
           "AC092066.6", "AL691452.1", "AP001257.1", "APOC1", "APOC1P1", "APOC2",
           "APOC4", "APOC4-APOC2", "APOE", "BCAM", "BCL3", "BIN1", "BLOC1S3",
           "CBLC", "CCDC89", "CD46P1", "CEACAM16", "CEACAM19", "CLASRP",
           "CLPTM1", "CR1", "CR1L", "CREBZF", "CTB-129P6.11", "CTB-129P6.4",
           "CTB-171A8.1", "DYRK3", "EED", "EXOC3L2", "GEMIN7", "HMHA1", "IGSF23",
           "MARK4", "MS4A2", "MS4A4A", "MS4A4E", "MS4A6A", "MS4A6E", "PFKFB2",
           "PICALM", "POLR2E", "PVR", "PVRL2", "RELB", "RP11-343H5.6",
           "RP11-730K11.1", "RP11-736I10.2", "RP11-78B10.2", "SETP17", 
           "SLC25A1P1", "SNORA70", "snoZ6", "SORL1", "TOMM40", "TRAPPC6A",
           "YOD1", "ZNF112")
shared_set = c("ABCA7", "AC005779.1", "AC005779.2", "AC006126.3", "AC006126.4",
               "AC011558.5", "AC138472.4", "AL355353.1", "AP001257.1", "APOC1",
               "APOC1P1", "APOC2", "APOC4", "APOC4-APOC2", "APOE", "BCAM", "BCL3",
               "BIN1", "BLOC1S3", "CBLC", "CCDC25", "CD2AP", "CEACAM22P",
               "CENPQ", "CLASRP", "CLPTM1", "CLU", "CNN2", "CTB-129P6.11",
               "CTB-129P6.4", "ELP3", "ESCO2", "EXOC3L2", "GPR115", "HMHA1",
               "MARK4", "MS4A2", "MS4A3", "MS4A4A", "MS4A4E", "MS4A6A", "MUT",
               "NUGGC", "PNOC", "POLR2E", "PPP1R37", "PVRL2", "RELB", "RNU6-1276P",
               "RNU6-2", "RNU6-611P", "RP1-251M9.3", "RP11-138I18.1",
               "RP11-138I18.2", "RP11-16P20.3", "RP11-385F7.1", "RP11-736I10.2",
               "RP11-812I20.2", "SLC24A4", "TNFRSF21", "TOMM40", "TRAPPC6A",
               "USP6NL", "Y_RNA", "ZNF112")