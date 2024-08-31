import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LinearRegression, LogisticRegression, ElasticNet, MultiTaskElasticNet
from sklearn.metrics import roc_auc_score, average_precision_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import warnings
# Silence the SettingWithCopyWarning
warnings.filterwarnings("ignore", category = pd.core.common.SettingWithCopyWarning)

# Load in cleaned data
df = pd.read_csv('/opt/genomics/IPHinvestigators/joyfu/dementia-epilepsy/data/modeling/sample_data_full.csv')

# Separate features (SNPs) and target outcomes
X = df.drop(columns = ['AD', 'LOE'])
y = df[['AD', 'LOE']]

# Placeholder for demographics and SNP columns
demographics = ['age_last_visit', 'female', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
SNP_columns = [col for col in df.columns if col.startswith('chr')]

# Generate a list of random seeds
np.random.seed(1314)  # Set seed for reproducibility
random_seeds = np.random.randint(0, 10000, size = 1000)  # Generate 1000 random seeds
# Placeholder lists for recording results
results = []
# Initialize an empty DataFrame to store the results
shared_df = pd.DataFrame()

# Loop over each random seed
for seed in random_seeds:
    X = df.drop(columns = ['AD', 'LOE'])
    y = df[['AD', 'LOE']]
    # Train-test split with the current seed
    y_stratify = y.apply(lambda row: '_'.join(row.values.astype(str)), axis = 1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, stratify = y_stratify, random_state = seed)
    # Standardize demographic features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train[demographics])
    X_test_scaled = scaler.transform(X_test[demographics])
    # Regress out demographics from AD and LOE using logistic regression
    log_reg_ad = LogisticRegression()
    log_reg_ad.fit(X_train_scaled, y_train['AD'])
    y_train_pred_ad = log_reg_ad.predict_proba(X_train_scaled)[:, 1]
    y_train_resid_ad = y_train['AD'] - y_train_pred_ad
    log_reg_loe = LogisticRegression()
    log_reg_loe.fit(X_train_scaled, y_train['LOE'])
    y_train_pred_loe = log_reg_loe.predict_proba(X_train_scaled)[:, 1]
    y_train_resid_loe = y_train['LOE'] - y_train_pred_loe
    # Store residuals
    y_train_resid = pd.DataFrame({
        'AD_residual': y_train_resid_ad,
        'LOE_residual': y_train_resid_loe
    })
    # Apply the same transformation to the test set
    y_test_pred_ad = log_reg_ad.predict_proba(X_test_scaled)[:, 1]
    y_test_resid_ad = y_test['AD'] - y_test_pred_ad
    y_test_pred_loe = log_reg_loe.predict_proba(X_test_scaled)[:, 1]
    y_test_resid_loe = y_test['LOE'] - y_test_pred_loe
    y_test_resid = pd.DataFrame({
        'AD_residual': y_test_resid_ad,
        'LOE_residual': y_test_resid_loe
    })
    
    # Isolate SNP data
    X_train_SNPs = X_train[SNP_columns]
    X_test_SNPs = X_test[SNP_columns]

    # Elastic Net separate
    elnet_ad = ElasticNet(alpha = 0.019, l1_ratio = 0.778) 
    elnet_ad.fit(X_train_SNPs, y_train_resid['AD_residual'])
    y_pred_ad_resid = elnet_ad.predict(X_test_SNPs)
    elnet_loe = ElasticNet(alpha = 0.019, l1_ratio = 0.778) # 0.019, 0.778 # 0.019, 1
    elnet_loe.fit(X_train_SNPs, y_train_resid['LOE_residual'])
    y_pred_loe_resid = elnet_loe.predict(X_test_SNPs)
    auc_ad_elnet = roc_auc_score(y_test['AD'], y_pred_ad_resid)
    auc_loe_elnet = roc_auc_score(y_test['LOE'], y_pred_loe_resid)
    auprc_ad_elnet = average_precision_score(y_test['AD'], y_pred_ad_resid)
    auprc_loe_elnet = average_precision_score(y_test['LOE'], y_pred_loe_resid)
    mse_ad_elnet = mean_squared_error(y_test['AD'], y_pred_ad_resid)
    mse_loe_elnet = mean_squared_error(y_test['LOE'], y_pred_loe_resid)
    coefficients_ad = pd.DataFrame(elnet_ad.coef_.T, index = SNP_columns, columns = ['AD_coef'])
    coefficients_loe = pd.DataFrame(elnet_loe.coef_.T, index = SNP_columns, columns = ['LOE_coef'])
    ad_specific_features = coefficients_ad[(coefficients_ad['AD_coef'] != 0) ]
    loe_specific_features = coefficients_loe[(coefficients_loe['LOE_coef'] != 0)]
    column1 = ad_specific_features.index
    column2 = loe_specific_features.index
    # Find the overlap between the two columns
    overlap = column1[column1.isin(column2)]

    # Multi-Task Elastic Net 
    multi_task_elnet = MultiTaskElasticNet(alpha = 0.019, l1_ratio = 0.667) # 0.019, 0.667  # 0.013, 0.889
    multi_task_elnet.fit(X_train_SNPs, y_train_resid)
    y_pred_resid_multi_task = multi_task_elnet.predict(X_test_SNPs)
    auc_ad_multi_task = roc_auc_score(y_test['AD'], y_pred_resid_multi_task[:, 0])
    auc_loe_multi_task = roc_auc_score(y_test['LOE'], y_pred_resid_multi_task[:, 1])
    auprc_ad_multi_task = average_precision_score(y_test['AD'], y_pred_resid_multi_task[:, 0])
    auprc_loe_multi_task = average_precision_score(y_test['LOE'], y_pred_resid_multi_task[:, 1])
    mse_ad_multi_task = mean_squared_error(y_test['AD'], y_pred_resid_multi_task[:, 0])
    mse_loe_multi_task = mean_squared_error(y_test['LOE'], y_pred_resid_multi_task[:, 1])
    coefficients = pd.DataFrame(multi_task_elnet.coef_.T, index = SNP_columns, columns = ['AD_coef', 'LOE_coef'])
    # Identify shared features (non-zero coefficients in both AD and LOE)
    shared_features = coefficients[(coefficients['AD_coef'] != 0) & (coefficients['LOE_coef'] != 0)]
    shared_features['SNP'] = shared_features.index
    shared_features['seed'] = seed
    shared_df = pd.concat([shared_df, shared_features], ignore_index = True)
    
    # Record the results
    results.append({
        'seed': seed,
        'model': 'Multi-Task Elastic Net',
        'auc_ad': auc_ad_multi_task,
        'auprc_ad': auprc_ad_multi_task,
        'mse_ad': mse_ad_multi_task,
        'auc_loe': auc_loe_multi_task,
        'auprc_loe': auprc_loe_multi_task,
        'mse_loe': mse_loe_multi_task,
        'overlap_single': None,
        
    })

    results.append({
        'seed': seed,
        'model': 'Separate Elastic Net (AD)',
        'auc_ad': auc_ad_elnet,
        'auprc_ad': auprc_ad_elnet,
        'mse_ad': mse_ad_elnet,
        'auc_loe': None,
        'auprc_loe': None,
        'mse_loe': None,
        'overlap_single': overlap
    })

    results.append({
        'seed': seed,
        'model': 'Separate Elastic Net (LOE)',
        'auc_ad': None,
        'auprc_ad': None,
        'mse_ad': None,
        'auc_loe': auc_loe_elnet,
        'auprc_loe': auprc_loe_elnet,
        'mse_loe': mse_loe_elnet,
        'overlap_single': overlap
    })
    
result_df = pd.DataFrame(results, columns=['seed', 'model', 'auc_ad', 'auprc_ad', 'mse_ad', 
                                           'auc_loe', 'auprc_loe', 'mse_loe', 'overlap_single'])

# Output the DataFrame to a CSV file
result_df.to_csv('/opt/genomics/IPHinvestigators/joyfu/dementia-epilepsy/output/result_elnet.csv', index = False)
shared_df.to_csv('/opt/genomics/IPHinvestigators/joyfu/dementia-epilepsy/output/shared_features_elnet.csv', index = False)