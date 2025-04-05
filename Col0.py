import numpy as np
import pandas as pd
import xgboost as xgb
import torch
import matplotlib.pyplot as plt
import random
from itertools import combinations

##########################################################
##                  Raw XGboosst
##########################################################
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, classification_report, confusion_matrix, roc_auc_score,
    precision_recall_curve, average_precision_score, f1_score
)

# Reproducibility
seed = 42
np.random.seed(seed)
torch.manual_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)

# Load data
data = pd.read_csv("Drought_category_prediction.csv")
y = data['Drought_tolerance'].astype(int)  # Class 1 = drought-tolerant
feature_columns = [col for col in data.columns if "_count" in col]
X = data[feature_columns]

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.2, random_state=42)

# Scale
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# XGBoost Model
xgb_model = xgb.XGBClassifier(objective='binary:logistic', use_label_encoder=False, eval_metric='logloss')

# Grid search parameters
param_grid = {
    'n_estimators': [100, 250],
    'max_depth': [2, 3],
    'learning_rate': [0.3],
    'gamma': [0],
    'subsample': [1],
    'colsample_bytree': [0.8],
    'min_child_weight': [1]
}

grid_search = GridSearchCV(
    estimator=xgb_model,
    param_grid=param_grid,
    scoring='roc_auc',
    cv=3,
    verbose=1,
    n_jobs=-1
)

grid_search.fit(X_train_scaled, y_train)
best_model = grid_search.best_estimator_

# Predict
y_pred_probs = best_model.predict_proba(X_test_scaled)[:, 1]
threshold = 0.45
y_pred = (y_pred_probs >= threshold).astype(int)

# Evaluation
conf_matrix = confusion_matrix(y_test, y_pred)
accuracy = accuracy_score(y_test, y_pred)
roc_auc = roc_auc_score(y_test, y_pred_probs)
print("Confusion Matrix:\n", conf_matrix)
print("Accuracy:", round(accuracy, 4))
print("ROC AUC:", round(roc_auc, 4))
print(classification_report(y_test, y_pred))

# Save predictions
prediction_output = pd.DataFrame({
    'sequence_name': data.loc[y_test.index, 'sequence_name'] if 'sequence_name' in data.columns else y_test.index,
    'true_label': y_test.values,
    'predicted_label': y_pred,
    'prob_Class0': 1 - y_pred_probs,
    'prob_Class1': y_pred_probs
})
prediction_output.to_csv("xgboost_predictions_drought_tolerance.csv", index=False)

# Feature Importance
importance_df = pd.DataFrame(
    best_model.feature_importances_,
    index=feature_columns,
    columns=['Importance']
).sort_values(by='Importance', ascending=False)

importance_df.head(20).plot(kind='bar', figsize=(10, 5), title='Top 20 Feature Importances')
plt.tight_layout()
plt.show()

# PR Curve
precision, recall, _ = precision_recall_curve(y_test, y_pred_probs)
ap = average_precision_score(y_test, y_pred_probs)
plt.figure(figsize=(8, 5))
plt.plot(recall, precision, marker='.')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title(f'Precision-Recall Curve (AP = {ap:.2f})')
plt.grid(True)
plt.tight_layout()
plt.show()


############precision recall
from sklearn.metrics import f1_score

thresholds = np.arange(0.1, 0.9, 0.01)
f1s = []
for t in thresholds:
    y_pred = (y_pred_probs >= t).astype(int)
    f1s.append(f1_score(y_test, y_pred))

plt.plot(thresholds, f1s)
plt.xlabel("Threshold")
plt.ylabel("Drought tolerance Class  F1-score")
plt.title("Drought tolerance Class  F1 vs. Threshold")
plt.grid(True)
plt.show()



#######################################################
##                      Hypertunning + scale        ##
#########################################################
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier

# Rough estimate from your class balance
scale_weight = (y_train == 0).sum() / (y_train == 1).sum()

# Define model with initial imbalance handling
xgb_model = XGBClassifier(
    objective='binary:logistic',
    use_label_encoder=False,
    eval_metric='logloss',
    scale_pos_weight=scale_weight
)

# Grid parameters (feel free to expand this later)
param_grid = {
    'n_estimators': [100, 200],
    'max_depth': [3, 5],
    'learning_rate': [0.01, 0.1],
    'gamma': [0, 1],
    'subsample': [0.7, 1.0],
    'colsample_bytree': [0.7, 1.0]
}

# Use F1-score for class 1
grid_search = GridSearchCV(
    estimator=xgb_model,
    param_grid=param_grid,
    scoring='f1',
    cv=3,
    verbose=2,
    n_jobs=-1
)

grid_search.fit(X_train_scaled, y_train)

# Final model
best_model = grid_search.best_estimator_
print("Best parameters:", grid_search.best_params_)






#############################################################
##                          Final                          ##
#############################################################
# Re-create the scaled training set with column names
# Use scaled data but keep column names
X_train_df = pd.DataFrame(X_train_scaled, columns=feature_columns)
X_test_df = pd.DataFrame(X_test_scaled, columns=feature_columns)




# Retrain with feature names intact
final_model = xgb.XGBClassifier(
    objective='binary:logistic',
    use_label_encoder=False,
    eval_metric='logloss',
    colsample_bytree=0.7,
    gamma=1,
    learning_rate=0.1,
    max_depth=5,
    n_estimators=200,
    subsample=0.7,
    scale_pos_weight=2.0,
    random_state=42
)

final_model.fit(X_train_df, y_train)

# And also use the DataFrame version for prediction
y_pred_prob = final_model.predict_proba(X_test_df)[:, 1]


# Apply best threshold from earlier
best_thresh = 0.40
y_pred = (y_pred_prob >= best_thresh).astype(int)

# Final evaluation
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, f1_score

print("📊 Final Evaluation")
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("Classification Report:\n", classification_report(y_test, y_pred))
print("ROC AUC:", roc_auc_score(y_test, y_pred_prob))
print("Class 1 F1-score:", f1_score(y_test, y_pred, pos_label=1))

#### Class 1 accuaracy
# From confusion matrix: TP and FN for class 1
# From confusion matrix: TP and FN for class 1
cm = confusion_matrix(y_test, y_pred)
TP = cm[1, 1]
FN = cm[1, 0]

# Class 1 accuracy (same as recall for class 1)
class_1_accuracy = TP / (TP + FN)
print(f"Class 1 Accuracy: {class_1_accuracy:.4f}")




################################################################
##                          Top 20 features                    ##     
#################################################################

# Get feature importances from the trained model
importance_df = pd.DataFrame({
    'Feature': final_model.get_booster().feature_names,
    'Importance': final_model.feature_importances_
}).sort_values(by='Importance', ascending=False)

top_20_feats = importance_df['Feature'].head(20).tolist()
print("Top 20 features:\n", top_20_feats)



##################################################################
##                          Shapley                             ##
##################################################################
import shap

# Build SHAP explainer using the DataFrame
explainer = shap.Explainer(final_model, X_train_df)

# Compute SHAP values on the test set
shap_values = explainer(X_test_df)

# Slice SHAP values and features for the top 20 only
top_20_index = [X_train_df.columns.get_loc(feat) for feat in top_20_feats]
shap_values_top20 = shap_values[:, top_20_index]
X_test_top20 = X_test_df[top_20_feats]
final_model.fit(X_train_df, y_train)
shap.summary_plot(shap_values_top20.values, features=X_test_top20, feature_names=top_20_feats)



################################################################
##
###################################################################
from itertools import combinations
from sklearn.metrics import f1_score, accuracy_score
import pandas as pd
import random

# Start from your top 20 features
top_feats = top_20_feats  # from SHAP or feature importance
results = []

# Define range of interaction sizes (2 to 6)
for k in range(2, 6):
    sampled_combos = random.sample(list(combinations(top_feats, k)), 100)  # Sample 100 combos per k

    for combo in sampled_combos:
        # Subset data
        X_train_sub = X_train_df[list(combo)]
        X_test_sub = X_test_df[list(combo)]

        # Retrain on subset
        model = xgb.XGBClassifier(
            objective='binary:logistic',
            use_label_encoder=False,
            eval_metric='logloss',
            scale_pos_weight=2.0,  # class imbalance
            max_depth=5,
            learning_rate=0.1,
            n_estimators=200,
            subsample=0.7,
            colsample_bytree=0.7,
            gamma=1,
            random_state=42
        )
        model.fit(X_train_sub, y_train)
        y_pred_prob = model.predict_proba(X_test_sub)[:, 1]
        y_pred = (y_pred_prob >= 0.35).astype(int)

        # Score
        f1 = f1_score(y_test, y_pred, pos_label=1)
        acc = accuracy_score(y_test, y_pred)
        roc = roc_auc_score(y_test, y_pred_prob)

        results.append({
            'features': combo,
            'f1_class1': f1,
            'accuracy': acc,
            'roc_auc': roc,
            'k': k
        })
# Create DataFrame
combo_df = pd.DataFrame(results)

# Top 20 combinations by F1
top_20_combos = combo_df.sort_values(by='f1_class1', ascending=False).head(20)
top_20_combos.to_csv("top20_motif_combos.csv", index=False)

# Plot
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 7))
bar_labels = [
    ' + '.join(feats) + f"\nF1={f1:.2f}"
    for feats, f1 in zip(top_20_combos['features'], top_20_combos['f1_class1'])
]

plt.barh(range(len(top_20_combos)), top_20_combos['f1_class1'], color='teal')
plt.yticks(ticks=range(len(top_20_combos)), labels=bar_labels, fontsize=9)
plt.xlabel("Class 1 F1-score")
plt.title("🔝 Top 20 Motif Combinations for Predicting Drought Tolerance")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.grid(axis='x', linestyle='--')
plt.show()
