import os
import math
import pandas as pd
import numpy as np
from joblib import dump
import matplotlib.pyplot as plt
import sklearn.preprocessing as preprocessing

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix, accuracy_score, roc_auc_score

from common.pharmacy_common import PharmacyCommon
from tabulate import tabulate

#Configs
project_path        = "/home/mrcong/Code/phamarcy_code/qsar_ml_find_hdac2_inhibitor"
data_source         = os.path.join(f"{project_path}/data_for_modeling/train_test_data", "new_HDAC2_train_test_data.xlsx")
#MACCS, Morgan2
fpts_name           = "Morgan2"
bits                = 2048
model_result_file   = os.path.join(f"{project_path}/results/KQDD/20240315/files", f"{fpts_name}_{bits}bits_model_result.xlsx")
roc_fig_path        = os.path.join(f"{project_path}/results/KQDD/20240315/files", f"{fpts_name}_{bits}bits_roc_fig.png")
model_saving_path   = f"{project_path}/results/KQDD/20240315/models/{fpts_name}_{bits}_bits"

#Create dir
if not os.path.exists(model_saving_path):
#printing if the path exists or not
    os.makedirs(model_saving_path)

print(f"[+] Start training model with fingerprint: {fpts_name} and bits {bits}")
#Init
common = PharmacyCommon()

#Import data
train_dataset = pd.read_excel(data_source, sheet_name='train_dataset')
validation_dataset = pd.read_excel(data_source, sheet_name='validation_dataset')
test_dataset = pd.read_excel(data_source, sheet_name='test_dataset')
print(len(train_dataset), len(test_dataset), len(validation_dataset))
"""
Ultilities
"""
def gen_fpts(data) -> np.ndarray:
    if fpts_name == "Morgan2":
        return common.gen_ecfp4_fpts(data=data, bits=bits)
    elif fpts_name == "MACCS":
        return common.gen_maccs_fpts(data=data)
    elif fpts_name == "ECFP6":
        return common.gen_ecfp6_fpts(data=data, bits=bits)
    elif fpts_name == "Rdkit":
        return common.gen_rdkit_fpts(data=data, bits=bits)

def model_evaluation_calculation(cm):
    tp = cm[0][0]; tn = cm[1][1]; fp = cm[0][1]; fn = cm[1][0]
    ac = (tp+tn)/(tp+tn+fp+fn)
    mcc = (tp*tn - fp*fn) / math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    precision = tp / (tp +fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)
    return ac, precision, recall, mcc, f1

def me_result(cm, model_name):
    cm_string = "Confusion matrix of " + model_name
    print(cm_string)
    print(cm)
    ac, se, sp, mcc, f1 = model_evaluation_calculation(cm)
    print("Comparision:")
    table = [['Model', 'AC', 'SE', 'SP', 'MCC', 'F1'], [model_name, ac, se, sp, mcc, f1]]
    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
""""""

#X data
X_Train = gen_fpts(train_dataset['SMILES'])
X_Test = gen_fpts(test_dataset['SMILES'])
X_Validation = gen_fpts(validation_dataset['SMILES'])

#y data
y_Train = np.array(train_dataset['Categories'])
y_Test = np.array(test_dataset['Categories'])
y_Validation = np.array(validation_dataset['Categories'])

#Original data
print("Original data:")
print(y_Train[0:5])
print(y_Test[0:5])
print(y_Validation[0:5])

#Label encoder
label_encoder = preprocessing.LabelEncoder()
y_Train = label_encoder.fit_transform(y_Train)
y_Test = label_encoder.transform(y_Test)
y_Validation = label_encoder.transform(y_Validation)
#Class encoded
print("Class encoded:")
print(list(label_encoder.classes_))
print(label_encoder.transform(label_encoder.classes_))
print("Encoded data:")
print(y_Train[0:5])
print(y_Test[0:5])
print(y_Validation[0:5])

#Model training
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
print("[+] Model training:")
rf = RandomForestClassifier(criterion='entropy', n_estimators=200, random_state=42)
rf.fit(X_Train, y_Train)

knn = KNeighborsClassifier(n_neighbors=5, metric='minkowski', p=2)
knn.fit(X_Train, y_Train)

svm = SVC(kernel='rbf', probability=True, random_state=42)
svm.fit(X_Train, y_Train)

xgb = XGBClassifier(objective='binary:logistic',n_estimators=200, tree_method="hist")
xgb.fit(X_Train, y_Train)

#Model evaluations
print("[+] Model evaluation")
X_Total = np.concatenate((X_Train, X_Validation), axis=0)
y_Total = np.concatenate((y_Train, y_Validation), axis=0)
cols = ['Model', 'Fingerprint', 
        '10-fold cross validation AC (only train and validation sets)', 'Test-set Accuracy', 'Test-set Precision', 'Test-set Recall', 'Test-set MCC', 'Test-set F1', 'Test-set AUC', 
        'Training error', 'Test error', 'Validation error',
        'Validation-set Accuracy', 'Validation-set Precision', 'Validation-set Recall', 'Validation-set MCC', 'Validation-set F1', 'Validation-set AUC']

model_result = None
cv = KFold(n_splits=10, random_state=1, shuffle=True)

for model_name in ['KNN', 'RF', 'SVM', 'XgBoost']:
    cv_scores = None
    y_pred_test = None
    y_pred_train = None
    y_pred_validation = None
    model = None
    
    if model_name == 'KNN':
        model = knn
    elif model_name == 'RF':
        model = rf
    elif model_name == 'SVM':
        model = svm
    elif model_name == 'XgBoost':
        model = xgb
    else:
        print("Error")
    
    #Cross validation
    cv_scores = cross_val_score(model, X_Total, y_Total, scoring='accuracy', cv=cv, n_jobs=-1)
    
    #Test set
    y_pred_test = model.predict(X_Test)
    y_proba_test = model.predict_proba(X_Test)[:, 1]
    auc_score_test = roc_auc_score(y_Test, y_proba_test)
    cm_test = confusion_matrix(y_Test, y_pred_test)
    test_ac, test_precision, test_recall, test_mcc, test_f1 = model_evaluation_calculation(cm_test)
    test_error = 1 - test_ac
    
    #validation set
    y_pred_validation = model.predict(X_Validation)
    y_proba_validation = model.predict_proba(X_Validation)[:, 1]
    auc_score_validation = roc_auc_score(y_Validation, y_proba_validation)
    cm_validation = confusion_matrix(y_Validation, y_pred_validation)
    validation_ac, validation_precision, validation_recall, validation_mcc, validation_f1 = model_evaluation_calculation(cm_validation)
    validation_error = 1 - validation_ac
    
    #train set
    y_pred_train = model.predict(X_Train)
    y_proba_train = model.predict_proba(X_Train)[:, 1]
    auc_score_train = roc_auc_score(y_Train, y_proba_train)
    cm_train = confusion_matrix(y_Train, y_pred_train)
    train_ac, _, _, _, _ = model_evaluation_calculation(cm_train)
    train_error = 1 - train_ac
    
    row_result = pd.DataFrame([[model_name, fpts_name, 
                                np.mean(cv_scores)*100, test_ac*100, test_precision*100, test_recall*100, test_mcc*100, test_f1*100, auc_score_test*100, 
                                train_error*100, test_error*100, validation_error*100,
                                validation_ac*100, validation_precision*100, validation_recall*100, validation_mcc*100, validation_f1*100, auc_score_validation*100]], 
                                columns=cols)
    if model_result is None:
        model_result = row_result
    else:
        model_result = pd.concat([model_result, row_result], ignore_index=True)

#Export to output file
print(f"[+] Export model result to {model_result_file}")
model_result.to_excel(model_result_file)

#Saving models
print(f"[+] Saving models")
dump(knn, os.path.join(model_saving_path, f"knn_{fpts_name}_{bits}_bits.joblib"))
dump(rf, os.path.join(model_saving_path, f"rf_{fpts_name}_{bits}_bits.joblib"))
dump(svm, os.path.join(model_saving_path, f"svm_{fpts_name}_{bits}_bits.joblib"))
xgb.save_model(os.path.join(model_saving_path, f"xgb_{fpts_name}_{bits}_bits.json"))


#Export ROC curves
from sklearn.metrics import RocCurveDisplay
fig, ax = plt.subplots(figsize=(6, 6))

RocCurveDisplay.from_estimator(
    estimator=xgb, 
    X=X_Test, 
    y=y_Test,
    name=f"ROC curve for XGBoost",
    color='red',
    ax=ax)

RocCurveDisplay.from_estimator(
    estimator=rf, 
    X=X_Test, 
    y=y_Test,
    name=f"ROC curve for RF",
    color='blue',
    linestyle="dashed",
    ax=ax)

RocCurveDisplay.from_estimator(
    estimator=knn, 
    X=X_Test, 
    y=y_Test,
    name=f"ROC curve for KNN",
    color='gray',
    linestyle='dotted',
    ax=ax)

RocCurveDisplay.from_estimator(
    estimator=svm, 
    X=X_Test, 
    y=y_Test,
    name=f"ROC curve for SVM",
    color='gray',
    linestyle='dotted',
    ax=ax)

plt.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
plt.axis("square")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"{fpts_name}_{bits} bits models ROC curves")
plt.legend()
print(f"[+] Save ROC curves in path: {roc_fig_path}")
plt.savefig(roc_fig_path)
plt.close()