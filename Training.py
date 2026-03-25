import pandas as pd

from sklearn.model_selection import GroupShuffleSplit
from ClassifierModels import *
from sklearn.preprocessing import LabelEncoder
import joblib
from sklearn.metrics import roc_curve, roc_auc_score

##data has been preprocessed with minma for all the gene models.

final_embeddings_family=pd.read_csv('../DataForModeling_withsynteny/Top100preprocessedfeatures_modeling.csv')

splitter = GroupShuffleSplit(test_size=.1, random_state = 0, n_splits=2)
split = splitter.split(final_embeddings_family, groups=final_embeddings_family['family'])
label_encoder = LabelEncoder()
for train_inds, test_inds in split:
    
    train = final_embeddings_family.iloc[train_inds]
    test = final_embeddings_family.iloc[test_inds]

finaltrain_feature=train.drop(['gene', 'family', 'Label'], axis=1)
finaltrain_response=train['Label']
finaltrain_response_ = label_encoder.fit_transform(finaltrain_response)

model=RF(features=finaltrain_feature, response=finaltrain_response_, rescale_type='none')

model.grid_search()

model=model.train_rf()

truelabels=label_encoder.fit_transform(test['Label'].values)
# print(test['Label'])

testfeature =test.drop(['gene', 'family', 'Label'], axis=1)

prob=model.predict_proba(testfeature)

fpr, tpr, thresholds = roc_curve(test['Label'], prob[:,1], pos_label=1)
roc_auc = roc_auc_score(test['Label'].values, prob[:,1]) 
print('ROC_AUC:', roc_auc)

joblib.dump(model, "SupplementalData/DataForModeling_withsynteny/myrandomForestmode_top100features_March13_2026.joblib")