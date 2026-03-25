import pandas as pd
import numpy as np
from ClassifierModels import *
import joblib
# import importlib

# importlib.reload(ClassifierModels)

rawfeatures=pd.read_csv("SupplementalData/DataForModeling_withsynteny/Preprocessedwithlabelsandfamily_2290_finalallgenesnopangenecount_39034genes_april29.csv")

loaded_rf = joblib.load("SupplementalData/DataForModeling_withsynteny/myrandomForestmode_2290_allfeatures_March13_2026.joblib")

importances=loaded_rf.feature_importances_
feature_imp_df = pd.DataFrame({'Feature': rawfeatures.columns[1:-3], 'Gini Importance': importances}).sort_values('Gini Importance', ascending=False) 

feature_imp_df.iloc[0:100,].to_csv('Allfeaturesimportance_2290features_syntenycorrected.csv', index=False)

predictiondata = pd.DataFrame({
    'GeneID':rawfeatures.gene.values,
    'Predicted_Probability': loaded_rf.predict_proba(np.array(rawfeatures.drop(['gene', 'Label', 'Part','family'], axis=1)))[:, 1],
    
})

predictiondata.to_csv('Allgenes_predictedprobabilities_2290features_syntenycorrected.csv', index=False)