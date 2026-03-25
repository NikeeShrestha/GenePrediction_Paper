import pandas as pd
import numpy as np
from ClassifierModels import *
import joblib
# import importlib

# importlib.reload(ClassifierModels)

##Arabidopsis gene prediction

rawfeatures=pd.read_csv("SupplementalData/DataForModeling_withsynteny/PreprocessedArabidopsis_top99features.csv")

loaded_rf = joblib.load("SupplementalData/DataForModeling_withsynteny/myrandomForestmode_top99features_March17_2026_rice_tair.joblib")

predictiondata = pd.DataFrame({
    'GeneID':rawfeatures.gene_id.values,
    'Predicted_Probability': loaded_rf.predict_proba(np.array(rawfeatures.drop(['gene_id'], axis=1)))[:, 1]
})

predictiondata.to_csv('Arabidopsisgene_withprediction.csv', index=False)

##Rice gene Prediction

rawfeatures=pd.read_csv("SupplementalData/DataForModeling_withsynteny/Preprocessedrice_top99features.csv")

loaded_rf = joblib.load("SupplementalData/DataForModeling_withsynteny/myrandomForestmode_top99features_March17_2026_rice_tair.joblib")

predictiondata = pd.DataFrame({
    'GeneID':rawfeatures.gene_id.values,
    'Predicted_Probability': loaded_rf.predict_proba(np.array(rawfeatures.drop(['gene_id'], axis=1)))[:, 1]
})

predictiondata.to_csv('Ricegene_withprediction.csv', index=False)