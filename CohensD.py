import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

final_embeddings_family=pd.read_csv("Continuousfeatures_forcohensd.csv")

features=final_embeddings_family.columns[:-1]
features=features[1:]

##PSE-c vs PAG
mannwhitneytest=[]
for i in features:
    tempdf=final_embeddings_family[['Label',i]]
    group1=final_embeddings_family[final_embeddings_family['Label']=='PSE-c'][i]
    group2=final_embeddings_family[final_embeddings_family['Label']=='PAG'][i]
    sd_y=np.std(group1,ddof=1)
    sd_n=np.std(group2,ddof=1)
    mean_y=np.mean(group1)
    mean_n=np.mean(group2)
    pval=mannwhitneyu(group1, group2)[1]
    stat=mannwhitneyu(group1, group2)[0]
    pooled_sd = np.sqrt(((len(group1) - 1) * sd_y ** 2 + (len(group2) - 1) * sd_n ** 2) / 
                        (len(group1) + len(group2) - 2))
    cohens_d = (mean_y - mean_n) / pooled_sd
    
    
    mannwhitneytest.append({
        'feature':i,
        'mean_PAG':mean_y,
        'mean_PSE-c':mean_n,
        'sd_PAG':sd_y,
        'sd_PSE-c':sd_n,
        'whitney_pval':pval,
        'whitney_stat':stat,
        'cohens_d':cohens_d
    })

mannwhitneytest=pd.DataFrame.from_dict(mannwhitneytest)
mannwhitneytest.to_csv('Mannwhitneytestforallfeatures_PSE-cvsPAG.csv', index=False)


##nonPAG_nonPSE-c vs PAG
mannwhitneytest=[]
for i in features:
    tempdf=final_embeddings_family[['Label',i]]
    group1=final_embeddings_family[final_embeddings_family['Label']=='nonPAG_nonPSE-c'][i]
    group2=final_embeddings_family[final_embeddings_family['Label']=='PAG'][i]
    sd_y=np.std(group1,ddof=1)
    sd_n=np.std(group2,ddof=1)
    mean_y=np.mean(group1)
    mean_n=np.mean(group2)
    pval=mannwhitneyu(group1, group2)[1]
    stat=mannwhitneyu(group1, group2)[0]
    pooled_sd = np.sqrt(((len(group1) - 1) * sd_y ** 2 + (len(group2) - 1) * sd_n ** 2) / 
                        (len(group1) + len(group2) - 2))
    cohens_d = (mean_y - mean_n) / pooled_sd
    
    
    mannwhitneytest.append({
        'feature':i,
        'mean_PAG':mean_y,
        'mean_nonPAG_nonPSE-c':mean_n,
        'sd_PAG':sd_y,
        'sd_nonPAG_nonPSE-c':sd_n,
        'whitney_pval':pval,
        'whitney_stat':stat,
        'cohens_d':cohens_d
    })

mannwhitneytest=pd.DataFrame.from_dict(mannwhitneytest)
mannwhitneytest.to_csv('Mannwhitneytestforallfeatures_nonPAG_nonPSE-cvsPAG.csv', index=False)

##nonPAG_nonPSE-c vs PSE-c
mannwhitneytest=[]
for i in features:
    tempdf=final_embeddings_family[['Label',i]]
    group1=final_embeddings_family[final_embeddings_family['Label']=='nonPAG_nonPSE-c'][i]
    group2=final_embeddings_family[final_embeddings_family['Label']=='PSE-c'][i]
    sd_y=np.std(group1,ddof=1)
    sd_n=np.std(group2,ddof=1)
    mean_y=np.mean(group1)
    mean_n=np.mean(group2)
    pval=mannwhitneyu(group1, group2)[1]
    stat=mannwhitneyu(group1, group2)[0]
    pooled_sd = np.sqrt(((len(group1) - 1) * sd_y ** 2 + (len(group2) - 1) * sd_n ** 2) / 
                        (len(group1) + len(group2) - 2))
    cohens_d = (mean_y - mean_n) / pooled_sd
    
    
    mannwhitneytest.append({
        'feature':i,
        'mean_nonPAG_nonPSE-c':mean_y,
        'mean_PSE-c':mean_n,
        'sd_nonPAG_nonPSE-c':sd_y,
        'sd_PSE-c':sd_n,
        'whitney_pval':pval,
        'whitney_stat':stat,
        'cohens_d':cohens_d
    })

mannwhitneytest=pd.DataFrame.from_dict(mannwhitneytest)
mannwhitneytest.to_csv('Mannwhitneytestforallfeatures_nonPAG_nonPSE-cvsPSE-c.csv', index=False)