from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

pd.option_context('mode.use_inf_as_na', True)
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
import numpy as np
import textwrap
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import ttest_ind, mannwhitneyu
import textwrap

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Nimbus Sans"
})

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Nimbus Sans",
    'font.size':18
})

fig1A=pd.read_csv('PredictionProbability_modelwith2675features_syntenycorrected_noscaffoldgenes.csv')
fig1A['Label']=fig1A['Label'].replace('Conditional Phenotype','Uncategorized gene models')
fig1A['Label']=fig1A['Label'].replace('Validated Gene Models','Phenotype associated genes')
fig1A['Label']=fig1A['Label'].replace('Non-Validated Gene Models','Loss of function tolerant genes')
fig1A['Label']=fig1A['Label'].replace('Rest of the Gene Models','Uncategorized gene models')

fig1A['Model']="2,675 Features"

# hueorder=['Loss of function tolerant genes','Phenotype associated genes', 'Uncategorized gene models']

fig1D=pd.read_csv('PredictionProbability_modelwith2290features_syntenycorrected_noscaffoldgenes.csv')
fig1D['Label']=fig1D['Label'].replace('Validated Gene Models','Phenotype associated genes')
fig1D['Label']=fig1D['Label'].replace('Non-Validated Gene Models','Loss of function tolerant genes')
fig1D['Label']=fig1D['Label'].replace('Rest of the Gene Models','Uncategorized gene models')
fig1D['Label']=fig1D['Label'].replace('Conditional Phenotype','Uncategorized gene models')

fig1D['Model']="2,290 Features"


fig1F=pd.read_csv('PredictionProbabilityallgenes_top100features_syntenycorrected_noscaffoldgenes.csv')


fig1F['Label']=fig1F['Label'].replace('Validated Gene Models','Phenotype associated genes')
fig1F['Label']=fig1F['Label'].replace('Non-Validated Gene Models','Loss of function tolerant genes')
fig1F['Label']=fig1F['Label'].replace('Rest of the Gene Models','Uncategorized gene models')
fig1F['Label']=fig1F['Label'].replace('Conditional Phenotype','Uncategorized gene models')

fig1F['Model']="Top 100 Features"

finaldf=pd.concat([fig1A, fig1D, fig1F], axis=0)


##TOp panel
#Panel A

xticklabels=['Premature stop \nenriched-common genes','Phenotype\nassociated genes', 'Uncategorized\ngene models']

fig1A=pd.read_csv('all_2675fetaures_ROCcurve_syntenycorrected.csv')
fig1B=pd.read_csv('all_2290fetaures_ROCcurve_syntenycorrcted.csv')
fig1C=pd.read_csv('top100features_ROCcurves_syntenycorrcted.csv')

fpr, tpr, thresholds = roc_curve(fig1A['trulabel'], fig1A['probability'], pos_label=1)
roc_auc = roc_auc_score(fig1A['trulabel'], fig1A['probability'])

fig=plt.figure(figsize=(10,10))

ax=fig.add_subplot(2,1,1)

# gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1])


# Plot the ROC curve 
plt.plot(fpr, tpr, label='2,675 Features (area = %0.2f)' % roc_auc, linewidth=2, color='blue') 

##########
fpr, tpr, thresholds = roc_curve(fig1B['trulabel'], fig1B['probability'], pos_label=1)
roc_auc = roc_auc_score(fig1B['trulabel'], fig1B['probability'])

# fig=plt.figure(figsize=(20,15))

# gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1])

# ax1 = fig.add_subplot(gs[0, 0]) 

# Plot the ROC curve 
plt.plot(fpr, tpr, label='2,290 Features (area = %0.2f)' % roc_auc, linewidth=2, color='green') 

##########
fpr, tpr, thresholds = roc_curve(fig1C['trulabel'], fig1C['probability'], pos_label=1)
roc_auc = roc_auc_score(fig1C['trulabel'], fig1C['probability'])

# fig=plt.figure(figsize=(20,15))

# gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1])

# Plot the ROC curve 
plt.plot(fpr, tpr, label='Top 100 Features (area = %0.2f)' % roc_auc, linewidth=2, color='red') 
# roc curve for tpr = fpr  
plt.plot([0, 1], [0, 1], 'k--', linewidth=2)


# plt.scatter(fpr[optimal_threshold_index], tpr[optimal_threshold_index], marker='s',
#             color='red', label='Optimal Threshold: {:.2f}'.format(optimal_threshold),
#            s=200, zorder=3)
plt.xlabel('False positive rate') 
plt.ylabel('True positive rate') 
# plt.title('ROC Curve') 
plt.legend(loc="lower right", fontsize=14, frameon=False)
# ax.set_yticklabels(ticks)
# ax.set_xticklabels(ticks)


ax=plt.gca()
yticks=ax.get_yticks()
# print(yticks)
ax.set_yticks(yticks[1:])
ax.set_yticklabels([str(round(x,2)) for x in yticks[1:]])

yticks=ax.get_xticks()
ax.set_xticks(yticks[1:])
ax.set_xticklabels([str(round(x,2)) for x in yticks[1:]])

# if i=='exonNum':
#     plt.ylim([-0.02, yticks[-1]])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

plt.ylim(-0.01,1.01)
plt.xlim(-0.01,1.01)

# plt.savefig(f'AllFigure2_ROCcurve.png', dpi=350, bbox_inches='tight')
ax.text(-0.1, 1.05, 'A', transform=ax.transAxes, fontweight='bold', va='top', ha='left')



# fig=plt.figure(figsize=(10,5))

##run bottom cells to generate finaldf

ax=fig.add_subplot(2,1,2)

sns.boxplot(x='Label',y='Predicted_Probability', data=finaldf, hue='Model', palette=['blue', 'green', 'red'])

xticklabels=['Phenotype\nassociated genes', 'Premature stop \nenriched-common genes','Uncategorized\ngene models']

ax=plt.gca()

ax.set_xticks(ax.get_xticks())
ax.set_xticklabels(xticklabels)

ax.set_yticks(ax.get_yticks())
ax.set_yticklabels([round(y,2) for y in ax.get_yticks()])

plt.ylim(0,1)

plt.ylabel('Predicted probability of \nbeing phenotype associated genes')

plt.legend().remove()

plt.xlabel('')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

ax.text(-0.1, 1.05, 'B', transform=ax.transAxes, fontweight='bold', va='top', ha='left')

plt.savefig(f'Figure2_merged_boxplot.svg', dpi=350, bbox_inches='tight')
plt.savefig(f'Figure2_merged_boxplot.png', dpi=350, bbox_inches='tight')
plt.show()