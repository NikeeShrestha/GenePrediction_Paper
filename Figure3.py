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
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import textwrap
import matplotlib.lines as mlines
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import numpy as np


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

fig=plt.figure(figsize=(20,12))

##panel 3A

fig3f=pd.read_csv('Figure3A.csv')

labels=['Uncategorized maize gene models','PAGs-maize', 'PSE-c-maize']


color={'PAGs-maize':'red',
      # 'Old PAGs':'red',
      'Uncategorized maize gene models':'green',
      'PSE-c-maize':'goldenrod'}

ax=fig.add_subplot(2,2,1)
for i in labels:
    subset = fig3f[fig3f['3A'] == i]
#     sns.histplot(subset['Altalelle'], stat='density',
#                  bins=100, alpha=0.8, label=f'{label} (n={len(subset)}/9939)', 
#                  kde=True, color=color[label],line_kws={'linewidth': 2})
    sns.kdeplot(subset['Predicted_Probability'],
                alpha=0.5, fill=True, color=color[i], label=f'{i} (n={len(subset):,})')
plt.legend(loc='upper left', frameon=False)


ax=plt.gca()
yticks=ax.get_yticks()

ax.set_yticks(yticks)

ax.set_yticklabels([round(x,2) for x in yticks])

plt.xlabel(None)

plt.xlabel('Predicted probability of being phenotype associated genes')

# # ax.get_legend().remove()


ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

yticks=ax.get_xticks()

ax.set_xticks(yticks)

ax.set_xticklabels([round(x,2) for x in yticks])

plt.xlim(0,1.01)
xticklabels = [textwrap.fill(label.get_text(), 15) for label in ax.get_xticklabels()]
ax.set_xticklabels(xticklabels)

ax.text(-0.1, 1.05,'A' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')

##3B

nulltempdata1=pd.read_csv('ricegenemodelswithprediction_syntenycorrected.csv')
# nulltempdata

x1=nulltempdata1[nulltempdata1['Label']=='validated']['Predicted_Probability']
x2=nulltempdata1[nulltempdata1['Label']=='rest']['Predicted_Probability']

print('rice:', mannwhitneyu(x1,x2))

nulltempdata1['Label']=nulltempdata1['Label'].replace('validated', 'PAGs-rice')
nulltempdata1['Label']=nulltempdata1['Label'].replace('rest', 'Uncategorized rice gene models')

color={'PAGs-rice':'maroon',
      # 'Old PAGs':'red',
      'Uncategorized rice gene models':'green',
      'Hold-out phenotype associated genes':'red'}


ax=fig.add_subplot(2,2,2)



for i in nulltempdata1['Label'].unique():
    subset = nulltempdata1[nulltempdata1['Label'] == i]
#     sns.histplot(subset['Altalelle'], stat='density',
#                  bins=100, alpha=0.8, label=f'{label} (n={len(subset)}/9939)', 
#                  kde=True, color=color[label],line_kws={'linewidth': 2})
    sns.kdeplot(subset['Predicted_Probability'],
                alpha=0.5, fill=True, color=color[i], label=f'{i} (n={len(subset):,})')
    
# plt.legend(loc='upper left', bbox_to_anchor=(-0.05,1.2), frameon=False)
plt.legend(loc='upper left',frameon=False)
ax=plt.gca()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

# plt.ylim(0,6)
yticks1=ax.get_yticks()
# print(yticks)
ax.set_yticks(yticks1)
ax.set_yticklabels([x for x in yticks1])

yticks=ax.get_xticks()
ax.set_xticks(yticks)
ax.set_xticklabels([str(round(x,2)) for x in yticks])

plt.xlim(0,1)

plt.xlabel('Predicted probability of being phenotype associated genes')
ax.text(-0.1, 1.05,'B' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')


##arabidopsis

nulltempdata2=pd.read_csv('Arabisswithprediction_syntenycorrected.csv')
# df2=pd.read_excel('Figure5/Arabidopsis/PAGS_arabidopsis_class.xlsx')
# df2=df2[df2['Class']=='ESN']
nulltempdata2['Label']=nulltempdata2['Label'].replace('validated', 'PAGs-arabidopsis')
nulltempdata2['Label']=nulltempdata2['Label'].replace('rest', 'Uncategorized arabidopsis gene models')

# nulltempdata2=pd.merge(df1, df2, left_on='GeneID', right_on='Gene', how='left')

# nulltempdata2['Class']=nulltempdata2['Class'].replace('MRP', 'Phenotype associated arabidopsis genes')
# nulltempdata2['Class']=nulltempdata2['Class'].replace('ESN', 'Phenotype associated arabidopsis genes')
# nulltempdata2['Class']=nulltempdata2['Class'].replace('CND', 'Phenotype associated arabidopsis genes')
# nulltempdata2['Class']=nulltempdata2['Class'].replace('CLB', 'Phenotype associated arabidopsis genes')

# nulltempdata2['Class']=nulltempdata2['Class'].fillna('rest')

# nulltempdata2['Class']=nulltempdata2['Class'].replace('rest', 'Uncategorized arabidopsis gene models')

x1=nulltempdata2[nulltempdata2['Label']=='PAGs-arabidopsis']['Predicted_Probability']
x2=nulltempdata2[nulltempdata2['Label']=='Uncategorized arabidopsis gene models']['Predicted_Probability']

print('arabidopsis:', mannwhitneyu(x1,x2))
# nulltempdata
# nulltempdata2['Label']=nulltempdata2['Label'].replace('validated', 'Phenotype associated arabidopsis genes')
# nulltempdata2['Label']=nulltempdata2['Label'].replace('rest', 'Uncategorized arabidopsis gene models')
# nulltempdata2['Label']=nulltempdata2['Class']
# print(nulltempdata1['Label'].unique())

color={'PAGs-arabidopsis':'maroon',
      # 'Old PAGs':'red',
      'Uncategorized arabidopsis gene models':'green',
      'Hold-out phenotype associated genes':'red'}

# fig=plt.figure(figsize=(15,6))
ax=fig.add_subplot(2,2,3)

for i in nulltempdata2['Label'].unique():
    subset = nulltempdata2[nulltempdata2['Label'] == i]
#     sns.histplot(subset['Altalelle'], stat='density',
#                  bins=100, alpha=0.8, label=f'{label} (n={len(subset)}/9939)', 
#                  kde=True, color=color[label],line_kws={'linewidth': 2})
    sns.kdeplot(subset['Predicted_Probability'],
                alpha=0.5, fill=True, color=color[i], label=f'{i} (n={len(subset):,})')
    
# plt.legend(loc='upper left', bbox_to_anchor=(-0.05,1.2), frameon=False)
plt.legend(loc='upper left',frameon=False)
ax=plt.gca()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

# yticks=ax.get_yticks()
# print(yticks)

yticks1=ax.get_yticks()
ax.set_yticks(yticks1)
ax.set_yticklabels([x for x in yticks1])


# plt.xlim(0,6)
yticks=ax.get_xticks()
ax.set_xticks(yticks)
ax.set_xticklabels([str(round(x,2)) for x in yticks])

plt.xlabel('Predicted probability of being phenotype associated genes')
ax.text(-0.1, 1.05,'C' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')
# plt.savefig('Figure5/Arabidopsis_rice.png', dpi=350, bbox_inches='tight')



##panel D

ax=fig.add_subplot(2,2,4)

fig3f=pd.read_csv('Figure3D.csv')


labels=['Uncategorized maize gene models','New PAGs-maize', 'Hold-out PAGs-maize',  'Hold-out PSE-c-maize', 'PSE-r-maize']

color={'Uncategorized maize gene models':'green',
      # 'Old PAGs':'red',
      'New PAGs-maize':'maroon',
      'Hold-out PAGs-maize':'red',
      'Hold-out PSE-c-maize': 'goldenrod',
      'PSE-r-maize': 'plum'}

for i in labels:
    subset = fig3f[fig3f['list'] == i]
#     sns.histplot(subset['Altalelle'], stat='density',
#                  bins=100, alpha=0.8, label=f'{label} (n={len(subset)}/9939)', 
#                  kde=True, color=color[label],line_kws={'linewidth': 2})
    sns.kdeplot(subset['Predicted_Probability'],
                alpha=0.5, fill=True, color=color[i], label=f'{i} (n={len(subset):,})')
plt.legend(loc='upper left', frameon=False)


ax=plt.gca()
yticks=ax.get_yticks()

ax.set_yticks(yticks)

ax.set_yticklabels([round(x,2) for x in yticks])

plt.xlabel(None)

plt.xlabel('Predicted probability of being phenotype associated genes')

# # ax.get_legend().remove()


ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

yticks=ax.get_xticks()

ax.set_xticks(yticks)

ax.set_xticklabels([round(x,2) for x in yticks])

plt.xlim(0,1.01)
xticklabels = [textwrap.fill(label.get_text(), 15) for label in ax.get_xticklabels()]
ax.set_xticklabels(xticklabels)

ax.text(-0.1, 1.05,'D' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')

plt.savefig(f'Figure3_rice_tair.svg', dpi=350, bbox_inches='tight')
plt.savefig(f'Figure3_rice_tair.png', dpi=350, bbox_inches='tight')

plt.show()