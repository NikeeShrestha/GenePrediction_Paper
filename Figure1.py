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

from scipy.stats import mannwhitneyu

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

##figure 1 A

df=pd.read_csv('../Nullgenes_prematurestopcodon.csv')
nullgenes=(df[df['Altalelle']>=0.25]['gene'])
fig = plt.figure(figsize=(15, 8))

# Set up GridSpec: 2 rows, 6 columns
gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1.5])

# Top two panels (each span 3 columns)
ax1 = fig.add_subplot(gs[0, 0:2]) 

sns.histplot(df['Altalelle'], alpha=0.5, color='purple')
plt.xlim(0,1)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
for spine in ax1.spines.values():
    spine.set_linewidth(4)
ax1.tick_params(width=4, length=8)

plt.ylim(0,250)
# plt.ylim(0.3,0.5)
yticks=ax1.get_yticks()
# print(yticks)
xticks=ax1.get_xticks()
ax1.set_yticks(yticks)
ax1.set_yticklabels([str(int(x)) for x in yticks])


xticks=ax1.get_xticks()
ax1.set_xticks(xticks)
ax1.set_xticklabels([str(round(x,2)) for x in xticks])

plt.axvline(x=0.25, color='r', linestyle='--', linewidth=2)

plt.xlabel('Alternate allele frequency')

# plt.text(0.26, 200, f'n = {len(set(nullgenes))}  gene models')
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontweight='bold', va='top', ha='left')

##panel B
fig1B=pd.read_csv('Figure1B.csv')
ax2 = fig.add_subplot(gs[0, 2:4])

colors={'Validated Gene Models':'Reds', 'Non-Validated Gene Models':'Blues',
        'Rest of the Gene Models':'green','Conditional Phenotype':'brown'}

line_colors = {
    'Validated Gene Models': 'darkred',
    'Non-Validated Gene Models': 'darkblue',
    'Rest of the Gene Models': 'darkgreen',
    'Conditional Phenotype': 'saddlebrown'
}

custom_legend_handles=[]
# labels=['Rest of the Gene Models', 'Conditional Phenotype', 'Validated Gene Models', 'Non-Validated Gene Models']
labels=['Validated Gene Models','Non-Validated Gene Models']

customlabels={'Validated Gene Models': 'Phenotype associated genes (PAGs)','Non-Validated Gene Models': 'Premature stop codon enriched-common genes (PSE-c)'}
for i in labels:
    sns.kdeplot(data=fig1B[fig1B['Label']==i],x='PC1',y='PC2',
                    fill=True,cmap=colors[i], alpha=0.8,label=i,levels=7,
                    linewidths=0, ax=ax2)
    sns.kdeplot(
    data=fig1B[fig1B['Label'] ==i],
    x='PC1', y='PC2',
    fill=False,
    color=line_colors[i],
    linewidths=2,
    levels=7,
    ax=ax2
    )
    
    facecolor_with_alpha = to_rgba(line_colors[i], alpha=0.5)

    custom_legend_handles.append(Line2D([0], [0], marker='s', color='w', markerfacecolor=facecolor_with_alpha, markersize=20, label=customlabels[i],alpha=1))


# plt.gca().set_aspect('auto')

yticks=ax2.get_yticks()

ax2.set_yticks(yticks)
ax2.set_yticklabels([str(round(x,2)) for x in yticks])


yticks=ax2.get_xticks()

ax2.set_xticks(yticks)
ax2.set_xticklabels([str(round(x,2)) for x in yticks])

ax2.set_xlim(-5,5)
ax2.set_ylim(-5,10)
# if i=='exonNum':
#     plt.ylim([-0.02, yticks[-1]])
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
for spine in ax2.spines.values():
    spine.set_linewidth(4)
ax2.tick_params(width=4, length=8)

labels = [h.get_label() for h in custom_legend_handles]

# 2. Define the new desired order
# For example, if you want to reorder manually:
new_order = [1, 0]  # <-- whatever order you want (indices)

# 3. Reorder both handles and labels
reordered_handles = [custom_legend_handles[i] for i in new_order]
reordered_labels = [labels[i] for i in new_order]

# 4. Now plot legend with reordered handles and labels
ax2.legend(handles=reordered_handles, labels=reordered_labels,
           loc='upper right', frameon=False)

# ax2.legend(handles=custom_legend_handles,loc='upper right',frameon=False)
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontweight='bold', va='top', ha='left')

##panelC
final_embeddings_family=pd.read_csv('Figure1C.csv', index_col=0)
final_embeddings_family=final_embeddings_family.sample(frac=1)
labels=[ 'Non-Validated Gene Models','Validated Gene Models']
final_embeddings_family = final_embeddings_family.rename(columns={'exonNum': 'Number of exons', 'UTR3length': '3\'UTR length', 
                                                                 'isoforms': 'Number of isoforms', 'GC': 'GC content'})
final_embeddings_familytrait=final_embeddings_family.columns[1:5]
num=0

plotnum=['C','D', 'E', 'F']

for i in final_embeddings_familytrait:
    ax=fig.add_subplot(gs[1, num])
    tempdf=final_embeddings_family[['Label', i]]
    # tempdf[i]=-np.log10(tempdf[i]+0.0001)
    tempdf = tempdf.copy()  # Make a proper copy if tempdf is a slice
    group1=tempdf[tempdf['Label']=='Validated Gene Models'][i]
    group2=tempdf[tempdf['Label']=='Non-Validated Gene Models'][i]
    pval=mannwhitneyu(group1, group2)[1]

    # tempdf.loc[:, i] = np.log10(tempdf[i].clip(lower=1))
    # tempdf.loc[:, i]=-np.log10(tempdf[i]+0.0001)
    sns.boxplot(data=tempdf, x='Label', y=i, palette=['blue', 'red'],
                linewidth=2,width=0.2, boxprops=dict(edgecolor="black"),
               whiskerprops=dict(color="black"),medianprops=dict(color="black"),capprops=dict(color="black"),
               flierprops=dict(markerfacecolor="black", markeredgecolor="black"), order=labels)
    
    # plt.ylabel(f'$log10$({i})')
    plt.ylabel(f'{i}')
    plt.xlabel(None)
    
    newticklabels = ['PSE-c','PAGs']
    wrapped_labels = [textwrap.fill(label, 15) for label in newticklabels] 
    
    plt.xticks(range(len(labels)), wrapped_labels)
    
    # print(yticks)
   
 
    ax=plt.gca()
    yticks=ax.get_yticks()
   
    if pval>0.05:
        print(pval)
        pvalue='n.s'
        text= f'{pvalue}'
    if pval<0.05 and pval>0.01:
        print(pval)
        pvalue='*'
        text= fr'{pvalue}'
    if 0.001 < pval < 0.01:
        print(pval)
        pvalue='**'
        text=f'{pvalue}'
    if pval<0.001:
        print(pval)
        pvalue='***'
        # text=r'$p\textrm{-}value < 0.001$'
        # r'$p\textrm{-}value < 0.001$'
        text=f'{pvalue}'

    
    if i=='Number of exons':
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(int(x)) for x in yticks])
        plt.ylim([-0.2, 60])
        plt.plot([0,1], [55,55],color='black', linewidth=2)
        plt.text(0.4, 56, text, fontsize=18)
    if i=='3\'UTR length':
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(int(x)) for x in yticks])
        plt.ylim([-100, yticks[-1]])
        plt.plot([0,1], [2700,2700],color='black', linewidth=2)
        # print(text)
        plt.text(0.4, 2750, text, fontsize=18)
    if i=='Number of isoforms':
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(int(x)) for x in yticks])
        plt.ylim([0.5, yticks[-1]])
        plt.plot([0,1], [9.5,9.5],color='black', linewidth=2)
        plt.text(0.4, 9.7, text, fontsize=18)
        
    if i=='GC content':
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{round(x,2)}\%' for x in yticks])
        plt.ylim([0.35, yticks[-1]])
        plt.plot([0,1], [0.85,0.85],color='black', linewidth=2)
        plt.text(0.4, 0.87, text, fontsize=18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(4)
    ax.tick_params(width=4, length=8)
    
        
    ax.text(-0.3, 1.05, plotnum[num], transform=ax.transAxes, fontweight='bold', va='top', ha='left')
    num+=1
    

plt.savefig(f'Figure1.svg', dpi=350, bbox_inches='tight')
plt.show()

