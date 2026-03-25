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
    'font.size':24
})

##figure 3

genexticklabels=['Bottom quartile gene models', 'Top quartile gene models']
#panel A

fig3A=pd.read_csv('../Figure2/PredictionProbability_100features_allgenes_nopangenecount_allgenes_top100features_nogeneexpressionsyntenycorrected.csv')
q1=fig3A['Predicted_Probability'].quantile(0.25)
q3=fig3A['Predicted_Probability'].quantile(0.75)

fig=plt.figure(figsize=(30, 15))

# fig = plt.figure(figsize=(10, 6))
gs = fig.add_gridspec(2, 6)

# gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1.5])

# ax=fig.add_subplot(3,2,1)
ax=fig.add_subplot(gs[0, :3]) 


ax = sns.histplot(fig3A['Predicted_Probability'], stat='density', bins=50, alpha=0.2, color='gray')

# Create a normalization based on the range of values
norm = mcolors.Normalize(vmin=fig3A['Predicted_Probability'].min(), 
                         vmax=fig3A['Predicted_Probability'].max())

# Choose a colormap
cmap = plt.get_cmap('viridis')

# # Loop over each bar (patch) in the histogram
# for patch in ax.patches:
#     # Compute the bin center
#     bin_center = patch.get_x() + patch.get_width() / 2
#     # Set the facecolor based on the bin center
#     patch.set_facecolor(cmap(norm(bin_center)))
      
for patch in ax.patches:
    bin_center = patch.get_x() + patch.get_width() / 2
    if bin_center <= q1:
        color = 'mediumorchid'
    elif bin_center < q3:
        color = 'gray'
    else:
        color = 'yellowgreen'
    patch.set_facecolor(color)
    
plt.xlim(0,1)

plt.axvline(x=q1, ls='--', color='black', linewidth=3)
plt.axvline(x=q3, ls='--', color='black', linewidth=3)

myax=plt.gca()
myax.spines['right'].set_visible(False)
myax.spines['top'].set_visible(False)
for spine in myax.spines.values():
    spine.set_linewidth(4)
myax.tick_params(width=4, length=8)

# plt.ylim(0,80)
# plt.ylim(0.3,0.5)
yticks=myax.get_yticks()
# print(yticks)
xticks=myax.get_xticks()
myax.set_yticks(yticks)
myax.set_yticklabels([str(round(x,2)) for x in yticks])

myax.set_xticks(xticks)
myax.set_xticklabels([str(round(x,2)) for x in xticks])


plt.text(0.05, 1.6, 'Bottom quartile\ngene models')
plt.text(0.8, 1.6, 'Top quartile\ngene models')

plt.xlabel('Predicted probability of being phenotype associated genes')
ax.text(-0.12, 1.05,'A' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')


#figure3B

fig3B=pd.read_csv('Twocategoriesgenes_frequencymarker_witheffect_onlyprimarytranscriptcds_syntenycorrected_2290_top100.csv')
fig3B['Category']=fig3B['Category'].replace('High Probability', 'Top quartile gene models')
fig3B['Category']=fig3B['Category'].replace('Low Probability', 'Bottom quartile gene models')

stopgaineddf=pd.read_csv('../Nullgenes_prematurestopcodon.csv')


size=fig3B[fig3B['Category']=='Top quartile predicted genes'].shape[0]
finaldf=fig3B
finaldf=pd.merge(stopgaineddf,finaldf, left_on='gene', right_on='GeneID')
# # sns.histplot(data=finaldf, x='Altalelle',stat='density',bins=100, alpha=0.2, hue='prob')
fig3bfortest=finaldf
# ax=fig.add_subplot(3,2,2)

ax=fig.add_subplot(gs[0, 3:])
# labels=['High Probability', 'Low Probability']
# color={'High Probability':'red', 'Low Probability':'blue'}
size={'Top quartile gene models':fig3B[fig3B['Category']=='Top quartile gene models'].shape[0], 
      'Bottom quartile gene models':fig3B[fig3B['Category']=='Bottom quartile gene models'].shape[0]}

sns.kdeplot(data=finaldf, x='Altalelle', hue='Category',
            alpha=0.9, fill=True, hue_order=['Top quartile gene models', 'Bottom quartile gene models'], palette=['yellowgreen', 'mediumorchid'])

ax=plt.gca()

subset1=finaldf[finaldf['Category'] == 'Bottom quartile gene models']
subset2=finaldf[finaldf['Category'] == 'Top quartile gene models']


custom_lines = [
      mlines.Line2D([], [], color='mediumorchid', linestyle='-', linewidth=2,
              label=f"Bottom quartile gene models(n={len(subset1)}/{size['Bottom quartile gene models']})"),
      mlines.Line2D([], [], color='yellowgreen', linestyle='-', linewidth=2,
              label=f"Top quartile gene models(n={len(subset2)}/{size['Top quartile gene models']})")
  ]

ax.legend(handles=custom_lines,  loc='upper right', frameon=False)

# ax=plt.gca()
yticks=ax.get_yticks()
# print(yticks)
ax.set_yticks(yticks)
ax.set_yticklabels([str(round(x,2)) for x in yticks])



yticks=ax.get_xticks()
ax.set_xticks(yticks)
ax.set_xticklabels([str(round(x,2)) for x in yticks])

# if i=='exonNum':
#     plt.ylim([-0.02, yticks[-1]])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)
plt.xlim(0,1)

plt.xlabel('Frequency of alternate allele causing PTC')
ax.text(-0.12, 1.05,'B' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')


## panel C
fig3D=pd.read_csv('markerfrequencytwocategories_noscaffgenes_syntenycorrected_2290_top100.csv')
# ax=fig.add_subplot(3,2,3)

ax=fig.add_subplot(gs[1, 0:2])
fig3D['Category'] = pd.Categorical(
    fig3D['Category'],
    categories=['Low Probability', 'High Probability'],
    ordered=True
)

custom_palette = {'high': 'darkslategray', 'moderate': 'darkturquoise', 'low': 'lightblue'}

sns.lineplot(data=fig3D, x='Category',
             y='value', hue='variable', err_style='bars', errorbar=('ci', 95),
             palette=custom_palette, err_kws={'capsize':5,'capthick':3, 'elinewidth':2, },
             linewidth=3, markersize=8,estimator=np.mean, ax=ax)
ax.set_xlim(-0.1, 1.1) 
# ax=plt.xticks([0,1],  rotation=0)
ax1=plt.gca()
yticks=ax1.get_yticks()

ax1.set_yticks(yticks)

ax1.set_yticklabels([round(x,2) for x in yticks])

ax1.set_xticklabels(['Bottom quartile\ngene models', 'Top quartile\ngene models'])

plt.xlabel(None)

plt.ylabel('Variant frequency per kilobase')

plt.legend(title='variant effect', loc='best', frameon=False)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
for spine in ax1.spines.values():
    spine.set_linewidth(4)
ax1.tick_params(width=4, length=8)
plt.ylim(0,7)
ax.text(-0.12, 1.05,'C' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')

##figure 3D
# ax=fig.add_subplot(3,2,4)

ax=fig.add_subplot(gs[1, 2:4])

fig3C=pd.read_csv('peercentmarkereffect_noscaffgenes_2290_top100_syntenycorrected_cds_primarytranscript.csv', index_col=0)
desired_order = ["high", "moderate", "low"] 
fig3C = fig3C.loc[['Low Probability', 'High Probability']]

fig3C.plot(kind='bar', stacked=True,
                color={'high': 'darkslategray', 'moderate': 'darkturquoise', 'low': 'lightblue' },
                width=0.3, ax=ax)
ax=plt.ylabel('Proportion of variants')
# plt.xlabel(''angle=0)
ax=plt.xticks([0,1],  rotation=0)
ax=plt.legend(title='variant effect', loc='center',frameon=False)

ax=plt.gca()
# ax.set_xticklabels(['Bottom quartile\ngene models', 'Top quartile\ngene models'])

ax=plt.gca()
yticks=ax.get_yticks()

ax.set_yticks(yticks[:-1])

ax.set_yticklabels([str(f'{int(x)}\%') for x in yticks[:-1]])

plt.xlabel(None)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

ax.set_xticks([])

# Define custom label positions and box styles
# xtick_labels = ['Low probability', 'High probability']
xtick_labels=genexticklabels
xtick_positions = [bar.get_x() + bar.get_width()/2 for bar in ax.containers[0]]

# Define box colors
xtick_colors = {
    'Bottom quartile gene models': 'mediumorchid',
    'Top quartile gene models': 'yellowgreen'
}

# Add text with colored box for each tick
for xpos, label in zip(xtick_positions, xtick_labels):
    ax.text(
        xpos, -3,  # y slightly below the x-axis
        label,
        ha='center', va='top',
        weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor=xtick_colors[label], edgecolor='black', linewidth=0.4)
    )
# ax.set_xticklabels(genexticklabels)
ax.text(-0.12, 1.05,'D' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')   

##panel F

df=pd.read_csv('genesnpscoremean2026.csv', index_col=0)
# quartile=fig3B[['GeneID','Category']]
# df=pd.merge(df,quartile, on='GeneID' )
# df=df.fillna(0)

df['Category']=df['Category'].replace('Low Probability','Bottom quartile gene models')
df['Category']=df['Category'].replace('High Probability','Top quartile gene models')

df_longer=pd.melt(df, value_vars=df.columns[[1,4,5,6,7,8]], value_name='shotscore', id_vars=['GeneID', 'Category'])
# ax=fig.add_subplot(3,2,5)

ax=fig.add_subplot(gs[1, 4:6])
# sns.barplot(x="day", y="total_bill", data=df, )
# plt.show()
sns.barplot(x="variable", y="shotscore", data=df_longer, hue='Category',errorbar=('se', 1),
           palette=['mediumorchid', 'yellowgreen'], hue_order=['Bottom quartile gene models', 'Top quartile gene models'],
           errwidth=3, errcolor='black', capsize=0.05, ax=ax)

ax=plt.gca()
yticks=ax.get_yticks()

ax.set_yticks(yticks)

ax.set_yticklabels([round(x,2) for x in yticks])

plt.xlabel(None)

plt.ylabel('Average absolute zero shot score')

plt.legend(title=None, frameon=False, loc='upper right')

# ax.get_legend().remove()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(4)
ax.tick_params(width=4, length=8)

# xticklabels=ax.get_xticklabels()

# xticklabels=['Gene\nbody'] + [x for x in xticklabels[1:]]
xticklabels=['gene\nbody', 'CDS', '5\'UTR', '3\'UTR', 'intron']

ax.set_xticklabels([x for x in xticklabels])
ax.text(-0.12, 1.05,'E' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')

plt.plot([-0.2,0.3], [df_longer[df_longer['variable']=='genenmean']['shotscore'].quantile(0.75)-0.02,
                      df_longer[df_longer['variable']=='genenmean']['shotscore'].quantile(0.75)-0.02], color='black')
plt.text(-0.18, df_longer[df_longer['variable']=='genenmean']['shotscore'].quantile(0.75)+0.02-0.02, '***', color='black')


plt.plot([0.7,1.2], [df_longer[df_longer['variable']=='CDSmean']['shotscore'].quantile(0.75)-0.02,
                      df_longer[df_longer['variable']=='CDSmean']['shotscore'].quantile(0.75)-0.02], color='black')
plt.text(0.75, df_longer[df_longer['variable']=='CDSmean']['shotscore'].quantile(0.75)+0.02-0.02, '***', color='black')


plt.plot([1.8,2.3], [df_longer[df_longer['variable']=='mUTR5mean']['shotscore'].quantile(0.75)+0.1,
                      df_longer[df_longer['variable']=='mUTR5mean']['shotscore'].quantile(0.75)+0.1], color='black')
plt.text(1.9, df_longer[df_longer['variable']=='mUTR5mean']['shotscore'].quantile(0.75)+0.15, '**', color='black')



plt.plot([2.8,3.3], [df_longer[df_longer['variable']=='mUTR3mean']['shotscore'].quantile(0.75)+0.25,
                      df_longer[df_longer['variable']=='mUTR3mean']['shotscore'].quantile(0.75)+0.25], color='black')
plt.text(2.9, df_longer[df_longer['variable']=='mUTR3mean']['shotscore'].quantile(0.75)+0.27, '**', color='black')


plt.plot([3.7,4.2], [df_longer[df_longer['variable']=='intronmean']['shotscore'].quantile(0.75),
                      df_longer[df_longer['variable']=='intronmean']['shotscore'].quantile(0.75)], color='black')
plt.text(3.8, df_longer[df_longer['variable']=='intronmean']['shotscore'].quantile(0.75)+0.02, '**', color='black')
ax.text(-0.12, 1.05,'E' ,transform=ax.transAxes, fontweight='bold', va='top', ha='left')



plt.savefig(f'Figure3_March19.svg', dpi=350, bbox_inches='tight')
plt.savefig(f'Figure3_March19.png', dpi=350, bbox_inches='tight')

plt.show()