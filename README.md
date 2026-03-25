# Gene Prediction – Code Description

Scripts for training, predicting, evaluating, and visualizing a Random Forest classifier that identifies phenotype-associated genes (PAGs) vs. premature stop codon-enriched common genes (PSE-c) in maize.

---

## Core Model Scripts

### `ClassifierModels.py`
Defines reusable classes for machine learning workflows.

- **`base`** – Parent class with feature preprocessing (standardization or min-max scaling), prediction, and K-Fold cross-validation with optional Observed vs. Predicted visualization.
- **`RF(base)`** – Random Forest subclass supporting hyperparameter grid search (`RandomizedSearchCV`), model training, K-Fold cross-validation, and Gini-based feature importance.

### `Training.py`
Trains the Random Forest classifier on the top 100 features. Splits data by gene family, runs grid search, trains the final model, evaluates using ROC-AUC, and saves the trained model as a `.joblib` file.

---

## Prediction Scripts

Load a pre-trained model and generate predicted probabilities (of being a PAG) for all gene models. Each script also exports feature importances ranked by Gini score. These scripts are related to Supplemental Dataset S3.

| Script | Feature Set | Output Files |
|---|---|---|
| `Prediction_Top100Features.py` | Top 100 | `Allgenes_predictedprobabilities_Top100features_syntenycorrected.csv`, `Allfeaturesimportance_100features_syntenycorrected.csv` |
| `Prediction_2290Features.py` | 2,290 features | `Allgenes_predictedprobabilities_2290features_syntenycorrected.csv`, `Allfeaturesimportance_2290features_syntenycorrected.csv` |
| `Prediction_2675Features.py` | 2,675 features | `Allgenes_predictedprobabilities_2675features_syntenycorrected.csv`, `Allfeaturesimportance_2675features_syntenycorrected.csv` |

### `Prediction_rice_arabidopsis.py`
Applies the maize-trained model to predict gene function in **rice** and **Arabidopsis**, producing per-species prediction CSVs. This script is related to Supplemental Dataset S5.

---

## Statistical Analysis

### `CohensD.py`
Performs pairwise feature-level comparisons across gene categories using **Mann-Whitney U test** and **Cohen's d** effect size. Runs three comparisons:
1. PSE-c vs. PAG
2. nonPAG_nonPSE-c vs. PAG
3. nonPAG_nonPSE-c vs. PSE-c

Each output CSV contains per-feature means, standard deviations, Mann-Whitney statistic, p-value, and Cohen's d. This script is related to Supplemental Dataset S2.

---

## Figure Scripts

### `Figure1.py`
Characterizes PAG vs. PSE-c gene models across six panels:
- **A** – Histogram of alternate allele frequency; threshold line marks PSE-c cutoff (≥0.25).
- **B** – 2D KDE of PCA scores (PC1 vs. PC2) showing separation between PAGs and PSE-c.
- **C–F** – Boxplots comparing number of exons, 3′ UTR length, number of isoforms, and GC content between gene categories, with Mann-Whitney significance annotations.

**Output**: `Figure1.svg`

### `Figure2.py`
Evaluates model performance across the three feature sets:
- **A** – Overlaid ROC curves for all three models with AUC.
- **B** – Boxplots of predicted probabilities by gene category across models.

**Outputs**: `Figure2_merged_boxplot.svg`, `Figure2_merged_boxplot.png`

### `Figure3.py`
Assesses cross-species transferability of the maize model:
- **A** – KDE of predicted probabilities for maize gene categories.
- **B** – Rice PAGs vs. uncategorized gene models.
- **C** – Arabidopsis PAGs vs. uncategorized gene models.
- **D** – Novel/hold-out PAG categories in maize (new PAGs, hold-out PAGs, hold-out PSE-c, PSE-r).

**Outputs**: `Figure3_rice_tair.svg`, `Figure3_rice_tair.png`

### `Figure4.py`
Compares functional genomic properties of top vs. bottom quartile predicted gene models:
- **A** – Histogram of all predicted probabilities colored by quartile.
- **B** – KDE of premature stop codon frequency by quartile.
- **C** – Variant frequency per kilobase (high/moderate/low effect) by quartile.
- **D** – Stacked bar chart of variant effect proportions by quartile.
- **E** – Average absolute zero-shot scores across genomic regions (gene body, CDS, 5′UTR, 3′UTR, intron) with significance annotations.

**Outputs**: `Figure4.svg`, `Figure4.png`

---

## Dependencies

| Package | Purpose |
|---|---|
| `pandas`, `numpy` | Data manipulation |
| `scikit-learn` | Random Forest, cross-validation, preprocessing |
| `scipy` | Statistical tests (Mann-Whitney U, Spearman/Pearson r) |
| `matplotlib`, `seaborn` | Plotting and figure generation |
| `joblib` | Saving/loading trained models |
