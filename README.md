# Multi-omic preterm birth machine learning analysis

## Working title
Immunological Maladaptation as a Predictor of Spontaneous Preterm Birth in Human Pregnancies

## Abstract
Dysregulated maternal immune adaptation during pregnancy plays a central role in the pathogenesis of spontaneous preterm birth (sPTB). However, predictive models of sPTB based on maternal immune features, ideally before clinical symptoms arise, remain limited. In a nested case-control study embedded within a population-based, low-risk pregnancy cohort, we decoded a pattern of abnormal immune response in mothers' blood that precedes sPTB by weeks to months. Prominent features include heightened sensitivity to adrenergic signals within myeloid and T cell populations in early pregnancy, followed by increased production of pro-inflammatory cytokines in the third trimester. CD4+ T cells exhibited gene expression patterns indicative of a Th17-skewed, neuroactive protein–responsive phenotype. Our study provides a multi-omics resource and a conceptual framework for identifying individuals at increased risk for sPTB with broad translational implications for advancing targeted preventive measures.


## Description
This repository contains code to perform machine learning analyses for the prediction of preterm birth. 

Blood samples were drawn from pregnant women at three time points during gestation (T1, T2, T3) and analysed with i) mass cytometry to profile immune cells, and ii) an aptamer-based platform to profile the proteome. In total, seven sets of measurements were derived. These are:

1. Frequency data (freq) representing the abundance of each of 27 immune cell populations; 
2. Intracellular cytokine expression data (ICS) at baseline or in response to stimulation; 
3. Endogenous phosphosignaling responses (US);
4. Phosphosignaling response after Isoproterenol (ISO) stimulation;
5. Phosphosignaling response after lipopolysaccharide (LPS) stimulation;
6. Phosphosignaling response after simultaneous PI and IL-2, IL-4, and IL-6 (PIIL246) stimulation;
7. Proteomics data

To predict preterm birth with multi-omic longitudinal data, we evaluated four configurations with machine learning. Configurations refer to i) the unit of observation for model training and evaluation and ii) whether preterm birth was modeled with classification (gestational age at birth < 37 weeks) or regression (the gestational age at birth in weeks). In the sample-level data configuration, each individual blood sample is the unit of observation. This configuration of the data addresses whether sets of immune variables can predict preterm birth at early time points during gestation. In the patient-level data configuration, the unit of observation is the patient and addresses whether the immune trajectory of each patient can predict preterm birth. Classification models were trained to predict preterm birth class labels (gestational age at birth < 37 weeks). 

Model training followed a nested LOOCV structure, with an outer loop used to evaluate the stacked generalization (SG) ensemble and an inner loop used to train the individual submodels. For each omic dataset, least absolute shrinkage and selection operator (Lasso) and Random Forest (RF) submodels were trained within the inner LOOCV loop. To integrate predictive signals across omics layers, we used stacked generalization (SG) in the outer LOOCV loop. In each outer iteration, submodel out-of-sample predictions were used as inputs for a second-layer Lasso model constrained to non-negative coefficients. This nested cross-validation design ensures that SG performance estimates are unbiased and avoids information leakage between submodels and the ensemble layer.


## File Guide


```shell
.
├── README.md
├── configs
│   ├── cfg_1.yaml
│   ├── cfg_2.yaml
│   ├── cfg_3.yaml
│   └── cfg_4.yaml
├── data
│   └── prep
│       └── prince
│           ├── 20240903_totalData.rds
│           └── 20240903_metadata_full.csv
├── pipelines
│   └── snakefile.smk
├── results
│   └── cfg_1
│       ├── submodels         # model predictions (.rds) per iteration
│       └── model_analysis    # summary metrics, feature importance, plots
│   ├── cfg_2
│   ├── cfg_3
│   └── cfg_4
└── src
    └── R
        ├── pipeline_utils.R
        ├── sub_model_train.R        
        ├── sg_model_train.R
        └── model_analysis.R
```


### Configurations

- cfg_1.yaml: patient-level classification analysis
- cfg_2.yaml: patient-level regression analysis
- cfg_3.yaml: sample-level classification analysis
- cfg_4.yaml: sample-level regression analysis

### Code

- `snakefile.smk`: This file runs the scripts below, training the submodels, the ensemble, and generating results.
- `pipeline_utils.R`: This file contains the functions for training all models.
- `sub_model_train.R`: This script trains the submodels with nested LOOCV.
- `sg_model_train.R`: This script trains the SG models with submodel predictions using nested LOOCV.
- `model_analysis.R`: This script analyses model performance and summarizes coefficients and importance scores. 

### Data

- `20240903_totalData.rds`: processed data in long-format from all patients, timepoints, and omics
- `20240903_metadata_full.csv`: sample and patient metadata
- The output of `sub_model_train.R` and `sg_model_train.R` are provided separately to reviewers so that `model_analysis.R` can be evaluated without running the model training scripts.

## Running the pipeline

For each configuration file (`config.yaml`), you can run the snakemake pipeline below:

```shell
cd /root/to/project/
snakemake --snakefile pipelines/snakefile.smk --cores 100 --config cfg=cfg_1
```
- The `--config` argument selects one of the configuration files in `configs/`
- We recommend providing a large `--cores` value so many submodels can be trained in parallel. The R scripts internally cap usage at (available cores - 1) to avoid oversubscribing the machine.
- Model training scripts auto-detect completed submodels and SG models, so re-running or resuming `snakefile.smk` will only train missing models before proceeding to downstream steps.



## Environment

Code was developed and tested with:

- R 4.4.1
- Snakemake 9.14.2
- Key R packages: tidymodels, glmnet, ranger, healthyR.ai, pROC, openxlsx


