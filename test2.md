# Building and evaluation of a PBPK model for raltegravir in adults



| Version     | 0.1                  |
| ----------- | -------------------- |
| OSP Version | 8.0                  |
| Author      | Ibrahim Ince (incei) |

# Table of Contents

- [1 Introduction](#1-introduction)
- [2 Methods](#2-methods)
  - [2.1 Modeling strategy](#2_1 Modeling strategy)
  - [2.2 Data used](#22-data-used)
  - [2.3 Model parameters and assumptions](#2 3-model-parameters-and-assumptions)
- [3 Results and Discussion](#3-results-and-discussion)
  - [3.1 Raltegravir final input parameters](#3.1-raltegravir-final-input-parameters)
  - [3.2 Raltegravir Diagnostics Plots](#32-raltegravir-diagnostics-plots)
  - [3.3: Raltegravir Concentration-Time profiles](#33-raltegravir-concentration-time-profiles)
- [4 Conclusion](#4-conclusion)
- [5 References](#5-references)
- [Introduction](#introduction)
  - [CYP3A4 DDI Network](#cyp3a4-ddi-network)
  - [Objective](#objective)
  - [Qualification Workflow](#qualification-workflow)
- [Qualification of Use Case CYP3A4-mediated DDI](#qualification-of-use-case-cyp3a4-mediated-ddi)
- [Appendix](#appendix)
  - [Evaluation of DDI network (coupled perpetrator / victim simulations)](#evaluation-of-ddi-network-coupled-perpetrator--victim-simulations)
    - [Itraconazole-midazolam DDI](#itraconazole-midazolam-ddi)
    - [Rifampicin-midazolam DDI](#rifampicin-midazolam-ddi)
  - [Evaluation of victim and perpetrator models](#evaluation-of-victim-and-perpetrator-models)
    - [Itraconazole Model](#itraconazole-model)
      - [Itraconazole input parameters](#itraconazole-input-parameters)
      - [Itraconazole model evaluation](#itraconazole-model-evaluation)
        - [Itraconazole concentration-time profiles](#itraconazole-concentration-time-profiles)
        - [Itraconazole overall goodness of fit](#itraconazole-overall-goodness-of-fit)
    - [Midazolam Model](#midazolam-model)
      - [Midazolam input parameters](#midazolam-input-parameters)
      - [Midazolam model evaluation](#midazolam-model-evaluation)
        - [Midazolam concentration-time profiles](#midazolam-concentration-time-profiles)
        - [Midazolam overall goodness of fit](#midazolam-overall-goodness-of-fit)
    - [Rifampicin Model](#rifampicin-model)
      - [Rifampicin input parameters](#rifampicin-input-parameters)
      - [Rifampicin model evaluation](#rifampicin-model-evaluation)
        - [Rifampicin concentration-time profiles](#rifampicin-concentration-time-profiles)
        - [Rifampicin overall goodness of fit](#rifampicin-overall-goodness-of-fit)
  - [Mathematical Implementation of drug-drug interactions](#mathematical-implementation-of-drug-drug-interactions)
  - [Open Systems Pharmacology Suite (OSPS) Introduction](#open-systems-pharmacology-suite-osps-introduction)

# 1 Introduction

The presented model building and evaluation report evaluates the performance of PBPK model for raltegravir in adults.

Raltegravir, sold under the brand name Isentress, is an antiretroviral medication used to treat HIV/AIDS by blocking the establishment of post-integration HIV latency. It is also used as part of post exposure prophylaxis to prevent HIV infection following potential exposure. Raltegravir is only taken orally and is mainly metabolized by UGT1A1 (~70%) [(Kassahun 2007](# 5 References)). The final raltegravir model feautures metabolism by UGT1A1 and to a minor extent by UGT1A9. Additionaly, there is excretion via glomerular filtration. The model adequately describes the pharmacokinetics of raltegravir in adults.

The raltegravir model is a whole-body PBPK model, allowing for dynamic translation between individuals with organs expressing UGT1A1. The raltegravir report demonstrates the level of confidence of the raltegravir PBPK model build with the OSP suite with regard to reliable predictions of Raltegravir PBPK adults during model-informed drug development. The presented raltegravir PBPK model as well as the respective evaluation plan and PBPK report are provided open-source and transparently documented (https://github.com/Incei/Raltegravir-Model).

# 2 Methods

## 2_1 Modeling strategy

The general concept of building a PBPK model has previously been described by Kuepfer et al. ([Kuepfer 2016](# 5 References)) Regarding the relevant anthropometric (height, weight) and physiological parameters (e.g. blood flows, organ volumes, binding protein concentrations, hematocrit, cardiac output) in adults was gathered from the literature and has been previously published ([PK-Sim Ontogeny Database Version 7.3](# 5 References)). The information was incorporated into PK-Sim® and was used as default values for the simulations in adults.

The  applied activity and variability of plasma proteins and active processes that are integrated into PK-Sim® are described in the publicly available ‘PK-Sim® Ontogeny Database Version 7.3 ([Schlender 2016](# 5 References)) or otherwise referenced for the specific process.

First, a base mean model was built using data from the single dose escalation study study to find an appropriate structure describing the PK of Raltegravir plasma. The mean PK model was developed using a typical European individual. Unknown parameters were identified using the Parameter Identification module provided in PK-Sim®. Structural model selection was mainly guided by visual inspection of the resulting description of data and biological plausibility.

Once the appropriate structural model was identified, additional parameters for different formulations were identified. 

A final PBPK model was established and simulations were compared to the reported data to evaluate model appropriateness and to assess model qualification, by means of diagnostics plots and predicted versus observed concentration-time profiles, of which the results support an adequate prediction of the PK in adults.

During model building, uncertainties in data-quality, as well as study differences may cause not being able to adequately describe the PK of all reported clinical study data. 

## 22 Data used

### 2.2.1	In vitro / physico-chemical data

A literature search was performed to collect available information on physical chemical properties of raltegravir. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit**    | **Raltegravir literature**                                  | **Description**                                 |
| :-------------- | ----------- | ----------------------------------------------------------- | ----------------------------------------------- |
| MW              | g/mol       | 586.2 ([drugbank.ca](# 5 References))                       | Molecular weight                                |
| pKa             |             | 7.67 ([Moss 2012](# 5 References))                          | Acid dissociation constant                      |
| Solubility (pH) | mg/L        | Reference pH-dependent table  ([Moss 2013](# 5 References)) | Solubility                                      |
| logP            |             | 0.58 ([Moss 2012](# 5 References))                          | Partition coefficient between octanol and water |
| fu              |             | 0.17 ([Laufer 2009](# 5 References))                        | Fraction unbound                                |
| Km UGT1A1       | µM          | 99 ([Kassahun 2007](# 5 References))                        | Michaelis-Menten constant                       |
| Vmax UGT1A1     | nmol/min/mg | 0.89 ([Kassahun 2007](# 5 References))                      | Maximum rate of reaction                        |
| Km UGT1A9       | µM          | 296 ([Kassahun 2007](# 5 References))                       | Michaelis-Menten constant                       |
| Vmax UGT1A9     | nmol/min/mg | 0.53 ([Kassahun 2007](# 5 References))                      | Maximum rate of reaction                        |

### 2.2.2	Clinical data

A literature search was performed to collect available clinical data on Raltegravir in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                      | Study description                                            |
| :------------------------------- | :----------------------------------------------------------- |
| [Iwamoto 2008](# 5 References)   | Single- and multiple-dose escalation study in healthy subjects |
| [Iwamoto 2009](# 5 References)   | Effects of ritonavir and efavirenz on safety, tolerability and pharmacokinetics of raltegravir |
| [Markowitz 2006](# 5 References) | Monotherapy, followed by a longer term combination therapy of raltegravir versus efavirenz |
| [Kassahun 2007](# 5 References)  | Pharmacokinetics study in healthy adults                     |
| [Rhee 2014](# 5 References)      | Pediatric formulation study in healthy adults                |
| [Wenning 2009](# 5 References)   | Effect of rifampin on the pharmacokinetics of raltegravir    |

## 2 3 Model parameters and assumptions

### 2.3.1	Absorption

As no intravenous data is currently available to study systemic clearance of raltegravir in vivo, only oral data was used for model building. For oral administration the following parameters play a role with regards to the absorption kinetics of a compound, which can be estimated with PBPK: solubility, lipophilicity and intestinal permeability. Moss et al. ([Moss 2013](#5 References)) published values for raltegravir solubility in population groups with very low pH, low pH, medium pH, high pH and very high pH, after a single 400 mg dose of raltegravir. For the raltegravir PBPK model we have applied the medium pH group for creating a pH dependent solubility profile throughout the intestinal tract. The lipophilicity as well as pKa of raltegravir was also published by Moss et al ([Moss 2012, Moss 2013](#5 References)) to be 0.58 (as log partition coefficient between octanol and water (pH 7) and 6.67 (acid), respectively. These values were applied and fixed into the raltegravir PBPK model, without further optimization. Regarding intestinal transcellular permeability (Pint), Moss et al (Moss 2012) reported a range of apical to basolateral apparent permeability in Caco-2 monolayer at different pH values. Using published functions Pint can be calculated from Caco-2 cell membrane permeability measurements (Parrot et al. ([Parrot 2002](#5 References)), Thelen et al. ([Thelen 2010](#5 References), Sun et al. ([Sun 2002](#5 References)) and Sjögren et al. ([Sjögren 2013](#5 References)). However as no reference/calibrator compound was available to correct for inter-study variability, these functions could not be applied, and it was decided to estimate the Pint from in vivo clinical data instead. Nevertheless, for plausibility check, a theoretical Pint was calculated using the aforementioned functions without correction, resulting in a range of Pint from 2.14E-04 to 1.47E-09 cm/min. The finally estimated (based on in vivo data) Pint falls within this range.

**Table 2.** Reported Caco-permeability and calculated theoretical effective permeability (intestinal transcellular permeability) values for raltegravir via different reported functions, lacking a reference compound for correcting inter-study variability.

| Reference publication of reported function | **Ph apical-basolateral** | **Peff (cm/min)** | **Reference compound available for correcting Inter study variability** |
| ------------------------------------------ | ------------------------- | ----------------- | ------------------------------------------------------------ |
| Raltegravir Caco permeability (Moss 2012 ) | 7.4                       | 6.60E-6           | -                                                            |
| Raltegravir Caco permeability (Moss 2012)  | 6.5                       | 9.20E-6           | -                                                            |
| Parrot 2002                                | 7.4                       | 2.14E-04          | Not available                                                |
| Thelen 2010                                | 7.4                       | 1.47E-09          | Not available                                                |
| Sjögren 2013                               | 7.4                       | 1.03606E-6        | Not available                                                |
| Sun 2002                                   | 7.4                       | 2.77E-4           | Not available                                                |
| Sun 2002                                   | 6.5                       | 2.86E-4           | Not available                                                |
| Simcyp (*)                                 | 6.5                       | 4.62E-4           | Not available                                                |

*Not published as paper, Simcyp applied an adapted version of Sun et al 2002 

### 2.3.2	Distribution

Laufer et al. ([Laufer 2009](#5 References)) published a fu in humans to be 0.17. Barau et al 2013 reported that raltegravir binds to serum albumin, and not alpha glycoprotein, and are built-in as such into the PBPK model.

After testing the available organ-plasma partition coefficient and cell permeability calculation methods built in PK-Sim, observed clinical data was best described by choosing the partition coefficient calculation by Rodgers and Rowland, and cell permeability calculation by PK-Sim standard. Specific organ permeability normalized to surface area was automatically calculated by PK-Sim.

### 2.3.3	Metabolism and Elimination

Kassahun et al. ([Kassahun 2007](#5 References)) studied the absorption, metabolism, and excretion of raltegravir in healthy volunteers after a single oral dose of 200 mg (200Ci) of [14C] raltegravir. Human liver microsomal incubations confirmed the dominat role of UGT metabolism for raltegravir. Additionally, data from incubations using cDNA-expressed UGTs indicate that the major mechanism of clearance of raltegravir in humans is UGT1A1-mediated glucuronidation. Raltegravir was in particular converted to M2 by UGT1A1 and 1A9. The apparent Km values for the glucuronidation of raltegravir by UGT1A1 and UGT1A9 were 99 (standard deviation (SD): 16) and 296 (SD: 55) µM, respectively. The corresponding Vmax values (nmol/min/mg) were 0.89 (SD: 0.05) for UGT1A1, and 0.53 (SD: 0.06) for UGT1A9.

Based on this information, the reported in vitro Km values for UGT1A1 and 1A9 in the model. Reported Vmax values were of the dimension nmol/min/mg protein and thus not directly transferable into the PBPK model. Therefore, a unique scaling factor fUGT on the in vitro Vmax values was estimated to match observed in vivo data, and keeping the relative relationship between those in vitro values (0.89 and 0.53 nmol/min/mg) for UGT1A1 and UGT1A9 fixed according to:

Vmax,UGT1A1 = fUGT * Vmax,in-vitro,UGT1A1 

Vmax,UGT1A9 = fUGT * Vmax,in-vitro,UGT1A9

It is especially important to fix the relative contribution of both enzymes as a ratio to ensure that, when scaling to other populations (e.g. children where both UGT’s undergo a different ontogeny pattern, or patients who have differently reduced amounts of UGT1A1 vs 1A9) the relative contributions can be adequately scaled. 
Note that the estimated scaling factor fUGT will be directly implemented into the final in vivo Vmax values (only Vmax,UGT1A1 and Vmax,UGT1A9 will be reported in section [3](## 3 Results and Discussion))

Finally, as ~9% of the dose is excreted in human urine as unchanged parent compound, GFR is introduced in the raltegravir PBPK model.

# 3 Results and Discussion

The PBPK model **raltegravir** was developed with clinical pharmacokinetic data covering 4 different oral formulation and a dose range of 10-1600mg, including single dose (SD) as well as multiple dose (MD) clinical data. 

As there were 4 different oral formulations available for model evaluation, all formulations require an estimation of the dissolution kinetics via a Weibull function. This function requires the estimation of 2 parameters, the dissolution time (time where 50% of the drug is dissolved), and dissolution shape (shape parameter of the Weibull function). Therefore, to minimize the amount of parameters for fitting, as a first step, the PK study data (lactose formulation) by Iwamoto et al. ([Iwamoto 2007](# 5 References)) was fitted which includes SD escalation and hast a broad dose-range (10mg-1600mg) to capture (non-) linearity. During the model-fitting, the following parameters were estimated (all other parameters were fixed to reported values):

- Vmax (as unique scaling factor fUGT, as described in section [2.3.3](## 2.3 Model parameters and assumptions)) 
- Weibul function paramters: Dissolution time and dissolution shape
- Specific intestinal permeability (transcellular)

The fit resulted in an adequate description of all data. As there is no iv data available, it was not possible to clearly distinguish between clearance and absorption, resulting in a considerable correlation between Vmax and dissolution shape (Weibull). An attempt to fix Vmax to reported in vitro values, and only estimating absorption (lipophilicity and intestinal transcellular permeability) resulted in an underprediction of the clearance, and clearly indicated a need for increase in clearance. As described above, no reported intestinal permeability was found other than caco-permeability. Caco-permeability could not be translated to effective intestinal permeability without a reference compound. Therefore it was decided to continue with the model where both Pint and Vmax were estimated.

As a second step, clincial study data for all other formulations summarised in section [2.2.2](## 2.2 Data used) were included for model fitting, including film-coated tablets (100-400mg MD, 200-400mg SD), chewable tablets (400mg fasted + fed) and oral granules in suspension (400mg). In this step, only the Weibull functions were estimated with all other parameters fixed based on the first step. Finally, as the Weibull functions were highly correlated (as expected), only dissolution shape was estimated as a last step. The model results show that the PBPK model of raltegravir adequately described the date for all formulations and doses available.

## 3.1 Raltegravir final input parameters

The compound parameter values of the final raltegravir PBPK model are illustrated below.



# Formulation: chewable tablet

Type: Weibull

## Parameters

| Name                             | Value                | Value Origin             |      |
| -------------------------------- | -------------------- | ------------------------ | ---- |
| Dissolution time (50% dissolved) | 1.0000049774E-05 min | Parameter Identification |      |
| Lag time                         | 0 min                |                          |      |
| Dissolution shape                | 0.050078869          | Parameter Identification |      |
| Use as suspension                | Yes                  |                          |      |

# Formulation: filmcoated tablet (original Merck formulation)

Type: Weibull

## Parameters

| Name                             | Value      | Value Origin             |      |
| -------------------------------- | ---------- | ------------------------ | ---- |
| Dissolution time (50% dissolved) | 500 min    | Parameter Identification |      |
| Lag time                         | 0 min      |                          |      |
| Dissolution shape                | 0.03536656 | Parameter Identification |      |
| Use as suspension                | Yes        |                          |      |

# Compound: Raltegravir

## Parameters

| Name                                             | Value                 | Value Origin                           | Alternative | Default |      |
| ------------------------------------------------ | --------------------- | -------------------------------------- | ----------- | ------- | ---- |
| Solubility table                                 | 40 mg/l               | Publication-In Vitro-Moss 2013 Table 2 | Moss 2013   | True    |      |
| Lipophilicity                                    | 0.58 Log Units        | Parameter Identification               | Moss 2012   | True    |      |
| Fraction unbound (plasma, reference value)       | 0.17                  | Publication-In Vitro-Laufer 2009       | Measurement | True    |      |
| Specific intestinal permeability (transcellular) | 2.8481843854E-07 cm/s | Parameter Identification               | Fit         | True    |      |
| F                                                | 1                     | Publication-Other-Drugbank.ca          |             |         |      |
| Is small molecule                                | Yes                   |                                        |             |         |      |
| Molecular weight                                 | 444.4163 g/mol        | Publication-Other-Drugbank.ca          |             |         |      |
| Plasma protein binding partner                   | Albumin               |                                        |             |         |      |

## Calculation methods

| Name                    | Value               |      |
| ----------------------- | ------------------- | ---- |
| Partition coefficients  | Rodgers and Rowland |      |
| Cellular permeabilities | PK-Sim Standard     |      |

## Processes

### Systemic Process: Glomerular Filtration-Kassahun 2007

Species: Human

#### Parameters

| Name         | Value | Value Origin                       |      |
| ------------ | ----: | ---------------------------------- | ---- |
| GFR fraction |     1 | Publication-In Vitro-Kassahun 2007 |      |

### Metabolizing Enzyme: UGT1A1-Kassahun 2007

Molecule: UGT1A1

#### Parameters

| Name                               | Value                                 | Value Origin                       |      |
| ---------------------------------- | ------------------------------------- | ---------------------------------- | ---- |
| In vitro Vmax for liver microsomes | 2.7351231632 nmol/min/mg mic. protein | Parameter Identification           |      |
| Km                                 | 99 µM                                 | Publication-In Vitro-Kassahun 2007 |      |

### Metabolizing Enzyme: UGT1A9-Kassahun 2007

Molecule: UGT1A9

#### Parameters

| Name                               | Value                                 | Value Origin                       |      |
| ---------------------------------- | ------------------------------------- | ---------------------------------- | ---- |
| In vitro Vmax for liver microsomes | 1.6287812095 nmol/min/mg mic. protein | Parameter Identification           |      |
| Km                                 | 296 µM                                | Publication-In Vitro-Kassahun 2007 |      |

# Formulation: Weibull (granules)

Type: Weibull

## Parameters

| Name                             | Value                | Value Origin             |      |
| -------------------------------- | -------------------- | ------------------------ | ---- |
| Dissolution time (50% dissolved) | 0.00010000047426 min | Parameter Identification |      |
| Lag time                         | 0 min                |                          |      |
| Dissolution shape                | 0.0654456264         | Parameter Identification |      |
| Use as suspension                | Yes                  |                          |      |

# Formulation: Weibull (lactose formulation)

Type: Weibull

## Parameters

| Name                             | Value              | Value Origin             |      |
| -------------------------------- | ------------------ | ------------------------ | ---- |
| Dissolution time (50% dissolved) | 2.30152527E-10 min | Parameter Identification |      |
| Lag time                         | 0 min              |                          |      |
| Dissolution shape                | 0.0389537131       | Parameter Identification |      |
| Use as suspension                | Yes                |                          |      |

## 3.2 Raltegravir Diagnostics Plots

Below you find the goodness-of-fit visual diagnostic plots for raltegravir PBPK model performance (observed versus individually simulated plasma concentration and weighted residuals versus time) of all data used for model building.

![001_plotGOFMergedPredictedVsObserved.png](images\003_3_Results_and_Discussion\002_3_2_Raltegravir_Diagnostics_Plots\001_plotGOFMergedPredictedVsObserved.png)

![002_plotGOFMergedResidualsOverTime.png](images\003_3_Results_and_Discussion\002_3_2_Raltegravir_Diagnostics_Plots\002_plotGOFMergedResidualsOverTime.png)

GMFE = 1.487148 

## 3.3: Raltegravir Concentration-Time profiles

Simulated versus observed plasma concentration-time profiles of all data are listed below.

![001_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\001_plotTimeProfile.png)

![002_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\002_plotTimeProfile.png)

![003_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\003_plotTimeProfile.png)

![004_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\004_plotTimeProfile.png)

![005_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\005_plotTimeProfile.png)

![006_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\006_plotTimeProfile.png)

![007_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\007_plotTimeProfile.png)

![008_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\008_plotTimeProfile.png)

![009_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\009_plotTimeProfile.png)

![010_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\010_plotTimeProfile.png)

![011_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\011_plotTimeProfile.png)

![012_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\012_plotTimeProfile.png)

![013_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\013_plotTimeProfile.png)

![014_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\014_plotTimeProfile.png)

![015_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\015_plotTimeProfile.png)

![016_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\016_plotTimeProfile.png)

![017_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\017_plotTimeProfile.png)

![018_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\018_plotTimeProfile.png)

![019_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\019_plotTimeProfile.png)

![020_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\020_plotTimeProfile.png)

![021_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\021_plotTimeProfile.png)

![022_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\022_plotTimeProfile.png)

![023_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\023_plotTimeProfile.png)

![024_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\024_plotTimeProfile.png)

![025_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\025_plotTimeProfile.png)

![026_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\026_plotTimeProfile.png)

![027_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\027_plotTimeProfile.png)

![028_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\028_plotTimeProfile.png)

![029_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\029_plotTimeProfile.png)

![030_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\030_plotTimeProfile.png)

![031_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\031_plotTimeProfile.png)

![032_plotTimeProfile.png](images\003_3_Results_and_Discussion\003_3_3__Raltegravir_Concentration-Time_profiles\032_plotTimeProfile.png)

# 4 Conclusion

The final raltegravir PBPK model applies metabolism by UGT1A1,  UGT1A9 and glomerular filtration and adequately describes the pharmacokinetics of raltegravir in adults receiving SD, MD of Raltegravir ranging from 10mg to 1600mg, including four different oral formulations. 

This model could be applied for the investigation of DDI, and translation to special populations such as pediatrics with regard to UGT1A1 and 1A9 metabolism.

# 5 References

**Iwamoto 2008** Iwamoto M, Wenning LA, Petry AS, Laethem M, De Smet M, Kost JT, Merschman SA, Strohmaier KM, Ramael S, Lasseter KC, Stone JA, Gottesdiener KM, Wagner JA. Safety, tolerability, and pharmacokinetics of raltegravir after single and multiple doses in healthy subjects. Clin Pharmacol Ther. 2008 Feb;83(2):293-9. Epub 2007 Aug 22.

**Iwamoto 2009** Iwamoto M, Wenning LA, Nguyen BY, Teppler H, Moreau AR, Rhodes RR, Hanley WD, Jin B, Harvey CM, Breidinger SA, Azrolan N, Farmer HF Jr, Isaacs RD, Chodakewitz JA, Stone JA, Wagner JA. Effects of omeprazole on plasma levels of raltegravir. Clin Infect Dis. 2009 Feb 15;48(4):489-92. doi: 10.1086/596503.

**Kassahun 2007** Kassahun K, McIntosh I, Cui D, Hreniuk D, Merschman S, Lasseter K, Azrolan N, Iwamoto M, Wagner JA, Wenning LA. Metabolism and disposition in humans of raltegravir (MK-0518), an anti-AIDS drug targeting the human immunodeficiency virus 1 integrase enzyme. Drug Metab Dispos. 2007 Sep;35(9):1657-63. Epub 2007 Jun 25.

**Kuepfer 2016** Kuepfer L, Niederalt C, Wendl T, Schlender JF, Willmann S, Lippert J, Block M, Eissing T, Teutonico D. Applied Concepts in PBPK Modeling: How to Build a PBPK/PD Model.CPT Pharmacometrics Syst Pharmacol. 2016 Oct;5(10):516-531. doi: 10.1002/psp4.12134. Epub 2016 Oct 19. 

**Wenning 2009** Larissa A. Wenning,, William D. Hanley, Diana M. Brainard, Amelia S. Petry, Kalyan Ghosh, Bo Jin, Eric Mangin, Thomas C. Marbury, Jolene K. Berg, Jeffrey A. Chodakewitz, Julie A. Stone,1 Keith M. Gottesdiener, John A. Wagner, and Marian Iwamoto. Effect of Rifampin, a Potent Inducer of Drug-Metabolizing Enzymes, on the Pharmacokinetics of Raltegravir. Antimicrob Agents Chemother. 2009 Jul; 53(7): 2852–2856.

**Laufer 2009** Laufer R, Paz OG, Di Marco A, Bonelli F, Monteagudo E, Summa V, Rowley M. Quantitative prediction of human clearance guiding the development of Raltegravir (MK-0518, isentress) and related HIV integrase inhibitors. Drug Metab Dispos. 2009 Apr;37(4):873-83. doi: 10.1124/dmd.108.023804. Epub 2009 Jan 14.

**Markowitz 2006** Markowitz M, Morales-Ramirez JO, Nguyen BY, Kovacs CM, Steigbigel RT, Cooper DA, Liporace R, Schwartz R, Isaacs R, Gilde LR, Wenning L, Zhao J, Teppler H. Antiretroviral activity, pharmacokinetics, and tolerability of MK-0518, a novel inhibitor of HIV-1 integrase, dosed as monotherapy for 10 days in treatment-naive HIV-1-infected individuals. J Acquir Immune Defic Syndr. 2006 Dec 15;43(5):509-15.

**Moss 2013** Moss DM, Siccardi M, Back DJ, Owen A. Predicting intestinal absorption of raltegravir using a population-based ADME simulation. J Antimicrob Chemother. 2013 Jul;68(7):1627-34. doi: 10.1093/jac/dkt084. Epub 2013 Mar 20.

**Moss 2012** Moss DM, Siccardi M, Murphy M, Piperakis MM, Khoo SH, Back DJ, Owen A. Divalent metals and pH alter raltegravir disposition in vitro. Antimicrob Agents Chemother. 2012 Jun;56(6):3020-6. doi: 10.1128/AAC.06407-11. Epub 2012 Mar 26.

**Parrott 2008** Parrott N, Lave T. Applications of physiologically based absorption models in drug discovery and development. Mol Pharm. 2008 Sep-Oct;5(5):760-75. doi: 10.1021/mp8000155. Epub 2008 Jun 12.

**PK-Sim Ontogeny Database Version 7.3** (https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/blob/38cf71b384cfc25cfa0ce4d2f3addfd32757e13b/PK-Sim%20Ontogeny%20Database%20Version%207.3.pdf)

**Rhee 2014** Rhee EG, Rizk ML, Brainard DM, Gendrano IN 3rd, Jin B, Wenning LA, Wagner JA, Iwamoto M. A pharmacokinetic comparison of adult and paediatric formulations of raltegravir in healthy adults. Antivir Ther. 2014;19(6):619-24. doi: 10.3851/IMP2765. Epub 2014 Mar 7.

**Schlender 2016** Schlender JF, Meyer M, Thelen K, Krauss M, Willmann S, Eissing T, Jaehde U. Development of a Whole-Body Physiologically Based Pharmacokinetic Approach to Assess the Pharmacokinetics of Drugs in Elderly Individuals. Clin Pharmacokinet. 2016 Dec;55(12):1573-1589. 

**Sjögren 2013** Sjögren E, Westergren J, Grant I, Hanisch G, Lindfors L, Lennernäs H, Abrahamsson B, Tannergren C. In silico predictions of gastrointestinal drug absorption in pharmaceutical product development: application of the mechanistic absorption model GI-Sim. Eur J Pharm Sci. 2013 Jul 16;49(4):679-98. doi: 10.1016/j.ejps.2013.05.019. Epub 2013 May 29.

**Sun 2002** Sun D, Lennernas H, Welage LS, Barnett JL, Landowski CP, Foster D, Fleisher D, Lee KD, Amidon GL. Comparison of human duodenum and Caco-2 gene expression profiles for 12,000 gene sequences tags and correlation with permeability of 26 drugs. Pharm Res. 2002 Oct;19(10):1400-16. 

**Thelen 2011** Thelen K, Coboeken K, Willmann S, Burghaus R, Dressman JB, Lippert J. Evolution of a detailed physiological model to simulate the gastrointestinal transit and absorption process in humans, part 1: oral solutions. J Pharm Sci. 2011 Dec;100(12):5324-45. doi: 10.1002/jps.22726. Epub 2011 Oct 12

# Introduction

## CYP3A4 DDI Network

## CYP3A4 DDI Network

To qualify the OSP suite for the prediction of the CYP3A4 DDI potential of new drugs, a set of verified PBPK models of index perpetrators, covering the range from strong CYP3A4 induction to strong inhibition, and respective CYP3A4 DDI victim drugs is specified to set up a CYP3A4-mediated DDI modeling network. 



The following perpetrator compounds were selected: 

- rifampicin (strong CYP3A4 inducer), 
- efavirenz (moderate CYP3A4 inducer), 
- fluvoxamine (weak/moderate CYP3A4 inhibitor), 
- erythromycin (moderate CYP3A4 inhibitor), 
- clarithromycin (strong CYP 3A4 inhibitor), 
- itraconazole including metabolites (strong CYP3A4 inhibitior). 



The following sensitive CYP3A4 substrates as victim drugs were selected:

- midazolam, 
- triazolam, 
- alprazolam,
- alfentanil. 



A detailed description of the performed workflow in terms of DDI network modeling has been previously published [[1](#reference)].

Figure 1 shows the prespecified and developed DDI modeling network of interacting perpetrator and victim drugs for the OSP suite qualification of predicting CYP3A4-mediated DDI.



**Figure** **1: CYP3A4 DDI modeling network**
![DDI CYP3A4 network](images/DDI_CYP3A4_Compound_Network.png)



### Reference

[1] [Hanke N, Frechen S, Moj D, Britz H, Eissing T, Wendl T, Lehr T. PBPK Models for CYP3A4 and P-gp DDI Prediction: A Modeling Network of Rifampicin, Itraconazole, Clarithromycin, Midazolam, Alfentanil, and Digoxin. CPT Pharmacometrics Syst Pharmacol. (2018)](https://doi.org/10.1002/psp4.12343)

## Objective

## Objective 

This qualification report evaluates the predictive performance of the open systems pharmacology (OSP) suite to predict cytochrome P450 3A4 (CYP3A4)-mediated drug-drug interactions. A PBPK model network of selected index CYP3A4 DDI perpetrators, covering the range from strong induction to strong inhibition, and respective sensitive index CYP3A4 victim drugs has been built and evaluated. All models are whole-body PBPK models, allowing for dynamic DDI assessment in organs expressing CYP3A4. 

The qualification report demonstrates the level of confidence of the OSP suite with regard to reliable PBPK predictions of CYP3A4-mediated DDIs by means of pre-specified qualification measures. The presented PBPK models as well as the respective qualification plan and qualification report are transparently documented and provided open-source [[1,2](#reference)]. 



### References

[1] [www.open-systems-pharmacology.org](http://www.open-systems-pharmacology.org)

[2] [INSERT LINK TO REPOSITORY](INSERT LINK TO REPOSITORY)





## Qualification Workflow

## Automatic (re)-qualification workflow 

OSP provides a dynamic landscape of model repositories and a database of observed clinical data [[1](#reference)]. Additionally, a technical framework to assess confidence of a specific intended use has been developed (qualification runner and reporting engine). This framework allows for an automatic (re)-qualification workflow of the OSP suite, comprising the following steps (Figure 1):

- PBPK model development and verification with observed data,
- Qualification plan generation,
- Qualification plan execution,
- Qualification report generation.

**Figure 1: OSP suite automatic (re)-qualification workflow**
![OSP qualification workflow](images/OSP_Qualification_Workflow_1.png)

In a first step, the respective qualification scenario is saved in a special qualification repository on GitHub [[1](#reference)]. This qualification scenario repository contains a detailed qualification plan that links and combines respective models and data to address the use case that shall be qualified. Therefore, the qualification plan consists of: 

- PK-Sim project files,
- Additional model building steps (if applicable),
- Description of potential cross-dependencies between PK-Sim project files (if applicable),
- Observed data (needed for model development and verification),
- Qualification scenario description text modules
- Detailed report settings to describe the generation of charts and qualification measures. 

PK-Sim projects, observed data sets, and qualification scenario text modules are deposited in distinct repositories and are referenced by the qualification plan (Figure 2).

**Figure 2: Qualification scenario repository landscape on GitHub**
![OSP qualification workflow detail](images/OSP_Qualification_Workflow_2.png)

In a second step the qualification runner [[2](#reference)] processes the qualification plan, i.e. all project parts are exported and prepared for the reporting engine [[3](#reference)]. The reporting engine provides a validated environment (currently implemented in MATLAB®, a transfer to R is in development) for model execution and finally generates the qualification report. This report contains the evaluation of the individual PBPK models with observed data (i.e. standard goodness of fit plots, visual predictive checks) and a comprehensive qualification of the specific use case assessing the predictive performance of the OSP suite by means of a predefined set of qualification measures and charts. 

The automated execution of the described workflow can be triggered to assess re-qualification in case new data, changes in model structure or parameterization, or new OSP suite releases arise.

### References

[1] [www.open-systems-pharmacology.org](http://www.open-systems-pharmacology.org)

[2] [https://github.com/Open-Systems-Pharmacology/QualificationRunner](https://github.com/Open-Systems-Pharmacology/QualificationRunner)

[3] [https://github.com/Open-Systems-Pharmacology/Reporting-Engine](https://github.com/Open-Systems-Pharmacology/Reporting-Engine)



# Qualification of Use Case CYP3A4-mediated DDI

|                Perpetrator |        Victim | Predicted AUC Ratio | Observed AUC Ratio | Pred/Obs AUC Ratio | Predicted CMAX Ratio | Observed CMAX Ratio | Pred/Obs CMAX Ratio |       Reference |
| -------------------------: | ------------: | ------------------: | -----------------: | -----------------: | -------------------: | ------------------: | ------------------: | --------------: |
|   Itraconazole, 200 mg, PO | Midazolam, IV |              2.3152 |             3.2258 |             0.7177 |               1.0117 |                   - |                   - |    Olkkola 1996 |
|   Itraconazole, 100 mg, PO | Midazolam, PO |              3.4697 |             5.7451 |            0.60394 |               1.9613 |              2.5588 |              0.7665 |     Ahonen 1995 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              4.9294 |               7.97 |            0.61849 |               2.2277 |                3.12 |             0.71402 |    Backman 1998 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              1.1177 |               2.63 |              0.425 |                1.074 |                1.92 |             0.55939 |    Backman 1998 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              5.2277 |               10.8 |            0.48404 |               2.2378 |                 3.4 |             0.65818 |    Olkkola 1994 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              4.5464 |                3.4 |             1.3372 |               2.1997 |                 1.8 |               1.222 |    Olkkola 1996 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              6.1785 |                6.6 |            0.93613 |               2.5665 |                 2.5 |              1.0266 |    Olkkola 1996 |
|    Itraconazole, 50 mg, PO | Midazolam, PO |              3.2322 |                  2 |             1.6161 |               2.2008 |                   - |                   - |  Templeton 2010 |
|   Itraconazole, 200 mg, PO | Midazolam, PO |              7.3895 |                4.7 |             1.5722 |               3.4395 |                   - |                   - |  Templeton 2010 |
|   Itraconazole, 400 mg, PO | Midazolam, PO |              9.2123 |                5.4 |              1.706 |               3.7303 |                   - |                   - |  Templeton 2010 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |              0.4937 |            0.45833 |             1.0772 |              0.75214 |                   - |                   - |     Gorski 2003 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |             0.47784 |            0.37931 |             1.2598 |              0.90627 |                   - |                   - |   Kharasch 1997 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |             0.48247 |            0.52113 |            0.92582 |              0.90861 |                1.01 |             0.89965 |   Kharasch 2004 |
|       Rifampicin, 5 mg, PO | Midazolam, IV |             0.79138 |               0.84 |            0.94212 |              0.98019 |              1.0323 |             0.94956 |   Kharasch 2011 |
|      Rifampicin, 10 mg, PO | Midazolam, IV |             0.70846 |               0.77 |            0.92007 |               0.9682 |              1.0645 |             0.90952 |   Kharasch 2011 |
|      Rifampicin, 25 mg, PO | Midazolam, IV |             0.61322 |               0.63 |            0.97337 |               0.9493 |             0.83871 |              1.1319 |   Kharasch 2011 |
|      Rifampicin, 75 mg, PO | Midazolam, IV |             0.53865 |                0.6 |            0.89776 |              0.92855 |              1.3226 |             0.70208 |   Kharasch 2011 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |             0.47327 |            0.65501 |            0.72254 |              0.90373 |               1.106 |             0.81715 |       Link 2008 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |             0.47882 |               0.51 |            0.93887 |              0.90704 |                   - |                   - | Phimmasone 2001 |
|     Rifampicin, 600 mg, PO | Midazolam, IV |              0.4937 |            0.57947 |            0.85199 |              0.75227 |                   - |                   - |     Szalat 2007 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |            0.038283 |              0.041 |            0.93374 |              0.09826 |            0.061818 |              1.5895 |    Backman 1996 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |            0.038597 |              0.023 |             1.6781 |             0.098261 |               0.054 |              1.8196 |    Backman 1998 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |             0.11838 |              0.132 |            0.89681 |               0.2238 |               0.202 |              1.1079 |    Backman 1998 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |             0.02853 |            0.12449 |            0.22917 |             0.070257 |             0.16957 |             0.41434 |      Chung 2006 |
|     Rifampicin, 450 mg, PO | Midazolam, PO |            0.038473 |            0.44118 |           0.087206 |               0.0827 |             0.22581 |             0.36624 |        Eap 2004 |
|     Rifampicin, 450 mg, PO | Midazolam, PO |            0.033981 |           0.052239 |            0.65049 |                0.085 |             0.11154 |             0.76207 |        Eap 2004 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |            0.047973 |           0.045349 |             1.0579 |              0.11693 |            0.067039 |              1.7442 |     Gorski 2003 |
|     Rifampicin, 300 mg, PO | Midazolam, PO |            0.043642 |           0.057161 |             0.7635 |              0.10789 |             0.12092 |             0.89225 |     Gurley 2006 |
|     Rifampicin, 300 mg, PO | Midazolam, PO |            0.043642 |           0.060317 |            0.72354 |              0.10789 |             0.10762 |              1.0024 |    Gurley 2008a |
|     Rifampicin, 600 mg, PO | Midazolam, PO |            0.035413 |           0.052632 |            0.67284 |             0.075883 |             0.10989 |             0.69054 |   Kharasch 2004 |
|       Rifampicin, 5 mg, PO | Midazolam, PO |             0.42184 |                0.8 |             0.5273 |              0.55551 |                 0.8 |             0.69438 |   Kharasch 2011 |
|      Rifampicin, 10 mg, PO | Midazolam, PO |             0.29076 |               0.68 |            0.42759 |              0.42549 |             0.93333 |             0.45588 |   Kharasch 2011 |
|      Rifampicin, 25 mg, PO | Midazolam, PO |             0.16247 |                0.4 |            0.40616 |              0.27343 |             0.50667 |             0.53966 |   Kharasch 2011 |
|      Rifampicin, 75 mg, PO | Midazolam, PO |            0.080359 |               0.25 |            0.32143 |              0.15436 |                0.34 |               0.454 |   Kharasch 2011 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |            0.030501 |           0.015549 |             1.9616 |             0.077383 |            0.034865 |              2.2195 |       Link 2008 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |             0.20754 |              0.123 |             1.6873 |              0.27316 |               0.162 |              1.6862 |    Reitman 2011 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |             0.34588 |              0.383 |            0.90307 |              0.49107 |               0.403 |              1.2185 |    Reitman 2011 |
|     Rifampicin, 600 mg, PO | Midazolam, PO |             0.93395 |              0.815 |             1.1459 |              0.93066 |               0.731 |              1.2731 |    Reitman 2011 |
|                        AUC |        Number |           Ratio [%] |                    |                    |                      |                     |                     |                 |
| -------------------------: |        -----: |           --------: |                    |                    |                      |                     |                     |                 |
|               Points total |            38 |                   - |                    |                    |                      |                     |                     |                 |
| Points within Guest et al. |            28 |             73.6842 |                    |                    |                      |                     |                     |                 |
|       Points within 2-fold |            31 |             81.5789 |                    |                    |                      |                     |                     |                 |

![003_plotDDIRatioAUCpredictedVsObserved.png](images\01_Qualification_of_Use_Case_CYP3A4-mediated_DDI\003_plotDDIRatioAUCpredictedVsObserved.png)

![004_plotDDIRatioAUCresidualsVsObserved.png](images\01_Qualification_of_Use_Case_CYP3A4-mediated_DDI\004_plotDDIRatioAUCresidualsVsObserved.png)

|                       CMAX | Number | Ratio [%] |
| -------------------------: | -----: | --------: |
|               Points total |     30 |         - |
| Points within Guest et al. |     16 |   53.3333 |
|       Points within 2-fold |     25 |   83.3333 |

![006_plotDDIRatioCMAXpredictedVsObserved.png](images\01_Qualification_of_Use_Case_CYP3A4-mediated_DDI\006_plotDDIRatioCMAXpredictedVsObserved.png)

![007_plotDDIRatioCMAXresidualsVsObserved.png](images\01_Qualification_of_Use_Case_CYP3A4-mediated_DDI\007_plotDDIRatioCMAXresidualsVsObserved.png)

# Appendix

## Evaluation of DDI network (coupled perpetrator / victim simulations)

## Evaluation of DDI network (coupled perpetrator / victim simulations)

The published DDI studies between the respective perpetrators and victim drugs were simulated and compared to observed data.

For further references please refer to the database for observed data [[1](#reference)] and respective model repositories [[2-3](#reference)].



### References

[1] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

[2] [https://github.com/Open-Systems-Pharmacology/Rifampicin-Midazolam-DDI](https://github.com/Open-Systems-Pharmacology/Rifampicin-Midazolam-DDI)

[3] [https://github.com/Open-Systems-Pharmacology/Itraconazole-Midazolam-DDI](https://github.com/Open-Systems-Pharmacology/Itraconazole-Midazolam-DDI)



### Itraconazole-midazolam DDI

### Itraconazole-midazolam DDI

The itraconazole-midazolam interaction was evaluated using 5 clinical DDI studies involving studies with intravenous and oral administration of midazolam.

For further references please refer to the respective model repository [[2](#reference)] and the database of observed data [[3](#reference)].



#### References

[1] [https://github.com/Open-Systems-Pharmacology/Itraconazole-Midazolam-DDI](https://github.com/Open-Systems-Pharmacology/Itraconazole-Midazolam-DDI)

[2] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)



![001_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\001_plotComparisonTimeProfile.png)

![002_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\002_plotComparisonTimeProfile.png)

![003_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\003_plotComparisonTimeProfile.png)

![004_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\004_plotComparisonTimeProfile.png)

![005_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\005_plotComparisonTimeProfile.png)

![006_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\006_plotComparisonTimeProfile.png)

![007_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Itraconazole-midazolam_DDI\007_plotComparisonTimeProfile.png)

### Rifampicin-midazolam DDI

### Rifampicin-midazolam DDI

The rifampicin-midazolam interaction was evaluated using 14 clinical DDI studies involving studies with intravenous and oral administration of midazolam.

For further references please refer to the respective model repository [[2](#reference)] and the database of observed data [[3](#reference)].



#### References

[1] [https://github.com/Open-Systems-Pharmacology/Rifampicin-Midazolam-DDI](https://github.com/Open-Systems-Pharmacology/Rifampicin-Midazolam-DDI)

[2] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

![001_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\001_plotComparisonTimeProfile.png)

![002_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\002_plotComparisonTimeProfile.png)

![003_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\003_plotComparisonTimeProfile.png)

![004_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\004_plotComparisonTimeProfile.png)

![005_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\005_plotComparisonTimeProfile.png)

![006_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\006_plotComparisonTimeProfile.png)

![007_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\007_plotComparisonTimeProfile.png)

![008_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\008_plotComparisonTimeProfile.png)

![009_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\009_plotComparisonTimeProfile.png)

![010_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\010_plotComparisonTimeProfile.png)

![011_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\011_plotComparisonTimeProfile.png)

![012_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\012_plotComparisonTimeProfile.png)

![013_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\013_plotComparisonTimeProfile.png)

![014_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\014_plotComparisonTimeProfile.png)

![015_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\015_plotComparisonTimeProfile.png)

![016_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\016_plotComparisonTimeProfile.png)

![017_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\017_plotComparisonTimeProfile.png)

![018_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\018_plotComparisonTimeProfile.png)

![019_plotComparisonTimeProfile.png](images\02_Appendix\Evaluation_of_DDI_network__coupled_perpetrator___victim_simulations_\Rifampicin-midazolam_DDI\019_plotComparisonTimeProfile.png)

## Evaluation of victim and perpetrator models

## Evaluation of victim and perpetrator models

The PBPK models of **midazolam**, **itraconazole**, and **rifampicin** were developed with clinical data covering a broad dosing range possible for intravenous as well as oral, single and multiple dose administration of all drugs. Plasma or serum concentrations, fractions excreted to urine or bile, and other clinical measurements were included for model development whenever available. 

For further references please refer to the database for observed data [[1](#reference)] and respective model repositories [[2-4](#reference)].



### References

[1] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

[2] [https://github.com/Open-Systems-Pharmacology/Midazolam-Model](https://github.com/Open-Systems-Pharmacology/Midazolam-Model)

[3] [https://github.com/Open-Systems-Pharmacology/Itraconazole-Model](https://github.com/Open-Systems-Pharmacology/Itraconazole-Model)

[4] [https://github.com/Open-Systems-Pharmacology/Rifampicin-Model](https://github.com/Open-Systems-Pharmacology/Rifampicin-Model)





### Itraconazole Model

### Itraconazole model

Itraconazole is a triazole agent prescribed for the treatment of fungal infections. It is predominantly metabolized by CYP3A4, resulting in the sequential formation of several metabolites, starting with the major metabolite hydroxy-itraconazole, followed by keto-itraconazole and N-desalkyl-itraconazole. All three metabolites are further metabolized by CYP3A4. Parent and metabolites are reported to competitively inhibit CYP3A4 and hence, are included in the PBPK model. The metabolites inhibit their own formation and itraconazole inhibits further conversion of its metabolites by CYP3A4. Itraconazole has been proposed as one of the most appropriate strong CYP3A4 inhibitors for clinical DDI studies, to replace the no longer recommended CYP3A4 inhibitor drug ketoconazole. 

The presented model is based on the published model by Hanke et al. 2018 [[1](#reference)]. It model was built and evaluated using multiple clinical studies, covering a dose range from 100 to 400 mg after intravenous and oral administration in different formulations, administered under fasted conditions or together with food . Although the plasma concentrations of keto-itraconazole and N-desalkyl-itraconazole are lower than those of itraconazole and hydroxy-itraconazole, N-desalkyl-itraconazole is reported to be a very potent inhibitor in vitro, and integration of the two further metabolites into the model with their inhibitory effects helped to describe the strong nonlinearity and plasma accumulation of itraconazole. The model applies sequential metabolism of itraconazole to hydroxy-itraconazole to keto-itraconazole to N-desalkyl-itraconazole by CYP3A4, including competitive inhibition of CYP3A4 by the parent drug and all three metabolites, and glomerular filtration. Additional features of the model represent P-gp inhibition.

Please note that in comparison to the published version by Hanke et al. 2018 [[1](https://github.com/sfrechen/Itraconazole-Model/blob/master/README.md#reference)], the dissolution and solubility has been optimized and updated for the  administration of itraconazole capsules in fasted state.

For further references please refer to the itraconazole model repository [[2](#reference)] and the database of observed data [[3](#reference)].



#### References

[1] [Hanke N, Frechen S, Moj D, Britz H, Eissing T, Wendl T, Lehr T. PBPK models for CYP3A4 and P-gp DDI prediction: a modeling network of rifampicin, itraconazole, clarithromycin, midazolam, alfentanil and digoxin. CPT: Pharmacometrics & Systems Pharmacology (2018)](https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1002/psp4.12343)

[2] [https://github.com/Open-Systems-Pharmacology/Itraconazole-Model](https://github.com/Open-Systems-Pharmacology/Itraconazole-Model)

[3] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

#### Itraconazole input parameters

### Input parameters

# Compound: Itraconazole

## Parameters

| Name                                             | Value                   | Value Origin                                                 | Alternative                           | Default |      |
| ------------------------------------------------ | ----------------------- | ------------------------------------------------------------ | ------------------------------------- | ------- | ---- |
| Solubility at reference pH                       | 8 mg/l                  | Publication-Taupitz et al. 2013                              | Solution fasted (Taupitz et al. 2013) | True    |      |
| Reference pH                                     | 6.5                     | Publication-Taupitz et al. 2013                              | Solution fasted (Taupitz et al. 2013) | True    |      |
| Solubility at reference pH                       | 1.58 mg/l               |                                                              | Solution fed                          | False   |      |
| Reference pH                                     | 6.5                     |                                                              | Solution fed                          | False   |      |
| Solubility at reference pH                       | 0.9728307177 mg/l       | Parameter Identification-Parameter Identification-Value updated from 'Capsule fasted' on 2019-05-15 12:25 | Capsule fasted                        | False   |      |
| Reference pH                                     | 6.5                     |                                                              | Capsule fasted                        | False   |      |
| Solubility at reference pH                       | 0.7 mg/l                |                                                              | Capsule fed                           | False   |      |
| Reference pH                                     | 6.5                     |                                                              | Capsule fed                           | False   |      |
| Lipophilicity                                    | 4.624 Log Units         | Parameter Identification-Fit                                 | Fit                                   | True    |      |
| Fraction unbound (plasma, reference value)       | 0.6016197247 %          |                                                              | Templeton, 2008                       | True    |      |
| Permeability                                     | 0.1111202419 cm/min     | Parameter Identification-Fit                                 | Fit                                   | False   |      |
| Specific intestinal permeability (transcellular) | 5.3261558344E-05 dm/min | Parameter Identification-Fit                                 | Fit                                   | True    |      |
| Cl                                               | 2                       |                                                              |                                       |         |      |
| Is small molecule                                | Yes                     |                                                              |                                       |         |      |
| Molecular weight                                 | 705.633 g/mol           |                                                              |                                       |         |      |
| Plasma protein binding partner                   | Albumin                 |                                                              |                                       |         |      |
| Enable supersaturation                           | No                      |                                                              |                                       |         |      |

## Calculation methods

| Name                    | Value               |      |
| ----------------------- | ------------------- | ---- |
| Partition coefficients  | Rodgers and Rowland |      |
| Cellular permeabilities | PK-Sim Standard     |      |

## Processes

### Metabolizing Enzyme: CYP3A4-Isoherranen 2004

Molecule: CYP3A4
Metabolite: Hydroxy-Itraconazole

#### Parameters

| Name                             | Value                          | Value Origin                 |      |
| -------------------------------- | ------------------------------ | ---------------------------- | ---- |
| In vitro Vmax/recombinant enzyme | 0.27 pmol/min/pmol rec. enzyme | Publication-Isoherranen 2004 |      |
| Km                               | 2.0688492598 nmol/l            | Publication-Isoherranen 2004 |      |
| kcat                             | 0.0402937875 1/min             | Unknown                      |      |

### Systemic Process: Glomerular Filtration-GFR

Species: Human

#### Parameters

| Name         | Value | Value Origin                 |      |
| ------------ | ----: | ---------------------------- | ---- |
| GFR fraction |     1 | Publication-Isoherranen 2004 |      |

### Inhibition: CYP3A4-Isoherranen, 2004

Molecule: CYP3A4

#### Parameters

| Name | Value      | Value Origin                               |      |
| ---- | ---------- | ------------------------------------------ | ---- |
| Ki   | 1.3 nmol/l | Parameter Identification-Isoherranen, 2004 |      |

### Inhibition: ABCB1-Shityakov 2014

Molecule: ABCB1

#### Parameters

| Name | Value        | Value Origin               |      |
| ---- | ------------ | -------------------------- | ---- |
| Ki   | 0.008 µmol/l | Publication-Shityakov 2014 |      |

#### Itraconazole model evaluation

### Model evaluation

##### Itraconazole concentration-time profiles

#### Concentration-time profiles

![001_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\001_plotTimeProfile.png)

![002_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\002_plotTimeProfile.png)

![003_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\003_plotTimeProfile.png)

![004_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\004_plotTimeProfile.png)

![005_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\005_plotTimeProfile.png)

![006_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\006_plotTimeProfile.png)

![007_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\007_plotTimeProfile.png)

![008_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\008_plotTimeProfile.png)

![009_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\009_plotTimeProfile.png)

![010_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\010_plotTimeProfile.png)

![011_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\011_plotTimeProfile.png)

![012_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\012_plotTimeProfile.png)

![013_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\013_plotTimeProfile.png)

![014_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\014_plotTimeProfile.png)

![015_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\015_plotTimeProfile.png)

![016_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\016_plotTimeProfile.png)

![017_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\017_plotTimeProfile.png)

![018_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_concentration-time_profiles\018_plotTimeProfile.png)

##### Itraconazole overall goodness of fit

#### Overall goodness of fit

![001_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\001_plotGOFMergedResiduals.png)

![002_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\002_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.505025 

![004_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\004_plotGOFMergedResiduals.png)

![005_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\005_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.503085 

![007_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\007_plotGOFMergedResiduals.png)

![008_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\008_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.658847 

![010_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\010_plotGOFMergedResiduals.png)

![011_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Itraconazole_Model\Itraconazole_model_evaluation\Itraconazole_overall_goodness_of_fit\011_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.529602 

### Midazolam Model

### Midazolam model

Midazolam is a widely-used sedative, approved as premedication before surgical interventions. It is almost exclusively metabolized by CYP3A4, turning it into a sensitive probe and victim drug for the investigation of in vivo CYP3A4 activity. Midazolam shows substantial first pass metabolism, resulting in a bioavailability of under 50%. Less than 1% of a midazolam dose is excreted unchanged in urine.

The herein presented model represents an update of the midazolam model publisdhed by Hanke et al. [[1](#reference)]. The model has been  developed using in particular published pharmacokinetic clinical data by Hohmann et al. 2015 [[2](#reference)], Hyland et al. 2009 [[3](#reference)] and Thummel et al. 1996 [[4](#reference)]. It has then been evaluated by comparing observed data to simulations of a large number of clinical studies covering a dose range from 0.05 mg/kg to 20 mg after intravenous and oral administrations. Furthermore, it has been evaluated within a CYP3A4 DDI modeling network as a victim drug. 

Model features include:

- CYP3A4 metabolism
- (direct) UGT1A4 metabolism
- excretion into urine via glomerular filtration
- a decrease in the permeability between the intracellular and interstitial space  (parameters "P (intracellular->interstitial)" and "P (interstitial->intracellular)") in intestinal mucosa to optimize quantitatively the extent of gut wall metabolism
- and binding to a hypothetical binding partner in the brain to optimize a late redistribution phase in midazolam plasma concentrations.

For further references please refer to the midazolam model repository [[5](#reference)] and the database of observed data [[6](#reference)].



#### References

[1] [Hanke N, Frechen S, Moj D, Britz H, Eissing T, Wendl T, Lehr T. PBPK models for CYP3A4 and P-gp DDI prediction: a modeling network of rifampicin, itraconazole, clarithromycin, midazolam, alfentanil and digoxin. CPT: Pharmacometrics & Systems Pharmacology (2018)](https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1002/psp4.12343)

[2] [Hohmann N, Kocheise F, Carls A, Burhenne J, Haefeli WE, Mikus G. Midazolam microdose to determine systemic and pre-systemic metabolic CYP3A activity in humans. Br J Clin Pharmacol (2015)](https://doi.org/10.1111/bcp.12502)

[3] [Hyland R, Osborne T, Payne A, Kempshall S, Logan YR, Ezzeddine K, Jones B. In vitro and in vivo glucuronidation of midazolam in humans. Br J Clin Pharmacol (2009)](https://doi.org/10.1111/j.1365-2125.2009.03386.x)

[4] [Thummel KE, O'Shea D, Paine MF, Shen DD, Kunze KL, Perkins JD, Wilkinson GR. Oral first-pass elimination of midazolam involves both gastrointestinal and hepatic CYP3A-mediated metabolism. Clin Pharmacol Ther (1996)](https://doi.org/10.1016/S0009-9236(96)90177-0)

[5] [https://github.com/Open-Systems-Pharmacology/Midazolam-Model](https://github.com/Open-Systems-Pharmacology/Midazolam-Model)

[6] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

#### Midazolam input parameters

### Input parameters

# Compound: Midazolam

## Parameters

| Name                                             | Value                   | Value Origin                                                 | Alternative                             | Default |      |
| ------------------------------------------------ | ----------------------- | ------------------------------------------------------------ | --------------------------------------- | ------- | ---- |
| Solubility at reference pH                       | 0.13 mg/ml              | Publication-Heikkinen 2012                                   | Aqueous solubility                      | False   |      |
| Reference pH                                     | 5                       | Publication-Heikkinen 2012                                   | Aqueous solubility                      | False   |      |
| Solubility at reference pH                       | 0.049 mg/ml             | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 | FaSSIF                                  | True    |      |
| Reference pH                                     | 6.5                     | Publication-Heikkinen 2012                                   | FaSSIF                                  | True    |      |
| Solubility at reference pH                       | 0.09 mg/ml              | Publication-Heikkinen 2012                                   | FeSSIF                                  | False   |      |
| Reference pH                                     | 5                       | Publication-Heikkinen 2012                                   | FeSSIF                                  | False   |      |
| Lipophilicity                                    | 2.8972038771 Log Units  | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 | Optimized                               | True    |      |
| Lipophilicity                                    | 3.53 Log Units          | Publication-Simcyp                                           | LogP (Simcyp)                           | False   |      |
| Lipophilicity                                    | 3 Log Units             | Publication-Dagenais 2009                                    | LogP (experimental) (Dagenais)          | False   |      |
| Lipophilicity                                    | 3.37 Log Units          | Publication-GastroPlus                                       | LogP (GastroPlus)                       | False   |      |
| Lipophilicity                                    | 3.1 Log Units           | Publication-Rodgers and Rowland                              | LogP (experimental) (Rodgers & Rowland) | False   |      |
| Fraction unbound (plasma, reference value)       | 0.032                   | Publication-Simcyp                                           | Simcyp                                  | False   |      |
| Fraction unbound (plasma, reference value)       | 0.031                   | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 | Gertz et al. 2010                       | True    |      |
| Fraction unbound (plasma, reference value)       | 0.022                   | Publication-Lown et al. 1995                                 | Lown et al. 1995                        | False   |      |
| Fraction unbound (plasma, reference value)       | 0.016                   |                                                              | Björkman et al. 2001 (men)              | False   |      |
| Fraction unbound (plasma, reference value)       | 0.02                    |                                                              | Björkman et al. 2001 (women)            | False   |      |
| Specific intestinal permeability (transcellular) | 0.00015549970673 cm/min | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 | Optimized                               | True    |      |
| Cl                                               | 1                       |                                                              |                                         |         |      |
| F                                                | 1                       |                                                              |                                         |         |      |
| Is small molecule                                | Yes                     |                                                              |                                         |         |      |
| Molecular weight                                 | 325.78 g/mol            |                                                              |                                         |         |      |
| Plasma protein binding partner                   | Albumin                 |                                                              |                                         |         |      |
| Enable supersaturation                           | No                      |                                                              |                                         |         |      |

## Calculation methods

| Name                    | Value               |      |
| ----------------------- | ------------------- | ---- |
| Partition coefficients  | Rodgers and Rowland |      |
| Cellular permeabilities | PK-Sim Standard     |      |

## Processes

### Specific Binding: GABRG2-Buhr 1997

Molecule: GABRG2

#### Parameters

| Name | Value      | Value Origin                                                 |      |
| ---- | ---------- | ------------------------------------------------------------ | ---- |
| koff | 1 1/min    | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |
| Kd   | 1.8 nmol/l |                                                              |      |

### Metabolizing Enzyme: CYP3A4-Patki et al. 2003

Molecule: CYP3A4
Metabolite: 1-OH-Midazolam

#### Parameters

| Name                               | Value                         | Value Origin |      |
| ---------------------------------- | ----------------------------- | -----------: | ---- |
| In vitro Vmax for liver microsomes | 0.18 nmol/min/mg mic. protein |              |      |
| Km                                 | 3.9 µmol/l                    |              |      |

### Metabolizing Enzyme: CYP3A4-Simcyp

Molecule: CYP3A4
Metabolite: 1-OH-Midazolam

#### Parameters

| Name                             | Value                          | Value Origin |      |
| -------------------------------- | ------------------------------ | -----------: | ---- |
| In vitro Vmax/recombinant enzyme | 2.16 pmol/min/pmol rec. enzyme |              |      |
| Km                               | 2.16 µmol/l                    |              |      |

### Metabolizing Enzyme: UGT1A4-Optimized

Molecule: UGT1A4
Metabolite: Midazolam-N-glucuronide

#### Parameters

| Name                                        | Value                        | Value Origin                                                 |      |
| ------------------------------------------- | ---------------------------- | ------------------------------------------------------------ | ---- |
| In vitro Vmax for liver microsomes          | 276 pmol/min/mg mic. protein |                                                              |      |
| Content of CYP proteins in liver microsomes | 58 pmol/mg mic. protein      | Publication-Achour 2014                                      |      |
| Km                                          | 37.8 µmol/l                  | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |
| kcat                                        | 3.5911771641 1/min           | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |

### Systemic Process: Glomerular Filtration-Optimized

Species: Human

#### Parameters

| Name         |        Value | Value Origin                                                 |      |
| ------------ | -----------: | ------------------------------------------------------------ | ---- |
| GFR fraction | 0.6401025724 | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |

### Metabolizing Enzyme: CYP3A4-Optimized

Molecule: CYP3A4

#### Parameters

| Name                               | Value                        | Value Origin                                                 |      |
| ---------------------------------- | ---------------------------- | ------------------------------------------------------------ | ---- |
| In vitro Vmax for liver microsomes | 850 pmol/min/mg mic. protein |                                                              |      |
| Km                                 | 4 µmol/l                     | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |
| kcat                               | 8.7607941215 1/min           | Parameter Identification-Parameter Identification-Value updated from 'PI Hohmann iv+po, Hyland feUr MDZG, Thummel feUr unchanged - Pint' on 2019-04-09 16:10 |      |

### Metabolizing Enzyme: CYP3A4-Galentin et al. 2004

Molecule: CYP3A4
Metabolite: 1-OH-Midazolam

#### Parameters

| Name                             | Value                          | Value Origin |      |
| -------------------------------- | ------------------------------ | -----------: | ---- |
| In vitro Vmax/recombinant enzyme | 1.96 pmol/min/pmol rec. enzyme |              |      |
| Km                               | 2.69 µmol/l                    |              |      |

### Metabolizing Enzyme: CYP3A4-Ito et al. 2003

Molecule: CYP3A4
Metabolite: 1-OH-Midazolam

#### Parameters

| Name                               | Value                         | Value Origin |      |
| ---------------------------------- | ----------------------------- | -----------: | ---- |
| In vitro Vmax for liver microsomes | 4.41 nmol/min/mg mic. protein |              |      |
| Km                                 | 3.8 µmol/l                    |              |      |

### Metabolizing Enzyme: CYP3A4-GastroPlus

Molecule: CYP3A4

#### Parameters

| Name                               | Value                        | Value Origin |      |
| ---------------------------------- | ---------------------------- | -----------: | ---- |
| In vitro Vmax for liver microsomes | 850 pmol/min/mg mic. protein |              |      |
| Km                                 | 4 µmol/l                     |              |      |

### Metabolizing Enzyme: UGT1A4-Klieber et al. 2008

Molecule: UGT1A4

#### Parameters

| Name                                        | Value                        | Value Origin            |      |
| ------------------------------------------- | ---------------------------- | ----------------------- | ---- |
| In vitro Vmax for liver microsomes          | 276 pmol/min/mg mic. protein |                         |      |
| Content of CYP proteins in liver microsomes | 58 pmol/mg mic. protein      | Publication-Achour 2014 |      |
| Km                                          | 37.8 µmol/l                  |                         |      |

#### Midazolam model evaluation

### Model evaluation

##### Midazolam concentration-time profiles

#### Concentration-time profiles

![001_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\001_plotTimeProfile.png)

![002_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\002_plotTimeProfile.png)

![003_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\003_plotTimeProfile.png)

![004_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\004_plotTimeProfile.png)

![005_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\005_plotTimeProfile.png)

![006_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\006_plotTimeProfile.png)

![007_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\007_plotTimeProfile.png)

![008_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\008_plotTimeProfile.png)

![009_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\009_plotTimeProfile.png)

![010_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\010_plotTimeProfile.png)

![011_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\011_plotTimeProfile.png)

![012_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\012_plotTimeProfile.png)

![013_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\013_plotTimeProfile.png)

![014_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\014_plotTimeProfile.png)

![015_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\015_plotTimeProfile.png)

![016_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\016_plotTimeProfile.png)

![017_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\017_plotTimeProfile.png)

![018_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\018_plotTimeProfile.png)

![019_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\019_plotTimeProfile.png)

![020_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\020_plotTimeProfile.png)

![021_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\021_plotTimeProfile.png)

![022_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\022_plotTimeProfile.png)

![023_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\023_plotTimeProfile.png)

![024_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\024_plotTimeProfile.png)

![025_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\025_plotTimeProfile.png)

![026_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\026_plotTimeProfile.png)

![027_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\027_plotTimeProfile.png)

![028_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\028_plotTimeProfile.png)

![029_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\029_plotTimeProfile.png)

![030_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\030_plotTimeProfile.png)

![031_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\031_plotTimeProfile.png)

![032_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\032_plotTimeProfile.png)

![033_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\033_plotTimeProfile.png)

![034_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_concentration-time_profiles\034_plotTimeProfile.png)

##### Midazolam overall goodness of fit

#### Overall goodness of fit

![001_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_overall_goodness_of_fit\001_plotGOFMergedResiduals.png)

![002_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Midazolam_Model\Midazolam_model_evaluation\Midazolam_overall_goodness_of_fit\002_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.457709 

### Rifampicin Model

### Rifampicin model

Rifampicin is an antibiotic used for the treatment of mycobacterium infections, including tuberculosis and leprosy. For the investigation of DDIs, rifampicin is an established potent inducer of multiple drug metabolizing enzymes (CYP3A4, CYP2B6, CYP2C8, CYP2C9, CYP2C19) and transporters (P-gp, MRP2, MRP3, MRP4, OATP1A2). In addition to its inducing capabilities, rifampicin also competitively inhibits enzymes and transporters like CYP3A4, P-gp, OATP1B1 and OATP1B3. This qualification report focuses on the CYP3A4  effects of rifampicin. While induction by rifampicin involves upregulation of enzymes and transporters and therefore takes several days to fully develop, competitive inhibition has an instantaneous effect and is strongest at the time of highest exposure to the inhibitor. As a consequence, the effects of rifampicin caused via competitive inhibition are most prominent 1 - 2 h after its oral administration and of relatively short duration. These complex and opposing actions of rifampicin demand careful consideration of the timing of interacting drugs during clinical studies and modeling. 

The presented model is based on the published model by Hanke et al. 2018 [[1](#reference)]. It was built and evaluated using several clinical studies, covering a dosing range from 300 to 600 mg after intravenous and oral administration. Integrating and testing processes that were described as vital to the pharmacokinetics of rifampicin resulted in a final model that applies transport by OATP1B1, metabolism by arylacetamide deacetylase (AADAC), transport by P-gp and glomerular filtration. Furthermore, auto-induction of OATP1B1, AADAC and P-gp expression has been incorporated . Studies that measured pharmacokinetic profiles of rifampicin at different days of a 600 mg po once daily regimen indicate that the clearance of rifampicin increases over time due to auto-induction of elimination processes. Therefore, induction of OATP1B1, AADAC and P-gp expression is included in the rifampicin model. The hypothesis of AADAC induction by rifampicin was based on the fact, that 

- rifampicin induces its own metabolism, and 
- B-esterases are inducible by rifampicin via PXR. 

For applications in the context of DDI, the model features CYP3A4 induction, transient competitive 
CYP3A4 inhibition. Additional features represent P-gp induction and inhibition.

For further references please refer to the rifampicin model repository [[2](#reference)] and the database of observed data [[3](#reference)].



#### References

[1] [Hanke N, Frechen S, Moj D, Britz H, Eissing T, Wendl T, Lehr T. PBPK models for CYP3A4 and P-gp DDI prediction: a modeling network of rifampicin, itraconazole, clarithromycin, midazolam, alfentanil and digoxin. CPT: Pharmacometrics & Systems Pharmacology (2018)](https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1002/psp4.12343)

[2] [https://github.com/Open-Systems-Pharmacology/Rifampicin-Model](https://github.com/Open-Systems-Pharmacology/Rifampicin-Model)

[3] [https://github.com/Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)

#### Rifampicin input parameters

### Input parameters

# Compound: Rifampicin

## Parameters

| Name                                             | Value           | Value Origin                                            | Alternative        | Default |      |
| ------------------------------------------------ | --------------- | ------------------------------------------------------- | ------------------ | ------- | ---- |
| Solubility at reference pH                       | 2800 mg/l       | Database-DrugBank                                       | Aqueous solubility | True    |      |
| Reference pH                                     | 7               | Database-DrugBank                                       | Aqueous solubility | True    |      |
| Solubility at reference pH                       | 450 mg/l        |                                                         | FaSSIF             | False   |      |
| Reference pH                                     | 6.5             |                                                         | FaSSIF             | False   |      |
| Solubility at reference pH                       | 1600 mg/l       |                                                         | FeSSIF             | False   |      |
| Reference pH                                     | 5               |                                                         | FeSSIF             | False   |      |
| Solubility at reference pH                       | 1211.88 mg/l    |                                                         | Optimized          | False   |      |
| Reference pH                                     | 7               |                                                         | Optimized          | False   |      |
| Lipophilicity                                    | 2.5 Log Units   | Publication-Parameter Identification-Hanke et al. 2018  | Optimized          | True    |      |
| Lipophilicity                                    | 2.7 Log Units   |                                                         | DrugBank logP      | False   |      |
| Fraction unbound (plasma, reference value)       | 17 %            | Publication-Other-Templeton 2011 (equilibrium dialysis) | Templeton 2011     | True    |      |
| Fraction unbound (plasma, reference value)       | 16 %            |                                                         | Baneyx 2013        | False   |      |
| Fraction unbound (plasma, reference value)       | 11 %            |                                                         | Boman 1974         | False   |      |
| Specific intestinal permeability (transcellular) | 1.24E-05 cm/min | Publication-Parameter Identification-Hanke et al. 2018  | Optimized          | True    |      |
| Is small molecule                                | Yes             |                                                         |                    |         |      |
| Molecular weight                                 | 822.94 g/mol    |                                                         |                    |         |      |
| Plasma protein binding partner                   | Albumin         |                                                         |                    |         |      |
| Enable supersaturation                           | No              |                                                         |                    |         |      |

## Calculation methods

| Name                    | Value               |      |
| ----------------------- | ------------------- | ---- |
| Partition coefficients  | Rodgers and Rowland |      |
| Cellular permeabilities | PK-Sim Standard     |      |

## Processes

### Metabolizing Enzyme: AADAC-Nakajima 2011

Molecule: AADAC

#### Parameters

| Name                 | Value          | Value Origin                                           |      |
| -------------------- | -------------- | ------------------------------------------------------ | ---- |
| Enzyme concentration | 1 µmol/l       |                                                        |      |
| Vmax                 | 6.5 µmol/l/min |                                                        |      |
| Km                   | 195.1 µmol/l   |                                                        |      |
| kcat                 | 9.865 1/min    | Publication-Parameter Identification-Hanke et al. 2018 |      |

### Transport Protein: P-gp-Collett 2004

Molecule: P-gp

#### Parameters

| Name                      | Value           | Value Origin                                           |      |
| ------------------------- | --------------- | ------------------------------------------------------ | ---- |
| Transporter concentration | 60 nmol/l       |                                                        |      |
| Vmax                      | 2.87 µmol/l/min |                                                        |      |
| Km                        | 55 µmol/l       |                                                        |      |
| kcat                      | 0.6088 1/min    | Publication-Parameter Identification-Hanke et al. 2018 |      |

### Transport Protein: OATP1B1-Tirona 2003

Molecule: OATP1B1

#### Parameters

| Name                      | Value            | Value Origin                                           |      |
| ------------------------- | ---------------- | ------------------------------------------------------ | ---- |
| Transporter concentration | 109.6 µmol/l     |                                                        |      |
| Vmax                      | 0.372 µmol/l/min |                                                        |      |
| Km                        | 1.5 µmol/l       |                                                        |      |
| kcat                      | 7.796 1/min      | Publication-Parameter Identification-Hanke et al. 2018 |      |

### Systemic Process: Glomerular Filtration-GFR

Species: Human

#### Parameters

| Name         | Value | Value Origin                             |      |
| ------------ | ----: | ---------------------------------------- | ---- |
| GFR fraction |     1 | Publication-Assumption-Hanke et al. 2018 |      |

### Inhibition: CYP3A4-Kajosaari 2005

Molecule: CYP3A4

#### Parameters

| Name | Value       | Value Origin                      |      |
| ---- | ----------- | --------------------------------- | ---- |
| Ki   | 18.5 µmol/l | Publication-Kajosaari et al. 2005 |      |

### Inhibition: P-gp-Reitman 2011

Molecule: P-gp

#### Parameters

| Name | Value      | Value Origin                                                 |      |
| ---- | ---------- | ------------------------------------------------------------ | ---- |
| Ki   | 169 µmol/l | Publication-Assumption-Reitman 2011 (IC50 = Ki (169 µM / (1+ (0.1 µM / 177 µM) ) |      |

### Induction: CYP3A4-Templeton 2011

Molecule: CYP3A4

#### Parameters

| Name | Value       | Value Origin                                       |      |
| ---- | ----------- | -------------------------------------------------- | ---- |
| EC50 | 0.34 µmol/l | Publication-Templeton 2011 (weighted mean for FHH) |      |
| Emax | 9           | Publication-Templeton 2011 (weighted mean for FHH) |      |

### Induction: P-gp-Greiner 1999

Molecule: P-gp

#### Parameters

| Name | Value       | Value Origin                               |      |
| ---- | ----------- | ------------------------------------------ | ---- |
| EC50 | 0.34 µmol/l | Publication-Assumption-Hanke et al. 2018   |      |
| Emax | 2.5         | Publication-Assumption-Greiner et al. 1999 |      |

### Induction: OATP1B1-Dixit 2007

Molecule: OATP1B1

#### Parameters

| Name | Value       | Value Origin                                           |      |
| ---- | ----------- | ------------------------------------------------------ | ---- |
| EC50 | 0.34 µmol/l | Publication-Assumption-Hanke et al. 2018               |      |
| Emax | 0.383       | Publication-Parameter Identification-Hanke et al. 2018 |      |

### Induction: AADAC-Assumed

Molecule: AADAC

#### Parameters

| Name | Value       | Value Origin                                           |      |
| ---- | ----------- | ------------------------------------------------------ | ---- |
| EC50 | 0.34 µmol/l | Publication-Assumption-Hanke et al. 2018               |      |
| Emax | 0.985       | Publication-Parameter Identification-Hanke et al. 2018 |      |

#### Rifampicin model evaluation

### Model evaluation

##### Rifampicin concentration-time profiles

#### Concentration-time profiles

![001_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\001_plotTimeProfile.png)

![002_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\002_plotTimeProfile.png)

![003_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\003_plotTimeProfile.png)

![004_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\004_plotTimeProfile.png)

![005_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\005_plotTimeProfile.png)

![006_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\006_plotTimeProfile.png)

![007_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\007_plotTimeProfile.png)

![008_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\008_plotTimeProfile.png)

![009_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\009_plotTimeProfile.png)

![010_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\010_plotTimeProfile.png)

![011_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\011_plotTimeProfile.png)

![012_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\012_plotTimeProfile.png)

![013_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\013_plotTimeProfile.png)

![014_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\014_plotTimeProfile.png)

![015_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\015_plotTimeProfile.png)

![016_plotTimeProfile.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_concentration-time_profiles\016_plotTimeProfile.png)

##### Rifampicin overall goodness of fit

#### Overall goodness of fit

![001_plotGOFMergedResiduals.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_overall_goodness_of_fit\001_plotGOFMergedResiduals.png)

![002_plotGOFMergedPredictedVsObserved.png](images\02_Appendix\Evaluation_of_victim_and_perpetrator_models\Rifampicin_Model\Rifampicin_model_evaluation\Rifampicin_overall_goodness_of_fit\002_plotGOFMergedPredictedVsObserved.png)

GMFE = 1.489023 

## Mathematical Implementation of drug-drug interactions

## Mathematical Implementation of drug-drug interactions



### DDI modeling: Competitive inhibition 

A detailed representation of the mathematical implementation of competitive enzyme inhibition  can be found in the OSP manual [[1](#references)].



### DDI modeling: Mechanism-based inhibition 

A detailed representation of the mathematical implementation of mechanism-based enzyme inhibition  can be found in the OSP manual [[2](#references)].



### DDI modeling: Induction 

A detailed representation of the mathematical implementation of enzyme induction can be found in the OSP manual [[3](#references)].

 

### References

[1] [https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#competitive-inhibition-simple-setting-with-one-inhibitor](https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#competitive-inhibition-simple-setting-with-one-inhibitor)

[2] [https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#irreversible-inhibition](https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#irreversible-inhibition)

[3] [https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#enzyme-induction](https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-compounds-defining-inhibition-induction-processes#enzyme-induction)



## Open Systems Pharmacology Suite (OSPS) Introduction

## Open Systems Pharmacology Suite (OSPS) Introduction

OSPS is a tool for PBPK modelling and simulation of drugs in laboratory animals and humans. PK-Sim and MoBi are part of the Open Systems Pharmacology Suite (OSPS) [[1](#references)]. PK-Sim® is based on a generic PBPK-model with 18 organs and tissues. Represented organs/tissues include arterial and venous blood, adipose tissue (separable adipose, excluding yellow marrow), brain, lung, bone (including yellow marrow), gonads, heart, kidneys, large intestine, liver, muscle, portal vein, pancreas, skin, small intestine, spleen and stomach, as shown in Figure 1.

Each organ consists of four sub-compartments namely the plasma, red blood cells (which together build the vascular space), interstitial space, and cellular space. Distribution between the plasma and red blood cells as well as between the interstitial and cellular compartments can be permeability-limited. In the brain, the permeation barrier is located between the vascular and the interstitial space. PK-Sim® estimates model parameters (intestinal permeability [[2](#references)], organ partition coefficients [[3,4](#references)], and permeabilities) from physico-chemical properties of compounds (molecular weight, pKa, ace/base properties) and the composition of each tissue compartment (lipids, water and proteins). Partition coefficients can be calculated using a variety of methods available in PK-Sim®, for example the internal PK-Sim® method [[3,4](#references)] or that of Rodgers and Rowland [[5-7](#references)]. 

Physiological databases included in the software incorporate the dependencies of organ weights, organ blood flows and gastrointestinal parameters (gastrointestinal length, radius of each section, intestinal surface area [[2]) with the user-defined body weight and height of the individual [[8](#references)]. Thereby, PK Sim® allows generating realistic virtual populations. For a detailed description of the PBPK model structure implemented in PK Sim®, see Willmann et al. [[2,4,8,9](#references)] or the Open Systems Pharmacology (OSP) Suite homepage (<http://www.open-systems-pharmacology.org/>).

 

**Figure** **1: Structure of the Whole Body PBPK Model integrated in PK-Sim®**

![generic PBPK model](images/PK-Sim_PBPK_generic_model_scheme.png)

### References

[1] [www.open-systems-pharmacology.org](http://www.open-systems-pharmacology.org/)

[2] [Willmann S, Schmitt W, Keldenich J, Lippert J, Dressman JB. A physiological model for the estimation of the fraction dose absorbed in humans.J Med Chem. 2004 Jul 29;47(16):4022-31.](https://www.ncbi.nlm.nih.gov/pubmed/15267240)

[3] [Haerter MW, K.J., Schmitt W, *Estimation of physicochemical and ADME parameters.* , in *Handbook of Combinatorial Chemistry: Drugs, Catalysts, Materials*, H.W. Nicolaou KC HR, Editor. 2002, Wiley VCH Verlag GmbH: Weinheim, Germany. p. 743-60.](https://onlinelibrary.wiley.com/doi/pdf/10.1002/3527603034.ch26)

[4] [Willmann S, Lippert J, Schmitt W. From physicochemistry to absorption and distribution: predictive mechanistic modelling and computational tools. Expert Opin Drug Metab Toxicol. 2005 Jun;1(1):159-68.](https://www.ncbi.nlm.nih.gov/pubmed/16922658)

[5] [Rodgers, T, D. Leahy, and M. Rowland. Physiologically based pharmacokinetic modeling 1: predicting the tissue distribution of moderate-to-strong bases. J Pharm Sci. 2005 Jun;94(6):1259-76.](https://www.ncbi.nlm.nih.gov/pubmed/15858854)

[6] [Rodgers T, Rowland M. Physiologically based pharmacokinetic modelling 2: predicting the tissue distribution of acids, very weak bases, neutrals and zwitterions. J Pharm Sci. 2006 Jun;95(6):1238-57.](https://www.ncbi.nlm.nih.gov/pubmed/16639716)

[7] [Rodgers T, Rowland M. Mechanistic approaches to volume of distribution predictions: understanding the processes. Pharm Res. 2007 May;24(5):918-33. Epub 2007 Mar 20.](https://www.ncbi.nlm.nih.gov/pubmed/17372687)

[8] [Willmann S, Höhn K, Edginton A, Sevestre M, Solodenko J, Weiss W, Lippert J, Schmitt W. Development of a physiology-based whole-body population model for assessing the influence of individual variability on the pharmacokinetics of drugs. J Pharmacokinet Pharmacodyn. 2007 Jun;34(3):401-31. Epub 2007 Mar 13.](https://www.ncbi.nlm.nih.gov/pubmed/17431751)

[9] [Willmann S, Lippert J, Sevestre M, Solodenko J, Schmitt W. PK-Sim®: a physiologically based pharmacokinetic ‘whole-body’ model. Biosilico 2003.1(4):121-24.](https://www.sciencedirect.com/science/article/pii/S1478538203023424?via%3Dihub)

