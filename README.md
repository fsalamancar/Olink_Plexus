# 🧬 PLEXUS Biomarker Discovery

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)  
![Status](https://img.shields.io/badge/status-Research%20Project-yellow)


This repository contains the analysis of the PLEXUS dataset aimed at identifying protein-based blood biomarkers associated with Inflammatory Bowel Disease (IBD), using Olink proteomics technology.

The project focuses on two main components:

**Identifying differential biomarkers between:**
- 🔹 Ulcerative Colitis (UC) and Crohn's Disease (CD)

**Identifying biomarkers across CD subtypes based on disease location:**
- 🔹 Ileal
- 🔹 Colonic
- 🔹 Ileocolonic

## 📁 Project Structure

The repository is divided into two main folders:

```
├── Plexus_CDvariants          # CD subtypes comparison (ileal, colonic, ileocolonic)
└── Plexus_UCvsCD              # UC vs CD comparison
```

Each folder follows a common structure:

```
├── 1.Datos_Raw                        # Raw proteomics data (Olink)
├── 2.Datos_Limpios_UCvsCD            # Clean data for UC vs CD
├── 2.1.Datos_Limpios_CD_Comparissons # Clean data for CD subtype comparisons
├── 3.Previsualizaciones              # Exploratory plots and quality control
└── Output
    ├── data                          # Processed outputs
    └── tests                         # Statistical test results
```

## ⚙️ Methodology

The analysis follows a three-stage pipeline:

### 1️⃣ Covariate Selection
📄 `Feature_selection.ipynb`  
Performs initial data cleaning and uses classification models to select relevant covariates for downstream modeling.

### 2️⃣ Data Cleaning
📄 `Data_Cleaning.R`  
Processes and prepares the dataset, ensuring quality and consistency for association analysis.

### 3️⃣ Logistic Regression Modeling
📄 `Modeling.R`  
Runs protein-by-protein logistic regression to identify associations with the selected phenotype, adjusted by the previously selected covariates.

## 📈 Report 
The report is located here:
[Canva](https://www.canva.com/design/DAGpUaEiEHM/yiK0eYn5EvBvasw39ivyOQ/edit)

## 🧪 Technology Stack

- 🔬 Olink Explore – Protein biomarker quantification platform  
- 🐍 Python & 📘 R – Data wrangling, visualization, and modeling  
- 🧠 Feature selection using classification models  
- 📊 Logistic regression for association analysis

## 💾 Data
The data is shared via request

## 📦 Requirements

All required packages and dependencies are listed in `requirements.txt`.

## 👤 Author

**Francisco Salamanca**  
Bioinformatician | MSc in Bioinformatics  
Universidad Nacional de Colombia | Institute of Clinical Molecular Biology

[GitHub](https://github.com/fsalamancar) • [Website](https://fsalamancar.github.io/) • [LinkedIn](https://www.linkedin.com/in/fjosesala/) • [IKMB](https://www.ikmb.uni-kiel.de/people/francisco-salamanca/)
