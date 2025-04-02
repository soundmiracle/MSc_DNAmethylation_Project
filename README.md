---

```
# 🧬 Genome-wide DNA Methylation Correlation Analysis

This repository contains R scripts and supporting files for my MSc project conducted at University College London. The project focused on analysing **long-range correlation structures of DNA methylation** across various tissues and ethnic groups.

---

## 📁 Project Structure

```
MSc_DNAmethylation_Project/
├── README.md
├── 01_preprocessing/
│   └── methylation_preprocess.R
├── 02_correlation/
│   └── correlation_calc.R
├── 03_HCR_detection/
│   └── hcr_analysis.R
├── 04_visualisation/
│   ├── heatmap.R
│   └── violin_plot.R
└── data/
    └── humanmethylation450_15017482_v1-2.csv.gz
```

---

## 📋 Project Summary

We aimed to quantify DNA methylation correlation patterns across the genome by:

- Aligning methylation data across 27 tissues/cell types.
- Computing pairwise CpG site correlations within 10kb genomic windows.
- Extracting **Highly Correlated Regions (HCRs)** of CpG sites.
- Comparing correlation conservation across datasets to infer biological consistency.

---

## 💻 How to Use

### Step 1: Preprocessing (Optional)
```r
source("01_preprocessing/methylation_preprocess.R")
```

### Step 2: Correlation Calculation
```bash
Rscript 02_correlation/correlation_calc.R /path/to/your/RDS/files/
```

This script loops through `.RDS` beta-value matrices and calculates pairwise CpG correlations within 10kb windows.

### Step 3: HCR Detection
```bash
Rscript 03_HCR_detection/hcr_analysis.R chromosome_number block_size
```
- `chromosome_number` — e.g. `"14"`
- `block_size` — in base pairs (e.g., `100000` or `1000000`)

### Step 4: Visualisation
```r
source("04_visualisation/heatmap.R")
source("04_visualisation/violin_plot.R")
```

Generates correlation heatmaps, violin plots of beta values, and inter-dataset comparison charts.

---

## 📦 Required R Packages

Make sure the following packages are installed:
- `ggplot2`
- `dplyr`
- `reshape2`
- `gtools`
- `gridExtra`
- `gplots`

Install them via:

```r
install.packages(c("ggplot2", "dplyr", "reshape2", "gtools", "gridExtra", "gplots"))
```

---

## 📈 Outputs

- `.RData` files: CpG-CpG correlation matrices
- `.txt` files: Summary stats per genomic block
- `.png` files: Heatmaps, density plots, violin plots
- `cor.datasets` list object: Collection of region-wise correlation matrices

---

## 📎 Notes

- Methylation matrices are `.RDS` format (CpG x sample).
- CpG manifest (`humanmethylation450_15017482_v1-2.csv.gz`) is required for annotation.
- All results are computed per chromosome and per block size (e.g. 100kb or 1Mb).
- Correlations are squared (r²) to assess strength, and average correlation per region is used as a measure of consistency.

---

## 👨‍🔬 Author

**Young Jun Kim**  
MSc in Genetics of Human Disease, University College London  
📧 [young.kim.22@ucl.ac.uk](mailto:young.kim.22@ucl.ac.uk)  
📞 +44 07377 613434  
🔗 [LinkedIn Profile](https://www.linkedin.com/in/youngjunkim95)

---

```

---
