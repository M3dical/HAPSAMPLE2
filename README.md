# HAP-SAMPLE2

Code for Data-based Resampling for Association Studies with Admixture

### Introduction

Given the advancements in statistical genetics, simulation methods are becoming increasingly crucial for enhancing genotype-phenotype association studies (Wright et al., 2007; Su et al., 2011). Simulated datasets are invaluable for conducting power analyses and evaluating various statistical learning approaches (Yuan et al., 2012). In response to these needs, we introduce HAP-SAMPLE2, an advanced iteration that builds upon the foundational methodologies of HAP-SAMPLE (Wright et al., 2007). HAP-SAMPLE2 represents a significant evolution in simulation technology, offering robust tools for generating datasets for admixed populations. By resampling from actual datasets and incorporating both common and rare variants, HAP-SAMPLE2 facilitates the creation of comprehensive phenotypic data suitable for case/control studies and quantitative trait analyses.
### Installation

Please install the following R packages. Version numbers for these packages were used in the HAP-SAMPLE2 vignette.

```bash
parallel      4.3.1
vcfR          1.15.0
Matrix        1.6-5
slam          0.1-50
matrixStats   1.0.0
doParallel    1.0.17
foreach       1.5.2
DirichletReg  0.7-1
LEA           3.12.2
ggplot2       3.4.2
```

### Example Usage

Please check the HAP-SAMPLE2 vignette in the /vignettes folder to run example simulations using the data in the /inst/extdata folder. Copy the files to your working directory and update the directories in the vignette.
