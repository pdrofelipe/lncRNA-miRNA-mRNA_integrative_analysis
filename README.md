# Construction of lncRNA–miRNA–mRNA ceRNA Regulatory Networks

This repository contains the R pipeline used to construct competitive endogenous RNA (ceRNA) networks based on experimentally validated interactions. The workflow integrates differentially expressed miRNAs (DE-miRNAs) with their target mRNAs and lncRNAs using the `multiMiR` package, leveraging the **miRTarBase** database.

---

## 📋 Overview & Methodology

The regulatory network is built through a three-phase data integration approach:

1. **miRNA $\rightarrow$ mRNA Interactions:** Identification of target mRNAs regulated by differentially expressed miRNAs (DE-miRNAs).
2. **lncRNA $\rightarrow$ miRNA Interactions:** Identification of lncRNAs acting as molecular sponges for the selected miRNAs based on official HGNC nomenclature patterns (e.g., `LINC`, `-AS`, `HG`, `MALAT1`, `NEAT1`).
3. **ceRNA Network Integration:** Merging lncRNA–miRNA and miRNA–mRNA interactions to construct full **lncRNA–miRNA–mRNA** regulatory axes.

---

## 🛠️ Prerequisites & Dependencies

Ensure you have R (version 4.0 or higher) installed. The script automatically manages dependency installation via `BiocManager`.

### Required R Packages
- **[multiMiR](https://bioconductor.org/packages/multiMiR/)**: Access to microRNA-target databases.
- **[tidyverse](https://www.tidyverse.org/)**: Data manipulation and transformation.

To install dependencies manually, run:
```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("multiMiR", "tidyverse"))
