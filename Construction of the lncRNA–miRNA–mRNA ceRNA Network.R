# Title: Construction of the lncRNA–miRNA–mRNA ceRNA Network
# Methodology: Integration of DE-miRNAs with validated mRNA/lncRNA targets 
# using multiMiR (mirtarbase)

# 1. Environment Setup -----------------------------------------------------

# Recommended: Use R Projects (.Rproj) instead of hardcoded setwd()
# setwd('C:/Users/use_your_own_pc')

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("multiMiR", "tidyverse"))

library(multiMiR)
library(tidyverse)

# 2. Data Import -----------------------------------------------------------

# Importing Differentially Expressed (DE) miRNA data
# Note: csv2 is used for Brazilian/European format (sep = ';', dec = ',')
degs_data <- read.csv2("DEMIRNA_INPUT.csv")
target_mirnas <- degs_data$Symbol

# 3. Phase 1: miRNA -> mRNA Interactions -----------------------------------
# Retrieving experimentally validated mRNA targets for the input miRNAs

res_mirna_deg <- get_multimir(
  org     = "hsa",
  target  = target_mirnas,
  table   = "mirtarbase", # Validated interactions only
  summary = TRUE
)

# Processing the miRNA-mRNA interaction table
df_mirna_deg <- res_mirna_deg@data %>%
  select(mature_mirna_id, target_symbol, database, experiment) %>%
  distinct()

# 4. Phase 2: lncRNA -> miRNA Interactions ---------------------------------
# Identifying lncRNAs that interact with the validated miRNAs

unique_mirnas <- unique(df_mirna_deg$mature_mirna_id)

res_mirna_lnc <- get_multimir(
  org     = "hsa",
  mirna   = unique_mirnas,
  table   = "mirtarbase",
  summary = TRUE
)

# Filtering for lncRNAs based on official nomenclature patterns
# Patterns include: LINC, Antisense (-AS), Host Genes (HG), and known lncRNAs
df_lnc_mirna <- res_mirna_lnc@data %>%
  filter(str_detect(target_symbol, 
                    "^LINC|^MIR\\d+HG|^LOC|OT$|-AS\\d+$|GAS5|MALAT1|NEAT1|XIST|HOTAIR")) %>%
  select(target_symbol, mature_mirna_id, database, experiment) %>%
  rename(lncRNA_symbol = target_symbol) %>%
  distinct()

# 5. Phase 3: ceRNA Network Integration ------------------------------------
# Merging lncRNA-miRNA and miRNA-mRNA data to establish the regulatory axes

ceRNA_network <- inner_join(df_lnc_mirna, df_mirna_deg, 
                            by = "mature_mirna_id", 
                            relationship = "many-to-many")

# 6. Exporting Results -----------------------------------------------------

write.csv(df_mirna_deg, "1_Validated_miRNA_mRNA.csv", row.names = FALSE)
write.csv(df_lnc_mirna, "2_Validated_lncRNA_miRNA.csv", row.names = FALSE)
write.csv(ceRNA_network, "3_Complete_ceRNA_Network.csv", row.names = FALSE)

message("Analysis complete. Files saved to directory.")