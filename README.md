# Devcomhomologene
Identification of ortholog ligand, receptor, TF or TG in other vertebrate species
To extend DEVCOM beyond human and mouse, we implemented DevComHomology, a module that projects DEVCOMdb L-R-TF-TG networks to additional vertebrates. DEVCOMdb is curated primarily in human and mouse; given a target species and its gene identifiers, DevComHomology normalizes genes to NCBI IDs and queries NCBI HomoloGene to retrieve orthologs. For each ligand, receptor, TF and TG in DEVCOMdb, the corresponding orthologs are substituted to generate species-specific CCC networks. Homologous interactions can currently be inferred for 20 commonly used vertebrate species, according to HomoloGene statistics (https://www.ncbi.nlm.nih.gov/homologene/statistics/). This strategy enables users to reuse human or mouse-based DEVCOMdb priors when applying DEVCOM in other vertebrate model systems.

```


library(dplyr)
library(tidyr)
library(homologene)
library(Devcomhomologene)

# Suppose the function LRI_homology_conversion() is defined in your environment,
# or source the file where you stored it:
# source("path/to/your/LRI_homology_conversion.R")

# Create a small example ligand-receptor CSV for demonstration
example_df <- data.frame(
  L1 = c("FGF2", "EGF", "VEGFA"),
  L2 = c("FGF1", NA, NA),
  L3 = NA,
  L4 = NA,
  R1 = c("FGFR1", "EGFR", "KDR"),
  R2 = NA,
  R3 = NA,
  R4 = NA,
  R5 = NA,
  stringsAsFactors = FALSE
)

# Save example to disk
example_input <- "example_lr_human.csv"
write.csv(example_df, example_input, row.names = FALSE)

# Call the conversion: mapping from human (9606) to mouse (10090)
output_file <- "example_lr_human_to_mouse.csv"
LRI_homology_conversion(
  source_file = example_input,
  inTax = 9606,
  outTax = 10090,
  output_file = output_file
)

# Check the output
cat("Converted CSV saved to:", output_file, "\n")
read <- read.csv(output_file, stringsAsFactors = FALSE)
print(read)

```
