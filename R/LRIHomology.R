#' @title LRIHomology conversion
#'
#' @description This package can convert ligand and receptor data related to human or mouse development into data for different species.
#' Convert Homology Data for Ligand and Receptor Analysis
#'
#' @param source_file The file path of the source ligand-receptor data.
#' @param inTax The Taxonomy ID of the source species (e.g., 9606 for human).
#' @param outTax The Taxonomy ID of the target species.
#' @param output_file The file path to save the converted data.
#' @return A CSV file with the converted homology data.
#'
#' @export LRI_homology_conversion
LRI_homology_conversion <- function(source_file, inTax, outTax, output_file) {

  # Helper to ensure required packages are installed and loaded
  ensure_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }

  # Ensure dplyr, tidyr, homologene are installed and loaded
  ensure_packages(c("dplyr", "tidyr", "homologene"))

  # Read input ligand-receptor table
  genelist <- read.csv(source_file, header = TRUE, stringsAsFactors = FALSE)

  # Branch logic based on source species homologene taxonomy id
  if (inTax == 9606) {  # 9606 is human
    # Columns expected in human L-R tables (adjust if your table uses different column names)
    columns <- c("L1", "L2", "L3", "L4", "R1", "R2", "R3", "R4", "R5")

    # Prepare a list to store homology lookup tables for each column
    homology_list <- list()

    # For each column, run homologene mapping and store results
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      # Rename columns to standard names: SourceGene, TargetGene
      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      homology_list[[col_name]] <- homology_data
    }

    # Function to replace a single gene by mapped target gene (returns NA if not mapped)
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && !is.na(gene) && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # Apply replacements rowwise: map each ligand/receptor column using the precomputed lookup
    expanded_data <- genelist %>%
      rowwise() %>%
      mutate(
        L1 = list(replace_gene(L1, homology_list[["L1"]])),
        L2 = list(replace_gene(L2, homology_list[["L2"]])),
        L3 = list(replace_gene(L3, homology_list[["L3"]])),
        L4 = list(replace_gene(L4, homology_list[["L4"]])),
        R1 = list(replace_gene(R1, homology_list[["R1"]])),
        R2 = list(replace_gene(R2, homology_list[["R2"]])),
        R3 = list(replace_gene(R3, homology_list[["R3"]])),
        R4 = list(replace_gene(R4, homology_list[["R4"]])),
        R5 = list(replace_gene(R5, homology_list[["R5"]]))
      ) %>%
      ungroup()

    # Unnest all list-columns, keep rows where ligand side and receptor side have at least one mapped gene
    final_data <- expanded_data %>%
      tidyr::unnest(cols = dplyr::everything()) %>%
      dplyr::filter(
        rowSums(!is.na(dplyr::select(., L1, L2, L3, L4))) > 0 &
          rowSums(!is.na(dplyr::select(., R1, R2, R3, R4, R5))) > 0
      ) %>%
      dplyr::mutate(across(.cols = dplyr::everything(), ~ ifelse(is.na(.), "", .)))

    # Save result to CSV
    write.csv(final_data, output_file, row.names = FALSE)

  } else if (inTax == 10090) {  # 10090 is mouse
    # Read again to be safe (could reuse genelist defined above)
    genelist <- read.csv(source_file, header = TRUE, stringsAsFactors = FALSE)

    # Define columns expected for mouse input (adjust if different)
    columns <- c("L1", "R1", "R2")

    homology_list <- list()

    # Compute homology mapping per column
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      homology_list[[col_name]] <- homology_data
    }

    # Replacement helper (returns NA if no mapping)
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && !is.na(gene) && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # Map genes rowwise using the computed homology tables
    expanded_data <- genelist %>%
      rowwise() %>%
      mutate(
        L1 = list(replace_gene(L1, homology_list[["L1"]])),
        R1 = list(replace_gene(R1, homology_list[["R1"]])),
        R2 = list(replace_gene(R2, homology_list[["R2"]]))
      ) %>%
      ungroup()

    # Keep rows where L side has at least one mapped gene and R side has at least one mapped gene
    final_data <- expanded_data %>%
      tidyr::unnest(cols = dplyr::everything()) %>%
      dplyr::filter(
        rowSums(!is.na(dplyr::select(., L1))) > 0 &
          rowSums(!is.na(dplyr::select(., R1, R2))) > 0
      ) %>%
      dplyr::mutate(across(.cols = dplyr::everything(), ~ ifelse(is.na(.), "", .)))

    write.csv(final_data, output_file, row.names = FALSE)

  } else {
    stop("Unsupported species for homologene conversion. Only human (9606) and mouse (10090) are supported.")
  }
}
