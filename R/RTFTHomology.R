#' @title RTFTHomology conversion
#'
#' @description This package provides homology analysis for human and mouse
#' receptors, transcription factors, and target genes.
#'
#' @param source_file The file path of the source R-TF/TF-T data.
#' @param inTax The Taxonomy ID of the source species (e.g., 9606 for human).
#' @param outTax The Taxonomy ID of the target species.
#' @param output_file The file path to save the converted data.
#' @return A CSV file with the converted homology data.
#'
#' @export RTFT_homology_conversion
RTFT_homology_conversion <- function(source_file, inTax, outTax, output_file) {

  # Helper function to ensure required packages are installed and loaded
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

  # Read input data
  genelist <- read.csv(source_file, header = TRUE)

  # Choose processing logic based on species homologene taxonomy ID
  if (inTax == 9606) {  # 9606 is human homologene ID
    # Homology analysis code for human
    columns <- c("from", "to")

    # Create an empty list to store homology mapping tables
    homology_list <- list()

    # Perform homology analysis for each column and store results
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      # Rename the column names
      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      # Store in list
      homology_list[[col_name]] <- homology_data
    }

    # Replacement function for gene mapping
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # Perform gene replacement row by row
    expanded_data <- genelist %>%
      rowwise() %>%
      mutate(
        from = list(replace_gene(from, homology_list[["from"]])),
        to   = list(replace_gene(to,   homology_list[["to"]]))
      ) %>%
      ungroup()

    # Generate final output table
    final_data <- expanded_data %>%
      unnest(cols = everything()) %>%
      filter(
        rowSums(!is.na(select(., from))) > 0 &
          rowSums(!is.na(select(., to))) > 0
      ) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

    # Save to CSV file
    write.csv(final_data, output_file, row.names = FALSE)

  } else if (inTax == 10090) {  # 10090 is mouse homologene ID
    # Homology analysis code for mouse
    genelist <- read.csv(source_file, header = TRUE)

    # Define column names
    columns <- c("from", "to")

    # Create list for storing homology mapping tables
    homology_list <- list()

    # Perform homology analysis and store results
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      # Rename the column names
      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      homology_list[[col_name]] <- homology_data
    }

    # Replacement function for gene mapping
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # Perform gene replacement
    expanded_data <- genelist %>%
      rowwise() %>%
      mutate(
        from = list(replace_gene(from, homology_list[[1]])),
        to   = list(replace_gene(to,   homology_list[[2]]))
      ) %>%
      ungroup()

    # Generate final output table
    final_data <- expanded_data %>%
      unnest(cols = everything()) %>%
      # Keep rows where both "from" and "to" have at least one non-NA value
      filter(
        rowSums(!is.na(select(., from))) > 0 &
          rowSums(!is.na(select(., to))) > 0
      ) %>%
      # Replace NA with blank while keeping original data types
      mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

    # Save to CSV file
    write.csv(final_data, output_file, row.names = FALSE)

  } else {
    stop("Unsupported species for homologene conversion. Only human (9606) and mouse (10090) are supported.")
  }
}
