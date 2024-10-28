#' @title LRIHomology conversion
#'
#' @description This package can convert ligand and receptor data related to human or mouse development into data for different species.
#' Convert Homology Data for Ligand and Receptor Analysis
#'
#' @param source_file The file path of the source ligand-receptor data.
#' @param source_tax The Taxonomy ID of the source species (e.g., 9606 for human).
#' @param target_tax The Taxonomy ID of the target species.
#' @param output_file The file path to save the converted data.
#' @return A CSV file with the converted homology data.

#' @export LRI_homology_conversion
LRI_homology_conversion <- function(source_file, inTax, outTax, output_file) {

  # 内嵌确保所需包加载的函数
  ensure_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }

  # 确保 dplyr, tidyr, homologene 包已安装并加载
  ensure_packages(c("dplyr", "tidyr", "homologene"))

  # 读取数据
  genelist <- read.csv(source_file, header = TRUE)

  # 根据输入物种的 homologene ID 来选择处理逻辑
  if (inTax == 9606) {  # 9606 是人类的 homologene ID
    # 人类的同源性分析代码
    columns <- c("L1", "L2", "L3", "L4", "R1", "R2", "R3", "R4", "R5")

    # 创建一个空的同源列表
    homology_list <- list()

    # 对每个列名进行同源分析，并存储结果
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      # 重命名列名
      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      # 将同源数据存储到列表
      homology_list[[col_name]] <- homology_data
    }

    # 替换基因的函数
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # 执行基因替换
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

    # 生成最终的表格
    final_data <- expanded_data %>%
      unnest(cols = everything()) %>%
      filter(
        rowSums(!is.na(select(., L1, L2, L3, L4))) > 0 &
          rowSums(!is.na(select(., R1, R2, R3, R4, R5))) > 0
      ) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

    # 保存为CSV文件
    write.csv(final_data, output_file, row.names = FALSE)

  } else if (inTax == 10090) {  # 10090 是小鼠的 homologene ID
    # 小鼠的同源性分析代码
    genelist <- read.csv(source_file, header = TRUE)

    # 定义列名
    columns <- c("L1", "R1", "R2")

    # 创建一个空的同源列表
    homology_list <- list()

    # 对每个列名进行同源分析，并存储结果
    for (col_name in columns) {
      gene_set <- as.character(genelist[[col_name]])
      homology_data <- homologene(gene_set, inTax = inTax, outTax = outTax)

      # 重命名列名
      homology_data <- homology_data %>%
        rename(SourceGene = colnames(.)[1], TargetGene = colnames(.)[2])

      # 将同源数据存储到列表
      homology_list[[col_name]] <- homology_data
    }

    # 替换基因的函数
    replace_gene <- function(gene, homology_df) {
      if (length(gene) == 1 && gene %in% homology_df$SourceGene) {
        return(homology_df %>% filter(SourceGene == gene) %>% pull(TargetGene))
      } else {
        return(NA)
      }
    }

    # 执行基因替换
    expanded_data <- genelist %>%
      rowwise() %>%
      mutate(
        L1 = list(replace_gene(L1, homology_list[[1]])),
        R1 = list(replace_gene(R1, homology_list[[2]])),
        R2 = list(replace_gene(R2, homology_list[[3]]))
      ) %>%
      ungroup()

    # 生成最终的表格
    final_data <- expanded_data %>%
      unnest(cols = everything()) %>%
      # 保留满足条件的行：L组和R组至少各有一个非NA值
      filter(
        rowSums(!is.na(select(., L1))) > 0 &
          rowSums(!is.na(select(., R1, R2))) > 0
      ) %>%
      # 替换 NA 为空白，同时确保保持原始数据类型
      mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

    # 保存为CSV文件
    write.csv(final_data, output_file, row.names = FALSE)
  } else {
    stop("Unsupported species for homologene conversion. Only human (9606) and mouse (10090) are supported.")
  }
}
