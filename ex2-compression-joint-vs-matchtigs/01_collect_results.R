#!/usr/bin/env Rscript

library(tidyverse)

get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)

    if (length(file_arg) == 0) {
        return(normalizePath(".", winslash = "/", mustWork = TRUE))
    }

    normalizePath(dirname(sub("^--file=", "", file_arg[1])),
                  winslash = "/",
                  mustWork = TRUE)
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."),
                           winslash = "/",
                           mustWork = TRUE)
results_dir <- file.path(repo_root, "002-compression-joint-vs-matchtigs", "results")
output_path <- file.path(script_dir, "final_results.tsv")

result_files <- list.files(results_dir,
                           pattern = "\\.tsv$",
                           full.names = TRUE,
                           recursive = TRUE)

if (length(result_files) == 0) {
    stop("No experiment result TSV files found in ", results_dir)
}

raw_results <- map_dfr(result_files, function(path) {
    rel_path <- sub(paste0("^", stringr::fixed(results_dir), "/"), "", path)
    dataset <- dirname(rel_path)
    k_value <- as.integer(sub("\\.tsv$", "", basename(path)))

    lines <- readr::read_lines(path)
    header <- str_split(lines[1], "\\t", simplify = TRUE)
    body <- str_split(lines[-1], "[[:space:]]+")
    parsed <- do.call(rbind, body) %>%
        as_tibble(.name_repair = "minimal")
    colnames(parsed) <- header

    parsed %>%
        mutate(across(-method, as.double)) %>%
        mutate(dataset = dataset, k = k_value, .before = 1)
})

kmer_counts <- raw_results %>%
    filter(method == "kmer_count") %>%
    transmute(dataset, k, kmer_count = superstring_std_basic)

final_results <- raw_results %>%
    filter(method != "kmer_count") %>%
    pivot_longer(cols = starts_with("superstring_") | starts_with("mask_"),
                 names_to = "metric",
                 values_to = "size_bytes") %>%
    tidyr::extract(metric,
                   into = c("component", "compressor_family", "compression_setting"),
                   regex = "^(superstring|mask)_(std|geco)_(basic|max)$",
                   remove = TRUE) %>%
    left_join(kmer_counts, by = c("dataset", "k")) %>%
    mutate(
        compressor = recode(compressor_family, std = "xz", geco = "geco3"),
        bits_per_kmer = size_bytes / kmer_count * 8,
        dataset_label = str_replace_all(dataset, "_", " ")
    ) %>%
    select(dataset,
           dataset_label,
           k,
           method,
           kmer_count,
           component,
           compressor,
           compression_setting,
           size_bytes,
           bits_per_kmer) %>%
    arrange(dataset,
            k,
            method,
            compressor,
            compression_setting,
            component)

write_tsv(final_results, output_path)
