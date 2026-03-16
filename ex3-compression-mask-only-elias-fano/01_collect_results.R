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
results_dir <- file.path(repo_root, "003-compression-mask-only-elias-fano", "results")
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
    transmute(dataset, k, kmer_count = mask)

lb_len <- raw_results %>%
    filter(method == "lowerbound_length") %>%
    transmute(dataset, k, lb_length = superstring)

lb_runs <- raw_results %>%
    filter(method == "lowerbound_matchtig_count") %>%
    transmute(dataset, k, lb_runs = mask)

final_results <- raw_results %>%
    filter(method != "kmer_count", method != "lowerbound_length", method != "lowerbound_matchtig_count") %>%
    left_join(kmer_counts, by = c("dataset", "k")) %>%
    left_join(lb_len, by = c("dataset", "k")) %>%
    left_join(lb_runs, by = c("dataset", "k")) %>%
    mutate(
        S_bits_per_kmer = superstring / kmer_count,
        M_bits_per_kmer = mask / kmer_count,
        dataset_label = str_replace_all(dataset, "_", " ")
    ) %>%
    mutate(lb = lb_length * 2 / kmer_count + lb_runs * log2(lb_length / lb_runs + 1) / kmer_count) %>%
    select(dataset,
           dataset_label,
           k,
           method,
           kmer_count,
           S_bits_per_kmer,
           M_bits_per_kmer,
           lb) %>%
    arrange(dataset,
            k,
            method)

write_tsv(final_results, output_path)
