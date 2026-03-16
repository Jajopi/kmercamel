#!/usr/bin/env Rscript

library(tidyverse)

df0 <- list.files("input",
                  "\\.tsv$",
                  full.names = TRUE,
                  recursive = TRUE) |>
    map_dfr( ~ {
        data <- read_tsv(.x) #col_names = c("alg", "len", "runs"), 
        rel_path <- sub("^input/", "", .x)
        subdir <- dirname(rel_path)
        filename <- basename(.x)
        filename_no_ext <- sub("\\.tsv$", "", filename)
        data %>%
            mutate(genome = subdir, k = filename_no_ext)
    })

# List all .tsv files in the input directory
tsv_files <- list.files(
    "input",
    pattern = "\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
)
print(tsv_files)


df1 <- df0 %>%
    mutate(alg = method, len = superstring, runs = mask) %>%
    filter(k != 21, genome != "N._gonorrhoeae_pangenome") %>%
    group_by(genome, k) %>%
    mutate(kmer_count = runs[alg == "kmer_count"]) %>%
    mutate(LB_runs = if (any(alg == "lowerbound_matchtig_count"))
        runs[alg == "lowerbound_matchtig_count"]
        else
            NA) %>%
    mutate(LB_length = len[alg == "lowerbound_length"]) %>%
    # lowerbound_length = length[alg == "cyclecover", genome == genome, k == k], lowerbound_matchtig_count = runs[alg == "cyclecover", genome == genome, k == k]) %>%
    arrange(genome, k) %>%
    filter(alg != "kmer_count" &
               alg != "lowerbound_length" &
               alg != "lowerbound_matchtig_count") %>%
    # mutate(
    #     S_alg = recode(
    #         alg,
    #         "greedy-maxone" = "greedy",
    #         "greedy-minone" = "greedy",
    #         "greedy-minrun" = "greedy",
    #         "matchtigs" = "matchtigs",
    #         "simplitigs" = "simplitigs",
    #         .default = "joint"
    #     )
    # ) %>%
    mutate(run_penalty = case_when(grepl("^joint-", alg) ~ as.numeric(sub("joint-", "", alg)), TRUE ~ NA_real_))

show(df1)

df2 <- df1 %>%
    mutate(
        genome = recode(
            genome,
            "E._coli_pangenome" = "E.coli pg.",
            "N._gonorrhoeae_pangenome" = "N.gono. pg.",
            "S._Pneumoniae_pangenome" = "S.pneumo. pg.",
            "SARS-CoV-2_pangenome" = "SC2 pg."
        )
    ) %>%
    mutate(alg = recode(
        alg,
        "greedy-maxone" = "greedyMS-maxone",
        "greedy-minone" = "greedyMS-minone",
        "greedy-minrun" = "greedyMS-minrun",
        "matchtigs" = "greedy matchtigs",
        "simplitigs" = "greedy simplitigs",
        "eulertigs" = "optimal simplitigs",
         .default = alg
    ))

df <- df2


df %>%
    write_tsv("final_data.tsv")


# df1 <- df0 %>%
#     mutate(
#         genome = recode(
#             genome,
#             "sars-cov-2_pangenome_gisaid.unitigs_k128" = "sc2-pg",
#             "spneumo_pangenome_RASE_db.unitigs_k128" = "sp-pg",
#             "ngono_pangenome_RASE_db.unitigs_k128" = "ng-pg",
#             "ecoli_pangenome_661k_all.unitigs_k128" = "ec-pg-all",
#             "ecoli_pangenome_661k_HQ.unitigs_k128" = "ec-pg-hq",
#             "human_genome_assembly_chm13.v2.0" = "hg-t2t",
#             "human_genome_illumina.unitigs_minfreq2_k32" = "hg-ilm",
#             "human_rnaseq_srx348811.unitigs_minfreq2_k32" = "rna-ilm",
#             "human_microbiome_illumina_srs063932.unitigs_minfreq2_k32" = "mtg-ilm",
#             "hprc_human_assemblies_HQ.unitigs_k64" = "hprc",
#             "minikraken4GB_k31" = "minikr4gib",
#             "minikraken8GB_k31" = "minikr8gib"
#         )
#     )

# df2 <- df1 %>%
#     select(-r, -o, -z, -pref) %>%
#     rename(dataset = genome) %>%
#     select(dataset, k, S_alg, kmer_count, everything())

# df3 <-
#     bind_rows(df2,
#               df2 %>% transmute(dataset, S_alg = "cyclecover", kmer_count=kmer_count, k, l = CC_LB)) %>%
#     select(-CC_LB) %>%
#     arrange(across(everything())) %>%
#     unique()
