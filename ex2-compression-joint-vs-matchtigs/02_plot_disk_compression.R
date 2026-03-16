#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 12
w <- 20
u <- "cm"
fn <- "fig_Pareto_front.pdf"

ps <- 3 # point size
fs <- 13 # font size


theme_set(theme_bw(base_size = fs))

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


nejm <- pal_nejm("default")(8)

get_method_color <- function(method) {
    dplyr::case_when(
        str_detect(method, "unitigs|simplitigs|eulertigs") ~ nejm[2],
        str_detect(method, "matchtigs") ~ nejm[2],
        str_detect(method, "greedy") ~ nejm[4],
        str_detect(method, "joint") ~ nejm[1],
        TRUE ~ NA_character_
    )
}

sorting_key_major <- function(method) {
    dplyr::case_when(
        str_detect(method, "unitigs") ~ -1,
        str_detect(method, "simplitigs|eulertigs") ~ 0,
        str_detect(method, "matchtigs") ~ 1,
        str_detect(method, "greedy") ~ 2,
        str_detect(method, "joint") ~ 4,
        TRUE ~ Inf
    )
}

sorting_key_minor <- function(method) {
    out <- rep(0, length(method))
    joint_idx <- str_detect(method, "^joint-")
    out[joint_idx] <- as.numeric(sub("^joint-", "", method[joint_idx]))
    out
}

script_dir <- get_script_dir()
results_path <- file.path(script_dir, "final_results.tsv")
plot_dir <- file.path(script_dir, "plots")

if (!file.exists(results_path)) {
    stop("Missing input file: ", results_path,
         ". Run 01_collect_results.R first.")
}

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

all_results <- read_tsv(results_path, show_col_types = FALSE) %>%
    mutate(
        dataset_label = factor(dataset_label, levels = c("S. Pneumoniae pangenome", "SARS-CoV-2 pangenome", "E. coli pangenome")),
        method_color = get_method_color(method),
        method_major = sorting_key_major(method),
        method_minor = sorting_key_minor(method),
        component = factor(component, levels = c("superstring", "mask")),
        compressor = factor(compressor, levels = c("xz", "geco3"))
    ) %>%
    filter(dataset != "N._gonorrhoeae_pangenome")

make_plot <- function(k_value, compression_setting_name) {
    selected <- all_results %>%
        filter(k == k_value,
               compression_setting == compression_setting_name) %>%
        arrange(dataset, method_major, method_minor, compressor, component)

    if (nrow(selected) == 0) {
        return(invisible(NULL))
    }

    method_levels <- selected %>%
        distinct(method, method_major, method_minor) %>%
        arrange(method_major, method_minor) %>%
        pull(method)

    method_positions <- tibble(method = method_levels,
                               method_index = seq_along(method_levels)) %>%
        mutate(method_label = recode(
            method,
            "matchtigs" = "greedy matchtigs",
            "eulertigs" = "optimal simplitigs",
            "greedy-minrun" = "greedy MS (minrun)",
            "greedy-minone" = "greedy MS (minone)",
            "greedy-maxone" = "greedy MS (maxone)",
            .default = method
        )) %>%
        mutate(method_label = if_else(
            str_detect(method, "^joint-"),
            paste0("Pareto-optimized MS (", sub("^joint-", "", method), ")"),
            method_label
        ))

    offsets_alphas <- tibble(
        compressor = factor(c("xz", "geco3"), levels = c("xz", "geco3")),
        offset = c(-0.225, 0.225),
        alpha = c(0.5, 1.0)
    )

    bar_data <- selected %>%
        left_join(method_positions, by = "method") %>%
        left_join(offsets_alphas, by = "compressor") %>%
        group_by(dataset, dataset_label, method, compressor) %>%
        arrange(component, .by_group = TRUE) %>%
        mutate(
            ymin = lag(cumsum(bits_per_kmer), default = 0),
            ymax = cumsum(bits_per_kmer),
            xmin = method_index + offset - 0.2,
            xmax = method_index + offset + 0.2,
            fill_color = if_else(component == "mask", "black", method_color)
        ) %>%
        ungroup()

    line_data <- bind_rows(
        selected %>%
            filter(method %in% c("matchtigs", "greedy-minrun")) %>%
            group_by(dataset, dataset_label, compressor, method) %>%
            summarise(total = sum(bits_per_kmer),
                      color = first(method_color),
                      .groups = "drop") %>%
            group_by(dataset, dataset_label, compressor) %>%
            slice_min(total, n = 1, with_ties = FALSE) %>%
            ungroup(),
        selected %>%
            filter(str_detect(method, "joint")) %>%
            group_by(dataset, dataset_label, compressor, method) %>%
            summarise(total = sum(bits_per_kmer),
                      color = first(method_color),
                      .groups = "drop") %>%
            group_by(dataset, dataset_label, compressor) %>%
            slice_min(total, n = 1, with_ties = FALSE) %>%
            ungroup()
    ) %>%
        left_join(offsets_alphas %>% select(compressor, alpha), by = "compressor") %>%
        filter(!is.na(total), total >= 0, total <= 2.5)

    # show(line_data)

    plot_obj <- ggplot() +
        geom_rect(data = bar_data,
                  aes(xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax,
                      fill = fill_color,
                      alpha = alpha),
                  color = NA,
                  ) +
        geom_hline(data = line_data,
                   aes(yintercept = total, color = color, alpha = alpha),
                   linetype = "dashed",
                   linewidth = 0.6,
                #    alpha = if (compression_setting == "basic") 0.4 else 1.0
                   ) +
        scale_fill_identity() +
        scale_alpha_identity(guide = "none") +
        scale_color_identity() +
        scale_x_continuous(breaks = method_positions$method_index,
                           labels = method_positions$method_label,
                           expand = expansion(mult = c(0.04, 0.04))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
        coord_cartesian(ylim = c(0, NA)) +
        facet_wrap(~ dataset_label, ncol = 3) +
        labs(x = NULL,
             y = "bits per distinct k-mer") +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.margin = margin(1.5, 3, 1.5, 10)
        )

    out_base <- file.path(plot_dir, paste0("fig_disk_compr-k_", k_value, "-", compression_setting_name))

    ggsave(paste0(out_base, ".pdf"),
           plot = plot_obj,
           width = w,
           height = h,
           units = u)
    # ggsave(paste0(out_base, ".png"),
    #        plot = plot_obj,
    #        width = 6.8,
    #        height = 6.5,
    #        units = "in",
    #        dpi = 300)
}

plot_jobs <- all_results %>%
    distinct(k) %>%
    arrange(k)

purrr::walk(plot_jobs$k, function(k) {
    make_plot(k, "basic")
    make_plot(k, "max")
})
