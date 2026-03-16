#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 12
w <- 20
u <- "cm"
fn <- "fig_Pareto_front.pdf"

ps <- 3 # point size
fs <- 13 # font size


# Algorithms to keep (and their display order)
# alg_levels <- c("global", "greedy-matchtigs", "simplitigs", "cyclecover")

my_shapes <- c(
    "Pareto-optimized MS" = 18,
    "greedyMS-maxone" = 6,
    "greedyMS-minone" = 24,
    "greedyMS-minrun" = 5,
    "greedy matchtigs" = 20,
    "optimal simplitigs" = 17
)

nejm <- pal_nejm("default")(8)

my_colors <- c(
    "Pareto-optimized MS" = nejm[1],
    "greedyMS-maxone" = nejm[4],
    "greedyMS-minone" = nejm[4],
    "greedyMS-minrun" = nejm[4],
    "greedy matchtigs" = nejm[2],
    "optimal simplitigs" = nejm[2]
)

#--------------------------

df0 <- read_tsv("final_data.tsv")

df1 <- df0 %>%
    filter(k != 21, genome != "ng-pang", alg != "unitigs", alg != "greedy simplitigs") %>%
    mutate(genome = str_replace(genome, "SC2 pg.", "SARS-CoV2 pg.")) %>%
    mutate(genome = factor(genome, levels = c("S.pneumo. pg.", "SARS-CoV2 pg.", "E.coli pg.")),
           type = if_else(str_starts(alg, "joint-"), "Pareto-optimized MS", alg))

df <- df1


all_types <- unique(df$type)
greedy_types <- sort(all_types[str_detect(all_types, "^greedyMS-")])
type_levels <- c("Pareto-optimized MS", greedy_types, "greedy matchtigs", "optimal simplitigs")
type_levels <- type_levels[type_levels %in% all_types]


#--------------------------

theme_set(theme_bw(base_size = fs))

ggplot(df,
       aes(
           x = runs / kmer_count,
           y = len / kmer_count,
           shape = type,
           color = type
       )) +
    facet_grid(genome ~ k,
               scales = "free_x",
               labeller = labeller(
                   k = function(x)
                       paste0("k = ", x)
               )) +
    geom_point(
        data = df %>% filter(type == "greedy matchtigs"),
        size = ps,
        stroke = 1.2
    ) +
    geom_point(
        data = df %>% filter(type != "greedy matchtigs"),
        size = ps,
        stroke = 1.2
    ) +
    # geom_line(
    #     data = df %>% filter(type == "joint"),
    #     aes(group = interaction(genome, k)),
    #     linetype = "solid",
    #     color = nejm[1],
    #     size = 0.3,
    #     show.legend = FALSE
    # ) +
    geom_vline(
        data = df %>% group_by(genome, k) %>% slice(1) %>% ungroup() %>% filter(!is.na(LB_runs)),
        aes(xintercept = LB_runs / kmer_count),
        linetype = "dashed",
        color = "black"
    ) +
    geom_hline(
        data = df %>% group_by(genome, k) %>% slice(1) %>% ungroup() %>% filter(!is.na(LB_length)),
        aes(yintercept = LB_length / kmer_count),
        linetype = "dashed",
        color = "black"
    ) +
    scale_x_continuous(
        "runs per dist. k-mer",
        limits = c(0, NA),
        breaks = scales::breaks_pretty(n = 4),
        labels = function(x) ifelse(x == 0, "0", sprintf("%g", x)),
        expand = expansion(mult = c(0.005, 0.015))
    ) +
    scale_y_continuous("chars per dist. k-mer",
                       limits = c(1, 2.1),
                       expand = c(0, 0)) +
    scale_shape_manual(values = my_shapes,
                       min = 1,
                       breaks = type_levels) +
    scale_color_manual(values = my_colors, breaks = type_levels) +
    labs(shape = NULL, color = NULL) +
    theme(
        legend.position = "top",
        legend.spacing.y = unit(0, "lines"),
        legend.key.height = unit(0.6, "lines"),
        legend.box.spacing = unit(0, "lines"),
        legend.margin = margin(0, 0, 5, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(color = "gray60", fill = NA),
        strip.text = element_text(margin = margin(3, 3, 3, 3)),
            plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    )

ggsave(fn,
       height = h,
       width = w,
       unit = u)

ggsave(
    sub("\\.pdf$", ".png", fn),
    height = h,
    width = w,
    units = u,
    dpi = 300
)
