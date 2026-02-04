# ~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                            #
#                                                                            #
#                                                                            #
#             Targeting hyperarousal to improve sleep:                       #
#   A network intervention analysis of a digital intervention for insomnia   #
#                                                                            #
#                                                                            #
#                     Analysis/code by Dr. Linda Betz                        #
#                                                                            #
#                                                                            #
#                                                                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(patchwork)
library(mgm)
library(qgraph)
library(emmeans)
library(bootnet)

# ----------- Central Definition of Item-Labels / Symptoms -----------

# Define order and names here
item_lookup <- c(
  "ISI_1" = "dis",
  "ISI_2" = "dms",
  "ISI_3" = "ema",
  "ISI_4" = "dsat",
  "ISI_5" = "idf",
  "ISI_6" = "qol",
  "ISI_7" = "wry",
  "PHQ9_1" = "loi",
  "PHQ9_2" = "dep",
  "PHQ9_4" = "fatg",
  "PHQ9_5" = "app",
  "PHQ9_6" = "wrth",
  "PHQ9_7" = "con",
  "PHQ9_8" = "mot",
  "PHQ9_9" = "sui",
  "GAD7_1" = "nerv",
  "GAD7_2" = "wctl",
  "GAD7_3" = "wmch",
  "GAD7_4" = "relx",
  "GAD7_5" = "rstl",
  "GAD7_6" = "irr",
  "GAD7_7" = "afrd"
)

node_labels <- c("treat", unname(item_lookup))

# Generate facet levels for ggplot
facet_levels <- paste0(rep(c("ISI", "PHQ9", "GAD7"), c(7, 8, 7)), ": ", item_lookup)

# Define node colors for consistency
node_cols <- c("#e6a8ff",
               rep("#d6c1f5", 7),
               rep("grey80", 8),
               rep("#fff9f9", 7))

# ----------- Load preprocessed data from somnovia RCT -----------

data_raw <- read_csv("somnovia_RCT/data/raw/data_raw.csv")  %>%
  rename(t0_GAD7_6 = to_GAD7_6) %>%
  select(-matches("PHQ9_3")) # exclude item 3 from PHQ9 (sleep item -> topological overlap)


# ----------- Visualization: Individual Symptom Trajectories -----------

plot_data <- data_raw %>%
  select(group, matches("ISI|PHQ9|GAD7"), -matches("sum")) %>%
  pivot_longer(-group, names_to = "raw_name", values_to = "value") %>%
  mutate(
    group = factor(
      if_else(group == 1, "somnovia", "control"),
      levels = c("control", "somnovia")
    ),
    time = toupper(str_extract(raw_name, "t0|t1|t2")),
    clean_item = str_remove(raw_name, "t0_|t1_|t2_"),
    quest = str_extract(clean_item, "ISI|PHQ9|GAD7"),
    node = item_lookup[clean_item],
    facet_label = factor(paste0(quest, ": ", node), levels = facet_levels)
  ) %>%
  group_by(group, facet_label, time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

ggplot(plot_data, aes(
  x = time,
  y = value,
  color = group,
  group = group
)) +
  facet_wrap( ~ facet_label, scales = "free_y", ncol = 4) +
  geom_line(alpha = 0.8) + geom_point() +
  scale_color_manual(values = c("control" = "grey69", "somnovia" = "#BB9DD9")) +
  theme_minimal(base_size = 10) +
  labs(title = "Changes in Individual Items", y = "Score Value", x = NULL)

ggsave(
  "somnovia_NIA/plots/Supplementary_Figure_1.png",
  width = 6,
  height = 6.5
)


# ----------- Group-Effect Estimates on Individual Items (ANCOVA) -----------

# To assess item-level intervention effects, ANCOVA models were fitted
# for all individuals items at T1 and T2

outcomes <- names(item_lookup)

## ----------- T1 -----------
ancova_models_t1 <- tibble(outcome = outcomes) %>%
  mutate(model = map(outcome, ~ lm(reformulate(
    c("group", paste0("t0_", .x)), paste0("t1_", .x)
  ), data = data_raw))) %>%
  mutate(
    emm = map(model, ~ emmeans(.x, "group")),
    cohen_d_effect = map2(emm, model, ~ eff_size(
      .x, sigma = sigma(.y), edf = df.residual(.y)
    ))
  ) %>%
  mutate(cohen_d_effect = map(
    cohen_d_effect,
    ~ as_tibble(.x) %>% select(effect.size) %>% as.numeric()
  ))

ancova_models_t1 %>%
  mutate(tidy = map(model, broom::tidy),
         glance = map(model, broom::glance)) %>%
  unnest(tidy) %>%
  filter(term == "group") %>%
  transmute(outcome,
            term,
            estimate = round(estimate, 2),
            p.value,
            cohen_d_effect = cohen_d_effect) %>%
  knitr::kable()

## ----------- T2 -----------
ancova_models_t2 <- tibble(outcome = outcomes) %>%
  mutate(model = map(outcome, ~ lm(reformulate(
    c("group", paste0("t0_", .x)), paste0("t2_", .x)
  ), data = data_raw))) %>%
  mutate(
    emm = map(model, ~ emmeans(.x, "group")),
    cohen_d_effect = map2(emm, model, ~ eff_size(
      .x, sigma = sigma(.y), edf = df.residual(.y)
    ))
  ) %>%
  mutate(cohen_d_effect = map(
    cohen_d_effect,
    ~ as_tibble(.x) %>% select(effect.size) %>% as.numeric()
  ))

ancova_models_t2 %>%
  mutate(tidy = map(model, broom::tidy),
         glance = map(model, broom::glance)) %>%
  unnest(tidy) %>%
  filter(term == "group") %>%
  transmute(outcome,
            term,
            estimate = round(estimate, 2),
            p.value,
            cohen_d_effect = cohen_d_effect) %>%
  knitr::kable()


# ----------- Transform data for NIA -----------

data_nia_t0 <- data_raw %>%
  select(group, matches("t0_ISI|t0_PHQ9|t0_GAD7")) %>%
  select(-matches("sum"))

nrow(data_nia_t0) # 290

data_nia_t1 <- data_raw %>%
  select(group, matches("t1_ISI|t1_PHQ9|t1_GAD7")) %>%
  select(-matches("sum")) %>%
  drop_na() # NIA doesn't allow missings

nrow(data_nia_t1) # 248

data_nia_t2 <- data_raw %>%
  select(group, matches("t2_ISI|t2_PHQ9|t2_GAD7")) %>%
  select(-matches("sum")) %>%
  drop_na() # NIA doesn't allow missings

nrow(data_nia_t2) # 244


# Remove Baseline Differences

# To focus on treatment effects, we regress T1/T2 scores on T0 scores
# and extract the residuals. These residuals represent follow-up scores
# adjusted for baseline severity.

data_nia_adjusted_t1 <- tibble(outcome = outcomes) %>%
  mutate(model = map(outcome, ~ lm(
    reformulate(c(paste0("t0_", .x)), paste0("t1_", .x)), data = data_raw %>% filter(!is.na(t1_GAD7_sum))
  ))) %>%
  mutate(residuals = map(model, ~ residuals(.x))) %>%
  select(-model) %>%
  pivot_wider(names_from = outcome,
              names_prefix = "t1_",
              values_from = residuals) %>%
  unnest(cols = everything()) %>%
  bind_cols(data_nia_t1 %>% select(group)) %>%
  select(group, everything())


data_nia_adjusted_t2 <- tibble(outcome = outcomes) %>%
  mutate(model = map(outcome, ~ lm(
    reformulate(c(paste0("t0_", .x)), paste0("t2_", .x)), data = data_raw %>% filter(!is.na(t2_GAD7_sum))
  ))) %>%
  mutate(residuals = map(model, ~ residuals(.x))) %>%
  select(-model) %>%
  pivot_wider(names_from = outcome,
              names_prefix = "t2_",
              values_from = residuals) %>%
  unnest(cols = everything()) %>%
  bind_cols(data_nia_t2 %>% select(group)) %>%
  select(group, everything())


# ----------- Network Intervention Analysis (MGM) -----------

fit_mgm_network <- function(data) {
  set.seed(1)
  estimateNetwork(
    as.matrix(data),
    default = "mgm",
    type = c("c", rep("g", 22)),
    # treatment categorical, symptoms continuous
    level = c(2, rep(1, 22)),
    criterion = "CV",
    rule = "AND",
    nFolds = 10,
    order = 2,
    binarySign = TRUE
  )
}

mgm_t0 <- fit_mgm_network(data_nia_t0)
mgm_t1 <- fit_mgm_network(data_nia_adjusted_t1)
mgm_t2 <- fit_mgm_network(data_nia_adjusted_t2)


# ----------- Visualization -----------

# Calculate average layout based on all three time points
tmp_graphs <- list(
  qgraph(mgm_t0$graph, DoNotPlot = TRUE),
  qgraph(mgm_t1$graph, DoNotPlot = TRUE),
  qgraph(mgm_t2$graph, DoNotPlot = TRUE)
)
average_layout <- averageLayout(tmp_graphs[[1]], tmp_graphs[[2]], tmp_graphs[[3]], repulsion = 1.07)

pdf(width = 8,
    height = 3.75,
    file = "somnovia_NIA/plots/Figure_1.pdf")
layout(matrix(c(1, 2, 3), 1, 3))

# Common plot settings
q_args <- list(
  color = node_cols,
  labels = node_labels,
  theme = "colorblind",
  layout = average_layout,
  shape = c("square", rep("circle", 22)),
  label.cex = 1.25,
  cut = 0,
  minimum = 0.01,
  edge.width = 0.5,
  negDashed = TRUE
)

do.call(qgraph, c(
  list(
    input = mgm_t0$graph,
    title = "a) Baseline",
    maximum = 0.2
  ),
  q_args
))
do.call(qgraph, c(
  list(
    input = mgm_t1$graph,
    title = "b) 3 months",
    maximum = 0.3
  ),
  q_args
))
do.call(qgraph, c(
  list(
    input = mgm_t2$graph,
    title = "c) 6 months",
    maximum = 0.2
  ),
  q_args
))

dev.off()


# ----------- Edge Accuracy/Stability (Bootstrapping) -----------

# The following function is an adaptation of the "plot.bootnet" function in bootnet.
# It isolates edges for a target node (here: "treat") to ensure readability and visibility,
# especially in graphs with many nodes.
#
# Edges are ordered by their proportion of non-zero appearances (-prop0).

plot_boot_edges <- function(boot_obj,
                            target_node = "treat",
                            plot_title = NULL) {
  sumTable <- summary(boot_obj, statistics = "edge") %>%
    mutate(
      node1 = str_replace_all(node1, "group", "treat"),
      node2 = str_replace_all(node2, "group", "treat")
    ) %>%
    filter(node1 == target_node | node2 == target_node) %>%
    mutate(
      # Use the global item_lookup to clean edge labels
      id_clean = str_replace_all(id, item_lookup),
      id_clean = str_replace_all(id_clean, "group", "treat"),
      id_clean = str_remove(id_clean, "t0_|t1_|t2_"),
      alpha_val = 0.25 + (1 - 0.25) * (1 - prop0)
    )
  
  ggplot(sumTable, aes(y = fct_reorder(id_clean, -prop0))) +
    geom_segment(
      aes(
        x = q2.5_non0,
        xend = q97.5_non0,
        yend = fct_reorder(id_clean, -prop0),
        alpha = alpha_val
      ),
      colour = "black",
      linewidth = 1
    ) +
    geom_point(aes(
      x = mean_non0,
      alpha = alpha_val,
      colour = "mean"
    ), size = 3) +
    geom_point(aes(x = sample, colour = "sample"), size = 3) +
    geom_label(aes(x = 0, label = format(round(prop0, 2), nsmall = 2)), size = 2.5, alpha = 0.8) +
    scale_color_manual(name = "",
                       values = c("mean" = "black", "sample" = "red")) +
    scale_alpha_identity() +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9),
      strip.text = element_blank()
    ) +
    labs(x = "\n Estimated Edge Weight", y = "Edges") +
    ggtitle(plot_title)
}


# Run Bootstrapping & Plotting

set.seed(1)
mgm_t0_boot <- bootnet(mgm_t0, nBoots = 500, statistics =
                         c("edge"))
mgm_t1_boot <- bootnet(mgm_t1, nBoots = 500, statistics =
                         c("edge"))
mgm_t2_boot <- bootnet(mgm_t2, nBoots = 500, statistics =
                         c("edge"))

plot_edges_boot_t0 <- plot_boot_edges(mgm_t0_boot, target_node = "treat", plot_title = "a) Baseline")
plot_edges_boot_t1 <- plot_boot_edges(mgm_t1_boot, target_node = "treat", plot_title = "b) 3 months")
plot_edges_boot_t2 <- plot_boot_edges(mgm_t2_boot, target_node = "treat", plot_title = "b) 6 months")


(
  plot_boot_edges(mgm_t0_boot, plot_title = "a) Baseline") +
    plot_boot_edges(mgm_t1_boot, plot_title = "b) 3 months") +
    plot_boot_edges(mgm_t2_boot, plot_title = "c) 6 months")
) +
  plot_layout(axes = "collect", guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "somnovia_NIA/plots/Supplementary_Figure_2.png",
  width = 9,
  height = 5
)
