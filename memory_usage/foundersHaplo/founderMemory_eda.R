## founderMemory_eda.R
## Exploratory visualisation of MapPop memory benchmark results
## Run AFTER founderMemory_benchmark.R has produced results/memory_results.rds

library(tidyverse)

# ── Load ──────────────────────────────────────────────────────────────────────
df <- readRDS("memory_usage/foundersHaplo/results/memory_results.rds")


# Pretty labels for facets
df <- df %>%
  mutate(
    segSites_label = paste0("segSites = ", segSites),
    nChr_label     = paste0("nChr = ", nChr),
    method         = factor(method, levels = c("quick", "macs"))
  )

# ── 1. Overview table ─────────────────────────────────────────────────────────
cat("\n── Memory summary (MB) ──\n")
df %>%
  arrange(method, nChr, segSites, nInd) %>%
  print(n = Inf)

# ── 2. Line plot: mem_mb vs nInd, faceted by segSites × nChr ─────────────────
p1 <- ggplot(df, aes(x = nInd, y = mem_mb, colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_continuous(labels = scales::label_number(suffix = " MB")) +
  facet_grid(nChr_label ~ segSites_label) +
  labs(
    title    = "MapPop memory usage by method and population size",
    subtitle = "Faceted by number of chromosomes and segregating sites",
    x        = "Number of individuals (log scale)",
    y        = "Memory (MB)",
    colour   = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("memory_usage/foundersHaplo/results/memory_results.png", p1, width = 10, height = 7, dpi = 150)

message("Saved: results/plot_mem_vs_nInd.png")

# ── 3. Bar chart: absolute memory for every scenario ─────────────────────────
p2 <- df %>%
  mutate(scenario = paste0("nChr=", nChr, "\nseg=", segSites)) %>%
  ggplot(aes(x = factor(nInd, labels = scales::label_number(big.mark=",")(nInd_range)),
             y = mem_mb, fill = method)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  scale_y_continuous(labels = scales::label_number(suffix = " MB")) +
  facet_wrap(~ scenario, ncol = 4) +
  labs(
    title  = "MapPop memory: quick vs macs",
    x      = "nInd",
    y      = "Memory (MB)",
    fill   = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

# pull nInd_range from data so script is self-contained
nInd_range <- sort(unique(df$nInd))

p2 <- df %>%
  mutate(scenario = paste0("nChr=", nChr, " | seg=", segSites),
         nInd_label = scales::label_number(big.mark = ",")(nInd)) %>%
  ggplot(aes(x = reorder(nInd_label, nInd), y = mem_mb, fill = method)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  scale_y_continuous(labels = scales::label_number(suffix = " MB")) +
  facet_wrap(~ scenario, ncol = 2) +
  labs(
    title = "MapPop memory: quick vs macs",
    x     = "nInd",
    y     = "Memory (MB)",
    fill  = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("memory_usage/foundersHaplo/results/plot_bar_comparison.png", p2, width = 10, height = 7, dpi = 150)
message("Saved: results/plot_bar_comparison.png")

# ── 4. Memory ratio quick / macs ──────────────────────────────────────────────
ratio_df <- df %>%
  select(method, nInd, nChr, segSites, mem_mb) %>%
  pivot_wider(names_from = method, values_from = mem_mb) %>%
  mutate(ratio_quick_over_macs = quick / macs)

cat("\n── quick / macs memory ratio ──\n")
print(ratio_df, n = Inf)

p3 <- ratio_df %>%
  mutate(nChr_label     = paste0("nChr = ", nChr),
         segSites_label = paste0("segSites = ", segSites)) %>%
  ggplot(aes(x = nInd, y = ratio_quick_over_macs)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#2c7fb8", linewidth = 0.8) +
  geom_point(size = 2.5, colour = "#2c7fb8") +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  facet_grid(nChr_label ~ segSites_label) +
  labs(
    title    = "Memory ratio: quick / macs",
    subtitle = "Values > 1 → quick uses more memory than macs",
    x        = "nInd (log scale)",
    y        = "Ratio (quick / macs)"
  ) +
  theme_bw(base_size = 12)

ggsave("memory_usage/foundersHaplo/results/plot_ratio.png", p3, width = 9, height = 6, dpi = 150)
message("Saved: results/plot_ratio.png")

# ── 5. Heatmap: memory landscape ─────────────────────────────────────────────
p4 <- df %>%
  ggplot(aes(x = factor(nInd), y = factor(segSites), fill = mem_mb)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mem_mb, 1)), size = 3) +
  scale_fill_viridis_c(option = "plasma", name = "MB", na.value = "grey80") +
  facet_grid(nChr_label ~ method) +
  scale_x_discrete(labels = function(x) scales::label_number(big.mark = ",")(as.numeric(x))) +
  labs(
    title = "Memory heatmap (MB)",
    x     = "nInd",
    y     = "segSites"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("memory_usage/foundersHaplo/results/plot_heatmap.png", p4, width = 10, height = 7, dpi = 150)
message("Saved: results/plot_heatmap.png")

message("\nAll EDA plots saved to results/")

