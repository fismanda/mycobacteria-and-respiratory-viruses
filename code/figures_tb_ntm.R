# ============================================================
# Figures for TB/NTM mycobacterial paper
# Ontario, Canada 2011-2024
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)   # for combining panels

# ============================================================
# FIGURE 1 — Time series of TB and NTM weekly incidence
# ============================================================
# Requires: weekly dataset exported from SAE
# Variables needed: week_date, tb_inc_overall, ntm_inc_overall
# (annualized incidence per 100,000, already calculated in Stata)

# Load data — adjust path as needed
# weekly <- read.csv("tb_weekly_time_series_small.csv")
# weekly$week_date <- as.Date(weekly$week_date, format = "%d%b%Y")

# Pandemic shading periods
pandemic_start <- as.Date("2020-03-15")
study_end      <- as.Date("2024-03-28")

fig1 <- ggplot(weekly, aes(x = week_date)) +

  # Pandemic shading
  annotate("rect",
           xmin = pandemic_start, xmax = study_end,
           ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.5) +

  # Series lines
  geom_line(aes(y = tb_inc_overall,  color = "Tuberculosis"),
            linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = ntm_inc_overall, color = "Non-tuberculous mycobacteria"),
            linewidth = 0.5, alpha = 0.7) +

  # Optional: 4-week rolling average smoother
  # geom_smooth(aes(y = tb_inc_overall,  color = "Tuberculosis"),
  #             method = "loess", span = 0.1, se = FALSE, linewidth = 1.2) +
  # geom_smooth(aes(y = ntm_inc_overall, color = "Non-tuberculous mycobacteria"),
  #             method = "loess", span = 0.1, se = FALSE, linewidth = 1.2) +

  # Pandemic annotation
  annotate("text",
           x = as.Date("2022-01-01"), y = Inf,
           label = "Pandemic era", vjust = 1.5,
           size = 3, color = "#856404") +

  scale_color_manual(
    values = c("Tuberculosis" = "#378ADD",
               "Non-tuberculous mycobacteria" = "#C0392B")
  ) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    name   = "Annualized incidence per 100,000",
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL, color = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = c(0.25, 0.92),
    legend.background = element_blank(),
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_text(size = 9),
    axis.text        = element_text(size = 9),
    plot.margin      = margin(8, 12, 4, 8)
  )

# ============================================================
# FIGURE 2 — Forest plot: pandemic IRR by outcome and ICH status
# ============================================================
# Update these values with final model results

forest_dat <- data.frame(
  outcome = c("TB","TB","TB","NTM","NTM","NTM"),
  stratum = c("Overall","Non-ICH","ICH",
              "Overall","Non-ICH","ICH"),
  irr     = c(1.206, 1.204, 1.194,
              0.677, 0.683, 0.633),
  lo      = c(1.077, 1.071, 0.752,
              0.615, 0.616, 0.487),
  hi      = c(1.351, 1.353, 1.896,
              0.745, 0.757, 0.821)
)

forest_dat <- forest_dat %>%
  mutate(
    ypos  = c(6, 5, 4, 2, 1, 0),
    color = c("#378ADD","#1D9E75","#D85A30",
              "#378ADD","#1D9E75","#D85A30"),
    sig   = ifelse(lo > 1 | hi < 1, "solid", "open"),
    label = factor(
      paste0(outcome, " \u2014 ", stratum),
      levels = rev(c("TB \u2014 Overall","TB \u2014 Non-ICH","TB \u2014 ICH",
                     "NTM \u2014 Overall","NTM \u2014 Non-ICH","NTM \u2014 ICH"))
    )
  )

# IRR label for right margin
forest_dat$irr_label <- sprintf("%.3f (%.3f\u2013%.3f)",
                                forest_dat$irr,
                                forest_dat$lo,
                                forest_dat$hi)

fig2 <- ggplot(forest_dat, aes(x = irr, y = ypos)) +

  # null line
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +

  # divider between TB and NTM sections
  geom_hline(yintercept = 3, linetype = "dotted",
             color = "grey70", linewidth = 0.5) +

  # CI lines
  geom_segment(aes(x = lo, xend = hi, y = ypos, yend = ypos,
                   color = color), linewidth = 1.0) +

  # CI caps
  geom_segment(aes(x = lo, xend = lo,
                   y = ypos - 0.18, yend = ypos + 0.18,
                   color = color), linewidth = 1.0) +
  geom_segment(aes(x = hi, xend = hi,
                   y = ypos - 0.18, yend = ypos + 0.18,
                   color = color), linewidth = 1.0) +

  # Points — open if CI crosses 1.0
  geom_point(aes(color = color,
                 fill  = ifelse(sig == "solid", color, "white"),
                 shape = sig),
             size = 3.5, stroke = 1.2) +

  # Group labels
  annotate("text", x = 0.38, y = 5.0,
           label = "Tuberculosis",
           hjust = 0, fontface = "bold", size = 3.5, color = "grey25") +
  annotate("text", x = 0.38, y = 1.0,
           label = "Non-tuberculous\nmycobacteria",
           hjust = 0, fontface = "bold", size = 3.5,
           color = "grey25", lineheight = 0.9) +

  # IRR values on right
  geom_text(aes(x = 2.25, label = irr_label),
            hjust = 0, size = 3.0, color = "grey30") +

  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_manual(values = c("solid" = 21, "open" = 21),
                     guide = "none") +
  scale_x_log10(
    breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0),
    labels = c("0.4","0.5","0.6","0.7","0.8","1.0","1.2","1.5","2.0"),
    limits = c(0.38, 3.5)
  ) +
  scale_y_continuous(
    breaks = forest_dat$ypos,
    labels = paste0("  ", forest_dat$stratum),
    limits = c(-0.8, 7.2)
  ) +

  # Legend
  annotate("point", x = 2.85, y = 6.8, color="#378ADD",
           size=3, shape=21, fill="#378ADD") +
  annotate("point", x = 2.85, y = 6.3, color="#1D9E75",
           size=3, shape=21, fill="#1D9E75") +
  annotate("point", x = 2.85, y = 5.8, color="#D85A30",
           size=3, shape=21, fill="#D85A30") +
  annotate("text", x=2.95, y=6.8, label="Overall",  hjust=0, size=3) +
  annotate("text", x=2.95, y=6.3, label="Non-ICH",  hjust=0, size=3) +
  annotate("text", x=2.95, y=5.8, label="ICH",      hjust=0, size=3) +
  annotate("text", x=2.25, y=7.0, label="IRR (95% CI)",
           hjust=0, size=3, color="grey50", fontface="italic") +

  labs(
    x       = "Incidence rate ratio (log scale)",
    y       = NULL,
    caption = paste0("Pandemic vs pre-pandemic era (ref). ",
                     "Open circles = 95% CI crosses 1.0.\n",
                     "ICH = immunocompromising host condition.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    plot.caption       = element_text(size = 8, color = "grey55"),
    plot.margin        = margin(8, 100, 8, 8)
  ) +
  coord_cartesian(clip = "off")

# ============================================================
# Save figures
# ============================================================

ggsave("Figure1_timeseries.pdf", fig1, width = 8, height = 4)
ggsave("Figure1_timeseries.png", fig1, width = 8, height = 4, dpi = 300)

ggsave("Figure2_forest.pdf", fig2, width = 8, height = 5)
ggsave("Figure2_forest.png", fig2, width = 8, height = 5, dpi = 300)

# Combined two-panel figure (optional)
fig_combined <- fig1 / fig2 +
  plot_annotation(tag_levels = 'A')
ggsave("Figure1_2_combined.pdf", fig_combined, width = 8, height = 9)
ggsave("Figure1_2_combined.png", fig_combined, width = 8, height = 9, dpi = 300)
