library(ggplot2)
library(dplyr)
library(patchwork)

strata <- c("Overall","Non-ICH","ICH",
            "Under-65","Over-65")

pan_dat <- data.frame(
  stratum = rep(strata, 2),
  outcome = rep(c("TB","NTM"), each=5),
  irr = c(1.206,1.204,1.194,1.416,0.950,
          0.677,0.683,0.633,0.657,0.689),
  lo  = c(1.077,1.071,0.752,1.223,0.794,
          0.615,0.616,0.487,0.563,0.610),
  hi  = c(1.351,1.353,1.896,1.639,1.136,
          0.745,0.757,0.821,0.768,0.778)) %>%
  mutate(
    stratum = factor(stratum, levels=rev(strata)),
    sig = ifelse(lo > 1 | hi < 1, "sig", "ns"),
    irr_label = sprintf("%.2f (%.2f-%.2f)", irr, lo, hi))

dlnm_dat <- data.frame(
  stratum = rep(strata, 2),
  outcome = rep(c("TB","NTM"), each=5),
  irr = c(1.110,1.113,1.039,1.150,1.044,
          1.006,1.014,0.945,0.973,1.033),
  lo  = c(1.055,1.056,0.859,1.084,0.967,
          0.953,0.960,0.839,0.905,0.974),
  hi  = c(1.168,1.173,1.257,1.221,1.126,
          1.062,1.071,1.062,1.047,1.095)) %>%
  mutate(
    stratum = factor(stratum, levels=rev(strata)),
    sig = ifelse(lo > 1 | hi < 1, "sig", "ns"),
    irr_label = sprintf("%.2f (%.2f-%.2f)", irr, lo, hi))

make_forest <- function(dat, title, xlab, xlim, xbreaks) {
  ggplot(dat, aes(x=irr, y=stratum,
                  color=outcome, shape=sig)) +
    geom_vline(xintercept=1, linetype="dashed",
               color="grey50", linewidth=0.5) +
    geom_errorbarh(aes(xmin=lo, xmax=hi),
                   height=0.25, linewidth=0.8,
                   position=position_dodge(0.6)) +
    geom_point(size=3.5, stroke=1.2,
               position=position_dodge(0.6)) +
    geom_text(aes(x=xlim[2]*0.82,
                  label=irr_label),
              position=position_dodge(0.6),
              hjust=0, size=2.8, color="grey30") +
    scale_color_manual(values=c(
      "TB"="#C0392B","NTM"="#2980B9")) +
    scale_shape_manual(values=c(
      "sig"=19,"ns"=21), guide="none") +
    scale_x_log10(breaks=xbreaks,
                  limits=xlim) +
    labs(x=xlab, y=NULL,
         title=title, color=NULL) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom",
          panel.grid.major.y=element_line(
            color="grey92"),
          panel.grid.minor=element_blank())
}

p_pan <- make_forest(pan_dat,
                     "A. Pandemic era vs pre-pandemic",
                     "Incidence rate ratio (log scale)",
                     c(0.45, 2.8),
                     c(0.5,0.7,1.0,1.2,1.5,2.0))

p_dlnm <- make_forest(dlnm_dat,
                      "B. SARS-CoV-2 cumulative effect (per SD)",
                      "Incidence rate ratio (log scale)",
                      c(0.78, 1.45),
                      c(0.8,0.9,1.0,1.1,1.2,1.3,1.4))

p_combined <- p_pan + p_dlnm +
  plot_layout(guides="collect") &
  theme(legend.position="bottom")

png("fig2_forest.png", width=10, height=5,
    units="in", res=300)
print(p_combined)
dev.off()

print(p_combined)
