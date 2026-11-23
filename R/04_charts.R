# 04_charts.R ---------------------------------------------------------------
# Makes the charts in outputs/charts from the results tables.
# RG Economics chart style: white background, Arial (or Liberation Sans, which looks the same, if Arial
# isn't installed), a bold title that states the finding, a black box round each chart and no grid lines.
# Navy (#051542) for the main line and cyan (#51d6eb) for the second. Lines are labelled at their right-hand
# end rather than with a legend.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales) })
dir.create("outputs/charts", recursive = TRUE, showWarnings = FALSE)

navy <- "#051542"; cyan <- "#51d6eb"; grey <- "#565656"; ink <- "#000000"
fam <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"

theme_rg <- function(base = 13) theme_classic(base_size = base, base_family = fam) + theme(
  plot.background  = element_rect(fill = "white", colour = NA),
  panel.background = element_rect(fill = "white", colour = NA),
  panel.border     = element_rect(colour = ink, fill = NA, linewidth = 0.7),
  panel.grid       = element_blank(),
  axis.line        = element_blank(),
  axis.ticks       = element_line(colour = ink, linewidth = 0.5),
  axis.text        = element_text(colour = ink),
  strip.background = element_blank(),
  strip.text       = element_text(face = "bold", hjust = 0, colour = ink, size = base),
  plot.title       = element_text(face = "bold", size = base * 1.6, colour = ink),
  plot.subtitle    = element_text(size = base * 1.05, colour = ink, margin = margin(b = 10)),
  plot.caption     = element_text(size = base * 0.75, colour = ink, hjust = 0, margin = margin(t = 10)),
  plot.title.position = "plot", plot.caption.position = "plot",
  legend.position = "none", panel.spacing.x = unit(1.2, "lines"))

bands <- read_csv("outputs/tables/scenario_bands.csv", show_col_types = FALSE) %>% mutate(year = as.integer(fy))
paths <- read_csv("assumptions/milburn_paths.csv", show_col_types = FALSE)
fy_lab <- function(y) paste0(substr(y, 3, 4), "/", substr(y + 1, 3, 4))
src_mc <- "Source: RG Economics modelling (not official figures) using ONS LFS/APS, ILO WP140, GO-Science AI Scenarios 2030 and Milburn Review data table 3.\nShaded bands: 10th-90th percentile of 2,000 model runs per scenario; lines: median run."
short_s <- c("S1 Slow Burn" = "Slow Burn", "S2 Open Frontier" = "Open Frontier", "S3 Augmented Growth" = "Augmented Growth",
             "S4 Transformation Economy" = "Transformation Economy", "S5 Take-Off" = "Take-Off")
mil <- paths %>% transmute(year = fy_start_year, `Milburn baseline` = baseline_neet_rate,
                           `Milburn with reforms` = reform_neet_rate) %>%
  pivot_longer(-year, names_to = "series", values_to = "v")

# ---- Chart 1: all five scenarios, with AI acting only through the jobs young people do -------------
b1 <- bands %>% filter(variant %in% c("base", "jobs_only")) %>%
  transmute(scenario = factor(short_s[scenario], levels = short_s), year,
            lo_n = 100 * no_reform_p10, md_n = 100 * no_reform_p50, hi_n = 100 * no_reform_p90,
            lo_r = 100 * package_p10,   md_r = 100 * package_p50,   hi_r = 100 * package_p90)
last1 <- factor("Take-Off", levels = short_s)
end1 <- bind_rows(
  b1 %>% filter(scenario == last1, year == 2036) %>% transmute(scenario, y = md_n, lab = "AI scenario,\nno reform", col = navy),
  b1 %>% filter(scenario == last1, year == 2036) %>% transmute(scenario, y = md_r, lab = "AI scenario,\nwith Milburn's\npolicies", col = cyan),
  tibble(scenario = last1, y = paths$baseline_neet_rate[paths$fy_start_year == 2036], lab = "Milburn\nbaseline", col = grey),
  tibble(scenario = last1, y = paths$reform_neet_rate[paths$fy_start_year == 2036], lab = "Milburn with\npolicies", col = grey))
p1 <- ggplot(b1, aes(year)) +
  geom_line(data = mil, mapping = aes(x = year, y = v, linetype = series), colour = grey, linewidth = 0.6, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo_n, ymax = hi_n), fill = navy, alpha = 0.18) +
  geom_ribbon(aes(ymin = lo_r, ymax = hi_r), fill = cyan, alpha = 0.30) +
  geom_line(aes(y = md_n), colour = navy, linewidth = 1.3) +
  geom_line(aes(y = md_r), colour = cyan, linewidth = 1.3) +
  geom_text(data = end1, aes(x = 2036.9, y = y, label = lab, colour = I(col)), hjust = 0, fontface = "bold",
            size = 3.4, lineheight = 0.9, family = fam, inherit.aes = FALSE) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = c(2026, 2031, 2036), labels = fy_lab, expand = expansion(mult = 0.04)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(5, 18), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(2026, 2036), clip = "off") +
  labs(title = "AI in young people's own jobs barely moves NEET in three scenarios, and adds about 2 points in the other two",
       subtitle = "Share of 16-24 year olds not in education, employment or training (NEET), UK. AI acts only through the jobs young people do (no extra economy-wide slump)",
       x = NULL, y = NULL, caption = src_mc) +
  theme_rg() + theme(plot.margin = margin(10, 120, 10, 10))
ggsave("outputs/charts/chart1_neet_ai_in_jobs.png", p1, width = 16, height = 6.8, dpi = 170, bg = "white")

# ---- Chart 2: Transformation Economy and Take-Off, three versions --------------------------------
vlab <- c(jobs_only = "AI in jobs only", slump_fade = "Plus slump, then\nrecovery (central)", slump_held = "Plus slump held\nto 2036 (upper\nbound)")
vcol <- c(jobs_only = cyan, slump_fade = grey, slump_held = navy)
b2 <- bands %>% filter(variant %in% names(vlab)) %>%
  transmute(scenario = factor(short_s[scenario], levels = short_s[4:5]), variant, year,
            lo = 100 * no_reform_p10, md = 100 * no_reform_p50, hi = 100 * no_reform_p90)
end2 <- b2 %>% filter(scenario == "Take-Off", year == 2036) %>%
  transmute(scenario, y = md, lab = vlab[variant], col = vcol[variant])
mil2 <- mil %>% filter(series == "Milburn baseline")
p2 <- ggplot(b2, aes(year, md, colour = variant, fill = variant)) +
  geom_line(data = mil2, mapping = aes(x = year, y = v), colour = grey, linetype = "22", linewidth = 0.6, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.3) +
  geom_text(data = end2, aes(x = 2036.9, y = y, label = lab), colour = end2$col, hjust = 0, fontface = "bold",
            size = 3.6, lineheight = 0.9, family = fam, inherit.aes = FALSE) +
  annotate("text", x = 2031, y = 12.25, label = "Milburn baseline", colour = grey, hjust = 0, size = 3.3, family = fam) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_colour_manual(values = vcol) + scale_fill_manual(values = vcol) +
  scale_x_continuous(breaks = c(2026, 2030, 2033, 2036), labels = fy_lab, expand = expansion(mult = 0.03)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(xlim = c(2026, 2036), clip = "off") +
  labs(title = "In the two disruptive scenarios, how deep the slump is and how long it lasts matter most",
       subtitle = "Youth NEET rate with no reform, UK. Slump: an extra fall in youth job openings from 2029, sized so youth NEET rises as much as in 2008-11\n(Take-Off: as much as up to the 2011 peak). Take-Off slump: 131 of 2,000 runs dropped because they needed openings to fall by more than 85%",
       x = NULL, y = NULL, caption = src_mc) +
  theme_rg() + theme(plot.margin = margin(10, 140, 10, 10))
ggsave("outputs/charts/chart2_disruptive_three_versions.png", p2, width = 15, height = 6.4, dpi = 170, bg = "white")

# ---- Chart 3: AI-proofing test, share of each policy's no-AI effect it keeps ------------------------------------------------------------
rob <- read_csv("outputs/tables/policy_robustness.csv", show_col_types = FALSE)
lev_order <- c("Full package", "Activation, health and welfare (supply)", "Subsidised youth hires (demand)",
               "Education and training places")
short_w <- sub("Transformation Economy", "Transformation\nEconomy", short_s); names(short_w) <- names(short_s)
col_lab <- function(s, v) paste0(short_w[s], ifelse(v == "base", "", paste0("\n", c(jobs_only = "AI in jobs only",
                                  slump_fade = "+ slump, recovery", slump_held = "+ slump held")[v])))
rob <- rob %>% mutate(col = col_lab(scenario, variant),
                      ord = match(scenario, names(short_s)) * 10 + match(variant, c("base", "jobs_only", "slump_fade", "slump_held")),
                      col = reorder(factor(col), ord),
                      lever = factor(lever, levels = rev(lev_order)),
                      lab = paste0(round(100 * p50), "%\n(", round(100 * p10), "-", round(100 * p90), ")"))
p3 <- ggplot(rob, aes(col, lever, fill = p50)) +
  geom_tile(colour = "white", linewidth = 2) +
  geom_text(aes(label = lab, colour = p50 < 0.72), size = 3.6, lineheight = 0.9, family = fam) +
  scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white"), guide = "none") +
  scale_fill_gradient2(low = navy, mid = "white", high = cyan, midpoint = 1, limits = c(0.5, 1.2), oob = squish,
                       labels = percent, name = "Share of no-AI\neffect kept\n(median run)") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Milburn's policies keep most of their punch, unless the slump is held in place all the way to 2036",
       subtitle = "Fall in NEET each policy delivers in 2036/37 as a share of what it delivers with no AI change. Pass mark: 80%. Brackets: 10th-90th percentile of runs",
       x = NULL, y = NULL,
       caption = "Source: RG Economics modelling (not official figures). Policies are modelled as fixed numbers of places and hires, sized to match Milburn's reform path with no AI change.") +
  theme_rg(12) + theme(axis.ticks = element_blank(), legend.position = "right", legend.title = element_text(size = 10),
                       legend.key.height = unit(1.6, "cm"))
ggsave("outputs/charts/chart3_policy_robustness.png", p3, width = 17, height = 6.6, dpi = 170, bg = "white")

# ---- Chart 5: what drives the result (one assumption changed at a time) ---------------------------
tor <- read_csv("outputs/tables/sensitivity_oneway.csv", show_col_types = FALSE)
nice <- c(A_2028 = "AI use by 2028", A_2030 = "AI use by 2030", A_2036 = "AI use by 2036",
          m_disp = "How hard AI hits exposed jobs", b_broaden = "AI spreading to more tasks",
          g_comp = "New tasks in jobs AI helps with", r_realloc = "Lost hiring that reappears elsewhere",
          tau_realloc = "How long that takes", D_demand = "Size of the economy-wide effect", D_temp = "2028 market correction",
          lag_sep = "Delay before job losses rise", G_anchor = "Today's AI effect on job openings",
          sigma_sep = "Job losses among existing workers", alpha = "How much hiring depends on openings")
lc1 <- function(x) ifelse(grepl("^AI", x), x, paste0(tolower(substr(x, 1, 1)), substring(x, 2)))
for (s in c("S4 Transformation Economy", "S5 Take-Off")) for (vr in c("slump_fade", "jobs_only")) {
  d <- tor %>% filter(scenario == s, variant == vr) %>% mutate(delta = 100 * (neet36 - central), label = nice[param]) %>%
    group_by(label) %>% mutate(span = max(abs(delta))) %>% ungroup() %>%
    filter(span > 0.005) %>% mutate(label = reorder(label, span))
  topp <- d %>% distinct(label, span) %>% arrange(desc(span)) %>% slice(1:2) %>% pull(label) %>% as.character()
  p4 <- ggplot(d, aes(delta, label, fill = side)) +
    geom_col(width = 0.6, position = "identity") +
    geom_vline(xintercept = 0, colour = ink, linewidth = 0.5) +
    scale_fill_manual(values = c(high = navy, low = cyan),
                      labels = c(high = "Assumption at high end of range", low = "Assumption at low end of range")) +
    scale_x_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), x, " pts")) +
    labs(title = paste0(short_s[s], c(slump_fade = " (slump, then recovery)", jobs_only = " (AI in jobs only)")[vr],
                        ": what drives the result"),
         subtitle = sprintf("Biggest two: %s; %s.\nChange in the 2036/37 youth NEET rate (no reform) when one assumption moves to the end of its range. Central = %.1f%%",
                            lc1(topp[1]), lc1(topp[2]), 100 * d$central[1]),
         x = NULL, y = NULL, caption = paste0("Source: RG Economics modelling (not official figures). One assumption at a time, others at central values", if (vr == "slump_fade") "; the slump is re-solved each time." else ".")) +
    theme_rg() + theme(axis.ticks.y = element_blank(), legend.position = "top", legend.justification = "left",
                       legend.title = element_blank())
  ggsave(sprintf("outputs/charts/chart5_sensitivity_%s_%s.png", substr(s, 1, 2), vr), p4, width = 14, height = 6.4, dpi = 170, bg = "white")
}

# ---- Chart 4: what the UK data shows so far, APS 2021-2025 ----------------------------------------
aps <- read_csv("outputs/tables/uk_aps_exposure_check.csv", show_col_types = FALSE) %>%
  pivot_longer(c(idx_young, idx_older), names_to = "age", values_to = "idx") %>%
  mutate(age = recode(age, idx_young = "Aged 16-24", idx_older = "Aged 25+"),
         group = factor(group, levels = c("Low exposure", "Some exposure", "High substitution exposure")))
lab5 <- aps %>% filter(year == 2025) %>% group_by(group) %>%
  mutate(vj = if (abs(diff(idx)) < 3) ifelse(idx == max(idx), -0.4, 1.4) else 0.5) %>% ungroup()
p5 <- ggplot(aps, aes(year, idx, colour = age)) +
  annotate("rect", xmin = 2020.8, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = 100, colour = grey, linewidth = 0.4, linetype = "22") +
  geom_line(linewidth = 1.3) + geom_point(size = 2.6) +
  geom_text(data = lab5, aes(x = 2025.12, label = paste0(age, ": ", round(idx)), vjust = vj), hjust = 0,
            fontface = "bold", size = 3.6, family = fam) +
  facet_wrap(~ group) +
  scale_colour_manual(values = c(`Aged 16-24` = navy, `Aged 25+` = cyan)) +
  scale_x_continuous(breaks = 2021:2025, limits = c(2020.8, 2026.4)) +
  labs(title = "Whether young workers in AI-exposed jobs are falling behind depends on the start year",
       subtitle = "UK employment by how exposed the job is to AI, 2022 = 100. Shaded: 2021, still affected by Covid. Descriptive only, survey data",
       x = NULL, y = NULL,
       caption = "Source: ONS ad hoc 3537 (APS, UK, Jan-Dec 2021-2025); ILO WP140; ONS SOC 2020 coding index. Occupation groups classed by RG Economics' mapping of ILO exposure grades.\nSuppressed cells (too few survey responses) counted as zero.") +
  theme_rg()
ggsave("outputs/charts/chart4_uk_exposure_check.png", p5, width = 13.5, height = 5.8, dpi = 170, bg = "white")
cat("Charts written.\n")
