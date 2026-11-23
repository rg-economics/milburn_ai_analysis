# 03_run.R ------------------------------------------------------------------
# Fits the model, runs the five GO-Science AI scenarios many times (Monte Carlo), tests whether
# Milburn's policies still work, runs the sensitivity checks and saves the results tables.
# Charts are made in 04_charts.R.
#
# Transformation Economy (S4) and Take-Off (S5) are each run three ways (see section 3):
#   jobs_only  : AI acts only through the jobs young people do
#   slump_fade : plus an extra fall in youth job openings from 2029 (the "slump"), sized in every run so
#                the youth NEET rate in 2030/31 rises as much as it did in 2008-11, then lifts by early 2037.
#                This is the central case
#   slump_held : as slump_fade, but the slump stays in place to 2036. This is the upper bound

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(purrr) })
source("R/02_model.R")
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- Read in the data and assumptions --------------------------------------------------------------------
gpd <- read_csv("assumptions/global.csv", show_col_types = FALSE)
gp  <- as.list(setNames(gpd$central, gpd$param))
scen_long <- read_csv("assumptions/scenarios.csv", show_col_types = FALSE)
scen_names <- unique(scen_long$scenario)
stocks <- read_csv("data/derived/stocks_2026q2.csv", show_col_types = FALSE)
haz    <- read_csv("data/derived/hazards_all_ages.csv", show_col_types = FALSE)
entry  <- read_csv("data/derived/entry_at_16.csv", show_col_types = FALSE)
q2324  <- read_csv("data/derived/neet_23_24.csv", show_col_types = FALSE)
occ    <- read_csv("data/derived/occupations_youth_exposure.csv", show_col_types = FALSE) %>%
  filter(!is.na(s)) %>% mutate(w_youth = w_youth / sum(w_youth))
paths  <- read_csv("assumptions/milburn_paths.csv", show_col_types = FALSE)
hist   <- read_csv("data/derived/neet_history.csv", show_col_types = FALSE)
targets <- read_csv("assumptions/scenario_targets.csv", show_col_types = FALSE)
pop    <- sum(stocks$level)
start_rate <- hist$neet_rate[hist$period == "Apr-Jun 2026"]
set.seed(gp$seed)

sp_of <- function(scn, which = "central")
  { d <- scen_long %>% filter(scenario == scn); as.list(setNames(d[[which]], d$param)) }
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
fyv <- function(run, y, col = "neet") fy_avg(run, col)[[as.character(y)]]

# ---- 1. Fit the model to today's data and to Milburn's baseline ------------------------------------------------------------
calibrate_all <- function(g) {
  cal <- calibrate_ss(stocks, haz, entry, q2324, g)
  tq  <- quarterly_target(paths, start_rate)
  b   <- calibrate_baseline(cal, tq)
  list(cal = cal, b = b, tq = tq)
}
CB <- calibrate_all(gp); cal <- CB$cal; b <- CB$b
base_run <- simulate(cal, b$phi, b$Sb, alpha = gp$alpha)
base_fy  <- fy_avg(base_run)
cat("Calibration fit:", signif(cal$fit, 3), "| fitted: ksep", round(cal$ksep, 2), "xexit", round(cal$xexit, 3),
    "kfind", round(cal$kfind, 2), "kui", round(cal$kui, 2), "\n")
write_csv(tibble(fy = paths$fy, milburn_baseline = paths$baseline_neet_rate,
                 model_no_ai = round(100 * as.numeric(base_fy), 2)), "outputs/tables/calibration_check.csv")
write_csv(tibble(param = c("ksep", "xexit", "kfind", "kui", "fit_objective"),
                 value = c(cal$ksep, cal$xexit, cal$kfind, cal$kui, cal$fit),
                 meaning = c("youth multiplier on all-age E->U and E->I hazards",
                             "quarterly exit hazard from education (all learners, not funded places)",
                             "youth multiplier on all-age U->E and I->E hazards",
                             "youth multiplier on all-age U->I hazard",
                             "sum of squared relative errors in steady-state stocks")),
          "outputs/tables/fitted_flow_parameters.csv")

# ---- 2. Milburn's policies, sized to match his reform path when there is no AI change ----
# There are three types: education and training places (educ), subsidised hires (demand) and
# activation, health and welfare support (supply). Each is a fixed number of places or hires per quarter.
# Step 1: size each policy so the 2036/37 fall in the NEET rate is split between them as assumed.
# Step 2: scale the whole package up over time so, with no AI change, it follows Milburn's reform
# path (his Table 3) in every year.
gap_fy <- paths$baseline_neet_rate - paths$reform_neet_rate
P0 <- approx(c(T0, paths$fy_start_year + 0.625), c(0, gap_fy) / max(gap_fy), xout = qtimes, rule = 2)$y
fy_q <- floor(qtimes - 0.25)
neet36 <- function(run) fyv(run, 2036)
mkpol <- function(Ke = 0, Ks = 0, Kd = 0, beta = gp$beta_takeup, P = P0, spell = gp$edp_spell_q)
  list(P = P, Ke = Ke, Ks = Ks, Kd = Kd, beta = beta, edp_spell = spell)
run_pol <- function(pol, V = rep(1, NQ), sep = rep(1, NQ), alpha) simulate(cal, b$phi, b$Sb, V, sep, alpha, pol)

calibrate_policy <- function(alpha, dshare, eshare = gp$policy_educ_share, spell = gp$edp_spell_q) {
  base36 <- neet36(simulate(cal, b$phi, b$Sb, alpha = alpha))
  target_total <- base36 - paths$reform_neet_rate[paths$fy_start_year == 2036] / 100
  tgt <- c(educ = eshare, demand = (1 - eshare) * dshare, supply = (1 - eshare) * (1 - dshare)) * target_total
  solve1 <- function(f, t) if (t <= 0) 0 else uniroot(function(z) (base36 - neet36(run_pol(f(z), alpha = alpha))) - t,
                                   c(0, 0.001), extendInt = "upX", tol = 1e-12)$root
  k <- c(solve1(function(z) mkpol(Ke = z, spell = spell), tgt[["educ"]]),
         solve1(function(z) mkpol(Kd = z, spell = spell), tgt[["demand"]]),
         solve1(function(z) mkpol(Ks = z, spell = spell), tgt[["supply"]]))
  names(k) <- c("educ", "demand", "supply")
  # Fit a smooth profile over time: one value per financial year with straight lines in between, chosen so
  # the package follows Milburn's reform path in every year. A penalty stops it jumping around year to year
  knots_t <- c(T0, paths$fy_start_year + 0.625)
  Pk <- function(v) approx(knots_t, c(0, v), xout = qtimes, rule = 2)$y
  tg_all <- paths$reform_neet_rate / 100
  obj <- function(v) {
    r <- fy_avg(run_pol(mkpol(k[["educ"]], k[["supply"]], k[["demand"]], P = Pk(v), spell = spell), alpha = alpha))
    1e4 * sum((as.numeric(r[as.character(paths$fy_start_year)]) - tg_all)^2) + 0.05 * sum(diff(v, differences = 2)^2)
  }
  v0 <- gap_fy / max(gap_fy)
  fit <- optim(v0, obj, method = "L-BFGS-B", lower = 0, upper = 5, control = list(maxit = 300))
  P <- Pk(fit$par)
  list(k = k, P = P, spell = spell)
}
# Education share: chosen so that, with no AI change, the package puts Milburn's 440,000 young people into
# work by 2036/37 (annex s.4 and s.8) and also cuts the NEET rate by his full amount.
# So matching the 440,000 is by design.
into_work_k <- function(kk, alpha) {
  r0 <- simulate(cal, b$phi, b$Sb, alpha = alpha)
  r1 <- run_pol(mkpol(kk$k[["educ"]], kk$k[["supply"]], kk$k[["demand"]], P = kk$P, spell = kk$spell), alpha = alpha)
  (fyv(r1, 2036, "E") - fyv(r0, 2036, "E")) * pop / 1000
}
es_fit <- uniroot(function(es) into_work_k(calibrate_policy(gp$alpha, gp$policy_demand_share, eshare = es), gp$alpha) -
                    gp$milburn_into_work_k, c(0.02, 0.7), tol = 1e-3)$root
cat("Education share fitted to Milburn's", gp$milburn_into_work_k, "k into work:", round(es_fit, 3), "\n")
gp$policy_educ_share <- es_fit
write_csv(tibble(param = "policy_educ_share_fitted", value = es_fit), "outputs/tables/fitted_education_share.csv")
alpha_grid <- seq(0.3, 0.7, by = 0.05); dshare_grid <- c(0.2, 0.3, 0.4)
# Fit the policies once for each value of alpha and each policy split used in the Monte Carlo
pol_cache <- expand_grid(alpha = alpha_grid, dshare = dshare_grid) %>% mutate(kk = map2(alpha, dshare, calibrate_policy))
get_kk <- function(a, d) pol_cache$kk[[which(abs(pol_cache$alpha - a) < 1e-9 & abs(pol_cache$dshare - d) < 1e-9)]]
kk_central <- get_kk(gp$alpha, gp$policy_demand_share)
cat("Policy volumes per quarter at full strength, people (central):",
    paste(names(kk_central$k), round(kk_central$k * pop), collapse = "; "), "\n")

levers <- list(package = function(kk, be) mkpol(kk$k[["educ"]], kk$k[["supply"]], kk$k[["demand"]], be, kk$P, kk$spell),
               educ    = function(kk, be) mkpol(Ke = kk$k[["educ"]], beta = be, P = kk$P, spell = kk$spell),
               demand  = function(kk, be) mkpol(Kd = kk$k[["demand"]], beta = be, P = kk$P, spell = kk$spell),
               supply  = function(kk, be) mkpol(Ks = kk$k[["supply"]], beta = be, P = kk$P, spell = kk$spell))

# Check: with no AI change, does the model follow Milburn's reform path (his Table 3)?
pkg0 <- fy_avg(run_pol(levers$package(kk_central, gp$beta_takeup), alpha = gp$alpha))
write_csv(tibble(fy = paths$fy, milburn_reform = paths$reform_neet_rate,
                 model_reform_no_ai = round(100 * as.numeric(pkg0[as.character(paths$fy_start_year)]), 2)),
          "outputs/tables/reform_path_check.csv")

# Run one scenario with no policies, with the full package and with each policy on its own.
# Returns yearly NEET rates and employment, plus job openings in 2030 and 2036
evaluate <- function(sp, g, alpha, dshare, beta = g$beta_takeup, kk = get_kk(alpha, dshare)) {
  ai <- ai_paths(sp, g, occ)
  nr <- simulate(cal, b$phi, b$Sb, ai$V, ai$sep, alpha)
  out <- list(no_reform = fy_avg(nr), emp_no_reform = fy_avg(nr, "E"))
  for (l in names(levers)) {
    r <- run_pol(levers[[l]](kk, beta), ai$V, ai$sep, alpha)
    out[[l]] <- fy_avg(r); out[[paste0("emp_", l)]] <- fy_avg(r, "E")
  }
  w36 <- fy_q == 2036
  out$V36 <- mean(ai$V[w36]); out$Vjobs36 <- mean(ai$V_jobs[w36]); out$L36 <- mean(ai$L[w36])
  out$V30 <- mean(ai$V[fy_q == 2030])
  out
}
# A scenario with no AI change at all, used as the comparison for the AI-proofing test
noai <- list(A_2028 = gp$A0, A_2030 = gp$A0, A_2036 = gp$A0, m_disp = 0, b_broaden = 0, g_comp = 0,
             r_realloc = 0, tau_realloc = 4, D_demand = 0, D_temp = 0, lag_sep = 4)
# Each policy's effect: how much it cuts the 2036/37 NEET rate compared with no policies
eff_of <- function(e) sapply(names(levers), function(l) e$no_reform[["2036"]] - e[[l]][["2036"]])
base_eff <- pol_cache %>% mutate(eff = map2(alpha, dshare, function(a, d) eff_of(evaluate(noai, gp, a, d))))
get_base_eff <- function(a, d) base_eff$eff[[which(abs(base_eff$alpha - a) < 1e-9 & abs(base_eff$dshare - d) < 1e-9)]]

# ---- 3. Scenario versions and the slump ----------------------------------------------
# The slump is an extra fall in youth job openings (on top of the scenario's own demand effect). It starts
# in 2029 and is sized so the youth NEET rate in 2030/31 matches the scenario's description. After mid-2030:
#   slump_fade : it lifts in a straight line by early 2037 (central case: a recession, then a recovery)
#   slump_held : it stays at its 2030 size until 2036 (upper bound)
# Runs where even an 85% extra fall in openings can't reach the target are dropped (counted in slump_size.csv).
neet_at <- function(sp, g, alpha, y, cl = cal, bl = b) {
  ai <- ai_paths(sp, g, occ); fyv(simulate(cl, bl$phi, bl$Sb, ai$V, ai$sep, alpha), y) }
# Find the size of slump that hits the 2030/31 target (NA if it would need more than slump_max)
solve_slump <- function(sp, g, alpha, scn, tr, cl = cal, bl = b) {
  tg <- targets[targets$scenario == scn, ]; sp$slump_w <- Inf
  f <- function(S) { sp$slump <- S; neet_at(sp, g, alpha, tg$target_fy_start, cl, bl) - tr }
  lo <- -gp$slump_max
  if (f(lo) < 0 || f(0.5) > 0) return(NA)
  uniroot(f, c(lo, 0.5), tol = 1e-7)$root
}
variants <- bind_rows(
  tibble(scenario = scen_names, variant = "base") %>% filter(!scenario %in% targets$scenario),
  expand_grid(scenario = targets$scenario, variant = c("jobs_only", "slump_fade", "slump_held")))
# Set up a scenario's assumptions for one version, sizing the slump if there is one
prep_sp <- function(sp, g, alpha, scn, variant, tr = NULL, cl = cal, bl = b) {
  sp$slump <- 0; sp$slump_w <- Inf
  if (grepl("^slump", variant)) {
    if (is.null(tr)) tr <- base_fy[["2030"]] + targets$target_change[targets$scenario == scn]
    S <- solve_slump(sp, g, alpha, scn, tr, cl, bl)
    if (is.na(S)) return(NULL)
    sp$slump <- S
    if (variant == "slump_fade")   sp$slump_w <- 2036.75 - 2030.5
  }
  sp
}

# Runs with every assumption at its central value
central_sp <- map(seq_len(nrow(variants)), function(i)
  prep_sp(sp_of(variants$scenario[i]), gp, gp$alpha, variants$scenario[i], variants$variant[i]))
central_eval <- map(seq_len(nrow(variants)), function(i) evaluate(central_sp[[i]], gp, gp$alpha, gp$policy_demand_share))
central <- map_dfr(seq_len(nrow(variants)), function(i) { e <- central_eval[[i]]
  tibble(scenario = variants$scenario[i], variant = variants$variant[i], fy = names(e$no_reform),
         no_reform = as.numeric(e$no_reform), package = as.numeric(e$package),
         emp_no_reform = as.numeric(e$emp_no_reform), emp_package = as.numeric(e$emp_package)) })
write_csv(central, "outputs/tables/central_paths.csv")
central_V <- map_dfr(seq_len(nrow(variants)), function(i) tibble(
  scenario = variants$scenario[i], variant = variants$variant[i], D_demand = central_sp[[i]]$D_demand,
  slump = central_sp[[i]]$slump, slump_unwind_years = central_sp[[i]]$slump_w,
  openings_2030_vs_noAI = central_eval[[i]]$V30, openings_2036_vs_noAI = central_eval[[i]]$V36,
  openings_jobs_channel_only_2036 = central_eval[[i]]$Vjobs36))
write_csv(central_V, "outputs/tables/central_openings_and_slump.csv")
print(central_V)

# ---- 4. Monte Carlo: run each scenario 2,000 times with different assumptions --------------------------------------------------------------------
# Pick a random value for each scenario assumption from within its range
draw_sp <- function(s) {
  d <- scen_long %>% filter(scenario == s)
  v <- runif(nrow(d), d$low, d$high); names(v) <- d$param
  v[c("tau_realloc", "lag_sep")] <- round(v[c("tau_realloc", "lag_sep")])
  v["tau_realloc"] <- max(1, v["tau_realloc"])
  as.list(v)
}
rng <- function(p) c(gpd$low[gpd$param == p], gpd$high[gpd$param == p])
# Pick random values for the assumptions that apply to every scenario
draw_g <- function() {
  g <- gp
  g$G_anchor  <- runif(1, rng("G_anchor")[1], rng("G_anchor")[2])
  g$sigma_sep <- runif(1, rng("sigma_sep")[1], rng("sigma_sep")[2])
  g$beta_takeup <- runif(1, rng("beta_takeup")[1], rng("beta_takeup")[2])
  g$alpha  <- sample(alpha_grid, 1)
  g$dshare <- sample(dshare_grid, 1)
  g
}
N <- as.integer(Sys.getenv("NDRAWS", gp$n_draws))
t_mc <- Sys.time()
mc <- map_dfr(seq_len(nrow(variants)), function(vi) {
  s <- variants$scenario[vi]; vr <- variants$variant[vi]
  set.seed(gp$seed + match(s, scen_names))   # all versions of a scenario use the same random draws, so they can be compared like for like
  map_dfr(seq_len(N), function(i) {
    sp <- draw_sp(s); g <- draw_g()
    sp <- prep_sp(sp, g, g$alpha, s, vr)
    if (is.null(sp)) return(NULL)
    e <- evaluate(sp, g, g$alpha, g$dshare, g$beta_takeup)
    # Share of each policy's no-AI effect it keeps in this run
    ret <- eff_of(e) / get_base_eff(g$alpha, g$dshare)
    tibble(scenario = s, variant = vr, draw = i, fy = names(e$no_reform), no_reform = as.numeric(e$no_reform),
           package = as.numeric(e$package), emp_no_reform = as.numeric(e$emp_no_reform),
           emp_package = as.numeric(e$emp_package),
           ret_package = ret[["package"]], ret_educ = ret[["educ"]], ret_demand = ret[["demand"]], ret_supply = ret[["supply"]],
           emp_effect_package = as.numeric(e$emp_package) - as.numeric(e$emp_no_reform),
           V36 = e$V36, Vjobs36 = e$Vjobs36, L36 = e$L36, V30 = e$V30, slump = sp$slump, slump_w = sp$slump_w,
           !!!sp[c("A_2036", "m_disp", "b_broaden", "g_comp", "r_realloc", "D_demand", "tau_realloc", "lag_sep")],
           G_anchor = g$G_anchor, sigma_sep = g$sigma_sep, beta_takeup = g$beta_takeup, alpha = g$alpha, dshare = g$dshare)
  })
})
cat("Monte Carlo:", N, "draws x", nrow(variants), "scenario variants in",
    round(difftime(Sys.time(), t_mc, units = "secs")), "s\n")
saveRDS(mc, "outputs/tables/mc_draws.rds")

# Median and 10th-90th percentile range across runs, for each year
bands <- mc %>% group_by(scenario, variant, fy) %>%
  summarise(across(c(no_reform, package, emp_no_reform, emp_package),
                   list(p10 = ~ q(.x, .1), p50 = ~ q(.x, .5), p90 = ~ q(.x, .9))), .groups = "drop")
write_csv(bands, "outputs/tables/scenario_bands.csv")

# Size of the slump across runs, and how many runs were dropped
slump <- mc %>% filter(fy == "2036") %>% group_by(scenario, variant) %>%
  summarise(runs_kept = n(), runs_dropped = N - n(),
            slump_p10 = q(slump, .1), slump_p50 = q(slump, .5), slump_p90 = q(slump, .9),
            unwind_years_p50 = q(slump_w, .5),
            openings_2030_p50 = q(V30, .5), openings_2036_p50 = q(V36, .5),
            openings_jobs_channel_2036_p50 = q(Vjobs36, .5), .groups = "drop")
write_csv(slump, "outputs/tables/slump_size.csv")
print(slump)

# ---- 5. AI-proofing test: does each policy keep at least 80% of its no-AI effect? ------------------------------------------------------------------
rob <- mc %>% filter(fy == "2036") %>%
  pivot_longer(starts_with("ret_"), names_to = "lever", values_to = "retained") %>%
  group_by(scenario, variant, lever) %>%
  summarise(p10 = q(retained, .1), p50 = q(retained, .5), p90 = q(retained, .9), .groups = "drop") %>%
  mutate(lever = recode(lever, ret_package = "Full package", ret_educ = "Education and training places",
                        ret_demand = "Subsidised youth hires (demand)", ret_supply = "Activation, health and welfare (supply)"))
rob_test <- bind_rows(
    rob %>% filter(variant %in% c("base", "jobs_only")) %>% mutate(set = "AI in jobs only"),
    rob %>% filter(variant %in% c("base", "slump_held")) %>% mutate(set = "With held slump in S4/S5"),
    rob %>% filter(variant %in% c("base", "slump_fade")) %>% mutate(set = "With fading slump in S4/S5")) %>%
  group_by(set, lever) %>%
  summarise(worst_p10 = min(p10), worst_p50 = min(p50), worst_scenario = scenario[which.min(p10)],
            passes_strict_p10 = all(p10 >= gp$pass_threshold), passes_median = all(p50 >= gp$pass_threshold),
            .groups = "drop")
write_csv(rob, "outputs/tables/policy_robustness.csv")
write_csv(rob_test, "outputs/tables/ai_proofing_test.csv")
print(rob %>% filter(scenario %in% targets$scenario), n = 40); print(rob_test, n = 20)

# How many extra young people the package gets into work by 2036/37 (thousands)
emp_eff <- mc %>% filter(fy == "2036") %>% group_by(scenario, variant) %>%
  summarise(emp_effect_p10 = q(emp_effect_package, .1), emp_effect_p50 = q(emp_effect_package, .5),
            emp_effect_p90 = q(emp_effect_package, .9), .groups = "drop") %>%
  mutate(across(starts_with("emp"), ~ round(.x * pop / 1000)))
e0 <- evaluate(noai, gp, gp$alpha, gp$policy_demand_share)
emp_eff <- bind_rows(tibble(scenario = "No AI change", variant = "none",
                            emp_effect_p50 = round((e0$emp_package[["2036"]] - e0$emp_no_reform[["2036"]]) * pop / 1000)), emp_eff)
write_csv(emp_eff, "outputs/tables/package_employment_effect_k.csv")
print(emp_eff)

# ---- 6. Breaking point: how big a permanent fall in youth job openings (from 2027) each policy can cope with ---------------------
be_c <- get_base_eff(gp$alpha, gp$policy_demand_share)
brk <- map_dfr(seq(0, 0.6, by = 0.05), function(cut) {
  V <- ifelse(qtimes >= 2027, 1 - cut, 1)
  nr <- neet36(simulate(cal, b$phi, b$Sb, V, rep(1, NQ), gp$alpha))
  map_dfr(c(0, 0.5, 1), function(beta) map_dfr(names(levers), function(l) {
    wp <- neet36(run_pol(levers[[l]](kk_central, beta), V, rep(1, NQ), gp$alpha))
    tibble(hiring_cut = cut, beta_takeup = beta, lever = l, neet36_no_reform = nr, retained = (nr - wp) / be_c[[l]])
  }))
})
write_csv(brk, "outputs/tables/policy_breaking_point.csv")
print(brk %>% filter(beta_takeup == 0.5) %>% select(hiring_cut, lever, retained) %>%
        pivot_wider(names_from = lever, values_from = retained), n = 20)
print(brk %>% distinct(hiring_cut, neet36_no_reform), n = 20)

# ---- 7. Headline results for 2036/37 ---------------------------------------------------------------------
head36 <- mc %>% filter(fy == "2036") %>% group_by(scenario, variant) %>%
  summarise(neet_no_reform_p50 = q(no_reform, .5), neet_no_reform_p10 = q(no_reform, .1),
            neet_no_reform_p90 = q(no_reform, .9),
            neet_package_p50 = q(package, .5), neet_package_p10 = q(package, .1), neet_package_p90 = q(package, .9),
            .groups = "drop") %>%
  mutate(across(starts_with("neet"), ~ round(100 * .x, 1)),
         extra_neet_vs_milburn_baseline_k = round((neet_no_reform_p50 / 100 - base_fy[["2036"]]) * pop / 1000),
         extra_neet_vs_milburn_reform_k   = round((neet_package_p50 / 100 - 0.068) * pop / 1000))
h30 <- mc %>% filter(fy == "2030") %>% group_by(scenario, variant) %>%
  summarise(neet30_no_reform_p50 = round(100 * q(no_reform, .5), 1), .groups = "drop")
head36 <- head36 %>% left_join(h30, by = c("scenario", "variant"))
write_csv(head36, "outputs/tables/headline_2036_37.csv")
print(head36, width = 220)

# ---- 8. Sensitivity checks -----------------------------------------------------------------------
# One at a time: move each assumption to the low and high end of its range, keeping the others at their
# central values. The slump is re-sized each time
tornado <- map_dfr(seq_len(nrow(variants)), function(vi) {
  s <- variants$scenario[vi]; vr <- variants$variant[vi]
  if (vr == "slump_held") return(NULL)
  c0 <- evaluate(prep_sp(sp_of(s), gp, gp$alpha, s, vr), gp, gp$alpha, gp$policy_demand_share)
  rows <- list()
  sc <- scen_long %>% filter(scenario == s, low != high)
  for (i in seq_len(nrow(sc))) for (side in c("low", "high")) {
    sp <- sp_of(s); sp[[sc$param[i]]] <- sc[[side]][i]
    if (sc$param[i] %in% c("tau_realloc", "lag_sep")) sp[[sc$param[i]]] <- round(sp[[sc$param[i]]])
    spp <- prep_sp(sp, gp, gp$alpha, s, vr); if (is.null(spp)) next
    e <- evaluate(spp, gp, gp$alpha, gp$policy_demand_share)
    rows[[length(rows) + 1]] <- tibble(scenario = s, variant = vr, param = sc$param[i], side = side, value = sc[[side]][i],
                                       neet36 = e$no_reform[["2036"]], central = c0$no_reform[["2036"]])
  }
  for (pn in c("G_anchor", "sigma_sep", "alpha")) for (side in c("low", "high")) {
    g <- gp; g[[pn]] <- gpd[[side]][gpd$param == pn]
    spp <- prep_sp(sp_of(s), g, g$alpha, s, vr); if (is.null(spp)) next
    e <- evaluate(spp, g, g$alpha, gp$policy_demand_share)
    rows[[length(rows) + 1]] <- tibble(scenario = s, variant = vr, param = pn, side = side, value = g[[pn]],
                                       neet36 = e$no_reform[["2036"]], central = c0$no_reform[["2036"]])
  }
  bind_rows(rows)
})
write_csv(tornado, "outputs/tables/sensitivity_oneway.csv")

# Which assumptions matter most across all Monte Carlo runs (rank correlation with the 2036/37 NEET rate)
gsa <- mc %>% filter(fy == "2036") %>% group_by(scenario, variant) %>%
  summarise(across(c(A_2036, m_disp, b_broaden, g_comp, r_realloc, D_demand, tau_realloc, lag_sep,
                     G_anchor, sigma_sep, alpha),
                   ~ suppressWarnings(cor(.x, no_reform, method = "spearman"))), .groups = "drop")
write_csv(gsa, "outputs/tables/sensitivity_rank_correlation.csv")

# How the policy split and the length of funded education places change the results.
# Education shares: 12% (if education delivered only the gap between 500,000 fewer NEETs and 440,000 more
# in work), 20% (the same gap using the annex's rounder 400,000) and 36% (education's share of the costs).
# Places last one year, two years or indefinitely. Run for S4 and S5 with a slump that lifts
split_sens <- map_dfr(c(0.12, 0.20, 0.36), function(es) map_dfr(c(4, 8, Inf), function(spq) {
  kk <- calibrate_policy(gp$alpha, gp$policy_demand_share, eshare = es, spell = spq)
  e0s <- evaluate(noai, gp, gp$alpha, gp$policy_demand_share, kk = kk); be <- eff_of(e0s)
  map_dfr(targets$scenario, function(s) {
    e <- evaluate(prep_sp(sp_of(s), gp, gp$alpha, s, "slump_fade"), gp, gp$alpha, gp$policy_demand_share, kk = kk)
    r <- eff_of(e) / be
    tibble(educ_share = es, edu_spell_quarters = spq, scenario = s, neet36_package = e$package[["2036"]],
           ret_package = r[["package"]], ret_educ = r[["educ"]], ret_demand = r[["demand"]], ret_supply = r[["supply"]],
           into_work_no_ai_k = round((e0s$emp_package[["2036"]] - e0s$emp_no_reform[["2036"]]) * pop / 1000))
  })
}))
write_csv(split_sens, "outputs/tables/sensitivity_lever_split_and_spell.csv")
print(split_sens, n = 30)

# How the fitted youth flow rates change the results: a lower and higher share of education leavers going
# straight into work, and a test that sets youth job finding to the all-age rate
calib_sens <- map_dfr(list(list(piE_leaver = 0.3), list(piE_leaver = 0.7), list(kfind_fixed = 1)), function(ch) {
  g <- gp; for (n in names(ch)) g[[n]] <- ch[[n]]
  cs <- calibrate_ss(stocks, haz, entry, q2324, g); bs <- calibrate_baseline(cs, CB$tq)
  map_dfr(scen_names, function(s) map_dfr(if (s %in% targets$scenario) c("jobs_only", "slump_fade", "slump_held") else "base", function(vr) {
    sp <- prep_sp(sp_of(s), gp, gp$alpha, s, vr, cl = cs, bl = bs)
    if (is.null(sp)) return(NULL)
    ai <- ai_paths(sp, gp, occ)
    r <- simulate(cs, bs$phi, bs$Sb, ai$V, ai$sep, gp$alpha)
    tibble(change = paste(names(ch), unlist(ch), sep = " = "), scenario = s, variant = vr,
           neet30 = fyv(r, 2030), neet36 = fyv(r, 2036), slump = sp$slump, slump_w = sp$slump_w,
           ksep = cs$ksep, xexit = cs$xexit, kfind = cs$kfind, kui = cs$kui, fit = cs$fit)
  }))
})
write_csv(calib_sens, "outputs/tables/sensitivity_calibration.csv")
print(calib_sens, n = 30)

# Two ways of sizing the slump. Central: as a rise in the youth NEET rate (the 2008-11 average rise of
# 1.7 points for S4 and the rise to the 2011 peak of 3.0 points for S5), added to the model's no-AI
# 2030/31 rate. Sensitivity: as the 2009-11 levels instead (15.6% and 16.9%)
base30 <- base_fy[["2030"]]
be_c0 <- get_base_eff(gp$alpha, gp$policy_demand_share)
target_sens <- map_dfr(list(list(lab = "change +1.7 / +3.0 pts (central)", S4 = NA, S5 = NA),
                            list(lab = "level 15.6% / 16.9%", S4 = 0.156, S5 = 0.169)), function(tt)
  map_dfr(targets$scenario, function(s) map_dfr(c("slump_fade", "slump_held"), function(vr) {
    tr <- tt[[substr(s, 1, 2)]]; if (is.na(tr)) tr <- NULL
    sp <- prep_sp(sp_of(s), gp, gp$alpha, s, vr, tr = tr); if (is.null(sp)) return(NULL)
    e <- evaluate(sp, gp, gp$alpha, gp$policy_demand_share)
    tibble(targets = tt$lab, scenario = s, variant = vr, target30 = if (is.null(tr)) NA else tr,
           neet30 = e$no_reform[["2030"]], neet36 = e$no_reform[["2036"]], package36 = e$package[["2036"]],
           ret_package = (eff_of(e) / be_c0)[["package"]], slump = sp$slump, openings_2030 = e$V30)
  })))
write_csv(target_sens, "outputs/tables/sensitivity_targets.csv")
print(target_sens)

# ---- 9. What the UK data shows so far: youth employment by AI exposure group (APS) ----
# Occupations are grouped by how much AI could replace their work: high (score 0.3 or more), some (0.05 to
# 0.3) and low. Shown from two start years: 2021 is still affected by Covid, and 2022 is the last full year
# before generative AI tools became widely used (from late 2022).
aps <- read_csv("data/derived/aps_soc3_by_year.csv", show_col_types = FALSE) %>% filter(!is.na(s))
aps_grp <- aps %>% mutate(group = case_when(s >= 0.3 ~ "High substitution exposure",
                                            s >= 0.05 ~ "Some exposure", TRUE ~ "Low exposure")) %>%
  group_by(year, group) %>% summarise(older = sum(all - y1624), y1624 = sum(y1624),
                                      youth_cells_suppressed = sum(y_suppressed), n_occ = n(), .groups = "drop")
aps_chk <- aps_grp %>% group_by(group) %>%
  mutate(idx_young = 100 * y1624 / y1624[year == 2022], idx_older = 100 * older / older[year == 2022],
         idx_young_2021 = 100 * y1624 / y1624[year == 2021], idx_older_2021 = 100 * older / older[year == 2021]) %>% ungroup()
write_csv(aps_chk, "outputs/tables/uk_aps_exposure_check.csv")
print(aps_chk, width = 200)

# Save the main inputs the charts need
saveRDS(list(base_fy = base_fy, paths = paths, hist = hist, pop = pop, cal = cal, kk = kk_central,
             variants = variants),
        "outputs/tables/run_meta.rds")
cat("Done.\n")
