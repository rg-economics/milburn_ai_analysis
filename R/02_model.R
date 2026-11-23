# 02_model.R ----------------------------------------------------------------
# The model of 16-24 year olds, with a channel for AI to cut hiring.
# Every young person is in one of four states (stored as shares of all 16-24 year olds):
#   E  = in work (and not in education)
#   Ed = in education or training
#   U  = NEET and unemployed (looking for work)
#   I  = NEET and inactive (not looking for work)
# Each step is one quarter. Time is written as a decimal year, so Apr-Jun 2026 = 2026.25.

T0 <- 2026.25; T1 <- 2037.00
qtimes <- seq(T0, T1, by = 0.25)
NQ <- length(qtimes)

# ---- How AI affects job openings ----------------------------------------------
adoption_path <- function(sp, A0) {
  # Share of firms using AI each quarter: straight lines between the scenario's values
  # for 2028, 2030 and 2036 (placed at mid-year)
  approx(x = c(T0, 2028.5, 2030.5, 2036.75, 2040),
         y = c(A0, sp$A_2028, sp$A_2030, sp$A_2036, sp$A_2036), xout = qtimes, rule = 2)$y
}

ai_paths <- function(sp, gp, occ) {
  # Works out, for each quarter:
  #   V   = youth job openings compared with a world with no AI change (1 = no change)
  #   sep = how much more often existing young workers lose their jobs (1 = no change)
  A  <- adoption_path(sp, gp$A0)
  dA <- pmax(0, A - gp$A0)
  d  <- gp$G_anchor / gp$A0                         # cut in openings for each extra unit of AI use
  s_eff <- pmin(1, occ$s + sp$b_broaden * occ$c)
  c_eff <- occ$c * (1 - sp$b_broaden)
  w <- occ$w_youth
  # Work out each occupation's cut and gain (rows = quarters, columns = occupations).
  # s_eff: how much AI replaces the work; c_eff: how much it helps with it.
  # b_broaden lets AI spread from helping with tasks to replacing them
  cut  <- outer(dA, s_eff) * d * sp$m_disp; cut[] <- pmin(cut, 0.95)  # cap any one cut at 95%
  gain <- outer(dA, c_eff) * sp$g_comp
  L    <- as.vector(cut %*% w)                      # share of youth hiring lost, weighted by where young people work
  G    <- as.vector(gain %*% w)                     # extra hiring where AI helps with the work or creates new tasks
  # Some of the lost hiring reappears in other, less exposed jobs, but only gradually
  R <- numeric(NQ); for (t in 2:NQ) R[t] <- R[t - 1] + (sp$r_realloc * L[t] - R[t - 1]) / sp$tau_realloc
  # Each scenario's wider effect on demand builds up as AI use grows to its 2030 level, then stays there.
  # (The GO-Science outlooks suggest unemployment may stop rising after 2030 as workers move into new jobs.
  # Any further damage still comes through the hiring cuts above.)
  dA_30 <- max(sp$A_2030 - gp$A0, 1e-6)
  build <- pmin(1, dA / dA_30)
  # The "slump" is an extra fall in youth openings on top of the scenario's own demand effect.
  # It is only used in Transformation Economy and Take-Off. It builds in a straight line from slump_start
  # (2029, when both scenarios say AI reaches human level) to mid-2030. After that it either stays
  # in place (slump_w = Inf) or lifts in a straight line over slump_w years.
  slump <- if (is.null(sp$slump)) 0 else sp$slump
  s0 <- if (is.null(gp$slump_start)) 2029 else gp$slump_start
  w_un <- if (is.null(sp$slump_w)) Inf else sp$slump_w
  sbuild <- pmin(1, pmax(0, (qtimes - s0) / (2030.5 - s0)))
  hold <- if (is.finite(w_un)) 1 - pmin(1, pmax(0, (qtimes - 2030.5) / w_un)) else rep(1, NQ)
  Dm <- (1 + sp$D_demand * build) * (1 + slump * sbuild * hold) + sp$D_temp * (qtimes >= 2028 & qtimes < 2030)
  V  <- pmax(0.05, (1 - L + G + R) * Dm)
  lag <- round(sp$lag_sep)
  Llag <- c(rep(0, lag), L)[1:NQ]
  sep <- 1 + gp$sigma_sep * Llag                     # job losses rise by sigma_sep times the hiring cut, after a delay
  V_jobs <- pmax(0.05, 1 - L + G + R)                # openings from AI in young people's own jobs only (no demand effect or slump)
  list(A = A, L = L, G = G, R = R, Dm = Dm, V = V, V_jobs = V_jobs, sep = sep)
}

# ---- Running the model forward ------------------------------------------------
# cal:    the fitted flow rates (from calibrate_ss)
# phi:    how easy it is to find a job each quarter, set so the no-AI run follows Milburn's baseline
# Sb:     number of job seekers each quarter in that baseline run (more seekers means more competition)
# V, sep: job openings and job losses from ai_paths (1 = no AI change)
# pol:    Milburn's policies (see 03_run.R). NULL means no policies
simulate <- function(cal, phi, Sb, V = rep(1, NQ), sep = rep(1, NQ), alpha = 0.5,
                     pol = NULL, x0 = cal$x0, keep = TRUE) {
  x <- matrix(NA_real_, NQ, 4, dimnames = list(NULL, c("E", "Ed", "U", "I")))
  x[1, ] <- x0
  S_rec <- numeric(NQ)
  h <- cal$h
  edp <- 0   # young people currently on a Milburn-funded education or training place (counted within Ed)
  edp_rec <- numeric(NQ)
  P <- if (is.null(pol)) rep(0, NQ) else pol$P
  for (t in 1:NQ) {
    E <- x[t, 1]; Ed <- x[t, 2]; U <- x[t, 3]; I <- x[t, 4]
    # Job seekers: the unemployed, plus inactive people and education leavers, weighted by how
    # likely they are to find work compared with the unemployed
    S  <- U + cal$wI * I + cal$wEd * Ed
    S_rec[t] <- S
    if (t == NQ) break
    # Job finding goes up with openings and down with more competition from other seekers
    # (a standard matching function). F = 1 means as easy as in the no-AI baseline
    F  <- phi[t] * V[t]^alpha * (S / Sb[t])^(-alpha)
    ue <- min(0.95, h$UE * F); ie <- min(0.95, h$IE * F)
    pE <- min(0.95, cal$piE * F)
    eu <- h$EU * cal$ksep * sep[t]; ei <- h$EI * cal$ksep * sep[t]
    iu <- h$IU
    ui <- h$UI
    rU <- cal$rU; rI <- cal$rI; rE <- cal$rE
    xo <- cal$xexit
    edp_rec[t] <- edp
    # Number of people moving between each pair of states this quarter
    fEU <- eu * E; fEI <- ei * E; fEEd <- rE * E
    fUE <- ue * U; fUI <- ui * U; fUEd <- rU * U
    fIE <- ie * I; fIU <- iu * I; fIEd <- rI * I
    # People on funded places leave when their place ends; everyone else leaves education at the usual rate
    sp_q <- if (is.null(pol$edp_spell)) Inf else pol$edp_spell
    lev_p <- if (is.finite(sp_q)) edp / sp_q else xo * edp
    lev <- xo * (Ed - edp) + lev_p
    edp <- edp - lev_p
    fEdE <- pE * lev; fEdU <- (1 - pE) * lev * cal$toU; fEdI <- (1 - pE) * lev * (1 - cal$toU)
    # New 16 year olds join and 25 year olds leave (the shares always add up to 1)
    inflow <- cal$entry / 36
    ow <- cal$omega * x[t, ]; out <- ow / sum(ow) / 36
    if (Ed > 0) edp <- edp * (1 - out[2] / Ed)        # funded learners also leave the model at 25
    nE  <- E  + fUE + fIE + fEdE - fEU - fEI - fEEd
    nEd <- Ed + fEEd + fUEd + fIEd - lev
    nU  <- U  + fEU + fIU + fEdU - fUE - fUI - fUEd
    nI  <- I  + fEI + fUI + fEdI - fIE - fIU - fIEd
    # Milburn's policies: fixed numbers of places and hires each quarter, scaled up over time by P.
    # Together they can never take more than half of all NEETs in one quarter
    if (!is.null(pol) && !is.null(pol$Kd)) {
      neet <- U + I
      relF <- F / phi[t]                                  # job finding compared with the no-AI baseline
      hD <- pol$Kd * P[t] * max(V[t], 0.05)^pol$beta      # subsidised hires: take-up falls as openings fall
      hS <- pol$Ks * P[t] * relF                          # activation into work: depends on job finding
      eE <- pol$Ke * P[t]                                 # education / training places: fixed
      tot <- hD + hS + eE
      if (tot > 0.5 * neet) { sc <- 0.5 * neet / tot; hD <- hD * sc; hS <- hS * sc; eE <- eE * sc }
      su <- if (neet > 0) U / neet else 0.5
      nE <- nE + hD + hS; nEd <- nEd + eE; edp <- edp + eE
      nU <- nU - su * (hD + hS + eE); nI <- nI - (1 - su) * (hD + hS + eE)
    }
    x[t + 1, ] <- c(nE, nEd, nU, nI) + inflow - out
  }
  if (!keep) return(x[NQ, ])
  data.frame(t = qtimes, x, neet = x[, "U"] + x[, "I"], S = S_rec, edp = edp_rec)
}

# ---- Fitting the model to today's data ----------------------------------------
# The ONS doesn't publish flow rates for young people, so we start from the all-age rates and fit four
# numbers so the model reproduces today's numbers of young people in each state
calibrate_ss <- function(stocks, haz, entry, q2324, gp) {
  x0 <- setNames(stocks$share, stocks$state)[c("E", "Ed", "U", "I")]
  h  <- list(EU = haz$EU, EI = haz$EI, UE = haz$UE, UI = haz$UI, IE = haz$IE, IU = haz$IU)
  q  <- c(E = 1 - q2324$neet_U - q2324$neet_I - gp$ed_share_23_24, Ed = gp$ed_share_23_24,
          U = q2324$neet_U, I = q2324$neet_I)
  base <- list(h = h, x0 = x0, wI = h$IE / h$UE, toU = gp$ed_exit_to_U,
               rU = gp$ret_ed_U, rI = gp$ret_ed_I, rE = gp$ret_ed_E,
               entry = setNames(entry$share, entry$state)[c("E", "Ed", "U", "I")],
               omega = q / x0)
  mk <- function(th) {
    cal <- base
    cal$ksep  <- exp(th[1])                 # young people's job loss rate compared with all ages
    cal$xexit <- plogis(th[2])              # share of people leaving education each quarter
    # young people's job finding rate compared with all ages (can be fixed at 1 as a sensitivity)
    cal$kfind <- if (is.null(gp$kfind_fixed) || is.na(gp$kfind_fixed)) exp(th[3]) else gp$kfind_fixed
    cal$kui   <- exp(th[4])                 # young people's unemployed-to-inactive rate compared with all ages
    cal$piE   <- gp$piE_leaver              # share of education leavers going straight into work
    cal$h$UE <- h$UE * cal$kfind; cal$h$IE <- h$IE * cal$kfind; cal$h$UI <- h$UI * cal$kui
    cal$wEd <- cal$xexit * cal$piE / cal$h$UE
    cal
  }
  ss <- function(cal) {  # run the model forward until it settles down, with no changes to openings
    x <- x0
    one <- rep(1, NQ)
    for (k in 1:30) x <- simulate(cal, one, rep(1e9, NQ), alpha = 0, x0 = x, keep = FALSE)
    x
  }
  # Fit to 2026. There are four numbers but only three independent targets, so a small penalty
  # keeps the multipliers close to 1 and picks between solutions that fit almost equally well
  obj <- function(th) { x <- ss(mk(th)); sum(((x - x0) / x0)^2) + 1e-4 * sum(th[c(1, 3, 4)]^2) }
  fit <- optim(c(0, qlogis(0.05), 0, 0), obj, control = list(maxit = 4000, reltol = 1e-14))
  cal <- mk(fit$par); cal$fit <- fit$value; cal$ss <- ss(cal)
  cal
}

calibrate_baseline <- function(cal, target_q, alpha = 0.5) {
  # Set how easy it is to find a job (phi), quarter by quarter, so the run with no AI and no policies
  # follows Milburn's baseline NEET rate exactly. This match is by design, not a test of the model.
  phi <- rep(1, NQ); Sb <- rep(1, NQ)
  x <- cal$x0
  for (t in 1:(NQ - 1)) {
    f <- function(p) {
      ph <- phi; ph[t] <- p
      xx <- simulate_step(cal, x, ph[t])
      (xx["U"] + xx["I"]) - target_q[t + 1]
    }
    phi[t] <- uniroot(f, c(0.05, 6), tol = 1e-10)$root
    S <- (x["U"] + cal$wI * x["I"] + cal$wEd * x["Ed"]); Sb[t] <- S
    x <- simulate_step(cal, x, phi[t])
  }
  Sb[NQ] <- x["U"] + cal$wI * x["I"] + cal$wEd * x["Ed"]
  phi[NQ] <- phi[NQ - 1]
  list(phi = phi, Sb = as.numeric(Sb))
}

simulate_step <- function(cal, x, phi_t) {
  # Runs one quarter with no AI change. Used only when fitting the baseline above
  ph <- rep(phi_t, NQ)
  S  <- x["U"] + cal$wI * x["I"] + cal$wEd * x["Ed"]
  out <- simulate(cal, ph, rep(S, NQ), alpha = 0, x0 = x, keep = TRUE)
  unlist(out[2, c("E", "Ed", "U", "I")])
}

fy_avg <- function(df, col = "neet") {
  fy <- floor(df$t - 0.25)       # quarters from Apr of year y to Jan of year y+1 count as financial year y
  tapply(df[[col]], fy, mean)
}

quarterly_target <- function(paths, start_rate) {
  # Turns Milburn's yearly baseline into a quarterly path. Each financial year's value sits at the middle
  # of that year, starting from the actual Apr-Jun 2026 rate
  approx(x = c(T0, paths$fy_start_year + 0.625), y = c(start_rate, paths$baseline_neet_rate) / 100,
         xout = qtimes, rule = 2)$y
}
