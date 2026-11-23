# 01_data.R -----------------------------------------------------------------
# Builds every input the model needs from the published files in data/raw.
# Everything here is public data (ONS, ILO and DWP's Milburn Review tables). No microdata.
# Results are saved to data/derived/*.csv

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(readr)
})

raw <- function(f) file.path("data", "raw", f)
der <- function(f) file.path("data", "derived", f)
num <- function(x) suppressWarnings(as.numeric(x))

# 1. How many young people are in each state, UK, Apr-Jun 2026 ------------------------------------------
neet <- read_excel(raw("neet.xlsx"), sheet = "People - SA", col_names = FALSE, .name_repair = "minimal")
neet <- neet[, 1:16]; names(neet) <- c("period", paste0("v", 1:15))
row  <- neet %>% filter(str_trim(period) == "Apr-Jun 2026")
stopifnot(nrow(row) == 1)
pop      <- num(row$v4)            # 16-24 population
neet_U   <- num(row$v2)            # NEET unemployed 16-24
neet_I   <- num(row$v3)            # NEET inactive 16-24
neet1617 <- c(U = num(row$v7), I = num(row$v8), pop = num(row$v9))

a06 <- read_csv(der("a06_youth.csv"), show_col_types = FALSE)  # made by 00_extract_a06.R
emp_not_fte <- a06 %>% filter(period == "Apr-Jun 2026") %>% pull(AGNJ)

stocks <- tibble(
  state = c("E", "Ed", "U", "I"),
  level = c(emp_not_fte, pop - emp_not_fte - neet_U - neet_I, neet_U, neet_I)
) %>% mutate(share = level / pop)
write_csv(stocks, der("stocks_2026q2.csv"))

# Where 16 year olds start out when they join the model (based on NEET rates for 16-17 year olds)
entry <- tibble(state = c("E", "Ed", "U", "I"),
                share = c(0, 1 - (neet1617["U"] + neet1617["I"]) / neet1617["pop"],
                          neet1617["U"] / neet1617["pop"], neet1617["I"] / neet1617["pop"]))
write_csv(entry, der("entry_at_16.csv"))

# Share of 23-24 year olds who are unemployed or inactive (average of latest four quarters, not seasonally adjusted).
# Used to work out who leaves the model when they turn 25
n2324 <- read_excel(raw("neet.xlsx"), sheet = "People 23-24 - NSA", col_names = FALSE, .name_repair = "minimal")
n2324 <- n2324[, 1:6]; names(n2324) <- c("period", "neet", "U", "I", "pop", "rate")
n2324 <- n2324 %>% mutate(across(-period, num)) %>% filter(!is.na(pop)) %>% tail(4)
write_csv(tibble(neet_U = mean(n2324$U / n2324$pop), neet_I = mean(n2324$I / n2324$pop)),
          der("neet_23_24.csv"))

# Past NEET rates, used for the starting point and the charts
hist <- neet %>% transmute(period = str_trim(period), neet_rate = num(v5)) %>%
  filter(!is.na(neet_rate), !is.na(period), grepl("^[A-Z][a-z]{2}-[A-Z][a-z]{2} [0-9]{4}", period))
write_csv(hist, der("neet_history.csv"))

# 2. How often people move between work, unemployment and inactivity each quarter -----
# All ages 16-64 (ONS X02). We take the average of the latest 8 quarters
x02 <- read_excel(raw("x02.xls"), sheet = "Labour market flows hazard rate",
                  col_names = FALSE, .name_repair = "minimal")[, c(1:4, 6:8, 10:12)]  # columns 5 and 9 are blank
names(x02) <- c("period", "EE", "EU", "EI", "UU", "UE", "UI", "II", "IE", "IU")
x02 <- x02 %>% mutate(across(-period, num)) %>% filter(!is.na(EE))
haz <- x02 %>% tail(8) %>% summarise(across(-period, ~ mean(.x) / 100))
haz$window <- paste(x02$period[nrow(x02) - 7], "to", x02$period[nrow(x02)])
write_csv(haz, der("hazards_all_ages.csv"))

# 3. AI exposure scores for each occupation: ILO Working Paper 140, Table A1 ---------------------
# Occupations here use the international ISCO-08 codes (4-digit)
txt <- readLines(raw("wp140.txt"))
# In the PDF text some exposure labels run onto two lines, so a row can start with just
# "Exposure" (the end of "Minimal Exposure"). We accept that too so no occupations are missed.
pat <- "^\\s*(Gradient \\d|Minimal Exposure|Not Exposed|Exposure)?\\s*(\\d{4})\\s{2,}(.+?)\\s{2,}(\\d?\\.\\d+|0|1)\\s+(\\d?\\.\\d+|0)\\s*$"
m <- str_match(txt, pat)
ilo <- tibble(isco4 = m[, 3], title = str_trim(m[, 4]), mu = num(m[, 5]), sd = num(m[, 6])) %>%
  filter(!is.na(isco4)) %>% distinct(isco4, .keep_all = TRUE) %>%
  mutate(gradient = case_when(            # the ILO's own rules for its four exposure grades (WP140 Table 5)
    mu >= 0.6 & mu - sd >= 0.5  ~ "G4",
    mu >= 0.5 & mu + sd >= 0.5  ~ "G3",
    mu >= 0.4 & mu + sd >= 0.5  ~ "G2",
    mu <  0.4 & mu + sd >= 0.5  ~ "G1",
    TRUE ~ "None"))
write_csv(ilo, der("ilo_isco4_scores.csv"))

# 4. Link UK occupation codes (SOC 2020) to the ILO's ISCO-08 codes -----------------
# The ONS coding index gives an ISCO code for each job title. Where one UK occupation covers
# several ISCO codes, each gets a weight based on how many job titles it has
idx <- read_excel(raw("socindex.xlsx"), sheet = "SOC2020 coding index")
cw <- idx %>% transmute(soc4 = as.character(SOC_2020), isco4 = as.character(`ISCO-08 code based on SOC2020`)) %>%
  filter(str_detect(soc4, "^\\d{4}$"), str_detect(isco4, "^\\d{4}$")) %>%
  count(soc4, isco4, name = "n_titles") %>%
  group_by(soc4) %>% mutate(w = n_titles / sum(n_titles)) %>% ungroup()

# If an ISCO occupation has no ILO score, use the average of its wider group (3-digit, then 2-digit)
ilo3 <- ilo %>% group_by(isco3 = substr(isco4, 1, 3)) %>% summarise(mu3 = mean(mu), sd3 = mean(sd))
ilo2 <- ilo %>% group_by(isco2 = substr(isco4, 1, 2)) %>% summarise(mu2 = mean(mu), sd2 = mean(sd))
g <- read_csv(file.path("assumptions", "global.csv"), show_col_types = FALSE)
gp <- setNames(g$central, g$param)
# Turn each ILO grade into two scores: how much AI could replace the work (s) and how much it
# could help with it (c). The weights are set in assumptions/global.csv
grad_w <- function(gr) tibble(
  s = case_when(gr == "G4" ~ gp["s_w_G4"], gr == "G3" ~ gp["s_w_G3"], TRUE ~ 0),
  c = case_when(gr %in% c("G1") ~ gp["c_w_G1"], gr == "G2" ~ gp["c_w_G2"], gr == "G3" ~ gp["c_w_G3"], TRUE ~ 0))

cw <- cw %>% left_join(ilo %>% select(isco4, mu, sd, gradient), by = "isco4") %>%
  mutate(isco3 = substr(isco4, 1, 3), isco2 = substr(isco4, 1, 2)) %>%
  left_join(ilo3, by = "isco3") %>% left_join(ilo2, by = "isco2") %>%
  mutate(imputed = is.na(mu),
         sd = coalesce(sd, sd3, sd2), mu = coalesce(mu, mu3, mu2),
         gradient = if_else(is.na(gradient),
           case_when(mu >= 0.6 & mu - sd >= 0.5 ~ "G4", mu >= 0.5 & mu + sd >= 0.5 ~ "G3",
                     mu >= 0.4 & mu + sd >= 0.5 ~ "G2", mu < 0.4 & mu + sd >= 0.5 ~ "G1", TRUE ~ "None"),
           gradient))
cw <- bind_cols(cw, grad_w(cw$gradient))
# Weighted average score for each UK occupation
soc4_exp <- cw %>% group_by(soc4) %>%
  summarise(mu = sum(w * mu), s = sum(w * s), c = sum(w * c), share_imputed = sum(w * imputed), .groups = "drop")

# 5. How many young people work in each occupation group (APS, ONS ad hoc 3537) ------------------
# Blank (suppressed) cells are counted as zero and flagged
read_aps <- function(sh) {
  d <- read_excel(raw("occ_age.xlsx"), sheet = sh, col_names = FALSE, .name_repair = "minimal")[, 2:7]
  names(d) <- c("label", "a1624", "a2534", "a3549", "a50p", "total")
  d %>% filter(str_detect(label, "^\\d{3} ")) %>%
    transmute(soc3 = substr(label, 1, 3), title = str_trim(substring(label, 5)),
              y1624 = coalesce(num(a1624), 0), y_suppressed = is.na(num(a1624)), all = rowSums(cbind(num(a1624), num(a2534), num(a3549), num(a50p))
                                                              , na.rm = TRUE),
              year = as.integer(paste0("20", substr(sh, 3, 4))))
}
aps <- bind_rows(lapply(paste0("JD", 21:25), read_aps))
# Average up to the 104 occupation groups (SOC 3-digit) and weight each by its share of youth employment in 2025
soc3_exp <- soc4_exp %>% mutate(soc3 = substr(soc4, 1, 3)) %>%
  group_by(soc3) %>% summarise(across(c(mu, s, c, share_imputed), mean), .groups = "drop")
occ <- aps %>% filter(year == 2025) %>% left_join(soc3_exp, by = "soc3") %>%
  mutate(w_youth = y1624 / sum(y1624))
write_csv(occ, der("occupations_youth_exposure.csv"))
write_csv(aps %>% left_join(soc3_exp, by = "soc3"), der("aps_soc3_by_year.csv"))

cat("Stocks 2026Q2:\n"); print(stocks)
cat("Hazards:\n"); print(haz)
cat(sprintf("ILO units parsed: %d | SOC3 groups: %d | youth-weighted s = %.3f, c = %.3f | imputed share = %.2f\n",
            nrow(ilo), nrow(occ), sum(occ$w_youth * occ$s, na.rm = TRUE),
            sum(occ$w_youth * occ$c, na.rm = TRUE), sum(occ$w_youth * occ$share_imputed, na.rm = TRUE)))
