# 00_extract_a06.R ----------------------------------------------------------
# Takes the 16-24 year old series we need out of ONS table A06 (seasonally adjusted).
# Each ONS series has a short code (its CDID), so we pick the columns using those codes.
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(readr) })
a <- read_excel(file.path("data", "raw", "a06.xls"), sheet = "People", col_names = FALSE, .name_repair = "minimal")
codes <- as.character(unlist(a[8, ]))
want <- c("JN5R", "JN62", "JN6B", "AGNJ", "AGOL", "AGPM")
out <- a[10:nrow(a), c(1, match(want, codes))]
names(out) <- c("period", want)
out <- out %>% filter(!is.na(JN5R), grepl("^[A-Z][a-z]{2}-[A-Z][a-z]{2} \\d{4}", period)) %>%
  mutate(across(-period, as.numeric))
dir.create(file.path("data", "derived"), showWarnings = FALSE, recursive = TRUE)
write_csv(out, file.path("data", "derived", "a06_youth.csv"))
