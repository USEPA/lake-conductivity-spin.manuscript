library(tidyverse)
library(xtable)

dat <- read_csv("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/pc1_documents/publications/spatialstat2023.si/output/cv_model_comp.csv")

dat <- dat[-41, -c(3, 5, 9, 13)]

dat <- dat |>
  mutate(formula = case_when(
    formula == "intercept" ~ "No",
    formula == "full" ~ "Yes"
  )) |>
  mutate(spcov_type = case_when(
    spcov_type == "exponential" ~ "Exp",
    spcov_type == "gaussian" ~ "Gau",
    spcov_type == "none" ~ "None"
  )) |>
  mutate(anisotropy = case_when(
    anisotropy ~ "Yes",
    !anisotropy ~ "No"
  )) |>
  mutate(partition_factor = case_when(
    partition_factor == "NULL" ~ "No",
    partition_factor == "~DSGN_CYCLE" ~ "Yes"
  ))

dat <- dat |>
  mutate(
    bias = format(round(bias, digits = 3)),
    mspe = str_c(round(mspe, digits = 3), " (", rank(mspe), ")"),
    r2 = str_c(round(r2, digits = 3), " (", rank(desc(r2)), ")"),
    cover = round(cover, digits = 3),
    time = round(time, digits = 2)
  )

dat <- dat |>
  select(formula, spcov_type, anisotropy, partition_factor, index, bias, mspe, r2, cover, time)

print(xtable(dat), include.rownames = FALSE)
