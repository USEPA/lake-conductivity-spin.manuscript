################################################################################
######### Load relevant software
################################################################################

library(tidyverse)
library(spmodel)
library(sf)
library(parallel)
library(here)

################################################################################
######### Set paths and read in data
################################################################################


nla <- read_csv(here("data", "nla_obs.csv"), guess_max = Inf)
outpath <- here("output")


################################################################################
######### Data cleaning
################################################################################

nla <- nla |>
  st_as_sf(coords = c("XCOORD", "YCOORD"), crs = "+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")

# epsg https://epsg.io/102003

nla <- nla |>
  mutate(
    DSGN_CYCLE = as.factor(DSGN_CYCLE),
    BinCropWs = if_else(PctCropWs > 0, 1, 0),
    BinUrb = if_else(PctUrb > 0, 1, 0)
  )

################################################################################
######### Finding local index separately for each year using kmeans
################################################################################

set.seed(10)
nfold <- 10
nla$fold <- sample(rep(seq_len(nfold), length.out = NROW(nla)))

################################################################################
######### Some visualizations
################################################################################

ggplot(nla, aes(x = DSGN_CYCLE, y = log(COND_RESULT))) +
  geom_point()

ggplot(nla, aes(x = Tmean8110Cat, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = Precip8110Ws, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = AREA_HA, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = DSGN_CYCLE, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = PctCropWs, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = PctUrb, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = CaOWs, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

ggplot(nla, aes(x = SWs, y = log(COND_RESULT))) +
  geom_point() + facet_grid(~ DSGN_CYCLE)

################################################################################
######### Creating CV functions
################################################################################


find_cv_all <- function(trial, arglist, fold, data) {

  set.seed(10)

  find_cv_ind <- function(arglist, fold, data) {

    arglist <- lapply(arglist, function(x) x[[1]])
    fold_vals <- unique(fold)
    cv_out <- lapply(fold_vals, function(x) {
      holdout <- fold == x
      find_cv_loss(holdout, arglist, data)
    })
    mspe <- mean(unlist(lapply(cv_out, function(x) x$spe)))
    mae <- median(unlist(lapply(cv_out, function(x) x$ae)))
    bias <- mean(unlist(lapply(cv_out, function(x) x$bias)))
    cor <- mean(unlist(lapply(cv_out, function(x) x$cor)))
    cover <- mean(unlist(lapply(cv_out, function(x) x$cover)))
    na_length <- sum(unlist(lapply(cv_out, function(x) x$na_length)))
    list(mspe = mspe, mae = mae, bias = bias, cor = cor, cover = cover, na_length = na_length)
  }

  find_cv_loss <- function(holdout, arglist, data) {

    data_train <- data[!holdout, , ]
    data_test <- data[holdout, , ]
    spmod <- splm(
      formula = arglist$formula,
      data = data[!holdout, , ],
      spcov_type = as.character(arglist$spcov_type),
      estmethod = as.character(arglist$estmethod),
      random = arglist$random,
      partition_factor = arglist$partition_factor,
      anisotropy = arglist$anisotropy,
      local = list(method = arglist$index, groups = 30)
      # local = list(index = arglist$index_val[!holdout])
    )

    cv_est_all <- predict(spmod, newdata = data_test, interval = "prediction")
    cv_est <- cv_est_all[, "fit", drop = TRUE]
    true_val <- log(data_test$COND_RESULT)
    # handle NA
    is_na_val <- is.na(cv_est)
    cv_est <- cv_est[!is_na_val]
    true_val <- true_val[!is_na_val]
    # continue
    error <- cv_est - true_val
    # coverage
    lwr <- cv_est_all[, "lwr", drop = TRUE]
    upr <- cv_est_all[, "upr", drop = TRUE]
    cover <- lwr <= true_val & true_val <= upr
    list(spe = error^2, ae = abs(error), bias = error, cor = cor(cv_est, true_val), cover = mean(cover), na_length = sum(is_na_val))
  }

  arglist <- arglist[trial, ]
  find_cv_ind(arglist, fold, data)
}

fit_all <- function(trial, arglist, data) {

  set.seed(10)

  arglist <- arglist[trial, ]
  arglist <- lapply(arglist, function(x) x[[1]])
  start_time <- proc.time()
  spmod <- splm(
    formula = arglist$formula,
    data = data,
    spcov_type = as.character(arglist$spcov_type),
    estmethod = as.character(arglist$estmethod),
    random = arglist$random,
    partition_factor = arglist$partition_factor,
    anisotropy = arglist$anisotropy,
    local = list(method = as.character(arglist$index), groups = 30)
  )
  end_time <- (proc.time() - start_time)["elapsed"]
  end_time
}

find_cv_time_spin <- function(trial, arglist, fold, data) {

  set.seed(10)

  find_cv_ind <- function(arglist, fold, data) {

    arglist <- lapply(arglist, function(x) x[[1]])
    fold_vals <- unique(fold)
    start_time <- proc.time()
    cv_out <- lapply(fold_vals, function(x) {
      holdout <- fold == x
      find_cv_loss(holdout, arglist, data)
    })
    end_time <- (proc.time() - start_time)["elapsed"]
    mspe <- mean(unlist(lapply(cv_out, function(x) x$spe)))
    mae <- median(unlist(lapply(cv_out, function(x) x$ae)))
    bias <- mean(unlist(lapply(cv_out, function(x) x$bias)))
    cor <- mean(unlist(lapply(cv_out, function(x) x$cor)))
    cover <- mean(unlist(lapply(cv_out, function(x) x$cover)))
    na_length <- sum(unlist(lapply(cv_out, function(x) x$na_length)))
    list(mspe = mspe, mae = mae, bias = bias, cor = cor, cover = cover, na_length = na_length, time = end_time)
  }

  find_cv_loss <- function(holdout, arglist, data) {

    data_train <- data[!holdout, , ]
    data_test <- data[holdout, , ]
    spmod <- splm(
      formula = arglist$formula,
      data = data[!holdout, , ],
      spcov_type = as.character(arglist$spcov_type),
      estmethod = as.character(arglist$estmethod),
      random = arglist$random,
      partition_factor = arglist$partition_factor,
      anisotropy = arglist$anisotropy,
      local = list(method = arglist$index, groups = 30)
    )

    cv_est_all <- predict(spmod, newdata = data_test, interval = "prediction")
    cv_est <- cv_est_all[, "fit", drop = TRUE]
    true_val <- log(data_test$COND_RESULT)
    # handle NA
    is_na_val <- is.na(cv_est)
    cv_est <- cv_est[!is_na_val]
    true_val <- true_val[!is_na_val]
    # continue
    error <- cv_est - true_val
    # coverage
    lwr <- cv_est_all[, "lwr", drop = TRUE]
    upr <- cv_est_all[, "upr", drop = TRUE]
    cover <- lwr <= true_val & true_val <= upr
    list(spe = error^2, ae = abs(error), bias = error, cor = cor(cv_est, true_val), cover = mean(cover), na_length = sum(is_na_val))
  }

  arglist <- arglist[trial, ]
  find_cv_ind(arglist, fold, data)
}

find_cv_time_nospin <- function(trial, arglist, fold, data) {

  set.seed(10)

  find_cv_ind <- function(arglist, fold, data) {

    arglist <- lapply(arglist, function(x) x[[1]])
    fold_vals <- unique(fold)
    start_time <- proc.time()
    cv_out <- lapply(fold_vals, function(x) {
      holdout <- fold == x
      find_cv_loss(holdout, arglist, data)
    })
    end_time <- (proc.time() - start_time)["elapsed"]
    mspe <- mean(unlist(lapply(cv_out, function(x) x$spe)))
    mae <- median(unlist(lapply(cv_out, function(x) x$ae)))
    bias <- mean(unlist(lapply(cv_out, function(x) x$bias)))
    cor <- mean(unlist(lapply(cv_out, function(x) x$cor)))
    cover <- mean(unlist(lapply(cv_out, function(x) x$cover)))
    na_length <- sum(unlist(lapply(cv_out, function(x) x$na_length)))
    list(mspe = mspe, mae = mae, bias = bias, cor = cor, cover = cover, na_length = na_length, time = end_time)
  }

  find_cv_loss <- function(holdout, arglist, data) {

    data_train <- data[!holdout, , ]
    data_test <- data[holdout, , ]
    spmod <- splm(
      formula = arglist$formula,
      data = data[!holdout, , ],
      spcov_type = as.character(arglist$spcov_type),
      estmethod = as.character(arglist$estmethod),
      random = arglist$random,
      partition_factor = arglist$partition_factor,
      anisotropy = arglist$anisotropy,
      local = FALSE
    )

    cv_est_all <- predict(spmod, newdata = data_test, interval = "prediction")
    cv_est <- cv_est_all[, "fit", drop = TRUE]
    true_val <- log(data_test$COND_RESULT)
    # handle NA
    is_na_val <- is.na(cv_est)
    cv_est <- cv_est[!is_na_val]
    true_val <- true_val[!is_na_val]
    # continue
    error <- cv_est - true_val
    # coverage
    lwr <- cv_est_all[, "lwr", drop = TRUE]
    upr <- cv_est_all[, "upr", drop = TRUE]
    cover <- lwr <= true_val & true_val <= upr
    list(spe = error^2, ae = abs(error), bias = error, cor = cor(cv_est, true_val), cover = mean(cover), na_length = sum(is_na_val))
  }

  arglist <- arglist[trial, ]
  find_cv_ind(arglist, fold, data)
}

write_preds <- function(nla_year, object, outpath, spin = TRUE) {


  nla_preds <- read_csv(here("data", paste0("nla_pred", nla_year, ".csv")), guess_max = Inf)

  nla_preds <- nla_preds %>% 
    mutate(
      BinCropWs = if_else(PctCropWs > 0, 1, 0),
      BinUrb = if_else(PctUrb > 0, 1, 0),
      DSGN_CYCLE = factor(nla_year),
      UNIQUE_ID = "1"
    ) |>
    select(AREA_HA, DSGN_CYCLE, Tmean8110Cat, Precip8110Ws, BinCropWs, PctCropWs,
           BinUrb, PctUrb, CaOWs, SWs, UNIQUE_ID, Longitude, Latitude) %>%
    filter(complete.cases(.)) |>
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) |>
    st_transform(crs = "+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")

  set.seed(10)
  start_time <- proc.time()
  object_preds <- augment(object, newdata = nla_preds, se_fit = TRUE, local = list(method = "all", parallel = TRUE))
  end_time <- (proc.time() - start_time)["elapsed"]

  object_preds_coords <- st_coordinates(object_preds)
  object_preds_df <- st_drop_geometry(object_preds)
  object_preds_df$XCOORD <- object_preds_coords[, 1]
  object_preds_df$YCOORD <- object_preds_coords[, 2]

  if (spin) {
    write_csv(object_preds_df, paste0(outpath, "preds_spin_", nla_year, ".csv"))
  } else {
    write_csv(object_preds_df, paste0(outpath, "preds_nospin_", nla_year, ".csv"))
  }
  data.frame(time = end_time, n_pred = NROW(nla_preds), nla_year = nla_year)
}

################################################################################
######### Model Comparison
################################################################################

################################################################################
######### Creating the Grid
################################################################################

formula_int <- log(COND_RESULT) ~ 1
formula_vars <- log(COND_RESULT) ~ AREA_HA + DSGN_CYCLE + Tmean8110Cat +
  Precip8110Ws + BinCropWs + BinCropWs:PctCropWs + BinUrb + BinUrb:PctUrb +
  CaOWs + SWs
# formula_vars <- log(COND_RESULT) ~ Tmean8110Cat +
#   Precip8110Ws + BinCropWs + BinCropWs:PctCropWs + BinUrb + BinUrb:PctUrb +
#   I(CaOWs + SWs)


formula <- list(formula_int, formula_vars)
spcov_type <- c("exponential", "gaussian", "none")
estmethod <- "reml"
anisotropy <- c(TRUE, FALSE)
# random <- list(~ UNIQUE_ID + DSGN_CYCLE)
random <- list(~ UNIQUE_ID)
part_none <- NULL
part <- ~ DSGN_CYCLE
partition_factor <- list(part_none, part)
index <- c("kmeans", "random")

cv_grid <- expand.grid(
  formula = formula,
  spcov_type = spcov_type,
  estmethod = estmethod,
  anisotropy = anisotropy,
  random = random,
  partition_factor = partition_factor,
  index = index
)

which_bad <- which(cv_grid$spcov_type == "none" & cv_grid$anisotropy)
cv_grid <- cv_grid[-which_bad, , drop = FALSE]
print(cv_grid)

numCores <- pmin(detectCores(), NROW(cv_grid))
cl <- makeCluster(numCores)
clusterEvalQ(cl, library(spmodel))
cv_vals <- parLapply(cl, seq_len(NROW(cv_grid)), safely(find_cv_all), cv_grid, nla$fold, nla)
stopCluster(cl)

cv_grid$mspe <- map_dbl(cv_vals, ~ .x$result$mspe)
cv_grid$mae <- map_dbl(cv_vals, ~ .x$result$mae)
cv_grid$bias <- map_dbl(cv_vals, ~ .x$result$bias)
cv_grid$r2 <- map_dbl(cv_vals, ~ .x$result$cor^2)
cv_grid$cover <- map_dbl(cv_vals, ~ .x$result$cover)
cv_grid$na_length <- map_dbl(cv_vals, ~ .x$result$na_length)
print(cv_grid)

# print_cv_grid <- cv_grid
# print_cv_grid$formula <- rep(c("intercept", "full"), length.out = NROW(cv_grid))
# print_cv_grid$random <- as.character(print_cv_grid$random)
# print_cv_grid$partition_factor <- as.character(print_cv_grid$partition_factor)
# write_csv(print_cv_grid, paste0(outpath, "cv_model_comp.csv"))

################################################################################
######### Timing All Models
################################################################################

fit_time <- lapply(seq_len(NROW(cv_grid)), safely(fit_all), cv_grid, nla)
cv_grid$time <- map_dbl(fit_time, ~ .x$result)

best_mspe <- which.min(cv_grid$mspe)
# add another row as best grid
cv_grid[NROW(cv_grid) + 1, ] <- cv_grid[best_mspe, ]
cv_grid[NROW(cv_grid), ]
cv_grid$index[NROW(cv_grid)] <- NA
cv_grid$mspe[NROW(cv_grid)] <- NA
cv_grid$mae[NROW(cv_grid)] <- NA
cv_grid$bias[NROW(cv_grid)] <- NA
cv_grid$r2[NROW(cv_grid)] <- NA
cv_grid$cover[NROW(cv_grid)] <- NA
cv_grid$na_length[NROW(cv_grid)] <- NA
start_time <- proc.time()
spmod_nospin <- splm(
  log(COND_RESULT) ~ AREA_HA + DSGN_CYCLE + Tmean8110Cat +
    Precip8110Ws + BinCropWs + BinCropWs:PctCropWs + BinUrb + BinUrb:PctUrb +
    CaOWs + SWs,
  data = nla,
  spcov_type = "exponential",
  estmethod = "reml",
  random = ~ UNIQUE_ID,
  partition_factor = NULL,
  anisotropy = TRUE,
  local = FALSE
)
end_time <- (proc.time() - start_time)["elapsed"]
cv_grid$time[NROW(cv_grid)] <- end_time

print_cv_grid <- cv_grid
print_cv_grid$formula <- rep(c("intercept", "full"), length.out = NROW(cv_grid))
print_cv_grid$random <- as.character(print_cv_grid$random)
print_cv_grid$partition_factor <- as.character(print_cv_grid$partition_factor)
write_csv(print_cv_grid, paste0(outpath, "cv_model_comp.csv"))

################################################################################
######### Compare Spatial Indexing to Full -- CV
################################################################################

################################################################################
######### Creating the Grid
################################################################################
best_grid <- cv_grid[which.min(cv_grid$mspe), ]
best_grid


cv_vals_spin <- lapply(1, safely(find_cv_time_spin), best_grid, nla$fold, nla)
best_grid_spin <- best_grid
best_grid_spin$mspe <- map_dbl(cv_vals_spin, ~ .x$result$mspe)
best_grid_spin$mae <- map_dbl(cv_vals_spin, ~ .x$result$mae)
best_grid_spin$bias <- map_dbl(cv_vals_spin, ~ .x$result$bias)
best_grid_spin$r2 <- map_dbl(cv_vals_spin, ~ .x$result$cor^2)
best_grid_spin$cover <- map_dbl(cv_vals_spin, ~ .x$result$cover)
best_grid_spin$time <- map_dbl(cv_vals_spin, ~ .x$result$time)
best_grid_spin$na_length <- map_dbl(cv_vals_spin, ~ .x$result$na_length)


cv_vals_nospin <- lapply(1, safely(find_cv_time_nospin), best_grid, nla$fold, nla)
best_grid_nospin <- best_grid
best_grid_nospin$mspe <- map_dbl(cv_vals_nospin, ~ .x$result$mspe)
best_grid_nospin$mae <- map_dbl(cv_vals_nospin, ~ .x$result$mae)
best_grid_nospin$bias <- map_dbl(cv_vals_nospin, ~ .x$result$bias)
best_grid_nospin$r2 <- map_dbl(cv_vals_nospin, ~ .x$result$cor^2)
best_grid_nospin$cover <- map_dbl(cv_vals_nospin, ~ .x$result$cover)
best_grid_nospin$time <- map_dbl(cv_vals_nospin, ~ .x$result$time)
best_grid_nospin$na_length <- map_dbl(cv_vals_nospin, ~ .x$result$na_length)

best_grid_comp <- rbind(best_grid_spin, best_grid_nospin)

print_best_grid <- best_grid_comp
print_best_grid$formula <- "full"
print_best_grid$random <- as.character(print_best_grid$random)
print_best_grid$partition_factor <- as.character(print_best_grid$partition_factor)
print_best_grid$index <- c("kmeans", "none")
write_csv(print_best_grid, paste0(outpath, "best_grid_comp.csv"))

################################################################################
######### Compare Spatial Indexing to Full -- All Data Fit
################################################################################

set.seed(10)
spmod_spin <- splm(
  log(COND_RESULT) ~ AREA_HA + DSGN_CYCLE + Tmean8110Cat +
    Precip8110Ws + BinCropWs + BinCropWs:PctCropWs + BinUrb + BinUrb:PctUrb +
    CaOWs + SWs,
  data = nla,
  spcov_type = "exponential",
  estmethod = "reml",
  random = ~ UNIQUE_ID,
  partition_factor = NULL,
  anisotropy = TRUE,
  local = list(method = "kmeans", groups = 30)
)

# fixed effects out
fe_spin <- tidy(spmod_spin) |>
  mutate(spin = TRUE)
fe_nospin <- tidy(spmod_nospin) |>
  mutate(spin = FALSE)
fe_out <- bind_rows(fe_spin, fe_nospin) |>
  arrange(term)
write_csv(fe_out, paste0(outpath, "fe_out.csv"))

# spatial covariance out
sp_spin <- tidy(spmod_spin, effects = "spcov") |>
  mutate(spin = TRUE)
sp_nospin <- tidy(spmod_nospin, effects = "spcov") |>
  mutate(spin = FALSE)
sp_out <- bind_rows(sp_spin, sp_nospin) |>
  arrange(term)
write_csv(sp_out, paste0(outpath, "sp_out.csv"))

# random covariance out
rand_spin <- tidy(spmod_spin, effects = "randcov") |>
  mutate(spin = TRUE)
rand_nospin <- tidy(spmod_nospin, effects = "randcov") |>
  mutate(spin = FALSE)
rand_out <- bind_rows(rand_spin, rand_nospin) |>
  arrange(term)
write_csv(rand_out, paste0(outpath, "rand_out.csv"))

# varcomp out
varcomp_spin <- varcomp(spmod_spin) |>
  mutate(spin = TRUE)
varcomp_nospin <- varcomp(spmod_nospin) |>
  mutate(spin = FALSE)
varcomp_out <- bind_rows(varcomp_spin, varcomp_nospin) |>
  arrange(varcomp)
write_csv(varcomp_out, paste0(outpath, "varcomp_out.csv"))

# augment out
spmod_spin_aug <- augment(spmod_spin)
spmod_spin_coords <- st_coordinates(spmod_spin_aug)
spmod_spin_aug <- st_drop_geometry(spmod_spin_aug)
spmod_spin_aug$XCOORD <- spmod_spin_coords[, 1]
spmod_spin_aug$YCOORD <- spmod_spin_coords[, 2]
names(spmod_spin_aug)
names(spmod_spin_aug)[1] <- "LogCn"
write_csv(spmod_spin_aug, paste0(outpath, "spmod_spin_aug.csv"))

################################################################################
######### Prediction
################################################################################

# calling for each year (SPIN)
pred_2007 <- write_preds(2007, spmod_spin, outpath, spin = TRUE)
pred_2012 <- write_preds(2012, spmod_spin, outpath, spin = TRUE)
pred_2017 <- write_preds(2017, spmod_spin, outpath, spin = TRUE)

pred_out <- bind_rows(pred_2007, pred_2012, pred_2017)
write_csv(pred_out, paste0(outpath, "pred_spin_out.csv"))

# calling for each year (no SPIN)
pred_nospin_2007 <- write_preds(2007, spmod_nospin, outpath, spin = FALSE)
pred_nospin_2012 <- write_preds(2012, spmod_nospin, outpath, spin = FALSE)
pred_nospin_2017 <- write_preds(2017, spmod_nospin, outpath, spin = FALSE)

pred_out <- bind_rows(pred_nospin_2007, pred_nospin_2012, pred_nospin_2017)
write_csv(pred_out, paste0(outpath, "pred_nospin_out.csv"))

################################################################################
######### Session Information
################################################################################

sessionInfo()
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
#
# Matrix products: default
#
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8
#
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] data.table_1.14.8 sf_1.0-13         spmodel_0.4.1     lubridate_1.9.2   forcats_1.0.0     stringr_1.5.0     dplyr_1.1.2
# [8] purrr_1.0.1       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.2     tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10        compiler_4.2.2     pillar_1.9.0       class_7.3-20       tools_4.2.2        bit_4.0.5          viridisLite_0.4.2
# [8] timechange_0.2.0   lifecycle_1.0.3    gtable_0.3.3       lattice_0.20-45    pkgconfig_2.0.3    rlang_1.1.1        Matrix_1.6-1
# [15] DBI_1.1.3          cli_3.6.1          rstudioapi_0.14    e1071_1.7-13       withr_2.5.0        hms_1.1.3          generics_0.1.3
# [22] vctrs_0.6.2        bit64_4.0.5        classInt_0.4-9     grid_4.2.2         tidyselect_1.2.0   glue_1.6.2         R6_2.5.1
# [29] fansi_1.0.4        vroom_1.6.3        farver_2.1.1       tzdb_0.4.0         magrittr_2.0.3     scales_1.2.1       units_0.8-2
# [36] colorspace_2.1-0   labeling_0.4.2     KernSmooth_2.23-20 utf8_1.2.3         stringi_1.7.12     proxy_0.4-27       munsell_0.5.0
# [43] crayon_1.5.2
