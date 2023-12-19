################################################################################
######### Load relevant software
################################################################################

library(tidyverse)
library(spmodel)
library(sf)
library(knitr)

path <- str_c(here(), "/")

LogCn_limits <- c(0, 11.1)

################################################################################
######### Figure
################################################################################

nla <- read_csv(str_c(path, "data/nla_obs.csv")) |>
  mutate(LogCn = log(COND_RESULT))

nla_2007 <- nla |>
  filter(DSGN_CYCLE == "2007")

nla_2012 <- nla |>
  filter(DSGN_CYCLE == "2012")

nla_2017 <- nla |>
  filter(DSGN_CYCLE == "2017")

######################################
######### Subfigure
######################################

logcn_hist_2007 <- ggplot(nla_2007, aes(x = LogCn)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  labs(x = "Log Conductivity (LogCn)", y = "Frequency", title = "2007") +
  lims(x = LogCn_limits) +
  theme_bw(base_size = 18) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_hist_2007.jpeg"), logcn_hist_2007, dpi = 300,
       width = 5.21, height = 5.27)

######################################
######### Subfigure
######################################

logcn_map_2007 <- ggplot(nla_2007, aes(x = XCOORD, y = YCOORD, color = LogCn)) +
  geom_point(size = 1.5) +
  labs(title = "2007") +
  scale_color_viridis_c(limits = LogCn_limits) +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_map_2007.jpeg"), logcn_map_2007, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/logcn_map_2007.jpeg"))

######################################
######### Subfigure
######################################

logcn_hist_2012 <- ggplot(nla_2012, aes(x = LogCn)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  labs(x = "Log Conductivity (LogCn)", y = "Frequency", title = "2012") +
  lims(x = LogCn_limits) +
  theme_bw(base_size = 18) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_hist_2012.jpeg"), logcn_hist_2012, dpi = 300,
       width = 5.21, height = 5.27)

######################################
######### Subfigure
######################################

logcn_map_2012 <- ggplot(nla_2012, aes(x = XCOORD, y = YCOORD, color = LogCn)) +
  geom_point(size = 1.5) +
  labs(title = "2012") +
  scale_color_viridis_c(limits = LogCn_limits) +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_map_2012.jpeg"), logcn_map_2012, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/logcn_map_2012.jpeg"))

######################################
######### Subfigure
######################################

logcn_hist_2017 <- ggplot(nla_2017, aes(x = LogCn)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  labs(x = "Log Conductivity (LogCn)", y = "Frequency", title = "2017") +
  lims(x = LogCn_limits) +
  theme_bw(base_size = 18) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_hist_2017.jpeg"), logcn_hist_2017, dpi = 300,
       width = 5.21, height = 5.27)

######################################
######### Subfigure
######################################

logcn_map_2017 <- ggplot(nla_2017, aes(x = XCOORD, y = YCOORD, color = LogCn)) +
  geom_point(size = 1.5) +
  labs(title = "2017") +
  scale_color_viridis_c(limits = LogCn_limits) +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/logcn_map_2017.jpeg"), logcn_map_2017, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/logcn_map_2017.jpeg"))


################################################################################
######### Figure
################################################################################

pred_2007 <- read_csv(str_c(path, "preds_spin_2007.csv")) |>
  mutate(
    .fitted = if_else(.fitted > 11.1, 11.1, .fitted),
    .fitted = if_else(.fitted < 0, 0, .fitted),
    .se.fit = if_else(.se.fit > 1.1, 1.1, .se.fit)
  )

pred_2012 <- read_csv(str_c(path, "preds_spin_2012.csv")) |>
  mutate(
    .fitted = if_else(.fitted > 11.1, 11.1, .fitted),
    .fitted = if_else(.fitted < 0, 0, .fitted),
    .se.fit = if_else(.se.fit > 1.1, 1.1, .se.fit)
  )

pred_2017 <- read_csv(str_c(path, "preds_spin_2017.csv")) |>
  mutate(
    .fitted = if_else(.fitted > 11.1, 11.1, .fitted),
    .fitted = if_else(.fitted < 0, 0, .fitted),
    .se.fit = if_else(.se.fit > 1.1, 1.1, .se.fit)
  )


######################################
######### Subfigure
######################################

preds_map_2007 <- ggplot(pred_2007, aes(x = XCOORD, y = YCOORD, color = .fitted)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(limits = LogCn_limits) +
  labs(color = "Preds ", title = "2007") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/preds_map_2007.jpeg"), preds_map_2007, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/preds_map_2007.jpeg"))

######################################
######### Subfigure
######################################

se_map_2007 <- ggplot(pred_2007, aes(x = XCOORD, y = YCOORD, color = .se.fit)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(option = "A") +
  labs(color = "StdErs", title = "2007") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/se_map_2007.jpeg"), se_map_2007, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/se_map_2007.jpeg"))

######################################
######### Subfigure
######################################

preds_map_2012 <- ggplot(pred_2012, aes(x = XCOORD, y = YCOORD, color = .fitted)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(limits = LogCn_limits) +
  labs(color = "Preds ", title = "2012") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/preds_map_2012.jpeg"), preds_map_2012, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/preds_map_2012.jpeg"))

######################################
######### Subfigure
######################################

se_map_2012 <- ggplot(pred_2012, aes(x = XCOORD, y = YCOORD, color = .se.fit)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(option = "A") +
  labs(color = "StdErs", title = "2012") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/se_map_2012.jpeg"), se_map_2012, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/se_map_2012.jpeg"))



######################################
######### Subfigure
######################################

preds_map_2017 <- ggplot(pred_2017, aes(x = XCOORD, y = YCOORD, color = .fitted)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(limits = LogCn_limits) +
  labs(color = "Preds ", title = "2017") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/preds_map_2017.jpeg"), preds_map_2017, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/preds_map_2017.jpeg"))

######################################
######### Subfigure
######################################

se_map_2017 <- ggplot(pred_2017, aes(x = XCOORD, y = YCOORD, color = .se.fit)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(option = "A") +
  labs(color = "StdErs", title = "2017") +
  coord_fixed() +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(str_c(path, "figures/se_map_2017.jpeg"), se_map_2017, dpi = 300,
       width = 5.21, height = 5.27)

plot_crop(str_c(path, "figures/se_map_2017.jpeg"))

################################################################################
######### Figure
################################################################################

######################################
######### Subfigure
######################################

sp_out <- read_csv(str_c(path, "sp_out.csv"))
sp_out
x_coords <- c(1, sqrt(2) / 2, 0, -sqrt(2)/2)
y_coords <- c(0, sqrt(2) / 2, 1, sqrt(2) / 2)
coords <- rbind(x_coords, y_coords)
coords
rotate <- 1.51
scale <- 0.325
rotate_clockwise_inv <- matrix(c(cos(rotate), -sin(rotate), sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
scale_yaxis_inv <- matrix(c(1, 0, 0, scale), nrow = 2, ncol = 2, byrow = TRUE)
new_coords <- (rotate_clockwise_inv %*% scale_yaxis_inv) %*% coords
new_coords
# find distances
sqrt(new_coords[1, ]^2 + new_coords[2, ]^2) * 327303

x_seq <- seq(-1, 1, length.out = 100)
y_seq <- seq(- 1, 1, length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
coords <- t(as.matrix(grid))
rotate <- 1.51
scale <- .325
range <- 1
rotate_clockwise_inv <- matrix(c(cos(rotate), -sin(rotate), sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
scale_yaxis_inv <- matrix(c(1, 0, 0, scale), nrow = 2, ncol = 2, byrow = TRUE)
new_coords <- (rotate_clockwise_inv %*% scale_yaxis_inv) %*% coords
grid$x_new <- new_coords[1, ]
grid$y_new <- new_coords[2, ]
dist <- sqrt(grid$x_new^2 + grid$y_new^2)
grid$cov <- exp(-dist/range)
ggplot(grid, aes(x = x_new, y = y_new, fill = cov)) +
  geom_tile() +
  scale_fill_viridis_c()

get_line <- function(range_frac, dir) {
  tibble::tibble(
    Distance = seq(0, 1, length.out = 1000),
    Correlation = exp(-Distance / range_frac),
    Direction = dir
  )
}


dat <- bind_rows(map2(c(1, 0.743, 0.325), c("N-S", "NE-SW", "E-W"), get_line))
anis1 <- ggplot(dat, aes(x = Distance, y = Correlation, color = Direction)) +
   geom_line(size = 2) +
   scale_color_viridis_d(end = 0.9) +
   theme_bw(base_size = 18) +
   scale_x_continuous(breaks = c(0, 1), labels = c("0", "327.30")) +
   labs(x = "Distance (km)")

ggsave(str_c(path, "figures/anis1.jpeg"), anis1, dpi = 300,
       width = 5.21, height = 5.27)


######################################
######### Subfigure
######################################
range <- 327303
get_ellipse <- function(r_frac) {

  range <- 327303
  r <- range * r_frac
  theta_seq <- seq(0, 2 * pi, length.out = 10000)
  x_orig <- r * cos(theta_seq)
  y_orig <- r * sin(theta_seq)
  coords <- rbind(x_orig, y_orig)
  rotate <- 1.51
  scale <- 0.325

  rotate_clockwise_inv <- matrix(c(cos(rotate), -sin(rotate), sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis_inv <- matrix(c(1, 0, 0, scale), nrow = 2, ncol = 2, byrow = TRUE)
  new_coords <- (rotate_clockwise_inv %*% scale_yaxis_inv) %*% coords
  x_new <- new_coords[1, ]
  y_new <- new_coords[2, ]
  dat <- data.frame(x = x_new, y = y_new, cor = round(exp(-r / range), 3))

}

anis_seq <- c(0.2, 0.6, 1)
ellipse_df <- bind_rows(map(anis_seq, get_ellipse))
anis2 <- ggplot(ellipse_df , aes(x = x, y = y, color = cor)) +
  geom_point(size = 1) +
  scale_color_viridis_c(name = "Cor", end = 0.9, option = "A") +
  labs(x = "Easting (km)", y = "Northing (km)") +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(-range, range), breaks = c(-range, 0, range), labels = c("327.30", "0", "327.30")) +
  scale_y_continuous(limits = c(-range, range), breaks = c(-range, 0, range), labels = c("327.30", "0", "327.30")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

ggsave(str_c(path, "figures/anis2.jpeg"), anis2, dpi = 300,
       width = 5.21, height = 5.27)


