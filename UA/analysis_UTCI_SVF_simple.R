#' Relationship between UTCI and sky view factor (SVF) - simplified analysis
#' Routes forum + park (10 sites), Berlin-Adlershof, 2026-05-29
#' Author: Adelina Petkova
#' Date: 09.2026
#'
#' Questions
#'   H1: Is UTCI related to SVF?
#'   H2: Does this relationship change over the day?
#'
#' Approach
#'   * All sites share the same weather. For every hour we therefore compute
#'     the UTCI anomaly: how much warmer (+) or cooler (-) a site was than
#'     the mean of all sites measured in that hour.
#'   * SVF has one value per site, so the site (n = 10) is the unit of the
#'     simple regressions.
#'
#'   1  descriptives
#'   2  H1: correlation + linear regression of site-mean UTCI anomaly on SVF
#'   3  H2: (a) the same regression separately for morning, midday, afternoon
#'          (b) linear mixed model with an SVF x period interaction, using
#'              all site-hours, site as random effect and AR(1) errors
#'          (c) descriptive figures: hourly SVF slopes and a heatmap
#'
#' Outputs are written to output_simple/ (csv tables + png figures).

#### settings ####
library(ggplot2)
library(dplyr)
library(lubridate)
library(nlme)   # mixed models (ships with R)

out_dir <- "output_simple"
dir.create(out_dir, showWarnings = FALSE)

tz_local <- "Europe/Berlin"        # UTC+2 (CEST) during the campaign

# Globe diameter [m] used for Tmrt. analysis_UA.R uses 0.15 m (ISO 7726
# standard globe). The Kestrel 5400 globe is a 1-inch (0.0254 m) globe, so
# check which value is correct for your instrument - it changes Tmrt a lot.
globe_D <- 0.0254
globe_eps <- 0.95

# UTCI is only defined for 10-m wind speeds >= 0.5 m/s. The Kestrel reports
# 0 m/s below its stall speed, so the 10-m wind is clipped at 0.5 m/s.
clip_wind_utci <- TRUE

# The campaign targeted daytime conditions; the few site-hours after sunset
# (21:13 CEST) are dropped. FALSE keeps all hours.
exclude_night <- TRUE

#### 0 data preparation ####
# Run the existing preprocessing scripts in separate environments so that
# they do not overwrite each other's objects.
env_park <- new.env(); env_forum <- new.env()
sys.source("analysis_UA.R",    envir = env_park)
sys.source("analysis_forum.R", envir = env_forum)
source("calc_UTCI.R")

# exact (mean) measurement time of every site-hour, from the filtered 10-s data
mean_times <- function(env) {
  env$filtered_data %>%
    mutate(hour = floor_date(Time, "hour")) %>%
    group_by(site.name, hour) %>%
    summarise(time_mid = mean(Time), .groups = "drop")
}

d <- bind_rows(
  env_park$model_data1  %>% left_join(mean_times(env_park),  by = c("site.name", "hour")) %>% mutate(route = "park"),
  env_forum$model_data2 %>% left_join(mean_times(env_forum), by = c("site.name", "hour")) %>% mutate(route = "forum")
)

calc_tmrt <- function(Tg, Ta, Va, D = globe_D, epsilon = globe_eps) {
  ((Tg + 273.15)^4 + (1.1e8 * Va^0.6 / (epsilon * D^0.4)) * (Tg - Ta))^0.25 - 273.15
}

# calc_utci() converts the 1.77-m Kestrel wind to 10 m internally; clipping is
# done on the measured wind so that the converted 10-m wind is >= 0.5 m/s.
ws_min <- 0.5 / (log(10 / 0.01) / log(1.77 / 0.01))

d <- d %>%
  mutate(
    site       = factor(site.name),
    route      = factor(route),
    local_time = with_tz(time_mid, tz_local),
    t_local    = hour(local_time) + minute(local_time) / 60,   # decimal local hour
    hour_f     = factor(format(with_tz(hour, tz_local), "%H")),
    period = case_when(
              t_local >= 5  & t_local < 11 ~ "morning",
              t_local >= 11 & t_local < 16 ~ "midday",
              t_local >= 16 & t_local < 21 ~ "afternoon",
              TRUE                         ~ "night"),
    period     = factor(period, levels = c("morning", "midday", "afternoon", "night")),
    Tmrt       = calc_tmrt(Globe.Temp, Temp, Wind.Speed),
    ws_utci    = if (clip_wind_utci) pmax(Wind.Speed, ws_min) else Wind.Speed,
    UTCI       = calc_utci(ta = Temp, tmrt = Tmrt, hur = Rel..Hum., ws = ws_utci)
  ) %>%
  filter(!is.na(UTCI), !is.na(svf)) %>%
  filter(!(exclude_night & period == "night")) %>%
  droplevels() %>%
  group_by(hour) %>%
  mutate(UTCI_anom = UTCI - mean(UTCI)) %>%   # UTCI anomaly (see header)
  ungroup()

cat("\n==== DATA ====\n")
cat("site-hours:", nrow(d), " sites:", nlevels(d$site), " hours:", nlevels(d$hour_f), "\n")
print(table(d$period))
write.csv(d %>% select(-hour_f), file.path(out_dir, "utci_svf_data.csv"), row.names = FALSE)

# plot style
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92")))
save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 200, bg = "white")
}
site_labels <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = site), size = 3.3, colour = "grey25",
                           box.padding = 0.4, min.segment.length = 0.3, seed = 1)
} else {
  geom_text(aes(label = site), size = 3.3, colour = "grey25", vjust = -0.9)
}

# simple regression of site-mean UTCI anomaly on SVF (n = 10 sites)
site_regression <- function(site_means, label) {
  ct  <- cor.test(site_means$svf, site_means$UTCI_anom)          # Pearson
  fit <- lm(UTCI_anom ~ svf, data = site_means)
  ci  <- confint(fit)["svf", ]
  data.frame(subset = label, n_sites = nrow(site_means),
             slope_per_0.1SVF = coef(fit)[["svf"]] / 10,
             lwr = ci[[1]] / 10, upr = ci[[2]] / 10,
             r = unname(ct$estimate), R2 = summary(fit)$r.squared,
             t = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value)
}

#### 1 descriptives ####
site_table <- d %>%
  group_by(route, site, svf) %>%
  summarise(n_hours = n(), UTCI_mean = mean(UTCI), UTCI_min = min(UTCI),
            UTCI_max = max(UTCI), UTCI_anom_mean = mean(UTCI_anom), .groups = "drop") %>%
  arrange(svf) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
cat("\n==== 1 SITE SUMMARY ====\n"); print(as.data.frame(site_table))
write.csv(site_table, file.path(out_dir, "01_site_summary.csv"), row.names = FALSE)

p1 <- ggplot(d, aes(local_time, UTCI, group = site, colour = svf)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
  scale_colour_gradient(low = "#9ec5f4", high = "#0d366b", name = "SVF") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours") +
  labs(x = "Local time (CEST)", y = "UTCI (\u00b0C)",
       title = "UTCI at the 10 sites, coloured by sky view factor")
save_plot(p1, "01_utci_diurnal_by_site.png")

#### 2 H1: is UTCI related to SVF? ####
cat("\n==== 2 H1: SITE-MEAN UTCI ANOMALY vs SVF (n = 10 sites) ====\n")
site_means <- d %>% group_by(site, route, svf) %>%
  summarise(UTCI_anom = mean(UTCI_anom), .groups = "drop")
h1 <- site_regression(site_means, "all measured hours")
print(h1, digits = 3)
write.csv(h1, file.path(out_dir, "02_H1_regression.csv"), row.names = FALSE)

p2 <- ggplot(site_means, aes(svf, UTCI_anom)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2a78d6", fill = "#2a78d6",
              alpha = 0.15, linewidth = 0.7) +
  geom_point(aes(shape = route), colour = "#0d366b", size = 3) +
  site_labels +
  scale_shape_manual(values = c(forum = 16, park = 17), name = "Route") +
  labs(x = "Sky view factor", y = "Mean UTCI anomaly (K)",
       title = "Site-mean UTCI anomaly vs sky view factor",
       subtitle = sprintf("n = 10 sites: r = %.2f, p = %.3f, slope = %+.2f K per 0.1 SVF",
                          h1$r, h1$p, h1$slope_per_0.1SVF))
save_plot(p2, "02_H1_site_mean_vs_svf.png", w = 7, h = 5)

#### 3 H2: does the relationship change over the day? ####
# (a) the H1 regression, separately for each period of the day
cat("\n==== 3a H2: REGRESSION PER PERIOD (n = 10 sites each) ====\n")
period_means <- d %>% group_by(period, site, route, svf) %>%
  summarise(UTCI_anom = mean(UTCI_anom), .groups = "drop")
h2_periods <- bind_rows(lapply(levels(d$period), function(p)
  site_regression(filter(period_means, period == p), p)))
print(h2_periods, digits = 3)
write.csv(h2_periods, file.path(out_dir, "03a_H2_regression_per_period.csv"), row.names = FALSE)

p3a <- ggplot(period_means, aes(svf, UTCI_anom)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2a78d6", fill = "#2a78d6",
              alpha = 0.15, linewidth = 0.6) +
  geom_point(colour = "#0d366b", size = 2.2) +
  facet_wrap(~ period, nrow = 1) +
  labs(x = "Sky view factor", y = "Mean UTCI anomaly (K)",
       title = "Site-mean UTCI anomaly vs SVF by period of the day (n = 10 sites)")
save_plot(p3a, "03a_H2_site_means_by_period.png", w = 10, h = 4)

# (b) linear mixed model with all site-hours
#   UTCI ~ hour + route + SVF x period, random intercept for site
#   * hour (factor) removes the weather shared by all sites
#   * route accounts for the two instruments
#   * site random effect: repeated measurements at the same site
#   * AR(1): consecutive hours at a site are more similar than distant hours
cat("\n==== 3b H2: MIXED MODEL WITH SVF x PERIOD INTERACTION ====\n")
ar1 <- corCAR1(form = ~ t_local | site)
m_const <- lme(UTCI ~ route + hour_f + svf,              random = ~ 1 | site,
               correlation = ar1, data = d, method = "ML")
m_int   <- lme(UTCI ~ route + hour_f + svf + svf:period, random = ~ 1 | site,
               correlation = ar1, data = d, method = "ML")
cat("Does the SVF effect differ between periods? (likelihood-ratio test)\n")
lrt <- anova(m_const, m_int)
print(lrt)

# SVF effect in each period (per 0.1 SVF), refitted with REML for estimation.
# The period main effects are contained in hour_f, so only period:svf is added.
m_periods <- lme(UTCI ~ route + hour_f + period:svf, random = ~ 1 | site,
                 correlation = ar1, data = d)
tt   <- summary(m_periods)$tTable
rows <- grep("svf", rownames(tt))
h2_mixed <- data.frame(
  period = factor(sub("period(.*):svf", "\\1", rownames(tt)[rows]), levels(d$period)),
  slope  = tt[rows, "Value"] / 10,
  lwr    = (tt[rows, "Value"] - qt(0.975, tt[rows, "DF"]) * tt[rows, "Std.Error"]) / 10,
  upr    = (tt[rows, "Value"] + qt(0.975, tt[rows, "DF"]) * tt[rows, "Std.Error"]) / 10,
  t = tt[rows, "t-value"], df = tt[rows, "DF"], p = tt[rows, "p-value"], row.names = NULL)
print(h2_mixed, digits = 3)
write.csv(h2_mixed, file.path(out_dir, "03b_H2_mixed_model_period_slopes.csv"), row.names = FALSE)
write.csv(data.frame(chisq = lrt$L.Ratio[2], df = diff(lrt$df), p = lrt$`p-value`[2],
                     AIC_constant = lrt$AIC[1], AIC_interaction = lrt$AIC[2]),
          file.path(out_dir, "03b_H2_mixed_model_interaction_test.csv"), row.names = FALSE)

p3b <- ggplot(h2_mixed, aes(period, slope)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.12, colour = "#2a78d6", linewidth = 0.7) +
  geom_point(colour = "#2a78d6", size = 3.5) +
  geom_text(aes(label = sprintf("p = %.3f", p)), nudge_x = 0.28, size = 3.5, colour = "grey25") +
  labs(x = NULL, y = "Change in UTCI per +0.1 SVF (K)",
       title = "SVF effect on UTCI by period of the day (mixed model, 95% CI)",
       subtitle = sprintf("Difference between periods: chi-square = %.2f, df = %d, p = %.3f",
                          lrt$L.Ratio[2], diff(lrt$df), lrt$`p-value`[2]))
save_plot(p3b, "03b_H2_mixed_model_period_slopes.png", w = 7, h = 4.5)

# basic model check: residuals should scatter evenly around 0
p3b_res <- ggplot(data.frame(fitted = fitted(m_periods),
                             resid = resid(m_periods, type = "normalized")),
                  aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(colour = "#2a78d6", size = 1.5) +
  labs(x = "Fitted UTCI (\u00b0C)", y = "Normalised residual",
       title = "Residuals of the mixed model")
save_plot(p3b_res, "03b_mixed_model_residuals.png", w = 6, h = 4)

# (c) descriptive figures: SVF slope for every hour, and a heatmap
hourly <- d %>%
  group_by(hour) %>%
  group_modify(function(g, key) {
    fit <- lm(UTCI ~ svf, data = g)
    ci  <- confint(fit)["svf", ]
    tibble(local_time = with_tz(key$hour, tz_local), n = nrow(g),
           slope = coef(fit)[["svf"]] / 10, lwr = ci[[1]] / 10, upr = ci[[2]] / 10)
  }) %>% ungroup()
write.csv(hourly, file.path(out_dir, "03c_hourly_slopes.csv"), row.names = FALSE)

p3c <- ggplot(hourly, aes(local_time, slope)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#2a78d6", alpha = 0.15) +
  geom_line(colour = "#2a78d6", linewidth = 0.7) +
  geom_point(colour = "#2a78d6", size = 2.2) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours") +
  labs(x = "Local time (CEST)", y = "Change in UTCI per +0.1 SVF (K)",
       title = "SVF slope of UTCI for each measurement hour (95% CI)")
save_plot(p3c, "03c_hourly_slopes.png")

heat_data <- d %>%
  mutate(hour_local = as.integer(format(with_tz(hour, tz_local), "%H")),
         site_lab = reorder(factor(sprintf("%s (%.2f)", site, svf)), svf))
lim <- max(abs(heat_data$UTCI_anom))
p3c_heat <- ggplot(heat_data, aes(hour_local, site_lab, fill = UTCI_anom)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_vline(xintercept = c(10.5, 15.5), colour = "grey30", linetype = "dashed", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948",
                       midpoint = 0, limits = c(-lim, lim), name = "UTCI anomaly (K)") +
  scale_x_continuous(breaks = seq(5, 23, 1), expand = c(0, 0)) +
  labs(x = "Local time (h, CEST)", y = "Site (SVF), ordered by SVF",
       title = "UTCI anomaly relative to the hourly mean of all sites",
       subtitle = paste("Red = warmer than average, blue = cooler; white = no measurement;",
                        "dashed lines = period boundaries", sep = "\n")) +
  theme(panel.grid.major = element_blank())
save_plot(p3c_heat, "03c_heatmap_utci_anomaly.png", w = 10, h = 5)

cat("\nAll tables and figures written to", normalizePath(out_dir), "\n")
