#' Statistical analysis: relationship between UTCI and sky view factor (SVF)
#' Routes forum + park (10 sites), Berlin-Adlershof, 2026-05-29
#' Author: Adelina Petkova
#' Date: 09.2026
#'
#' Hypothesis
#'   H1: UTCI depends on SVF (open sites are thermally more stressful).
#'   H2: The strength / sign of that relationship changes over the day
#'       (expected: positive during high sun, ~0 or negative at night).
#'
#' Design problems the statistics must respect
#'   * SVF is a SITE property -> only 10 independent SVF values. The ~190
#'     site-hours are repeated measures, not independent replicates
#'     (pseudo-replication). Site is therefore a random effect, and the
#'     test of SVF has at most 10 - (number of site-level terms) df.
#'   * All sites share the same weather -> the diurnal cycle has to be
#'     removed (hour fixed effects / smooth of time) before SVF is tested.
#'   * Each route was measured with its own Kestrel -> route (= instrument)
#'     enters as a fixed effect.
#'   * Sites are visited in a fixed order within each hour, so a site is
#'     always measured at the same minute. The GAM uses the exact
#'     measurement time to account for this.
#'
#' Analyses
#'   0  data preparation (re-uses analysis_UA.R and analysis_forum.R)
#'   1  descriptives
#'   2  H1 - is there a relationship at all?
#'        2a site-level correlation + permutation test (n = 10)
#'        2b linear mixed model UTCI ~ SVF + route + hour + (1 | site)
#'   3  H2 - how does it change over the day?
#'        3a hour-by-hour cross-sectional regressions (n = 10 per hour)
#'        3b mixed model with SVF x period interaction
#'        3c GAM with a time-varying SVF coefficient
#'   4  mechanism: same hourly regressions for Tmrt, Ta and wind speed
#'   5  diagnostics & leave-one-site-out sensitivity
#'
#' Outputs are written to output/ (csv tables + png figures).

#### settings ####
library(ggplot2)
library(dplyr)
library(lubridate)
library(nlme)   # mixed models (ships with R)
library(mgcv)   # GAM (ships with R)

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)

tz_local <- "Europe/Berlin"        # UTC+2 (CEST) during the campaign
lat <- 52.433; lon <- 13.528       # centre of the two routes (protocol)

# Globe diameter [m] used for Tmrt. analysis_UA.R uses 0.15 m (ISO 7726
# standard globe). The Kestrel 5400 globe is a 1-inch (0.0254 m) globe, so
# check which value is correct for your instrument - it changes Tmrt a lot.
globe_D <- 0.15
globe_eps <- 0.95

# UTCI is only defined for 10-m wind speeds >= 0.5 m/s. The Kestrel reports
# 0 m/s below its stall speed, so the 10-m wind is clipped at 0.5 m/s.
clip_wind_utci <- TRUE

set.seed(42)

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
    summarise(time_mid = mean(Time), n_raw = n(), .groups = "drop")
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
    period     = case_when(t_local  <  5   ~ "night",
                           t_local  <  11  ~ "morning",
                           t_local  <  15  ~ "midday",
                           TRUE            ~ "afternoon"),
    period     = factor(period, levels = c("morning", "midday", "afternoon", "night")),
    Tmrt       = calc_tmrt(Globe.Temp, Temp, Wind.Speed),
    ws_utci    = if (clip_wind_utci) pmax(Wind.Speed, ws_min) else Wind.Speed,
    UTCI       = calc_utci(ta = Temp, tmrt = Tmrt, hur = Rel..Hum., ws = ws_utci),
    D_Tmrt     = Tmrt - Temp
  ) %>%
  filter(!is.na(UTCI), !is.na(svf)) %>%
  group_by(hour) %>%
  # anomaly relative to the mean of all sites measured in the same hour:
  # removes the common weather signal
  mutate(UTCI_anom = UTCI - mean(UTCI), Tmrt_anom = Tmrt - mean(Tmrt),
         Ta_anom = Temp - mean(Temp), ws_anom = Wind.Speed - mean(Wind.Speed)) %>%
  ungroup()

cat("\n==== DATA ====\n")
cat("site-hours:", nrow(d), " sites:", nlevels(d$site),
    " hours:", nlevels(d$hour_f), "\n")
cat("UTCI outside validity (Tmrt - Ta not in [-30, 70] K):",
    sum(d$D_Tmrt < -30 | d$D_Tmrt > 70), "\n")
cat("Wind speed clipped for UTCI:", sum(d$Wind.Speed < ws_min), "of", nrow(d), "\n")
print(table(d$period))

write.csv(d %>% select(-hour_f), file.path(out_dir, "utci_svf_data.csv"), row.names = FALSE)

# SVF colour scale: sequential blue, light = enclosed, dark = open sky
svf_scale <- scale_colour_gradient(low = "#9ec5f4", high = "#0d366b", name = "SVF")
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92")))
save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 200, bg = "white")
}

#### 1 descriptives ####
site_table <- d %>%
  group_by(route, site, svf) %>%
  summarise(n = n(),
            UTCI_mean = mean(UTCI), UTCI_min = min(UTCI), UTCI_max = max(UTCI),
            UTCI_anom_mean = mean(UTCI_anom),
            UTCI_anom_day = mean(UTCI_anom[period != "night"]),
            Tmrt_mean = mean(Tmrt), Ta_mean = mean(Temp), ws_mean = mean(Wind.Speed),
            .groups = "drop") %>%
  arrange(svf) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
cat("\n==== 1 SITE SUMMARY ====\n"); print(as.data.frame(site_table))
write.csv(site_table, file.path(out_dir, "01_site_summary.csv"), row.names = FALSE)

p1 <- ggplot(d, aes(local_time, UTCI, group = site, colour = svf)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.2) + svf_scale +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours") +
  labs(x = "Local time (CEST)", y = "UTCI (\u00b0C)",
       title = "UTCI at the 10 sites, coloured by sky view factor")
save_plot(p1, "01_utci_diurnal_by_site.png")

#### 2 H1: is there a relationship at all? ####
cat("\n==== 2a SITE-LEVEL CORRELATION (n = 10 sites) ====\n")
# Honest unit of analysis for a site property: one value per site.
perm_test_r <- function(x, y, n_perm = 9999) {
  r_obs  <- cor(x, y)
  r_perm <- replicate(n_perm, cor(x, sample(y)))
  c(r = r_obs, p_perm = (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1))
}
# mean UTCI anomaly (to the hourly all-site mean) per site, for the whole
# campaign, daytime only and for each period of the day
site_period <- bind_rows(
  d %>% mutate(subset = "all hours"),
  d %>% filter(period != "night") %>% mutate(subset = "daytime"),
  d %>% mutate(subset = as.character(period))
) %>%
  group_by(subset, site, svf) %>%
  summarise(UTCI_anom = mean(UTCI_anom), .groups = "drop")
subsets <- c("all hours", "daytime", levels(d$period))
site_level <- bind_rows(lapply(subsets, function(v) {
  g  <- filter(site_period, subset == v)
  y  <- g$UTCI_anom; x <- g$svf
  pt <- perm_test_r(x, y)
  sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  fit <- lm(y ~ x)
  data.frame(subset = v, n_sites = nrow(g),
             slope_per_0.1SVF = unname(coef(fit)[2]) / 10,
             pearson_r = pt[["r"]], p_permutation = pt[["p_perm"]],
             spearman_rho = unname(sp$estimate), p_spearman = sp$p.value)
}))
print(site_level, digits = 3)
write.csv(site_level, file.path(out_dir, "02a_site_level_correlation.csv"), row.names = FALSE)

p2a <- ggplot(filter(site_period, subset %in% levels(d$period)) %>%
                mutate(subset = factor(subset, levels(d$period))),
              aes(svf, UTCI_anom)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2a78d6", fill = "#2a78d6",
              alpha = 0.15, linewidth = 0.6) +
  geom_point(colour = "#0d366b", size = 2.2) +
  facet_wrap(~ subset, nrow = 1) +
  labs(x = "Sky view factor", y = "Mean UTCI anomaly (K)",
       title = "Site-mean UTCI anomaly vs SVF by period of the day (n = 10 sites)")
save_plot(p2a, "02a_site_level_by_period.png", w = 11, h = 4)

cat("\n==== 2b LINEAR MIXED MODEL: UTCI ~ SVF + route + hour + (1|site), AR(1) ====\n")
# hour_f absorbs the shared weather; the SVF effect is tested against the
# between-site variance (denominator df = sites - site-level terms).
# Residuals of consecutive hours at one site are correlated -> continuous-time
# AR(1) error structure (see section 5 for the comparison without it).
ar1 <- corCAR1(form = ~ t_local | site)
m_null <- lme(UTCI ~ route + hour_f,       random = ~ 1 | site, correlation = ar1,
              data = d, method = "ML")
m_svf  <- lme(UTCI ~ route + hour_f + svf, random = ~ 1 | site, correlation = ar1,
              data = d, method = "ML")
m_svf_reml <- update(m_svf, method = "REML")
print(summary(m_svf_reml)$tTable["svf", , drop = FALSE])
cat("Likelihood-ratio test for SVF:\n"); print(anova(m_null, m_svf))
ci_svf <- intervals(m_svf_reml, which = "fixed")$fixed["svf", ]
cat(sprintf("SVF effect: %.2f K per 0.1 SVF (95%% CI %.2f to %.2f)\n",
            ci_svf[2] / 10, ci_svf[1] / 10, ci_svf[3] / 10))

#### 3 H2: how does it change over the day? ####
cat("\n==== 3a HOUR-BY-HOUR REGRESSIONS (cross-sectional, n ~ 10 per hour) ====\n")
hourly_reg <- function(data, response) {
  data %>%
    group_by(hour) %>%
    group_modify(function(g, key) {
      fit <- lm(reformulate("svf", response), data = g)
      ci  <- confint(fit)["svf", ]
      sp  <- suppressWarnings(cor.test(g$svf, g[[response]], method = "spearman"))
      tibble(n = nrow(g),
             local_time = with_tz(key$hour, tz_local) + 30 * 60,
             sun_elev = mean(g$sun_elev),
             slope = coef(fit)[["svf"]] / 10,          # per 0.1 SVF
             lwr = ci[[1]] / 10, upr = ci[[2]] / 10,
             r = cor(g$svf, g[[response]]),
             p = summary(fit)$coefficients["svf", 4],
             rho = unname(sp$estimate), p_spearman = sp$p.value)
    }) %>%
    ungroup() %>%
    mutate(response = response,
           p_BH = p.adjust(p, "BH"))          # multiple testing over 19 hours
}
hr_utci <- hourly_reg(d, "UTCI")
print(as.data.frame(hr_utci %>% transmute(local = format(local_time, "%H:%M"), n,
                                          sun_elev = round(sun_elev, 1),
                                          slope = round(slope, 2), lwr = round(lwr, 2),
                                          upr = round(upr, 2), r = round(r, 2),
                                          p = signif(p, 2), p_BH = signif(p_BH, 2),
                                          rho = round(rho, 2))))
write.csv(hr_utci, file.path(out_dir, "03a_hourly_regressions_UTCI.csv"), row.names = FALSE)

p3a <- ggplot(hr_utci, aes(local_time, slope)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#2a78d6", alpha = 0.15) +
  geom_line(colour = "#2a78d6", linewidth = 0.7) +
  geom_point(aes(shape = p_BH < 0.05), colour = "#2a78d6", size = 2.5) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("not significant", "p(BH) < 0.05"), name = NULL) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours") +
  labs(x = "Local time (CEST)", y = "Change in UTCI per +0.1 SVF (K)",
       title = "Hourly cross-sectional slope of UTCI on SVF (95% CI)")
save_plot(p3a, "03a_hourly_slope_UTCI_SVF.png")

p3a_sc <- ggplot(d, aes(svf, UTCI_anom)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2a78d6", fill = "#2a78d6",
              alpha = 0.15, linewidth = 0.6) +
  geom_point(aes(shape = route), size = 1.8, colour = "#0d366b") +
  facet_wrap(~ format(with_tz(hour, tz_local), "%H:00")) +
  labs(x = "Sky view factor", y = "UTCI anomaly to hourly mean (K)",
       title = "UTCI vs SVF for each measurement hour")
save_plot(p3a_sc, "03a_scatter_by_hour.png", w = 11, h = 8)

cat("\n==== 3b MIXED MODEL WITH SVF x PERIOD INTERACTION ====\n")
m_add <- lme(UTCI ~ route + hour_f + svf,              random = ~ 1 | site,
             correlation = ar1, data = d, method = "ML")
m_int <- lme(UTCI ~ route + hour_f + svf + svf:period, random = ~ 1 | site,
             correlation = ar1, data = d, method = "ML")
cat("Likelihood-ratio test: does the SVF effect differ between periods?\n")
print(anova(m_add, m_int))

# period-specific SVF slopes (per 0.1 SVF) from the interaction model.
# Note: nlme uses within-site df for these terms, which is optimistic for a
# site-level predictor -> compare with the site-level tests in 2a.
m_int_reml <- lme(UTCI ~ route + hour_f + period:svf, random = ~ 1 | site,
                  correlation = ar1, data = d)
tt <- summary(m_int_reml)$tTable
rows <- grep("svf", rownames(tt))
period_slopes <- data.frame(
  period = sub("period(.*):svf", "\\1", rownames(tt)[rows]),
  slope  = tt[rows, "Value"] / 10,
  lwr    = (tt[rows, "Value"] - qt(0.975, tt[rows, "DF"]) * tt[rows, "Std.Error"]) / 10,
  upr    = (tt[rows, "Value"] + qt(0.975, tt[rows, "DF"]) * tt[rows, "Std.Error"]) / 10,
  df     = tt[rows, "DF"], p = tt[rows, "p-value"], row.names = NULL)
print(period_slopes, digits = 3)
write.csv(period_slopes, file.path(out_dir, "03b_period_slopes_mixed_model.csv"), row.names = FALSE)

cat("\n==== 3c GAM WITH TIME-VARYING SVF EFFECT ====\n")
# UTCI = route + f0(time) + f1(time) * SVF + site random effect
# f0 = shared diurnal cycle, f1 = SVF effect as a smooth function of time
g_const <- gam(UTCI ~ route + s(t_local, k = 10) + svf + s(site, bs = "re"),
               data = d, method = "REML")
g_vary  <- gam(UTCI ~ route + s(t_local, k = 10) + s(t_local, by = svf, k = 10) +
                 s(site, bs = "re"), data = d, method = "REML")
print(summary(g_vary))
cat("AIC constant SVF effect:", round(AIC(g_const), 1),
    "| AIC time-varying SVF effect:", round(AIC(g_vary), 1), "\n")
print(anova(g_const, g_vary, test = "F"))
# The GAM above assumes independent residuals (optimistic p-values). Refit the
# same model with AR(1) errors as a check.
g_vary_ar <- gamm(UTCI ~ route + s(t_local, k = 10) + s(t_local, by = svf, k = 10),
                  random = list(site = ~ 1), correlation = ar1, data = d, method = "REML")
cat("Time-varying SVF term with AR(1) errors:\n")
print(summary(g_vary_ar$gam)$s.table)

# f1(t) with 95 % CI
nd <- data.frame(t_local = seq(min(d$t_local), max(d$t_local), length.out = 200),
                 svf = 1, route = levels(d$route)[1], site = levels(d$site)[1])
pr <- predict(g_vary, nd, type = "terms", se.fit = TRUE)
term <- "s(t_local):svf"
gam_curve <- data.frame(t_local = nd$t_local,
                        slope = pr$fit[, term] / 10,
                        lwr = (pr$fit[, term] - 1.96 * pr$se.fit[, term]) / 10,
                        upr = (pr$fit[, term] + 1.96 * pr$se.fit[, term]) / 10)
write.csv(gam_curve, file.path(out_dir, "03c_gam_svf_effect_curve.csv"), row.names = FALSE)

p3c <- ggplot(gam_curve, aes(t_local, slope)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#2a78d6", alpha = 0.15) +
  geom_line(colour = "#2a78d6", linewidth = 0.8) +
  geom_point(data = hr_utci %>% mutate(t_local = hour(local_time) + minute(local_time) / 60),
             colour = "grey40", size = 1.8) +
  scale_x_continuous(breaks = seq(4, 24, 2)) +
  labs(x = "Local time (h, CEST)", y = "Change in UTCI per +0.1 SVF (K)",
       title = "Time-varying SVF effect on UTCI (GAM, 95% CI)",
       subtitle = "Line: GAM estimate; grey points: hour-by-hour regressions")
save_plot(p3c, "03c_gam_time_varying_svf_effect.png")

#### 4 mechanism: which UTCI input does SVF act on? ####
cat("\n==== 4 HOURLY SVF SLOPES FOR UTCI INPUTS ====\n")
mech <- bind_rows(hourly_reg(d, "UTCI"), hourly_reg(d, "Tmrt"),
                  hourly_reg(d, "Temp"), hourly_reg(d, "Wind.Speed")) %>%
  mutate(response = recode(response, UTCI = "UTCI (K)", Tmrt = "Tmrt (K)",
                           Temp = "Air temperature (K)", Wind.Speed = "Wind speed (m/s)"))
write.csv(mech, file.path(out_dir, "04_hourly_regressions_all_variables.csv"), row.names = FALSE)
print(as.data.frame(mech %>% group_by(response) %>%
                      summarise(mean_slope = mean(slope), n_sig_BH = sum(p_BH < 0.05))))

p4 <- ggplot(mech, aes(local_time, slope)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#2a78d6", alpha = 0.15) +
  geom_line(colour = "#2a78d6", linewidth = 0.6) +
  geom_point(aes(shape = p_BH < 0.05), colour = "#2a78d6", size = 2) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("not significant", "p(BH) < 0.05"), name = NULL) +
  facet_wrap(~ response, scales = "free_y") +
  scale_x_datetime(date_labels = "%H", date_breaks = "2 hours") +
  labs(x = "Local time (h, CEST)", y = "Change per +0.1 SVF",
       title = "Hourly SVF slopes for UTCI and its input variables")
save_plot(p4, "04_mechanism_hourly_slopes.png", w = 10, h = 6)

#### 5 diagnostics & sensitivity ####
cat("\n==== 5 DIAGNOSTICS ====\n")
res <- data.frame(fitted = fitted(m_int_reml), resid = resid(m_int_reml, type = "normalized"),
                  site = d$site, t_local = d$t_local)
cat("Shapiro-Wilk on normalised residuals: p =",
    signif(shapiro.test(res$resid)$p.value, 3), "\n")
# autocorrelation of the model WITHOUT the AR(1) term, to justify using it
m_iid <- update(m_int_reml, correlation = NULL)
res_iid <- data.frame(resid = resid(m_iid, type = "normalized"), site = d$site,
                      t_local = d$t_local)
lag1 <- res_iid %>% arrange(site, t_local) %>% group_by(site) %>%
  summarise(r1 = cor(resid[-n()], resid[-1]), .groups = "drop")
cat("Mean within-site lag-1 residual autocorrelation (independent errors):",
    round(mean(lag1$r1), 2), "\n")
cat("AIC independent errors / AR(1):", round(AIC(m_iid), 1), "/",
    round(AIC(m_int_reml), 1), "\n")
cat("Period slopes assuming independent errors (per 0.1 SVF), for comparison:\n")
tt_ar <- summary(m_iid)$tTable; rows_ar <- grep("svf", rownames(tt_ar))
print(data.frame(period = sub("period(.*):svf", "\\1", rownames(tt_ar)[rows_ar]),
                 slope = tt_ar[rows_ar, "Value"] / 10, p = tt_ar[rows_ar, "p-value"],
                 row.names = NULL), digits = 3)

p5 <- ggplot(res, aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(colour = "#2a78d6", size = 1.5) +
  labs(x = "Fitted UTCI (\u00b0C)", y = "Normalised residual",
       title = "Residuals of the SVF x period mixed model")
save_plot(p5, "05_residuals.png", w = 6, h = 4)

# leave-one-site-out: does a single site drive the result?
loso <- bind_rows(lapply(levels(d$site), function(s) {
  fit <- lme(UTCI ~ route + hour_f + period:svf, random = ~ 1 | site,
             correlation = ar1, data = droplevels(filter(d, site != s)))
  tt <- summary(fit)$tTable; rows <- grep("svf", rownames(tt))
  data.frame(left_out = s, period = sub("period(.*):svf", "\\1", rownames(tt)[rows]),
             slope = tt[rows, "Value"] / 10, p = tt[rows, "p-value"], row.names = NULL)
}))
cat("\nLeave-one-site-out range of period slopes (per 0.1 SVF):\n")
print(as.data.frame(loso %>% group_by(period) %>%
                      summarise(min_slope = min(slope), max_slope = max(slope),
                                n_p_below_0.05 = sum(p < 0.05))), digits = 3)
write.csv(loso, file.path(out_dir, "05_leave_one_site_out.csv"), row.names = FALSE)

cat("\nAll tables and figures written to", normalizePath(out_dir), "\n")
