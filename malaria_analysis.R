# ============================================================================
# malaria_analysis.R
# Full epidemiological analysis of the Kenya malaria dataset (Malaria Atlas
# Project surveys) for the Power BI-style dashboard. Mirrors the methodology
# in Calvin_Malaria.Project.Rmd and writes every statistic to malaria_stats.json
# so the dashboard generator consumes genuine R outputs.
#
# Run:  Rscript malaria_analysis.R
# ============================================================================

suppressMessages({
  library(jsonlite)
  library(dplyr)
})

# ---- 1. Load & clean (identical rules to the R project) --------------------
raw <- read.csv("kenya_malaria_raw.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8")

clean <- raw %>%
  filter(!is.na(pr), !is.na(examined), !is.na(positive), !is.na(year_start),
         examined > 0, positive <= examined) %>%
  mutate(
    year       = as.integer(year_start),
    prev_pct   = pr * 100,
    setting    = ifelse(rural_urban %in% c("RURAL", "URBAN", "PERI_URBAN"),
                        rural_urban, "Not recorded"),
    method     = ifelse(is.na(method) | method == "", "Not recorded", method)
  )

n_raw     <- nrow(raw)
n_surveys <- nrow(clean)
n_exam    <- sum(clean$examined)
n_pos     <- sum(clean$positive)
overall   <- n_pos / n_exam * 100
n_coords  <- sum(!is.na(clean$latitude))

# ---- 2. Helper: weighted prevalence over a subset --------------------------
wpr <- function(d) sum(d$positive) / sum(d$examined) * 100

# ---- 3. Temporal -----------------------------------------------------------
yearly <- clean %>%
  group_by(year) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined),
            mean_prev = mean(prev_pct), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

decades <- data.frame(
  label = c("1985\u201389", "1990\u201399", "2000\u201304", "2005\u201309", "2010\u201320"),
  lo    = c(1985, 1990, 2000, 2005, 2010),
  hi    = c(1989, 1999, 2004, 2009, 2020)
)
decade_rows <- lapply(seq_len(nrow(decades)), function(i) {
  d <- clean[clean$year >= decades$lo[i] & clean$year <= decades$hi[i], ]
  list(label = decades$label[i], surveys = nrow(d), positive = sum(d$positive),
       examined = sum(d$examined), weighted_prev = round(wpr(d), 2))
})

# Month-of-year seasonality (month the survey fieldwork started)
month_rows <- clean %>%
  filter(!is.na(month_start)) %>%
  group_by(month = month_start) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

# ---- 4. Geography ----------------------------------------------------------
zone_of <- function(lon, lat) {
  ifelse(is.na(lon) | is.na(lat), "No coordinates",
    ifelse(lon < 36.0, "Western (Lake Victoria basin)",
      ifelse(lon > 38.5 & lat < 0.5, "Coast",
        ifelse(lat > 1.0 & lon > 38.0, "North & East", "Central / highlands"))))
}
clean$zone <- zone_of(clean$longitude, clean$latitude)

zone_rows <- clean %>%
  group_by(zone) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

# Latitude bands: altitude proxy via latitude (highlands = higher latitude inland)
lat_rows <- clean %>%
  filter(!is.na(latitude)) %>%
  mutate(band = cut(latitude, breaks = c(-6, -1, 0.5, 2.5, 6),
                    labels = c("Coastal plain (<1\u00b0S)", "Lowland (1\u00b0S\u20130.5\u00b0N)",
                               "Highlands (0.5\u20132.5\u00b0N)", "Far north (>2.5\u00b0N)"))) %>%
  group_by(band) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

hotspots <- clean %>%
  arrange(desc(prev_pct)) %>%
  slice_head(n = 10) %>%
  select(site_name, year, setting, examined, positive, prev_pct, method)

# ---- 5. Breakdowns ---------------------------------------------------------
setting_rows <- clean %>%
  group_by(setting) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined),
            mean_prev = mean(prev_pct), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

method_rows <- clean %>%
  group_by(method) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined),
            mean_prev = mean(prev_pct), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

age_rows <- clean %>%
  mutate(age_group = case_when(
    is.na(lower_age) | is.na(upper_age) ~ "Not recorded",
    lower_age == 0 & upper_age <= 6 ~ "Under 6 (0\u20135 yrs)",
    lower_age < 15 & upper_age <= 18 ~ "Children / school-age",
    lower_age >= 15 ~ "Adults (15+)",
    TRUE ~ "All / mixed ages")) %>%
  group_by(age_group) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

rdt_rows <- clean %>%
  filter(!is.na(rdt_type), rdt_type != "") %>%
  count(rdt_type, name = "surveys") %>%
  arrange(desc(surveys))

# Trend by setting & by method (weighted prevalence per year)
setting_year <- clean %>%
  filter(setting %in% c("RURAL", "URBAN")) %>%
  group_by(year, setting) %>%
  summarise(positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

method_year <- clean %>%
  filter(method %in% c("Microscopy", "RDT")) %>%
  group_by(year, method) %>%
  summarise(positive = sum(positive), examined = sum(examined), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

# ---- 6. Sample-size reliability --------------------------------------------
size_bins <- clean %>%
  mutate(bin = cut(examined, breaks = c(0, 20, 50, 100, 250, 500, 1e9),
                   labels = c("<20", "20\u201349", "50\u201399", "100\u2013249",
                              "250\u2013499", "500+"))) %>%
  group_by(bin) %>%
  summarise(surveys = n(), positive = sum(positive), examined = sum(examined),
            mean_prev = mean(prev_pct), .groups = "drop") %>%
  mutate(weighted_prev = positive / examined * 100)

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
sum_stats <- list(
  surveys = n_surveys, examined = n_exam, positive = n_pos,
  weighted_prev = round(overall, 1),
  median_exam = q(clean$examined, 0.5),
  iqr_exam = c(q(clean$examined, 0.25), q(clean$examined, 0.75)),
  median_prev = round(q(clean$prev_pct, 0.5), 1),
  iqr_prev = round(c(q(clean$prev_pct, 0.25), q(clean$prev_pct, 0.75)), 1),
  largest_site = clean$site_name[which.max(clean$examined)],
  largest_exam = max(clean$examined),
  largest_year = clean$year[which.max(clean$examined)]
)

# Boxplot stats for method / setting
box_stats <- function(d, by) {
  lapply(split(d, d[[by]]), function(g) {
    list(group = g[[by]][1], n = nrow(g), min = round(min(g$prev_pct), 1),
         q1 = round(q(g$prev_pct, 0.25), 1), median = round(q(g$prev_pct, 0.5), 1),
         q3 = round(q(g$prev_pct, 0.75), 1), max = round(max(g$prev_pct), 1))
  })
}
box_method <- unname(box_stats(clean, "method"))
box_setting <- unname(box_stats(clean, "setting"))

# ---- 7. Logistic regression (the real R glm) -------------------------------
model_data <- clean %>%
  filter(setting %in% c("RURAL", "URBAN")) %>%
  mutate(setting = factor(setting, levels = c("RURAL", "URBAN")),
         method = factor(method, levels = c("Microscopy", "RDT")))

fit <- glm(cbind(positive, examined - positive) ~ year + setting + method,
           family = binomial, data = model_data)

or_ci <- exp(cbind(OR = coef(fit), confint.default(fit)))
pvals <- summary(fit)$coefficients[, 4]
model_rows <- list(
  list(term = "Survey year (per +1 year)", or = round(or_ci["year", 1], 3),
       lo = round(or_ci["year", 2], 3), hi = round(or_ci["year", 3], 3),
       p = pvals["year"]),
  list(term = "Setting: Urban vs Rural", or = round(or_ci["settingURBAN", 1], 3),
       lo = round(or_ci["settingURBAN", 2], 3), hi = round(or_ci["settingURBAN", 3], 3),
       p = pvals["settingURBAN"]),
  list(term = "Method: RDT vs Microscopy", or = round(or_ci["methodRDT", 1], 3),
       lo = round(or_ci["methodRDT", 2], 3), hi = round(or_ci["methodRDT", 3], 3),
       p = pvals["methodRDT"])
)

# ---- 8. Assemble + write JSON ----------------------------------------------
out <- list(
  kpi = list(n_surveys = n_surveys, n_raw = n_raw, n_dropped = n_raw - n_surveys,
             examined = n_exam, positive = n_pos, pr = round(overall, 1),
             year_min = min(clean$year), year_max = max(clean$year),
             n_coords = n_coords),
  yearly = yearly,
  decades = decade_rows,
  months = month_rows,
  zones = zone_rows,
  lat_bands = lat_rows,
  settings = setting_rows,
  methods = method_rows,
  ages = age_rows,
  rdt = rdt_rows,
  setting_year = setting_year,
  method_year = method_year,
  hotspots = hotspots,
  size_bins = size_bins,
  sum_stats = sum_stats,
  box_method = box_method,
  box_setting = box_setting,
  model = model_rows,
  model_extra = list(n = nrow(model_data), deviance = round(deviance(fit), 0),
                     aic = round(AIC(fit), 0), year_pct = round((1 - exp(coef(fit)["year"])) * 100, 1))
)

writeLines(toJSON(out, auto_unbox = TRUE, na = "null", digits = 6),
           "malaria_stats.json")

cat(sprintf("R analysis complete: %d surveys | %s examined | %s positive | %.1f%% prevalence\n",
            n_surveys, format(n_exam, big.mark = ","), format(n_pos, big.mark = ","), overall))
cat(sprintf("Model: year OR %.3f | urban OR %.3f | RDT OR %.3f | n = %d\n",
            exp(coef(fit)["year"]), exp(coef(fit)["settingURBAN"]),
            exp(coef(fit)["methodRDT"]), nrow(model_data)))
cat("Wrote malaria_stats.json\n")
