# Kenya Malaria Prevalence Dashboard

An interactive, Power BI-style dashboard of *Plasmodium falciparum* malaria
prevalence in Kenya, built from 2,149 geo-referenced community and
school-based surveys (Malaria Atlas Project, 1985–2020). Every statistic on
the dashboard is computed in **R** — cleaning, summaries, and a binomial
logistic regression — and rendered as a self-contained HTML dashboard that
works offline in any browser.

> **Live demo:** https://calvinokoth9528-cloud.github.io/kenya-malaria-dashboard/

## Visuals

![Dashboard overview](docs/screenshots/overview.png)

![Statistical model section](docs/screenshots/model.png)

## What the dashboard shows

| Section | Content |
|---|---|
| **KPI cards** | 2,149 surveys · 186,765 individuals examined · 38,161 infections · 20.4% overall prevalence |
| **Temporal trends** | Weighted prevalence by year (fell from ~46% in 1985 to ~8% by 2020) and by decade |
| **Seasonality & trends** | Prevalence by month of year, rural vs urban by year, microscopy vs RDT by year |
| **Geography** | Interactive map of all survey sites with a Kenya outline, plus prevalence by region and latitude band |
| **Breakdowns** | Setting, diagnostic method, age group, RDT types, hotspots |
| **Reliability** | Sample-size vs prevalence, per-survey statistics |
| **Statistical model** | Binomial logistic regression: year OR **0.930** (~7% decline/year), urban OR **0.483**, RDT OR **0.501** |
| **Data explorer** | Searchable table of all 2,149 surveys |

Power BI-style slicers (setting, diagnostic method, survey year) re-filter
every chart and the data table instantly.

## Dataset

- **Source:** Malaria Atlas Project (`malariaAtlas` R package), Kenya surveys
- **Scope:** 3,607 raw survey records; 2,149 retained after cleaning
  (1,458 removed for missing or invalid parasitological data)
- **Recorded fields:** site, coordinates, year and month, age range,
  individuals examined, positive results, parasite rate, diagnostic method,
  RDT type, species (all *P. falciparum*)

## Methodology (R)

`malaria_analysis.R` (R 4.6.1, `tidyverse`) performs the full analysis and
writes `malaria_stats.json`:

1. **Cleaning** — drop missing/invalid records; derive prevalence and
   survey-level variables (identical rules to the analysis report).
2. **Summaries** — weighted prevalence (infections ÷ examined) by year,
   decade, month, zone, setting, method, age group, latitude band, and
   sample-size band.
3. **Modelling** — binomial logistic regression
   `cbind(positive, examined - positive) ~ year + setting + method`
   (RURAL and Microscopy as reference levels), reported as adjusted odds
   ratios with 95% confidence intervals.

## Run it locally

Open `index.html` in any browser — no server or internet required (Chart.js
is vendored).

## Reproduce the analysis

```bash
# 1. Compute all statistics in R (writes malaria_stats.json)
Rscript malaria_analysis.R

# 2. Rebuild the dashboard from the R outputs (writes index.html)
python build_malaria_dashboard.py index.html
```

## Repository structure

```
├── index.html                 # Interactive dashboard (self-contained)
├── malaria_analysis.R         # R analysis -> malaria_stats.json
├── malaria_stats.json         # All statistics computed by R
├── build_malaria_dashboard.py # Python generator (reads R outputs + CSV)
├── kenya_malaria_raw.csv      # Source survey data (Malaria Atlas Project)
├── chart.umd.min.js           # Vendored Chart.js (offline support)
├── _kenya_outline.json        # Kenya boundary for the map
└── docs/screenshots/          # Dashboard visuals
```

## Skills demonstrated

- R data analysis (data cleaning, tidyverse summarisation, binomial
  logistic regression with `glm()`)
- Epidemiological interpretation (spatial, temporal and seasonal trends)
- Web dashboard development (HTML, JavaScript, Chart.js)
- Power BI-style interactive visual design with cross-filtering
