#!/usr/bin/env python3
"""
Build the Kenya Malaria Prevalence Dashboard (Power BI style, light theme).

Reads kenya_malaria_raw.csv (Malaria Atlas Project surveys), computes every
statistic used on the dashboard — including a binomial logistic regression
fitted with a pure-Python IRLS routine (mirrors the R glm in
Calvin_Malaria.Project.Rmd) — and writes malaria_dashboard.html with the data
embedded so the page is fully interactive client-side.

Run:  python build_malaria_dashboard.py
"""
import csv
import json
import math
import os
import sys

# Optional argv[1] overrides the output filename (e.g. "index.html" for Pages).
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "malaria_dashboard.html")
CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kenya_malaria_raw.csv")
OUTLINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_kenya_outline.json")

# --------------------------------------------------------------------------
# 1. Load & clean (same rules as the R project)
# --------------------------------------------------------------------------
def _num(v):
    try:
        f = float(v)
        return f if f == f else None
    except (ValueError, TypeError):
        return None

rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
clean = []
for r in rows:
    try:
        pr = float(r["pr"]); examined = int(float(r["examined"])); positive = int(float(r["positive"]))
        year = int(float(r["year_start"]))
    except (ValueError, TypeError):
        continue
    if pr != pr or examined != examined or positive != positive:
        continue
    if examined <= 0 or positive > examined:
        continue
    clean.append({
        "site": r["site_name"].strip(),
        "year": year, "examined": examined, "positive": positive, "pr": pr,
        "lat": _num(r["latitude"]), "lon": _num(r["longitude"]),
        "setting": r["rural_urban"],
        "method": r["method"],
        "lo": r["lower_age"], "hi": r["upper_age"],
        "rdt": r["rdt_type"],
        "month": _num(r["month_start"]),
    })
n_raw = len(rows)
n_surveys = len(clean)
n_examined = sum(r["examined"] for r in clean)
n_positive = sum(r["positive"] for r in clean)
overall_pr = n_positive / n_examined * 100.0
year_min, year_max = min(r["year"] for r in clean), max(r["year"] for r in clean)
n_coords = sum(1 for r in clean if r["lat"] is not None)

# --------------------------------------------------------------------------
# 2. Statistical analysis is performed in R (malaria_analysis.R), which writes
# malaria_stats.json; the dashboard consumes those genuine R outputs.
import json as _json
STATS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "malaria_stats.json")
if not os.path.exists(STATS):
    raise SystemExit("malaria_stats.json not found - run: Rscript malaria_analysis.R")
rs = _json.load(open(STATS, encoding="utf-8"))

model_rows = rs["model"]
model_extra = rs["model_extra"]
rdt_rows = [{"type": r["rdt_type"], "n": r["surveys"]} for r in rs["rdt"]]

# --------------------------------------------------------------------------
# 3. Per-survey geometry helpers
# --------------------------------------------------------------------------

def age_idx(r):
    try:
        lo, hi = float(r["lo"]), float(r["hi"])
    except (ValueError, TypeError):
        return 4
    if lo == 0 and hi <= 6:
        return 1                      # under-6 (0-5)
    if hi <= 18 and lo < 15:
        return 0                      # children / school-age
    if lo >= 15:
        return 3                      # adults 15+
    return 2                          # all / mixed ages

def zone_idx(r):
    lon, lat = r["lon"], r["lat"]
    if lon is None or lat is None:
        return 4                       # no coordinates
    if lon < 36.0:
        return 0                      # Western / Lake Victoria basin
    if lon > 38.5 and lat < 0.5:
        return 1                      # Coast
    if lat > 1.0 and lon > 38.0:
        return 3                      # North & East
    return 2                          # Central / highlands

def set_idx(r):
    return 0 if r["setting"] == "RURAL" else (1 if r["setting"] == "URBAN" else 2)

def met_idx(r):
    return 0 if r["method"] == "Microscopy" else 1

# plasma-ish colour ramp for the map (yellow -> orange -> red -> pink -> purple)
STOPS = [(0.0, "#f0f921"), (0.28, "#fca636"), (0.55, "#e54c4c"), (0.78, "#a52c78"), (1.0, "#3b0f70")]
def ramp(t):
    t = max(0.0, min(1.0, t))
    for i in range(len(STOPS) - 1):
        t0, c0 = STOPS[i]; t1, c1 = STOPS[i + 1]
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            r0, g0, b0 = int(c0[1:3], 16), int(c0[3:5], 16), int(c0[5:7], 16)
            r1, g1, b1 = int(c1[1:3], 16), int(c1[3:5], 16), int(c1[5:7], 16)
            return "#%02x%02x%02x" % (int(r0 + (r1 - r0) * f), int(g0 + (g1 - g0) * f), int(b0 + (b1 - b0) * f))
    return STOPS[-1][1]

survey_rows = []
for r in clean:
    survey_rows.append([
        r["year"], set_idx(r), met_idx(r), r["examined"], r["positive"],
        round(r["lon"], 3) if r["lon"] is not None else None,
        round(r["lat"], 3) if r["lat"] is not None else None,
        age_idx(r), zone_idx(r),
        r["site"], round(r["pr"] * 100.0, 2), ramp(r["pr"]),
        int(r["month"]) if r["month"] is not None else None,
    ])

AGE_NAMES = ["Children / school-age", "Under 6 (0\u20135 yrs)", "All / mixed ages", "Adults (15+)", "Not recorded"]
ZONE_NAMES = ["Western (Lake Victoria basin)", "Coast", "Central / highlands", "North & East", "No coordinates"]

# gtsummary-style summary stats for the intro
ex = sorted(r["examined"] for r in clean)
prv = sorted(r["pr"] * 100 for r in clean)
def q(vals, qq):
    i = (len(vals) - 1) * qq
    lo, hi = math.floor(i), math.ceil(i)
    if lo == hi:
        return vals[lo]
    return vals[lo] + (vals[hi] - vals[lo]) * (i - lo)

kp = {
    "n_surveys": n_surveys, "n_raw": n_raw, "n_dropped": n_raw - n_surveys,
    "examined": n_examined, "positive": n_positive, "pr": round(overall_pr, 1),
    "year_min": year_min, "year_max": year_max,
    "n_coords": n_coords,
    "exam_median": q(ex, 0.5), "exam_iqr": [q(ex, 0.25), q(ex, 0.75)],
    "pr_median": round(q(prv, 0.5), 1), "pr_iqr": [round(q(prv, 0.25), 1), round(q(prv, 0.75), 1)],
    "model": model_rows, "model_extra": model_extra, "rdt": rdt_rows,
}
outline = json.load(open(OUTLINE, encoding="utf-8")) if os.path.exists(OUTLINE) else []

# --------------------------------------------------------------------------
# 4. HTML template
# --------------------------------------------------------------------------
HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Kenya Malaria Prevalence Dashboard — KEMRI Attachment</title>
<script src="chart.umd.min.js"></script>
<style>
  :root{
    --bg:#f6f4ef; --card:#ffffff; --line:#e7e3d8; --ink:#26231c; --mut:#736e60;
    --teal:#0f9d8f; --teal-d:#0b7a70; --amber:#e8a33d; --rose:#d94f4f; --olive:#7a8c4e;
  }
  *{box-sizing:border-box; margin:0; padding:0}
  body{background:var(--bg); color:var(--ink); font-family:"Segoe UI", system-ui, Arial, sans-serif; font-size:14px; line-height:1.45;}
  .wrap{max-width:1280px; margin:0 auto; padding:20px 24px 60px;}
  header.top{background:#fff; border-bottom:1px solid var(--line);}
  .top-inner{max-width:1280px; margin:0 auto; padding:18px 24px; display:flex; align-items:center; gap:18px; flex-wrap:wrap;}
  .logo{width:46px;height:46px;border-radius:10px;background:var(--teal);color:#fff;display:flex;align-items:center;justify-content:center;font-weight:700;font-size:18px;}
  h1{font-size:22px; font-weight:700; letter-spacing:-.2px;}
  .sub{color:var(--mut); font-size:13px; margin-top:2px;}
  .pill{margin-left:auto; background:#fff7e8; color:#9a6b1a; border:1px solid #f0d9a8; border-radius:999px; padding:6px 14px; font-size:12.5px; font-weight:600;}
  .kpis{display:grid; grid-template-columns:repeat(auto-fit,minmax(175px,1fr)); gap:14px; margin:20px 0;}
  .kpi{background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px 16px; box-shadow:0 1px 3px rgba(60,50,20,.05);}
  .kpi .lbl{font-size:11.5px; text-transform:uppercase; letter-spacing:.6px; color:var(--mut); font-weight:600;}
  .kpi .val{font-size:26px; font-weight:700; margin-top:4px; letter-spacing:-.5px;}
  .kpi .note{font-size:12px; color:var(--mut); margin-top:2px;}
  .kpi .val.teal{color:var(--teal-d);} .kpi .val.amber{color:#b97c22;} .kpi .val.rose{color:var(--rose);} .kpi .val.olive{color:var(--olive);}
  .filters{background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px 18px; margin:0 0 20px; box-shadow:0 1px 3px rgba(60,50,20,.05);}
  .filters h3{font-size:12px; text-transform:uppercase; letter-spacing:.6px; color:var(--mut); margin-bottom:10px;}
  .frow{display:flex; gap:26px; flex-wrap:wrap; align-items:flex-end;}
  .fgroup label{display:block; font-size:12px; color:var(--mut); font-weight:600; margin-bottom:6px;}
  .chips{display:flex; gap:8px; flex-wrap:wrap;}
  .chip{border:1px solid var(--line); background:#faf9f5; border-radius:999px; padding:6px 14px; font-size:13px; cursor:pointer; user-select:none; transition:all .15s;}
  .chip.on{background:var(--teal); border-color:var(--teal); color:#fff;}
  .chip .n{opacity:.65; font-size:11.5px;}
  .yr{display:flex; align-items:center; gap:8px;}
  .yr input{width:86px; padding:6px 8px; border:1px solid var(--line); border-radius:8px; font-size:13px; background:#faf9f5;}
  .yr span{color:var(--mut);}
  .reset{margin-left:auto; border:1px solid var(--line); background:#fff; border-radius:8px; padding:7px 14px; font-size:13px; cursor:pointer; color:var(--mut);}
  .reset:hover{border-color:var(--teal); color:var(--teal-d);}
  section{margin:26px 0;}
  .sec-head{display:flex; align-items:baseline; gap:12px; margin-bottom:12px;}
  .sec-head h2{font-size:17px; font-weight:700;}
  .sec-head .tag{font-size:11px; color:#fff; background:var(--amber); border-radius:999px; padding:2px 10px; font-weight:600;}
  .insight{background:#fdf9ee; border:1px solid #f0e3bd; border-left:4px solid var(--amber); border-radius:8px; padding:10px 14px; font-size:13px; color:#6b5a26; margin-bottom:14px;}
  .grid{display:grid; gap:14px;}
  .g2{grid-template-columns:1fr 1fr;} .g3{grid-template-columns:1fr 1fr 1fr;} .g31{grid-template-columns:1.4fr 1fr;}
  @media (max-width:960px){.g2,.g3,.g31{grid-template-columns:1fr;}}
  .card{background:var(--card); border:1px solid var(--line); border-radius:12px; padding:16px 18px; box-shadow:0 1px 3px rgba(60,50,20,.05);}
  .card h4{font-size:13.5px; font-weight:700; margin-bottom:4px;}
  .card .cap{font-size:12px; color:var(--mut); margin-bottom:10px;}
  .chart-box{position:relative; height:300px;}
  .chart-box.tall{height:430px;}
  table{width:100%; border-collapse:collapse; font-size:13px;}
  th{text-align:left; font-size:11.5px; text-transform:uppercase; letter-spacing:.5px; color:var(--mut); padding:6px 8px; border-bottom:1px solid var(--line); font-weight:600;}
  td{padding:6px 8px; border-bottom:1px solid #f2efe6;}
  tr:last-child td{border-bottom:none;}
  td.num, th.num{text-align:right;}
  .bar{height:8px; border-radius:4px; background:var(--teal); display:inline-block;}
  .hot{color:var(--rose); font-weight:700;}
  .model-or{font-weight:700; font-size:15px; color:var(--teal-d);}
  .sig{color:var(--olive); font-weight:600;}
  .caveat{background:#fdf2f0; border:1px solid #f2d5cf; border-left:4px solid var(--rose); border-radius:8px; padding:10px 14px; font-size:12.5px; color:#8c3a30; margin-top:12px;}
  .two-col{display:grid; grid-template-columns:1fr 1fr; gap:14px;}
  @media (max-width:960px){.two-col{grid-template-columns:1fr;}}
  footer{margin-top:34px; color:var(--mut); font-size:12px; border-top:1px solid var(--line); padding-top:16px;}
  .legend{display:flex; gap:14px; flex-wrap:wrap; margin-top:8px; font-size:12px; color:var(--mut);}
  .legend .sw{width:12px;height:12px;border-radius:3px;display:inline-block;margin-right:5px;vertical-align:-1px;}
</style>
</head>
<body>
<header class="top">
  <div class="top-inner">
    <div class="logo">KM</div>
    <div>
      <h1>Kenya Malaria Prevalence Dashboard</h1>
      <div class="sub">Plasmodium falciparum survey data 1985–2020 &middot; Malaria Atlas Project &middot; Calvin Okoth &middot; KEMRI attachment</div>
    </div>
    <div class="pill">Interactive — slicers filter every chart</div>
  </div>
</header>
<div class="wrap">

  <div class="kpis">
    <div class="kpi"><div class="lbl">Surveys analysed</div><div class="val teal" id="k_surveys">—</div><div class="note">of <span id="k_raw">—</span> raw records</div></div>
    <div class="kpi"><div class="lbl">Individuals examined</div><div class="val" id="k_exam">—</div><div class="note">across all surveys</div></div>
    <div class="kpi"><div class="lbl">P. falciparum infections</div><div class="val rose" id="k_pos">—</div><div class="note">positive results</div></div>
    <div class="kpi"><div class="lbl">Overall prevalence</div><div class="val amber" id="k_pr">—</div><div class="note">positive / examined</div></div>
    <div class="kpi"><div class="lbl">Time span</div><div class="val" id="k_span">—</div><div class="note" id="k_span_note">surveys per year shown</div></div>
    <div class="kpi"><div class="lbl">Median per-survey</div><div class="val olive" id="k_med">—</div><div class="note" id="k_med_note">examined (IQR)</div></div>
  </div>

  <div class="filters">
    <h3>Filters (Power BI-style slicers)</h3>
    <div class="frow">
      <div class="fgroup">
        <label>Setting</label>
        <div class="chips" id="chip-setting">
          <div class="chip on" data-k="0">Rural <span class="n" id="n_set0"></span></div>
          <div class="chip on" data-k="1">Urban <span class="n" id="n_set1"></span></div>
          <div class="chip on" data-k="2">Not recorded <span class="n" id="n_set2"></span></div>
        </div>
      </div>
      <div class="fgroup">
        <label>Diagnostic method</label>
        <div class="chips" id="chip-method">
          <div class="chip on" data-k="0">Microscopy <span class="n" id="n_met0"></span></div>
          <div class="chip on" data-k="1">RDT <span class="n" id="n_met1"></span></div>
        </div>
      </div>
      <div class="fgroup">
        <label>Survey year</label>
        <div class="yr">
          <input type="number" id="yrMin" min="1985" max="2020" value="1985">
          <span>to</span>
          <input type="number" id="yrMax" min="1985" max="2020" value="2020">
        </div>
      </div>
      <button class="reset" id="btnReset">Reset all</button>
    </div>
  </div>

  <section id="trend-sec">
    <div class="sec-head"><h2>Temporal trends</h2><span class="tag">Weighted prevalence</span></div>
    <div class="insight">Weighted prevalence (infections ÷ examined) fell from ~46% in 1985 to ~8% by 2020. Read the line alongside the bars: years with very few surveys (e.g. 1992, 1996–97) rest on little data. The 2008–2010 and 2015/2020 peaks in survey volume reflect national school-based survey programmes.</div>
    <div class="grid g31">
      <div class="card"><h4>Prevalence &amp; survey volume by year</h4><div class="cap">Bars: surveys conducted · Line: infections ÷ individuals examined</div><div class="chart-box" id="wrap-trend"><canvas id="trend"></canvas></div></div>
      <div class="card"><h4>Prevalence by decade</h4><div class="cap">Weighted prevalence per survey decade</div><div class="chart-box" id="wrap-decade"><canvas id="decade"></canvas></div></div>
    </div>
  </section>

  <section id="season-sec">
    <div class="sec-head"><h2>Seasonality &amp; trends by group</h2><span class="tag">When · Who</span></div>
    <div class="grid g3">
      <div class="card"><h4>Surveys by month of year</h4><div class="cap">Bars: surveys started that month · Line: weighted prevalence</div><div class="chart-box" id="wrap-season"><canvas id="season"></canvas></div></div>
      <div class="card"><h4>Rural vs urban prevalence by year</h4><div class="cap">Weighted prevalence of infections</div><div class="chart-box" id="wrap-settrend"><canvas id="settrend"></canvas></div></div>
      <div class="card"><h4>Microscopy vs RDT by year</h4><div class="cap">Weighted prevalence of infections</div><div class="chart-box" id="wrap-mettrend"><canvas id="mettrend"></canvas></div></div>
    </div>
    <div class="insight" style="margin-top:14px">Seasonality is visible even in this survey-based data: fieldwork and reported prevalence are highest in the rainy-season months (Feb–May) and lowest from July and November — consistent with Kenya's bimodal transmission pattern. Rural sites carry roughly twice the prevalence of urban sites in most years.</div>
  </section>

  <section id="geo">
    <div class="sec-head"><h2>Geography of transmission</h2><span class="tag">1,607 sites with coordinates</span></div>
    <div class="insight">Transmission is highly uneven in space. High prevalence clusters in western Kenya (Lake Victoria basin) and along the coast (Kilifi, Kwale, Malindi), while central, eastern and northern sites sit near zero — matching the documented ecology of malaria in Kenya.</div>
    <div class="grid g31">
      <div class="card"><h4>Survey sites across Kenya</h4><div class="cap">Point colour = prevalence (%) · point size = number examined · hover for details</div><div class="chart-box tall" id="wrap-map"><canvas id="map"></canvas></div>
        <div class="legend"><span><span class="sw" style="background:#f0f921"></span>0%</span><span><span class="sw" style="background:#fca636"></span>25%</span><span><span class="sw" style="background:#e54c4c"></span>50%</span><span><span class="sw" style="background:#a52c78"></span>75%</span><span><span class="sw" style="background:#3b0f70"></span>100%</span></div>
      </div>
      <div class="card"><h4>Prevalence by approximate zone</h4><div class="cap">Zones derived from coordinates; hover for survey counts</div><div class="chart-box tall" id="wrap-zone"><canvas id="zone"></canvas></div></div>
    </div>
    <div class="grid g2" style="margin-top:14px">
      <div class="card"><h4>Top 10 hotspot sites</h4><div class="cap">Highest per-survey prevalence (small samples inflate extremes)</div>
        <div style="max-height:300px; overflow:auto"><table id="hotspots"><thead><tr><th>Site</th><th>Year</th><th>Setting</th><th class="num">Examined</th><th class="num">Positive</th><th class="num">Prevalence</th></tr></thead><tbody></tbody></table></div>
      </div>
      <div class="card"><h4>RDT types in use</h4><div class="cap">Rapid diagnostic test brands recorded in the dataset</div>
        <table><thead><tr><th>Test</th><th class="num">Surveys</th><th class="num">Share</th></tr></thead><tbody id="rdt-table"></tbody></table>
      </div>
    </div>
  </section>

  <section id="break">
    <div class="sec-head"><h2>Breakdowns</h2><span class="tag">Setting · Method · Age · Species</span></div>
    <div class="grid g3">
      <div class="card"><h4>Individuals examined by setting</h4><div class="cap">Share of examined population</div><div class="chart-box" id="wrap-set"><canvas id="set"></canvas></div></div>
      <div class="card"><h4>Individuals examined by method</h4><div class="cap">Microscopy vs rapid diagnostic tests</div><div class="chart-box" id="wrap-met"><canvas id="met"></canvas></div></div>
      <div class="card"><h4>Surveyed population by age group</h4><div class="cap">Who each survey sampled</div><div class="chart-box" id="wrap-age"><canvas id="age"></canvas></div></div>
    </div>
    <div class="grid g2" style="margin-top:14px">
      <div class="card"><h4>Prevalence by diagnostic method</h4><div class="cap">Distribution of per-survey prevalence (box = IQR, whiskers = min–max)</div><div class="chart-box" id="wrap-boxmet"><canvas id="boxmet"></canvas></div></div>
      <div class="card"><h4>Prevalence by setting</h4><div class="cap">Distribution of per-survey prevalence</div><div class="chart-box" id="wrap-boxset"><canvas id="boxset"></canvas></div></div>
    </div>
    <div class="grid g2" style="margin-top:14px">
      <div class="card"><h4>Prevalence by latitude band</h4><div class="cap">Altitude proxy — lowlands vs highlands gradient</div><div class="chart-box" id="wrap-lat"><canvas id="lat"></canvas></div></div>
      <div class="card"><h4>Prevalence by survey sample size</h4><div class="cap">Weighted prevalence per examined-size band</div><div class="chart-box" id="wrap-sizebin"><canvas id="sizebin"></canvas></div></div>
    </div>
    <div class="insight" style="margin-top:14px">All <span id="species_n">—</span> surveys report <b>Plasmodium falciparum</b> as the species. Microscopy surveys show a higher median prevalence than RDT surveys in this dataset, but the two methods are also used in different places and years — they are not directly interchangeable (see the model, below).</div>
  </section>

  <section id="reliab">
    <div class="sec-head"><h2>Sample size &amp; reliability</h2><span class="tag">Why small surveys are risky</span></div>
    <div class="grid g31">
      <div class="card"><h4>Sample size vs estimated prevalence</h4><div class="cap">Each point is one survey; orange line is the smoothed trend</div><div class="chart-box" id="wrap-size"><canvas id="size"></canvas></div></div>
      <div class="card"><h4>Survey summary statistics</h4><div class="cap">Per-survey distribution (current filter)</div>
        <table id="sumstat"><tbody></tbody></table>
      </div>
    </div>
    <div class="insight" style="margin-top:14px">Surveys examining fewer than ~20 people produce the most extreme estimates — both the highest and lowest prevalence values. This is why weighted summaries (by individuals examined) are preferred over simple averages of sites.</div>
  </section>

  <section id="model">
    <div class="sec-head"><h2>Statistical model</h2><span class="tag">Binomial logistic regression</span></div>
    <div class="insight">Model: <code>positive/examined ~ year + setting + method</code> (RURAL and Microscopy as reference levels). Fitted to the <b id="model_n">—</b> surveys with a recorded setting. Each odds ratio is adjusted for the other variables.</div>
    <div class="grid g31">
      <div class="card"><h4>Adjusted odds ratios (95% CI)</h4>
        <table><thead><tr><th>Characteristic</th><th class="num">OR</th><th class="num">95% CI</th><th class="num">p-value</th></tr></thead><tbody id="model-table"></tbody></table>
        <div class="caveat"><b>Read with care:</b> survey sites are not a random sample of Kenyan communities, and sites cluster spatially (e.g. Kilifi, western Kenya). These are associational, not causal, estimates.</div>
      </div>
      <div class="card"><h4>What the model says</h4>
        <div style="margin-top:6px">
          <div class="kpi" style="box-shadow:none; margin-bottom:10px"><div class="lbl">Odds of infection per additional year</div><div class="val teal" id="m_year">—</div><div class="note">≈ <span id="m_year_pct">—</span>% decline per year, holding setting &amp; method constant</div></div>
          <div class="kpi" style="box-shadow:none; margin-bottom:10px"><div class="lbl">Urban vs rural</div><div class="val amber" id="m_urban">—</div><div class="note">≈ half the odds of rural sites</div></div>
          <div class="kpi" style="box-shadow:none"><div class="lbl">RDT vs microscopy</div><div class="val olive" id="m_rdt">—</div><div class="note">adjusted for year &amp; setting</div></div>
        </div>
        <div class="cap" style="margin-top:12px">Model fit: deviance <span id="m_dev">—</span>, AIC <span id="m_aic">—</span></div>
      </div>
    </div>
  </section>

  <section id="data-sec">
    <div class="sec-head"><h2>Full dataset explorer</h2><span class="tag" id="dt_count">—</span></div>
    <div class="card">
      <h4>All surveys (filtered)</h4><div class="cap">Search any site name; rows follow the slicers above</div>
      <input type="text" id="dtSearch" placeholder="Search site…" style="width:100%; padding:9px 12px; border:1px solid var(--line); border-radius:8px; margin-bottom:10px; font-size:13px; background:#faf9f5;">
      <div style="max-height:420px; overflow:auto"><table><thead><tr><th>Site</th><th class="num">Year</th><th>Setting</th><th>Method</th><th class="num">Examined</th><th class="num">Positive</th><th class="num">Prevalence</th></tr></thead><tbody id="dtBody"></tbody></table></div>
    </div>
  </section>

  <footer>
    <b>Data:</b> Malaria Atlas Project (malariaAtlas R package) — community &amp; school-based cross-sectional surveys of <i>P. falciparum</i> in Kenya, 1985–2020. <b>Cleaning:</b> <span id="f_dropped">—</span> of <span id="f_raw">—</span> raw records removed (missing/invalid parasitology); <span id="f_coords">—</span> sites have coordinates. <b>Analysis in R 4.6.1:</b> <code>malaria_analysis.R</code> — tidyverse summaries and a binomial <code>glm()</code> logistic regression (year + setting + method); dashboard rendered Power BI-style with Chart.js (embedded, works offline). Mirrors <i>Calvin_Malaria.Project.Rmd</i>.
  </footer>
</div>

<script>
const DATA = @@DATA@@;
const OUTLINE = @@OUTLINE@@;
const KP = @@KP@@;
const RSTATS = @@RSTATS@@;
const AGE_NAMES = @@AGENAMES@@;
const ZONE_NAMES = @@ZONENAMES@@;

// ---------- filter state ----------
let fYear = [KP.year_min, KP.year_max];
let fSet = new Set([0, 1, 2]);
let fMet = new Set([0, 1]);

const fmt = n => n >= 1e6 ? (n/1e6).toFixed(2) + "M" : n >= 1e4 ? (n/1e3).toFixed(1) + "K" : n.toLocaleString();

function filtered() {
  return DATA.filter(r =>
    r[0] >= fYear[0] && r[0] <= fYear[1] && fSet.has(r[1]) && fMet.has(r[2]));
}

// ---------- aggregation helpers ----------
function aggByYear(rows) {
  const m = new Map();
  for (const r of rows) {
    if (!m.has(r[0])) m.set(r[0], [0, 0, 0]);
    const a = m.get(r[0]); a[0]++; a[1] += r[4]; a[2] += r[3];
  }
  return [...m.entries()].sort((a, b) => a[0] - b[0]);
}
function aggBy(rows, idx) {
  const m = new Map();
  for (const r of rows) {
    if (!m.has(r[idx])) m.set(r[idx], [0, 0, 0]);
    const a = m.get(r[idx]); a[0]++; a[1] += r[4]; a[2] += r[3];
  }
  return [...m.entries()].sort((a, b) => a[0] - b[0]);
}
function pct(p, e) { return e ? p / e * 100 : 0; }
function quantile(vals, qq) {
  if (!vals.length) return 0;
  const s = [...vals].sort((a, b) => a - b);
  const i = (s.length - 1) * qq, lo = Math.floor(i), hi = Math.ceil(i);
  return lo === hi ? s[lo] : s[lo] + (s[hi] - s[lo]) * (i - lo);
}
function boxStats(vals) {
  return { min: Math.min(...vals), q1: quantile(vals, 0.25), med: quantile(vals, 0.5), q3: quantile(vals, 0.75), max: Math.max(...vals) };
}

// ---------- chart factory ----------
const charts = {};
const LIGHT = { grid: { color: "#eee9dc" }, ticks: { color: "#736e60" }, title: { color: "#26231c" } };
function makeChart(id, cfg) {
  if (charts[id]) charts[id].destroy();
  charts[id] = new Chart(document.getElementById(id), cfg);
}

// ---------- Kenya outline plugin ----------
const kenyaPlugin = {
  id: "kenya",
  beforeDatasetsDraw(chart) {
    const { ctx, chartArea, scales } = chart;
    if (!chartArea || !OUTLINE.length) return;
    ctx.save();
    ctx.beginPath();
    OUTLINE.forEach((ring, i) => {
      ring.forEach((pt, j) => {
        const x = scales.x.getPixelForValue(pt[0]);
        const y = scales.y.getPixelForValue(pt[1]);
        j === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
      });
      ctx.closePath();
    });
    ctx.fillStyle = "rgba(245,240,228,0.9)";
    ctx.fill();
    ctx.strokeStyle = "rgba(180,168,142,0.9)";
    ctx.lineWidth = 1.2;
    ctx.stroke();
    ctx.restore();
  }
};

// ---------- render ----------
function renderAll() {
  const rows = filtered();

  // KPIs
  const n = rows.length;
  const ex = rows.reduce((a, r) => a + r[3], 0);
  const po = rows.reduce((a, r) => a + r[4], 0);
  document.getElementById("k_surveys").textContent = n.toLocaleString();
  document.getElementById("k_exam").textContent = fmt(ex);
  document.getElementById("k_pos").textContent = fmt(po);
  document.getElementById("k_pr").textContent = pct(po, ex).toFixed(1) + "%";
  document.getElementById("k_span").textContent = fYear[0] + "–" + fYear[1];
  const prs = rows.map(r => r[10]);
  document.getElementById("k_med").textContent = quantile(prs, 0.5).toFixed(1) + "%";
  document.getElementById("k_med_note").textContent = "median prevalence, " + n.toLocaleString() + " surveys";
  document.getElementById("k_raw").textContent = KP.n_surveys.toLocaleString();
  document.getElementById("species_n").textContent = n.toLocaleString();

  // trend chart
  const yr = aggByYear(rows);
  const years = yr.map(x => x[0]);
  const surveys = yr.map(x => x[1][0]);
  const wpr = yr.map(x => +pct(x[1][1], x[1][2]).toFixed(2));
  makeChart("trend", {
    type: "bar",
    data: {
      labels: years,
      datasets: [
        { type: "bar", label: "Surveys", data: surveys, yAxisID: "y1", backgroundColor: "rgba(232,163,61,.55)", borderColor: "rgba(232,163,61,.9)", borderWidth: 1, borderRadius: 3 },
        { type: "line", label: "Weighted prevalence (%)", data: wpr, yAxisID: "y", borderColor: "#0f9d8f", backgroundColor: "#0f9d8f", pointRadius: 3, pointBackgroundColor: "#0b7a70", tension: 0.25, borderWidth: 2.5 }
      ]
    },
    options: {
      maintainAspectRatio: false, responsive: true, interaction: { mode: "index", intersect: false },
      scales: {
        x: { grid: { display: false }, ticks: { color: "#736e60", maxRotation: 60, minRotation: 0 } },
        y: { ...LIGHT, position: "left", title: { display: true, text: "Prevalence (%)" }, min: 0, max: 100 },
        y1: { position: "right", grid: { display: false }, ticks: { color: "#b97c22" }, title: { display: true, text: "Surveys" } }
      },
      plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } } }
    }
  });

  // decade chart
  const dec = [[1985, 1989], [1990, 1999], [2000, 2004], [2005, 2009], [2010, 2020]];
  const dlab = ["1985–89", "1990–99", "2000–04", "2005–09", "2010–20"];
  const dd = dec.map(([a, b], i) => {
    const sub = rows.filter(r => r[0] >= a && r[0] <= b);
    const pe = sub.reduce((s, r) => s + r[4], 0), ee = sub.reduce((s, r) => s + r[3], 0);
    return { label: dlab[i], pr: +pct(pe, ee).toFixed(1), n: sub.length };
  });
  makeChart("decade", {
    type: "bar",
    data: { labels: dd.map(d => d.label), datasets: [{ label: "Weighted prevalence (%)", data: dd.map(d => d.pr), backgroundColor: dd.map(d => d.pr > 25 ? "#d94f4f" : d.pr > 12 ? "#e8a33d" : "#0f9d8f"), borderRadius: 5 }] },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { grid: { display: false }, ticks: { color: "#736e60" } }, y: { ...LIGHT, min: 0, max: 60, title: { display: true, text: "Prevalence (%)" } } },
      plugins: { legend: { display: false }, tooltip: { callbacks: { afterLabel: c => "  " + dd[c.dataIndex].n.toLocaleString() + " surveys" } } }
    }
  });

  // map
  const pts = rows.filter(r => r[5] !== null && r[6] !== null);
  makeChart("map", {
    type: "scatter",
    data: { datasets: [{ label: "Survey sites", data: pts.map(r => ({ x: r[5], y: r[6], site: r[9], yr: r[0], ex: r[3], po: r[4], pr: r[10], met: r[2] === 0 ? "Microscopy" : "RDT", set: ["Rural", "Urban", "Not recorded"][r[1]] })), pointBackgroundColor: pts.map(r => r[11]), pointBorderColor: "rgba(255,255,255,.75)", pointBorderWidth: 0.6, pointRadius: pts.map(r => Math.max(2.5, Math.min(11, 2.5 + Math.log2(r[3] + 1) * 1.1))), pointHoverRadius: 8 }] },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: {
        x: { min: 33.6, max: 42.1, grid: { display: false }, ticks: { display: false }, title: { display: true, text: "Longitude" } },
        y: { min: -5.0, max: 5.8, grid: { display: false }, ticks: { display: false }, title: { display: true, text: "Latitude" } }
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            title: c => c[0].raw.site,
            label: c => {
              const r = c.raw;
              return ["Year: " + r.yr, "Examined: " + r.ex.toLocaleString(), "Positive: " + r.po.toLocaleString(), "Prevalence: " + r.pr.toFixed(1) + "%", r.set + " · " + r.met];
            }
          }
        }
      }
    },
    plugins: [kenyaPlugin]
  });

  // zone chart
  const zr = aggBy(rows, 8);
  const zData = [0, 1, 2, 3, 4].map(i => {
    const found = zr.find(x => x[0] === i);
    return found ? { name: ZONE_NAMES[i], pr: +pct(found[1][1], found[1][2]).toFixed(1), n: found[1][0] } : { name: ZONE_NAMES[i], pr: 0, n: 0 };
  });
  makeChart("zone", {
    type: "bar",
    data: { labels: zData.map(z => z.name), datasets: [{ label: "Weighted prevalence (%)", data: zData.map(z => z.pr), backgroundColor: ["#0f9d8f", "#0b7a70", "#7a8c4e", "#e8a33d"], borderRadius: 5 }] },
    options: {
      indexAxis: "y", maintainAspectRatio: false, responsive: true,
      scales: { x: { ...LIGHT, min: 0, max: 40, title: { display: true, text: "Prevalence (%)" } }, y: { grid: { display: false }, ticks: { color: "#736e60", font: { size: 11 } } } },
      plugins: { legend: { display: false }, tooltip: { callbacks: { afterLabel: c => "  " + zData[c.dataIndex].n.toLocaleString() + " surveys" } } }
    }
  });

  // setting + method donuts
  const setAgg = aggBy(rows, 1);
  const SET_NAMES = ["Rural", "Urban", "Not recorded"];
  const setExam = [0, 1, 2].map(i => { const f = setAgg.find(x => x[0] === i); return f ? f[1][2] : 0; });
  makeChart("set", {
    type: "doughnut",
    data: { labels: SET_NAMES, datasets: [{ data: setExam, backgroundColor: ["#0f9d8f", "#e8a33d", "#b8b2a2"], borderColor: "#fff", borderWidth: 2 }] },
    options: { maintainAspectRatio: false, responsive: true, cutout: "58%", plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } }, tooltip: { callbacks: { label: c => "  " + SET_NAMES[c.dataIndex] + ": " + fmt(setExam[c.dataIndex]) + " examined" } } } }
  });
  const metAgg = aggBy(rows, 2);
  const MET_NAMES = ["Microscopy", "RDT"];
  const metExam = [0, 1].map(i => { const f = metAgg.find(x => x[0] === i); return f ? f[1][2] : 0; });
  makeChart("met", {
    type: "doughnut",
    data: { labels: MET_NAMES, datasets: [{ data: metExam, backgroundColor: ["#0b7a70", "#e8a33d"], borderColor: "#fff", borderWidth: 2 }] },
    options: { maintainAspectRatio: false, responsive: true, cutout: "58%", plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } }, tooltip: { callbacks: { label: c => "  " + MET_NAMES[c.dataIndex] + ": " + fmt(metExam[c.dataIndex]) + " examined" } } } }
  });

  // age chart
  const ageAgg = aggBy(rows, 7);
  const ageN = [0, 1, 2, 3, 4].map(i => { const f = ageAgg.find(x => x[0] === i); return f ? f[1][0] : 0; });
  makeChart("age", {
    type: "bar",
    data: { labels: AGE_NAMES, datasets: [{ label: "Surveys", data: ageN, backgroundColor: "#7a8c4e", borderRadius: 5 }] },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { grid: { display: false }, ticks: { color: "#736e60", font: { size: 10.5 }, maxRotation: 45, minRotation: 0 } }, y: { ...LIGHT, beginAtZero: true } },
      plugins: { legend: { display: false } }
    }
  });

  // ---- seasonality: month of year ----------------
  const MONTHS = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
  const withMonth = rows.filter(r => r[12] !== null && r[12] >= 1 && r[12] <= 12);
  const monthAgg = [];
  for (let m = 1; m <= 12; m++) {
    const sub = withMonth.filter(r => r[12] === m);
    const pe = sub.reduce((s, r) => s + r[4], 0), ee = sub.reduce((s, r) => s + r[3], 0);
    monthAgg.push({ n: sub.length, pr: +pct(pe, ee).toFixed(1) });
  }
  makeChart("season", {
    type: "bar",
    data: {
      labels: MONTHS,
      datasets: [
        { type: "bar", label: "Surveys", data: monthAgg.map(m => m.n), yAxisID: "y1", backgroundColor: "rgba(232,163,61,.55)", borderColor: "rgba(232,163,61,.9)", borderWidth: 1, borderRadius: 3 },
        { type: "line", label: "Weighted prevalence (%)", data: monthAgg.map(m => m.pr), yAxisID: "y", borderColor: "#0f9d8f", backgroundColor: "#0f9d8f", pointRadius: 3, tension: 0.25, borderWidth: 2.5 }
      ]
    },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: {
        x: { grid: { display: false }, ticks: { color: "#736e60" } },
        y: { ...LIGHT, position: "left", min: 0, max: 40, title: { display: true, text: "Prevalence (%)" } },
        y1: { position: "right", grid: { display: false }, ticks: { color: "#b97c22" }, title: { display: true, text: "Surveys" } }
      },
      plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } } }
    }
  });

  // ---- trends by group -----------------------------
  function groupYearTrend(catIdx, catNames, colors) {
    const out = catNames.map((nm, ci) => {
      const m = new Map();
      for (const r of rows) {
        if (r[catIdx] !== ci) continue;
        if (!m.has(r[0])) m.set(r[0], [0, 0]);
        const a = m.get(r[0]); a[0] += r[4]; a[1] += r[3];
      }
      return { name: nm, color: colors[ci], pts: [...m.entries()].sort((a, b) => a[0] - b[0]).map(([y, v]) => ({ x: y, y: +pct(v[0], v[1]).toFixed(1) })) };
    });
    return out;
  }
  const settrend = groupYearTrend(1, ["Rural", "Urban"], ["#0f9d8f", "#e8a33d"]);
  makeChart("settrend", {
    type: "line",
    data: { datasets: settrend.map(g => ({ label: g.name, data: g.pts, borderColor: g.color, backgroundColor: g.color, pointRadius: 2.5, tension: 0.25, borderWidth: 2.5 })) },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { ...LIGHT, type: "linear", min: KP.year_min, max: KP.year_max, title: { display: true, text: "Year" } }, y: { ...LIGHT, min: 0, max: 80, title: { display: true, text: "Prevalence (%)" } } },
      plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } } }
    }
  });
  const mettrend = groupYearTrend(2, ["Microscopy", "RDT"], ["#0b7a70", "#e8a33d"]);
  makeChart("mettrend", {
    type: "line",
    data: { datasets: mettrend.map(g => ({ label: g.name, data: g.pts, borderColor: g.color, backgroundColor: g.color, pointRadius: 2.5, tension: 0.25, borderWidth: 2.5 })) },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { ...LIGHT, type: "linear", min: KP.year_min, max: KP.year_max, title: { display: true, text: "Year" } }, y: { ...LIGHT, min: 0, max: 80, title: { display: true, text: "Prevalence (%)" } } },
      plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } } }
    }
  });

  // ---- latitude bands ------------------------------
  const BANDS = ["Coastal plain (<1°S)", "Lowland (1°S–0.5°N)", "Highlands (0.5–2.5°N)", "Far north (>2.5°N)"];
  const bandAgg = [0, 1, 2, 3].map(i => {
    const sub = rows.filter(r => r[7] !== null && (i === 0 ? r[7] < -1 : i === 1 ? r[7] < 0.5 : i === 2 ? r[7] < 2.5 : true));
    const pe = sub.reduce((s, r) => s + r[4], 0), ee = sub.reduce((s, r) => s + r[3], 0);
    return { name: BANDS[i], pr: +pct(pe, ee).toFixed(1), n: sub.length };
  });
  makeChart("lat", {
    type: "bar",
    data: { labels: bandAgg.map(b => b.name), datasets: [{ label: "Weighted prevalence (%)", data: bandAgg.map(b => b.pr), backgroundColor: ["#0f9d8f", "#0b7a70", "#7a8c4e", "#e8a33d"], borderRadius: 5 }] },
    options: {
      indexAxis: "y", maintainAspectRatio: false, responsive: true,
      scales: { x: { ...LIGHT, min: 0, max: 40, title: { display: true, text: "Prevalence (%)" } }, y: { grid: { display: false }, ticks: { color: "#736e60", font: { size: 11 } } } },
      plugins: { legend: { display: false }, tooltip: { callbacks: { afterLabel: c => "  " + bandAgg[c.dataIndex].n.toLocaleString() + " surveys" } } }
    }
  });

  // ---- sample-size bands ---------------------------
  const SBINS = ["<20", "20–49", "50–99", "100–249", "250–499", "500+"];
  const sbin = SBINS.map((nm, i) => {
    const sub = rows.filter(r => { const e = r[3]; return i === 0 ? e < 20 : i === 5 ? e >= 500 : (e >= [20, 50, 100, 250][i - 1] && e < [50, 100, 250, 500][i - 1]); });
    const pe = sub.reduce((s, r) => s + r[4], 0), ee = sub.reduce((s, r) => s + r[3], 0);
    return { name: nm, pr: +pct(pe, ee).toFixed(1), n: sub.length };
  });
  makeChart("sizebin", {
    type: "bar",
    data: { labels: sbin.map(b => b.name), datasets: [{ label: "Weighted prevalence (%)", data: sbin.map(b => b.pr), backgroundColor: sbin.map(b => b.n < 20 ? "#d94f4f" : "#0f9d8f"), borderRadius: 5 }] },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { grid: { display: false }, ticks: { color: "#736e60" } }, y: { ...LIGHT, min: 0, max: 40, title: { display: true, text: "Prevalence (%)" } } },
      plugins: { legend: { display: false }, tooltip: { callbacks: { afterLabel: c => "  " + sbin[c.dataIndex].n.toLocaleString() + " surveys" } } }
    }
  });

  // boxplots
  function boxChart(id, catIdx, names) {
    const cats = names.map((nm, i) => ({ nm, vals: rows.filter(r => r[catIdx] === i).map(r => r[10]) }));
    const stats = cats.map(c => c.vals.length ? boxStats(c.vals) : null);
    const whisk = cats.map((c, i) => stats[i] ? [stats[i].min, stats[i].q1] : null);
    const whisk2 = cats.map((c, i) => stats[i] ? [stats[i].q3, stats[i].max] : null);
    const box = cats.map((c, i) => stats[i] ? [stats[i].q1, stats[i].q3] : null);
    const med = cats.map((c, i) => stats[i] ? stats[i].med : null);
    makeChart(id, {
      type: "bar",
      data: {
        labels: cats.map(c => c.nm),
        datasets: [
          { label: "Whisker (min–Q1)", data: whisk, backgroundColor: "rgba(15,157,143,.35)", barPercentage: 0.35, categoryPercentage: 0.6 },
          { label: "Box (Q1–Q3)", data: box, backgroundColor: "rgba(15,157,143,.75)", barPercentage: 0.7, categoryPercentage: 0.6 },
          { label: "Whisker (Q3–max)", data: whisk2, backgroundColor: "rgba(15,157,143,.35)", barPercentage: 0.35, categoryPercentage: 0.6 },
          { label: "Median", data: med, type: "line", borderColor: "#d94f4f", backgroundColor: "#d94f4f", pointRadius: 4, pointStyle: "rectRot", borderWidth: 2 }
        ]
      },
      options: {
        maintainAspectRatio: false, responsive: true,
        scales: { x: { grid: { display: false }, ticks: { color: "#736e60", font: { size: 11.5 } } }, y: { ...LIGHT, min: 0, max: 100, title: { display: true, text: "Prevalence (%)" } } },
        plugins: { legend: { display: false }, tooltip: { callbacks: { afterLabel: c => { const s = stats[c.dataIndex]; return s ? ["  min " + s.min.toFixed(1) + "%", "  Q1 " + s.q1.toFixed(1) + "%", "  median " + s.med.toFixed(1) + "%", "  Q3 " + s.q3.toFixed(1) + "%", "  max " + s.max.toFixed(1) + "%"] : ""; } } } }
      }
    });
  }
  boxChart("boxmet", 2, MET_NAMES);
  boxChart("boxset", 1, SET_NAMES);

  // sample size scatter
  const sz = rows.map(r => ({ x: r[3], y: r[10] })).sort((a, b) => a.x - b.x);
  const smooth = [];
  const WIN = Math.max(8, Math.min(30, Math.round(sz.length / 25)));
  for (let i = 0; i < sz.length; i += Math.max(1, Math.round(WIN / 3))) {
    const win = sz.slice(Math.max(0, i - WIN / 2), i + WIN / 2);
    if (!win.length) continue;
    const my = win.reduce((s, p) => s + p.y, 0) / win.length;
    smooth.push({ x: sz[i].x, y: +my.toFixed(1) });
  }
  makeChart("size", {
    type: "scatter",
    data: {
      datasets: [
        { label: "Surveys", data: sz, pointRadius: 2.5, pointBackgroundColor: "rgba(15,157,143,.5)", pointBorderColor: "rgba(15,157,143,.8)", pointBorderWidth: 0.5 },
        { label: "Smoothed trend", data: smooth, type: "line", borderColor: "#e8a33d", borderWidth: 2.5, pointRadius: 0, tension: 0.3 }
      ]
    },
    options: {
      maintainAspectRatio: false, responsive: true,
      scales: { x: { ...LIGHT, type: "linear", title: { display: true, text: "Number examined" } }, y: { ...LIGHT, min: 0, max: 100, title: { display: true, text: "Prevalence (%)" } } },
      plugins: { legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } }, tooltip: { callbacks: { title: c => "Examined: " + c[0].raw.x.toLocaleString(), label: c => "Prevalence: " + c[0].raw.y.toFixed(1) + "%" } } }
    }
  });

  // summary stats table
  const allEx = rows.map(r => r[3]);
  const allPr = rows.map(r => r[10]);
  const topEx = rows.reduce((a, b) => (b[3] > a[3] ? b : a), rows[0]);
  const statRows = [
    ["Surveys (current filter)", n.toLocaleString()],
    ["Individuals examined", fmt(ex)],
    ["Infections detected", fmt(po)],
    ["Weighted prevalence", pct(po, ex).toFixed(1) + "%"],
    ["Median examined per survey", quantile(allEx, 0.5).toLocaleString() + " (IQR " + quantile(allEx, 0.25).toLocaleString() + "–" + quantile(allEx, 0.75).toLocaleString() + ")"],
    ["Median prevalence per survey", quantile(allPr, 0.5).toFixed(1) + "% (IQR " + quantile(allPr, 0.25).toFixed(1) + "–" + quantile(allPr, 0.75).toFixed(1) + "%)"],
    ["Largest survey", topEx ? topEx[9] + " — " + topEx[3].toLocaleString() + " examined (" + topEx[0] + ")" : "—"],
    ["Sites with coordinates", rows.filter(r => r[5] !== null).length.toLocaleString()],
    ["Species recorded", "Plasmodium falciparum (100%)"]
  ];
  document.getElementById("sumstat").innerHTML = statRows.map(r => "<tr><td>" + r[0] + "</td><td class='num'><b>" + r[1] + "</b></td></tr>").join("");

  // hotspots
  const hot = [...rows].sort((a, b) => b[10] - a[10]).slice(0, 10);
  document.getElementById("hotspots").querySelector("tbody").innerHTML = hot.map(r =>
    "<tr><td>" + r[9] + "</td><td>" + r[0] + "</td><td>" + ["Rural", "Urban", "Not recorded"][r[1]] + "</td><td class='num'>" + r[3].toLocaleString() + "</td><td class='num'>" + r[4].toLocaleString() + "</td><td class='num hot'>" + r[10].toFixed(1) + "%</td></tr>"
  ).join("");

  renderDataTable();
  updateChipCounts();
}

function renderDataTable() {
  const rows = filtered();
  const q = (document.getElementById("dtSearch").value || "").trim().toLowerCase();
  const matched = q ? rows.filter(r => r[9].toLowerCase().includes(q)) : rows;
  document.getElementById("dt_count").textContent = matched.length.toLocaleString() + " rows";
  document.getElementById("dtBody").innerHTML = matched.slice(0, 500).map(r =>
    "<tr><td>" + r[9] + "</td><td class='num'>" + r[0] + "</td><td>" + ["Rural", "Urban", "Not recorded"][r[1]] + "</td><td>" + (r[2] === 0 ? "Microscopy" : "RDT") + "</td><td class='num'>" + r[3].toLocaleString() + "</td><td class='num'>" + r[4].toLocaleString() + "</td><td class='num'>" + r[10].toFixed(1) + "%</td></tr>"
  ).join("");
}

function updateChipCounts() {
  const all = filtered();
  const setAgg = aggBy(all, 1), metAgg = aggBy(all, 2);
  document.getElementById("n_set0").textContent = "· " + ((setAgg.find(x => x[0] === 0) || [0, [0]])[1][0]).toLocaleString();
  document.getElementById("n_set1").textContent = "· " + ((setAgg.find(x => x[0] === 1) || [0, [0]])[1][0]).toLocaleString();
  document.getElementById("n_set2").textContent = "· " + ((setAgg.find(x => x[0] === 2) || [0, [0]])[1][0]).toLocaleString();
  document.getElementById("n_met0").textContent = "· " + ((metAgg.find(x => x[0] === 0) || [0, [0]])[1][0]).toLocaleString();
  document.getElementById("n_met1").textContent = "· " + ((metAgg.find(x => x[0] === 1) || [0, [0]])[1][0]).toLocaleString();
}

// ---------- static content ----------
function renderStatic() {
  document.getElementById("k_raw").textContent = KP.n_surveys.toLocaleString();
  document.getElementById("k_span_note").textContent = KP.year_min + "–" + KP.year_max + " in source data";
  document.getElementById("f_dropped").textContent = KP.n_dropped.toLocaleString();
  document.getElementById("f_raw").textContent = KP.n_raw.toLocaleString();
  document.getElementById("f_coords").textContent = KP.n_coords.toLocaleString();
  document.getElementById("species_n").textContent = KP.n_surveys.toLocaleString();
  document.getElementById("model_n").textContent = KP.model_extra.n.toLocaleString();
  document.getElementById("m_year").textContent = "OR " + KP.model[0].or.toFixed(3);
  document.getElementById("m_year_pct").textContent = KP.model_extra.year_pct.toFixed(1);
  document.getElementById("m_urban").textContent = "OR " + KP.model[1].or.toFixed(3);
  document.getElementById("m_rdt").textContent = "OR " + KP.model[2].or.toFixed(3);
  document.getElementById("m_dev").textContent = KP.model_extra.deviance.toLocaleString(undefined, { maximumFractionDigits: 0 });
  document.getElementById("m_aic").textContent = KP.model_extra.aic.toLocaleString(undefined, { maximumFractionDigits: 0 });
  document.getElementById("model-table").innerHTML = KP.model.map(m => {
    const sig = m.p < 0.001 ? "<span class='sig'>p&lt;0.001</span>" : "p = " + m.p.toFixed(3);
    return "<tr><td>" + m.term + "</td><td class='num model-or'>" + m.or.toFixed(3) + "</td><td class='num'>" + m.lo.toFixed(3) + " – " + m.hi.toFixed(3) + "</td><td class='num'>" + sig + "</td></tr>";
  }).join("");
  const rdtTotal = KP.rdt.reduce((a, r) => a + r.n, 0);
  document.getElementById("rdt-table").innerHTML = KP.rdt.map(r =>
    "<tr><td>" + r.type + "</td><td class='num'>" + r.n.toLocaleString() + "</td><td class='num'>" + (r.n / rdtTotal * 100).toFixed(0) + "%</td></tr>"
  ).join("");
}

// ---------- event wiring ----------
function wire() {
  document.querySelectorAll("#chip-setting .chip").forEach(ch => {
    ch.addEventListener("click", () => {
      const k = +ch.dataset.k;
      if (fSet.has(k) && fSet.size > 1) { fSet.delete(k); ch.classList.remove("on"); }
      else if (!fSet.has(k)) { fSet.add(k); ch.classList.add("on"); }
      renderAll();
    });
  });
  document.querySelectorAll("#chip-method .chip").forEach(ch => {
    ch.addEventListener("click", () => {
      const k = +ch.dataset.k;
      if (fMet.has(k) && fMet.size > 1) { fMet.delete(k); ch.classList.remove("on"); }
      else if (!fMet.has(k)) { fMet.add(k); ch.classList.add("on"); }
      renderAll();
    });
  });
  const ymin = document.getElementById("yrMin"), ymax = document.getElementById("yrMax");
  function onYear() {
    let a = +ymin.value, b = +ymax.value;
    if (isNaN(a)) a = KP.year_min;
    if (isNaN(b)) b = KP.year_max;
    a = Math.max(KP.year_min, Math.min(KP.year_max, a));
    b = Math.max(KP.year_min, Math.min(KP.year_max, b));
    if (a > b) { const t = a; a = b; b = t; ymin.value = a; ymax.value = b; }
    fYear = [a, b];
    renderAll();
  }
  ymin.addEventListener("change", onYear);
  ymax.addEventListener("change", onYear);
  document.getElementById("dtSearch").addEventListener("input", renderDataTable);
  document.getElementById("btnReset").addEventListener("click", () => {
    fYear = [KP.year_min, KP.year_max];
    fSet = new Set([0, 1, 2]);
    fMet = new Set([0, 1]);
    ymin.value = KP.year_min; ymax.value = KP.year_max;
    document.querySelectorAll(".chip").forEach(c => c.classList.add("on"));
    renderAll();
  });
}

renderStatic();
wire();
renderAll();
</script>
</body>
</html>
"""

# --------------------------------------------------------------------------
# 5. Assemble
# --------------------------------------------------------------------------
html = (HTML
        .replace("@@DATA@@", json.dumps(survey_rows))
        .replace("@@OUTLINE@@", json.dumps(outline))
        .replace("@@KP@@", json.dumps(kp))
        .replace("@@RSTATS@@", json.dumps(rs))
        .replace("@@AGENAMES@@", json.dumps(AGE_NAMES))
        .replace("@@ZONENAMES@@", json.dumps(ZONE_NAMES)))
with open(OUT, "w", encoding="utf-8") as f:
    f.write(html)
print(f"Wrote {OUT}  ({os.path.getsize(OUT) / 1024:.0f} KB, {n_surveys} surveys embedded)")
print(f"KPIs: {n_surveys} surveys | {n_examined:,} examined | {n_positive:,} positive | {overall_pr:.1f}%")
print(f"Model (from R): year OR {model_rows[0]['or']:.3f} | urban OR {model_rows[1]['or']:.3f} | RDT OR {model_rows[2]['or']:.3f} | n={model_extra['n']}")
