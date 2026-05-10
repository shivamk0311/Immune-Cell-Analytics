from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.io as pio
from flask import Flask, render_template, request
from scipy import stats

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

app = Flask(__name__, template_folder=str(ROOT / 'templates'))

def conn():
    if not DB_PATH.exists():
        raise RuntimeError(f"{DB_PATH} not found. Run 'python load_data.py' first")
    return sqlite3.connect(DB_PATH)

def fig_html(fig, div_id):
    return pio.to_html(fig, include_plotlyjs=True, full_html=False, div_id=div_id)

def distinct(col_table, col):
    with conn() as c:
        rows = c.execute(
            f"SELECT DISTINCT {col} FROM {col_table} WHERE {col} IS NOT NULL ORDER BY {col};"
        ).fetchall()
    
    return [r[0] for r in rows]


@app.route("/")
def index():
    with conn() as c:
        counts = {
            t: c.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
            for t in ("projects", "subjects", "samples", "cell_counts")
        }
    return render_template("index.html", counts = counts)

@app.route("/frequencies")
def frequencies():
    sample_type = request.args.get("sample_type", "")
    condition = request.args.get("condition", "")

    where, params = ["1=1"], []
    if sample_type:
        where.append("s.sample_type = ?"); params.append(sample_type)
    
    if condition:
        where.append("sub.condition = ?"); params.append(condition)
    
    sql = f"""
    WITH totals AS (
        SELECT sample_id, SUM(count) AS total_count FROM cell_counts GROUP BY sample_id
    )
    SELECT cc.sample_id AS sample,
           t.total_count,
           cc.population,
           cc.count,
           ROUND(100.0 * cc.count / t.total_count, 2) AS percentage
    FROM cell_counts cc
    JOIN samples  s   ON s.sample_id  = cc.sample_id
    JOIN subjects sub ON sub.subject_id = s.subject_id
    JOIN totals   t   ON t.sample_id  = cc.sample_id
    WHERE {' AND '.join(where)}
    ORDER BY cc.sample_id, cc.population
    LIMIT 5000;
    """

    with conn() as c:
        df = pd.read_sql_query(sql, c, params=params)

    fig = px.box(
        df, x="population", y="percentage", points=False,
        title="Population relative frequency distribution"
            + (f" — {sample_type}" if sample_type else "")
            + (f" / {condition}" if condition else ""),
    )
    fig.update_layout(yaxis_title="relative frequency (%)", height=420)

    return render_template(
        "frequencies.html",
        table_html=df.head(200).to_html(classes="data", index=False, border=0),
        rows=len(df),
        plot=fig_html(fig, "freq-plot"),
        sample_types=distinct("samples", "sample_type"),
        conditions=distinct("subjects", "condition"),
        selected_sample_type=sample_type,
        selected_condition=condition,
    )

@app.route("/responders")
def responders():
    treatment = request.args.get("treatment", "miraclib")
    condition = request.args.get("condition", "melanoma")
    sample_type = request.args.get("sample_type", "PBMC")

    sql = """
    WITH totals AS (
        SELECT sample_id, SUM(count) AS total_count FROM cell_counts GROUP BY sample_id
    )
    SELECT s.sample_id      AS sample,
           sub.response      AS response,
           cc.population     AS population,
           100.0 * cc.count / t.total_count AS percentage
    FROM samples s
    JOIN subjects sub  ON sub.subject_id = s.subject_id
    JOIN cell_counts cc ON cc.sample_id  = s.sample_id
    JOIN totals t       ON t.sample_id   = s.sample_id
    WHERE sub.condition = ?
      AND sub.treatment = ?
      AND s.sample_type = ?
      AND sub.response IN ('yes', 'no');
    """
    with conn() as c:
        df = pd.read_sql_query(sql, c, params=(condition, treatment, sample_type))

    if df.empty:
        return render_template(
            "responders.html",
            empty=True, plot=None, stats_table=None,
            treatments=distinct("subjects", "treatment"),
            conditions=distinct("subjects", "condition"),
            sample_types=distinct("samples", "sample_type"),
            selected_treatment=treatment, selected_condition=condition,
            selected_sample_type=sample_type,
        )

    fig = px.box(
        df, x="population", y="percentage", color="response",
        category_orders={"population": POPULATIONS, "response": ["yes", "no"]},
        title=f"{condition} + {treatment} + {sample_type} — responders vs non-responders",
        points="all",
    )
    fig.update_layout(yaxis_title="relative frequency (%)", height=480, boxmode="group")

    rows = []
    for pop in POPULATIONS:
        sub = df[df.population == pop]
        r = sub.loc[sub.response == "yes", "percentage"].to_numpy()
        n = sub.loc[sub.response == "no",  "percentage"].to_numpy()
        if len(r) == 0 or len(n) == 0:
            rows.append({"population": pop, "n_yes": len(r), "n_no": len(n),
                         "mean_yes": None, "mean_no": None,
                         "mean_diff": None, "p_value": None, "significant": False})
            continue
        u, p = stats.mannwhitneyu(r, n, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_yes": int(len(r)), "n_no": int(len(n)),
            "mean_yes": round(float(r.mean()), 3),
            "mean_no":  round(float(n.mean()), 3),
            "mean_diff": round(float(r.mean() - n.mean()), 3),
            "p_value":  f"{p:.3g}",
            "significant": bool(p < 0.05),
        })
    stats_df = pd.DataFrame(rows)

    return render_template(
        "responders.html",
        empty=False,
        plot=fig_html(fig, "resp-plot"),
        stats_table=stats_df.to_html(classes="data", index=False, border=0),
        treatments=distinct("subjects", "treatment"),
        conditions=distinct("subjects", "condition"),
        sample_types=distinct("samples", "sample_type"),
        selected_treatment=treatment,
        selected_condition=condition,
        selected_sample_type=sample_type,
    )

@app.route("/subset")
def subset():
    treatment = request.args.get("treatment", "miraclib")
    condition = request.args.get("condition", "melanoma")
    sample_type = request.args.get("sample_type", "PBMC")
    time_point = int(request.args.get("time", 0))

    base = """
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE sub.condition = ?
      AND sub.treatment = ?
      AND s.sample_type = ?
      AND s.time_from_treatment_start = ?
    """
    params = (condition, treatment, sample_type, time_point)

    with conn() as c:
        samples = pd.read_sql_query(
            f"SELECT s.sample_id, sub.subject_id, sub.project_id, sub.response, sub.sex {base};",
            c, params=params,
        )
        avg_b = pd.read_sql_query(f"""
            SELECT AVG(cc.count) AS avg_b_cell, COUNT(*) AS n
            FROM samples s
            JOIN subjects sub  ON sub.subject_id = s.subject_id
            JOIN cell_counts cc ON cc.sample_id  = s.sample_id
            WHERE sub.condition = ? AND sub.treatment = ? AND s.sample_type = ?
              AND s.time_from_treatment_start = ?
              AND sub.sex = 'M' AND sub.response = 'yes'
              AND cc.population = 'b_cell';
        """, c, params=params)

    subj = samples.drop_duplicates("subject_id")
    by_project = samples.groupby("project_id").size().rename("n_samples").reset_index()
    by_response = subj.groupby("response", dropna=False).size().rename("n_subjects").reset_index()
    by_sex = subj.groupby("sex", dropna=False).size().rename("n_subjects").reset_index()

    avg_b_value = (round(float(avg_b.iloc[0]["avg_b_cell"]), 2)
                   if avg_b.iloc[0]["n"] else None)

    return render_template(
        "subset.html",
        treatments=distinct("subjects", "treatment"),
        conditions=distinct("subjects", "condition"),
        sample_types=distinct("samples", "sample_type"),
        selected_treatment=treatment,
        selected_condition=condition,
        selected_sample_type=sample_type,
        selected_time=time_point,
        n_samples=len(samples), n_subjects=len(subj),
        by_project=by_project.to_html(classes="data", index=False, border=0),
        by_response=by_response.to_html(classes="data", index=False, border=0),
        by_sex=by_sex.to_html(classes="data", index=False, border=0),
        avg_b_value=avg_b_value,
        avg_b_n=int(avg_b.iloc[0]["n"]),
    )

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)