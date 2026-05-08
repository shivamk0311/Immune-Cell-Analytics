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

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

app = Flask(__name__, template_folder=str(ROOT / 'templates'))

def conn():
    if not DB_PATH.exists():
        raise RuntimeError(f"{DB_PATH} not found. Run 'python load_data.py' first")
    return sqlite3.connect(DB_PATH)

def fig_html(fig, div_id):
    return pio.to_html(fig, include_plotlyjs=False, full_html=False, div_id=div_id)

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


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)