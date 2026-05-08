from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent 
DB_PATH = ROOT / "loblaw.db"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

def connect():
    if not DB_PATH.exists():
        raise SystemExit(f"{DB_PATH} not found — run `python load_data.py` first.")
    return sqlite3.connect(DB_PATH)

def frequencies(conn):
    sql = """
    WITH totals AS (
        SELECT sample_id, sum(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    ) 

    SELECT cc.sample_id AS sample,
            t.total_count AS total_count,
            cc.population AS population,
            cc.count AS count,
            ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
    FROM cell_counts cc
    JOIN totals t ON t.sample_id = cc.sample_id
    ORDER BY cc.sample_id, cc.population;
    """
    
    df = pd.read_sql_query(sql, conn)
    df.to_csv(OUT/"frequencies.csv", index=False)

    return df 


def responder_analysis(conn):

    sql = """
    WITH totals AS (
        SELECT sample_id, sum(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )

    SELECT s.sample_id AS sample,
           sub.response AS response,
           cc.population AS population,
           100 * cc.count/t.total_count AS percentage
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject.id
    JOIN cell_counts cc ON cc.sample_id = s.sample_id
    JOIN totals t ON t.sample_id = s.sample_id
    WHERE sub.condition = 'melanoma'
    AND   sub.treatment = 'miraclib'
    AND   s.sample_type = 'PBMC'
    AND   sub.response IN ('yes', 'no');
    """

    df = pd.read_sql_query(sql, conn)
    df.to_csv(OUT/ "responder_frequencies.csv", index=False)

    rows = []
    for pop in POPULATIONS:
        sub = df[df.population == pop]
        responders = sub.loc[sub.response == 'yes', "percentage"].to_numpy()
        non_responders = sub.loc[sub.response == 'no', "percentage"].to_numpy()

        if len(responders) == 0 or len(non_responders) == 0:
            rows.append((pop, len(responders), len(non_responders), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')))
            continue
        
        u, p = stats.mannwhitneyu(responders, non_responders, alternative = 'two-sided')

        rows.append((
            pop,
            len(responders),
            len(non_responders),
            float(responders.mean()),
            float(non_responders.mean()),
            float( responders.mean() - non_responders.mean()),
            float(u),
            float(p),
        ))

        stats_df = pd.DataFrame(rows, columns = [
            "population", "n_responders", "n_non_responders", "mean_pct_responders", "mean_pct_non_responders", "mean_diff", "mannwhitney_u", "p_value",
        ])

        stats_df["p_value_bonferroni"] = (stats_df["p_value"] * len(POPULATIONS)).clip(upper=1.0)
        stats_df["significant_alpha_0.05"] = stats_df["p_value"] < 0.05
        stats_df.to_csv(OUT / "responder_stats.csv", index=False)

        fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 4), sharey=False)
        for ax, pop in zip(axes, POPULATIONS):
            sub = df[df.population == pop]
            data = [
                sub.loc[sub.response == "yes", "percentage"].to_numpy(),
                sub.loc[sub.response == "no",  "percentage"].to_numpy(),
            ]
            ax.boxplot(data, showmeans=True)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["yes", "no"])
            ax.set_title(pop)
            ax.set_xlabel("response")
            ax.set_ylabel("relative frequency (%)")
            p = stats_df.loc[stats_df.population == pop, "p_value"].iloc[0]
            ax.text(0.5, 0.95, f"p = {p:.3g}", transform=ax.transAxes,
                    ha="center", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.8))
        fig.suptitle("Melanoma + miraclib + PBMC — responders (yes) vs non-responders (no)")
        fig.tight_layout()
        plot_path = OUT / "responder_boxplots.png"
        fig.savefig(plot_path, dpi=130)
        plt.close(fig)

        return df, stats_df, plot_path



