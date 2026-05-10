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
    JOIN subjects sub ON sub.subject_id = s.subject_id
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



def part4_subset(conn):
    base_filter = """
        sub.condition = 'melanoma'
    AND sub.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
    """

    samples = pd.read_sql_query(f"""
        SELECT s.sample_id, sub.subject_id, sub.project_id,
               sub.response, sub.sex, sub.age
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE {base_filter};
    """, conn)
    samples.to_csv(OUT / "subset_samples.csv", index=False)

    by_project = (samples.groupby("project_id")
                          .size().rename("n_samples").reset_index()
                          .sort_values("project_id"))
    by_project.to_csv(OUT / "subset_by_project.csv", index=False)

    
    subj = samples.drop_duplicates("subject_id")
    by_response = (subj.groupby("response", dropna=False)
                       .size().rename("n_subjects").reset_index())
    by_response.to_csv(OUT / "subset_by_response.csv", index=False)

    by_sex = (subj.groupby("sex", dropna=False)
                  .size().rename("n_subjects").reset_index())
    by_sex.to_csv(OUT / "subset_by_sex.csv", index=False)

    avg_b = pd.read_sql_query(f"""
        SELECT AVG(cc.count) AS avg_b_cell, COUNT(*) AS n
        FROM samples s
        JOIN subjects sub  ON sub.subject_id = s.subject_id
        JOIN cell_counts cc ON cc.sample_id  = s.sample_id
        WHERE {base_filter}
          AND sub.sex      = 'M'
          AND sub.response = 'yes'
          AND cc.population = 'b_cell';
    """, conn)
    avg = float(avg_b.iloc[0]["avg_b_cell"]) if avg_b.iloc[0]["n"] else float("nan")
    n = int(avg_b.iloc[0]["n"])

    summary = {
        "n_samples_total": int(len(samples)),
        "n_subjects_total": int(subj.shape[0]),
        "by_project": by_project.to_dict(orient="records"),
        "by_response": by_response.to_dict(orient="records"),
        "by_sex": by_sex.to_dict(orient="records"),
        "avg_b_cell_male_responder_t0": round(avg, 2) if n else None,
        "n_male_responder_t0_samples": n,
    }
    pd.Series(
        {"avg_b_cell_male_responder_t0": summary["avg_b_cell_male_responder_t0"],
         "n_male_responder_t0_samples":  summary["n_male_responder_t0_samples"]}
    ).to_csv(OUT / "subset_avg_b_cell.csv", header=["value"])
    return summary


def main() -> None:
    conn = connect()
    try:
        freq = frequencies(conn)
        print(f"[part 2] frequencies: {len(freq):,} rows -> outputs/frequencies.csv")

        _, stats_df, plot_path = responder_analysis(conn)
        print(f"[part 3] stats:")
        for _, r in stats_df.iterrows():
            sig = "*" if r["significant_alpha_0.05"] else " "
            print(f"   {sig} {r['population']:11s}  p={r['p_value']:.4g}  "
                  f"diff={r['mean_diff']:+.2f}%")
        print(f"[part 3] boxplots -> {plot_path.relative_to(ROOT)}")

        s = part4_subset(conn)
        print(f"[part 4] baseline samples: {s['n_samples_total']} "
              f"({s['n_subjects_total']} subjects)")
        print(f"[part 4] by project: {s['by_project']}")
        print(f"[part 4] by response: {s['by_response']}")
        print(f"[part 4] by sex: {s['by_sex']}")
        print(f"[part 4] avg b_cell (male, responder, t=0) = "
              f"{s['avg_b_cell_male_responder_t0']}  (n={s['n_male_responder_t0_samples']})")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
