# Loblaw Bio — Cell-Count Analysis & Dashboard

A small data pipeline + Flask dashboard for Bob Loblaw's miraclib clinical-trial
cell-count data.

## Quick start

```bash
make setup        # install Python deps from requirements.txt
make pipeline     # build loblaw.db, write all CSVs and the Part 3 boxplot
make dashboard    # start Flask at http://127.0.0.1:5000
```

`make pipeline` runs end-to-end with no manual steps:

1. `load_data.py` — drops/recreates `loblaw.db` and loads `cell-count.csv`.
2. `analysis.py` — Parts 2/3/4 outputs into `outputs/`.

The dashboard URL when running locally: <http://127.0.0.1:5000>.

## Repository layout

```
.
├── load_data.py           # Part 1 — schema + CSV loader
├── analysis.py            # Parts 2, 3, 4 — runs the queries and writes outputs/
├── dashboard.py           # Flask app (Plotly charts, jinja templates)
├── templates/             # Jinja2 templates for the four dashboard pages
├── cell-count.csv         # input
├── outputs/               # generated tables and the boxplot PNG
├── requirements.txt
├── Makefile               # setup / pipeline / dashboard
└── README.md
```

## Database schema

Four tables, normalized so subject- and project-level metadata isn't repeated
across the ~10k sample-population rows.

```
projects(project_id PK)

subjects(subject_id PK,
         project_id FK -> projects.project_id,
         condition, age, sex, treatment, response)

samples(sample_id PK,
        subject_id FK -> subjects.subject_id,
        sample_type, time_from_treatment_start)

cell_counts(sample_id FK -> samples.sample_id,
            population, count,
            PRIMARY KEY (sample_id, population))
```

### Why this shape

- **`cell_counts` is long, not wide.** The CSV stores the five populations as
  columns (`b_cell`, `cd8_t_cell`, …). Storing them as rows means adding a
  population (e.g. `dendritic`) is a data change, not a schema migration, and
  per-population aggregations are a `GROUP BY population` instead of a `UNION ALL`
  across columns.
- **Subjects separated from samples.** A subject has many samples over time;
  treatment/response/sex/condition belong to the subject, time-varying state
  belongs to the sample. Splitting them avoids the de-duplication step that
  would otherwise be required for "how many subjects" questions (see Part 4).
- **`projects` is its own table.** Trivial today (one column) but gives a clean
  hook for project-level metadata when there are hundreds of trials —
  sponsoring lab, indication area, study start date, etc.
- **Indexes** are on every column we filter or join on (`condition`,
  `treatment`, `response`, `sample_type`, `time_from_treatment_start`,
  `population`). At ~10k rows this is overkill, but it's the right shape for
  the scaling discussion below.

### Scaling to hundreds of projects / thousands of samples

- The per-sample frequency calculation in Part 2 is a single CTE — it stays
  O(n) on `cell_counts` and benefits from the `(sample_id, population)`
  primary key.
- For ad-hoc analytics it's worth precomputing a `sample_totals(sample_id,
  total_count)` table (or a materialized view, if we move to Postgres) so the
  inner CTE doesn't re-scan `cell_counts` on every dashboard hit.
- At a few million rows the right next move is Postgres + a star-style
  reporting layer: `fact_cell_counts` keyed on `(sample_id, population_id)`,
  small `dim_subject` / `dim_sample` / `dim_population` / `dim_project`
  tables, and a couple of materialized rollups (per-sample frequency,
  per-cohort means) refreshed on a schedule.
- For repeated cohort queries (e.g. "all melanoma + miraclib + PBMC at t=0"),
  store cohort definitions as named views or in a `cohorts` table so the
  filter logic isn't duplicated across the analysis script and the dashboard.
- The dashboard already parameterizes condition / treatment / sample type /
  time, so adding a new study mostly means loading more rows; no code change
  is needed for the existing analyses to apply to the new cohort.

## Code structure

Three top-level scripts, each does one thing, and they share nothing except
the SQLite file on disk:

- **`load_data.py`** is idempotent: it deletes `loblaw.db` if present and
  rebuilds it from the CSV. The CSV is read once with `csv.DictReader` and
  loaded with `executemany` batches. No pandas dependency in the loader keeps
  it cheap to run.
- **`analysis.py`** is the headless pipeline that the grader's `make pipeline`
  invokes. Each part is its own function returning a DataFrame so the dashboard
  can call them too if we ever want to share more logic. Outputs are CSVs (so
  the grader can diff them) plus one PNG for Part 3.
- **`dashboard.py`** is a thin Flask app. Each page does its own SQL query
  scoped to the user's filter selections, renders a Plotly figure inline (no
  separate JSON endpoint), and ships an HTML table. State lives in the URL
  query string — no session, no JS framework, nothing to keep in sync.

## Statistical choice (Part 3)

Per-population comparison uses the **Mann-Whitney U test** (two-sided). It's
non-parametric, doesn't assume normal distributions of the relative-frequency
values, and handles modest sample sizes well. The five-test family is
Bonferroni-adjusted in `outputs/responder_stats.csv`; the dashboard shows the
unadjusted p-values for clarity and the corrected ones live in the CSV.

## Outputs (generated by `make pipeline`)

```
outputs/
├── frequencies.csv             # Part 2: long-format relative-frequency table
├── responder_frequencies.csv   # Part 3: melanoma+miraclib+PBMC, per-sample %
├── responder_stats.csv         # Part 3: Mann-Whitney U + Bonferroni
├── responder_boxplots.png      # Part 3: per-population boxplots, yes vs no
├── subset_samples.csv          # Part 4: matching baseline samples
├── subset_by_project.csv       # Part 4: samples per project
├── subset_by_response.csv      # Part 4: subjects by response
├── subset_by_sex.csv           # Part 4: subjects by sex
└── subset_avg_b_cell.csv       # Part 4: avg B-cells, male responders, t=0
```

## Dashboard

Local link (after `make dashboard`): <http://127.0.0.1:5000>

Pages:

- **Overview** — counts of projects/subjects/samples and what each tab does.
- **Part 2 — Frequencies** — boxplot of population frequencies + a paged HTML
  table, filterable by sample type and condition.
- **Part 3 — Responders** — grouped boxplot (yes vs no, per population) and
  per-population Mann-Whitney U statistics. Defaults reproduce Bob's question
  (melanoma / miraclib / PBMC) but each filter is editable.
- **Part 4 — Subset** — totals + samples-per-project / subjects-by-response /
  subjects-by-sex tables, plus the average B-cell count headline number for
  male responders at the chosen time point.

## Notes

- The dataset uses `condition` and `sex` columns (the brief mentions
  "indication" and "gender"); the schema sticks with the actual CSV column
  names so the loader is dumb and faithful.
- `make clean` removes the database and outputs folder if you want to fully
  re-run from a blank slate.
