from __future__ import annotations
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "loblaw.db"

POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte" )

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition TEXT,
    age INTEGER,
    sex TEXT,
    treatment TEXT,
    response TEXT
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT,
    time_from_treatment_start INTEGER
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    population TEXT NOT NULL,
    count INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_condition ON subjects(condition);
CREATE INDEX idx_subjects_treatment ON subjects(treatment);
CREATE INDEX idx_subjects_response ON subjects(response);
CREATE INDEX idx_subjects_subject ON subjects(subject_id);
CREATE INDEX idx_subjects_time ON samples(time_from_treatment_start);
CREATE INDEX idx_subjects_type ON samples(sample_type);
CREATE INDEX idx_subjects_pop ON cell_counts(population);

"""


def init_db (db_path):
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    conn.executescript(SCHEMA)
    return conn 


def _to_int(s):
    s = (s or "").strip()
    return int(s) if s else None 

def load_csv(conn, csv_path):
    projects = set()
    subjects = {}
    samples = []
    counts = []

    with csv_path.open(newline = "") as f:
        for row in csv.DictReader(f):
            projects.add(row["project"])
            subjects.setdefault(row["subject"],
                                (
                                    row["subject"],
                                    row["project"],
                                    row["condition"] or None,
                                    _to_int(row["age"]),
                                    row["sex"] or None,
                                    row["treatment"] or None, 
                                    row["response"] or None,
                                ),
            )

            samples.append(
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"] or None,
                    _to_int(row["time_from_treatment_start"]),
                )
            )

            for pop in POPULATIONS:
                counts.append((row["sample"], pop, int(row[pop])))
    
    curr = conn.cursor()
    curr.executemany("INSERT INTO projects VALUES (?)", [(p,) for p in sorted(projects)])
    curr.executemany(
        "INSERT INTO subjects VALUES (?,?,?,?,?,?,?)",
        list(subjects.values()),
    )
    curr.executemany("INSERT INTO samples VALUES (?,?,?,?)", samples)
    curr.executemany("INSERT INTO cell_counts VALUES (?,?,?)", counts)
    conn.commit()

    return {
        "projects": len(projects),
        "subjects": len(subjects),
        "samples": len(samples),
        "cell_counts": len(counts),
    }

def main():
    if not CSV_PATH.exists():
        raise SystemExit(f"cell-count.csv not found at {CSV_PATH}")

    conn = init_db(DB_PATH)

    try:
        stats = load_csv(conn, CSV_PATH)
    finally:
        conn.close()
    
    print(f"Wrote {DB_PATH.name}:")
    for table, n in stats.items():
        print(f"  {table:12s}  {n:>8d} rows")

if __name__ == "__main__":
    main()