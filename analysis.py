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

