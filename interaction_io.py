"""Read precomputed interaction CSVs without importing Schrödinger."""
from pathlib import Path
import csv

import numpy as np
import pandas as pd


def load_interaction(f_name, hits=None, ignore=None):
    path = Path(f_name).with_suffix(".interaction")
    if path.exists():
        with path.open(newline="") as source:
            header = next(csv.reader(source), [])
        if len(header) != len(set(header)):
            raise ValueError("Interaction CSV needs unique feature names.")
        frame = pd.read_csv(path, index_col=0)
    else:
        # Conversion remains an optional, separate Schrödinger operation.
        try:
            from read_interaction import read_interaction
        except ImportError as exc:
            raise FileNotFoundError(
                f"Precomputed interaction file not found: {path}. "
                "Create it in the separate Schrödinger environment first."
            ) from exc
        raw = read_interaction(f_name, hits)
        frame = pd.DataFrame(raw[1:], columns=raw[0]).set_index(raw[0][0])

    if frame.empty or "ishit" not in frame or "docking_score" not in frame:
        raise ValueError("Interaction CSV needs compound rows, ishit and docking_score.")
    if frame.columns.has_duplicates or frame.index.hasnans or frame.index.has_duplicates:
        raise ValueError("Interaction CSV needs unique feature names and compound titles.")
    # Preserve supplied feature order, with docking_score last as expected by the model.
    names = [name for name in frame.columns if name not in ("ishit", "docking_score")]
    if not names:
        raise ValueError("Interaction CSV has no residue interaction features.")
    names.append("docking_score")
    if ignore is not None:
        ignored = Path(ignore).read_text().splitlines()
        frame = frame.loc[~frame.index.astype(str).isin(ignored)]
    if frame.empty:
        raise ValueError("No compounds remain after filtering.")
    labels = pd.to_numeric(frame["ishit"], errors="raise").to_numpy(dtype=np.float64)
    values = frame[names].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.float64)
    if not np.isfinite(values).all() or not np.isfinite(labels).all():
        raise ValueError("Interaction features and labels must be finite numeric values.")
    return frame.index.astype(str).to_numpy(), labels, names, values
