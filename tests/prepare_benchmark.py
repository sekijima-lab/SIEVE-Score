"""Prepare numeric-only inputs from all rows of the supplied MUV CSVs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from interaction_io import load_interaction


def prepare(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for target in ('CAG', 'FXI', 'HIV'):
        source = ROOT / 'dataset/muv' / target / 'glide-dock_best_pv.interaction'
        names, labels, features, X = load_interaction(source)
        X = X.astype(np.float32)  # Both forest releases internally use float32.
        y = (labels > 0).astype(np.int64)
        np.savez_compressed(output / (target + '.npz'), X=X, y=y,
                            names=names.astype('U'))
        manifest.append(dict(target=target, path=str(source.relative_to(ROOT)),
            sha256=hashlib.sha256(source.read_bytes()).hexdigest(), rows=len(y),
            positives=int(y.sum()), features=X.shape[1],
            scope='All supplied compounds; only the CSV header is excluded.'))
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(prepare(parser.parse_args().output), indent=2))
