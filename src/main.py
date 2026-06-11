import os
import csv
import argparse

from simulation import sweep, verify_distance_preservation
from threshold import SamplePoint
from visualization import render_all, render_diagrams
from rich.status import Status

CSV_FIELDS = ["code", "eta", "d", "p", "pl", "sigma", "errors", "shots"]

def save_samples(path: str, rows: list[dict]) -> None:
    with Status("Saving samples", spinner="arc"):
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        print(f"☑ Samples saved to {path}")

def row_to_pair(row: dict) -> tuple[str, SamplePoint]:
    sp = SamplePoint(
        eta=float(row["eta"]),
        distance=int(row["d"]),
        physical_error_rate=float(row["p"]),
        logical_error_rate=float(row["pl"]),
        standard_deviation=float(row["sigma"]),
        logical_errors=int(row["errors"]),
        shots=int(row["shots"]),
    )
    return row["code"], sp

def load_samples(path: str) -> list[tuple[str, SamplePoint]]:    
    with open(path, newline="") as f:
        return [row_to_pair(row) for row in csv.DictReader(f)]

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--max-shots", type=int, default=2_000_000)
    ap.add_argument("--target-errors", type=int, default=200)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--from-data", default=None)
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()

    os.makedirs(f"{args.outdir}", exist_ok=True)
    figures_outdir = f"{args.outdir}/figures"
    diagrams_outdir = f"{args.outdir}/diagrams"
    
    if args.from_data is not None:
        # Load from CSV and plot
        pairs = load_samples(args.from_data)
        render_all(pairs, figures_outdir)
    else:
        # Verify, sweep, save result and plot
        verify_distance_preservation()
        rows = sweep(
            max_shots=args.max_shots,
            target_errors=args.target_errors,
            num_workers=args.workers,
        )
        save_samples(f"{args.outdir}/samples.csv", rows)
        render_all([row_to_pair(row) for row in rows], figures_outdir)

    # Generate lattice and circuit diagrams
    render_diagrams(diagrams_outdir)
