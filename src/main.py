import os
import argparse

import pandas as pd
import pandera.pandas as pa
from simulation import sweep, verify_distance_preservation
from visualization import render_figures, render_diagrams, print_summary

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--max-shots", type=int, default=1_000_000)
    ap.add_argument("--target-errors", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=max(os.cpu_count() - 4, 1))
    ap.add_argument("--from-data", default=None)
    return ap.parse_args()

sample_schema = pa.DataFrameSchema({
    "code": pa.Column(str),
    "eta": pa.Column(float),
    "d": pa.Column(int),
    "p": pa.Column(float),
    "pl": pa.Column(float),
    "sigma": pa.Column(float),
    "errors": pa.Column(int),
    "shots": pa.Column(int),
})

if __name__ == "__main__":
    args = parse_args()

    os.makedirs(f"{args.outdir}", exist_ok=True)
    figures_outdir = f"{args.outdir}/figures"
    diagrams_outdir = f"{args.outdir}/diagrams"

    df: pd.DataFrame = None
    
    if args.from_data is not None:
        # Load from CSV and plot
        df = pd.read_csv(args.from_data)
        df = sample_schema.validate(df)
    else:
        # Verify, sweep, save result and plot
        verify_distance_preservation()
        df = sweep(
            max_shots=args.max_shots,
            target_errors=args.target_errors,
            num_workers=args.workers,
        )
        df.to_csv(f"{args.outdir}/samples.csv", index=False)
        print(f"☑ Samples saved to {args.outdir}/samples.csv")

    render_figures(df, figures_outdir)
    render_diagrams(diagrams_outdir)
    print_summary(df)
