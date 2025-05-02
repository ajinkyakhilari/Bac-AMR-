#!/usr/bin/env python3
"""
depth_plot.py: Generate a per-base depth file from a BAM and create per-contig depth plots.
Usage:
  python scripts/depth_plot.py -i mapping/sample.sorted.bam \
                              -o depth/sample.depth.txt \
                              -p depth/plots/
"""
import argparse
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Generate depth file and plots per contig from a BAM file"
    )
    parser.add_argument(
        "-i", "--bam", required=True,
        help="Input sorted BAM file"
    )
    parser.add_argument(
        "-o", "--depth", required=True,
        help="Output depth text file (tsv with contig, pos, depth)"
    )
    parser.add_argument(
        "-p", "--plot-dir", required=True,
        help="Directory where per-contig PNG plots will be saved"
    )
    args = parser.parse_args()

    # ensure output directories exist
    os.makedirs(os.path.dirname(args.depth), exist_ok=True)
    os.makedirs(args.plot_dir, exist_ok=True)

    # run samtools depth
    depth_cmd = ["samtools", "depth", "-a", args.bam]
    with open(args.depth, "w") as out_f:
        subprocess.check_call(depth_cmd, stdout=out_f)

    # read depth file
    df = pd.read_csv(
        args.depth,
        sep="\t",
        names=["contig", "pos", "depth"]
    )

    # plot each contig
    for contig, grp in df.groupby("contig"):
        fig, ax = plt.subplots()
        ax.plot(grp.pos, grp.depth)
        ax.set(
            title=f"{os.path.basename(args.bam)} – {contig}",
            xlabel="Position",
            ylabel="Depth"
        )
        out_png = os.path.join(args.plot_dir, f"{contig}.png")
        fig.savefig(out_png)
        plt.close(fig)

if __name__ == "__main__":
    main()
