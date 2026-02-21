#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- CONFIG ----------------
NULL_CSV = "external/pcrEfficiency/Results/null_scores_Geller.csv"
NULL_COL = "avg_efficiency"

METHOD_SCORES = {
    "PD-single-LPath": 1.910455,
    "Álvarez-Rodríguez et al.": 1.861111,
}

OUT_FIG = "Results/null_distribution_CVB3.png"
# ----------------------------------------


def load_null_scores(csv_path, col):
    vals = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                vals.append(float(r[col]))
            except:
                continue
    return np.array(vals, dtype=float)


def empirical_p_value(null, score):
    return (1.0 + np.sum(null >= score)) / float(1.0 + len(null))


def plot_null_with_methods(null, method_scores, out_path):

    # Larger, paper-appropriate size
    fig, ax = plt.subplots(figsize=(7, 5))

    # Histogram
    ax.hist(
        null,
        bins=25,
        density=True,
        color="#bdbdbd",
        edgecolor="black",
        linewidth=0.7,
        alpha=0.85,
        label="Random design"
    )

    colors = {
        "PD-single-LPath": "#d62728",
        "Álvarez-Rodríguez et al.": "#1f77b4",
    }

    # Vertical lines
    for method, score in method_scores.items():
        pval = empirical_p_value(null, score)

        ax.axvline(
            score,
            color=colors[method],
            linestyle="--",
            linewidth=2.5,
            label=method
        )

        ymax = ax.get_ylim()[1]
        ax.text(
            score+0.008,
            ymax * 1.2,
            "p = %.3g" % pval,
            rotation=90,
            color=colors[method],
            fontsize=11,
            ha="center",
            va="top",
            fontweight="bold"
        )

    # Labels
    ax.set_xlabel("Average predicted PCR efficiency", fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.set_ylim(0,11)

    ax.tick_params(labelsize=12)

    # Legend inside top-left, non-overlapping
    ax.legend(
        loc="upper left",
        fontsize=11,
        frameon=True,
        facecolor="white",
        edgecolor="black"
    )

    fig.tight_layout()

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():

    null = load_null_scores(NULL_CSV, NULL_COL)

    plot_null_with_methods(null, METHOD_SCORES, OUT_FIG)

    print("Saved figure to:", OUT_FIG)


if __name__ == "__main__":
    main()