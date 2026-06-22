#!/usr/bin/env python3
"""Compute cell type similarity matrix based on second-best label proportions."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Paths
tuning_dir = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning/",
)

probability_csv = (
    tuning_dir / "label_spreading_diagnosis_knn_alpha_0.2_n_5_probabilities.csv"
)

output_dir = tuning_dir / "celltype_second_best_ambiguity"
output_dir.mkdir(parents=True, exist_ok=True)


def load_probability_matrix(
    probability_csv: Path,
    high_conf_threshold: float = 0.8,
    margin_threshold: float = 0.8,
) -> pd.Dataframe:
    """Load probability matrix and check that probabilities sum to 1 per cell."""
    # Load probability matrix
    prob_df = pd.read_csv(
        probability_csv,
        index_col=0,
    )

    prob_df = prob_df.apply(pd.to_numeric, errors="coerce")
    prob_df = prob_df.dropna(axis=1, how="all")
    prob_df = prob_df.fillna(0)

    # check that probabilities sum to 1 per cell
    row_sums = prob_df.sum(axis=1)
    if np.allclose(row_sums.median(), 1.0, atol=0.05):
        min_possible_margin = 2 * high_conf_threshold - 1

        print("\nProbabilities appear to approximately sum to 1 per cell.")
        print(
            f"With top_score >= {high_conf_threshold}, "
            f"the smallest possible margin is approximately {min_possible_margin:.3f}.",
        )

        if margin_threshold < min_possible_margin:
            msg = (
                f"Therefore, no high-confidence cells can also have "
                f"margin <= {margin_threshold}."
            )
            raise ValueError(msg)
    else:
        msg = "Probabilities do not appear to sum to 1 per cell."
        raise ValueError(msg)

    return prob_df


def get_top_2_labels(
    prob_df: pd.DataFrame,
    high_conf_threshold: float = 0.8,
    second_score_threshold: float = 0.2,
    margin_threshold: float = 0.8,
) -> pd.DataFrame:
    """Get top and second-best label per cell, with flags for ambiguous cells."""
    # Get top and second-best label per cell
    labels = np.array(prob_df.columns)

    scores = prob_df.to_numpy(dtype=float)
    scores_for_ranking = np.nan_to_num(scores, nan=-np.inf)

    order = np.argsort(scores_for_ranking, axis=1)[:, ::-1]

    top_idx = order[:, 0]
    second_idx = order[:, 1]

    top_label = labels[top_idx]
    second_label = labels[second_idx]

    top_score = scores[np.arange(scores.shape[0]), top_idx]
    second_score = scores[np.arange(scores.shape[0]), second_idx]

    assignment_df = pd.DataFrame(
        {
            "top_label": top_label,
            "top_score": top_score,
            "second_label": second_label,
            "second_score": second_score,
        },
        index=prob_df.index,
    )

    assignment_df["margin"] = assignment_df["top_score"] - assignment_df["second_score"]

    assignment_df.to_csv(output_dir / "celltype_top_second_assignment_all_cells.csv")

    assignment_df["ambiguous_low_intermediate"] = (
        (assignment_df["top_score"] < high_conf_threshold)
        & (assignment_df["second_score"] >= second_score_threshold)
        & (assignment_df["margin"] <= margin_threshold)
    )

    assignment_df["high_confidence"] = assignment_df["top_score"] >= high_conf_threshold

    assignment_df["high_confidence_with_credible_second"] = (
        assignment_df["top_score"] >= high_conf_threshold
    ) & (assignment_df["second_score"] >= second_score_threshold)

    assignment_df.to_csv(output_dir / "celltype_top_second_assignment_with_flags.csv")

    return assignment_df


def _plot_second_best_heatmap(props, output_file, title):
    """Plot second-best proportion heatmap."""
    plot_df = props.copy()

    # Drop labels absent from this subset
    plot_df = plot_df.loc[
        plot_df.sum(axis=1) > 0,
        plot_df.sum(axis=0) > 0,
    ]

    if plot_df.empty:
        print(f"Skipping empty heatmap: {output_file.name}")
        return

    plt.figure(figsize=(12, 10))

    sns.heatmap(
        plot_df,
        cmap="viridis",
        square=True,
        linewidths=0.5,
        cbar_kws={"label": "Proportion where column label is second-best"},
    )

    plt.title(title)
    plt.xlabel("Second-best label")
    plt.ylabel("Highest-confidence label")
    plt.tight_layout()

    plt.savefig(
        output_file,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )

    plt.close()


def _build_pair_review_table(counts, props, global_spearman):
    """
    Build one row per unordered label pair.

    a_top_b_second_prop:
        Among cells whose top label is A, proportion with B as second-best.

    b_top_a_second_prop:
        Among cells whose top label is B, proportion with A as second-best.

    reciprocal_ambiguity_score:
        min(A -> B, B -> A)
    """
    records = []

    labels_present = sorted(set(counts.index).union(set(counts.columns)))

    for i, label_a in enumerate(labels_present):
        for label_b in labels_present[i + 1 :]:
            if label_a == label_b:
                continue

            a_to_b = (
                props.loc[label_a, label_b]
                if label_a in props.index and label_b in props.columns
                else 0.0
            )

            b_to_a = (
                props.loc[label_b, label_a]
                if label_b in props.index and label_a in props.columns
                else 0.0
            )

            a_to_b_count = (
                counts.loc[label_a, label_b]
                if label_a in counts.index and label_b in counts.columns
                else 0
            )

            b_to_a_count = (
                counts.loc[label_b, label_a]
                if label_b in counts.index and label_a in counts.columns
                else 0
            )

            n_a = counts.loc[label_a].sum() if label_a in counts.index else 0

            n_b = counts.loc[label_b].sum() if label_b in counts.index else 0

            corr = (
                global_spearman.loc[label_a, label_b]
                if label_a in global_spearman.index
                and label_b in global_spearman.columns
                else np.nan
            )

            records.append(
                {
                    "label_a": label_a,
                    "label_b": label_b,
                    "n_a_cells": n_a,
                    "n_b_cells": n_b,
                    "a_top_b_second_count": a_to_b_count,
                    "b_top_a_second_count": b_to_a_count,
                    "a_top_b_second_prop": a_to_b,
                    "b_top_a_second_prop": b_to_a,
                    "reciprocal_ambiguity_score": min(a_to_b, b_to_a),
                    "mean_ambiguity_score": (a_to_b + b_to_a) / 2,
                    "global_spearman": corr,
                },
            )

    review_df = pd.DataFrame.from_records(records)

    if review_df.empty:
        return review_df

    return review_df.sort_values(
        [
            "reciprocal_ambiguity_score",
            "mean_ambiguity_score",
            "global_spearman",
        ],
        ascending=False,
    )


def run_second_best_analysis(
    prob_df,
    assignment_df,
    subset_mask,
    analysis_name,
    title_prefix,
    review_cutoff=0.1,
):
    """Generate counts, proportions, heatmap, and review tables."""
    subset_df = assignment_df[subset_mask].copy()

    subset_df.to_csv(output_dir / f"{analysis_name}_cell_assignments.csv")

    print(f"\n{analysis_name}")
    print(f"Cells in subset: {subset_df.shape[0]:,}")

    expected_review_cols = [
        "label_a",
        "label_b",
        "n_a_cells",
        "n_b_cells",
        "a_top_b_second_count",
        "b_top_a_second_count",
        "a_top_b_second_prop",
        "b_top_a_second_prop",
        "reciprocal_ambiguity_score",
        "mean_ambiguity_score",
        "global_spearman",
    ]

    if subset_df.empty:
        print("No cells in this subset.")

        empty_review = pd.DataFrame(columns=expected_review_cols)

        empty_review.to_csv(
            output_dir / f"{analysis_name}_pair_review.csv",
            index=False,
        )

        return empty_review

    counts = pd.crosstab(
        subset_df["top_label"],
        subset_df["second_label"],
    )

    all_labels = list(prob_df.columns)

    counts = counts.reindex(
        index=all_labels,
        columns=all_labels,
        fill_value=0,
    )

    counts.to_csv(output_dir / f"{analysis_name}_second_best_counts.csv")

    row_totals = counts.sum(axis=1).replace(0, np.nan)

    props = counts.div(
        row_totals,
        axis=0,
    ).fillna(0)

    props.to_csv(output_dir / f"{analysis_name}_second_best_proportions.csv")

    _plot_second_best_heatmap(
        props=props,
        output_file=output_dir / f"{analysis_name}_second_best_heatmap.pdf",
        title=f"{title_prefix}\nP(second-best label | top label)",
    )

    # Global Spearman correlation between confidence profiles
    global_spearman = prob_df.corr(method="spearman")
    global_spearman.to_csv(output_dir / "global_celltype_spearman_correlation.csv")

    review_df = _build_pair_review_table(
        counts=counts,
        props=props,
        global_spearman=global_spearman,
    )

    review_df.to_csv(
        output_dir / f"{analysis_name}_pair_review.csv",
        index=False,
    )

    if not review_df.empty:
        candidates = review_df[
            review_df["reciprocal_ambiguity_score"] >= review_cutoff
        ].copy()

        candidates.to_csv(
            output_dir
            / f"{analysis_name}_curation_candidates_reciprocal_ge_{review_cutoff}.csv",
            index=False,
        )

    return None


def main() -> None:
    """Run second-best label analysis for ambiguous and high-confidence cells."""
    prob_df = load_probability_matrix(probability_csv)

    assignment_df = get_top_2_labels(prob_df)

    # Ambiguous low/intermediate-confidence cells
    run_second_best_analysis(
        assignment_df=assignment_df,
        subset_mask=assignment_df["ambiguous_low_intermediate"],
        analysis_name="ambiguous_low_intermediate",
        title_prefix=(
            "Second-best labels among ambiguous low/intermediate-confidence cells"
        ),
    )

    # High-confidence cells with credible second-best label
    run_second_best_analysis(
        assignment_df=assignment_df,
        subset_mask=assignment_df["high_confidence_with_credible_second"],
        analysis_name="high_confidence_credible_second",
        title_prefix=(
            "Second-best labels among high-confidence cells with credible alternatives"
        ),
    )

    # High-confidence cells only, regardless of second-best score
    run_second_best_analysis(
        assignment_df=assignment_df,
        subset_mask=assignment_df["high_confidence"],
        analysis_name="high_confidence_all_second",
        title_prefix=("Second-best labels among all high-confidence cells"),
    )


if __name__ == "__main__":
    main()
