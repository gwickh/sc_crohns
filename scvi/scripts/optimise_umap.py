from __future__ import annotations

import inspect
import math
import re
import time
from itertools import product
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from scipy import sparse
from sklearn.manifold import trustworthiness
from sklearn.metrics import silhouette_score

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------


def safe_name(value: object) -> str:
    """Make a string safe for use in filenames."""
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    return text.strip("_")


def dense_array(matrix) -> np.ndarray:
    """Convert a sparse or dense matrix to a dense NumPy array."""
    if sparse.issparse(matrix):
        return matrix.toarray()

    return np.asarray(matrix)


def run_neighbors(
    adata: ad.AnnData,
    n_neighbors: int,
    use_rep: str,
    metric: str,
    random_seed: int,
) -> None:
    """
    Run Scanpy neighbours while supporting Scanpy versions that use either
    `rng` or `random_state`.
    """
    kwargs = {
        "n_neighbors": n_neighbors,
        "use_rep": use_rep,
        "metric": metric,
    }

    parameters = inspect.signature(sc.pp.neighbors).parameters

    if "rng" in parameters:
        kwargs["rng"] = random_seed
    else:
        kwargs["random_state"] = random_seed

    sc.pp.neighbors(adata, **kwargs)


def run_umap(
    adata: ad.AnnData,
    min_dist: float,
    spread: float,
    random_seed: int,
) -> None:
    """
    Run Scanpy UMAP while supporting Scanpy versions that use either
    `rng` or `random_state`.
    """
    kwargs = {
        "min_dist": min_dist,
        "spread": spread,
        "init_pos": "spectral",
    }

    parameters = inspect.signature(sc.tl.umap).parameters

    if "rng" in parameters:
        kwargs["rng"] = random_seed
    else:
        kwargs["random_state"] = random_seed

    sc.tl.umap(adata, **kwargs)


def make_colour_dictionary(labels: list[str]) -> dict[str, str]:
    """Assign one consistent colour to every curated cell-type label."""
    if len(labels) <= len(sc.pl.palettes.default_102):
        colours = list(sc.pl.palettes.default_102[: len(labels)])
    else:
        cmap = plt.colormaps["hsv"].resampled(len(labels))
        colours = [cmap(index) for index in range(len(labels))]

    return dict(zip(labels, colours, strict=True))


def calculate_silhouette(
    embedding: np.ndarray,
    labels: np.ndarray,
) -> float:
    """
    Calculate a silhouette score when it is mathematically possible.

    This score is reported as a diagnostic but is not used to select the
    best UMAP, because optimising against curated labels can make the
    visualisation look artificially well separated.
    """
    labels = np.asarray(labels, dtype=str)
    n_groups = len(np.unique(labels))

    if n_groups < 2 or n_groups >= len(labels):
        return np.nan

    try:
        return float(silhouette_score(embedding, labels))
    except ValueError:
        return np.nan


def plot_umap_by_diagnosis(
    embedding: np.ndarray,
    obs: pd.DataFrame,
    diagnosis_column: str,
    label_column: str,
    diagnosis_order: list[str],
    colour_dictionary: dict[str, str],
    title: str,
    output_path: Path,
    random_seed: int,
    point_size: float,
    alpha: float,
    dpi: int,
) -> None:
    """
    Make one UMAP panel per diagnosis.

    All panels use the same axis limits and the same cell-type colours.
    """
    n_panels = len(diagnosis_order)

    figure_width = max(6.0, 5.0 * n_panels + 4.0)
    figure_height = 6.0

    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(figure_width, figure_height),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    axes = axes.ravel()
    rng = np.random.default_rng(random_seed)

    # Use identical limits for every diagnosis panel.
    x_min, x_max = embedding[:, 0].min(), embedding[:, 0].max()
    y_min, y_max = embedding[:, 1].min(), embedding[:, 1].max()

    x_padding = max((x_max - x_min) * 0.03, 0.1)
    y_padding = max((y_max - y_min) * 0.03, 0.1)

    for axis, diagnosis in zip(axes, diagnosis_order, strict=True):
        panel_mask = obs[diagnosis_column].astype(str).to_numpy() == diagnosis
        panel_indices = np.flatnonzero(panel_mask)

        # Shuffle plotting order so that one cell type does not systematically
        # cover all the others.
        panel_indices = rng.permutation(panel_indices)

        panel_labels = obs.iloc[panel_indices][label_column].astype(str).to_numpy()

        panel_colours = [colour_dictionary[label] for label in panel_labels]

        axis.scatter(
            embedding[panel_indices, 0],
            embedding[panel_indices, 1],
            c=panel_colours,
            s=point_size,
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

        axis.set_title(
            f"{diagnosis}\n{len(panel_indices):,} cells",
            fontsize=12,
        )

        axis.set_xlim(x_min - x_padding, x_max + x_padding)
        axis.set_ylim(y_min - y_padding, y_max + y_padding)
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.set_xticks([])
        axis.set_yticks([])

        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.2)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            label=label,
            markerfacecolor=colour_dictionary[label],
            markeredgecolor="none",
            markersize=6,
        )
        for label in colour_dictionary
    ]

    legend_columns = max(
        1,
        math.ceil(len(legend_handles) / 25),
    )

    fig.legend(
        handles=legend_handles,
        title=label_column,
        loc="center left",
        bbox_to_anchor=(0.99, 0.5),
        frameon=False,
        ncol=legend_columns,
        fontsize=8,
        title_fontsize=9,
    )

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 0.84, 0.94))

    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)


# ---------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------


def main() -> int:
    # -----------------------------------------------------------------
    # Input and output settings
    # -----------------------------------------------------------------

    input_path = Path(
        "project-area/data/crohns_scrnaseq/10c_14n_analysis/adata_final.h5ad",
    )

    output_dir = Path(
        "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/umap_parameter_optimisation"
    )

    diagnosis_column = "diagnosis"
    label_column = "final_curated_labels"
    use_rep = "X_embeddings"

    # Diagnosis panel order. Any additional diagnosis levels will be
    # appended after these.
    preferred_diagnosis_order = [
        "Normal",
        "Crohn's disease",
    ]

    # -----------------------------------------------------------------
    # UMAP parameter grid
    # -----------------------------------------------------------------

    n_neighbors_values = [
        50,
    ]

    min_dist_values = [
        0.5,
    ]

    spread_values = [
        1.0,
    ]

    distance_metrics = [
        "euclidean",
        # "cosine",
    ]

    # -----------------------------------------------------------------
    # Plot settings
    # -----------------------------------------------------------------

    point_size = 1.5
    alpha = 0.7
    plot_dpi = 300
    save_pdf = False

    # -----------------------------------------------------------------
    # Scoring settings
    # -----------------------------------------------------------------

    # Calculating trustworthiness on all ~180,000 cells would be expensive.
    # The same random subset is used for every parameter combination.
    score_sample_size = 2_000
    trustworthiness_n_neighbors = 15
    random_seed = 42

    # If True, write a copy of the AnnData object containing the chosen
    # UMAP in adata.obsm["X_umap_optimised"].
    write_updated_h5ad = False

    # -----------------------------------------------------------------
    # Load data
    # -----------------------------------------------------------------

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {input_path}")
    adata = sc.read_h5ad(input_path)

    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"Available obsm entries: {list(adata.obsm.keys())}")

    required_columns = {
        diagnosis_column,
        label_column,
    }

    missing_columns = required_columns.difference(adata.obs.columns)

    if missing_columns:
        raise KeyError(
            "Missing required adata.obs columns: "
            f"{sorted(missing_columns)}\n"
            f"Available columns: {list(adata.obs.columns)}",
        )

    if use_rep not in adata.obsm:
        raise KeyError(
            f"Representation {use_rep!r} was not found in adata.obsm.\n"
            f"Available representations: {list(adata.obsm.keys())}\n"
            "Change `use_rep` in main() to the correct representation.",
        )

    # Remove cells without a diagnosis.
    diagnosis_valid = adata.obs[diagnosis_column].notna() & adata.obs[
        diagnosis_column
    ].astype(str).str.strip().ne("")

    if not diagnosis_valid.all():
        n_removed = int((~diagnosis_valid).sum())
        print(f"Removing {n_removed:,} cells without a diagnosis.")
        adata = adata[diagnosis_valid].copy()

    # Retain cells without curated labels as an explicit group.
    adata.obs[label_column] = (
        adata.obs[label_column]
        .astype("string")
        .fillna("Unassigned")
        .replace("", "Unassigned")
    )

    adata.obs[diagnosis_column] = adata.obs[diagnosis_column].astype(str)

    observed_diagnoses = sorted(adata.obs[diagnosis_column].unique().tolist())

    diagnosis_order = [
        diagnosis
        for diagnosis in preferred_diagnosis_order
        if diagnosis in observed_diagnoses
    ]

    diagnosis_order.extend(
        diagnosis
        for diagnosis in observed_diagnoses
        if diagnosis not in diagnosis_order
    )

    cell_type_labels = sorted(adata.obs[label_column].astype(str).unique().tolist())

    colour_dictionary = make_colour_dictionary(cell_type_labels)

    print(f"Diagnosis order: {diagnosis_order}")
    print(f"Number of cell types: {len(cell_type_labels)}")
    print(f"Using representation: {use_rep}")

    # Save the cell-type colour mapping used in all plots.
    pd.DataFrame(
        {
            label_column: list(colour_dictionary.keys()),
            "colour": list(colour_dictionary.values()),
        },
    ).to_csv(
        output_dir / "cell_type_colour_mapping.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Select a fixed subset for scoring all parameter combinations
    # -----------------------------------------------------------------

    rng = np.random.default_rng(random_seed)

    n_score_cells = min(score_sample_size, adata.n_obs)

    score_indices = np.sort(
        rng.choice(
            adata.n_obs,
            size=n_score_cells,
            replace=False,
        ),
    )

    high_dimensional_sample = dense_array(adata.obsm[use_rep][score_indices])

    score_obs = adata.obs.iloc[score_indices]

    print(
        f"Using {n_score_cells:,} cells to calculate "
        "trustworthiness and silhouette scores.",
    )

    # -----------------------------------------------------------------
    # Run parameter grid
    # -----------------------------------------------------------------

    results = []
    best_trustworthiness = -np.inf
    best_embedding = None
    best_parameters = None

    total_combinations = (
        len(distance_metrics)
        * len(n_neighbors_values)
        * len(min_dist_values)
        * len(spread_values)
    )

    combination_number = 0

    for metric, n_neighbors in product(
        distance_metrics,
        n_neighbors_values,
    ):
        print(
            f"\nComputing neighbour graph: metric={metric}, n_neighbors={n_neighbors}",
        )

        neighbour_start = time.perf_counter()

        run_neighbors(
            adata=adata,
            n_neighbors=n_neighbors,
            use_rep=use_rep,
            metric=metric,
            random_seed=random_seed,
        )

        neighbour_seconds = time.perf_counter() - neighbour_start

        # Reuse this neighbour graph for all min_dist/spread combinations.
        for min_dist, spread in product(
            min_dist_values,
            spread_values,
        ):
            combination_number += 1

            if min_dist > spread:
                print(
                    f"Skipping min_dist={min_dist}, spread={spread}: "
                    "min_dist should not be greater than spread.",
                )
                continue

            print(
                f"\n[{combination_number}/{total_combinations}] "
                f"metric={metric}, "
                f"n_neighbors={n_neighbors}, "
                f"min_dist={min_dist}, "
                f"spread={spread}",
            )

            umap_start = time.perf_counter()

            run_umap(
                adata=adata,
                min_dist=min_dist,
                spread=spread,
                random_seed=random_seed,
            )

            umap_seconds = time.perf_counter() - umap_start

            embedding = np.asarray(
                adata.obsm["X_umap"],
                dtype=np.float64,
            ).copy()

            embedding_sample = embedding[score_indices]

            # Trustworthiness measures how well local neighbourhoods in
            # the latent space are retained in the two-dimensional UMAP.
            score_k = min(
                trustworthiness_n_neighbors,
                max(1, (n_score_cells - 1) // 2),
            )

            trustworthiness_score = trustworthiness(
                high_dimensional_sample,
                embedding_sample,
                n_neighbors=score_k,
                metric=metric,
            )

            cell_type_silhouette = calculate_silhouette(
                embedding_sample,
                score_obs[label_column].astype(str).to_numpy(),
            )

            diagnosis_silhouette = calculate_silhouette(
                embedding_sample,
                score_obs[diagnosis_column].astype(str).to_numpy(),
            )

            parameter_name = (
                f"metric_{safe_name(metric)}"
                f"__neighbors_{n_neighbors}"
                f"__min_dist_{min_dist:g}"
                f"__spread_{spread:g}"
            )

            png_path = output_dir / f"{parameter_name}.png"

            plot_title = (
                f"UMAP: {use_rep}\n"
                f"metric={metric}, "
                f"n_neighbors={n_neighbors}, "
                f"min_dist={min_dist:g}, "
                f"spread={spread:g}"
            )

            plot_umap_by_diagnosis(
                embedding=embedding,
                obs=adata.obs,
                diagnosis_column=diagnosis_column,
                label_column=label_column,
                diagnosis_order=diagnosis_order,
                colour_dictionary=colour_dictionary,
                title=plot_title,
                output_path=png_path,
                random_seed=random_seed,
                point_size=point_size,
                alpha=alpha,
                dpi=plot_dpi,
            )

            if save_pdf:
                plot_umap_by_diagnosis(
                    embedding=embedding,
                    obs=adata.obs,
                    diagnosis_column=diagnosis_column,
                    label_column=label_column,
                    diagnosis_order=diagnosis_order,
                    colour_dictionary=colour_dictionary,
                    title=plot_title,
                    output_path=output_dir / f"{parameter_name}.pdf",
                    random_seed=random_seed,
                    point_size=point_size,
                    alpha=alpha,
                    dpi=plot_dpi,
                )

            results.append(
                {
                    "metric": metric,
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist,
                    "spread": spread,
                    "trustworthiness": trustworthiness_score,
                    "cell_type_silhouette": cell_type_silhouette,
                    "diagnosis_silhouette": diagnosis_silhouette,
                    "neighbor_graph_seconds": neighbour_seconds,
                    "umap_seconds": umap_seconds,
                    "plot_file": png_path.name,
                },
            )

            print(
                f"Trustworthiness: {trustworthiness_score:.4f}; "
                f"cell-type silhouette: {cell_type_silhouette:.4f}; "
                f"diagnosis silhouette: {diagnosis_silhouette:.4f}",
            )

            if trustworthiness_score > best_trustworthiness:
                best_trustworthiness = trustworthiness_score
                best_embedding = embedding.copy()
                best_parameters = {
                    "metric": metric,
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist,
                    "spread": spread,
                    "trustworthiness": trustworthiness_score,
                    "cell_type_silhouette": cell_type_silhouette,
                    "diagnosis_silhouette": diagnosis_silhouette,
                }

    if not results or best_embedding is None:
        raise RuntimeError("No UMAP parameter combinations were completed.")

    # -----------------------------------------------------------------
    # Save ranked results
    # -----------------------------------------------------------------

    scores = pd.DataFrame(results).sort_values(
        "trustworthiness",
        ascending=False,
    )

    scores.insert(
        0,
        "rank",
        np.arange(1, len(scores) + 1),
    )

    scores.to_csv(
        output_dir / "umap_parameter_scores.csv",
        index=False,
    )

    # Save coordinates from the best-scoring UMAP.
    best_coordinates = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "UMAP_1": best_embedding[:, 0],
            "UMAP_2": best_embedding[:, 1],
            diagnosis_column: adata.obs[diagnosis_column].astype(str).to_numpy(),
            label_column: adata.obs[label_column].astype(str).to_numpy(),
        },
    )

    best_coordinates.to_csv(
        output_dir / "best_umap_coordinates.csv",
        index=False,
    )

    best_parameters_df = pd.DataFrame([best_parameters])

    best_parameters_df.to_csv(
        output_dir / "best_umap_parameters.csv",
        index=False,
    )

    # Also produce a clearly named copy of the best plot.
    best_plot_title = (
        f"Best UMAP ranked by trustworthiness: {use_rep}\n"
        f"metric={best_parameters['metric']}, "
        f"n_neighbors={best_parameters['n_neighbors']}, "
        f"min_dist={best_parameters['min_dist']:g}, "
        f"spread={best_parameters['spread']:g}, "
        f"trustworthiness={best_parameters['trustworthiness']:.4f}"
    )

    plot_umap_by_diagnosis(
        embedding=best_embedding,
        obs=adata.obs,
        diagnosis_column=diagnosis_column,
        label_column=label_column,
        diagnosis_order=diagnosis_order,
        colour_dictionary=colour_dictionary,
        title=best_plot_title,
        output_path=output_dir / "best_umap_by_diagnosis.png",
        random_seed=random_seed,
        point_size=point_size,
        alpha=alpha,
        dpi=plot_dpi,
    )

    # Optionally save the best UMAP inside a new AnnData file.
    if write_updated_h5ad:
        adata.obsm["X_umap_optimised"] = best_embedding

        adata.uns["umap_optimisation"] = {
            key: (value.item() if isinstance(value, np.generic) else value)
            for key, value in best_parameters.items()
        }

        adata.write_h5ad(output_dir / "adata_with_optimised_umap.h5ad")

    print("\nBest parameters:")
    for key, value in best_parameters.items():
        print(f"  {key}: {value}")

    print(f"\nResults written to: {output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
