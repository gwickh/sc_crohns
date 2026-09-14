#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy import sparse

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

# Plot colours
BLUE = "#2166AC"
RED = "#D73027"
GREY = "#BDBDBD"


def safe_filename(value: str) -> str:
    """Convert a cell-type name into a safe filename."""
    cleaned = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value),
    ).strip("_")

    return cleaned or "cell_type"


def make_unique_filenames(values: list[str]) -> dict[str, str]:
    """Create unique filenames if two labels produce the same safe name."""
    output = {}
    used = set()

    for value in values:
        base = safe_filename(value)
        candidate = base
        suffix = 2

        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1

        output[value] = candidate
        used.add(candidate)

    return output


def get_raw_counts(
    adata: ad.AnnData,
    counts_source: str,
) -> tuple[object, pd.Index]:
    """
    Retrieve raw counts and corresponding gene identifiers.

    counts_source can be:
        "counts" or another adata.layers key
        "X"
        "raw"
    """
    if counts_source == "X":
        matrix = adata.X
        gene_ids = adata.var_names

    elif counts_source == "raw":
        if adata.raw is None:
            raise ValueError(
                "counts_source='raw' was selected, but adata.raw is absent.",
            )

        matrix = adata.raw.X
        gene_ids = adata.raw.var_names

    else:
        if counts_source not in adata.layers:
            raise KeyError(
                f"Count layer {counts_source!r} was not found. "
                f"Available layers: {list(adata.layers.keys())}",
            )

        matrix = adata.layers[counts_source]
        gene_ids = adata.var_names

    if matrix.shape[0] != adata.n_obs:
        raise ValueError("The selected count matrix does not align with adata.obs.")

    if matrix.shape[1] != len(gene_ids):
        raise ValueError("The selected count matrix does not align with its gene IDs.")

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()

    if values.size and np.any(~np.isfinite(values)):
        raise ValueError("The selected count matrix contains non-finite values.")

    if values.size and np.any(values < 0):
        raise ValueError("The selected count matrix contains negative values.")

    if values.size and np.any(np.abs(values - np.rint(values)) > 1e-8):
        raise ValueError(
            "The selected count matrix contains non-integer values. "
            "PyDESeq2 requires raw counts, not normalized expression.",
        )

    return matrix, gene_ids.copy()


def normalise_feature_names(feature_ids: pd.Index) -> pd.Index:
    """
    Harmonise platform-specific feature names.

    Examples:
        A1BG_hg38 -> A1BG
        ENSG00000121410.12 -> ENSG00000121410
        ENSG00000121410.12_hg38 -> ENSG00000121410

    Only a terminal ``_hg38`` tag and a terminal dot followed entirely by
    digits are removed. Internal dots and text are retained.

    """
    names = pd.Index(feature_ids.astype(str).str.strip())

    # Remove the Parse genome tag first so a version immediately before it
    # becomes terminal and can then be removed safely.
    names = names.str.replace(
        r"(?i)_hg38$",
        "",
        regex=True,
    )

    # Remove Ensembl-style terminal versions such as .1 or .12, if present.
    names = names.str.replace(
        r"\.\d+$",
        "",
        regex=True,
    )

    # Also handle the uncommon ordering "GENE_hg38.1": removing the version
    # above exposes the genome suffix, which can now be removed.
    names = names.str.replace(
        r"(?i)_hg38$",
        "",
        regex=True,
    )

    invalid_names = names.isna() | (names == "")

    if invalid_names.any():
        bad = feature_ids[invalid_names].astype(str).tolist()[:10]
        raise ValueError(
            "Feature-name harmonisation produced missing or empty names. "
            f"Examples of original values: {bad}",
        )

    return pd.Index(names, name="gene_name")


def harmonise_count_features(
    matrix: object,
    feature_ids: pd.Index,
) -> tuple[object, pd.Index, pd.DataFrame]:
    """
    Harmonise feature names and sum count columns that map to the same gene.

    Summation is essential: simply dropping duplicated columns would discard
    counts from whichever platform used the other naming convention.
    """
    original_ids = pd.Index(feature_ids.astype(str), name="original_feature_id")
    harmonised_ids = normalise_feature_names(original_ids)

    if matrix.shape[1] != len(original_ids):
        raise ValueError(
            "The count matrix does not align with the supplied feature names.",
        )

    # factorize returns one integer group for every original feature column.
    feature_codes, unique_features = pd.factorize(
        harmonised_ids,
        sort=False,
    )

    unique_features = pd.Index(
        unique_features.astype(str),
        name="gene_name",
    )

    n_original = len(original_ids)
    n_harmonised = len(unique_features)

    # This sparse matrix maps original feature columns to harmonised genes.
    # Multiplication sums all columns belonging to the same gene.
    column_aggregator = sparse.csr_matrix(
        (
            np.ones(n_original, dtype=np.int8),
            (
                np.arange(n_original, dtype=np.int64),
                feature_codes,
            ),
        ),
        shape=(n_original, n_harmonised),
    )

    if sparse.issparse(matrix):
        harmonised_matrix = (matrix.tocsr() @ column_aggregator).tocsr()
        harmonised_matrix.eliminate_zeros()
    else:
        dense_matrix = np.asarray(matrix)
        harmonised_matrix = np.zeros(
            (dense_matrix.shape[0], n_harmonised),
            dtype=dense_matrix.dtype,
        )

        # Dense raw matrices of this size are unusual, but this avoids making
        # a dense feature-by-feature aggregation matrix.
        for original_position, harmonised_position in enumerate(feature_codes):
            harmonised_matrix[:, harmonised_position] += dense_matrix[
                :,
                original_position,
            ]

    audit = pd.DataFrame(
        {
            "original_feature_id": original_ids.astype(str),
            "harmonised_feature_id": harmonised_ids.astype(str),
        },
    )

    group_sizes = audit.groupby(
        "harmonised_feature_id",
        sort=False,
    )["original_feature_id"].transform("size")

    audit["n_original_features_for_gene"] = group_sizes.to_numpy()
    audit["was_harmonised"] = (
        audit["original_feature_id"] != audit["harmonised_feature_id"]
    )
    audit["was_collapsed"] = group_sizes.gt(1).to_numpy()

    n_hg38 = original_ids.str.contains(r"(?i)_hg38$", regex=True).sum()
    n_versioned = original_ids.str.contains(r"\.\d+(?:_hg38)?$", regex=True).sum()
    n_collapsed_columns = n_original - n_harmonised
    n_duplicate_genes = int(
        audit.loc[audit["was_collapsed"], "harmonised_feature_id"].nunique(),
    )

    print("\nFeature harmonisation:")
    print(f"  Original feature columns: {n_original:,}")
    print(f"  Features ending in _hg38: {n_hg38:,}")
    print(f"  Features with terminal numeric versions: {n_versioned:,}")
    print(f"  Unique harmonised genes: {n_harmonised:,}")
    print(f"  Duplicate columns summed: {n_collapsed_columns:,}")
    print(f"  Genes represented by multiple columns: {n_duplicate_genes:,}")

    if not unique_features.is_unique:
        raise RuntimeError(
            "Feature harmonisation failed to produce unique gene names.",
        )

    values = (
        harmonised_matrix.data
        if sparse.issparse(harmonised_matrix)
        else np.asarray(harmonised_matrix).ravel()
    )

    if values.size and np.any(np.abs(values - np.rint(values)) > 1e-8):
        raise ValueError(
            "Harmonised counts are not integer-valued.",
        )

    return harmonised_matrix, unique_features, audit


def write_harmonised_anndata(
    source_adata: ad.AnnData,
    harmonised_counts: object,
    harmonised_gene_names: pd.Index,
    feature_audit: pd.DataFrame,
    output_path: Path,
    original_counts_source: str,
) -> None:
    """
    Write a new AnnData containing the harmonised raw-count matrix.

    Cell metadata and observation-level embeddings are retained. Existing
    feature-level annotations, layers and varm entries are not copied because
    they no longer align after feature columns have been collapsed.
    """
    feature_summary = (
        feature_audit.groupby(
            "harmonised_feature_id",
            sort=False,
        )
        .agg(
            n_original_features=("original_feature_id", "size"),
            original_feature_ids=(
                "original_feature_id",
                lambda values: "|".join(values.astype(str)),
            ),
        )
        .reindex(harmonised_gene_names.astype(str))
    )

    feature_summary.index = pd.Index(
        harmonised_gene_names.astype(str),
        name="gene_name",
    )

    # X contains the harmonised raw integer counts. Do not also copy it into a
    # counts layer or .raw, because that would approximately double file size.
    harmonised_adata = ad.AnnData(
        X=harmonised_counts,
        obs=source_adata.obs.copy(),
        var=feature_summary,
    )

    # These arrays are aligned only to cells, so their dimensions remain valid.
    for key, values in source_adata.obsm.items():
        harmonised_adata.obsm[key] = values.copy()

    harmonised_adata.uns["feature_harmonisation"] = {
        "original_counts_source": str(original_counts_source),
        "count_storage": "X",
        "removed_terminal_suffix": "_hg38",
        "removed_terminal_version_pattern": r"\.\d+$",
        "duplicate_handling": "sum raw-count columns",
        "n_original_features": len(feature_audit),
        "n_harmonised_features": len(harmonised_gene_names),
        "n_collapsed_columns": int(len(feature_audit) - len(harmonised_gene_names)),
    }

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(f"Writing harmonised AnnData: {output_path}")

    harmonised_adata.write_h5ad(
        output_path,
        compression="gzip",
    )

    print(
        "Wrote harmonised AnnData with "
        f"{harmonised_adata.n_obs:,} cells and "
        f"{harmonised_adata.n_vars:,} unique genes.",
    )


def choose_contrast(
    diagnosis_levels: list[str],
    case_group: str | None,
    reference_group: str | None,
) -> tuple[str, str]:
    """Validate or automatically identify the diagnosis contrast."""
    if (case_group is None) != (reference_group is None):
        raise ValueError(
            "Set both case_group and reference_group, or leave both as None.",
        )

    if case_group is not None and reference_group is not None:
        if case_group not in diagnosis_levels:
            raise ValueError(
                f"Case {case_group!r} is absent. "
                f"Available diagnoses: {diagnosis_levels}",
            )

        if reference_group not in diagnosis_levels:
            raise ValueError(
                f"Reference {reference_group!r} is absent. "
                f"Available diagnoses: {diagnosis_levels}",
            )

        if case_group == reference_group:
            raise ValueError("case_group and reference_group must differ.")

        return case_group, reference_group

    if len(diagnosis_levels) != 2:
        raise ValueError(
            "Automatic contrast selection requires exactly two "
            f"diagnoses. Found: {diagnosis_levels}",
        )

    reference_pattern = re.compile(
        r"normal|healthy|control|^hc$",
        re.IGNORECASE,
    )

    reference_matches = [
        diagnosis
        for diagnosis in diagnosis_levels
        if reference_pattern.search(diagnosis)
    ]

    if len(reference_matches) != 1:
        raise ValueError(
            "Could not identify the reference diagnosis automatically. "
            f"Available diagnoses: {diagnosis_levels}. "
            "Set case_group and reference_group in main().",
        )

    reference_group = reference_matches[0]
    case_group = next(
        diagnosis for diagnosis in diagnosis_levels if diagnosis != reference_group
    )

    return case_group, reference_group


def aggregate_counts_by_sample(
    cell_counts: object,
    sample_ids: pd.Series,
    gene_ids: pd.Index,
) -> pd.DataFrame:
    """Sum cell-by-gene counts into sample-by-gene counts."""
    samples = sorted(sample_ids.astype(str).unique().tolist())

    sample_to_position = {sample: position for position, sample in enumerate(samples)}

    sample_codes = (
        sample_ids.astype(str).map(sample_to_position).to_numpy(dtype=np.int64)
    )

    aggregation_matrix = sparse.csr_matrix(
        (
            np.ones(len(sample_codes), dtype=np.int64),
            (
                sample_codes,
                np.arange(
                    len(sample_codes),
                    dtype=np.int64,
                ),
            ),
        ),
        shape=(len(samples), len(sample_codes)),
    )

    pseudobulk = aggregation_matrix @ cell_counts

    if sparse.issparse(pseudobulk):
        pseudobulk = pseudobulk.toarray()

    pseudobulk = np.asarray(pseudobulk)

    return pd.DataFrame(
        np.rint(pseudobulk).astype(np.int64),
        index=pd.Index(samples, name="sample_id"),
        columns=gene_ids,
    )


def classify_genes(
    results: pd.DataFrame,
    q_threshold: float,
    log2fc_threshold: float,
) -> pd.Series:
    """Classify genes according to the requested thresholds."""
    significant = results["q_value"].lt(q_threshold)

    downregulated = significant & results["log2FC"].lt(-log2fc_threshold)

    upregulated = significant & results["log2FC"].gt(log2fc_threshold)

    return pd.Series(
        np.select(
            [downregulated, upregulated],
            ["Down", "Up"],
            default="Not significant",
        ),
        index=results.index,
        dtype="string",
    )


def calculate_minus_log10_q(
    q_values: pd.Series,
) -> np.ndarray:
    """
    Convert q-values to -log10(q).

    Genes with missing q-values from independent filtering are treated as
    non-significant and plotted at zero.
    """
    values = q_values.fillna(1.0).clip(
        lower=np.finfo(float).tiny,
        upper=1.0,
    )

    return -np.log10(values.to_numpy(dtype=float))


def save_cell_type_volcano(
    results: pd.DataFrame,
    cell_type: str,
    case_group: str,
    reference_group: str,
    q_threshold: float,
    log2fc_threshold: float,
    output_prefix: Path,
) -> None:
    """Save one volcano plot for one cell type."""
    y_values = calculate_minus_log10_q(results["q_value"])

    fig, ax = plt.subplots(figsize=(7, 6))

    plot_groups = [
        ("Not significant", GREY),
        ("Down", BLUE),
        ("Up", RED),
    ]

    for label, colour in plot_groups:
        mask = results["significance"].eq(label).to_numpy()

        ax.scatter(
            results.loc[mask, "log2FC"],
            y_values[mask],
            s=10,
            alpha=0.65,
            color=colour,
            edgecolors="none",
            label=f"{label} (n={mask.sum():,})",
            rasterized=True,
        )

    ax.axvline(
        -log2fc_threshold,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.axvline(
        log2fc_threshold,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.axhline(
        -np.log10(q_threshold),
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.set_title(
        cell_type,
        fontweight="bold",
    )

    ax.set_xlabel(f"log2FC ({case_group} / {reference_group})")

    ax.set_ylabel(r"$-\log_{10}(q)$")

    ax.legend(
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()

    fig.savefig(
        output_prefix.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    fig.savefig(
        output_prefix.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_combined_volcano(
    results: pd.DataFrame,
    case_group: str,
    reference_group: str,
    q_threshold: float,
    log2fc_threshold: float,
    output_prefix: Path,
) -> None:
    """
    Save the combined volcano plot.

    Significant genes are coloured by cell type. All other genes are grey.
    """
    y_values = calculate_minus_log10_q(results["q_value"])

    fig, ax = plt.subplots(figsize=(11, 7))

    nonsignificant = results["significance"].eq("Not significant").to_numpy()

    ax.scatter(
        results.loc[nonsignificant, "log2FC"],
        y_values[nonsignificant],
        s=6,
        alpha=0.35,
        color=GREY,
        edgecolors="none",
        label="Not significant",
        rasterized=True,
    )

    cell_types = sorted(results["cell_type"].unique().tolist())

    if len(cell_types) <= 20:
        colour_map = plt.colormaps["tab20"].resampled(max(len(cell_types), 1))
    else:
        colour_map = plt.colormaps["turbo"].resampled(len(cell_types))

    for index, cell_type in enumerate(cell_types):
        significant = (
            results["cell_type"].eq(cell_type)
            & results["significance"].ne("Not significant")
        ).to_numpy()

        if not significant.any():
            continue

        ax.scatter(
            results.loc[significant, "log2FC"],
            y_values[significant],
            s=9,
            alpha=0.75,
            color=[colour_map(index)],
            edgecolors="none",
            label=(f"{cell_type} (n={significant.sum():,})"),
            rasterized=True,
        )

    ax.axvline(
        -log2fc_threshold,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.axvline(
        log2fc_threshold,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.axhline(
        -np.log10(q_threshold),
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.set_title(
        "Differential expression across cell types",
        fontweight="bold",
    )

    ax.set_xlabel(f"log2FC ({case_group} / {reference_group})")

    ax.set_ylabel(r"$-\log_{10}(q)$")

    ax.legend(
        title="Significant cell type",
        frameon=False,
        fontsize=7,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    fig.tight_layout()

    fig.savefig(
        output_prefix.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    fig.savefig(
        output_prefix.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def main() -> int:
    adata_path = (
        Path.home() / "project-area/data/crohns_scrnaseq/"
        "10c_14n_analysis/scvi_tools_output/"
        "adata_low_conf_hc_celltype_anchor_profiles.h5ad"
    )

    output_dir = (
        Path.home() / "project-area/data/crohns_scrnaseq/"
        "10c_14n_analysis/"
        "pseudobulk_dea"
    )

    # A new count-focused AnnData will be written here immediately after
    # feature names are harmonised and duplicate columns are summed.
    harmonised_h5ad_name = "adata_final.h5ad"

    # AnnData metadata columns.
    cell_type_column = "final_curated_labels"

    sample_column = "sample_id"

    diagnosis_column = "diagnosis"

    # Raw-count source:
    counts_source = "X"

    # Set these to exact values from adata.obs["diagnosis"].
    case_group = "Crohn's disease"
    reference_group = "Normal"

    # Filtering and significance settings.
    min_cells_per_sample = 10

    min_replicates_per_group = 2

    min_total_gene_count = 10

    q_threshold = 0.05

    log2fc_threshold = 1.0

    # PyDESeq2 settings.
    n_cpus = 4

    quiet = False

    # ================================================================

    adata_path = adata_path.expanduser()

    output_dir = output_dir.expanduser()

    harmonised_h5ad_path = output_dir / harmonised_h5ad_name

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if min_cells_per_sample < 1:
        raise ValueError("min_cells_per_sample must be at least 1.")

    if min_replicates_per_group < 2:
        raise ValueError("min_replicates_per_group must be at least 2.")

    if not 0 < q_threshold < 1:
        raise ValueError("q_threshold must lie between 0 and 1.")

    if log2fc_threshold < 0:
        raise ValueError("log2fc_threshold must be non-negative.")

    # Load AnnData
    adata = ad.read_h5ad(adata_path)

    keep_cells = (
        adata.obs["final_curated_labels"]
        .astype("string")
        .str.strip()
        .ne("low_conf")
        .fillna(False)
    )

    adata = adata[keep_cells].copy()

    # Check that sample_id exists
    if "sample_id" not in adata.obs.columns:
        raise KeyError(
            f"'sample_id' is absent. Available columns: {adata.obs.columns.tolist()}",
        )

    sample_ids = adata.obs["sample_id"].astype("string")

    # List unique sample IDs
    unique_sample_ids = sorted(sample_ids.dropna().unique().tolist())

    print("Unique sample_id values:")
    for sample_id in unique_sample_ids:
        print(sample_id)

    # Assign platform
    adata.obs["platform"] = np.where(
        sample_ids.str.contains(".h5ad", regex=False, na=False),
        "Parse",
        "10X",
    )

    adata.obs["platform"] = pd.Categorical(
        adata.obs["platform"],
        categories=["10X", "Parse"],
    )

    # Verify the assignment
    print("\nPlatform counts:")
    print(adata.obs["platform"].value_counts(dropna=False))

    print("\nSamples and assigned platforms:")
    print(
        adata.obs[["sample_id", "platform"]]
        .drop_duplicates()
        .sort_values(["platform", "sample_id"])
        .to_string(index=False),
    )

    # Obtain raw counts.
    counts, gene_ids = get_raw_counts(
        adata,
        counts_source,
    )

    # Harmonise feature names before any cell filtering or pseudobulk
    # aggregation. Terminal numeric versions and the Parse-specific _hg38
    # suffix are removed, and any resulting duplicate count columns are
    # summed so that both platforms contribute to one gene.
    counts, gene_ids, feature_harmonisation = harmonise_count_features(
        matrix=counts,
        feature_ids=gene_ids,
    )

    feature_harmonisation.to_csv(
        output_dir / "feature_harmonisation_audit.csv",
        index=False,
    )

    duplicate_feature_groups = feature_harmonisation.loc[
        feature_harmonisation["was_collapsed"]
    ].copy()

    duplicate_feature_groups.to_csv(
        output_dir / "collapsed_feature_groups.csv",
        index=False,
    )

    print(
        "Saved feature harmonisation audit to "
        f"{output_dir / 'feature_harmonisation_audit.csv'}",
    )

    write_harmonised_anndata(
        source_adata=adata,
        harmonised_counts=counts,
        harmonised_gene_names=gene_ids,
        feature_audit=feature_harmonisation,
        output_path=harmonised_h5ad_path,
        original_counts_source=counts_source,
    )

    platform_column = "platform"

    required_columns = [
        cell_type_column,
        sample_column,
        diagnosis_column,
        platform_column,
    ]

    missing_columns = [column for column in required_columns if column not in adata.obs]

    if missing_columns:
        raise KeyError(
            f"Missing adata.obs columns: {missing_columns}. "
            f"Available columns: {adata.obs.columns.tolist()}",
        )

    metadata = adata.obs[required_columns].copy()

    for column in required_columns:
        metadata[column] = metadata[column].astype("string").str.strip()

    # Remove cells with incomplete metadata.
    valid_metadata = metadata.notna().all(axis=1) & metadata.ne("").all(axis=1)

    if (~valid_metadata).any():
        print(f"Removing {(~valid_metadata).sum():,} cells with incomplete metadata.")

    valid_positions = np.flatnonzero(valid_metadata.to_numpy())

    metadata = metadata.loc[valid_metadata].copy()

    counts = counts[
        valid_positions,
        :,
    ]

    # Confirm that each donor has exactly one diagnosis.
    diagnoses_per_sample = metadata.groupby(
        sample_column,
        observed=True,
    )[diagnosis_column].nunique()

    invalid_samples = diagnoses_per_sample.index[diagnoses_per_sample.ne(1)].tolist()

    if invalid_samples:
        raise ValueError(
            "Each sample must have exactly one diagnosis. "
            f"Invalid samples: {invalid_samples}",
        )

    diagnosis_levels = sorted(metadata[diagnosis_column].unique().tolist())

    print(f"Diagnosis levels: {diagnosis_levels}")

    case_group, reference_group = choose_contrast(
        diagnosis_levels=diagnosis_levels,
        case_group=case_group,
        reference_group=reference_group,
    )

    print(f"Testing {case_group!r} versus {reference_group!r}.")

    print(f"Positive log2FC means higher expression in {case_group!r}.")

    # Keep only the selected diagnosis contrast.
    contrast_mask = metadata[diagnosis_column].isin(
        [
            reference_group,
            case_group,
        ],
    )

    contrast_positions = np.flatnonzero(contrast_mask.to_numpy())

    metadata = metadata.loc[contrast_mask].copy()

    counts = counts[
        contrast_positions,
        :,
    ]

    cell_types = sorted(metadata[cell_type_column].unique().tolist())

    safe_cell_type_names = make_unique_filenames(cell_types)

    inference = DefaultInference(n_cpus=n_cpus)

    all_results = []

    skipped_cell_types = []

    # ================================================================
    # RUN ONE PSEUDOBULK DEA PER CELL TYPE
    # ================================================================

    for cell_type in cell_types:
        print(f"\nProcessing cell type: {cell_type}")

        cell_type_mask = metadata[cell_type_column].eq(cell_type)

        cell_positions = np.flatnonzero(cell_type_mask.to_numpy())

        cell_metadata = metadata.loc[cell_type_mask].copy()

        cell_counts = counts[
            cell_positions,
            :,
        ]

        # Count the number of cells from each donor.
        cells_per_sample = cell_metadata[sample_column].value_counts()

        retained_samples = cells_per_sample.index[
            cells_per_sample.ge(min_cells_per_sample)
        ]

        retain_cell_mask = cell_metadata[sample_column].isin(retained_samples)

        retained_positions = np.flatnonzero(retain_cell_mask.to_numpy())

        cell_metadata = cell_metadata.loc[retain_cell_mask].copy()

        cell_counts = cell_counts[
            retained_positions,
            :,
        ]

        def skip(reason: str) -> None:
            print(f"Skipping {cell_type}: {reason}")

            skipped_cell_types.append(
                {
                    "cell_type": cell_type,
                    "reason": reason,
                },
            )

        if cell_metadata.empty:
            skip("No sample met min_cells_per_sample.")
            continue

        sample_metadata = (
            cell_metadata[
                [
                    sample_column,
                    diagnosis_column,
                    platform_column,
                ]
            ]
            .drop_duplicates()
            .set_index(sample_column)
        )

        pseudobulk_counts = aggregate_counts_by_sample(
            cell_counts=cell_counts,
            sample_ids=cell_metadata[sample_column],
            gene_ids=gene_ids,
        )

        sample_metadata = sample_metadata.loc[pseudobulk_counts.index].copy()

        replicate_counts = sample_metadata[diagnosis_column].value_counts()

        case_replicates = int(
            replicate_counts.get(
                case_group,
                0,
            ),
        )

        reference_replicates = int(
            replicate_counts.get(
                reference_group,
                0,
            ),
        )

        print(
            "Retained biological samples: "
            f"{case_group}={case_replicates}, "
            f"{reference_group}={reference_replicates}",
        )

        if (
            min(
                case_replicates,
                reference_replicates,
            )
            < min_replicates_per_group
        ):
            skip(
                "Insufficient biological replicates: "
                f"{case_group}={case_replicates}, "
                f"{reference_group}={reference_replicates}",
            )
            continue

        # Remove genes with almost no pseudobulk counts.
        expressed_genes = pseudobulk_counts.sum(axis=0).ge(min_total_gene_count)

        pseudobulk_counts = pseudobulk_counts.loc[
            :,
            expressed_genes,
        ]

        if pseudobulk_counts.shape[1] == 0:
            skip("No genes passed the total-count filter.")
            continue

        # PyDESeq2 metadata.
        deseq_metadata = pd.DataFrame(
            {
                "platform": pd.Categorical(sample_metadata[platform_column]),
                "diagnosis": pd.Categorical(
                    sample_metadata[diagnosis_column],
                    categories=[
                        reference_group,
                        case_group,
                    ],
                ),
            },
            index=pseudobulk_counts.index,
        )

        try:
            dds = DeseqDataSet(
                counts=pseudobulk_counts,
                metadata=deseq_metadata,
                design="~platform + diagnosis",
                refit_cooks=True,
                inference=inference,
                quiet=quiet,
            )

            dds.deseq2()

            statistics = DeseqStats(
                dds,
                contrast=[
                    "diagnosis",
                    case_group,
                    reference_group,
                ],
                alpha=q_threshold,
                inference=inference,
                quiet=quiet,
            )

            statistics.summary()

        except Exception as error:
            skip(f"PyDESeq2 failed: {type(error).__name__}: {error}")
            continue

        results = statistics.results_df.copy()

        results.index.name = "feature_id"

        results = results.reset_index().rename(
            columns={
                "baseMean": "base_mean",
                "log2FoldChange": "log2FC",
                "lfcSE": "log2FC_SE",
                "stat": "wald_statistic",
                "pvalue": "p_value",
                "padj": "q_value",
            },
        )

        # Features have already been harmonised to unique gene names before
        # pseudobulk aggregation, so no Ensembl dictionary lookup is needed.
        results["gene_name"] = results["feature_id"].astype("string")

        results.insert(
            1,
            "cell_type",
            cell_type,
        )

        results.insert(
            2,
            "case",
            case_group,
        )

        results.insert(
            3,
            "reference",
            reference_group,
        )

        results.insert(
            4,
            "n_case_samples",
            case_replicates,
        )

        results.insert(
            5,
            "n_reference_samples",
            reference_replicates,
        )

        results["significance"] = classify_genes(
            results=results,
            q_threshold=q_threshold,
            log2fc_threshold=log2fc_threshold,
        )

        # Put key gene identifiers first.
        preferred_columns = [
            "feature_id",
            "gene_name",
            "cell_type",
            "case",
            "reference",
            "n_case_samples",
            "n_reference_samples",
            "base_mean",
            "log2FC",
            "log2FC_SE",
            "wald_statistic",
            "p_value",
            "q_value",
            "significance",
        ]

        remaining_columns = [
            column for column in results.columns if column not in preferred_columns
        ]

        results = results[preferred_columns + remaining_columns]

        results = results.sort_values(
            [
                "q_value",
                "p_value",
            ],
            na_position="last",
        ).reset_index(drop=True)

        safe_name = safe_cell_type_names[cell_type]

        results.to_csv(
            output_dir / f"{safe_name}_DEA.csv",
            index=False,
        )

        save_cell_type_volcano(
            results=results,
            cell_type=cell_type,
            case_group=case_group,
            reference_group=reference_group,
            q_threshold=q_threshold,
            log2fc_threshold=log2fc_threshold,
            output_prefix=(output_dir / f"{safe_name}_volcano"),
        )

        all_results.append(results)

    # Save information about analyses that could not be performed.
    if skipped_cell_types:
        pd.DataFrame(skipped_cell_types).to_csv(
            output_dir / "skipped_cell_types.csv",
            index=False,
        )

    if not all_results:
        raise RuntimeError(
            "No cell type produced DEA results. See skipped_cell_types.csv.",
        )

    # ================================================================
    # COMBINED RESULTS
    # ================================================================

    combined_results = pd.concat(
        all_results,
        ignore_index=True,
    )

    combined_results.to_csv(
        output_dir / "combined_DEA.csv",
        index=False,
    )

    save_combined_volcano(
        results=combined_results,
        case_group=case_group,
        reference_group=reference_group,
        q_threshold=q_threshold,
        log2fc_threshold=log2fc_threshold,
        output_prefix=(output_dir / "combined_volcano"),
    )

    settings = [
        f"input={adata_path}",
        f"counts_source={counts_source}",
        f"harmonised_h5ad={harmonised_h5ad_path}",
        "feature_harmonisation=strip terminal numeric version and terminal _hg38; sum duplicates",
        f"cell_type_column={cell_type_column}",
        f"sample_column={sample_column}",
        f"diagnosis_column={diagnosis_column}",
        f"case={case_group}",
        f"reference={reference_group}",
        f"min_cells_per_sample={min_cells_per_sample}",
        f"min_replicates_per_group={min_replicates_per_group}",
        f"min_total_gene_count={min_total_gene_count}",
        f"q_threshold={q_threshold}",
        f"log2fc_threshold={log2fc_threshold}",
        f"python={sys.version.replace(chr(10), ' ')}",
    ]

    (output_dir / "run_settings.txt").write_text(
        "\n".join(settings) + "\n",
        encoding="utf-8",
    )

    print(f"\nDEA complete. Results written to {output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
