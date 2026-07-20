#!/usr/bin/env python3
"""Pseudobulk T cells and perform DEA to identify if PIM3 is a DEG."""

import re
import warnings
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy import sparse
from sklearn.decomposition import PCA

mpl.use("Agg")

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True


# Paths
reference_path = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy/adata_umap.h5ad",
)

annotated_path = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/"
    "scvi_tools_output/query_concat_curated.h5ad",
)

output_dir = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/pseudobulk_dea/")

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

# Settings
category_key = "category"
t_cell_category = "T cells"

sample_key = "sample_id"
condition_key = "Diagnosis"

# PIM3 Ensembl gene ID
TARGET_ENSEMBL_ID = "ENSG00000198355"


# Helper functions
def strip_ensembl_version(value) -> str:
    """Strip version off ensembl IDs."""
    if pd.isna(value):
        return ""

    return re.sub(
        r"\.\d+$",
        "",
        str(value).strip(),
    )


def normalise_ensembl_ids(values) -> pd.Index:
    """Return version-free Ensembl identifiers."""
    return pd.Index(
        [strip_ensembl_version(value) for value in values],
        dtype="object",
    )


def check_integer_counts(
    x,
    matrix_name: str,
) -> None:
    """Check that an expression matrix contains non-negative integer counts."""
    values = x.data if sparse.issparse(x) else np.asarray(x).ravel()

    if values.size == 0:
        msg = "{matrix_name} is empty."
        raise ValueError(msg)

    if not np.all(np.isfinite(values)):
        msg = f"{matrix_name} contains non-finite values."
        raise ValueError(msg)

    if np.min(values) < 0:
        msg = f"{matrix_name} contains negative values."
        raise ValueError(msg)

    if not np.allclose(
        values,
        np.rint(values),
        atol=1e-6,
    ):
        msg = f"{matrix_name} does not appear to contain raw integer \
            counts. Do not use normalised or log-transformed values \
            for pseudobulk DESeq2."
        raise ValueError(msg)


def get_ensembl_ids(
    var: pd.DataFrame,
    var_names: pd.Index,
    target_id: str | None = None,
    source_name: str = "matrix",
) -> tuple[pd.Index, str]:
    """Find Ensembl IDs in var."""
    if "ensembl_id" not in var.columns:
        msg = f"""
        ensembl_id column is absent from var.
        Available columns: {list(var.columns)}
        """
        raise KeyError(msg)

    ids = normalise_ensembl_ids(var["ensembl_id"])

    ensembl_fraction = pd.Series(ids).astype(str).str.startswith("ENSG").mean()

    if ensembl_fraction > 0.5:
        "ensembl_id".append(
            (
                f"var[{column!r}]",
                ids,
            ),
        )

    var_name_ids = normalise_ensembl_ids(var_names)

    var_name_ensembl_fraction = (
        pd.Series(var_name_ids).astype(str).str.startswith("ENSG").mean()
    )

    if var_name_ensembl_fraction > 0.5:
        "ensembl_id".append(
            (
                "var_names",
                var_name_ids,
            ),
        )

    if not "ensembl_id":
        msg = (
            f"No Ensembl-ID source was found for {source_name}. \
            Available np.var columns: {list(np.var.columns)}",
        )
        raise KeyError(msg)

    if target_id is None:
        source_label, ids = "ensembl_id"[0]

        print(f"{source_name}: using {source_label} for Ensembl IDs.")

        return ids, source_label

    target_id = strip_ensembl_version(target_id)

    for source_label, ids in "ensembl_id":
        if target_id in set(ids):
            print(f"{source_name}: found {target_id} using {source_label}.")

            return ids, source_label

    searched_sources = [label for label, _ in "ensembl_id"]

    msg = (
        f"{target_id} was not found in {source_name}.\
        Ensembl-ID sources searched: {searched_sources}",
    )

    raise KeyError(msg)


def collapse_duplicate_ensembl_ids(
    arr,
    ensembl_ids: pd.Index,
) -> tuple:
    """Collapse duplicated Ensembl features by summing their counts."""
    ensembl_ids = normalise_ensembl_ids(ensembl_ids)

    valid = pd.Series(ensembl_ids).astype(str).str.startswith("ENSG").to_numpy()

    n_invalid = int((~valid).sum())

    if n_invalid > 0:
        print(f"Dropping {n_invalid:,} features without valid ENSG identifiers.")

        arr = arr[:, valid]
        ensembl_ids = ensembl_ids[valid]

    if len(ensembl_ids) == 0:
        msg = "No valid Ensembl gene IDs remain."
        raise ValueError(msg)

    unique_ids = pd.Index(
        pd.unique(ensembl_ids),
        dtype="object",
        name="feature_id",
    )

    n_duplicates = len(ensembl_ids) - len(unique_ids)

    if n_duplicates == 0:
        print("No duplicated Ensembl IDs were found.")

        return arr, unique_ids

    print(f"Collapsing {n_duplicates:,} duplicated Ensembl features by summing counts.")

    gene_codes = pd.Categorical(
        ensembl_ids,
        categories=unique_ids,
        ordered=True,
    ).codes

    feature_to_gene = sparse.csr_matrix(
        (
            np.ones(
                len(ensembl_ids),
                dtype=np.int64,
            ),
            (
                np.arange(len(ensembl_ids)),
                gene_codes,
            ),
        ),
        shape=(
            len(ensembl_ids),
            len(unique_ids),
        ),
    )

    if sparse.issparse(arr):
        collapsed_arr = arr.tocsr() @ feature_to_gene

        collapsed_arr = collapsed_arr.tocsr()

    else:
        collapsed_arr = sparse.csr_matrix(np.asarray(arr)) @ feature_to_gene

        collapsed_arr = collapsed_arr.tocsr()

    return collapsed_arr, unique_ids


# =====================================================================
# Load objects
# =====================================================================

print(f"Loading reference object:\n{reference_path}")

reference = ad.read_h5ad(reference_path)

print(f"Reference object: {reference.n_obs:,} cells * {reference.n_vars:,} variables")

print(f"\nLoading annotated object:\n{annotated_path}")

annotated = ad.read_h5ad(annotated_path)

print(f"Annotated object: {annotated.n_obs:,} cells * {annotated.n_vars:,} variables")


# =====================================================================
# Validate cell identifiers
# =====================================================================

if not reference.obs_names.is_unique:
    raise ValueError("reference.obs_names contains duplicated cell identifiers.")

if not annotated.obs_names.is_unique:
    raise ValueError("annotated.obs_names contains duplicated cell identifiers.")

if category_key not in annotated.obs.columns:
    raise KeyError(
        f"{category_key!r} is absent from annotated.obs.\n"
        f"Available columns: {list(annotated.obs.columns)}",
    )


# =====================================================================
# Select curated T cells
# =====================================================================
# Selecting category == "T cells" automatically excludes cells whose
# broad annotation is "Unknown".

category_values = annotated.obs[category_key].astype("string")

annotated_tcell_ids = annotated.obs_names[
    category_values.eq(t_cell_category).fillna(False)
]

print(
    f"\nCells with {category_key} == {t_cell_category!r}: {len(annotated_tcell_ids):,}",
)

if len(annotated_tcell_ids) == 0:
    raise ValueError(f"No cells had {category_key} == {t_cell_category!r}.")


# =====================================================================
# Match T cells to reference object
# =====================================================================

matched_tcell_ids = annotated_tcell_ids[annotated_tcell_ids.isin(reference.obs_names)]

unmatched_tcell_ids = annotated_tcell_ids[
    ~annotated_tcell_ids.isin(reference.obs_names)
]

matched_fraction = len(matched_tcell_ids) / len(annotated_tcell_ids)

print(f"Matched T cells: {len(matched_tcell_ids):,}")

print(f"Unmatched T cells: {len(unmatched_tcell_ids):,}")

print(f"Fraction recovered: {matched_fraction:.2%}")

if len(matched_tcell_ids) == 0:
    print("\nExample reference cell IDs:")

    print(reference.obs_names[:10].tolist())

    print("\nExample annotated cell IDs:")

    print(annotated.obs_names[:10].tolist())

    raise ValueError(
        "No T-cell identifiers matched between the reference and annotated objects.",
    )

if matched_fraction < 0.95:
    print("\nExample unmatched T-cell identifiers:")

    print(unmatched_tcell_ids[:20].tolist())

    raise ValueError("Fewer than 95% of annotated T cells were recovered.")


# =====================================================================
# Choose raw-count matrix using Ensembl IDs only
# =====================================================================

if "counts" in reference.layers:
    count_ids, id_source = get_ensembl_ids(
        var=reference.var,
        var_names=reference.var_names,
        target_id="ensembl_id",
        source_name='reference.layers["counts"]',
    )

    check_integer_counts(
        reference.layers["counts"],
        'reference.layers["counts"]',
    )

    print('Using reference.layers["counts"] as the raw-count matrix.')

count_source = {
    "arr": reference.layers["counts"],
    "ensembl_ids": count_ids,
    "source": 'reference.layers["counts"]',
    "id_source": id_source,
}

count_arr = count_source["arr"]
count_ensembl_ids = count_source["ensembl_ids"]

print(f"\nSelected count source: {count_source['source']}")

print(f"Selected Ensembl-ID source: {count_source['id_source']}")


# =====================================================================
# Extract T-cell count rows
# =====================================================================

reference_row_indices = reference.obs_names.get_indexer(matched_tcell_ids)

if np.any(reference_row_indices < 0):
    raise RuntimeError("Some matched cells could not be indexed in the reference.")

tcell_count_X = count_arr[reference_row_indices, :]


# =====================================================================
# Collapse duplicated Ensembl identifiers
# =====================================================================

tcell_count_X, unique_ensembl_ids = collapse_duplicate_ensembl_ids(
    X=tcell_count_X,
    ensembl_ids=count_ensembl_ids,
)

if TARGET_ENSEMBL_ID not in unique_ensembl_ids:
    raise KeyError(
        f"{TARGET_ENSEMBL_ID} was lost while standardising "
        "the Ensembl feature identifiers.",
    )


# =====================================================================
# Construct T-cell observation metadata
# =====================================================================

tcell_obs = reference.obs.loc[matched_tcell_ids].copy()


transfer_columns = list(annotated.obs.columns)

for column in transfer_columns:
    tcell_obs[column] = annotated.obs.loc[
        matched_tcell_ids,
        column,
    ]


# =====================================================================
# Construct Ensembl-only var DataFrame
# =====================================================================

tcell_var = pd.DataFrame(
    {
        "ensembl_id": (unique_ensembl_ids.astype(str)),
    },
    index=pd.Index(
        unique_ensembl_ids.astype(str),
        name="feature_id",
    ),
)

# The index name is deliberately "feature_id", not "ensembl_id to prevent the AnnData
# HDF5 writing error caused by an index and column sharing the same name while containing
# different values.

# Construct T-cell AnnData

t_cells = ad.AnnData(
    X=tcell_count_X,
    obs=tcell_obs,
    var=tcell_var,
)

if not t_cells.var_names.is_unique:
    raise ValueError("Ensembl feature identifiers are still duplicated.")

# Store the same raw integer counts explicitly for pseudobulk DEA.
t_cells.layers["counts"] = t_cells.X.copy()

t_cells.uns["count_source"] = count_source["source"]

t_cells.uns["ensembl_id_source"] = count_source["id_source"]

t_cells.uns["target_ensembl_id"] = TARGET_ENSEMBL_ID


# =====================================================================
# Validate annotations and DEA metadata
# =====================================================================

unexpected_categories = set(t_cells.obs[category_key].dropna().astype(str).unique()) - {
    t_cell_category,
}

if unexpected_categories:
    raise ValueError(
        "The T-cell object contains unexpected categories: "
        f"{sorted(unexpected_categories)}",
    )

for required_column in [
    sample_key,
    condition_key,
]:
    if required_column not in t_cells.obs.columns:
        raise KeyError(
            f"{required_column!r} is required for pseudobulk DEA "
            "but is absent from t_cells.obs.",
        )

print("\nBiological samples per diagnosis:")

sample_metadata = t_cells.obs[
    [
        sample_key,
        condition_key,
    ]
].drop_duplicates()

print(sample_metadata[condition_key].value_counts(dropna=False))

print("\nT-cell counts per biological sample:")

print(
    t_cells.obs.groupby(
        [
            condition_key,
            sample_key,
        ],
        observed=True,
    )
    .size()
    .sort_values(ascending=False),
)


# =====================================================================
# Target Ensembl-ID diagnostics
# =====================================================================

target_index = t_cells.var_names.get_loc(TARGET_ENSEMBL_ID)

target_counts = t_cells.layers["counts"][
    :,
    target_index,
]

if sparse.issparse(target_counts):
    total_target_counts = int(target_counts.sum())

    target_positive_cells = int(target_counts.getnnz())

else:
    target_counts = np.asarray(target_counts).ravel()

    total_target_counts = int(target_counts.sum())

    target_positive_cells = int((target_counts > 0).sum())

print(f"\nTarget Ensembl ID: {TARGET_ENSEMBL_ID}")

print(f"Total target counts in curated T cells: {total_target_counts:,}")

print(
    f"T cells with detectable target counts: "
    f"{target_positive_cells:,} / "
    f"{t_cells.n_obs:,}",
)


# =====================================================================
# Final write checks
# =====================================================================

print(f"\nvar index name: {t_cells.var.index.name!r}")

print(
    "var_names equal ensembl_id column: "
    f"{np.array_equal(t_cells.var_names.astype(str), t_cells.var['ensembl_id'].astype(str).to_numpy())}",
)

check_integer_counts(
    t_cells.layers["counts"],
    't_cells.layers["counts"]',
)


# =====================================================================
# AnnData settings
# =====================================================================

counts_layer = "counts"

sample_key = "sample_id"

condition_key = "Diagnosis"

category_key = "category"
t_cell_category = "T cells"


# =====================================================================
# Model settings
# =====================================================================

# Add sample-level covariates when appropriate, for example:
#
covariates = ["platform"]

alpha = 0.05
lfc_cutoff = 1.0

n_cpus = 8

size_factors_fit_type = "poscounts"


# =====================================================================
# Filtering settings
# =====================================================================

min_cells_per_sample = 20

min_total_gene_count = 10
min_count_in_sample = 5
min_samples_expressed = 3

# The target Ensembl feature will be retained for testing whenever it has at least one count,
# even if it fails the generic gene filter.
force_include_target_gene = True

top_n_heatmap = 30


# =====================================================================
# Helper functions
# =====================================================================


def normalise_condition(value) -> str | None:
    """Map diagnosis values onto Crohns and Normal."""
    if pd.isna(value):
        return None

    value = str(value).strip().lower()
    value = value.replace("’", "'")

    if "crohn" in value:
        return "Crohns"

    if value in "normal":
        return "Normal"

    return None


def check_integer_counts(X, matrix_name: str) -> None:
    """Confirm that the selected expression matrix contains raw counts."""
    if sparse.issparse(X):
        values = X.data
    else:
        values = np.asarray(X).ravel()

    if values.size == 0:
        raise ValueError(f"{matrix_name} is empty.")

    if values.size > 1_000_000:
        rng = np.random.default_rng(0)
        values = rng.choice(
            values,
            size=1_000_000,
            replace=False,
        )

    if not np.all(np.isfinite(values)):
        raise ValueError(f"{matrix_name} contains non-finite values.")

    if np.min(values) < 0:
        raise ValueError(f"{matrix_name} contains negative values.")

    if not np.allclose(
        values,
        np.rint(values),
        atol=1e-6,
    ):
        raise ValueError(
            f"{matrix_name} does not appear to contain raw integer "
            "counts. Do not use log-normalised expression for DESeq2.",
        )


def canonicalise_ensembl_id(value) -> str:
    """
    Remove an optional Ensembl version suffix.

    Example:
        ENSG00000198355.4 -> ENSG00000198355

    """
    if pd.isna(value):
        return ""

    return re.sub(
        r"\.\d+$",
        "",
        str(value).strip(),
    )


def find_target_feature_by_ensembl_id(
    adata: ad.AnnData,
    target_ensembl_id: str,
) -> str:
    """
    Find a matrix feature using Ensembl IDs in adata.var_names only.

    No gene-symbol or gene-name columns are searched.
    """
    target = canonicalise_ensembl_id(target_ensembl_id)

    canonical_var_names = pd.Index(
        [canonicalise_ensembl_id(feature) for feature in adata.var_names],
    )

    matches = np.flatnonzero(canonical_var_names == target)

    if len(matches) == 0:
        raise KeyError(
            f"{target} was not found in adata.var_names.\n"
            f"First five var_names: "
            f"{adata.var_names[:5].astype(str).tolist()}",
        )

    if len(matches) > 1:
        matched_features = [str(adata.var_names[position]) for position in matches]

        raise ValueError(
            f"{target} matched multiple matrix features: "
            f"{matched_features}. The object should contain one "
            "feature per canonical Ensembl gene ID.",
        )

    return str(adata.var_names[matches[0]])


def validate_sample_metadata(
    obs: pd.DataFrame,
    sample_column: str,
    metadata_columns: list[str],
) -> pd.DataFrame:
    """
    Construct one metadata row per sample and confirm that each sample
    has exactly one value for every model variable.
    """
    records = []

    for sample, sample_obs in obs.groupby(
        sample_column,
        observed=True,
        sort=False,
    ):
        record = {
            sample_column: str(sample),
        }

        for column in metadata_columns:
            values = sample_obs[column].dropna().astype(str).unique()

            if len(values) == 0:
                record[column] = np.nan

            elif len(values) == 1:
                record[column] = values[0]

            else:
                raise ValueError(
                    f"Sample {sample!r} has multiple values for "
                    f"{column!r}: {values.tolist()}",
                )

        records.append(record)

    return pd.DataFrame(records).set_index(sample_column)


def aggregate_pseudobulk_counts(
    adata: ad.AnnData,
    layer: str,
    sample_column: str,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Sum raw cell-level counts within each biological sample.
    """
    X = adata.layers[layer]

    if sparse.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)

    sample_values = adata.obs[sample_column].astype(str)
    sample_categories = pd.Categorical(sample_values)

    if np.any(sample_categories.codes < 0):
        raise ValueError(f"{sample_column!r} contains missing sample identifiers.")

    n_samples = len(sample_categories.categories)
    n_cells = adata.n_obs

    # Sample-by-cell indicator matrix.
    aggregation_matrix = sparse.csr_matrix(
        (
            np.ones(n_cells, dtype=np.int64),
            (
                sample_categories.codes,
                np.arange(n_cells),
            ),
        ),
        shape=(n_samples, n_cells),
    )

    pseudobulk = aggregation_matrix @ X

    if sparse.issparse(pseudobulk):
        pseudobulk = pseudobulk.toarray()

    pseudobulk = np.asarray(pseudobulk)

    if not np.allclose(
        pseudobulk,
        np.rint(pseudobulk),
        atol=1e-6,
    ):
        raise ValueError("Pseudobulk values are not integers. Check the counts layer.")

    counts = pd.DataFrame(
        np.rint(pseudobulk).astype(np.int64),
        index=sample_categories.categories.astype(str),
        columns=adata.var_names.astype(str),
    )

    n_cells_per_sample = sample_values.value_counts().reindex(counts.index).astype(int)

    return counts, n_cells_per_sample


def check_design_rank(
    metadata: pd.DataFrame,
    variables: list[str],
) -> None:
    """
    Detect obvious rank deficiency/confounding in a main-effects design.
    """
    design_check = pd.get_dummies(
        metadata[variables],
        drop_first=True,
        dtype=float,
    )

    design_matrix = np.column_stack(
        [
            np.ones(metadata.shape[0]),
            design_check.to_numpy(dtype=float),
        ],
    )

    rank = np.linalg.matrix_rank(design_matrix)

    if rank < design_matrix.shape[1]:
        raise ValueError(
            "The model design is rank deficient. Diagnosis may be "
            "completely confounded with a covariate such as platform, "
            "or two covariates may encode the same information.",
        )


def save_figure(
    fig,
    output_stem: Path,
) -> None:
    """Save a Matplotlib figure as PDF and PNG."""
    fig.savefig(
        output_stem.with_suffix(".pdf"),
        bbox_inches="tight",
        pad_inches=0.2,
    )

    fig.savefig(
        output_stem.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
    )

    plt.close(fig)


# =====================================================================
# Load and validate AnnData
# =====================================================================


t_cells.var_names_make_unique()

sum_pim3 = t_cells.layers["counts"][:, t_cells.var_names.get_loc(TARGET_ENSEMBL_ID)]
sum_pim3 = (
    sum_pim3.toarray().ravel()
    if sparse.issparse(sum_pim3)
    else np.asarray(sum_pim3).ravel()
)

sample_counts = (
    pd.Series(sum_pim3, index=t_cells.obs_names)
    .groupby(t_cells.obs[sample_key].astype(str))
    .sum()
)

keep_samples = sample_counts[sample_counts > 1].index

t_cells = t_cells[t_cells.obs[sample_key].astype(str).isin(keep_samples)].copy()

print(f"Input object: {t_cells.n_obs:,} cells × {t_cells.n_vars:,} genes")

required_obs_columns = [
    sample_key,
    condition_key,
] + covariates

missing_obs_columns = [
    column for column in required_obs_columns if column not in t_cells.obs.columns
]

if missing_obs_columns:
    raise KeyError(
        f"Missing t_cells.obs columns: {missing_obs_columns}\n"
        f"Available columns: {list(t_cells.obs.columns)}",
    )

if counts_layer not in t_cells.layers:
    raise KeyError(
        f"Raw-count layer {counts_layer!r} was not found.\n"
        f"Available layers: {list(t_cells.layers.keys())}",
    )

check_integer_counts(
    t_cells.layers[counts_layer],
    f't_cells.layers["{counts_layer}"]',
)


# =====================================================================
# Safeguard: retain only curated T cells
# =====================================================================

if category_key in t_cells.obs.columns:
    category_mask = (
        t_cells.obs[category_key].astype("string").eq(t_cell_category).fillna(False)
    )

    unexpected = int((~category_mask).sum())

    if unexpected > 0:
        print(
            f"Removing {unexpected:,} cells whose "
            f"{category_key!r} is not {t_cell_category!r}.",
        )

        t_cells = t_cells[category_mask].copy()

if t_cells.n_obs == 0:
    raise ValueError("No curated T cells remain.")


# =====================================================================
# Locate target Ensembl ID
# =====================================================================

target_feature = find_target_feature_by_ensembl_id(
    adata=t_cells,
    target_ensembl_id=TARGET_ENSEMBL_ID,
)

target_features = [target_feature]

print(f"Target Ensembl ID: {TARGET_ENSEMBL_ID}")
print(f"Matched matrix feature: {target_feature}")


# =====================================================================
# Create sample-level pseudobulk counts
# =====================================================================

pseudobulk_counts, n_cells_per_sample = aggregate_pseudobulk_counts(
    adata=t_cells,
    layer=counts_layer,
    sample_column=sample_key,
)

metadata_columns = [
    condition_key,
] + covariates

sample_metadata = validate_sample_metadata(
    obs=t_cells.obs,
    sample_column=sample_key,
    metadata_columns=metadata_columns,
)

sample_metadata["n_cells"] = n_cells_per_sample.reindex(sample_metadata.index)

sample_metadata["condition"] = sample_metadata[condition_key].map(normalise_condition)

unmapped_samples = sample_metadata[sample_metadata["condition"].isna()]

if not unmapped_samples.empty:
    print("\nDropping samples with unrecognised diagnosis labels:")

    print(unmapped_samples[[condition_key]])

sample_metadata = sample_metadata[
    sample_metadata["condition"].isin(["Normal", "Crohns"])
].copy()

sample_metadata = sample_metadata[
    sample_metadata["n_cells"] >= min_cells_per_sample
].copy()

pseudobulk_counts = pseudobulk_counts.loc[sample_metadata.index].copy()

if pseudobulk_counts.empty:
    raise ValueError("No samples passed the condition and minimum-cell filters.")

sample_metadata["condition"] = pd.Categorical(
    sample_metadata["condition"],
    categories=["Normal", "Crohns"],
    ordered=True,
)

for covariate in covariates:
    if not pd.api.types.is_numeric_dtype(sample_metadata[covariate]):
        sample_metadata[covariate] = sample_metadata[covariate].astype("category")


# =====================================================================
# Sample checks
# =====================================================================

condition_sample_counts = (
    sample_metadata["condition"]
    .value_counts()
    .reindex(["Normal", "Crohns"])
    .fillna(0)
    .astype(int)
)

print("\nSamples retained per condition:")
print(condition_sample_counts)

if condition_sample_counts["Normal"] < 2 or condition_sample_counts["Crohns"] < 2:
    raise ValueError(
        "At least two independent biological samples are required "
        "in both Normal and Crohn's groups.",
    )

if condition_sample_counts["Normal"] < 3 or condition_sample_counts["Crohns"] < 3:
    warnings.warn(
        "One condition has fewer than three biological samples. "
        "Interpret inferential statistics cautiously.",
    )

if "platform" in sample_metadata.columns:
    platform_table = pd.crosstab(
        sample_metadata["platform"],
        sample_metadata["condition"],
        dropna=False,
    )

    platform_table.to_csv(output_dir / "platform_by_condition.csv")

    print("\nPlatform * condition:")
    print(platform_table)

    if "platform" not in covariates:
        warnings.warn(
            "'platform' is present but is not included in the model. "
            "Consider adding it to covariates if it is not confounded "
            "with condition.",
        )


# =====================================================================
# Save raw pseudobulk inputs
# =====================================================================

sample_metadata.to_csv(output_dir / "sample_metadata.csv")

pseudobulk_counts.to_csv(output_dir / "raw_pseudobulk_counts_all_genes.csv")


# =====================================================================
# Gene filtering
# =====================================================================

gene_filter_qc = pd.DataFrame(index=pseudobulk_counts.columns)

gene_filter_qc.index.name = "feature_id"


gene_filter_qc["total_count"] = pseudobulk_counts.sum(axis=0)

gene_filter_qc["n_samples_count_ge_threshold"] = pseudobulk_counts.ge(
    min_count_in_sample,
).sum(axis=0)

gene_filter_qc["passed_default_filter"] = (
    gene_filter_qc["total_count"] >= min_total_gene_count
) & (gene_filter_qc["n_samples_count_ge_threshold"] >= min_samples_expressed)

gene_filter_qc["forced_target_gene"] = False

keep_genes = gene_filter_qc["passed_default_filter"].copy()

if force_include_target_gene:
    for feature in target_features:
        if (
            feature in gene_filter_qc.index
            and gene_filter_qc.loc[feature, "total_count"] > 0
        ):
            if not keep_genes.loc[feature]:
                print(
                    f"Forcing {TARGET_ENSEMBL_ID} feature {feature!r} "
                    "through the generic expression filter.",
                )

                gene_filter_qc.loc[
                    feature,
                    "forced_target_gene",
                ] = True

            keep_genes.loc[feature] = True

gene_filter_qc["retained_for_DEA"] = keep_genes

gene_filter_qc.to_csv(output_dir / "gene_filter_QC.csv")

counts_filtered = pseudobulk_counts.loc[
    :,
    keep_genes,
].copy()

print(
    f"\nGenes retained for DEA: "
    f"{counts_filtered.shape[1]:,} / "
    f"{pseudobulk_counts.shape[1]:,}",
)

if counts_filtered.shape[1] == 0:
    raise ValueError("No genes passed the expression filter.")

target_tested_features = [
    feature for feature in target_features if feature in counts_filtered.columns
]

if not target_tested_features:
    warnings.warn(
        f"{TARGET_ENSEMBL_ID} has no detectable pseudobulk counts and cannot be tested.",
    )


# =====================================================================
# Model design
# =====================================================================

design_variables = covariates + ["condition"]

model_metadata = sample_metadata[design_variables].copy()

missing_model_values = model_metadata.isna().any(axis=1)

if missing_model_values.any():
    dropped_samples = model_metadata.index[missing_model_values].tolist()

    print("\nDropping samples with missing model covariates:")
    print(dropped_samples)

    model_metadata = model_metadata.loc[~missing_model_values].copy()

    sample_metadata = sample_metadata.loc[model_metadata.index].copy()

    counts_filtered = counts_filtered.loc[model_metadata.index].copy()

    pseudobulk_counts = pseudobulk_counts.loc[model_metadata.index].copy()

check_design_rank(
    metadata=model_metadata,
    variables=design_variables,
)

design_formula = "~ " + " + ".join(design_variables)

print(f"\nDESeq2 design: {design_formula}")
print("Contrast: Crohns versus Normal (positive log2FC = higher in Crohn's)")


# =====================================================================
# Fit PyDESeq2
# =====================================================================

inference = DefaultInference(n_cpus=n_cpus)

dds = DeseqDataSet(
    counts=counts_filtered,
    metadata=model_metadata,
    design=design_formula,
    refit_cooks=True,
    size_factors_fit_type=size_factors_fit_type,
    inference=inference,
)

dds.deseq2()

stat_res = DeseqStats(
    dds,
    contrast=[
        "condition",
        "Crohns",
        "Normal",
    ],
    alpha=alpha,
    inference=inference,
)

stat_res.summary()


# =====================================================================
# Extract DEA results
# =====================================================================

results = stat_res.results_df.copy()

results.index = results.index.astype(str)
results.index.name = "feature_id"

results = results.reset_index()


results["significant_FDR"] = results["padj"].notna() & (results["padj"] < alpha)

results["significant_FDR_and_LFC"] = results["significant_FDR"] & (
    results["log2FoldChange"].abs() >= lfc_cutoff
)

results["direction"] = np.select(
    [
        (results["significant_FDR_and_LFC"] & (results["log2FoldChange"] > 0)),
        (results["significant_FDR_and_LFC"] & (results["log2FoldChange"] < 0)),
    ],
    [
        "Upregulated",
        "Downregulated",
    ],
    default="Not significant",
)

results = results.sort_values(
    [
        "padj",
        "pvalue",
    ],
    na_position="last",
)

results.to_csv(
    output_dir / "Crohns_vs_Normal_DESeq2_all_genes.csv",
    index=False,
)

significant_results = results[results["significant_FDR_and_LFC"]].copy()

significant_results.to_csv(
    output_dir / "Crohns_vs_Normal_DESeq2_significant.csv",
    index=False,
)


# =====================================================================
# Target Ensembl-ID statistics
# =====================================================================

target_results = results[results["feature_id"].eq(target_feature)].copy()

target_results.to_csv(
    output_dir / f"{TARGET_ENSEMBL_ID}_DESeq2_statistics.csv",
    index=False,
)

print(f"\n{TARGET_ENSEMBL_ID} DEA statistics:")

if target_results.empty:
    print(f"{TARGET_ENSEMBL_ID} was not tested. Inspect gene_filter_QC.csv.")

else:
    columns_to_print = [
        "feature_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "direction",
    ]

    print(target_results[columns_to_print].to_string(index=False))


# =====================================================================
# Normalised counts
# =====================================================================

normalised_counts = pd.DataFrame(
    np.asarray(dds.layers["normed_counts"]),
    index=dds.obs_names.astype(str),
    columns=dds.var_names.astype(str),
)

normalised_counts.to_csv(output_dir / "DESeq2_normalised_pseudobulk_counts.csv")


# =====================================================================
# Variance-stabilised counts for PCA and heatmap
# =====================================================================

try:
    dds.vst(use_design=False)

    transformed_counts = pd.DataFrame(
        np.asarray(dds.layers["vst_counts"]),
        index=dds.obs_names.astype(str),
        columns=dds.var_names.astype(str),
    )

    transformation_name = "DESeq2 VST"

except Exception as error:
    warnings.warn(f"VST failed ({error}). Using log2 normalised counts instead.")

    transformed_counts = np.log2(normalised_counts + 1)

    transformation_name = "log2 normalised count + 1"

transformed_counts.to_csv(output_dir / "transformed_counts_for_visualisation.csv")


# =====================================================================
# Sample statistics
# =====================================================================

sample_metadata["raw_library_size"] = pseudobulk_counts.sum(axis=1)

sample_metadata["DESeq2_size_factor"] = pd.Series(
    np.asarray(dds.obs["size_factors"]),
    index=dds.obs_names.astype(str),
).reindex(sample_metadata.index)

sample_metadata.to_csv(output_dir / "sample_metadata_with_QC.csv")


# =====================================================================
# Summary statistics
# =====================================================================

n_up = int((significant_results["log2FoldChange"] > 0).sum())

n_down = int((significant_results["log2FoldChange"] < 0).sum())

# =====================================================================
# Plot 1: number of T cells per sample
# =====================================================================

cell_count_plot = (
    sample_metadata.reset_index()
    .rename(
        columns={
            sample_key: "sample",
        },
    )
    .sort_values(
        [
            "condition",
            "n_cells",
        ],
    )
)

fig, ax = plt.subplots(
    figsize=(
        max(9, 0.4 * cell_count_plot.shape[0]),
        5,
    ),
)

sns.barplot(
    data=cell_count_plot,
    x="sample",
    y="n_cells",
    hue="condition",
    dodge=False,
    ax=ax,
)

ax.set_title("Curated T cells contributing to each pseudobulk sample")

ax.set_xlabel("Biological sample")
ax.set_ylabel("Number of T cells")

ax.tick_params(
    axis="x",
    rotation=90,
)

ax.legend(
    title="Condition",
    frameon=False,
)

fig.tight_layout()

save_figure(
    fig,
    output_dir / "sample_T_cell_counts",
)


# =====================================================================
# Plot 2: pseudobulk library sizes
# =====================================================================

library_plot = (
    sample_metadata.reset_index()
    .rename(
        columns={
            sample_key: "sample",
        },
    )
    .sort_values(
        [
            "condition",
            "raw_library_size",
        ],
    )
)

fig, ax = plt.subplots(
    figsize=(
        max(9, 0.4 * library_plot.shape[0]),
        5,
    ),
)

sns.barplot(
    data=library_plot,
    x="sample",
    y="raw_library_size",
    hue="condition",
    dodge=False,
    ax=ax,
)

ax.set_yscale("log")

ax.set_title("Pseudobulk T-cell library size by sample")

ax.set_xlabel("Biological sample")
ax.set_ylabel("Total raw pseudobulk counts")

ax.tick_params(
    axis="x",
    rotation=90,
)

ax.legend(
    title="Condition",
    frameon=False,
)

fig.tight_layout()

save_figure(
    fig,
    output_dir / "pseudobulk_library_sizes",
)


# =====================================================================
# Plot 3: PCA
# =====================================================================

n_components = min(
    2,
    transformed_counts.shape[0] - 1,
    transformed_counts.shape[1],
)

if n_components == 2:
    pca = PCA(
        n_components=2,
        random_state=0,
    )

    pca_coordinates = pca.fit_transform(transformed_counts)

    pca_df = sample_metadata.copy()

    pca_df["PC1"] = pca_coordinates[:, 0]
    pca_df["PC2"] = pca_coordinates[:, 1]

    pca_df.to_csv(output_dir / "PCA_sample_coordinates.csv")

    fig, ax = plt.subplots(figsize=(7, 6))

    scatter_arguments = {
        "data": pca_df.reset_index(),
        "x": "PC1",
        "y": "PC2",
        "hue": "condition",
        "s": 90,
        "ax": ax,
    }

    if "platform" in pca_df.columns:
        scatter_arguments["style"] = "platform"

    sns.scatterplot(**scatter_arguments)

    for sample, row in pca_df.iterrows():
        ax.annotate(
            str(sample),
            (
                row["PC1"],
                row["PC2"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    explained = pca.explained_variance_ratio_ * 100

    ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)")

    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")

    ax.set_title(f"PCA of T-cell pseudobulk samples\n{transformation_name}")

    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    fig.tight_layout()

    save_figure(
        fig,
        output_dir / "pseudobulk_PCA",
    )


# =====================================================================
# Plot 4: volcano plot
# =====================================================================

volcano_df = results.dropna(
    subset=[
        "log2FoldChange",
        "pvalue",
    ],
).copy()

volcano_df["minus_log10_padj"] = -np.log10(
    volcano_df["padj"].fillna(1).clip(lower=1e-300),
)

fig, ax = plt.subplots(figsize=(8, 7))

colour_map = {
    "Upregulated": "red",
    "Downregulated": "blue",
    "Not significant": "grey",
}

for direction in ["Not significant", "Downregulated", "Upregulated"]:
    group = volcano_df[volcano_df["direction"].eq(direction)]

    ax.scatter(
        group["log2FoldChange"],
        group["minus_log10_padj"],
        c=colour_map[direction],
        s=12,
        alpha=0.6,
        label=direction,
    )


ax.axvline(
    -lfc_cutoff,
    linestyle="--",
    linewidth=0.8,
)

ax.axvline(
    lfc_cutoff,
    linestyle="--",
    linewidth=0.8,
)

ax.axhline(
    -np.log10(alpha),
    linestyle="--",
    linewidth=0.8,
)

target_volcano = volcano_df[volcano_df["feature_id"].eq(target_feature)]

for _, row in target_volcano.iterrows():
    ax.scatter(
        row["log2FoldChange"],
        row["minus_log10_padj"],
        c=colour_map[row["direction"]],
        s=180,
        marker=".",
        edgecolor="black",
        linewidth=0.8,
        zorder=10,
    )

    ax.annotate(
        "PIM3",
        (
            row["log2FoldChange"],
            row["minus_log10_padj"],
        ),
        xytext=(7, 7),
        textcoords="offset points",
        fontweight="bold",
    )

ax.set_title("Differential expression in T cells")

ax.set_xlabel("log2 fold change")

ax.set_ylabel("-log10 adjusted p-value")

ax.legend(
    frameon=False,
)

fig.tight_layout()

save_figure(
    fig,
    output_dir / "volcano_Crohns_vs_Normal",
)


# =====================================================================
# Plot 5: MA plot
# =====================================================================

ma_df = results.dropna(
    subset=[
        "baseMean",
        "log2FoldChange",
    ],
).copy()

fig, ax = plt.subplots(figsize=(8, 6))

not_sig = ma_df[~ma_df["significant_FDR_and_LFC"]]

sig = ma_df[ma_df["significant_FDR_and_LFC"]]

ax.scatter(
    np.log10(not_sig["baseMean"] + 1),
    not_sig["log2FoldChange"],
    s=10,
    c="grey",
    alpha=0.4,
    label="Not significant",
)

ax.scatter(
    np.log10(sig["baseMean"] + 1),
    sig["log2FoldChange"],
    s=14,
    c="red",
    alpha=0.75,
    label="Significant DEG",
)


ax.axhline(
    0,
    linewidth=0.8,
)

ax.axhline(
    lfc_cutoff,
    linestyle="--",
    linewidth=0.8,
)

ax.axhline(
    -lfc_cutoff,
    linestyle="--",
    linewidth=0.8,
)

target_ma = ma_df[ma_df["feature_id"].eq(target_feature)]

for _, row in target_ma.iterrows():
    x_value = np.log10(row["baseMean"] + 1)

    ax.scatter(
        x_value,
        row["log2FoldChange"],
        s=180,
        c=colour_map[row["direction"]],
        marker=".",
        edgecolor="black",
        linewidth=0.8,
        zorder=10,
    )

    ax.annotate(
        "PIM3",
        (
            x_value,
            row["log2FoldChange"],
        ),
        xytext=(7, 7),
        textcoords="offset points",
        fontweight="bold",
    )

ax.set_title("MA plot")

ax.set_xlabel("log10p1 normalised mean expression")

ax.set_ylabel("log2 fold change")

ax.legend(
    frameon=False,
)

fig.tight_layout()

save_figure(
    fig,
    output_dir / "MA_Crohns_vs_Normal",
)


# =====================================================================
# Plot 6: sample-level target expression
# =====================================================================

target_normalised_features = [
    feature for feature in target_features if feature in normalised_counts.columns
]

if target_normalised_features:
    target_expression = sample_metadata.copy()

    target_expression["target_normalised_count"] = normalised_counts[
        target_normalised_features
    ].sum(axis=1)

    target_expression["log2_target_normalised_count_plus1"] = np.log2(
        target_expression["target_normalised_count"] + 1,
    )

    target_raw_features = [
        feature for feature in target_features if feature in pseudobulk_counts.columns
    ]

    target_expression["target_raw_count"] = pseudobulk_counts[target_raw_features].sum(
        axis=1,
    )

    target_expression.to_csv(
        output_dir / f"{TARGET_ENSEMBL_ID}_sample_level_expression.csv",
    )

    target_condition_summary = target_expression.groupby(
        "condition",
        observed=True,
    ).agg(
        n_samples=(
            "target_normalised_count",
            "size",
        ),
        mean_normalised_count=(
            "target_normalised_count",
            "mean",
        ),
        median_normalised_count=(
            "target_normalised_count",
            "median",
        ),
        standard_deviation=(
            "target_normalised_count",
            "std",
        ),
        mean_log2_normalised_count=(
            "log2_target_normalised_count_plus1",
            "mean",
        ),
    )

    target_condition_summary.to_csv(
        output_dir / f"{TARGET_ENSEMBL_ID}_expression_summary_by_condition.csv",
    )

    fig, ax = plt.subplots(figsize=(6, 5.5))

    sns.boxplot(
        data=target_expression,
        x="condition",
        y="log2_target_normalised_count_plus1",
        order=[
            "Normal",
            "Crohns",
        ],
        showfliers=False,
        ax=ax,
    )

    sns.stripplot(
        data=target_expression,
        x="condition",
        y="log2_target_normalised_count_plus1",
        order=[
            "Normal",
            "Crohns",
        ],
        jitter=0.16,
        size=7,
        ax=ax,
    )

    ax.set_title("PIM3 expression in curated T cells")

    ax.set_xlabel("Condition")

    ax.set_ylabel("log2(normalised pseudobulked count)")

    fig.tight_layout()

    save_figure(
        fig,
        output_dir / f"{TARGET_ENSEMBL_ID}_expression_Crohns_vs_Normal",
    )


# =====================================================================
# Plot 7: heatmap of leading DE genes
# =====================================================================

ranked_for_heatmap = results.dropna(
    subset=[
        "padj",
        "log2FoldChange",
    ],
).sort_values(
    [
        "significant_FDR_and_LFC",
        "padj",
        "pvalue",
    ],
    ascending=[
        False,
        True,
        True,
    ],
)

heatmap_features = (
    ranked_for_heatmap["feature_id"].drop_duplicates().head(top_n_heatmap).tolist()
)

heatmap_features = [
    feature for feature in heatmap_features if feature in transformed_counts.columns
]

if len(heatmap_features) >= 2:
    heatmap_data = transformed_counts[heatmap_features].T

    gene_means = heatmap_data.mean(axis=1)

    gene_stds = heatmap_data.std(
        axis=1,
        ddof=0,
    ).replace(0, np.nan)

    heatmap_z = heatmap_data.sub(gene_means, axis=0).div(gene_stds, axis=0).fillna(0)

    ordered_samples = (
        sample_metadata.reset_index()
        .sort_values(
            [
                "condition",
                sample_key,
            ],
        )[sample_key]
        .astype(str)
        .tolist()
    )

    heatmap_z = heatmap_z[ordered_samples]

    # Heatmap rows are labelled with Ensembl feature IDs only.
    heatmap_z.index = heatmap_z.index.astype(str)

    heatmap_z.to_csv(output_dir / "top_DE_genes_heatmap_Z_scores.csv")

    fig_width = max(
        9,
        0.45 * len(ordered_samples),
    )

    fig_height = max(
        8,
        0.28 * len(heatmap_features),
    )

    fig, ax = plt.subplots(
        figsize=(
            fig_width,
            fig_height,
        ),
    )

    sns.heatmap(
        heatmap_z,
        cmap="vlag",
        center=0,
        linewidths=0.2,
        cbar_kws={"label": (f"Row Z-score of {transformation_name}")},
        ax=ax,
    )

    n_normal = int(
        (
            sample_metadata.loc[
                ordered_samples,
                "condition",
            ]
            == "Normal"
        ).sum(),
    )

    if 0 < n_normal < len(ordered_samples):
        ax.axvline(
            n_normal,
            linewidth=1.5,
        )

    ax.set_title("Leading Crohn's versus Normal T-cell DE genes")

    ax.set_xlabel("Biological sample")

    ax.set_ylabel("Gene")

    ax.tick_params(
        axis="x",
        rotation=90,
    )

    ax.tick_params(
        axis="y",
        labelsize=7,
    )

    fig.tight_layout()

    save_figure(
        fig,
        output_dir / "top_DE_genes_heatmap",
    )
