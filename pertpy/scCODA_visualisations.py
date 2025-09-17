import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label" 
OUTPATH = PATH + "/../../pertpy_output"

ref_name = "Mesenchymal"

# Load data
df = pd.read_csv(os.path.join(OUTPATH, "category_log2fc_df.csv"))
df = df.sort_values("log2FC_mean").reset_index(drop=True)
y = np.arange(len(df))

# Colors (alphabetical order)
cats = df["Cell Type"].astype(str).values
unique_cats = sorted(set(cats))
cmap = plt.get_cmap("tab20")
color_map = {c: cmap(i % cmap.N) for i, c in enumerate(unique_cats)}
colors = [color_map[c] for c in cats]

# Point sizes from inclusion probability
if "Inclusion probability" in df.columns:
    incl = df["Inclusion probability"].to_numpy()
    sizes = 300 * incl  
else:
    incl = np.ones(len(df))
    sizes = np.full(len(df), 50)

# Error bars
x = df["log2FC_mean"].to_numpy()
x_lo = df["log2FC_HDI_low"].to_numpy()
x_hi = df["log2FC_HDI_high"].to_numpy()

# Figure
fig, ax = plt.subplots(figsize=(8, max(5, 0.3*len(df))))

for i in range(len(df)):
    ax.errorbar(
        x[i], y[i],
        xerr=[[x[i] - x_lo[i]], [x_hi[i] - x[i]]],
        fmt='o',
        elinewidth=1.2, capsize=3,
        color=colors[i], ecolor=colors[i],
        markersize=np.sqrt(sizes[i]) 
    )

ax.axvline(0, ls='--', lw=1, color='gray')

# Cosmetics
ax.set_yticks(y)
ax.set_yticklabels(df["Cell Type"])
ax.set_xlabel("Log2 fold change (HDI)")
ax.set_title("Differential composition")
ax.grid(axis='x', alpha=0.2)

fig.tight_layout()

# Save
out_pdf = os.path.join(OUTPATH, "category_log2FC.pdf")
fig.savefig(out_pdf, bbox_inches="tight")


# Load data
df = pd.read_csv(os.path.join(OUTPATH, "category_effect_df.csv"))
df = df.sort_values("Final Parameter").reset_index(drop=True)
y = np.arange(len(df))

# Colors (alphabetical order)
cats = df["Cell Type"].astype(str).values
unique_cats = sorted(set(cats))
cmap = plt.get_cmap("tab20")
color_map = {c: cmap(i % cmap.N) for i, c in enumerate(unique_cats)}
colors = [color_map[c] for c in cats]

# Point sizes from inclusion probability
if "Inclusion probability" in df.columns:
    incl = df["Inclusion probability"].to_numpy()
    sizes = 300 * incl   # scale factor for visibility
else:
    incl = np.ones(len(df))
    sizes = np.full(len(df), 50)

# Error bars
x = df["Final Parameter"].to_numpy()
x_lo = df["HDI 3%"].to_numpy()
x_hi = df["HDI 97%"].to_numpy()

# Figure
fig, ax = plt.subplots(figsize=(8, max(5, 0.3*len(df))))

for i in range(len(df)):
    ax.errorbar(
        x[i], y[i],
        xerr=[[x[i] - x_lo[i]], [x_hi[i] - x[i]]],
        fmt='o',
        elinewidth=1.2, capsize=3,
        color=colors[i], ecolor=colors[i],
        markersize=np.sqrt(sizes[i]) 
    )

ax.axvline(0, ls='--', lw=1, color='gray')

# Cosmetics
ax.set_yticks(y)
ax.set_yticklabels(df["Cell Type"])
ax.set_xlabel(r"$\tilde{\beta}$ (HDI)")
ax.set_title("Differential composition")
ax.grid(axis='x', alpha=0.2)

ref_y = df.index[df["Cell Type"] == ref_name][0]  # row index in df
ax.text(
    0, ref_y, "*", color="black", fontsize=20,
    va="center", ha="center", fontweight="bold"
)

# Legend for inclusion probability (sizes)
for prob in np.linspace(0.1, 1.0, 10):
    plt.scatter([], [], s=300*prob, c="gray", alpha=0.6, label=f"{prob:.2f}")

leg2 = ax.legend(title="Inclusion probability", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

fig.tight_layout()

# Save
out_pdf = os.path.join(OUTPATH, "category_coef.pdf")
fig.savefig(out_pdf, bbox_inches="tight")



ref_name = "Mast cell"

# Load data
df = pd.read_csv(os.path.join(OUTPATH, "curated_log2fc_df.csv"))
df = df.sort_values("log2FC_mean").reset_index(drop=True)
y = np.arange(len(df))

# Colors (alphabetical order)
cats = df["Cell Type"].astype(str).values
unique_cats = sorted(set(cats))
cmap = sns.color_palette("hls", len(unique_cats))
color_map = {ct: cmap[i] for i, ct in enumerate(unique_cats)}
colors = [color_map[c] for c in cats]

# Point sizes from inclusion probability
if "Inclusion probability" in df.columns:
    incl = df["Inclusion probability"].to_numpy()
    sizes = 300 * incl  
else:
    incl = np.ones(len(df))
    sizes = np.full(len(df), 50)

# Error bars
x = df["log2FC_mean"].to_numpy()
x_lo = df["log2FC_HDI_low"].to_numpy()
x_hi = df["log2FC_HDI_high"].to_numpy()

# Figure
fig, ax = plt.subplots(figsize=(8, max(5, 0.3*len(df))))

for i in range(len(df)):
    ax.errorbar(
        x[i], y[i],
        xerr=[[x[i] - x_lo[i]], [x_hi[i] - x[i]]],
        fmt='o',
        elinewidth=1.2, capsize=3,
        color=colors[i], ecolor=colors[i],
        markersize=np.sqrt(sizes[i]) 
    )

ax.axvline(0, ls='--', lw=1, color='gray')

# Cosmetics
ax.set_yticks(y)
ax.set_yticklabels(df["Cell Type"])

ax.set_ylim(-0.5, len(df) - 0.5)
ax.margins(y=0)  

ax.set_xlabel("Log2 fold change (HDI)")
ax.set_title("Differential composition")
ax.grid(axis='x', alpha=0.2)

fig.tight_layout()

# Save
out_pdf = os.path.join(OUTPATH, "curated_log2FC.pdf")
fig.savefig(out_pdf, bbox_inches="tight")


# Load data
df = pd.read_csv(os.path.join(OUTPATH, "curated_effect_df.csv"))
df = df.sort_values("Final Parameter").reset_index(drop=True)
y = np.arange(len(df))

# Colors (alphabetical order)
cats = df["Cell Type"].astype(str).values
unique_cats = sorted(set(cats))
cmap = sns.color_palette("hls", len(unique_cats))
color_map = {ct: cmap[i] for i, ct in enumerate(unique_cats)}
colors = [color_map[c] for c in cats]

# Point sizes from inclusion probability
if "Inclusion probability" in df.columns:
    incl = df["Inclusion probability"].to_numpy()
    sizes = 300 * incl   # scale factor for visibility
else:
    incl = np.ones(len(df))
    sizes = np.full(len(df), 50)

# Error bars
x = df["Final Parameter"].to_numpy()
x_lo = df["HDI 3%"].to_numpy()
x_hi = df["HDI 97%"].to_numpy()

# Figure
fig, ax = plt.subplots(figsize=(8, max(5, 0.3*len(df))))

for i in range(len(df)):
    ax.errorbar(
        x[i], y[i],
        xerr=[[x[i] - x_lo[i]], [x_hi[i] - x[i]]],
        fmt='o',
        elinewidth=1.2, capsize=3,
        color=colors[i], ecolor=colors[i],
        markersize=np.sqrt(sizes[i]) 
    )

ax.axvline(0, ls='--', lw=1, color='gray')

# Cosmetics
ax.set_yticks(y)
ax.set_yticklabels(df["Cell Type"])
ax.set_ylim(-0.5, len(df) - 0.5)
ax.margins(y=0)  

ax.set_xlabel(r"$\tilde{\beta}$ (HDI)")
ax.set_title("Differential composition")
ax.grid(axis='x', alpha=0.2)

ref_y = df.index[df["Cell Type"] == ref_name][0]  # row index in df
ax.text(
    0, ref_y, "*", color="black", fontsize=20,
    va="center", ha="center", fontweight="bold"
)

# Legend for inclusion probability (sizes)
for prob in np.linspace(0.1, 1.0, 10):
    plt.scatter([], [], s=300*prob, c="gray", alpha=0.6, label=f"{prob:.2f}")

leg2 = ax.legend(title="Inclusion probability", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

fig.tight_layout()

# Save
out_pdf = os.path.join(OUTPATH, "curated_coef.pdf")
fig.savefig(out_pdf, bbox_inches="tight")