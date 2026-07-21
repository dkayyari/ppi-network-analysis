"""
Part 1: Degree, Clustering Coefficient, Average Clustering, Scale-Free Test
"""

import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats

# ── File Paths (relative to repo root) ──────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
PPI_FILE = os.path.join(REPO_ROOT, "data", "Human-PPI.txt")
OUT_DIR  = os.path.join(REPO_ROOT, "results", "part1_network_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Dark theme ───────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor":  "#0d1117",
    "axes.facecolor":    "#0d1117",
    "axes.edgecolor":    "#444",
    "axes.labelcolor":   "#cccccc",
    "xtick.color":       "#cccccc",
    "ytick.color":       "#cccccc",
    "text.color":        "#cccccc",
    "grid.color":        "#222",
    "grid.linestyle":    "--",
    "grid.alpha":        0.4,
    "font.family":       "DejaVu Sans",
})

# ── 1. Build Graph ───────────────────────────────────────
print("Loading PPI network...")
edges = pd.read_csv(PPI_FILE, sep="\t", skiprows=1, header=None, names=["A", "B"])
G = nx.Graph()
G.add_edges_from(zip(edges["A"], edges["B"]))
G.remove_edges_from(nx.selfloop_edges(G))
print(f"  Nodes: {G.number_of_nodes():,}  Edges: {G.number_of_edges():,}")

# ── 2. Degree ────────────────────────────────────────────
degree_dict = dict(G.degree())
df_degree = (pd.DataFrame(degree_dict.items(), columns=["Protein", "Degree"])
               .sort_values("Degree", ascending=False))
df_degree.to_csv(os.path.join(OUT_DIR, "node_degrees.csv"), index=False)

# ── 3 & 4. Clustering ────────────────────────────────────
print("Calculating clustering coefficients...")
cc_dict = nx.clustering(G)
df_cc = (pd.DataFrame(cc_dict.items(), columns=["Protein", "ClusteringCoeff"])
           .sort_values("ClusteringCoeff", ascending=False))
df_cc.to_csv(os.path.join(OUT_DIR, "node_clustering.csv"), index=False)
avg_cc = nx.average_clustering(G)

# ── 5. Degree distribution data ──────────────────────────
degrees    = np.array([d for _, d in G.degree()])
deg_vals, deg_counts = np.unique(degrees[degrees > 0], return_counts=True)
pk         = deg_counts / deg_counts.sum()   # probability

# Power-law fit on log-log
mask = (deg_vals > 0) & (pk > 0)
log_k  = np.log10(deg_vals[mask].astype(float))
log_pk = np.log10(pk[mask])
slope, intercept, r, _, _ = stats.linregress(log_k, log_pk)
gamma  = -slope
fit_line = 10 ** (intercept + slope * log_k)

# Degree vs clustering
deg_arr = np.array([degree_dict[n] for n in G.nodes()])
cc_arr  = np.array([cc_dict[n]     for n in G.nodes()])

# Bin average for degree vs clustering
log_bins = np.logspace(np.log10(max(deg_arr.min(),1)), np.log10(deg_arr.max()), 30)
bin_idx  = np.digitize(deg_arr, log_bins)
bin_deg, bin_cc = [], []
for i in range(1, len(log_bins)):
    sel = cc_arr[bin_idx == i]
    if len(sel) > 0:
        bin_deg.append(log_bins[i])
        bin_cc.append(sel.mean())

# Network stats
n_components = nx.number_connected_components(G)
largest_cc   = max(nx.connected_components(G), key=len)

# ── FIGURE ───────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10), facecolor="#0d1117")
fig.suptitle("Human Protein–Protein Interaction Network Analysis",
             fontsize=16, fontweight="bold", color="white", y=0.98)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32,
                       left=0.07, right=0.97, top=0.93, bottom=0.08)

# ── Panel A: Linear degree distribution ─────────────────
ax_a = fig.add_subplot(gs[0, 0])
ax_a.bar(deg_vals, pk, width=1, color="#3a7ebf", alpha=0.75)
ax_a.set_xlabel("Degree k", fontsize=11)
ax_a.set_ylabel("P(k)  — Fraction of nodes", fontsize=11)
ax_a.set_title("A.  Degree Distribution  (linear scale)", fontsize=12,
               color="white", pad=8)
ax_a.grid(True)

# ── Panel B: Log-log + power-law ────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
ax_b.scatter(deg_vals[mask], pk[mask], s=22, color="#5ba3e0", alpha=0.8,
             label="Observed P(k)", zorder=3)
ax_b.plot(deg_vals[mask], fit_line, color="orange", lw=2,
          label=f"Power-law fit\n\u03b3 = {gamma:.2f}, $R^2$ = {r**2:.3f}")
ax_b.set_xscale("log"); ax_b.set_yscale("log")
ax_b.set_xlabel("Degree k (log scale)", fontsize=11)
ax_b.set_ylabel("P(k) (log scale)", fontsize=11)
ax_b.set_title("B.  Log-Log Degree Distribution + Power-Law Fit", fontsize=12,
               color="white", pad=8)
ax_b.legend(fontsize=9, facecolor="#1a1a2e", edgecolor="#444", labelcolor="white")
ax_b.grid(True)

# ── Panel C: Degree vs Clustering ────────────────────────
ax_c = fig.add_subplot(gs[1, 0])
ax_c.scatter(deg_arr, cc_arr, s=4, color="#2ecc71", alpha=0.25,
             label="Individual nodes", zorder=2)
ax_c.plot(bin_deg, bin_cc, color="orange", lw=2.5,
          label="Bin average", zorder=3)
ax_c.set_xscale("log")
ax_c.set_xlabel("Degree k (log scale)", fontsize=11)
ax_c.set_ylabel("Clustering Coefficient C(k)", fontsize=11)
ax_c.set_title("C.  Degree vs. Clustering Coefficient", fontsize=12,
               color="white", pad=8)
ax_c.legend(fontsize=9, facecolor="#1a1a2e", edgecolor="#444", labelcolor="white")
ax_c.grid(True)

# ── Panel D: Network Summary table ───────────────────────
ax_d = fig.add_subplot(gs[1, 1])
ax_d.axis("off")
ax_d.set_title("D.  Network Summary", fontsize=12, color="white", pad=8)

rows = [
    ("Nodes (proteins)",        f"{G.number_of_nodes():,}"),
    ("Edges (interactions)",    f"{G.number_of_edges():,}"),
    ("Connected components",    f"{n_components}"),
    ("Largest component",       f"{len(largest_cc):,} nodes"),
    ("Max degree",              f"{deg_arr.max()}"),
    ("Mean degree",             f"{deg_arr.mean():.2f}"),
    ("Median degree",           f"{int(np.median(deg_arr))}"),
    ("Avg. clustering coeff.",  f"{avg_cc:.6f}"),
    ("Power-law exponent \u03b3", f"{gamma:.4f}"),
    ("R\u00b2 (log-log fit)",   f"{r**2:.4f}"),
    ("Scale-free?",             f"Approx  (\u03b3 = {gamma:.2f})"),
]

y = 0.92
for label, value in rows:
    ax_d.text(0.02, y, label,  transform=ax_d.transAxes,
              fontsize=10.5, color="#888888", va="top")
    ax_d.text(0.98, y, value, transform=ax_d.transAxes,
              fontsize=10.5, color="white", va="top", ha="right", fontweight="bold")
    ax_d.plot([0.02, 0.98], [y - 0.01, y - 0.01], color='#333', lw=0.5,
             transform=ax_d.transAxes, clip_on=False)
    y -= 0.082

plt.savefig(os.path.join(OUT_DIR, "part1_network_analysis.png"), dpi=150,
            facecolor="#0d1117", bbox_inches="tight")
plt.close()
print("  Saved: part1_network_analysis.png")

# ── Network Summary text ─────────────────────────────────
summary = (
    f"Human PPI Network Summary\n"
    f"==========================\n"
    f"Nodes             : {G.number_of_nodes():,}\n"
    f"Edges             : {G.number_of_edges():,}\n"
    f"Connected comps   : {n_components}\n"
    f"Largest component : {len(largest_cc):,} nodes\n"
    f"Avg degree        : {deg_arr.mean():.2f}\n"
    f"Max degree        : {deg_arr.max()}\n"
    f"Median degree     : {int(np.median(deg_arr))}\n"
    f"Avg clustering    : {avg_cc:.6f}\n"
    f"Power-law gamma   : {gamma:.4f}  (R2={r**2:.4f})\n"
    f"Scale-free?       : {'Likely YES' if 2 <= gamma <= 3 else 'Approx (gamma = ' + str(round(gamma,2)) + ')'}\n"
)
with open(os.path.join(OUT_DIR, "network_summary.txt"), "w", encoding="utf-8") as f:
    f.write(summary)

print("\n" + summary)
print(f"All Part 1 outputs saved to: {OUT_DIR}")