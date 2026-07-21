"""
Part 2: Shortest Path Lengths & Wilcoxon Test
"""

import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats
from itertools import combinations

# ── File Paths (relative to repo root) ──────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
PPI_FILE   = os.path.join(REPO_ROOT, "data", "Human-PPI.txt")
LIST1_FILE = os.path.join(REPO_ROOT, "data", "protein-list1.txt")
LIST2_FILE = os.path.join(REPO_ROOT, "data", "protein-list2.txt")
OUT_DIR    = os.path.join(REPO_ROOT, "results", "part2_shortest_paths")
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

# ── Build Graph ──────────────────────────────────────────
print("Loading PPI network...")
edges = pd.read_csv(PPI_FILE, sep="\t", skiprows=1, header=None, names=["A", "B"])
G = nx.Graph()
G.add_edges_from(zip(edges["A"], edges["B"]))
G.remove_edges_from(nx.selfloop_edges(G))
print(f"  Nodes: {G.number_of_nodes():,}  Edges: {G.number_of_edges():,}")

# ── Load Protein Lists ───────────────────────────────────
with open(LIST1_FILE) as f:
    list1 = [p.strip() for p in f if p.strip()]
with open(LIST2_FILE) as f:
    list2 = [p.strip() for p in f if p.strip()]

# ── Shortest Path Function ───────────────────────────────
def get_shortest_paths(G, proteins, label):
    in_net    = [p for p in proteins if p in G]
    not_found = [p for p in proteins if p not in G]
    print(f"\n{label}: {len(proteins)} proteins | {len(in_net)} in network | {len(not_found)} missing")
    if not_found:
        print(f"  Missing: {not_found}")
    paths = []
    for p1, p2 in combinations(in_net, 2):
        try:
            paths.append({"Protein1": p1, "Protein2": p2,
                          "ShortestPath": nx.shortest_path_length(G, p1, p2)})
        except nx.NetworkXNoPath:
            paths.append({"Protein1": p1, "Protein2": p2, "ShortestPath": np.nan})
    return pd.DataFrame(paths), in_net, not_found

print("\nCalculating shortest paths...")
df1, in1, miss1 = get_shortest_paths(G, list1, "Protein-List1")
df2, in2, miss2 = get_shortest_paths(G, list2, "Protein-List2")

df1.to_csv(os.path.join(OUT_DIR, "shortest_paths_list1.csv"), index=False)
df2.to_csv(os.path.join(OUT_DIR, "shortest_paths_list2.csv"), index=False)

sp1 = df1["ShortestPath"].dropna().values
sp2 = df2["ShortestPath"].dropna().values

# ── Wilcoxon Rank-Sum Test ───────────────────────────────
stat, p = stats.mannwhitneyu(sp1, sp2, alternative="two-sided")
sig = p < 0.05

result_txt = (
    f"Wilcoxon Rank-Sum (Mann-Whitney U) Test\n"
    f"========================================\n"
    f"List1: n={len(sp1)}, mean={sp1.mean():.3f}, median={np.median(sp1):.1f}\n"
    f"List2: n={len(sp2)}, mean={sp2.mean():.3f}, median={np.median(sp2):.1f}\n"
    f"U statistic : {stat:.2f}\n"
    f"p-value     : {p:.4e}\n"
    f"Significant : {'YES' if sig else 'NO'} (alpha=0.05)\n"
)
print("\n" + result_txt)
with open(os.path.join(OUT_DIR, "wilcoxon_test_result.txt"), "w", encoding="utf-8") as f:
    f.write(result_txt)

# ── FIGURE ───────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10), facecolor="#0d1117")
fig.suptitle("Shortest Path Length Analysis — Protein Set Comparison",
             fontsize=16, fontweight="bold", color="white", y=0.98)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32,
                       left=0.07, right=0.97, top=0.93, bottom=0.08)

bins = range(1, int(max(sp1.max(), sp2.max())) + 2)

# ── Panel A: Overlapping histogram ───────────────────────
ax_a = fig.add_subplot(gs[0, 0])
ax_a.hist(sp1, bins=bins, density=True, alpha=0.6, color="#5ba3e0",
          label=f"List1 (n={len(sp1)})", edgecolor="#0d1117")
ax_a.hist(sp2, bins=bins, density=True, alpha=0.6, color="#e05b7a",
          label=f"List2 (n={len(sp2)})", edgecolor="#0d1117")
ax_a.set_xlabel("Shortest Path Length", fontsize=11)
ax_a.set_ylabel("Density", fontsize=11)
ax_a.set_title("A.  Path Length Distribution (Overlapping)", fontsize=12,
               color="white", pad=8)
ax_a.legend(fontsize=9, facecolor="#1a1a2e", edgecolor="#444", labelcolor="white")
ax_a.grid(True)

# ── Panel B: Side-by-side bar chart ──────────────────────
ax_b = fig.add_subplot(gs[0, 1])
unique_lengths = sorted(set(sp1.astype(int)) | set(sp2.astype(int)))
c1 = [np.sum(sp1 == l) / len(sp1) for l in unique_lengths]
c2 = [np.sum(sp2 == l) / len(sp2) for l in unique_lengths]
x  = np.array(unique_lengths)
w  = 0.35
ax_b.bar(x - w/2, c1, width=w, color="#5ba3e0", alpha=0.85, label=f"List1")
ax_b.bar(x + w/2, c2, width=w, color="#e05b7a", alpha=0.85, label=f"List2")
ax_b.set_xlabel("Shortest Path Length", fontsize=11)
ax_b.set_ylabel("Fraction of pairs", fontsize=11)
ax_b.set_title("B.  Path Length Distribution (Side-by-Side)", fontsize=12,
               color="white", pad=8)
ax_b.set_xticks(x)
ax_b.legend(fontsize=9, facecolor="#1a1a2e", edgecolor="#444", labelcolor="white")
ax_b.grid(True)

# ── Panel C: Boxplot ─────────────────────────────────────
ax_c = fig.add_subplot(gs[1, 0])
bp = ax_c.boxplot([sp1, sp2], patch_artist=True, widths=0.5,
                  medianprops=dict(color="white", lw=2),
                  whiskerprops=dict(color="#888"),
                  capprops=dict(color="#888"),
                  flierprops=dict(marker="o", color="#888", alpha=0.3, markersize=3))
bp["boxes"][0].set_facecolor("#5ba3e0"); bp["boxes"][0].set_alpha(0.75)
bp["boxes"][1].set_facecolor("#e05b7a"); bp["boxes"][1].set_alpha(0.75)
ax_c.set_xticks([1, 2])
ax_c.set_xticklabels([f"List1\n(n={len(sp1)})", f"List2\n(n={len(sp2)})"], fontsize=10)
ax_c.set_ylabel("Shortest Path Length", fontsize=11)
ax_c.set_title("C.  Boxplot Comparison", fontsize=12, color="white", pad=8)
ax_c.text(0.5, 0.95, f"p = {p:.4e}  ({'*' if sig else 'ns'})",
          transform=ax_c.transAxes, ha="center", va="top",
          fontsize=11, color="orange" if sig else "#aaa")
ax_c.grid(True, axis="y")

# ── Panel D: Summary table ────────────────────────────────
ax_d = fig.add_subplot(gs[1, 1])
ax_d.axis("off")
ax_d.set_title("D.  Wilcoxon Test Summary", fontsize=12, color="white", pad=8)

rows = [
    ("List1 proteins",        f"{len(list1)}"),
    ("List1 in network",      f"{len(in1)}"),
    ("List1 missing",         f"{len(miss1)}  {miss1 if miss1 else ''}"),
    ("List2 proteins",        f"{len(list2)}"),
    ("List2 in network",      f"{len(in2)}"),
    ("List2 missing",         f"{len(miss2)}  {miss2 if miss2 else ''}"),
    ("List1 mean path",       f"{sp1.mean():.3f}"),
    ("List2 mean path",       f"{sp2.mean():.3f}"),
    ("List1 median path",     f"{np.median(sp1):.1f}"),
    ("List2 median path",     f"{np.median(sp2):.1f}"),
    ("U statistic",           f"{stat:.2f}"),
    ("p-value",               f"{p:.4e}"),
    ("Significant (a=0.05)",  "YES" if sig else "NO"),
]

y = 0.95
for label, value in rows:
    ax_d.text(0.02, y, label,  transform=ax_d.transAxes,
              fontsize=9.5, color="#888888", va="top")
    color = "orange" if label in ("p-value", "Significant (a=0.05)") else "white"
    ax_d.text(0.98, y, value, transform=ax_d.transAxes,
              fontsize=9.5, color=color, va="top", ha="right", fontweight="bold")
    ax_d.plot([0.02, 0.98], [y - 0.008, y - 0.008], color='#333', lw=0.5,
             transform=ax_d.transAxes, clip_on=False)
    y -= 0.068

plt.savefig(os.path.join(OUT_DIR, "part2_shortest_paths.png"), dpi=150,
            facecolor="#0d1117", bbox_inches="tight")
plt.close()
print("  Saved: part2_shortest_paths.png")
print(f"\nAll Part 2 outputs saved to: {OUT_DIR}")