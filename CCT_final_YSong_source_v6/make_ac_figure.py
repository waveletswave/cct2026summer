"""A_c resolution-portability figure. All values from standalone_CA/docs/validation_log.md."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

INK, SEC, MUT = "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
BLUE, BLUE_L, RED = "#2a78d6", "#cde2fb", "#d03b3b"

band28, band10 = (0.040, 0.079), (0.050, 0.10)
overlap, chosen, naive10 = (0.050, 0.079), 0.04757, 0.006

fig, ax = plt.subplots(figsize=(9.6, 3.6), dpi=200)
fig.patch.set_facecolor("white"); ax.set_facecolor("white")

ax.axvspan(*overlap, color=BLUE_L, alpha=0.55, lw=0, zorder=0)
ax.text(sum(overlap)/2 + 0.0045, 1.55, "bands overlap\n0.050 to 0.079 km²", ha="center", va="center",
        fontsize=9.5, color=SEC, linespacing=1.3)

for band, y in ((band28, 1.0), (band10, 0.0)):
    ax.plot(band, [y, y], color=BLUE, lw=9, solid_capstyle="round", zorder=3)
    ax.plot(band[0], y, "o", color=BLUE, ms=11, mec="white", mew=2, zorder=4)
ax.text(band28[1] + 0.0025, 1.0, "passing band (|t| < 2)", va="center", fontsize=9.5, color=SEC)
ax.annotate("objective A$_c$ = first sustained pass\n(0.0396 km² @ 28 m,  0.0500 km² @ 10 m)",
            xy=(band28[0] - 0.0012, 1.06), xytext=(0.0145, 1.75), fontsize=9.5, color=SEC,
            ha="center", linespacing=1.35,
            arrowprops=dict(arrowstyle="-", color=MUT, lw=0.9, shrinkB=3,
                            connectionstyle="arc3,rad=-0.18"))

ax.axvline(chosen, color=INK, lw=1.4, ls=(0, (4, 3)), zorder=5)
ax.text(0.0495, 2.62, "chosen A$_c$ = 0.0476 km²", ha="left", va="center",
        fontsize=11, color=INK, fontweight="bold")
ax.text(0.0495, 2.22, "60 cells @ 28 m  =  476 cells @ 10 m", ha="left", va="center",
        fontsize=9.5, color=SEC)

ax.plot(naive10, 0.0, "x", color=RED, ms=11, mew=3, zorder=5)
ax.text(naive10 - 0.004, -0.78, "hold the cell count instead:  60 cells @ 10 m = 0.006 km², outside every passing band",
        ha="left", va="center", fontsize=9.5, color=RED)
arrow = FancyArrowPatch((chosen - 0.002, 0.0), (naive10 + 0.002, 0.0),
                        arrowstyle="-|>", mutation_scale=14, color=RED, lw=1.2, ls=":", zorder=2)
ax.add_patch(arrow)

ax.set_yticks([0.0, 1.0], labels=["10 m grid", "28 m grid"], fontsize=11.5, color=INK)
ax.set_xlim(0, 0.112); ax.set_ylim(-1.15, 2.95)
ax.set_xlabel("channel-initiation support area  A$_c$  (km²)", fontsize=11, color=SEC)
ax.xaxis.set_tick_params(labelsize=10, colors=MUT)
ax.grid(axis="x", color=GRID, lw=0.8, zorder=-1)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(BASE)
ax.tick_params(axis="y", length=0)

fig.tight_layout()
fig.savefig("assets/ac_resolution_bands.png", bbox_inches="tight", facecolor="white")
print("saved")
