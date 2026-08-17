"""Synthetic-terrain teaching assets: D8 flow, accumulation, threshold sweep.
Pure numpy/scipy — no GIS deps. Terrain is labeled synthetic on every frame.
"""
import heapq
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

rng = np.random.default_rng(11)
INK, SEC, MUT = "#0b0b0b", "#52514e", "#898781"
BLUE, BLUE_D, BLUE_L = "#2a78d6", "#0d366b", "#cde2fb"
CELL = 30.0  # m, nominal
OUT = "assets"

# ---------------- terrain ----------------
N = 210
def octave(sigma):
    f = gaussian_filter(rng.random((N, N)), sigma)
    return (f - f.min()) / (f.max() - f.min())
base = 1.0 * octave(34) + 0.55 * octave(14) + 0.22 * octave(6)
base = (base - base.min()) / (base.max() - base.min())
yy, xx = np.mgrid[0:N, 0:N]
tilt = (N - yy) / N * 0.16            # gentle regional slope toward the south edge
dem = base + tilt
dem = (dem - dem.min()) / (dem.max() - dem.min())
dem = 200 + 600 * dem                  # 200–800 m
dem = dem[5:-5, 5:-5]                  # trim border artifacts
NR, NC = dem.shape

# ---------------- priority-flood depression filling ----------------
filled = dem.copy()
visited = np.zeros_like(dem, bool)
heap = []
EPS = 1e-3
for i in range(NR):
    for j in (0, NC - 1):
        heapq.heappush(heap, (dem[i, j], i, j)); visited[i, j] = True
for j in range(NC):
    for i in (0, NR - 1):
        if not visited[i, j]:
            heapq.heappush(heap, (dem[i, j], i, j)); visited[i, j] = True
NB = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
while heap:
    z, i, j = heapq.heappop(heap)
    for di, dj in NB:
        ni, nj = i + di, j + dj
        if 0 <= ni < NR and 0 <= nj < NC and not visited[ni, nj]:
            visited[ni, nj] = True
            filled[ni, nj] = max(dem[ni, nj], z + EPS)
            heapq.heappush(heap, (filled[ni, nj], ni, nj))

# ---------------- D8 receivers + accumulation ----------------
filled = filled + rng.random(filled.shape) * 5e-3   # break ties on filled flats
dist = np.array([np.hypot(di, dj) for di, dj in NB])
recv = -np.ones((NR, NC), int)               # flat index of receiver
for i in range(NR):
    for j in range(NC):
        best, bi = 0.0, -1
        for k, (di, dj) in enumerate(NB):
            ni, nj = i + di, j + dj
            if 0 <= ni < NR and 0 <= nj < NC:
                s = (filled[i, j] - filled[ni, nj]) / dist[k]
                if s > best:
                    best, bi = s, ni * NC + nj
        recv[i, j] = bi                       # -1 = edge outlet

order = np.argsort(filled, axis=None)[::-1]   # high -> low
acc = np.ones((NR, NC))
fr = acc.ravel(); rr = recv.ravel()
for idx in order:
    r = rr[idx]
    if r >= 0:
        fr[r] += fr[idx]
acc = fr.reshape(NR, NC)

# ---------------- rendering helpers ----------------
ls = LightSource(azdeg=315, altdeg=45)
hs = ls.hillshade(dem, vert_exag=3, dx=CELL, dy=CELL)   # raw DEM: no flat "lakes"
area_km2 = NR * NC * CELL * CELL / 1e6

def canvas():
    # FIXED canvas (no tight bbox) so every frame in a set has identical pixel dims
    fig, ax = plt.subplots(figsize=(6.4, 7.05), dpi=170)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.862, bottom=0.015)
    fig.patch.set_facecolor("white")
    ax.imshow(hs, cmap="gray", vmin=0.1, vmax=1.05, interpolation="bilinear")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#e1e0d9")
    ax.text(0.012, 0.014, "synthetic terrain, for illustration", transform=ax.transAxes,
            fontsize=8, color=MUT, ha="left", va="bottom")
    return fig, ax

def save(fig, name):
    fig.savefig(f"{OUT}/{name}", facecolor="white")
    plt.close(fig)
    print(name)

def title(ax, main, sub=None):
    ax.set_title(main, fontsize=14.5, color=INK, pad=30 if sub else 12,
                 fontweight="bold", loc="left")
    if sub:
        ax.text(0, 1.018, sub, transform=ax.transAxes, fontsize=10.5, color=SEC, va="bottom")

# ---- P2 frame a: terrain only ----
fig, ax = canvas()
title(ax, "Terrain is just numbers on a grid")
save(fig, "p2_a_terrain.png")

# ---- P2 frame b: virtual rain ----
fig, ax = canvas()
title(ax, "Drop one unit of rain on every cell")
pts = rng.integers(0, [NR, NC], size=(1400, 2))
ax.scatter(pts[:, 1], pts[:, 0], s=5, c=BLUE, alpha=0.75, lw=0)
save(fig, "p2_b_rain.png")

# ---- P2 frame c: accumulation ----
fig, ax = canvas()
title(ax, "Count the drops passing through each cell", "flow accumulation: the valleys emerge on their own")
la = np.log10(acc)
blues = LinearSegmentedColormap.from_list("b", ["#9ec5f4", "#2a78d6", "#0d366b"])
show = np.ma.masked_where(la < np.log10(40), la)   # hide hillslopes, light the valleys
ax.imshow(show, cmap=blues, interpolation="nearest", vmin=np.log10(40), vmax=la.max())
save(fig, "p2_c_acc.png")

# ---- P3: threshold sweep (sparse -> dense) ----
ac_list = [2.0, 1.0, 0.5, 0.25, 0.125, 0.06, 0.03, 0.015]  # km^2
for n, ac_km2 in enumerate(ac_list):
    thr_cells = ac_km2 * 1e6 / (CELL * CELL)
    stream = acc >= thr_cells
    dd = stream.sum() * CELL / 1000 / area_km2
    fig, ax = canvas()
    title(ax, f"A stream begins where enough land drains into it",
          f"support area  A$_c$ = {ac_km2:g} km²      drainage density {dd:.2f} km/km²")
    m = np.ma.masked_where(~stream, np.log10(acc))
    ax.imshow(m, cmap=blues, interpolation="nearest")
    save(fig, f"p3_sweep_{n}.png")

# ---- P4 topology: the outlet's watershed + its directed network ----
thr_cells = 0.125 * 1e6 / (CELL * CELL)
stream = acc >= thr_cells
out_flat = np.argmax(np.where(recv.ravel() < 0, acc.ravel(), 0))
oi, oj = divmod(out_flat, NC)

# terminal outlet of every cell (receivers are always lower -> ascending order works)
term = np.full(NR * NC, -1, int)
asc = np.argsort(filled, axis=None)
for idx in asc:
    r = rr[idx] if False else recv.ravel()[idx]
    term[idx] = idx if r < 0 else term[r]
basin = (term.reshape(NR, NC) == out_flat)

fig, ax = canvas()
title(ax, "Then make it a network: directed, ordered, one outlet",
      "the watershed = every cell whose water reaches this outlet")
ax.imshow(np.ma.masked_where(basin, np.ones_like(acc)), cmap="gray",
          vmin=0, vmax=1, alpha=0.62, interpolation="nearest")   # dim outside basin
ax.contour(basin.astype(float), levels=[0.5], colors="#b84a1e", linewidths=1.6)
m = np.ma.masked_where(~(stream & basin), np.log10(acc))
ax.imshow(m, cmap=blues, interpolation="nearest")
ax.plot(oj, oi, marker="*", ms=23, color="#eb6834", mec="white", mew=1.5, zorder=6)
ax.annotate("outlet", xy=(oj, oi), xytext=(oj + 16, oi + 26), fontsize=13,
            color="#b84a1e", fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color="#b84a1e", lw=1.4))
save(fig, "p4_topology.png")

# ---------------- P1: tiny-grid D8 explainer (3 frames) ----------------
def tiny_canvas():
    fig, ax = plt.subplots(figsize=(6.4, 5.9), dpi=170)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.875, bottom=0.03)
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    ax.set_xlim(-0.55, 4.55); ax.set_ylim(4.55, -0.55)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    return fig, ax

g = np.array([[86, 78, 72, 70, 74],
              [80, 69, 62, 61, 66],
              [74, 64, 55, 52, 58],
              [70, 60, 49, 43, 50],
              [67, 57, 46, 38, 44]], float)

def draw_grid(ax, highlight=None):
    for i in range(5):
        for j in range(5):
            hl = highlight == (i, j)
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, lw=1.2,
                          ec="#c3c2b7", fc="#fff8e1" if hl else "#f7f7f5", zorder=1))
            ax.text(j, i - 0.18, f"{g[i,j]:.0f}", ha="center", va="center",
                    fontsize=15, color=INK, fontweight="bold", zorder=3)
            ax.text(j, i + 0.24, "m", ha="center", va="center", fontsize=8.5,
                    color=MUT, zorder=3)

def steepest(i, j):
    best, bk = 0.0, None
    for k, (di, dj) in enumerate(NB):
        ni, nj = i + di, j + dj
        if 0 <= ni < 5 and 0 <= nj < 5:
            s = (g[i, j] - g[ni, nj]) / dist[k]
            if s > best:
                best, bk = s, (di, dj)
    return bk

fig, ax = tiny_canvas()
ax.set_title("A DEM: elevation on a grid", fontsize=13, color=INK, pad=12, fontweight="bold", loc="left")
draw_grid(ax)
save(fig, "p1_a_grid.png")

fig, ax = tiny_canvas()
ax.set_title("Each cell asks its 8 neighbors: who is steepest downhill?",
             fontsize=12, color=INK, pad=12, fontweight="bold", loc="left")
ci, cj = 2, 2
draw_grid(ax, highlight=(ci, cj))
bk = steepest(ci, cj)
for di, dj in NB:
    tgt = (ci + di, cj + dj)
    if 0 <= tgt[0] < 5 and 0 <= tgt[1] < 5:
        is_best = (di, dj) == bk
        ax.annotate("", xy=(cj + dj * 0.72, ci + di * 0.72), xytext=(cj + dj * 0.3, ci + di * 0.3),
                    arrowprops=dict(arrowstyle="-|>", lw=3 if is_best else 1.2,
                                    color=BLUE_D if is_best else "#c3c2b7"), zorder=4)
save(fig, "p1_b_choice.png")

fig, ax = tiny_canvas()
ax.set_title("Water flows downhill: every cell points to one neighbor", fontsize=12,
             color=INK, pad=12, fontweight="bold", loc="left")
draw_grid(ax)
for i in range(5):
    for j in range(5):
        bk = steepest(i, j)
        if bk:
            di, dj = bk
            ax.annotate("", xy=(j + dj * 0.66, i + di * 0.66), xytext=(j + dj * 0.22, i + di * 0.22),
                        arrowprops=dict(arrowstyle="-|>", lw=2.2, color=BLUE), zorder=4)
save(fig, "p1_c_arrows.png")

print(f"domain {NR}x{NC} cells at {CELL:.0f} m = {area_km2:.1f} km2; max acc {acc.max():.0f} cells")
