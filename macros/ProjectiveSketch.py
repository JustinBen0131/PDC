#!/usr/bin/env python3
# ProjectiveSketch.py
#
# Remove ALL visible gaps between pairs by solving the exact “domino seam”
# analytically. We anchor each pair at the bottom‑right (BR) corner of its
# right rectangle on a common baseline. For consecutive pairs j → j+1, we
# compute the horizontal shift so that the entire left edge of pair j+1
# lies *on or slightly left of* the entire right edge of pair j across the
# full vertical overlap. This guarantees touching (with a tiny controlled
# overlap) even when the tilts differ.

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ---------------------------- Canvas -----------------------------------
WIDTH, HEIGHT = 1039, 123
DPI = 200
MARGIN = 10

# ----------------------- Pair & Slat Settings --------------------------
N_PAIRS      = 24
FIRST_VERT   = 3

RECT_W       = 32.0   # identical rectangle width (unscaled units)
RECT_H       = 130.0  # identical rectangle height (unscaled units)

MIN_TILT_DEG = 1.0     # slightly smaller initial tilt
MAX_TILT_DEG = 46.0    # keep final attitude
EASE_POW     = 2.8     # gentler ramp so step-to-step tilt looks less abrupt

EDGE_COLOR     = (0.35, 0.35, 0.35)
LW_BASE        = 1.15   # thinner strokes to match the reference
LW_LEFT_EDGE   = 1.8    # still slightly emphasized
RAIL_LW        = 2.6
RAIL_VMARGIN   = 6.0

# -------------------------- Geometry Utils -----------------------------
def rot_ccw(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s],
                     [s,  c]])

def build_pair_from_BR(anchor_br, tilt_deg_cw_from_vertical):
    """
    Build one block pair given BR of the RIGHT rectangle.
    Returns rectA(left), rectB(right) as [TL, TR, BR, BL] in y-up.
    """
    theta = -np.deg2rad(tilt_deg_cw_from_vertical)  # CW from vertical
    R = rot_ccw(theta)
    u = R @ np.array([1.0, 0.0])   # width direction
    v = R @ np.array([0.0, 1.0])   # height direction (up)

    # Right rectangle (B), anchored at BR
    BR_B = anchor_br
    BL_B = BR_B - RECT_W * u
    TR_B = BR_B + RECT_H * v
    TL_B = TR_B - RECT_W * u

    # Left rectangle (A) immediately to the left
    BR_A = BL_B
    BL_A = BR_A - RECT_W * u
    TR_A = BR_A + RECT_H * v
    TL_A = TR_A - RECT_W * u

    rectA = np.vstack([TL_A, TR_A, BR_A, BL_A])
    rectB = np.vstack([TL_B, TR_B, BR_B, BL_B])
    return rectA, rectB

# ----------------------- Tilt Schedule ---------------------------------
tilts = []
for j in range(N_PAIRS):
    if j < FIRST_VERT:
        tilts.append(0.0)
    else:
        t = (j - FIRST_VERT) / max(1, (N_PAIRS - FIRST_VERT - 1))
        tilts.append(MIN_TILT_DEG + (MAX_TILT_DEG - MIN_TILT_DEG) * (t ** EASE_POW))

# ----------------------- Anchor Placement (seam‑tight) ------------------
# Analytical seam condition:
#   Right edge of pair j is a line through BR_j with slope tan(theta_j).
#   Left edge of pair j+1 is a line through BL_A of j+1 with slope tan(theta_{j+1}).
#   To eliminate wedges for the full vertical overlap, we choose x_{j+1} so that
#   for all y in [0, y_overlap], x_left_{j+1}(y) <= x_right_j(y) - delta,
#   where y_overlap = RECT_H * min(cos(theta_j), cos(theta_{j+1})) and
#   delta is a small positive overlap in *world* units.
#
# With BR anchors (x_j at y=0), this yields:
#   dx_touch = 2*RECT_W*cos(theta_{j+1}) - (tan(theta_{j+1}) - tan(theta_j)) * y_overlap
#   x_{j+1}  = x_j + dx_touch - delta

def compute_anchors(tilts_deg, overlap_world):
    """
    Balanced seam placement with BR anchors.
    We enforce equality at the *midpoint* of the common vertical span between
    the two edges so the apparent overlap stays uniform across the row.
    A small negative `overlap_world` still crushes antialias slivers.
    """
    thetas = -np.deg2rad(np.asarray(tilts_deg))  # CCW angles
    anchors = np.zeros(len(thetas), dtype=float)

    for j in range(len(thetas) - 1):
        t0 = thetas[j]
        t1 = thetas[j + 1]

        cos1 = np.cos(t1)
        sin1 = np.sin(t1)
        tan0 = np.tan(t0)
        tan1 = np.tan(t1)

        # Common vertical overlap of right edge (pair j, y in [0, RECT_H])
        # and left edge (pair j+1, y in [y_bl, y_bl + RECT_H*cos1]) in y-up space.
        y_bl   = max(0.0, -2.0 * RECT_W * sin1)
        y_top1 = y_bl + RECT_H * cos1
        y_low  = y_bl
        y_high = min(RECT_H, y_top1)

        # Degenerate overlap: fall back to touching by width only.
        if y_high <= y_low + 1e-9:
            sec1 = 1.0 / max(abs(cos1), 1e-9)
            dx   = 2.0 * RECT_W * sec1 - overlap_world
        else:
            # Balanced contact: use the midpoint of the valid overlap band.
            y_star = 0.5 * (y_low + y_high)
            sec1   = 1.0 / max(abs(cos1), 1e-9)
            dx     = 2.0 * RECT_W * sec1 + (tan1 - tan0) * y_star - overlap_world

        anchors[j + 1] = anchors[j] + dx

    return anchors



# Pass 1 (delta = 0) → estimate scale so we can express pixel overlap in world units
anchors_pre = compute_anchors(tilts, overlap_world=0.0)
BOTTOM_Y = 0.0
pairs_yup_pre = [build_pair_from_BR(np.array([ax, BOTTOM_Y]), a)
                 for ax, a in zip(anchors_pre, tilts)]
pre_all = np.vstack([np.vstack([rA, rB]) for rA, rB in pairs_yup_pre])
min_x, min_y = pre_all[:, 0].min(), pre_all[:, 1].min()
max_x, max_y = pre_all[:, 0].max(), pre_all[:, 1].max()
span_x, span_y = max_x - min_x, max_y - min_y
scale = min((WIDTH - 2*MARGIN) / span_x, (HEIGHT - 2*MARGIN) / span_y)

# Choose a seam overlap in *pixels*, then convert to world units.
# We use half a stroke plus a small constant to crush antialias slivers everywhere.
HALF_STROKE_PX = (LW_BASE * DPI / 72.0) * 0.5
BASE_EPS_PX    = 0.6
OVERLAP_WORLD  = (HALF_STROKE_PX + BASE_EPS_PX) / scale

# Pass 2 (with overlap) → final anchors and geometry
anchors_x = compute_anchors(tilts, overlap_world=OVERLAP_WORLD)
pairs_yup = [build_pair_from_BR(np.array([ax, BOTTOM_Y]), a)
             for ax, a in zip(anchors_x, tilts)]

# -------------------------- Fit To Canvas ------------------------------
all_pts = np.vstack([np.vstack([rA, rB]) for rA, rB in pairs_yup])
min_x, min_y = all_pts[:, 0].min(), all_pts[:, 1].min()
max_x, max_y = all_pts[:, 0].max(), all_pts[:, 1].max()
span_x, span_y = max_x - min_x, max_y - min_y

scale = min((WIDTH - 2*MARGIN) / span_x, (HEIGHT - 2*MARGIN) / span_y)

def to_screen(p_yup):
    x = (p_yup[0] - min_x) * scale + MARGIN
    y = (max_y - p_yup[1]) * scale + MARGIN  # y-down
    return np.array([x, y])

pairs_px = []
for rA, rB in pairs_yup:
    pairs_px.append((
        np.vstack([to_screen(pt) for pt in rA]),
        np.vstack([to_screen(pt) for pt in rB]),
    ))

# Rails vertical span (for the left rail only)
all_px = np.vstack([np.vstack([rA, rB]) for rA, rB in pairs_px])
ymin_px, ymax_px = all_px[:, 1].min() - RAIL_VMARGIN, all_px[:, 1].max() + RAIL_VMARGIN
left_edge_x  = pairs_px[0][0][0, 0]  # TL of first pair’s left rectangle

# ------------------------------ Draw -----------------------------------
fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI)
ax = plt.axes([0, 0, 1, 1])
ax.set_xlim(0, WIDTH)
ax.set_ylim(HEIGHT, 0)
ax.axis("off")

# Background
bg = Path([(0, 0), (WIDTH, 0), (WIDTH, HEIGHT), (0, HEIGHT), (0, 0)],
          [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
ax.add_patch(PathPatch(bg, facecolor="white", edgecolor="none"))

# Draw RIGHT → LEFT so the LEFT (earlier) pair ends up on top wherever there is overlap
for idx in reversed(range(len(pairs_px))):
    rA, rB = pairs_px[idx]
    for rect in (rA, rB):
        verts = np.vstack([rect, rect[0]])
        codes = [Path.MOVETO] + [Path.LINETO] * 4
        ax.add_patch(PathPatch(
            Path(verts, codes),
            facecolor="none",
            edgecolor=EDGE_COLOR,
            linewidth=LW_BASE,
            capstyle="round",
            joinstyle="round",
            antialiased=False,   # keep crisp joins; overlap removes slivers
        ))

# Draw left edges from RIGHT → LEFT and give them a very high zorder so they sit on top
for idx in reversed(range(len(pairs_px))):
    rA, _ = pairs_px[idx]
    TL, BL = rA[0], rA[3]
    ax.plot([TL[0], BL[0]], [TL[1], BL[1]],
            color=EDGE_COLOR, linewidth=LW_LEFT_EDGE, solid_capstyle="round", zorder=1000)


# Left rail + stubs
ax.plot([left_edge_x, left_edge_x], [ymin_px, ymax_px],
        color=EDGE_COLOR, linewidth=RAIL_LW, solid_capstyle="round")
first_rA = pairs_px[0][0]
TL0, BL0 = first_rA[0], first_rA[3]
ax.plot([left_edge_x, TL0[0]], [TL0[1], TL0[1]], color=EDGE_COLOR, linewidth=LW_BASE)
ax.plot([left_edge_x, BL0[0]], [BL0[1], BL0[1]], color=EDGE_COLOR, linewidth=LW_BASE)

# Save
out_path = "/Users/patsfan753/Desktop/ProjectiveSketch.png"
os.makedirs(os.path.dirname(out_path), exist_ok=True)
fig.savefig(out_path, dpi=DPI, bbox_inches="tight", pad_inches=0)
