#!/usr/bin/env python3
"""
phi_unwrap_visual.py
~~~~~~~~~~~~~~~~~~~~

Slide-quality visualisation of the EMCal **φ-seam un-wrapping** procedure.

The EMCal barrel has 256 fine-φ towers; tower 255 is next to tower 0.
A naïve arithmetic mean across the seam is therefore wrong by ~½ barrel.
The *un-wrapping* algorithm

  1. moves every hit into a continuous window centred on the seed tower,
  2. takes the energy-weighted centroid in that window, then
  3. refolds the answer back to the detector range 0…255.

This script shows the process on a toy 4-tower “cluster” that straddles
the seam.

Quick start
-----------
$ ./phi_unwrap_visual.py                      # writes ~/Desktop/phi_unwrap_demo.png
$ ./phi_unwrap_visual.py --png /tmp/demo.png  # custom output
"""
from __future__ import annotations

# ── std / third-party ──────────────────────────────────────────────
from pathlib import Path
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle

# ── EMCal constants ────────────────────────────────────────────────
NPHI    = 256            # number of fine-φ towers
PATCH_R = 2              # half-width of the zoom-in patch (towers)
FIGSIZE = (12, 5.0)      # in inches
DPI     = 150
CMAP    = "viridis"

# --- make "value 0" show up in grey instead of viridis-purple ------------
vir      = mpl.colormaps[CMAP].resampled(256)  # works with all 3.5+ releases
newcols  = vir(np.linspace(0, 1, 256))
newcols[0] = np.array([0.85, 0.85, 0.85, 1.0]) # light opaque grey RGBA
CMAP_GRID = mpl.colors.ListedColormap(newcols) # use this in imshow()

OUT_DEF = Path.home() / "Desktop" / "phi_unwrap_demo.png"


mpl.rcParams.update({
    "savefig.dpi":   DPI,
    "font.size":     12,
    "axes.labelsize":12,
    "axes.titlesize":13,
    "xtick.labelsize":11,
    "ytick.labelsize":11,
    "legend.fontsize":11,
    "axes.linewidth":1.3,
})

# ╔════════════════════════════════════════════════════════════════╗
# ║  Seam-un-wrapping helper                                      ║
# ╚════════════════════════════════════════════════════════════════╝
def phi_centroid_unwrapped(
    phi_indices: np.ndarray,
    energies:    np.ndarray
) -> tuple[np.ndarray, float]:
    """
    Shift all φ indices into a seed-centred ±½-barrel window, compute the
    energy-weighted centroid there, then fold it back to 0…255.
    """
    seed = int(phi_indices[0])

    phi_unw = phi_indices.astype(float).copy()
    dphi    = phi_unw - seed
    phi_unw[dphi < -NPHI/2] += NPHI
    phi_unw[dphi >  NPHI/2] -= NPHI

    phi_cent_unw = float(np.dot(phi_unw, energies) / energies.sum())
    phi_cent_det = (phi_cent_unw + NPHI) % NPHI     # detector range

    return phi_unw, phi_cent_det


# ═══════════════════════════════════════════════════════════════════
#  utilities for the figure
# ═══════════════════════════════════════════════════════════════════
def build_patch(
    eta:        np.ndarray,   # tower η indices
    phi_det:    np.ndarray,   # tower φ indices *in detector units*
    energy:     np.ndarray,
    centre_eta: float,
    centre_phi: float
):
    """
    Return a PATCH×PATCH sub-matrix around the centroid.
    The φ values are first made continuous (no 256 → 0 jump) so that
    imshow draws a proper rectangle grid.
    """
    # ── continuous φ window centred on the seed tower ───────────────
    seed_phi = int(phi_det[0])
    phi_cont = phi_det.astype(float).copy()

    dphi = phi_cont - seed_phi
    phi_cont[dphi < -NPHI/2] += NPHI   # e.g. −2  → 254
    phi_cont[dphi >  NPHI/2] -= NPHI   # e.g. 257 →   1

    # ── choose integer window ───────────────────────────────────────
    eta0  = int(round(centre_eta))
    phi0  = int(round(centre_phi))
    rows  = np.arange(eta0 - PATCH_R, eta0 + PATCH_R + 1)
    cols  = np.arange(phi0 - PATCH_R, phi0 + PATCH_R + 1)
    patch = np.zeros((rows.size, cols.size))

    for η, φ, w in zip(eta, phi_cont, energy):
        if rows[0] <= η <= rows[-1] and cols[0] <= φ <= cols[-1]:
            patch[int(η - rows[0]), int(φ - cols[0])] = w

    seed_r = int(eta[0]      - rows[0])
    seed_c = int(phi_cont[0] - cols[0])

    return patch, rows, cols, seed_r, seed_c


def annotate_patch(
    ax, *, patch, rows, cols,
    seed_r, seed_c,
    cog_eta, cog_phi,
    title: str,
    seam_band: bool = False,
    centroid_fmt: str | None = None,
    colour: str = "black",
):
    """Draw coloured squares, optional seam band and centroid marker."""
    im = ax.imshow(
        patch,
        origin="lower",
        cmap=CMAP_GRID,                      # <--- changed
        extent=[cols[0] - .5, cols[-1] + .5, rows[0] - .5, rows[-1] + .5],
        vmin=0,
        vmax=patch.max()
    )

    if seam_band:
        ax.axvspan(-0.5, 0.5, facecolor="lightgrey", alpha=.35, zorder=0)
        ax.axvline( 0.0, color="dimgrey", lw=1.6, ls="--", zorder=1)
        ax.text(0, rows[0]-0.8, "φ seam",
                ha="center", va="top", fontsize=11, color="dimgrey")

    ax.add_patch(
        Rectangle((cols[seed_c]-.5, rows[seed_r]-.5), 1, 1,
                  edgecolor="white", facecolor="none", lw=2.0, zorder=2)
    )

    for (r, c), val in np.ndenumerate(patch):
        if val == 0:        # skip empty squares
            continue
        ax.text(
            c+cols[0], r+rows[0], f"{val:.0f}",
            ha="center", va="center", weight="bold",
            color="black", fontsize=11,
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            zorder=3
        )

    marker = dict(marker="x", ms=13, mew=2.2, color=colour) \
             if colour == "black" else dict(marker="*", ms=13, color=colour)
    ax.plot(cog_phi, cog_eta, **marker, zorder=4)

    if centroid_fmt:
        ax.text(cog_phi+0.35*np.sign(marker["ms"]), cog_eta+0.55,
                centroid_fmt, ha="left" if colour!="black" else "right",
                va="bottom", color=colour, fontsize=10)

    ax.set_title(title, pad=8)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(cols[0], cols[-1]+1))
    ax.set_yticks(np.arange(rows[0], rows[-1]+1))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.grid(color="k", lw=0.7, ls=":", alpha=0.25)   # darker, crisper grid
    return im


# ═══════════════════════════════════════════════════════════════════
#  main visualisation
# ═══════════════════════════════════════════════════════════════════
def visualise(cluster: np.ndarray, out_png: Path | str = OUT_DEF) -> None:
    """
    Build the 2-panel slide that explains φ seam un-wrapping.
    """
    out_png = Path(out_png).expanduser().resolve()
    out_png.parent.mkdir(parents=True, exist_ok=True)

    # unpack cluster rows = (η, φ, E)
    eta, phi_raw, E = cluster.T
    phi_raw = phi_raw.astype(int)

    # raw φ in a seam-centred frame (−128 … +127)
    phi_raw_rel = ((phi_raw + NPHI//2) % NPHI) - NPHI//2

    # energy-weighted centroids in matching frames
    cog_eta     = float(np.dot(eta,        E) / E.sum())
    cog_phi_raw = float(np.dot(phi_raw_rel, E) / E.sum())

    # unwrapped indices and centroid
    phi_unw, phi_cent = phi_centroid_unwrapped(phi_raw, E)
    cog_phi_unw_rel   = float(np.dot(phi_unw - phi_unw[0], E) / E.sum())

    # patches
    patch_raw, rows, cols_raw, sr_raw, sc_raw = build_patch(
        eta, phi_raw_rel, E, cog_eta, cog_phi_raw
    )
    patch_unw, _, cols_unw, sr_unw, sc_unw = build_patch(
        eta, phi_unw - phi_unw[0], E, cog_eta, cog_phi_unw_rel
    )

    # canvas
    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=FIGSIZE, sharey=True, constrained_layout=True
    )

    # left panel – raw
    annotate_patch(
        axL, patch=patch_raw, rows=rows, cols=cols_raw,
        seed_r=sr_raw, seed_c=sc_raw,
        cog_eta=cog_eta, cog_phi=cog_phi_raw,
        title="Raw tower indices\n(discontinuous across φ seam)",
        seam_band=True,
        colour="black"
    )

    # ---- SINGLE x-axis for RAW panel ----
    axL.set_xticks(cols_raw)
    axL.set_xticklabels(cols_raw.astype(int))
    axL.set_xlabel("raw φ index")

    # right panel – unwrapped
    imR = annotate_patch(
        axR, patch=patch_unw, rows=rows, cols=cols_unw,
        seed_r=sr_unw, seed_c=sc_unw,
        cog_eta=cog_eta, cog_phi=cog_phi_unw_rel,
        title="After seam un-wrapping\n(continuous local φ)",
        seam_band=False,
        colour="red"
    )

    # ---- SINGLE x-axis for UNWRAPPED panel ----
    axR.set_xticks(cols_unw)
    axR.set_xticklabels(cols_unw.astype(int))
    axR.set_xlabel("unwrapped φ index (seed-centred)")

    # curved arrow
    axL.annotate(
        "unwrap →",
        xy=(cog_phi_unw_rel, cog_eta+.15),
        xytext=(cog_phi_raw-.5, cog_eta+.8),
        arrowprops=dict(arrowstyle="-|>",
                        connectionstyle="arc3,rad=-0.30",
                        lw=2.2, color="black"),
        ha="center", va="bottom", fontsize=12
    )

    # colour-bar
    cax = fig.add_axes([0.46, 0.18, 0.02, 0.64])   # [L, B, W, H]
    cb  = fig.colorbar(imR, cax=cax)
    cb.set_label("Tower energy (arb.)", labelpad=6)
    cb.ax.yaxis.set_offset_position('right')
    cb.ax.tick_params(length=0)

    # ---------- caption: explicit weighted-average formulas (TOP) ----------
    def _num_str(e_arr, phi_arr):
        """Return '8×(-2) + 10×(-1) + …' ready for math-text."""
        return r"\;+\;".join(fr"{int(e)}\times({int(p)})"
                             for e, p in zip(e_arr, phi_arr))

    Ei_sum = E.sum()

    # raw (discontinuous) ---------------------------------------------------
    num_raw   = _num_str(E, phi_raw_rel)
    phi_raw_w = np.dot(E, phi_raw_rel) / Ei_sum

    # unwrapped (continuous) -----------------------------------------------
    local_unw = phi_unw - phi_unw[0]          # 0,1,2,3 …
    num_unw   = _num_str(E, local_unw)
    phi_unw_w = np.dot(E, local_unw) / Ei_sum


    # helper – two lines, "=" perfectly aligned -----------------------------
    def _two_line_eq(left, frac_str, value, x, y):
        """
          line 1: <left>=<fraction>
          line 2:        =<value>
        """
        # ── first line ───────────────────────────────────────────────────
        line1 = rf"${left}=\dfrac{{{frac_str}}}{{{int(Ei_sum)}}}$"
        fig.text(x, y, line1, ha="left", va="top",
                 fontsize=6,
                 path_effects=[pe.withStroke(linewidth=3, foreground='white')])

        # ── measure horizontal position of the " = " column ─────────────
        ghost = fig.text(x, y, rf"${left}$", ha="left", va="top",
                         fontsize=6, alpha=0.0)          # invisible probe
        fig.canvas.draw()                                # flush layout
        renderer = fig.canvas.get_renderer()
        w_pix = ghost.get_window_extent(renderer=renderer).width
        ghost.remove()

        dx = w_pix / (fig.get_figwidth() * fig.dpi)      # → figure-fraction

        # ── second line (numeric value) ─────────────────────────────────
        eq_val = rf"$={value:.2f}$"                      # starts with "="
        fig.text(x + dx,            # shift by width of <left>
                 y - 0.048,         # vertical gap
                 eq_val,
                 ha="left", va="top", fontsize=6,
                 path_effects=[pe.withStroke(linewidth=3,
                                             foreground='white')])


    # left (raw) – slight right shift so it sits inside grey box
    _two_line_eq(r"\phi_{\text{weighted, naïve}}",
                 num_raw,  phi_raw_w,
                 0.10, 0.89)

    # right (unwrapped)
    _two_line_eq(r"\phi_{\text{weighted, unwrapped}}",
                 num_unw,  phi_unw_w,
                 0.60, 0.89)


    # save
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"✅  Saved illustration →  {out_png}")


# ═══════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════
def _parse_args():
    ap = argparse.ArgumentParser(description="Visualise EMCal φ-seam un-wrapping.")
    ap.add_argument("--png", default=str(OUT_DEF), metavar="FILE",
                    help=f"Output PNG file (default: {OUT_DEF})")
    return ap.parse_args()


def main():
    # toy cluster that straddles the φ = 0 seam
    toy_cluster = np.array([
        (40, 254,  8.0),   # η, φ, E
        (40, 255, 10.0),
        (41,   0, 12.0),
        (41,   1,  6.0),
    ])

    args = _parse_args()
    visualise(toy_cluster, args.png)


if __name__ == "__main__":
    main()
