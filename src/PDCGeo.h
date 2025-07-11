/***************************************************************************************************
 *  PDCGeo.h  –  header‑only, zero‑overhead geometry helpers for the sPHENIX EMCal
 *
 *  Rationale
 *  =========
 *  Precise photon reconstruction in the sPHENIX EMCal requires
 *  corrections that are *local* to a 2×2 “super‑tower” block and *global* in η/φ.  This file
 *  encapsulates the purely geometrical parts of that task:
 *
 *    1.  **Symmetric folding helpers** (`foldOnce`, `foldAndStep`) that map any local coordinate
 *        into the canonical range (‑0.5 … +1.5] while keeping the *coarse* block index consistent.
 *
 *    2.  **Ash distortion** tools (`phi::undoAsh`, `eta::undoAsh`) – the analytical inverse of the
 *        “Ash‑b” radial electric‑field distortion in drift chambers.  For φ we additionally apply a
 *        small edge fix (needed right at the barrel seams), whereas for η the plain inverse is
 *        sufficient.
 *
 *    3.  **Block‑to‑global converters** that turn a (coarse‑block index, local offset) pair into
 *        global η or φ, optionally corrected by a *run‑time* offset functor (used for the measured
 *        rigid barrel tilt).
 *
 *    4.  **Fine‑index helpers** – map the same (blk,loc) pair to the EMCal’s native 256 × 96 tower
 *        grid, taking care of periodicity in φ and hard clamping in η.
 *
 *    5.  **Local coordinate toolkit** (`computeLocal`) – a *single‑pass* tower scanner that
 *        computes   ΣE, ΣE·η, ΣE·φ   and returns both the block address **and** the folded local
 *        offsets.  This fully replaces legacy getBlockCord + getAvg{Eta,‑Phi}.
 *
 *  Design choices
 *  --------------
 *  •  **Header‑only & constexpr** – all helpers are inlined, so there is zero binary overhead.
 *  •  **No external state** – the code depends solely on its arguments.  Global parameters such as
 *     tower pitch, detector size, … are compile‑time constants defined up‑front.
 *  •  **C++17‑only** – uses nothing beyond language / STL facilities available in C++17 (the sPHENIX
 *     default).  ROOT is needed only for `TVector2::Phi_mpi_pi`.
 *
 *  Public interface (stable)
 *  -------------------------
 *      •  foldOnce(u)
 *      •  phi::undoAsh(u,b)
 *      •  eta::undoAsh(u,b[,pitch])
 *      •  phi::blockToGlobal(blk,loc[,offsetFn])
 *      •  eta::blockToGlobal(blk,loc[,offsetFn])
 *      •  finePhiIdx(blk,loc) / fineEtaIdx(blk,loc)
 *      •  struct LocalCoord  +  computeLocal(twrEta,twrPhi,twrE)
 ***************************************************************************************************/
#ifndef PDC_GEO_H
#define PDC_GEO_H   1

/* ─────────────────────────  Standard / ROOT includes  ───────────────────────── */
#include <cmath>        // fabs, sinh, asinh, copysign, fmod, …
#include <algorithm>    // clamp, max_element
#include <limits>       // quiet_NaN
#include <type_traits>  // std::true_type / std::false_type
#include <vector>
#include <TVector2.h>   // TVector2::Phi_mpi_pi

namespace PDC::Geo
{
/* ═══════════════════════════════════════════════════════════════════════
 *  1.  Global detector constants
 *      – fixed for the sPHENIX EMCal barrel (CEMC).
 * ═════════════════════════════════════════════════════════════════════ */
inline constexpr int   kFinePerBlock   = 2;          ///< #fine towers per coarse block side
inline constexpr int   kFinePhiBins    = 256;        ///< Full EMCal barrel in φ
inline constexpr int   kFineEtaBins    =  96;        ///< Barrel height in η
inline constexpr int   kCoarsePhiBins  = kFinePhiBins / kFinePerBlock;   ///< 128
inline constexpr int   kCoarseEtaBins  = kFineEtaBins / kFinePerBlock;   ///<  48
inline constexpr float kRadPerFine     = 2.F * static_cast<float>(M_PI) / kFinePhiBins; ///< Δφ per fine tower
inline constexpr float kEtaMin         = -1.1F;      ///< Centre of η‑row 0
inline constexpr float kDEtaPerFine    =  2.2F / kFineEtaBins;           ///< Δη per fine tower

/* ═══════════════════════════════════════════════════════════════════════
 *  2.  Low‑level helpers
 * ═════════════════════════════════════════════════════════════════════ */
/** Symmetric fold of a *local* coordinate
 *  Maps any value into (‑0.5 … +1.5], preserving relative ordering.
 *  This is the canonical range used throughout the clusteriser (legacy),
 *  hence we keep it as the public contract. */
[[nodiscard]] inline constexpr float
foldOnce(float u) noexcept
{
    if (u <= -0.5F || u > 1.5F)
    {
        u = std::fmod(u + 2.F, 2.F);        // bring into 0 … <2
        if (u > 1.5F) u -= 2.F;             // final range  (‑0.5 , +1.5]
    }
    return u;
}

/** Fold + coarse‑index bookkeeping
 *  *If* the coordinate jumps across a 2‑tower border we increment (or
 *  decrement) the coarse index so that (blk,loc) stays a *bona‑fide*
 *  representation of the same physical position. */
template<bool Periodic>
inline constexpr void
foldAndStep(float& loc, int& coarse, int nCoarse)
{
    if (loc > -0.5F && loc < 1.5F) return;      // nothing to do

    const bool crossedLeft = (loc < -0.5F);     // left ↔ decreasing φ

    loc = std::fmod(loc + 2.F, 2.F);
    if (loc >= 1.5F) loc -= 2.F;

    coarse += crossedLeft ? -1 : +1;

    /* φ wraps around, η clamps at physical ends */
    if constexpr (Periodic)
    {
        if (coarse < 0)               coarse += nCoarse;
        else if (coarse >= nCoarse)   coarse -= nCoarse;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  3.  Ash‑distortion inverse
 *
 *  undoAshGeneric  →   measured  ↦  true
 * ═════════════════════════════════════════════════════════════════════ */
namespace _impl
{
    /* Small φ edge correction to smoothly glue the two halves of the barrel
     * together.  Tuned on simulation, negligible in the central region. */
    struct PhiEdgeFix
    {
        [[nodiscard]] static constexpr float apply(float corr, float b) noexcept
        {
            constexpr float Aref = 0.012F, bref = 0.18F;
            float A = std::clamp(Aref * (b / bref) * (b / bref), 0.008F, 0.016F);
            const float s = static_cast<float>(M_PI) * (corr - 0.5F);
            return corr - A * std::sin(s) * (1.F - 0.25F * std::cos(s));
        }
    };
    /* No edge correction needed for η */
    struct NoEdgeFix
    {
        [[nodiscard]] static constexpr float apply(float corr, float) noexcept
        { return corr; }
    };

    /* Generic inverse – edge policy selected at compile time */
    template<class EdgePolicy>
    [[nodiscard]] inline constexpr float
    undoAshGeneric(float loc, float b,
                   float towerPitch = 1.F)              // φ: pitch=1, η: pitch≈Δη
    {
        if (std::fabs(b) < 1e-9F) return loc * towerPitch;   // b≈0 → identity
        loc = foldOnce(loc);

        const bool  right = (loc > 0.5F);                   // right half‑tower
        const float xMeas = right ? loc - 1.F : loc;        // map to (‑0.5 … +0.5]

        const double S    = std::sinh(1. / (2. * b));
        float xTrue       = static_cast<float>(b * std::asinh(2. * S * xMeas));

        float corr = (right ? xTrue + 1.F : xTrue) * towerPitch;
        corr       = EdgePolicy::apply(corr, b);

        /* Guard: keep within barrel (rarely needed) */
        if (std::fabs(corr) > 1.6F)
            corr = std::copysign(1.6F, corr);
        return corr;
    }
} // namespace _impl

/* Public wrappers – nothing fancy */
namespace phi
{
    [[nodiscard]] inline constexpr float
    undoAsh(float loc, float b) noexcept
    { return _impl::undoAshGeneric<_impl::PhiEdgeFix>(loc, b); }
}

namespace eta
{
    [[nodiscard]] inline constexpr float
    undoAsh(float loc, float b, float towerPitch = 1.F) noexcept
    { return _impl::undoAshGeneric<_impl::NoEdgeFix>(loc, b, towerPitch); }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  5.  (blk,loc) → global helpers
 *
 *  Converts a *coarse* block index plus *local* offset into a detector‑
 *  fixed global coordinate.
 * ═════════════════════════════════════════════════════════════════════ */
namespace _impl
{
    template<bool IsPhi, class OffsetProvider>
    [[nodiscard]] inline constexpr float
    blockToGlobalGeneric(int blk, float loc,
                         const OffsetProvider& off) noexcept
    {
        loc = foldOnce(loc);

        if constexpr (IsPhi)
        {
            if (blk < 0 || blk >= kCoarsePhiBins || !std::isfinite(loc))
                return std::numeric_limits<float>::quiet_NaN();

            const float fine = static_cast<float>(blk) * kFinePerBlock
                               + loc + 0.5F;
            return TVector2::Phi_mpi_pi(fine * kRadPerFine + off());
        }
        else    // η
        {
            if (blk < 0 || blk >= kCoarseEtaBins || !std::isfinite(loc))
                return std::numeric_limits<float>::quiet_NaN();

            const float fine = static_cast<float>(blk) * kFinePerBlock
                               + loc + 0.5F;
            return kEtaMin + fine * kDEtaPerFine + off();
        }
    }
} // namespace _impl

/* Thin, user‑friendly wrappers */
namespace phi
{
    template<class OffsetProvider = decltype([](){return 0.F;})>
    [[nodiscard]] inline constexpr float
    blockToGlobal(int blk, float loc,
                  const OffsetProvider& off = [](){return 0.F;}) noexcept
    { return _impl::blockToGlobalGeneric<true>(blk, loc, off); }
}

namespace eta
{
    template<class OffsetProvider = decltype([](){return 0.F;})>
    [[nodiscard]] inline constexpr float
    blockToGlobal(int blk, float loc,
                  const OffsetProvider& off = [](){return 0.F;}) noexcept
    { return _impl::blockToGlobalGeneric<false>(blk, loc, off); }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  6.  Fine‑index helpers
 *
 *  Return the *tower‑grid* index corresponding to (blk,loc).  φ wraps,
 *  η clamps at the physical ends of the barrel.
 * ═════════════════════════════════════════════════════════════════════ */
template<int Nfine, bool Periodic>
[[nodiscard]] inline constexpr int
fineIdx(int blk, float loc) noexcept
{
    loc = foldOnce(loc);
    int idx = static_cast<int>(std::floor(blk * kFinePerBlock + loc + 0.5F));
    if constexpr (Periodic)
    {
        if (idx < 0)            idx += Nfine;
        else if (idx >= Nfine)  idx -= Nfine;
        return idx;
    }
    else
    {
        return std::clamp(idx, 0, Nfine - 1);
    }
}
[[nodiscard]] inline constexpr int
finePhiIdx(int blk, float loc) noexcept
{ return fineIdx<kFinePhiBins, true >(blk, loc); }

[[nodiscard]] inline constexpr int
fineEtaIdx(int blk, float loc) noexcept
{ return fineIdx<kFineEtaBins, false>(blk, loc); }

/* ═══════════════════════════════════════════════════════════════════════
 *  7.  Local 2×2‑block coordinate toolkit
 * ═════════════════════════════════════════════════════════════════════ */
struct LocalCoord
{
    int   blkEta  = 0;          ///< coarse η index  (0 … 47)
    int   blkPhi  = 0;          ///< coarse φ index  (0 … 127)
    float locEta  = 0.F;        ///< local η in (‑0.5 … +1.5]
    float locPhi  = 0.F;        ///< local φ in (‑0.5 … +1.5]
};

namespace _impl
{
    /* Internal helper – accumulate ΣE, ΣE·η, reference φ‑tower, … */
    struct ScanOut
    {
        double sumE      = 0.0;
        double sumEeta   = 0.0;
        int    refPhiBin = 0;        ///< fine‑bin with the highest E
        double refE      = 0.0;
    };

    template<class IntVec, class FloatVec>
    [[nodiscard]] inline constexpr ScanOut
    scanTowers(const IntVec& eta, const IntVec& phi, const FloatVec& E) noexcept
    {
        ScanOut out{};
        const std::size_t nT = E.size();
        for (std::size_t i = 0; i < nT; ++i)
        {
            const double Ei = E[i];
            out.sumE    += Ei;
            out.sumEeta += Ei * (eta[i] + 0.5);

            if (Ei > out.refE) { out.refE = Ei; out.refPhiBin = phi[i]; }
        }
        return out;
    }
} // namespace _impl

/** Compute the *canonical* local coordinate of a tower cluster. */
template<class IntVec, class FloatVec>
[[nodiscard]] inline LocalCoord
computeLocal(const IntVec& towerEta,
             const IntVec& towerPhi,
             const FloatVec& towerE) noexcept
{
    /* ── 0) Quick validation ───────────────────────────────────────── */
    const std::size_t nT = towerE.size();
    if (nT == 0 || towerEta.size() != nT || towerPhi.size() != nT)
        return {};                                // invalid → zero‑init

    /* ── 1) ΣE, ΣE·η and ref‑φ tower (single pass) ────────────────── */
    const auto scan = _impl::scanTowers(towerEta, towerPhi, towerE);

    if (scan.sumE < 1e-12) return {};            // all towers zero

    const float etaFine = static_cast<float>(scan.sumEeta / scan.sumE);  // 0 … 96

    /* ── 2) Unwrap φ around the reference tower and accumulate ΣE·φ ─ */
    constexpr int   Nphi = kFinePhiBins;
    constexpr float Half = Nphi / 2.F;

    double sumEphi = 0.0;
    for (std::size_t i = 0; i < nT; ++i)
    {
        int bin = towerPhi[i];
        int d   = bin - scan.refPhiBin;
        if (d < -Half) bin += Nphi;              // unwrap left
        if (d >  Half) bin -= Nphi;              // unwrap right

        sumEphi += towerE[i] * bin;
    }
    float phiFine = static_cast<float>(sumEphi / scan.sumE);
    phiFine       = std::fmod(phiFine + Nphi, Nphi);      // back to [0,256)

    /* ── 3) Coarse indices + local offsets ─────────────────────────── */
    LocalCoord c;
    c.blkEta = static_cast<int>(std::floor(etaFine)) / kFinePerBlock;
    c.blkPhi = static_cast<int>(std::floor(phiFine)) / kFinePerBlock;
    c.locEta = etaFine - c.blkEta * kFinePerBlock;
    c.locPhi = phiFine - c.blkPhi * kFinePerBlock;

    /* ── 4) Canonical symmetric fold (exactly once per axis) ───────── */
    foldAndStep<false>(c.locEta, c.blkEta, kCoarseEtaBins);
    foldAndStep< true>(c.locPhi, c.blkPhi, kCoarsePhiBins);

    /* Final hard clamp in η (defensive, should *never* trigger) */
    if      (c.blkEta < 0)               { c.blkEta = 0;                   c.locEta = -0.499F; }
    else if (c.blkEta >= kCoarseEtaBins) { c.blkEta = kCoarseEtaBins - 1;  c.locEta =  1.499F; }

    return c;
}

} // namespace PDC::Geo
#endif /* PDC_GEO_H */
