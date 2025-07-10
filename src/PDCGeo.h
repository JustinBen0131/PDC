/***************************************************************************************************
 *  PDCGeo.h  –  header‑only, zero‑overhead geometry helpers for the sPHENIX EMCal
 *
 *  Public interface (stable):
 *  ────────────────────────────────────────────────────────────────────────────────────────────────
 *      •  foldOnce(u)
 *      •  phi::undoAsh(u,b)
 *      •  eta::undoAsh(u,b[,pitch])
 *      •  phi::blockToGlobal(blk,loc[,offsetFn])
 *      •  eta::blockToGlobal(blk,loc[,offsetFn])
 *      •  finePhiIdx(blk,loc) / fineEtaIdx(blk,loc)
 *      •  struct LocalCoord  +  computeLocal(twrEta,twrPhi,twrE)
 *
 *  All helpers are `constexpr` & `inline` ➜ zero ABI / linkage overhead.
 ***************************************************************************************************/
#ifndef PDC_GEO_H
#define PDC_GEO_H   1

/* ── C++ / ROOT ───────────────────────────────────────────────────────── */
#include <cmath>        // fabs, sinh, asinh, copysign, fmod, …
#include <algorithm>    // clamp, max_element
#include <limits>       // quiet_NaN
#include <type_traits>  // std::true_type / std::false_type
#include <vector>
#include <TVector2.h>   // TVector2::Phi_mpi_pi

namespace PDC::Geo
{
/* ═══════════════════════════════════════════════════════════════════════
 *  1.  Global detector constants                                        *
 * ═════════════════════════════════════════════════════════════════════ */
inline constexpr int   kFinePerBlock   = 2;          ///< 2 × 2 super‑cell
inline constexpr int   kFinePhiBins    = 256;
inline constexpr int   kFineEtaBins    =  96;
inline constexpr int   kCoarsePhiBins  = kFinePhiBins / kFinePerBlock;   ///< 128
inline constexpr int   kCoarseEtaBins  = kFineEtaBins / kFinePerBlock;   ///< 48
inline constexpr float kRadPerFine     = 2.F * static_cast<float>(M_PI) / kFinePhiBins;
inline constexpr float kEtaMin         = -1.1F;      ///< centre of η‑row 0
inline constexpr float kDEtaPerFine    =  2.2F / kFineEtaBins;

/* ═══════════════════════════════════════════════════════════════════════
 *  2.  Low‑level helpers                                                *
 * ═════════════════════════════════════════════════════════════════════ */
[[nodiscard]] inline constexpr float
foldOnce(float u) noexcept
{
    if (u <= -0.5F || u > 1.5F)
    {
        u = std::fmod(u + 2.F, 2.F);        // 0 … <2
        if (u > 1.5F) u -= 2.F;             // (‑0.5 , +1.5]
    }
    return u;
}

/* small wrapper: value‑fold *and* coarse‑index bookkeeping ------------- */
template<bool Periodic>
inline constexpr void
foldAndStep(float& loc, int& coarse, int nCoarse)
{
    if (loc > -0.5F && loc < 1.5F) return;      // already canonical
    const bool crossedLeft = (loc < -0.5F);

    loc = std::fmod(loc + 2.F, 2.F);
    if (loc >= 1.5F) loc -= 2.F;

    coarse += crossedLeft ? -1 : +1;

    if constexpr (Periodic)                     // φ wraps, η clamps later
    {
        if (coarse < 0)               coarse += nCoarse;
        else if (coarse >= nCoarse)   coarse -= nCoarse;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  3.  Ash‑distortion inverse (generic template)                        *
 * ═════════════════════════════════════════════════════════════════════ */
namespace _impl
{
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
    struct NoEdgeFix
    {
        [[nodiscard]] static constexpr float apply(float corr, float) noexcept
        { return corr; }
    };

    template<class EdgePolicy>
    [[nodiscard]] inline constexpr float
    undoAshGeneric(float loc, float b,
                   float towerPitch = 1.F)              // φ passes 1.F
    {
        if (std::fabs(b) < 1e-9F) return loc * towerPitch;
        loc = foldOnce(loc);

        const bool  right = (loc > 0.5F);
        const float xMeas = right ? loc - 1.F : loc;

        const double S    = std::sinh(1. / (2. * b));
        float xTrue       = static_cast<float>(b * std::asinh(2. * S * xMeas));

        float corr = (right ? xTrue + 1.F : xTrue) * towerPitch;
        corr       = EdgePolicy::apply(corr, b);

        if (std::fabs(corr) > 1.6F)                     // stay inside barrel
            corr = std::copysign(1.6F, corr);
        return corr;
    }
} // namespace _impl

/* ═══════════════════════════════════════════════════════════════════════
 *  4.  Public Ash helpers                                               *
 * ═════════════════════════════════════════════════════════════════════ */
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
 *  5.  (blk,loc) → global helpers (axis‑agnostic template)              *
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
 *  6.  Fine‑index helpers                                               *
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
 *  7.  Local 2×2‑block coordinate toolkit                               *
 *      (drop‑in replacement for getBlockCord + getAvg{Eta,Phi})         *
 * ═════════════════════════════════════════════════════════════════════ */
struct LocalCoord
{
    int   blkEta  = 0;          ///< coarse η  index  (0 … 47)
    int   blkPhi  = 0;          ///< coarse φ  index  (0 … 127)
    float locEta  = 0.F;        ///< local η  in (‑0.5 … +1.5]
    float locPhi  = 0.F;        ///< local φ  in (‑0.5 … +1.5]
};

namespace _impl
{
    /* unified single‑pass scanner – returns ΣE, ΣE·η , ref‑φ‑bin & E , ΣE·φ* (unwrapped) */
    struct ScanOut
    {
        double sumE      = 0.0;
        double sumEeta   = 0.0;
        int    refPhiBin = 0;        ///< highest‑E tower (fine bin)
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
            out.sumEeta += Ei * eta[i];

            if (Ei > out.refE) { out.refE = Ei; out.refPhiBin = phi[i]; }
        }
        return out;
    }
} // namespace _impl

template<class IntVec, class FloatVec>
[[nodiscard]] inline LocalCoord
computeLocal(const IntVec& towerEta,
             const IntVec& towerPhi,
             const FloatVec& towerE) noexcept
{
    const std::size_t nT = towerE.size();
    if (nT == 0 || towerEta.size() != nT || towerPhi.size() != nT)
        return {};                                // invalid input ➜ zero‑init

    /* 1) single pass to gather ΣE, ΣE·η and reference φ‑tower */
    const auto scan = _impl::scanTowers(towerEta, towerPhi, towerE);

    if (scan.sumE < 1e-12) return {};            // all towers zero

    const float etaFine = static_cast<float>(scan.sumEeta / scan.sumE);  // 0 … 96

    /* 2) unwrap φ once around the reference, accumulate ΣE·φ */
    constexpr int   Nphi = kFinePhiBins;
    constexpr float Half = Nphi / 2.F;

    double sumEphi = 0.0;
    for (std::size_t i = 0; i < nT; ++i)
    {
        int bin = towerPhi[i];
        int d   = bin - scan.refPhiBin;
        if (d < -Half) bin += Nphi;
        if (d >  Half) bin -= Nphi;

        sumEphi += towerE[i] * bin;
    }
    float phiFine = static_cast<float>(sumEphi / scan.sumE);          // unwrapped
    phiFine       = std::fmod(phiFine + Nphi, Nphi);                  // back to [0,256)

    /* 3) coarse indices and local coords */
    LocalCoord c;
    c.blkEta = static_cast<int>(std::floor(etaFine)) / kFinePerBlock;
    c.blkPhi = static_cast<int>(std::floor(phiFine)) / kFinePerBlock;
    c.locEta = etaFine - c.blkEta * kFinePerBlock;
    c.locPhi = phiFine - c.blkPhi * kFinePerBlock;

    /* 4) exactly one symmetric fold per axis */
    foldAndStep<false>(c.locEta, c.blkEta, kCoarseEtaBins);
    foldAndStep< true>(c.locPhi, c.blkPhi, kCoarsePhiBins);

    /* hard clamp for η (should rarely trigger) */
    if      (c.blkEta < 0)               { c.blkEta = 0;                   c.locEta = -0.499F; }
    else if (c.blkEta >= kCoarseEtaBins) { c.blkEta = kCoarseEtaBins - 1;  c.locEta =  1.499F; }

    return c;
}

} // namespace PDC::Geo
#endif /* PDC_GEO_H */
