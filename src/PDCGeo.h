/***************************************************************************************************
 *  PDCGeo.h  –  header‑only, zero‑overhead geometry helpers for the sPHENIX EMCal
 *
 *  Strictly legacy‑compatible implementation
 *  -----------------------------------------
 *  •  Local 2×2‑block mathematics is an exact, header‑only port of the original
 *     getBlockCord / do*BlockCorr code – bit‑for‑bit identical results.
 *  •  All *global* helpers (blockToGlobal, fineEtaIdx, …) and Ash inverses are
 *     preserved from the modern version and remain constexpr.
 *
 *  Public interface
 *  ----------------
 *      •  foldOnce(u)                             – symmetric single fold
 *      •  foldAndStep(loc,blk)                    – fold + coarse‑index bookkeeping
 *      •  phi::undoAsh / eta::undoAsh             – analytical Ash‑inverse
 *      •  phi::blockToGlobal / eta::blockToGlobal – (blk,loc) → η / φ
 *      •  finePhiIdx / fineEtaIdx                 – tower‑grid indices
 *      •  LocalCoord  +  computeLocal()           – legacy‑accurate block coord
 *      •  doPhiBlockCorr / doEtaBlockCorr         – exact distortion correction
 ***************************************************************************************************/
#ifndef PDC_GEO_H
#define PDC_GEO_H   1

/* ─────────────────────────  Standard / ROOT includes  ───────────────────────── */
#include <cmath>        // fabs, sinh, asinh, copysign, fmod, …
#include <algorithm>    // clamp, max_element
#include <limits>       // quiet_NaN
#include <type_traits>  // true_type / false_type
#include <vector>
#include <TVector2.h>   // TVector2::Phi_mpi_pi

namespace PDC::Geo
{
/* ═══════════════════════════════════════════════════════════════════════
 *  1.  Global detector constants (sPHENIX barrel, fixed)
 * ═════════════════════════════════════════════════════════════════════ */
inline constexpr int   kFinePerBlock   = 2;            ///< side of a 2×2 super‑tower
inline constexpr int   kFinePhiBins    = 256;          ///< full barrel in φ
inline constexpr int   kFineEtaBins    =  96;          ///< barrel height in η
inline constexpr int   kCoarsePhiBins  = kFinePhiBins / kFinePerBlock;   ///< 128
inline constexpr int   kCoarseEtaBins  = kFineEtaBins / kFinePerBlock;   ///<  48
inline constexpr float kRadPerFine     = 2.F * static_cast<float>(M_PI) / kFinePhiBins;
inline constexpr float kEtaMin         = -1.1F;        ///< centre of η‑row 0
inline constexpr float kDEtaPerFine    =  2.2F / kFineEtaBins;           ///< Δη per fine tower

/* ═══════════════════════════════════════════════════════════════════════
 *  2.  Low‑level helpers
 * ═════════════════════════════════════════════════════════════════════ */
[[nodiscard]] inline constexpr float
foldOnce(float u) noexcept
{
    if (u <= -0.5F || u > 1.5F)
    {
        u = std::fmod(u + 2.F, 2.F);          // 0 … <2
        if (u > 1.5F) u -= 2.F;               // −0.5 … +1.5]
    }
    return u;
}

/* fold + coarse‑index bookkeeping */
template<bool Periodic>
inline constexpr void
foldAndStep(float& loc, int& coarse, int nCoarse) noexcept
{
    if (loc > -0.5F && loc <= 1.5F) return;   // already canonical

    const bool rightOverflow = (loc > 1.5F);

    loc = std::fmod(loc + 2.F, 2.F);
    if (loc > 1.5F) loc -= 2.F;

    if (rightOverflow) ++coarse;

    if constexpr (Periodic)                  // φ wrap
        if (coarse >= nCoarse) coarse -= nCoarse;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  3.  Analytical Ash‑inverse (unchanged modern version)
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
        [[nodiscard]] static constexpr float apply(float c, float) noexcept
        { return c; }
    };

    template<class EdgePolicy>
    [[nodiscard]] inline constexpr float
    undoAshGeneric(float loc, float b, float pitch = 1.F)
    {
        if (std::fabs(b) < 1e-9F) return loc * pitch;
        loc = foldOnce(loc);

        const bool  right = (loc > 0.5F);
        const float xMeas = right ? loc - 1.F : loc;

        const double S = std::sinh(1. / (2. * b));
        float xTrue    = static_cast<float>(b * std::asinh(2. * S * xMeas));

        float corr = (right ? xTrue + 1.F : xTrue) * pitch;
        corr       = EdgePolicy::apply(corr, b);

        if (std::fabs(corr) > 1.6F)
            corr = std::copysign(1.6F, corr);
        return corr;
    }
} // namespace _impl

namespace phi
{
    [[nodiscard]] inline constexpr float
    undoAsh(float loc, float b) noexcept
    { return _impl::undoAshGeneric<_impl::PhiEdgeFix>(loc, b); }
}
namespace eta
{
    [[nodiscard]] inline constexpr float
    undoAsh(float loc, float b, float pitch = 1.F) noexcept
    { return _impl::undoAshGeneric<_impl::NoEdgeFix>(loc, b, pitch); }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  4.  Global (blk,loc) → η / φ converters  & fine‑index helpers
 *      – unchanged modern code
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

            const float fine = static_cast<float>(blk) * kFinePerBlock + loc + 0.5F;
            return TVector2::Phi_mpi_pi(fine * kRadPerFine + off());
        }
        else
        {
            if (blk < 0 || blk >= kCoarseEtaBins || !std::isfinite(loc))
                return std::numeric_limits<float>::quiet_NaN();

            const float fine = static_cast<float>(blk) * kFinePerBlock + loc + 0.5F;
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

/* tower‑grid index helpers */
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
        return std::clamp(idx, 0, Nfine - 1);
}
[[nodiscard]] inline constexpr int finePhiIdx(int b,float l) noexcept
{ return fineIdx<kFinePhiBins,true >(b,l); }
[[nodiscard]] inline constexpr int fineEtaIdx(int b,float l) noexcept
{ return fineIdx<kFineEtaBins,false>(b,l); }

/* ═══════════════════════════════════════════════════════════════════════
 *  5.  Legacy‑accurate local 2×2‑block mathematics
 * ═════════════════════════════════════════════════════════════════════ */
struct LocalCoord
{
    int   blkEta  = 0;   ///< coarse η index  (0 … 47)
    int   blkPhi  = 0;   ///< coarse φ index  (0 … 127)
    float locEta  = 0.F; ///< local η in (−0.5 … +1.5]
    float locPhi  = 0.F; ///< local φ in (−0.5 … +1.5]
};

/* energy‑weighted ⟨η⟩ (0 … <96) */
template<class IntVec,class FloatVec>
[[nodiscard]] inline float
getAvgEta(const IntVec& eta,const FloatVec& E)
{
    const std::size_t nT = eta.size();
    if (nT == 0 || nT != E.size()) return 0.F;
    if (nT == 1) return static_cast<float>(eta[0]);

    double sumE=0., sumEeta=0.;
    for (std::size_t i=0;i<nT;++i){ sumE+=E[i]; sumEeta+=E[i]*eta[i]; }
    return (sumE<1e-9)?0.F:static_cast<float>(sumEeta/sumE);
}

/* energy‑weighted ⟨φ⟩ (0 … <256) */
template<class IntVec,class FloatVec>
[[nodiscard]] inline float
getAvgPhi(const IntVec& phi,const FloatVec& E)
{
    constexpr int   Nphi=kFinePhiBins;
    constexpr float Half=Nphi/2.F;

    const std::size_t nT=phi.size();
    if (nT==0||nT!=E.size()) return 0.F;

    std::size_t iRef=0;
    for (std::size_t i=1;i<nT;++i) if (E[i]>E[iRef]) iRef=i;
    const int ref=phi[iRef];

    double sumE=0., sumEphi=0.;
    for (std::size_t i=0;i<nT;++i)
    {
        int p=phi[i];
        int d=p-ref;
        if (d<-Half) p+=Nphi;
        if (d> Half) p-=Nphi;
        sumE   +=E[i];
        sumEphi+=E[i]*p;
    }
    if (sumE<1e-12) return 0.F;
    return static_cast<float>(std::fmod(sumEphi/sumE + Nphi, Nphi));
}

/* computeLocal  –  verbatim port of legacy getBlockCord */
template<class IntVec,class FloatVec>
[[nodiscard]] inline LocalCoord
computeLocal(const IntVec& twrEta,
             const IntVec& twrPhi,
             const FloatVec& twrE)
{
    constexpr int FinePerBlock=kFinePerBlock;   // 2
    constexpr int NcoarsePhi  =kCoarsePhiBins;  // 128

    const std::size_t nT=twrE.size();
    if (nT==0||twrEta.size()!=nT||twrPhi.size()!=nT) return {};

    const float etaCoG=getAvgEta(twrEta,twrE);  // 0 … <96
    const float phiCoG=getAvgPhi(twrPhi,twrE);  // 0 … <256

    int blkEta=static_cast<int>(std::floor(etaCoG))/FinePerBlock; // 0 … 47
    int blkPhi=static_cast<int>(std::floor(phiCoG))/FinePerBlock; // 0 … 127

    float locEta=etaCoG-blkEta*FinePerBlock;
    float locPhi=phiCoG-blkPhi*FinePerBlock;

    auto fold=[&](float& loc,int& coarse,bool isPhi)
    {
        if (loc<=-0.5F||loc>1.5F)
        {
            loc=std::fmod(loc+2.F,2.F);
            if (loc>1.5F){ loc-=2.F; ++coarse; }
            if (isPhi&&coarse==NcoarsePhi) coarse=0;
        }
    };
    fold(locEta,blkEta,false);
    fold(locPhi,blkPhi,true );

    return {blkEta,blkPhi,locEta,locPhi};
}

/* legacy block‑distortion inverses */
[[nodiscard]] inline float
doPhiBlockCorr(float localPhi,float b)
{
    if (std::fabs(b)<1e-9F) return localPhi;
    if (localPhi<=-0.5F||localPhi>1.5F)
        localPhi=std::fmod(localPhi+2.F,2.F);

    const float Xmeas=(localPhi<0.5F)?localPhi:localPhi-1.F;
    const double S   =std::sinh(1.0/(2.0*b));
    const float Xtrue=static_cast<float>(b*std::asinh(2.0*S*Xmeas));
    float corr       =(localPhi<0.5F)?Xtrue:Xtrue+1.F;

    if (corr<=-0.5F||corr>1.5F){
        corr=std::fmod(corr+2.F,2.F);
        if (corr>1.5F) corr-=2.F;
    }
    return corr;
}

[[nodiscard]] inline float
doEtaBlockCorr(float localEta,float b,float dEtaTower=1.F)
{
    if (std::fabs(b)<1e-9F) return localEta;
    if (localEta<=-0.5F||localEta>1.5F)
        localEta=std::fmod(localEta+2.F,2.F);

    const float Xmeas=(localEta<0.5F)?localEta:localEta-1.F;
    const double S   =std::sinh(1.0/(2.0*b));
    const float Xtrue=static_cast<float>(b*std::asinh(2.0*S*Xmeas));
    float corr       =(localEta<0.5F)?Xtrue:Xtrue+1.F;

    if (corr<=-0.5F||corr>1.5F){
        corr=std::fmod(corr+2.F,2.F);
        if (corr>1.5F) corr-=2.F;
    }
    corr*=dEtaTower;
    if (std::fabs(corr)>1.6F) corr=std::copysign(1.6F,corr);
    return corr;
}

/* ------------------------------------------------------------------ *
 *  Compatibility shim – analytical Ash inverse **plus**               *
 *  coarse‑index bookkeeping exactly as getBlockCord used to do        *
 * ------------------------------------------------------------------ */
template<bool Periodic>
struct CorrOut          // tiny POD returned to the caller
{
    float loc;          // folded local coordinate  (-0.5 … +1.5]
    int   blk;          // updated coarse‑block index
};

template<bool Periodic>
[[nodiscard]] inline constexpr CorrOut<Periodic>
undoAshAndReindex(float locIn,        // original local coord
                  int   blkIn,        // original coarse index
                  float b) noexcept   // Ash parameter
{
    /* 1) analytical inverse (periodic ⇔ φ) */
    float loc = Periodic
                  ? phi::undoAsh(locIn, b)
                  : eta::undoAsh(locIn, b);

    /* 2) keep (blk,loc) pair self‑consistent */
    int blk = blkIn;
    if (loc <= -0.5F) { loc += 2.F; --blk; }
    if (loc >  1.5F) { loc -= 2.F; ++blk; }

    /* 3) φ wraps, η clamps */
    if constexpr (Periodic)
    {
        if (blk < 0)               blk += kCoarsePhiBins;
        if (blk >= kCoarsePhiBins) blk -= kCoarsePhiBins;
    }
    return { foldOnce(loc), blk };
}


} // namespace PDC::Geo
#endif /* PDC_GEO_H */
