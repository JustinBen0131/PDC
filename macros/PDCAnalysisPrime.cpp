#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH2F.h>
#include <regex>
#include <TH3F.h>
#include <TString.h>
#include <TVector2.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <unordered_map>
#include <TMath.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// ========================= Global switches & energy bins ====================

static bool isFirstPass = false;  // <-- set by the user request

// Binning-mode enum used by the playground helpers
enum class EBinningMode { kRange, kDiscrete };

// Tiny shim: your new file defines EnsureDir() but the playground
// helpers call ensureDir(). Make both resolve to the same thing.
static inline void ensureDir(const char* p)               { gSystem->mkdir(p,true); }
static inline void ensureDir(const std::string& p)        { gSystem->mkdir(p.c_str(), true); }
static inline void EnsureDir(const std::string& p)        { gSystem->mkdir(p.c_str(), true); } // you already had this one

// Marker/color tables used by some summary plots
static const Color_t kCol[5] = { kBlue+1, kBlue+1, kRed+1, kRed+1, kRed+1 };
static const Style_t kMk [5] = { 20,      24,      20,     24,     21      };


/**
 * \brief 8 energy bins: E_edges[] = {2,4,6,8,10,12,15,20,30}
 * => N_E=8
 */
constexpr double E_edges[] = {2,4,6,8,10,12,15,20,30};
constexpr int    N_E       = (sizeof(E_edges)/sizeof(E_edges[0])) - 1;

// ============================== Small helpers ===============================

struct BRes {
  double val  = std::numeric_limits<double>::quiet_NaN();  // best-fit b
  double err  = 0.0;                                       // σ_b
  double norm = std::numeric_limits<double>::quiet_NaN();  // best-fit Norm
};

// Variant keys and pretty labels for legends/titles
static const std::map<std::string, std::string> kEtaPretty = {
  {"originalEta", "no #eta-dep"},
  {"fullEta",     "|#eta| #leq 1.10"},
  {"etaCore",     "|#eta| #leq 0.20"},
  {"etaMid",      "0.20 < |#eta| #leq 0.70"},
  {"etaEdge",     "0.70 < |#eta| #leq 1.10"}
};

// Pretty labels for |z_vtx| slices (COARSE) and fixed draw order
static const std::map<std::string, std::string> kZPretty = {
  {"z00to10", "|z_{vtx}| 0-10 cm"},
  {"z10to20", "|z_{vtx}| 10-20 cm"},
  {"z20to30", "|z_{vtx}| 20-30 cm"},
  {"z30to45", "|z_{vtx}| 30-45 cm"},
  {"z45to60", "|z_{vtx}| 45-60 cm"}
};
static const std::vector<std::string> kZOrder = {
  "z00to10","z10to20","z20to30","z30to45","z45to60"
};

// Pretty labels for |z_vtx| slices (FINE 0-10 cm: 0-2,2-4,4-6,6-8,8-10) and fixed draw order
static const std::map<std::string, std::string> kZFinePretty = {
  {"z00to02", "|z_{vtx}| 0-2 cm"},
  {"z02to04", "|z_{vtx}| 2-4 cm"},
  {"z04to06", "|z_{vtx}| 4-6 cm"},
  {"z06to08", "|z_{vtx}| 6-8 cm"},
  {"z08to10", "|z_{vtx}| 8-10 cm"}
};
static const std::vector<std::string> kZFineOrder = {
  "z00to02","z02to04","z04to06","z06to08","z08to10"
};

static inline std::vector<std::pair<double,double>> MakeEnergySlices()
{
  std::vector<std::pair<double,double>> v; v.reserve(N_E);
  for (int i=0;i<N_E;++i) v.emplace_back(E_edges[i], E_edges[i+1]);
  return v;
}

static inline std::vector<double> MakeEnergyCenters()
{
  std::vector<double> c; c.reserve(N_E);
  for (int i=0;i<N_E;++i) c.push_back(0.5*(E_edges[i]+E_edges[i+1]));
  return c;
}

// Try a list of names and return the first TH3F found (or nullptr).
static TH3F* GetTH3FByNames(TFile* f, const std::vector<TString>& names)
{
  for (const auto& n : names)
  {
    TH3F* h = dynamic_cast<TH3F*>( f->Get(n) );
    if (h) return h;
  }
  return nullptr;
}


// ==========================================================================
// π0 mass study (standalone). Call from PDCAnalysisPrime() like:
//   RunPi0MassAnalysis(inFilePath.c_str(), outBaseDir.c_str());
// ==========================================================================

void RunPi0MassAnalysis(const char* inFilePath, const char* outBaseDir)
{
  // ---- basic guards --------------------------------------------------------
  if (!inFilePath || !outBaseDir) {
    std::cerr << "[Pi0MassAna] bad args\n";
    return;
  }
  std::unique_ptr<TFile> fin( TFile::Open(inFilePath, "READ") );
  if (!fin || fin->IsZombie()) {
    std::cerr << "[Pi0MassAna][ERROR] cannot open " << inFilePath << "\n";
    return;
  }
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  // ---- variants & cosmetics -----------------------------------------------
  static const char* kVarTags[]   = { "CLUSRAW","CLUSCP","EAgeom","EAetaE","EAEonly","EAmix","PDCraw","PDCcorr" };
  static const char* kVarPretty[] = {
    "No Correction",                          // CLUSRAW
    "Legacy Constant b-corr",                 // CLUSCP
    "Energy/|z_{vtx}|-Dep b-correction",      // EAgeom  <-- UPDATED LABEL
    "Energy/|#eta|-Dep b-correction",         // EAetaE
    "Energy-Dep b-correction",                // EAEonly
    "EA-mix",                                 // EAmix
    "PDC-raw",                                // PDCraw
    "PDC-corr"                                // PDCcorr
  };
  static const Color_t kVarColor[] = { kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kGray+2 };
  constexpr int NVAR = sizeof(kVarTags)/sizeof(kVarTags[0]);
    // Saving-name alias: map 'EAgeom' -> 'EAzE' for all output paths & filenames
    auto kSaveTag = [&](int v)->const char* {
      return (std::string(kVarTags[v])=="EAgeom") ? "EAzE" : kVarTags[v];
  };

  struct ViewCfg { const char* key; const char* sfx; const char* pretty; };
  static const ViewCfg kViews[] = {
    { "originalEta", "",        "no #eta-dep" },
    { "fullEta",     "fullEta", "|#eta| #leq 1.10" },
    { "etaCore",     "etaCore", "|#eta| #leq 0.20" },
    { "etaMid",      "etaMid",  "0.20 < |#eta| #leq 0.70" },
    { "etaEdge",     "etaEdge", "0.70 < |#eta| #leq 1.10" }
  };
  constexpr int N_VIEW = sizeof(kViews)/sizeof(kViews[0]);

  // pretty energy labels / file-safe labels
  auto eLabel  = [&](int i)->TString { return Form("E_%g_to_%g", E_edges[i], E_edges[i+1]); };
  auto ePretty = [&](int i)->TString { return Form("%g #leq E < %g", E_edges[i], E_edges[i+1]); };
  auto eCtrOf  = [&](int i)->double  { return 0.5*(E_edges[i] + E_edges[i+1]); };

  // ---------- output root (ONLY view folders at top level) -----------------
  const TString baseOut = TString(outBaseDir) + "/pi0MassAna";
  gSystem->mkdir(baseOut, true);

  // ---------- robust histogram name candidates -----------------------------
  auto nameCandsView = [&](const char* tag, double elo, double ehi, const char* sfx) {
    std::vector<TString> v;
    const char* fmts[] = { "%g","%.0f","%.1f","%.2f","%.3f" };
    if (!sfx || !*sfx) {
      // base (no η suffix): h_m_pi0_<tag>_<elo>_<ehi>
      for (auto f1 : fmts) for (auto f2 : fmts)
        v.emplace_back( Form("h_m_pi0_%s_%s_%s", tag, Form(f1, elo), Form(f2, ehi)) );
    } else {
      // Try BOTH historical orderings:
      //   A) h_m_pi0_<tag>_<elo>_<ehi>_<sfx>
      //   B) h_m_pi0_<tag>_<sfx>_<elo>_<ehi>   (booking style you showed)
      for (auto f1 : fmts) for (auto f2 : fmts) {
        v.emplace_back( Form("h_m_pi0_%s_%s_%s_%s", tag, Form(f1, elo), Form(f2, ehi), sfx) );
        v.emplace_back( Form("h_m_pi0_%s_%s_%s_%s", tag, sfx,            Form(f1, elo), Form(f2, ehi)) );
      }
    }
    return v;
  };

  // Depth-first recursive lookup through subdirectories
  auto getHistByCand = [&](const std::vector<TString>& cand)->TH1*
  {
    std::function<TH1*(TDirectory*, const TString&)> findIn =
      [&](TDirectory* dir, const TString& nm)->TH1*
    {
      if (!dir) return nullptr;

      // direct lookup first
      if (auto* obj = dir->Get(nm)) {
        if (auto* h = dynamic_cast<TH1*>(obj)) return h;
      }

      // then recurse
      TIter next(dir->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(next())) {
        const char* cls = key->GetClassName();
        if (strcmp(cls,"TDirectoryFile")!=0 && strcmp(cls,"TDirectory")!=0) continue;
        TDirectory* sub = dynamic_cast<TDirectory*>(dir->Get(key->GetName()));
        if (!sub) continue;
        if (auto* h = findIn(sub, nm)) return h;
      }
      return nullptr;
    };

    for (const auto& n : cand) {
      if (TH1* h = findIn(fin.get(), n)) return h;
    }
    return nullptr;
  };

  // ---- fit containers ------------------------------------------------------
  struct FitRes { double mean=std::numeric_limits<double>::quiet_NaN(),
                  meanErr=0.0, sig=std::numeric_limits<double>::quiet_NaN(),
                  sigErr=0.0; bool ok=false; };
  struct FitModel {
    double par[7]  = {0,0,0, 0,0,0,0};   // gaus(0..2) + pol3(3..6)
    double perr[7] = {0,0,0, 0,0,0,0};
    FitRes summary;
    bool   have=false;
  };

  // Per-view containers: [view][variant][Ebin]
  std::vector<std::vector<std::vector<FitModel>>> modelV(
    N_VIEW, std::vector<std::vector<FitModel>>(NVAR, std::vector<FitModel>(N_E)));
  std::vector<std::vector<std::vector<FitRes>>> summaryV(
    N_VIEW, std::vector<std::vector<FitRes>>(NVAR, std::vector<FitRes>(N_E)));
  std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>> hKeepV(N_VIEW);
  for (int iv=0; iv<N_VIEW; ++iv) {
    hKeepV[iv].resize(NVAR);
    for (int v=0; v<NVAR; ++v) hKeepV[iv][v].resize(N_E);
  }

  // ---- one-time fit helper -------------------------------------------------
  auto fitOnce_gaus_pol3 = [&](TH1* h, int iv, int v, int i, const TString& saveDir)->FitModel
  {
    FitModel M; if (!h || h->GetEntries()<10) return M;
    const double xlo = std::max(0.040, h->GetXaxis()->GetXmin());
    const double xhi = std::min(0.900, h->GetXaxis()->GetXmax());
    const double mu0 = 0.135, sg0 = 0.020;
    int bmu = h->GetXaxis()->FindBin(mu0);
    double A0 = std::max(h->GetBinContent(bmu), 1.0);

    std::unique_ptr<TF1> f( new TF1("fGPol3","gaus(0)+pol3(3)", xlo, xhi) );
    f->SetParameters(A0, mu0, sg0, 0,0,0,0);
    f->SetParLimits(2, 0.008, 0.060);
    f->SetLineColor(kVarColor[v]);
    f->SetNpx(1000);

    TFitResultPtr fr = h->Fit(f.get(), "QSNR");
    if (fr.Get() && fr->IsValid()) {
      for (int p=0;p<7;++p) { M.par[p]=f->GetParameter(p); M.perr[p]=f->GetParError(p); }
      M.summary.mean=M.par[1]; M.summary.meanErr=M.perr[1];
      M.summary.sig =M.par[2]; M.summary.sigErr =M.perr[2];
      M.summary.ok = std::isfinite(M.summary.mean) && std::isfinite(M.summary.sig);
      M.have = true;
    }

    // save PNG for this [view, variant, E]
    {
      TCanvas c("cFit","Invariant Mass Distribution",900,650);
      c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);

      const TString title = Form("Invariant Mass Distribution (%s)   (%s)   [%s]",
                                 kVarPretty[v], ePretty(i).Data(), kViews[iv].pretty);
      h->SetTitle(title);
      h->GetXaxis()->SetTitle("M_{#gamma#gamma}  [GeV]");
      h->GetYaxis()->SetTitle("Counts");

      h->SetLineColor(kBlack); h->SetMarkerStyle(20); h->SetMarkerSize(0.9);
      h->Draw("E1");

      f->SetLineColor(kRed+1); f->SetLineWidth(3); f->Draw("SAME");

      TF1 fp(Form("fp_%d_%d_%d",iv,v,i),"gaus",xlo,xhi);
      fp.SetParameters(M.par[0], M.par[1], M.par[2]);
      fp.SetLineColor(kBlue+1); fp.SetLineStyle(2); fp.SetLineWidth(2); fp.Draw("SAME");

      TF1 fb(Form("fb_%d_%d_%d",iv,v,i),"pol3",xlo,xhi);
      for(int p=0;p<4;++p) fb.SetParameter(p,M.par[3+p]);
      fb.SetLineColor(kGreen+2); fb.SetLineStyle(2); fb.SetLineWidth(2); fb.Draw("SAME");

      TLegend L(0.58,0.70,0.88,0.88); L.SetBorderSize(0); L.SetFillStyle(0);
      L.AddEntry(h,"data","lep"); L.AddEntry(&fp,"#pi^{0} (gaus)","l");
      L.AddEntry(&fb,"bkg (pol3)","l"); L.AddEntry(f.get(),"total","l"); L.Draw();

      if (M.summary.ok) {
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
        tx.DrawLatex(0.14,0.86,Form("#mu = %.4f #pm %.4f GeV", M.summary.mean, M.summary.meanErr));
        tx.DrawLatex(0.14,0.80,Form("#sigma = %.4f #pm %.4f GeV", M.summary.sig, M.summary.sigErr));
      }
      const TString outP = saveDir + "/" + TString::Format("pi0Mass_%s_%s.png", kSaveTag(v), eLabel(i).Data());

      c.Modified(); c.Update(); c.SaveAs(outP);
    }

    return M;
  };

  // -------------------------- helpers for overlays/summaries ----------------
  auto drawTotalModel = [&](int iv, int v, int i)
  {
    const auto &M = modelV[iv][v][i];
    if (!M.have) return;
    const double xlo = 0.040, xhi = 0.900;

    TF1 *fTot = new TF1(Form("mTot_%d_%d_%d",iv,v,i), "gaus(0)+pol3(3)", xlo, xhi);
    for (int p=0; p<7; ++p) fTot->SetParameter(p, M.par[p]);
    fTot->SetNpx(1000);
    fTot->SetLineColor(kVarColor[v]);
    fTot->SetLineWidth(3);
    fTot->Draw("SAME");
  };

  auto makeGE = [](const std::vector<double>& xs, const std::vector<double>& ys,
                   const std::vector<double>& exs, const std::vector<double>& eys,
                   Color_t col, int mstyle)->std::unique_ptr<TGraphErrors>
  {
    std::vector<double> X,Y,EX,EY;
    for (size_t k=0;k<xs.size();++k) {
      if (!std::isfinite(ys[k])) continue;
      X.push_back(xs[k]); Y.push_back(ys[k]);
      EX.push_back(k<exs.size()?exs[k]:0.0);
      EY.push_back(k<eys.size()?eys[k]:0.0);
    }
    auto g = std::make_unique<TGraphErrors>((int)X.size(),
               X.empty()?nullptr:&X[0], Y.empty()?nullptr:&Y[0],
               EX.empty()?nullptr:&EX[0], EY.empty()?nullptr:&EY[0]);
    g->SetMarkerStyle(mstyle); g->SetMarkerSize(1.2);
    g->SetMarkerColor(col); g->SetLineColor(col);
    return g;
  };

  // view-scoped overlay: use only this view’s histograms/models/summaries
  auto overlay2 = [&](int iv, int vA, int vR, const TString& outDir){
    for (int i=0;i<N_E;++i) {
      if (!hKeepV[iv][vA][i] || !hKeepV[iv][vR][i]) continue;
      TCanvas c("c2","overlay2",900,650); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);
      TH1 *hR=hKeepV[iv][vR][i].get(), *hA=hKeepV[iv][vA][i].get();
      hR->SetTitle(Form("%s vs %s   (%s)   [%s]",
                        kVarPretty[vA], kVarPretty[vR], ePretty(i).Data(), kViews[iv].pretty));
      auto maxWithErr = [](TH1* h)->double {
        if (!h) return 0.0;
        double m = 0.0;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
          const double y = h->GetBinContent(b) + h->GetBinError(b);
          if (y > m) m = y;
        }
        return m;
      };
      const double yMax = std::max(maxWithErr(hR), maxWithErr(hA));
      hR->SetMinimum(0.0);
      hR->SetMaximum(yMax * 1.08);
      hR->Draw("E1"); hA->Draw("E1 SAME");

      drawTotalModel(iv,vR,i);
      drawTotalModel(iv,vA,i);

      TLegend L(0.16, 0.5, 0.48, 0.62);
      L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextSize(0.03);
      L.AddEntry(hR, kVarPretty[vR], "p");
      L.AddEntry(hA, kVarPretty[vA], "p");
      L.Draw();

      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
      double y=0.87, dy=0.047;
      const auto &RR=summaryV[iv][vR][i], &RA=summaryV[iv][vA][i];
      tx.SetTextColor(kVarColor[vR]);
      if (RR.ok) tx.DrawLatex(0.14,y,Form("%s:  #mu=%.4f#pm%.4f,  #sigma=%.4f#pm%.4f",
                          kVarPretty[vR],RR.mean,RR.meanErr,RR.sig,RR.sigErr));
      else tx.DrawLatex(0.14,y,Form("%s:  fit n/a",kVarPretty[vR]));
      y-=dy;
      tx.SetTextColor(kVarColor[vA]);
      if (RA.ok) tx.DrawLatex(0.14,y,Form("%s:  #mu=%.4f#pm%.4f,  #sigma=%.4f#pm%.4f",
                          kVarPretty[vA],RA.mean,RA.meanErr,RA.sig,RA.sigErr));
      else tx.DrawLatex(0.14,y,Form("%s:  fit n/a",kVarPretty[vA]));
      tx.SetTextColor(kBlack);

      const TString outP = outDir + "/" + TString::Format("overlay_%s_vs_%s_%s.png",
                                        kSaveTag(vA), kSaveTag(vR), eLabel(i).Data());

      c.Modified(); c.Update(); c.SaveAs(outP);
    }
  };

  auto overlay3 = [&](int iv, int vB, int vC, int vR, const TString& outDir){
    for (int i=0;i<N_E;++i) {
      if (!hKeepV[iv][vR][i] || !hKeepV[iv][vB][i] || !hKeepV[iv][vC][i]) continue;
      TCanvas c("c3","overlay3",900,650); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);
      TH1* hR=hKeepV[iv][vR][i].get(); TH1* hB=hKeepV[iv][vB][i].get(); TH1* hC=hKeepV[iv][vC][i].get();
      hR->SetTitle(Form("%s vs %s vs %s   (%s)   [%s]",
                        kVarPretty[vR], kVarPretty[vB], kVarPretty[vC], ePretty(i).Data(), kViews[iv].pretty));
      auto maxWithErr = [](TH1* h)->double {
        if (!h) return 0.0;
        double m = 0.0;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
          const double y = h->GetBinContent(b) + h->GetBinError(b);
          if (y > m) m = y;
        }
        return m;
      };
      const double yMax = std::max(maxWithErr(hR), std::max(maxWithErr(hB), maxWithErr(hC)));
      hR->SetMinimum(0.0);
      hR->SetMaximum(yMax * 1.35);
      hR->Draw("E1"); hB->Draw("E1 SAME"); hC->Draw("E1 SAME");

      drawTotalModel(iv,vR,i);
      drawTotalModel(iv,vB,i);
      drawTotalModel(iv,vC,i);

        // Legend without error bars in entries: use "p" (markers only)
        TLegend L(0.58,0.64,0.88,0.89);
        L.SetBorderSize(0);
        L.SetFillStyle(0);
        L.AddEntry(hR, kVarPretty[vR], "p");
        L.AddEntry(hB, kVarPretty[vB], "p");
        L.AddEntry(hC, kVarPretty[vC], "p");
        L.Draw();

        // Per-plot text
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.024);
        double y=0.87, dy=0.047;
        const auto &RR=summaryV[iv][vR][i], &RB=summaryV[iv][vB][i], &RC=summaryV[iv][vC][i];
        tx.SetTextColor(kVarColor[vR]);
        if (RR.ok) tx.DrawLatex(0.14,y,Form("%s: #mu=%.4f#pm%.4f, #sigma=%.4f#pm%.4f",
                            kVarPretty[vR],RR.mean,RR.meanErr,RR.sig,RR.sigErr));
        else       tx.DrawLatex(0.14,y,Form("%s: fit n/a",kVarPretty[vR]));
        y -= dy;

        tx.SetTextColor(kVarColor[vB]);
        if (RB.ok) tx.DrawLatex(0.14,y,Form("%s: #mu=%.4f#pm%.4f, #sigma=%.4f#pm%.4f",
                            kVarPretty[vB],RB.mean,RB.meanErr,RB.sig,RB.sigErr));
        else       tx.DrawLatex(0.14,y,Form("%s: fit n/a",kVarPretty[vB]));
        y -= dy;

        tx.SetTextColor(kVarColor[vC]);
        if (RC.ok) tx.DrawLatex(0.14,y,Form("%s: #mu=%.4f#pm%.4f, #sigma=%.4f#pm%.4f",
                            kVarPretty[vC],RC.mean,RC.meanErr,RC.sig,RC.sigErr));
        else       tx.DrawLatex(0.14,y,Form("%s: fit n/a",kVarPretty[vC]));
        tx.SetTextColor(kBlack);

        // --------- Clean ANSI percentage-only table to terminal for CLUSRAW_CLUSCP_CLUSEAonly ----------
        // Prints Δ% resolution (σ/μ) vs RAW for: Legacy Constant b-corr (CP) and Energy-Dep b-correction (E-only)
        if (outDir.Contains("CLUSRAW_CLUSCP_CLUSEAonly")) {
          // Header once per view (on the first bin)
          if (i == 0) {
            std::cout << "\n\033[1m[Overlay] Δ% resolution vs RAW (σ/μ)  —  View: "
                      << kViews[iv].pretty << "\033[0m\n";
            std::cout << "\033[1m"
                      << std::left  << std::setw(16) << "E-bin"
                      << " | " << std::setw(16) << "Legacy Const-b"
                      << " | " << std::setw(16) << "Energy-Dep b"
                      << "\033[0m\n";
            std::cout << std::string(16 + 3 + 16 + 3 + 16, '-') << "\n";
          }

          auto pct = [](double num)->std::string {
            if (!std::isfinite(num)) return std::string("      n/a      ");
            std::ostringstream os; os.setf(std::ios::fixed);
            os << std::setprecision(2) << std::showpos << num << "%" << std::noshowpos;
            return os.str();
          };

          // Percent change in RESOLUTION (sigma/mu) vs RAW
          auto safe_res = [](const FitRes& R)->double {
            if (!R.ok || !std::isfinite(R.sig) || !std::isfinite(R.mean) || R.sig<=0.0 || R.mean<=0.0)
              return std::numeric_limits<double>::quiet_NaN();
            return R.sig / R.mean;
          };

          const double resRAW = safe_res(RR);  // RAW = No Correction
          const double resCP  = safe_res(RB);  // Legacy Constant b-corr
          const double resEO  = safe_res(RC);  // Energy-Dep b-correction

          double dCP = std::numeric_limits<double>::quiet_NaN();
          double dEO = std::numeric_limits<double>::quiet_NaN();
          if (std::isfinite(resRAW) && resRAW > 0.0 && std::isfinite(resCP))
            dCP = 100.0 * (resCP - resRAW) / resRAW;
          if (std::isfinite(resRAW) && resRAW > 0.0 && std::isfinite(resEO))
            dEO = 100.0 * (resEO - resRAW) / resRAW;

          std::cout << std::left << std::setw(16) << ePretty(i).Data()
                    << " | " << std::setw(16) << pct(dCP)
                    << " | " << std::setw(16) << pct(dEO)
                    << "\n";
        }


        const TString outP = outDir + "/" + TString::Format("overlay3_%s_%s_%s_%s.png",
                                       kSaveTag(vR), kSaveTag(vB), kSaveTag(vC), eLabel(i).Data());

        c.Modified(); c.Update(); c.SaveAs(outP);
    }
  };

  auto twoSeries = [&](int iv, int vB, const char* prettyB, const TString& outSumDir, const TString& baseName)
  {
    // requires CLUSRAW as the baseline series
    int vRAW=-1; for (int v=0; v<NVAR; ++v) if (std::string(kVarTags[v])=="CLUSRAW") vRAW=v;
    if (vRAW<0 || vB<0) return;

    std::vector<double> x,ex(N_E,0.0), yR,eyR,yB,eyB, yRs,eyRs,yBs,eyBs;
    for (int i=0;i<N_E;++i) {
      x.push_back(eCtrOf(i));
      const auto &RR=summaryV[iv][vRAW][i], &RB=summaryV[iv][vB][i];
      if (RR.ok) { yR.push_back(RR.mean); eyR.push_back(RR.meanErr); yRs.push_back(RR.sig); eyRs.push_back(RR.sigErr); }
      else       { yR.push_back(NAN);     eyR.push_back(0.0);        yRs.push_back(NAN);   eyRs.push_back(0.0); }
      if (RB.ok) { yB.push_back(RB.mean); eyB.push_back(RB.meanErr); yBs.push_back(RB.sig); eyBs.push_back(RB.sigErr); }
      else       { yB.push_back(NAN);     eyB.push_back(0.0);        yBs.push_back(NAN);   eyBs.push_back(0.0); }
    }

    auto gR_mean = makeGE(x,yR,ex,eyR,kBlack,20);
    auto gB_mean = makeGE(x,yB,ex,eyB,kRed+1,24);
    auto gR_sig  = makeGE(x,yRs,ex,eyRs,kBlack,20);
    auto gB_sig  = makeGE(x,yBs,ex,eyBs,kRed+1,24);

    if (gR_mean->GetN()==0 && gB_mean->GetN()==0 &&
        gR_sig->GetN()==0  && gB_sig->GetN()==0) {
      std::cout << "[Summaries][" << kViews[iv].key << "] No valid points for " << baseName << " → skipping\n";
      return;
    }

    // mean vs E
    {
      TCanvas c("cMean","mean vs E",900,650); c.SetLeftMargin(0.14); c.SetBottomMargin(0.14);
      double yMin=+1e30,yMax=-1e30,X,Y;
      for (int k=0;k<gR_mean->GetN();++k){ gR_mean->GetPoint(k,X,Y); yMin=std::min(yMin,Y); yMax=std::max(yMax,Y); }
      for (int k=0;k<gB_mean->GetN();++k){ gB_mean->GetPoint(k,X,Y); yMin=std::min(yMin,Y); yMax=std::max(yMax,Y); }
      if (!std::isfinite(yMin) || !std::isfinite(yMax)) { yMin=0.10; yMax=0.16; }
      const double pad=0.15*(yMax-yMin);

      TH2F fr("frM","",100,E_edges[0]-0.5,E_edges[N_E]-0.5,100,yMin-pad,yMax+pad);
      fr.SetTitle(Form("#mu_{#pi^{0}}(E): %s vs %s   [%s];E [GeV];#mu_{#pi^{0}} [GeV]",
                       "No Correction", prettyB, kViews[iv].pretty));
      fr.SetStats(0); fr.Draw();
      gR_mean->Draw("P SAME"); gB_mean->Draw("P SAME");
      TLegend L(0.60,0.74,0.88,0.89); L.SetBorderSize(0); L.SetFillStyle(0);
      L.AddEntry(gR_mean.get(),"No Correction","p");
      L.AddEntry(gB_mean.get(),prettyB,"p"); L.Draw();
      c.SaveAs(outSumDir + "/" + baseName + "_Mean_vs_E.png");
    }
    // sigma vs E
    {
      TCanvas c("cSig","sigma vs E",900,650); c.SetLeftMargin(0.14); c.SetBottomMargin(0.14);
      double yMin=+1e30,yMax=-1e30,X,Y;
      for (int k=0;k<gR_sig->GetN();++k){ gR_sig->GetPoint(k,X,Y); yMin=std::min(yMin,Y); yMax=std::max(yMax,Y); }
      for (int k=0;k<gB_sig->GetN();++k){ gB_sig->GetPoint(k,X,Y); yMin=std::min(yMin,Y); yMax=std::max(yMax,Y); }
      if (!std::isfinite(yMin) || !std::isfinite(yMax)) { yMin=0.010; yMax=0.030; }
      const double pad=0.20*(yMax-yMin);

      TH2F fr("frS","",100,E_edges[0]-0.5,E_edges[N_E]-0.5,100,std::max(0.0,yMin-pad),yMax+pad);
      fr.SetTitle(Form("#sigma_{#pi^{0}}(E): %s vs %s   [%s];E [GeV];#sigma_{#pi^{0}} [GeV]",
                       "No Correction", prettyB, kViews[iv].pretty));
      fr.SetStats(0); fr.Draw();
      gR_sig->Draw("P SAME"); gB_sig->Draw("P SAME");
      TLegend L(0.60,0.74,0.88,0.89); L.SetBorderSize(0); L.SetFillStyle(0);
      L.AddEntry(gR_sig.get(),"No Correction","p");
      L.AddEntry(gB_sig.get(),prettyB,"p"); L.Draw();
      c.SaveAs(outSumDir + "/" + baseName + "_Sigma_vs_E.png");
    }
  };

  // -------------------------- make per-view folder skeletons ----------------
  // (Only views exist at top level; each view contains everything else.)
  std::map<TString, std::map<TString,TString>> overlayDirs; // overlayDirs[viewKey][name] = path
  for (int iv=0; iv<N_VIEW; ++iv) {
    const TString vroot = baseOut + "/" + kViews[iv].key;
    gSystem->mkdir(vroot, true);
    gSystem->mkdir(vroot + "/Summaries", true);
    gSystem->mkdir(vroot + "/zAnalysis", true);
    // per-variant folders
    for (int v=0; v<NVAR; ++v) gSystem->mkdir(vroot + "/" + TString(kSaveTag(v)), true);

      const std::vector<TString> names = {
        "CLUSRAW_EAonly", "CLUSRAW_EAetaE", "CLUSRAW_CLUSCP",
        "CLUSRAW_CLUSCP_CLUSEAonly", "CLUSRAW_CLUSCP_CLUSEAetaE",
        "CLUSCP_EAonly", "CLUSCP_EAetaE",
        "EAzE_EAetaE_EAonly"  // NEW: z_vtx-dep vs |eta|-dep vs E-only overlays
      };

    for (const auto& nm : names) {
      const TString p = vroot + "/" + nm;
      gSystem->mkdir(p, true);
      overlayDirs[kViews[iv].key][nm] = p;
    }
  }

  // -------------------------- entries table & CSV (once) --------------------
  {
    // Map tags to indices
    std::map<std::string,int> tag2idx;
    for (int v=0; v<NVAR; ++v) tag2idx[ kVarTags[v] ] = v;

    auto viewIndex = [](const std::string& s)->int {
      if (s=="fullEta") return 1;
      if (s=="etaCore") return 2;
      if (s=="etaMid")  return 3;
      if (s=="etaEdge") return 4;
      return 0; // originalEta/base case
    };

    auto isNumber = [](const std::string& s)->bool {
      char* end=nullptr;
      errno = 0;
      std::strtod(s.c_str(), &end);
      return end && *end=='\0' && errno==0;
    };

    auto findEbin = [&](double lo, double hi)->int {
      for (int i=0; i<N_E; ++i) {
        if (std::fabs(lo - E_edges[i])   < 1e-6 &&
            std::fabs(hi - E_edges[i+1]) < 1e-6) return i;
      }
      return -1;
    };

    // Accumulators per view
    std::vector<std::vector<std::vector<Long64_t>>> cnt(N_VIEW, std::vector<std::vector<Long64_t>>(NVAR, std::vector<Long64_t>(N_E,0)));

    // One recursive walk across the whole file
    std::function<void(TDirectory*)> walk = [&](TDirectory* dir)
    {
      if (!dir) return;
      TIter next(dir->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(next()))
      {
        const char* cls = key->GetClassName();
        const TString nm = key->GetName();

        if (!strcmp(cls,"TDirectoryFile") || !strcmp(cls,"TDirectory")) {
          if (auto* sub = dynamic_cast<TDirectory*>(dir->Get(nm))) walk(sub);
          continue;
        }
        if (!nm.BeginsWith("h_m_pi0_")) continue;

        std::unique_ptr<TObject> obj( dir->Get(nm) );
        TH1* h = dynamic_cast<TH1*>(obj.get());
        if (!h) continue;

        TString rest = nm; rest.ReplaceAll("h_m_pi0_","");
        std::vector<std::string> tok;
        std::unique_ptr<TObjArray> arr( rest.Tokenize("_") );
        if (arr) {
          for (int it=0; it<arr->GetEntriesFast(); ++it) {
            if (auto* os = dynamic_cast<TObjString*>(arr->At(it)))
              tok.emplace_back( os->GetString().Data() );
          }
        }
        if (tok.size() < 3) continue;

        auto itV = tag2idx.find(tok[0]);
        if (itV == tag2idx.end()) continue;
        int vIdx = itV->second;

        // supported shapes:
        //  TAG, Elo, Ehi
        //  TAG, Elo, Ehi, VIEW
        //  TAG, VIEW, Elo, Ehi
        double lo=0, hi=0; int iE=-1; int iView=0;

        if (tok.size()==3 && isNumber(tok[1]) && isNumber(tok[2])) {
          lo = std::stod(tok[1]); hi = std::stod(tok[2]);
          iE = findEbin(lo,hi);
          iView = 0; // originalEta
        } else if (tok.size()==4) {
          if (isNumber(tok[1]) && isNumber(tok[2])) {
            lo = std::stod(tok[1]); hi = std::stod(tok[2]);
            iE = findEbin(lo,hi);
            iView = viewIndex(tok[3]);
          } else {
            iView = viewIndex(tok[1]);
            if (!isNumber(tok[2]) || !isNumber(tok[3])) continue;
            lo = std::stod(tok[2]); hi = std::stod(tok[3]);
            iE = findEbin(lo,hi);
          }
        }
        if (iE>=0) cnt[iView][vIdx][iE] += static_cast<Long64_t>( h->GetEntries() );
      }
    };
    walk(fin.get());

    auto ebinPretty = [&](int i)->TString { return Form("%g-%g", E_edges[i], E_edges[i+1]); };

    std::cout << "\n======== π0 mass histogram entry counts (per view) ========\n";
    std::cout << "view, variant, E-bin, entries\n";

    // write CSV under originalEta/Summaries so top level stays view-only
    const TString csvPath = baseOut + "/originalEta/Summaries/pi0Mass_hist_counts.csv";
    gSystem->mkdir((baseOut + "/originalEta/Summaries"), true);
    std::ofstream csv(csvPath.Data());
    csv << "view,variant,Elo,Ehi,entries\n";

    for (int iv=0; iv<N_VIEW; ++iv)
      for (int v=0; v<NVAR; ++v)
        for (int i=0; i<N_E; ++i) {
          std::cout << kViews[iv].key << ", " << kVarPretty[v] << ", " << ebinPretty(i)
                    << ", " << cnt[iv][v][i] << "\n";
          csv << kViews[iv].key << ","
              << kVarPretty[v] << ","
              << E_edges[i] << "," << E_edges[i+1] << ","
              << cnt[iv][v][i] << "\n";
        }

    csv.close();
    std::cout << "---------------------------------------------------------------------\n";
    std::cout << "[counts] CSV written to: " << csvPath << "\n\n";
  }

  // -------------------------- PASS: fits & overlays per view ----------------
  for (int iv=0; iv<N_VIEW; ++iv) {
    const TString vroot   = baseOut + "/" + kViews[iv].key;
    const TString outSum  = vroot + "/Summaries";
    const TString outZ    = vroot + "/zAnalysis";

    // 1) per-variant fits saved in  <view>/<variant>/
    for (int v=0; v<NVAR; ++v) {
      const TString saveDir = vroot + "/" + TString(kSaveTag(v));
      for (int i=0; i<N_E; ++i) {
        // get this view’s histogram for [variant, E bin]
        std::vector<TString> cands = nameCandsView(kVarTags[v], E_edges[i], E_edges[i+1], kViews[iv].sfx);
        TH1* hOrig = getHistByCand(cands);
        if (!hOrig || hOrig->GetEntries()<=0 || hOrig->Integral()<=0.0) continue;

        // Fit image (detached clone)
        std::unique_ptr<TH1> hFit( static_cast<TH1*>(hOrig->Clone(
                                   Form("h_%s_%s_%s", kVarTags[v], eLabel(i).Data(), kViews[iv].key))) );
        hFit->SetDirectory(nullptr);
        hFit->GetXaxis()->SetTitle("M_{#gamma#gamma}  [GeV]");
        hFit->GetYaxis()->SetTitle("Counts");

        modelV[iv][v][i]   = fitOnce_gaus_pol3(hFit.get(), iv, v, i, saveDir);
        summaryV[iv][v][i] = modelV[iv][v][i].summary;

        // raw copy for overlays
        std::unique_ptr<TH1> hCopy( static_cast<TH1*>(hOrig->Clone(
                                     Form("hkeep_%s_%s_%s", kVarTags[v], eLabel(i).Data(), kViews[iv].key))) );
        hCopy->SetDirectory(nullptr);
        hCopy->SetLineColor(kVarColor[v]);
        hCopy->SetMarkerColor(kVarColor[v]);
        hCopy->SetMarkerStyle( v==0 ? 20 : 24 );
        hKeepV[iv][v][i] = std::move(hCopy);
      }
    }

    // 2) per-view overlays into view-scoped folders
    int vRAW=-1, vCP=-1, vEAetaE=-1, vEAonly=-1, vEAgeom=-1;
    for (int v=0; v<NVAR; ++v) {
      if (std::string(kVarTags[v])=="CLUSRAW") vRAW=v;
      if (std::string(kVarTags[v])=="CLUSCP")  vCP=v;
      if (std::string(kVarTags[v])=="EAetaE")  vEAetaE=v;
      if (std::string(kVarTags[v])=="EAEonly") vEAonly=v;
      if (std::string(kVarTags[v])=="EAgeom")  vEAgeom=v;
    }
      if (vRAW>=0 && vEAonly>=0) overlay2(iv, vEAonly, vRAW, overlayDirs[kViews[iv].key]["CLUSRAW_EAonly"]);
      if (vRAW>=0 && vEAetaE>=0) overlay2(iv, vEAetaE, vRAW, overlayDirs[kViews[iv].key]["CLUSRAW_EAetaE"]);
      if (vRAW>=0 && vCP>=0)     overlay2(iv, vCP,     vRAW, overlayDirs[kViews[iv].key]["CLUSRAW_CLUSCP"]);

      if (vCP>=0 && vEAonly>=0)  overlay2(iv, vEAonly, vCP,  overlayDirs[kViews[iv].key]["CLUSCP_EAonly"]);
      if (vCP>=0 && vEAetaE>=0)  overlay2(iv, vEAetaE, vCP,  overlayDirs[kViews[iv].key]["CLUSCP_EAetaE"]);

      if (vRAW>=0 && vCP>=0 && vEAonly>=0) overlay3(iv, vCP, vEAonly, vRAW, overlayDirs[kViews[iv].key]["CLUSRAW_CLUSCP_CLUSEAonly"]);
      if (vRAW>=0 && vCP>=0 && vEAetaE>=0) overlay3(iv, vCP, vEAetaE, vRAW, overlayDirs[kViews[iv].key]["CLUSRAW_CLUSCP_CLUSEAetaE"]);

      // NEW: triple overlay of z_vtx-dep, |eta|-dep, and E-only (no CLUSRAW/CLUSCP)
      if (vEAgeom>=0 && vEAetaE>=0 && vEAonly>=0)
        overlay3(iv, vEAetaE, vEAonly, vEAgeom, overlayDirs[kViews[iv].key]["EAzE_EAetaE_EAonly"]);


    // 3) per-view mean/sigma summary curves into <view>/Summaries/
    if (vCP>=0)     twoSeries(iv, vCP,     "Legacy Constant b-corr",         outSum, "CLUSRAW_vs_CLUSCP");
    if (vEAgeom>=0) twoSeries(iv, vEAgeom, "Energy/|z_{vtx}|-Dep b-correction", outSum, "CLUSRAW_vs_CLUSCPEA"); // label updated, filename unchanged

    // 4) z-correlation saved into <view>/zAnalysis/ (same plot for every view)
    {
      TH2F* h2 = dynamic_cast<TH2F*>( fin->Get("h2_truthReco_vz") );
      TH1*  hReco1D = dynamic_cast<TH1*>( fin->Get("h_reco_vz") );

      Long64_t nRecoEntries = 0;
      if (hReco1D) {
        nRecoEntries = static_cast<Long64_t>( hReco1D->GetEntries() );
      } else if (h2) {
        std::unique_ptr<TH1> projY(h2->ProjectionY("_py_tmp"));
        nRecoEntries = projY ? static_cast<Long64_t>(projY->GetEntries()) : 0;
      }

      if (h2 && h2->GetEntries() > 0) {
        TCanvas cZ("cZ","Truth vs Reco Vertex Z", 900, 750);
        cZ.SetLeftMargin(0.12); cZ.SetBottomMargin(0.12);

        h2->SetTitle(Form("Truth vs Reco Vertex Z   [%s];z_{truth} (cm);z_{reco} (cm)", kViews[iv].pretty));
        h2->Draw("COLZ");

        const double xlo = h2->GetXaxis()->GetXmin();
        const double xhi = h2->GetXaxis()->GetXmax();
        TLine diag(xlo, xlo, xhi, xhi);
        diag.SetLineColor(kRed+1);
        diag.SetLineStyle(2);
        diag.SetLineWidth(2);
        diag.Draw("SAME");

        const double rho = h2->GetCorrelationFactor();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.DrawLatex(0.16, 0.92, Form("#rho = %.3f   entries = %.0f   (reco 1D entries = %lld)",
                                      rho, h2->GetEntries(), static_cast<long long>(nRecoEntries)));

        cZ.SaveAs( (outZ + "/truth_vs_reco_vz.png").Data() );
      } else {
        std::cout << "[zAnalysis] WARNING: h2_truthReco_vz not found or empty in input file.\n";
      }
    }
  } // end per-view loop

  // -------------------------- Sigma overlay across variants (all E bins) --------------------------
  {
    // Required variant indices (tags from kVarTags)
    int vRAW   = -1; // "CLUSRAW"
    int vCP    = -1; // "CLUSCP"
    int vEAone = -1; // "EAEonly"
    int vEAeta = -1; // "EAetaE"
    int vEAgeo = -1; // "EAgeom"  <-- NEW: z_vtx dep
    for (int v = 0; v < NVAR; ++v) {
      const std::string tag = kVarTags[v];
      if      (tag == "CLUSRAW") vRAW   = v;
      else if (tag == "CLUSCP")  vCP    = v;
      else if (tag == "EAEonly") vEAone = v;
      else if (tag == "EAetaE")  vEAeta = v;
      else if (tag == "EAgeom")  vEAgeo = v;
    }

    if (vRAW < 0 || vCP < 0 || vEAone < 0 || vEAeta < 0 || vEAgeo < 0) {
      std::cerr << "[SigmaOverlay] Missing one of required variants (CLUSRAW, CLUSCP, EAEonly, EAetaE, EAgeom) — skipping overlay.\n";
    } else if (N_E <= 0) {
      std::cerr << "[SigmaOverlay] No energy bins available — skipping overlay.\n";
    } else {
      // Loop over ALL energy bins
      for (int iE = 0; iE < N_E; ++iE)
      {
        // X positions for the five categories and their display labels
        std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
        const char* xlabels[5] = {
          kVarPretty[vRAW],
          kVarPretty[vCP],
          kVarPretty[vEAone],
          kVarPretty[vEAeta],
          kVarPretty[vEAgeo]
        };

        // Helper to build a TGraphErrors for a given view index, with an x-offset to separate series within each category
        auto makeViewGraph = [&](int iview, double xoff, Style_t mstyle, Color_t col)->std::unique_ptr<TGraphErrors>
        {
          std::vector<double> y, ey;
          auto push = [&](const FitRes& R) {
            if (R.ok && std::isfinite(R.sig)) { y.push_back(R.sig); ey.push_back(R.sigErr); }
            else                               { y.push_back(std::numeric_limits<double>::quiet_NaN()); ey.push_back(0.0); }
          };

          // Order matches x: CLUSRAW, CLUSCP, EAEonly, EAetaE, EAgeom
          push(summaryV[iview][vRAW][iE]);
          push(summaryV[iview][vCP][iE]);
          push(summaryV[iview][vEAone][iE]);
          push(summaryV[iview][vEAeta][iE]);
          push(summaryV[iview][vEAgeo][iE]);

          // Filter NaNs while preserving x order; apply per-series x offset
          std::vector<double> xf, yf, exf, eyf;
          for (size_t k = 0; k < y.size(); ++k) {
            if (!std::isfinite(y[k])) continue;
            xf.push_back(x[k] + xoff);
            yf.push_back(y[k]);
            exf.push_back(0.0);
            eyf.push_back(ey[k]);
          }

          auto g = std::make_unique<TGraphErrors>(
            (int)xf.size(),
            xf.empty()?nullptr:&xf[0],
            yf.empty()?nullptr:&yf[0],
            exf.empty()?nullptr:&exf[0],
            eyf.empty()?nullptr:&eyf[0]
          );
          g->SetMarkerStyle(mstyle);
          g->SetMarkerSize(1.5);
          g->SetMarkerColor(col);
          g->SetLineColor(col);
          g->SetLineWidth(2);
          return g;
        };

        // Build graphs with small per-series x offsets (bin-centered at 1..5)
        // Order (left→right): etaCore (−0.18), etaMid (−0.06), etaEdge (+0.06), fullEta (+0.18)
        // kViews[]: 0:originalEta, 1:fullEta, 2:etaCore, 3:etaMid, 4:etaEdge
        auto g_core = makeViewGraph(2, -0.18, 20, kBlue+1);
        auto g_mid  = makeViewGraph(3, -0.06, 21, kGreen+2);
        auto g_edge = makeViewGraph(4,  0.06, 22, kMagenta+1);
        auto g_full = makeViewGraph(1,  0.18, 24, kRed+1);

        // Determine dynamic y-range from the available points, including error bars
        double ymin = +1e30, ymax = -1e30, X, Y;
        auto upd = [&](TGraphErrors* g) {
          if (!g) return;
          for (int k = 0; k < g->GetN(); ++k) {
            g->GetPoint(k, X, Y);
            const double eY = g->GetErrorY(k);
            ymin = std::min(ymin, Y - eY);
            ymax = std::max(ymax, Y + eY);
          }
        };
        upd(g_edge.get()); upd(g_mid.get()); upd(g_core.get()); upd(g_full.get());
        if (!(std::isfinite(ymin) && std::isfinite(ymax))) { ymin = 0.010; ymax = 0.040; }
        double span = ymax - ymin;
        if (span <= 0) { span = 1e-4; }
        const double padBottom = std::max(0.06 * span, 5e-4);
        const double padTop    = std::max(0.10 * span, 8e-4);
        ymin = std::max(0.0, ymin - padBottom);
        ymax = ymax + padTop;

        // Frame with categorical x-axis (5 bins)
        TCanvas cS(Form("cSigmaOverlay_%02d", iE), "Gaussian sigma (E bin)", 1000, 700);
        cS.SetLeftMargin(0.18);
        cS.SetBottomMargin(0.18);

        TH1F fr("fr_sigma", "", 5, 0.5, 5.5);
        fr.SetStats(0);
        fr.GetYaxis()->SetTitle("Gaussian #sigma  [GeV]");
        fr.GetYaxis()->SetLabelSize(0.025);
        fr.GetYaxis()->SetRangeUser(ymin, ymax);
        fr.GetXaxis()->SetTitle("");
        fr.GetXaxis()->SetTitleSize(0.0);   // hide x-axis title
        fr.GetXaxis()->SetLabelSize(0.032); // smaller category labels
        fr.GetYaxis()->SetTitleOffset(1.55);
        for (int ib = 1; ib <= 5; ++ib) fr.GetXaxis()->SetBinLabel(ib, xlabels[ib-1]);
        fr.GetXaxis()->CenterLabels(true);
        fr.Draw();

        // Draw graphs (left→right): Core, Mid, Edge, Full
        if (g_core && g_core->GetN() > 0) g_core->Draw("P SAME");
        if (g_mid  && g_mid->GetN()  > 0) g_mid->Draw("P SAME");
        if (g_edge && g_edge->GetN() > 0) g_edge->Draw("P SAME");
        if (g_full && g_full->GetN() > 0) g_full->Draw("P SAME");

        // Legend in numerical order (Core, Mid, Edge, Full)
        TLegend L(0.62, 0.70, 0.88, 0.88); L.SetBorderSize(0); L.SetFillStyle(0);
        L.AddEntry(g_core.get(), kViews[2].pretty, "p");
        L.AddEntry(g_mid.get(),  kViews[3].pretty, "p");
        L.AddEntry(g_edge.get(), kViews[4].pretty, "p");
        L.AddEntry(g_full.get(), kViews[1].pretty, "p");
        L.Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.DrawLatex(0.14, 0.93, Form("#sigma_{#pi^{0}} (%s) across variants", ePretty(iE).Data()));

        // Save the full (all-η) version for this energy bin
        {
          const TString outP = Form("%s/pi0Sigma_Ebin_%02d.png", baseOut.Data(), iE);
          cS.Modified(); cS.Update(); cS.SaveAs(outP);
          std::cout << "[SigmaOverlay] Wrote: " << outP << "\n";
        }

        // ALSO: save a second version that skips η-core entirely (still 5 categories)
        {
          double ymin2 = +1e30, ymax2 = -1e30, X2, Y2;
          auto upd2 = [&](TGraphErrors* g) {
            if (!g) return;
            for (int k = 0; k < g->GetN(); ++k) {
              g->GetPoint(k, X2, Y2);
              const double eY = g->GetErrorY(k);
              ymin2 = std::min(ymin2, Y2 - eY);
              ymax2 = std::max(ymax2, Y2 + eY);
            }
          };
          upd2(g_edge.get()); upd2(g_mid.get()); upd2(g_full.get());
          if (!(std::isfinite(ymin2) && std::isfinite(ymax2))) { ymin2 = 0.010; ymax2 = 0.040; }
          double span2 = ymax2 - ymin2;
          if (span2 <= 0) { span2 = 1e-4; }
          const double pad2b = std::max(0.06 * span2, 5e-4);
          const double pad2t = std::max(0.10 * span2, 8e-4);
          ymin2 = std::max(0.0, ymin2 - pad2b);
          ymax2 = ymax2 + pad2t;

          TCanvas cS2(Form("cSigmaOverlay_noCore_%02d", iE), "Gaussian sigma (no etaCore)", 1000, 700);
          cS2.SetLeftMargin(0.12);
          cS2.SetBottomMargin(0.18);

          TH1F fr2("fr_sigma_noCore", "", 5, 0.5, 5.5);
          fr2.SetStats(0);
          fr2.GetYaxis()->SetTitle("Gaussian #sigma  [GeV]");
          fr2.GetYaxis()->SetRangeUser(ymin2, ymax2);
          fr2.GetXaxis()->SetTitle("");
          fr2.GetXaxis()->SetTitleSize(0.0);
          fr2.GetXaxis()->SetLabelSize(0.032);
          fr2.GetYaxis()->SetTitleOffset(1.1);
          for (int ib = 1; ib <= 5; ++ib) fr2.GetXaxis()->SetBinLabel(ib, xlabels[ib-1]);
          fr2.GetXaxis()->CenterLabels(true);
          fr2.Draw();

          if (g_mid  && g_mid->GetN()  > 0) g_mid->Draw("P SAME");
          if (g_edge && g_edge->GetN() > 0) g_edge->Draw("P SAME");
          if (g_full && g_full->GetN() > 0) g_full->Draw("P SAME");

          TLegend L2(0.62, 0.70, 0.88, 0.88); L2.SetBorderSize(0); L2.SetFillStyle(0);
          L2.AddEntry(g_mid.get(),  kViews[3].pretty, "p");
          L2.AddEntry(g_edge.get(), kViews[4].pretty, "p");
          L2.AddEntry(g_full.get(), kViews[1].pretty, "p");
          L2.Draw();

          TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.035);
          tx2.DrawLatex(0.14, 0.93, Form("#sigma_{#pi^{0}} (%s) across variants  —  no |#eta| #leq 0.20", ePretty(iE).Data()));

          const TString outP2 = Form("%s/pi0Sigma_Ebin_%02d_noCore.png", baseOut.Data(), iE);
          cS2.Modified(); cS2.Update(); cS2.SaveAs(outP2);
          std::cout << "[SigmaOverlay] Wrote: " << outP2 << "\n";
        }
      } // for iE
    }
  }

    // -------------------------- Resolution overlay (σ/μ) across variants (all E bins) --------------------------
    {
      int vRAW   = -1, vCP=-1, vEAone=-1, vEAeta=-1, vEAgeo=-1;
      for (int v = 0; v < NVAR; ++v) {
        const std::string tag = kVarTags[v];
        if      (tag == "CLUSRAW") vRAW   = v;
        else if (tag == "CLUSCP")  vCP    = v;
        else if (tag == "EAEonly") vEAone = v;
        else if (tag == "EAetaE")  vEAeta = v;
        else if (tag == "EAgeom")  vEAgeo = v;
      }
      if (vRAW < 0 || vCP < 0 || vEAone < 0 || vEAeta < 0 || vEAgeo < 0) {
        std::cerr << "[ResolutionOverlay] Missing one of required variants (CLUSRAW, CLUSCP, EAEonly, EAetaE, EAgeom) — skipping overlay.\n";
      } else if (N_E <= 0) {
        std::cerr << "[ResolutionOverlay] No energy bins available — skipping overlay.\n";
      } else {
        for (int iE = 0; iE < N_E; ++iE) {
          // Category x-positions (bin centers)
          std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};

          // Local label shortener: ONLY used to build x-axis bin labels for this block.
          auto shortLabel = [&](const char* s)->const char* {
            if (!s) return "";
            if (strcmp(s, "No Correction") == 0)                      return "Raw";
            if (strcmp(s, "Legacy Constant b-corr") == 0)             return "Const-b";
            if (strcmp(s, "Energy-Dep b-correction") == 0)            return "Energy-Dep b";
            if (strcmp(s, "Energy/|#eta|-Dep b-correction") == 0)     return "Energy/|#eta|-Dep b";
            if (strcmp(s, "Energy/|z_{vtx}|-Dep b-correction") == 0)  return "Energy/|z_{vtx}|-Dep b";
            return s; // fallback unchanged
          };

          // Build the x-axis labels with the short mapping
          const char* xlabels_short[5] = {
            shortLabel(kVarPretty[vRAW]),
            shortLabel(kVarPretty[vCP]),
            shortLabel(kVarPretty[vEAone]),
            shortLabel(kVarPretty[vEAeta]),
            shortLabel(kVarPretty[vEAgeo])
          };

          auto makeViewGraphRes = [&](int iview, double xoff, Style_t mstyle, Color_t col)->std::unique_ptr<TGraphErrors>
          {
            std::vector<double> y, ey;
            auto push = [&](const FitRes& R) {
              if (R.ok && std::isfinite(R.sig) && std::isfinite(R.mean) && R.sig>0.0 && R.mean>0.0) {
                const double res = R.sig / R.mean;
                const double err = res * std::sqrt( std::pow(R.sigErr / R.sig,  2.0) +
                                                    std::pow(R.meanErr / R.mean, 2.0) );
                y.push_back(res); ey.push_back(err);
              } else {
                y.push_back(std::numeric_limits<double>::quiet_NaN()); ey.push_back(0.0);
              }
            };

            // Order matches x: CLUSRAW, CLUSCP, EAEonly, EAetaE, EAgeom
            push(summaryV[iview][vRAW][iE]);
            push(summaryV[iview][vCP][iE]);
            push(summaryV[iview][vEAone][iE]);
            push(summaryV[iview][vEAeta][iE]);
            push(summaryV[iview][vEAgeo][iE]);

            // Filter NaNs while preserving x order; apply per-series x offset
            std::vector<double> xf, yf, exf, eyf;
            for (size_t k = 0; k < y.size(); ++k) {
              if (!std::isfinite(y[k])) continue;
              xf.push_back(x[k] + xoff);
              yf.push_back(y[k]);
              exf.push_back(0.0);
              eyf.push_back(ey[k]);
            }

            auto g = std::make_unique<TGraphErrors>(
              (int)xf.size(),
              xf.empty()?nullptr:&xf[0],
              yf.empty()?nullptr:&yf[0],
              exf.empty()?nullptr:&exf[0],
              eyf.empty()?nullptr:&eyf[0]
            );
            g->SetMarkerStyle(mstyle);
            g->SetMarkerSize(1.5);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            g->SetLineWidth(2);
            return g;
          };

          auto r_core = makeViewGraphRes(2, -0.18, 20, kBlue+1);
          auto r_mid  = makeViewGraphRes(3, -0.06, 21, kGreen+2);
          auto r_edge = makeViewGraphRes(4,  0.06, 22, kMagenta+1);
          auto r_full = makeViewGraphRes(1,  0.18, 24, kRed+1);

          // Dynamic y-range
          double rymin = +1e30, rymax = -1e30, RX, RY;
          auto rupd = [&](TGraphErrors* g) {
            if (!g) return;
            for (int k = 0; k < g->GetN(); ++k) {
              g->GetPoint(k, RX, RY);
              const double eY = g->GetErrorY(k);
              rymin = std::min(rymin, RY - eY);
              rymax = std::max(rymax, RY + eY);
            }
          };
          rupd(r_core.get()); rupd(r_mid.get()); rupd(r_edge.get()); rupd(r_full.get());
          if (!(std::isfinite(rymin) && std::isfinite(rymax))) { rymin = 0.10; rymax = 0.20; }
          double rspan = rymax - rymin; if (rspan <= 0) rspan = 1e-3;
          const double rpadB = 0.18 * rspan;
          const double rpadT = 0.5 * rspan;
          rymin = std::max(0.0, rymin - rpadB);
          rymax = rymax + rpadT;

          // Frame and draw
          TCanvas cR(Form("cResolutionOverlay_%02d", iE), "Gaussian resolution (E bin)", 1000, 700);
          cR.SetLeftMargin(0.18);
          cR.SetBottomMargin(0.18);

          TH1F frr("fr_resolution", "", 5, 0.5, 5.5);
          frr.SetStats(0);
          frr.GetYaxis()->SetTitle("Gaussian resolution  #sigma / #mu");
          frr.GetYaxis()->SetLabelSize(0.025);
          frr.GetYaxis()->SetRangeUser(rymin, rymax);

          // X-axis cosmetics + mapped labels
          frr.GetXaxis()->SetTitle("");
          frr.GetXaxis()->SetTitleSize(0.0);
          frr.GetXaxis()->SetLabelSize(0.035);
          frr.GetXaxis()->SetLabelOffset(0.010); // a bit closer to ticks
          frr.GetXaxis()->CenterLabels(true);
          frr.LabelsOption("h","X");             // force horizontal bin labels

          frr.GetYaxis()->SetTitleOffset(1.55);

          for (int ib = 1; ib <= 5; ++ib)
            frr.GetXaxis()->SetBinLabel(ib, xlabels_short[ib-1]);

          frr.Draw();

          if (r_core && r_core->GetN() > 0) r_core->Draw("P SAME");
          if (r_mid  && r_mid->GetN()  > 0) r_mid->Draw("P SAME");
          if (r_edge && r_edge->GetN() > 0) r_edge->Draw("P SAME");
          if (r_full && r_full->GetN() > 0) r_full->Draw("P SAME");

          TLegend LR(0.62, 0.70, 0.88, 0.88); LR.SetBorderSize(0); LR.SetFillStyle(0);
          LR.AddEntry(r_core.get(), kViews[2].pretty, "p");
          LR.AddEntry(r_mid.get(),  kViews[3].pretty, "p");
          LR.AddEntry(r_edge.get(), kViews[4].pretty, "p");
          LR.AddEntry(r_full.get(), kViews[1].pretty, "p");
          LR.Draw();

          TLatex txR; txR.SetNDC(); txR.SetTextSize(0.035);
          txR.DrawLatex(0.14, 0.93, Form("#sigma_{#pi^{0}} / #mu_{#pi^{0}} (%s) across variants", ePretty(iE).Data()));

          const TString outR = Form("%s/pi0Resolution_Ebin_%02d.png", baseOut.Data(), iE);
          cR.Modified(); cR.Update(); cR.SaveAs(outR);
          std::cout << "[ResolutionOverlay] Wrote: " << outR << "\n";

          // No-Core version (exclude etaCore) — still 5 categories on x
          {
            double r2min = +1e30, r2max = -1e30, XX, YY;
            auto rupd2 = [&](TGraphErrors* g) {
              if (!g) return;
              for (int k = 0; k < g->GetN(); ++k) {
                g->GetPoint(k, XX, YY);
                const double eY = g->GetErrorY(k);
                r2min = std::min(r2min, YY - eY);
                r2max = std::max(r2max, YY + eY);
              }
            };
            rupd2(r_mid.get()); rupd2(r_edge.get()); rupd2(r_full.get());
            if (!(std::isfinite(r2min) && std::isfinite(r2max))) { r2min = 0.10; r2max = 0.20; }
            double rsp2 = r2max - r2min; if (rsp2 <= 0) rsp2 = 1e-3;
            const double r2padB = 0.06 * rsp2;
            const double r2padT = 0.10 * rsp2;
            r2min = std::max(0.0, r2min - r2padB);
            r2max = r2max + r2padT;

            TCanvas cR2(Form("cResolutionOverlay_noCore_%02d", iE), "Gaussian resolution (no etaCore)", 1000, 700);
            cR2.SetLeftMargin(0.18);
            cR2.SetBottomMargin(0.18);

            TH1F frr2("fr_resolution_noCore", "", 5, 0.5, 5.5);
            frr2.SetStats(0);
            frr2.GetYaxis()->SetTitle("Gaussian resolution  #sigma / #mu");
            frr2.GetYaxis()->SetRangeUser(r2min, r2max);

            // X-axis cosmetics + mapped labels (same mapping as above)
            frr2.GetXaxis()->SetTitle("");
            frr2.GetXaxis()->SetTitleSize(0.0);
            frr2.GetXaxis()->SetLabelSize(0.032);
            frr2.GetXaxis()->SetLabelOffset(0.010);
            frr2.GetXaxis()->CenterLabels(true);
            frr2.LabelsOption("h","X");

            frr2.GetYaxis()->SetTitleOffset(1.55);

            for (int ib = 1; ib <= 5; ++ib)
              frr2.GetXaxis()->SetBinLabel(ib, xlabels_short[ib-1]);

            frr2.Draw();

            if (r_mid  && r_mid->GetN()  > 0) r_mid->Draw("P SAME");
            if (r_edge && r_edge->GetN() > 0) r_edge->Draw("P SAME");
            if (r_full && r_full->GetN() > 0) r_full->Draw("P SAME");

            TLegend LR2(0.62, 0.70, 0.88, 0.88); LR2.SetBorderSize(0); LR2.SetFillStyle(0);
            LR2.AddEntry(r_mid.get(),  kViews[3].pretty, "p");
            LR2.AddEntry(r_edge.get(), kViews[4].pretty, "p");
            LR2.AddEntry(r_full.get(), kViews[1].pretty, "p");
            LR2.Draw();

            TLatex txR2; txR2.SetNDC(); txR2.SetTextSize(0.035);
            txR2.DrawLatex(0.14, 0.93, Form("#sigma_{#pi^{0}} / #mu_{#pi^{0}} (%s) across variants  —  no |#eta| #leq 0.20", ePretty(iE).Data()));

              const TString outR = Form("%s/pi0Resolution_Ebin_%02d.png", baseOut.Data(), iE);
              cR.Modified(); cR.Update(); cR.SaveAs(outR);
              std::cout << "[ResolutionOverlay] Wrote: " << outR << "\n";

              /* ---------- ADD: well-formatted summary table (per η-view) ---------- */
              auto shortVarName = [&](int v)->const char* {
                if      (v == vRAW)   return "Raw";
                else if (v == vCP)    return "Const-b";
                else if (v == vEAone) return "Energy-Dep b";
                else if (v == vEAeta) return "Energy/|#eta|-Dep b";
                else if (v == vEAgeo) return "Energy/|z_{vtx}|-Dep b";
                return kVarPretty[v];
              };
              auto resOf = [&](int iview, int v)->double {
                const FitRes& R = summaryV[iview][v][iE];
                if (!R.ok || !std::isfinite(R.sig) || !std::isfinite(R.mean) || R.sig<=0.0 || R.mean<=0.0)
                  return std::numeric_limits<double>::quiet_NaN();
                return R.sig / R.mean;
              };
              auto nEventsOf = [&](int iview, int v)->Long64_t {
                if (iview < 0 || iview >= N_VIEW) return 0;
                if (v < 0 || v >= NVAR) return 0;
                if (!hKeepV[iview][v][iE]) return 0;
                return static_cast<Long64_t>( hKeepV[iview][v][iE]->GetEntries() );
              };
              auto fmtD = [](double x, int p=5)->std::string{
                if (!std::isfinite(x)) return std::string("   n/a   ");
                std::ostringstream os; os.setf(std::ios::fixed); os<<std::setprecision(p)<<x; return os.str();
              };
              auto fmtPct = [](double x)->std::string{
                if (!std::isfinite(x)) return std::string("   n/a   ");
                std::ostringstream os; os.setf(std::ios::fixed); os<<std::showpos<<std::setprecision(2)<<x<<"%"<<std::noshowpos; return os.str();
              };

              const int W_ETA = 18, W_VAR = 20, W_N = 10, W_RES = 12, W_D = 12;

              // print one table per η-view (in the same order as the plot: Core, Mid, Edge, Full)
              std::vector<int> ivList = {2, 3, 4, 1}; // etaCore, etaMid, etaEdge, fullEta

              std::cout << "\n\033[1m[ResolutionOverlay] Summary (σ/μ)   E-bin = ["
                        << std::fixed << std::setprecision(0) << E_edges[iE] << "," << E_edges[iE+1]
                        << ") GeV\033[0m\n";
              std::cout << std::left << std::setw(W_ETA) << "η-range"
                        << " | " << std::setw(W_VAR) << "Variant"
                        << " | " << std::setw(W_N)   << "nEvents"
                        << " | " << std::setw(W_RES) << "RES(σ/μ)"
                        << " | " << std::setw(W_D)   << "Δ vs RAW"
                        << "\n";
              std::cout << std::string(W_ETA+3 + W_VAR+3 + W_N+3 + W_RES+3 + W_D, '-') << "\n";

              for (int iview : ivList) {
                // precompute RAW for this η view
                const double resRAW = resOf(iview, vRAW);

                // print a block header line for the η view
                std::cout << std::left << std::setw(W_ETA) << kViews[iview].pretty
                          << " | " << std::setw(W_VAR) << shortVarName(vRAW)
                          << " | " << std::setw(W_N)   << nEventsOf(iview, vRAW)
                          << " | " << std::setw(W_RES) << fmtD(resRAW, 5)
                          << " | " << std::setw(W_D)   << fmtPct(0.0)
                          << "\n";

                // the rest of the variants in the fixed order
                const int vList[4] = {vCP, vEAone, vEAeta, vEAgeo};
                for (int vv : vList) {
                  const double r = resOf(iview, vv);
                  const double d = (std::isfinite(resRAW) && resRAW>0.0 && std::isfinite(r))
                                   ? ( (r/resRAW - 1.0) * 100.0 ) : std::numeric_limits<double>::quiet_NaN();

                  std::cout << std::left << std::setw(W_ETA) << ""
                            << " | " << std::setw(W_VAR) << shortVarName(vv)
                            << " | " << std::setw(W_N)   << nEventsOf(iview, vv)
                            << " | " << std::setw(W_RES) << fmtD(r, 5)
                            << " | " << std::setw(W_D)   << fmtPct(d)
                            << "\n";
                }
                // separator between η-views
                std::cout << std::string(W_ETA+3 + W_VAR+3 + W_N+3 + W_RES+3 + W_D, '-') << "\n";
              }
              /* ---------- END add ---------- */

          }
        } // iE
      }
    }

  // -------------------------- Sigma RATIO vs E across η-views (money plots) --------------------------
  {
    // Find indices for RAW baseline and requested variants
    int vRAW = -1, vCP = -1, vEAone = -1, vEAeta = -1, vEAgeo = -1;
    for (int v = 0; v < NVAR; ++v) {
      const std::string tag = kVarTags[v];
      if      (tag == "CLUSRAW") vRAW   = v;
      else if (tag == "CLUSCP")  vCP    = v;
      else if (tag == "EAEonly") vEAone = v;
      else if (tag == "EAetaE")  vEAeta = v;
      else if (tag == "EAgeom")  vEAgeo = v;
    }

    if (vRAW < 0) {
      std::cerr << "[SigmaRatio] CLUSRAW index not found — skipping ratio plots.\n";
    } else {
      gSystem->mkdir(baseOut + "/SigmaRatios", true);

      // Build sigma(variant)/sigma(RAW) vs E with proper uncertainty propagation
      auto makeRatioGraph = [&](int iview, int vVar, Color_t col, Style_t mstyle)->std::unique_ptr<TGraphErrors>
      {
        std::vector<double> x, y, ex, ey;
        x.reserve(N_E); y.reserve(N_E); ex.assign(N_E, 0.0); ey.reserve(N_E);

        for (int i = 0; i < N_E; ++i) {
          const FitRes& Rs = summaryV[iview][vVar][i]; // variant
          const FitRes& Rr = summaryV[iview][vRAW][i]; // RAW baseline
          if (!Rs.ok || !Rr.ok) continue;
          if (!(std::isfinite(Rs.sig) && std::isfinite(Rr.sig))) continue;
          if (Rs.sig <= 0.0 || Rr.sig <= 0.0) continue;

          const double ratio = Rs.sig / Rr.sig;
          const double err   = ratio * std::sqrt( std::pow(Rs.sigErr / Rs.sig, 2.0) +
                                                  std::pow(Rr.sigErr / Rr.sig, 2.0) );
          x.push_back( eCtrOf(i) );
          y.push_back( ratio );
          ey.push_back( err );
        }

        auto g = std::make_unique<TGraphErrors>(
          (int)x.size(),
          x.empty()?nullptr:&x[0],
          y.empty()?nullptr:&y[0],
          ex.empty()?nullptr:&ex[0],
          ey.empty()?nullptr:&ey[0]
        );
        g->SetMarkerStyle(mstyle);
        g->SetMarkerSize(1.2);
        g->SetMarkerColor(col);
        g->SetLineColor(col);
        g->SetLineWidth(2);
        return g;
      };

      struct ViewStyle { int iv; Color_t col; Style_t ms; };
      // Order: Core, Mid, Edge, Full (η-dependent views only)
      std::vector<ViewStyle> views = {
        {2, kBlue+1,   20},  // etaCore
        {3, kGreen+2,  21},  // etaMid
        {4, kMagenta+1,22},  // etaEdge
        {1, kRed+1,    24}   // fullEta
      };

      auto moneyPlot = [&](int vVar, const char* varPretty, const char* varTag)
      {
        std::vector<std::unique_ptr<TGraphErrors>> gs;
        double ymin = +1e30, ymax = -1e30, X, Y;
        int Ntot = 0;

        for (const auto& vw : views) {
          auto g = makeRatioGraph(vw.iv, vVar, vw.col, vw.ms);
          Ntot += g->GetN();
          if (g->GetN() > 0) {
            for (int k = 0; k < g->GetN(); ++k) {
              g->GetPoint(k, X, Y);
              const double eY = g->GetErrorY(k);
              ymin = std::min(ymin, Y - eY);
              ymax = std::max(ymax, Y + eY);
            }
          }
          gs.emplace_back(std::move(g));
        }

        if (Ntot == 0) {
          std::cerr << "[SigmaRatio] No valid points for variant " << varTag << " — skipping.\n";
          return;
        }

        if (!(std::isfinite(ymin) && std::isfinite(ymax))) { ymin = 0.8; ymax = 1.2; }
        double span = ymax - ymin; if (span <= 0) span = 0.1;
        const double pad = 0.15 * span;
        ymin = std::max(0.0, ymin - pad);
        ymax = ymax + pad;

        TCanvas c("cRatio", "sigma ratio vs E", 950, 700);
        c.SetLeftMargin(0.14);
        c.SetBottomMargin(0.14);

        TH2F fr("fr_ratio", "", 100, E_edges[0]-0.5, E_edges[N_E]-0.5, 100, ymin, ymax);
        fr.SetTitle(Form("#sigma_{#pi^{0}} ratio: %s / No Correction;E [GeV];#sigma_{%s}/#sigma_{RAW}", varPretty, "#pi^{0}"));
        fr.SetStats(0);
        fr.Draw();

        // Reference line at ratio = 1
        TLine Lref(E_edges[0]-0.5, 1.0, E_edges[N_E]-0.5, 1.0);
        Lref.SetLineColor(kGray+2);
        Lref.SetLineStyle(2);
        Lref.Draw("SAME");

        for (const auto& g : gs) { if (g && g->GetN() > 0) g->Draw("P SAME"); }

        TLegend LR(0.62, 0.70, 0.88, 0.88); LR.SetBorderSize(0); LR.SetFillStyle(0);
        LR.AddEntry(gs[0].get(), kViews[2].pretty, "p"); // etaCore
        LR.AddEntry(gs[1].get(), kViews[3].pretty, "p"); // etaMid
        LR.AddEntry(gs[2].get(), kViews[4].pretty, "p"); // etaEdge
        LR.AddEntry(gs[3].get(), kViews[1].pretty, "p"); // fullEta
        LR.Draw();

        const TString outP = baseOut + "/SigmaRatios/sigmaRatio_" + TString(varTag) + "_over_CLUSRAW_vsE.png";
        c.Modified(); c.Update(); c.SaveAs(outP);
        std::cout << "[SigmaRatio] Wrote: " << outP << "\n";
      };

      // Produce the requested money plots:
      if (vCP    >= 0) moneyPlot(vCP,    "Legacy Constant b-corr",          "CLUSCP");
      if (vEAone >= 0) moneyPlot(vEAone, "Energy-Dep b-correction",         "EAEonly");
      if (vEAeta >= 0) moneyPlot(vEAeta, "Energy/|#eta|-Dep b-correction",  "EAetaE");
      if (vEAgeo >= 0) moneyPlot(vEAgeo, "Energy/|z_{vtx}|-Dep b-correction","EAgeom");

      // -------------------------- TERMINAL SUMMARY TABLE --------------------------
      // Print per-η-view tables of RMS percent differences vs CLUSRAW
      // Δ% = 100 * (σ_variant - σ_RAW) / σ_RAW   (one row per energy bin)
      auto pctStr = [](double num)->std::string {
        if (!std::isfinite(num)) return std::string("   n/a ");
        std::ostringstream os; os.setf(std::ios::fixed); os<<std::setprecision(2)<<std::setw(7)<<num;
        return os.str();
      };

      struct VarDef { int idx; const char* tag; };
      std::vector<VarDef> vars;
      if (vCP    >= 0) vars.push_back({vCP,    "CLUSCP "});
      if (vEAone >= 0) vars.push_back({vEAone, "EAEonly"});
      if (vEAeta >= 0) vars.push_back({vEAeta, "EAetaE "});
      if (vEAgeo >= 0) vars.push_back({vEAgeo, "EAgeom "});  // NEW in summary table

      // Views in numerical order: Core, Mid, Edge, Full
      std::vector<int> viewOrder = {2, 3, 4, 1};

      std::cout << "\n";
      std::cout << "╔═════════════════════════════════════════════════════════════════════════════════════╗\n";
      std::cout << "║   RMS percent differences vs CLUSRAW  (per η view, per E-bin)                      ║\n";
      std::cout << "║   Δ% = 100 × (σ_variant − σ_RAW) / σ_RAW                                           ║\n";
      std::cout << "╚═════════════════════════════════════════════════════════════════════════════════════╝\n";

      for (int iv : viewOrder)
      {
        if (vRAW < 0 || vars.empty()) continue;

        std::cout << "\n";
        std::cout << "View: " << kViews[iv].pretty << "   (key=" << kViews[iv].key << ")\n";
        std::cout << "─────────────────────────────────────────────────────────────────────────────────────\n";
        // Header
        std::cout << std::left << std::setw(16) << "E-bin"
                  << std::right;
        for (const auto& vd : vars) {
          std::cout << " | " << std::setw(12) << vd.tag;
        }
        std::cout << "\n";

        std::cout << std::string(16 + (int)vars.size() * (3 + 12), '-') << "\n";

        // Row per energy bin
        for (int i = 0; i < N_E; ++i)
        {
          const FitRes& Rr = summaryV[iv][vRAW][i];
          std::cout << std::left << std::setw(16)
                    << Form("[%.0f,%.0f) GeV", E_edges[i], E_edges[i+1])
                    << std::right;

          for (const auto& vd : vars)
          {
            const FitRes& Rs = summaryV[iv][vd.idx][i];
            double pct = std::numeric_limits<double>::quiet_NaN();
            if (Rr.ok && Rs.ok && std::isfinite(Rr.sig) && std::isfinite(Rs.sig) && Rr.sig > 0.0)
              pct = 100.0 * (Rs.sig - Rr.sig) / Rr.sig;

            std::cout << " | " << std::setw(12) << pctStr(pct);
          }
          std::cout << "\n";
        }

        // Per-view averages (unweighted over bins with valid entries)
        std::cout << std::string(16 + (int)vars.size() * (3 + 12), '-') << "\n";
        std::cout << std::left << std::setw(16) << "Avg Δ% (bins)"
                  << std::right;
        for (const auto& vd : vars)
        {
          double sum = 0.0; int n=0;
          for (int i = 0; i < N_E; ++i) {
            const FitRes& Rr = summaryV[iv][vRAW][i];
            const FitRes& Rs = summaryV[iv][vd.idx][i];
            if (Rr.ok && Rs.ok && std::isfinite(Rr.sig) && std::isfinite(Rs.sig) && Rr.sig > 0.0) {
              sum += 100.0 * (Rs.sig - Rr.sig) / Rr.sig;
              ++n;
            }
          }
          const double avg = (n>0) ? (sum / n) : std::numeric_limits<double>::quiet_NaN();
          std::cout << " | " << std::setw(12) << pctStr(avg);
        }
        std::cout << "\n";
      }
      std::cout << std::endl;
      // --------------------------------------------------------------------
    }
  }

    // -------------------------- NEW: resolutionRatios (originalEta only, 8/10/15 GeV variants) --------------------------
    {
      // Variant indices
      int vRAW=-1, vCP=-1, vEAgeo=-1, vEAeta=-1, vEAone=-1;
      for (int v=0; v<NVAR; ++v) {
        std::string t = kVarTags[v];
        if      (t=="CLUSRAW") vRAW  = v;
        else if (t=="CLUSCP")  vCP   = v;
        else if (t=="EAgeom")  vEAgeo= v; // |z_vtx| dep
        else if (t=="EAetaE")  vEAeta= v; // |eta|+E dep
        else if (t=="EAEonly") vEAone= v; // E-only
      }

      auto makeOneResolutionRatioPlot = [&](double E_hi_cap, const char* suffix)
      {
        // proceed only if all 5 are available
        if (!(vRAW>=0 && vCP>=0 && vEAgeo>=0 && vEAeta>=0 && vEAone>=0)) return;

        const int iv = 0; // originalEta

        // Determine max energy bin to include:
        // prefer exact [E_hi_cap-3, E_hi_cap) if present; else all bins with E_hi <= E_hi_cap
        int iMax = -1, iExact = -1;
        for (int i=0; i<N_E; ++i) {
          if (std::fabs(E_edges[i+1] - E_hi_cap) < 1e-6) iExact = i;                // exact terminal bin ends at cap
          if (E_edges[i+1] <= E_hi_cap + 1e-9)   iMax   = i;                        // <= cap
        }
        if (iExact >= 0) iMax = iExact;
        if (iMax < 0) return; // nothing within cap

        gSystem->mkdir(baseOut + "/resolutionRatios", true);

        // Per-series X offsets inside each E-bin (fraction of bin width):
        // RAW << CP < Center (E-only) < EA(|eta|+E) << EA(|z_vtx|+E)
        auto xOffsetFrac = [&](int vIdx)->double {
          if      (vIdx==vRAW)   return -0.16;
          else if (vIdx==vCP)    return -0.08;
          else if (vIdx==vEAone) return  0.00;
          else if (vIdx==vEAeta) return  0.08;
          else if (vIdx==vEAgeo) return  0.16;
          return 0.0;
        };

        // Build Gaussian resolution points (sigma/mu) for a variant, with per-bin X offsets
        auto makeResGraph = [&](int vIdx, Color_t col, Style_t mk)->std::unique_ptr<TGraphErrors>{
          std::vector<double> x, y, ex, ey;
          x.reserve(iMax+1); y.reserve(iMax+1); ex.clear(); ey.reserve(iMax+1);
          for (int i=0; i<=iMax; ++i) {
            const auto &R = summaryV[iv][vIdx][i];
            if (R.ok && std::isfinite(R.sig) && std::isfinite(R.mean) && R.sig>0.0 && R.mean>0.0) {
              const double res  = R.sig / R.mean;
              const double dres = res * std::sqrt( std::pow(R.sigErr/std::max(R.sig,1e-12),2.0) +
                                                   std::pow(R.meanErr/std::max(R.mean,1e-12),2.0) );
              const double eLo  = E_edges[i];
              const double eHi  = E_edges[i+1];
              const double eCtr = 0.5*(eLo + eHi);
              const double eWid = (eHi - eLo);
              const double xOff = xOffsetFrac(vIdx) * eWid;
              x.push_back( eCtr + xOff );
              y.push_back( res );
              ey.push_back( dres );
            }
          }
          ex.assign(x.size(), 0.0);
          auto g = std::make_unique<TGraphErrors>((int)x.size(),
                             x.empty()?nullptr:&x[0],
                             y.empty()?nullptr:&y[0],
                             ex.empty()?nullptr:&ex[0],
                             ey.empty()?nullptr:&ey[0]);
          g->SetMarkerStyle(mk); g->SetMarkerSize(1.2);
          g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
          return g;
        };

        auto gRAW   = makeResGraph(vRAW  , kVarColor[vRAW  ], 20);  // No Correction
        auto gCP    = makeResGraph(vCP   , kVarColor[vCP   ], 20);  // Legacy Constant b-corr
        auto gZDEP  = makeResGraph(vEAgeo, kVarColor[vEAgeo], 20);  // Energy/|z_{vtx}|-Dep
        auto gETADE = makeResGraph(vEAeta, kVarColor[vEAeta], 20);  // Energy/|#eta|-Dep
        auto gEONLY = makeResGraph(vEAone, kVarColor[vEAone], 20);  // Energy-Dep

        // Only draw if we have RAW and at least one corrected curve
        if (!(gRAW && (gCP || gZDEP || gETADE || gEONLY))) return;

        // X & Y ranges (use bin edges up to iMax; show last edge at cap)
        const double xl = E_edges[0];
        const double xr = E_edges[iMax+1];
        const double xr_tick = xr + 1e-6; // ensure terminal tick/label at edge is drawn

        auto scanY = [&](const TGraphErrors* G, double& lo, double& hi){
          if (!G) return;
          for (int k=0;k<G->GetN();++k){ double X,Y; G->GetPoint(k,X,Y); const double e=G->GetErrorY(k);
            lo = std::min(lo, Y - e); hi = std::max(hi, Y + e); }
        };
        double ylo=+1e30, yhi=-1e30;
        scanY(gRAW.get(),   ylo,yhi);
        scanY(gCP.get(),    ylo,yhi);
        scanY(gZDEP.get(),  ylo,yhi);
        scanY(gETADE.get(), ylo,yhi);
        scanY(gEONLY.get(), ylo,yhi);
        if (!(ylo<yhi)) { ylo=0.05; yhi=0.35; }
        // Asymmetric padding: modest lower pad, larger upper headroom
        const double span      = (yhi - ylo);
        const double padBottom = 0.10 * span;  // keep lower buffer small
        const double padTop    = 0.35 * span;  // increase upper buffer
        ylo = std::max(0.0, ylo - padBottom);
        yhi = yhi + padTop;

        // Canvas & pads  (reduced inter-pad gap; bottom pad slightly taller)
        TCanvas c("cResTopBot","Resolution & Ratios (originalEta)",1000,850);
        TPad topP("topP","",0,0.36,1,1);  topP.SetBottomMargin(0.008); topP.SetLeftMargin(0.15); topP.SetRightMargin(0.06); topP.Draw();
        TPad botP("botP","",0,0.00,1,0.36); botP.SetTopMargin(0.012);  botP.SetLeftMargin(0.15); botP.SetRightMargin(0.06); botP.SetBottomMargin(0.20); botP.Draw();

        // ---------------- TOP panel ----------------
        topP.cd();
        TH2F frU("frU","",100,xl,xr_tick,100,ylo,yhi);
        // No x-axis labels/ticks/title on the TOP panel
        frU.SetTitle("; ;RES = #sigma_{Gaus}/#mu_{Gaus}");
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        // Encourage integer ticks including terminal edge
          // Ticks: for the 8 GeV cap, force integer ticks (2,3,...,8); else keep optimized
          if (std::fabs(E_hi_cap - 8.0) < 1e-6) {
            // Range [2,8] → 6 primary intervals → labels at each integer
            frU.GetXaxis()->SetNdivisions(6, kFALSE);
          } else {
            frU.GetXaxis()->SetNdivisions(515, kTRUE);
          }

        // Make TOP-panel Y title larger and a bit closer
        frU.GetYaxis()->SetTitleSize(0.039);
        frU.GetYaxis()->SetTitleOffset(1.15);
        frU.Draw();

        if (gRAW)   gRAW->Draw("P SAME");
        if (gCP)    gCP->Draw("P SAME");
        if (gZDEP)  gZDEP->Draw("P SAME");
        if (gETADE) gETADE->Draw("P SAME");
        if (gEONLY) gEONLY->Draw("P SAME");

        // Requested top-left, two-column legend with the exact ordering/layout
        TLegend legL(0.16,0.74,0.46,0.90); legL.SetBorderSize(0); legL.SetFillStyle(0); legL.SetTextSize(0.037);
        TLegend legR(0.48,0.74,0.92,0.90); legR.SetBorderSize(0); legR.SetFillStyle(0); legR.SetTextSize(0.037);
        if (gRAW)   legL.AddEntry(gRAW.get(),   kVarPretty[vRAW]  , "p"); // No Correction
        if (gCP)    legL.AddEntry(gCP.get(),    kVarPretty[vCP]   , "p"); // Legacy Constant b-corr
        if (gEONLY) legR.AddEntry(gEONLY.get(), kVarPretty[vEAone], "p"); // Energy-Dep b-correction
        if (gETADE) legR.AddEntry(gETADE.get(), kVarPretty[vEAeta], "p"); // Energy/|#eta|-Dep b-corr
        if (gZDEP)  legR.AddEntry(gZDEP.get(),  kVarPretty[vEAgeo], "p"); // Energy/|z_vtx|-Dep b-corr
        legL.Draw(); legR.Draw();

        // ---------------- BOTTOM panel ----------------
        botP.cd();

        auto makeResAt = [&](int vIdx, int i, double& res, double& err)->bool{
          const auto &R = summaryV[iv][vIdx][i];
          if (R.ok && std::isfinite(R.sig) && std::isfinite(R.mean) && R.sig>0.0 && R.mean>0.0) {
            res  = R.sig / R.mean;
            err  = res * std::sqrt( std::pow(R.sigErr/std::max(R.sig,1e-12),2.0) +
                                    std::pow(R.meanErr/std::max(R.mean,1e-12),2.0) );
            return true;
          }
          return false;
        };

        auto makeRatioGraphToRAW = [&](int vIdx, Color_t col, Style_t mk)->std::unique_ptr<TGraphErrors>{
          // Bottom-panel bin offsets: center intentionally left empty.
          auto bottomOffsetFrac = [&](int vv)->double {
            if      (vv == vCP)    return -0.12;
            else if (vv == vEAone) return -0.04;
            else if (vv == vEAeta) return  0.04;
            else if (vv == vEAgeo) return  0.12;
            return 0.0;
          };

          std::vector<double> x,y,ex,ey;
          for (int i=0; i<=iMax; ++i) {
            double rv, dv, rr, dr;
            if (makeResAt(vIdx, i, rv, dv) && makeResAt(vRAW, i, rr, dr) && rv>0.0 && rr>0.0) {
              const double ratio = rv / rr;
              const double erel2 = std::pow(dv/std::max(rv,1e-12),2.0) + std::pow(dr/std::max(rr,1e-12),2.0);
              const double eLo   = E_edges[i], eHi = E_edges[i+1];
              const double eCtr  = 0.5*(eLo + eHi);
              const double eWid  = (eHi - eLo);
              const double xOff  = bottomOffsetFrac(vIdx) * eWid;
              x.push_back(eCtr + xOff);
              y.push_back(ratio);
              ey.push_back(ratio * std::sqrt(erel2));
            }
          }
          ex.assign(x.size(), 0.0);
          auto g = std::make_unique<TGraphErrors>((int)x.size(),
                         x.empty()?nullptr:&x[0],
                         y.empty()?nullptr:&y[0],
                         ex.empty()?nullptr:&ex[0],
                         ey.empty()?nullptr:&ey[0]);
          g->SetMarkerStyle(mk); g->SetMarkerSize(1.2);
          g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
          return g;
        };

        // Bottom ratio panel — all solid circles (style 20)
        auto rCP    = makeRatioGraphToRAW(vCP   , kVarColor[vCP]   , 20);
        auto rZDEP  = makeRatioGraphToRAW(vEAgeo, kVarColor[vEAgeo], 20);
        auto rETADE = makeRatioGraphToRAW(vEAeta, kVarColor[vEAeta], 20);
        auto rEONLY = makeRatioGraphToRAW(vEAone, kVarColor[vEAone], 20);

        double rlo=+1e30, rhi=-1e30;
        auto scanRY = [&](const TGraphErrors* G){
          if(!G) return;
          for(int k=0;k<G->GetN();++k){ double X,Y; G->GetPoint(k,X,Y);
            const double e=G->GetErrorY(k); rlo=std::min(rlo,Y-e); rhi=std::max(rhi,Y+e); }
        };
        scanRY(rCP.get()); scanRY(rZDEP.get()); scanRY(rETADE.get()); scanRY(rEONLY.get());
        if (!(rlo<rhi)) { rlo=0.9; rhi=1.1; }
        const double rpad = 0.06*(rhi-rlo);
        rlo -= rpad; rhi += rpad;

        TH2F frL("frL","",100,xl,xr_tick,100,rlo,rhi);
        frL.SetTitle(";E  [GeV];RES(Corr)/RES(Raw)");
        // Y-axis cosmetics
        frL.GetYaxis()->SetTitleSize(0.065);
        frL.GetYaxis()->SetTitleOffset(0.7);
        frL.GetYaxis()->SetLabelSize(0.050);
        // X-axis cosmetics
        frL.GetXaxis()->SetLabelSize(0.055);
        frL.GetXaxis()->SetTitleSize(0.052);
        frL.GetXaxis()->SetTitleOffset(1.05);
        frL.GetXaxis()->SetTickLength(0.045);
          if (std::fabs(E_hi_cap - 8.0) < 1e-6) {
            frL.GetXaxis()->SetNdivisions(6, kFALSE);  // integer ticks: 2..8
          } else {
            frL.GetXaxis()->SetNdivisions(515, kTRUE);
          }

        frL.Draw();

        TLine l1(xl,1.0,xr,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

        if (rCP)    rCP->Draw("P SAME");
        if (rETADE) rETADE->Draw("P SAME");
        if (rZDEP)  rZDEP->Draw("P SAME");
        if (rEONLY) rEONLY->Draw("P SAME");

        const TString outPNG = baseOut + "/resolutionRatios/resolutionRatios_originalEta_upto"
                               + TString::Format("%.0fGeV", E_hi_cap) + (suffix?TString(suffix):"") + ".png";
        c.Modified(); c.Update(); c.SaveAs(outPNG);

        // ---------- TEXT TABLES (TOP: absolute RES; BOTTOM: ratios to RAW) ----------
        auto fmt   = [](double v, int p=5)->std::string {
          if (!std::isfinite(v)) return std::string("   n/a   ");
          std::ostringstream os; os.setf(std::ios::fixed); os << std::setprecision(p) << v;
          return os.str();
        };
        auto pct   = [](double v)->std::string {
          if (!std::isfinite(v)) return std::string("    n/a    ");
          std::ostringstream os; os.setf(std::ios::fixed); os << std::setprecision(2) << v*100.0 << "%";
          return os.str();
        };
        auto safe_res_tbl = [&](int vIdx, int i)->double {
          const auto &R = summaryV[iv][vIdx][i];
          if (!R.ok || !std::isfinite(R.sig) || !std::isfinite(R.mean) || R.sig<=0.0 || R.mean<=0.0)
            return std::numeric_limits<double>::quiet_NaN();
          return R.sig / R.mean;
        };

        std::cout << "\n\033[1m[resolutionRatios] Tables (view=" << kViews[iv].pretty
                  << ", E up to " << E_edges[iMax+1] << " GeV)\033[0m\n";

        // TOP: absolute RES
        {
          std::cout << "\n\033[1mTOP (absolute RES = σ/μ)\033[0m\n";
          std::cout << std::left << std::setw(12) << "E-bin"
                    << " | " << std::setw(10) << "RAW"
                    << " | " << std::setw(10) << "Legacy"
                    << " | " << std::setw(10) << "E-only"
                    << " | " << std::setw(10) << "|η|-dep"
                    << " | " << std::setw(10) << "|z_vtx|-dep"
                    << "\n";
          std::cout << std::string(12+3+10+3+10+3+10+3+10+3+12, '-') << "\n";

          for (int i = 0; i <= iMax; ++i) {
            const double Rraw  = safe_res_tbl(vRAW  , i);
            const double Rcp   = safe_res_tbl(vCP   , i);
            const double Reo   = safe_res_tbl(vEAone, i);
            const double Reta  = safe_res_tbl(vEAeta, i);
            const double Rzdep = safe_res_tbl(vEAgeo, i);

            std::ostringstream ebin; ebin.setf(std::ios::fixed);
            ebin << std::setprecision(0) << E_edges[i] << "-" << E_edges[i+1];

            std::cout << std::left << std::setw(12) << ebin.str()
                      << " | " << std::setw(10) << fmt(Rraw)
                      << " | " << std::setw(10) << fmt(Rcp)
                      << " | " << std::setw(10) << fmt(Reo)
                      << " | " << std::setw(10) << fmt(Reta)
                      << " | " << std::setw(10) << fmt(Rzdep)
                      << "\n";
          }
        }

        // BOTTOM: ΔRES vs RAW
        {
          constexpr int W_EBIN  = 10;
          constexpr int W_COL   = 12;
          constexpr int W_WIN   = 22;

          auto colorPadPct = [&](double val, int width)->std::string {
            if (!std::isfinite(val)) {
              const std::string raw = "   n/a   ";
              return raw + std::string(std::max(0, width - (int)raw.size()), ' ');
            }
            const char* col = (val < 0.0 ? "\033[32m" : (val > 0.0 ? "\033[31m" : "\033[37m"));
            std::ostringstream os; os.setf(std::ios::fixed);
            os << std::showpos << std::setprecision(2) << val << "%" << std::noshowpos;
            const std::string raw = os.str();
            const int pad = std::max(0, width - (int)raw.size());
            return std::string(col) + raw + "\033[0m" + std::string(pad, ' ');
          };
          auto pctVal = [&](double Rraw, double Rv)->double {
            if (!std::isfinite(Rraw) || Rraw<=0.0 || !std::isfinite(Rv)) return std::numeric_limits<double>::quiet_NaN();
            return (Rv / Rraw - 1.0) * 100.0;
          };
          auto padTo = [&](std::string s, int w)->std::string {
            if ((int)s.size() < w) s += std::string(w - (int)s.size(), ' ');
            else if ((int)s.size() > w) s = s.substr(0, w);
            return s;
          };

          std::cout << "\n\033[1mBOTTOM — ΔRES vs RAW (percent; negative = better)\033[0m\n";
          std::cout << std::left << std::setw(W_EBIN) << "E-bin"
                    << " | " << std::setw(W_COL) << "ΔLegacy"
                    << " | " << std::setw(W_COL) << "ΔE-only"
                    << " | " << std::setw(W_COL) << "Δ|η|-dep"
                    << " | " << std::setw(W_COL) << "Δ|z_vtx|-dep"
                    << " | " << std::setw(W_WIN) << "Winner"
                    << "\n";
          std::cout << std::string(W_EBIN+3 + W_COL+3 + W_COL+3 + W_COL+3 + W_COL+3 + W_WIN, '-') << "\n";

          for (int i = 0; i <= iMax; ++i) {
            const double Rraw  = safe_res_tbl(vRAW  , i);
            const double Rcp   = safe_res_tbl(vCP   , i);
            const double Reo   = safe_res_tbl(vEAone, i);
            const double Reta  = safe_res_tbl(vEAeta, i);
            const double Rzdep = safe_res_tbl(vEAgeo, i);

            const double d_cp   = pctVal(Rraw, Rcp);
            const double d_eo   = pctVal(Rraw, Reo);
            const double d_eta  = pctVal(Rraw, Reta);
            const double d_zdep = pctVal(Rraw, Rzdep);

            struct Entry { const char* label; double val; };
            std::vector<Entry> cand = {
              {"Legacy",      d_cp  },
              {"E-only",      d_eo  },
              {"|η|-dep",     d_eta },
              {"|z_vtx|-dep", d_zdep}
            };
            cand.erase(std::remove_if(cand.begin(), cand.end(),
                       [](const Entry& e){ return !std::isfinite(e.val); }), cand.end());

            std::string winner = "n/a";
            if (!cand.empty()){
              const auto best = *std::min_element(cand.begin(), cand.end(),
                                [](const Entry& a, const Entry& b){ return a.val < b.val; });
              std::ostringstream w; w.setf(std::ios::fixed);
              w << best.label << " (" << std::showpos << std::setprecision(2) << best.val << "%" << std::noshowpos << ")";
              winner = w.str();
            }
            winner = padTo(winner, W_WIN);

            std::ostringstream ebin; ebin.setf(std::ios::fixed);
            ebin << std::setprecision(0) << E_edges[i] << "-" << E_edges[i+1];

            std::cout << std::left << std::setw(W_EBIN) << ebin.str()
                      << " | " << colorPadPct(d_cp  , W_COL)
                      << " | " << colorPadPct(d_eo  , W_COL)
                      << " | " << colorPadPct(d_eta , W_COL)
                      << " | " << colorPadPct(d_zdep, W_COL)
                      << " | " << winner
                      << "\n";
          }
          std::cout << std::endl;
        }
        // ---------- END TEXT TABLES ----------
      }; // end lambda

      // Generate the three parallel plots: up to 15 GeV (current), up to 10 GeV, up to 8 GeV
      makeOneResolutionRatioPlot(15.0, "");                 // existing (2..15)
      makeOneResolutionRatioPlot(10.0, "_upto10GeV");       // new (2..10)
      makeOneResolutionRatioPlot(8.0 , "_upto8GeV");        // new (2..8)
    }


    std::cout << "[Pi0MassAna] DONE. Outputs → " << baseOut << "\n";
}





static void SaveTH3Lego(TH3* h3, const char* outPNG, const char* etaLabelPretty)
{
  if (!h3) return;
  gStyle->SetOptStat(0);

  TCanvas c("c_lego","",1000,800);
  c.cd();
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin (0.12);
  gPad->SetBottomMargin(0.12);

  // Match the “good” plot: show only the high-E window (20-30 GeV)
  const double zMin = 20.0 + 1e-6;
  const double zMax = 30.0 - 1e-6;
  const int    zLo  = h3->GetZaxis()->FindBin(zMin);
  const int    zHi  = h3->GetZaxis()->FindBin(zMax);
  h3->GetZaxis()->SetRange(zLo, zHi);

  // Axis titles to match previous output
  h3->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
  h3->GetYaxis()->SetTitle("block #varphi_{local, 2#times 2}");
  h3->GetZaxis()->SetTitle("Cluster E  [GeV]");

  // Keep original title and append the η-range tag
  TString baseTitle = h3->GetTitle();
  h3->SetTitle(Form("%s  (#scale[0.8]{%s})", baseTitle.Data(), etaLabelPretty));

  // Draw with the same view settings as before
  h3->Draw("LEGO2 Z");
  gPad->SetTheta(26);
  gPad->SetPhi(38);

  c.SaveAs(outPNG);
}

// ============================== asinh model =================================
/**
 * \brief asinhModel
 * Y(X) = Norm × [ 2·b / sqrt(1 + 4·X² · sinh²(1/(2·b)) ) ]
 */
double asinhModel(double* x, double* par)
{
  double Norm = par[0];
  double b    = par[1];

  if (b < 1e-9) return 1e-12; // prevent division-by-zero anomalies

  double X   = x[0];
  double sh  = TMath::SinH(1.0/(2.0*b));
  double arg = 1.0 + 4.0*X*X*sh*sh;

  return Norm * ( (2.0*b) / TMath::Sqrt(arg) );
}

static BRes FitAsinh1D(TH1D* h, double xmin=-0.5, double xmax=0.5)
{
  BRes best;
  if (!h || h->GetEntries() <= 0) return best;

  TF1 trial("asinh", asinhModel, xmin, xmax, 2);
  trial.SetParNames("Norm", "bVal");
  trial.SetParLimits(0, 1e-9, 1e12);
  trial.SetParLimits(1, 1e-5, 2.0);

  const double peak = std::max(1.0, h->GetMaximum());
  const double bSeeds[] = {0.10, 0.14, 0.17, 0.20};

  double bestChi2 = 1e50;

  for (double b0 : bSeeds)
  {
    const double norm0 = peak / (2.0 * b0);
    trial.SetParameters(norm0, b0);
    TFitResultPtr r = h->Fit(&trial,"RQL0S","", xmin, xmax);
    if (!r.Get() || !r->IsValid()) continue;

    if (r->Chi2() < bestChi2)
    {
      bestChi2  = r->Chi2();
      best.val  = trial.GetParameter(1);   // b
      best.err  = trial.GetParError(1);
      best.norm = trial.GetParameter(0);   // FITTED Norm (use this to draw)
    }
  }
  return best;
}

// ========================= 2D BlockEta vs BlockPhi table ====================

static void Plot2DBlockEtaPhi(TH3F* hUnc3D,
                              TH3F* hCor3D,
                              bool  firstPass,
                              const std::vector<std::pair<double,double>>& eEdges,
                              const char* outDir,
                              const char* etaPretty)
{
  if (!hUnc3D) return;

  gStyle->SetOptStat(0);
  const int N = static_cast<int>(eEdges.size());

  auto tuneAxes = [](TH2* h)
  {
    const double tSize = 0.045, lSize = 0.035;
    h->GetXaxis()->SetTitleSize(tSize);
    h->GetYaxis()->SetTitleSize(tSize);
    h->GetZaxis()->SetTitleSize(tSize);
    h->GetXaxis()->SetLabelSize(lSize);
    h->GetYaxis()->SetLabelSize(lSize);
    h->GetZaxis()->SetLabelSize(lSize);
    h->GetXaxis()->SetTitleOffset(1.10);
    h->GetYaxis()->SetTitleOffset(1.35);
    h->GetZaxis()->SetTitleOffset(1.80);
  };

  // Uncorrected table
  TCanvas c2D_unc("c2D_unc","Uncorrected 2D block coords vs E",2400,1200);
  c2D_unc.Divide(4,2);
    for (int p=1; p<=8; ++p) {
      TPad *pad = (TPad*) c2D_unc.cd(p);
      pad->SetRightMargin(0.18);
      pad->SetLeftMargin (0.12);
      pad->SetBottomMargin(0.12);
      pad->SetTopMargin   (0.16);  // space for header
    }


  for (int i=0;i<N;++i)
  {
    const double eLo = eEdges[i].first;
    const double eHi = eEdges[i].second;
    const TString hdrTxt = Form("UNCORR:  E = %.1f - %.1f  GeV", eLo, eHi);

    int zLo = std::max(1, hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
    int zHi = std::min(hUnc3D->GetNbinsZ(),
                       hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));
    hUnc3D->GetZaxis()->SetRange(zLo,zHi);

    TH2D* h2 = (TH2D*) hUnc3D->Project3D("xy");
    if (!h2) continue;
    h2->SetDirectory(nullptr);
    h2->SetName(Form("h2_unc_xy_Ebin%d", i));
      h2->SetTitle("");  // avoid drawing a histogram title; use TLatex header instead
      h2->GetZaxis()->SetTitle("Entries");
      h2->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
      h2->GetYaxis()->SetTitle("block #varphi_{local, 2#times 2}");
      tuneAxes(h2);

      c2D_unc.cd(i+1);
      h2->Draw("COLZ");

      TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
      // left-top header; sits inside the pad, clear of the color bar
      tl.SetTextSize(0.050); tl.SetTextAlign(13);
      tl.DrawLatex(0.14, 0.96, Form("UNCORR:  E = %.1f - %.1f GeV   %s", eLo, eHi, etaPretty));


  }

  TString out2D_unc = Form("%s/BlockCoord2D_E_unc.png", outDir);
  c2D_unc.SaveAs(out2D_unc);

  // Corrected table only if requested and available
  if (!firstPass && hCor3D)
  {
    TCanvas c2D_cor("c2D_cor","Corrected 2D block coords vs E",2400,1200);
    c2D_cor.Divide(4,2);
      for (int p=1; p<=8; ++p) {
        TPad *pad = (TPad*) c2D_cor.cd(p);
        pad->SetRightMargin(0.18);
        pad->SetLeftMargin (0.12);
        pad->SetBottomMargin(0.12);
        pad->SetTopMargin   (0.16);  // space for header
      }


    for (int i=0;i<N;++i)
    {
      const double eLo = eEdges[i].first;
      const double eHi = eEdges[i].second;
      const TString hdrTxt = Form("CORR:  E = %.1f - %.1f  GeV", eLo, eHi);

      int zLo = std::max(1, hCor3D->GetZaxis()->FindBin(eLo + 1e-9));
      int zHi = std::min(hCor3D->GetNbinsZ(),
                         hCor3D->GetZaxis()->FindBin(eHi - 1e-9));
      hCor3D->GetZaxis()->SetRange(zLo,zHi);

      TH2D* h2 = (TH2D*) hCor3D->Project3D("xy");
      if (!h2) continue;
      h2->SetDirectory(nullptr);
      h2->SetName(Form("h2_cor_xy_Ebin%d", i));
        h2->SetTitle("");  // hide histogram title
        h2->GetZaxis()->SetTitle("Entries");
        h2->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
        h2->GetYaxis()->SetTitle("block #varphi_{local, 2#times 2}");
        tuneAxes(h2);

        c2D_cor.cd(i+1);
        h2->Draw("COLZ");

        TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
        tl.SetTextSize(0.050); tl.SetTextAlign(13);
        tl.DrawLatex(0.14, 0.96, Form("CORR:    E = %.1f - %.1f GeV   %s", eLo, eHi, etaPretty));

    }

    TString out2D_cor = Form("%s/BlockCoord2D_E_cor.png", outDir);
    c2D_cor.SaveAs(out2D_cor);
  }
}

// ====================== Overlay of uncorrected φ & η ========================

static void OverlayUncorrPhiEta(TH3F*  hUnc3D,
                                const std::vector<std::pair<double,double>>& eEdges,
                                const char* outDir,
                                const char* etaPretty)
{
  if (!hUnc3D) return;

  gStyle->SetOptStat(0);
  const int N = static_cast<int>(eEdges.size());

  TCanvas cOv("LocalPhiEtaOverlay_4by2",
              "Uncorr #varphi (red) and #eta (blue) overlays", 1600, 1200);
  cOv.Divide(4,2);

  for (int i=0;i<N;++i)
  {
    const double eLo = eEdges[i].first;
    const double eHi = eEdges[i].second;

    cOv.cd(i+1);

    const int zLo = std::max(1, hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
    const int zHi = std::min(hUnc3D->GetNbinsZ(), hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));

    // φ (Y-projection)
    hUnc3D->GetZaxis()->SetRange(zLo, zHi);
    hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX());

    TH1D* hPhi = static_cast<TH1D*>( hUnc3D->Project3D("y") );
    if (!hPhi) continue;
    hPhi->SetDirectory(nullptr);
    hPhi->SetName(Form("hPhi_E%d", i));
    hPhi->SetMarkerStyle(21);
    hPhi->SetMarkerColor(kRed+1);
    hPhi->SetLineColor  (kRed+1);

    // η (X-projection)
    hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());

    TH1D* hEta = static_cast<TH1D*>( hUnc3D->Project3D("x") );
    if (!hEta) { delete hPhi; continue; }
    hEta->SetDirectory(nullptr);
    hEta->SetName(Form("hEta_E%d", i));
    hEta->SetMarkerStyle(20);
    hEta->SetMarkerColor(kBlue + 1);
    hEta->SetLineColor  (kBlue + 1);

      // Fits (for display) + draw best-fit curves in matching colors
      BRes bPhi = FitAsinh1D(hPhi);
      BRes bEta = FitAsinh1D(hEta);

      const double ymax = std::max(hPhi->GetMaximum(), hEta->GetMaximum());
      hEta->GetYaxis()->SetRangeUser(0.0, 1.30 * std::max(1.0, ymax));

      TString sliceTitle = Form("%.1f < E < %.1f  GeV   (#scale[0.8]{%s})",
                                eLo, eHi, etaPretty);
      hEta->SetTitle(sliceTitle);
      hEta->GetXaxis()->SetTitle("#phi_{CG}/#eta_{CG}");

      // draw data
      hEta->Draw("E");
      hPhi->Draw("E SAME");

      // keep fit functions alive for the whole canvas lifetime
      static std::vector<std::unique_ptr<TF1>> s_keepFits;
      s_keepFits.reserve(2*eEdges.size());  // once per call is fine

      // raw handles for legend entries
      TF1* fEtaPtr = nullptr;
      TF1* fPhiPtr = nullptr;

      // η fit (blue)
      if (std::isfinite(bEta.val)) {
          auto fEta = std::make_unique<TF1>(Form("fEta_E%d", i), asinhModel, -0.5, 0.5, 2);
          if (std::isfinite(bEta.val) && std::isfinite(bEta.norm) && bEta.val > 1e-9) {
              fEta->SetParameters(bEta.norm, bEta.val);  // use FITTED Norm, not a heuristic
          }
        fEta->SetLineColor(kBlue+1);
        fEta->SetLineWidth(2);
        fEta->SetNpx(600);
        hEta->GetListOfFunctions()->Add(fEta.get());        // attach to hist for repaints
        fEta->Draw("L SAME");
        fEtaPtr = fEta.get();                                // keep raw handle for legend
        s_keepFits.emplace_back(std::move(fEta));            // own lifetime
      }

      // φ fit (red)
      if (std::isfinite(bPhi.val)) {
          auto fPhi = std::make_unique<TF1>(Form("fPhi_E%d", i), asinhModel, -0.5, 0.5, 2);
          if (std::isfinite(bPhi.val) && std::isfinite(bPhi.norm) && bPhi.val > 1e-9) {
              fPhi->SetParameters(bPhi.norm, bPhi.val);  // use FITTED Norm
          }
        fPhi->SetLineColor(kRed+1);
        fPhi->SetLineWidth(2);
        fPhi->SetNpx(600);
        hPhi->GetListOfFunctions()->Add(fPhi.get());
        fPhi->Draw("L SAME");
        fPhiPtr = fPhi.get();                                // raw handle for legend
        s_keepFits.emplace_back(std::move(fPhi));
      }

      gPad->Modified();
      gPad->Update();

      // legend: show only data (no fit curves)
      TLegend* leg = new TLegend(0.18,0.82,0.55,0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.038);
      leg->AddEntry(hPhi, "local #varphi", "lp");
      leg->AddEntry(hEta, "local #eta",    "lp");
      // no entries for fit lines
      leg->Draw();

      // print b ± σ_b
      TLatex tl; tl.SetNDC(); tl.SetTextSize(0.040); tl.SetTextAlign(32);
      tl.DrawLatex(0.88, 0.88, Form("b_{#varphi}=%.3g#pm%.2g", bPhi.val, bPhi.err));
      tl.DrawLatex(0.88, 0.84, Form("b_{#eta}=%.3g#pm%.2g",     bEta.val, bEta.err));
  }

  TString outName = Form("%s/LocalPhiEtaOverlay_4by2.png", outDir);
  cOv.SaveAs(outName);
}

// ============== Compute bφ(E) & bη(E) and save a combined plot ==============

struct BVecs { std::vector<double> ecenters, bphi, beta, bphiErr, betaErr; };

static BVecs MakeBvaluesVsEnergyPlot(TH3F* hUnc3D,
                                     const std::vector<std::pair<double,double>>& eEdges,
                                     const char* outDir,
                                     const char* etaPretty)
{
  BVecs R;
  if (!hUnc3D) return R;

  gStyle->SetOptStat(0);
  const int N = static_cast<int>(eEdges.size());

  R.ecenters.reserve(N);
  R.bphi    .reserve(N);
  R.beta    .reserve(N);
  R.bphiErr .reserve(N);
  R.betaErr .reserve(N);

  for (int i=0;i<N;++i)
  {
    const double eLo = eEdges[i].first;
    const double eHi = eEdges[i].second;

    const int zLo = std::max(1, hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
    const int zHi = std::min(hUnc3D->GetNbinsZ(), hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));

    // φ (Y-projection)
    hUnc3D->GetZaxis()->SetRange(zLo, zHi);
    hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX());
    TH1D* hPhi = static_cast<TH1D*>( hUnc3D->Project3D("y") );
    if (!hPhi) continue;
    hPhi->SetDirectory(nullptr);

    // η (X-projection)
    hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());
    TH1D* hEta = static_cast<TH1D*>( hUnc3D->Project3D("x") );
    if (!hEta) { delete hPhi; continue; }
    hEta->SetDirectory(nullptr);

    const BRes bPhi = FitAsinh1D(hPhi);
    const BRes bEta = FitAsinh1D(hEta);

    R.ecenters.push_back(0.5*(eLo+eHi));
    R.bphi.push_back(bPhi.val);
    R.beta.push_back(bEta.val);
    R.bphiErr.push_back(bPhi.err);
    R.betaErr.push_back(bEta.err);

    delete hPhi;
    delete hEta;
  }

  // Build scatter plot
  TGraphErrors* gPhi = new TGraphErrors(R.ecenters.size(), &R.ecenters[0], &R.bphi[0],
                                        nullptr, &R.bphiErr[0]);
  TGraphErrors* gEta = new TGraphErrors(R.ecenters.size(), &R.ecenters[0], &R.beta[0],
                                        nullptr, &R.betaErr[0]);

  gPhi->SetMarkerStyle(21);  gPhi->SetMarkerColor(kRed+1);
  gPhi->SetLineColor  (kRed+1);
  gEta->SetMarkerStyle(20);  gEta->SetMarkerColor(kBlue+1);
  gEta->SetLineColor  (kBlue+1);

  double ymin = +1e9, ymax = -1e9;
  for (size_t i=0;i<R.bphi.size();++i) { ymin = std::min(ymin, R.bphi[i]); ymax = std::max(ymax, R.bphi[i]); }
  for (size_t i=0;i<R.beta.size(); ++i) { ymin = std::min(ymin, R.beta[i]); ymax = std::max(ymax, R.beta[i]); }
  if (!std::isfinite(ymin) || !std::isfinite(ymax) || R.bphi.empty())
  { ymin = 0.0; ymax = 1.0; }
  const double yLo = ymin - 0.15*std::fabs(ymin);
  const double yHi = ymax + 0.25*std::fabs(ymax);

  const double xLo = E_edges[0]-0.5, xHi = E_edges[N]-0.5;

  TH2F* frame = new TH2F("bFrame", Form("best-fit  b  vs  E  (#scale[0.8]{%s});E  [GeV];b", etaPretty),
                         100, xLo, xHi, 100, yLo, yHi);
  frame->SetStats(0);

  TCanvas cB("bValues_vs_E","b values vs energy", 900, 700);
  frame->Draw();
  gPhi->Draw("P SAME");
  gEta->Draw("P SAME");

  TLegend leg(0.20,0.78,0.52,0.90);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
  leg.AddEntry(gPhi,"b_{#varphi}","lp");
  leg.AddEntry(gEta,"b_{#eta}","lp");
  leg.Draw();

  TString outName = Form("%s/bValues_vs_E.png", outDir);
  cB.SaveAs(outName);

  delete frame;
  delete gPhi;
  delete gEta;

  return R;
}

// ==================== (second-pass only) full fit canvases ==================

static void FitLocalPhiEta(TH3F*  hUnc3D,
                           TH3F*  hCor3D,
                           bool   isFirstPassFlag,
                           const std::vector<std::pair<double,double>>& eEdges,
                           const char* outDir)
{
  // This routine exists ONLY for the second pass (it generates
  // LocalPhiFits_4by2.png and LocalEtaFits_4by2.png with corrected overlays).
  if (isFirstPassFlag || !hUnc3D) return;

  // ----------( B )  Local-φ ------------------------------------------------
  TCanvas cFitsPhi("LocalPhiFits_4by2","Local #phi fits",1600,1200);
  cFitsPhi.Divide(4,2);

  const int N_E = static_cast<int>(eEdges.size());

  for (int i = 0; i < N_E; ++i)
  {
    double eLo = eEdges[i].first , eHi = eEdges[i].second;
    cFitsPhi.cd(i+1);

    int zLo = std::max(1 , hUnc3D->GetZaxis()->FindBin(eLo+1e-9));
    int zHi = std::min(hUnc3D->GetNbinsZ(),
                       hUnc3D->GetZaxis()->FindBin(eHi-1e-9));

    hUnc3D->GetZaxis()->SetRange(zLo,zHi);
    hUnc3D->GetXaxis()->SetRange(1,hUnc3D->GetNbinsX());

    std::unique_ptr<TH1D> hUnc( static_cast<TH1D*>(hUnc3D->Project3D("y")) );
    if(!hUnc){ continue; }
    hUnc->SetDirectory(nullptr);
    hUnc->SetName (Form("hUncPhi_E%d",i));
    hUnc->SetTitle(Form("Local #phi : E=[%.1f,%.1f)",eLo,eHi));
    hUnc->GetXaxis()->SetTitle("block #phi_{local, 2#times 2}");
    hUnc->GetYaxis()->SetTitle("counts");
    hUnc->SetMarkerStyle(20);
    hUnc->Draw("E");

    // Fit
    BRes best = FitAsinh1D(hUnc.get());
    TF1* bestFit = nullptr;
    if (std::isfinite(best.val))
    {
      bestFit = new TF1("bestF_phi", asinhModel, -0.5, 0.5, 2);
      bestFit->SetParameters(hUnc->GetMaximum(), best.val);
      bestFit->SetLineColor(kBlue+1);
      bestFit->SetLineWidth(2);
      bestFit->Draw("SAME");
    }

    // corrected overlay if provided
    TH1D* hCorRaw = nullptr;
    if (hCor3D){
      hCor3D->GetZaxis()->SetRange(zLo,zHi);
      hCor3D->GetXaxis()->SetRange(1,hCor3D->GetNbinsX());
      std::unique_ptr<TH1D> hCor( static_cast<TH1D*>(hCor3D->Project3D("y")) );
      if(hCor){
        hCor->SetDirectory(nullptr);
        hCor->SetName(Form("hCorPhi_E%d",i));
        hCor->SetMarkerStyle(21); hCor->SetMarkerColor(kRed);
        hCor->SetLineColor(kRed); hCor->Draw("SAME E");
        hCorRaw=hCor.release();
      }
    }

    TLegend lg(0.54,0.75,0.88,0.84); lg.SetBorderSize(0);
    lg.SetFillStyle(0); lg.SetTextSize(0.035);
    lg.AddEntry(hUnc.get(), "Uncorrected", "lp");
    if(hCorRaw) lg.AddEntry(hCorRaw, "Corrected", "lp");
    if(bestFit) lg.AddEntry(bestFit, Form("asinh fit (b=%.3g)", best.val), "l");
    lg.Draw();
  }

  cFitsPhi.SaveAs(Form("%s/LocalPhiFits_4by2.png",outDir));

  // ----------( C )  Local-η  --------------------------------------------
  TCanvas cFitsEta("LocalEtaFits_4by2","Local #eta fits",1600,1200);
  cFitsEta.Divide(4,2);

  for (int i = 0; i < N_E; ++i)
  {
    double eLo = eEdges[i].first , eHi = eEdges[i].second;
    cFitsEta.cd(i+1);

    int zLo = std::max(1 , hUnc3D->GetZaxis()->FindBin(eLo+1e-9));
    int zHi = std::min(hUnc3D->GetNbinsZ(),
                       hUnc3D->GetZaxis()->FindBin(eHi-1e-9));

    hUnc3D->GetZaxis()->SetRange(zLo,zHi);
    hUnc3D->GetYaxis()->SetRange(1,hUnc3D->GetNbinsY());

    std::unique_ptr<TH1D> hUnc( static_cast<TH1D*>(hUnc3D->Project3D("x")) );
    if(!hUnc){ continue; }
    hUnc->SetDirectory(nullptr);
    hUnc->SetName (Form("hUncEta_E%d",i));
    hUnc->SetTitle(Form("Local #eta : E=[%.1f,%.1f)",eLo,eHi));
    hUnc->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    hUnc->GetYaxis()->SetTitle("counts");
    hUnc->SetMarkerStyle(20);
    hUnc->Draw("E");

    // Fit
    BRes best = FitAsinh1D(hUnc.get());
    TF1* bestFit = nullptr;
    if (std::isfinite(best.val))
    {
      bestFit = new TF1("bestF_eta", asinhModel, -0.5, 0.5, 2);
      bestFit->SetParameters(hUnc->GetMaximum(), best.val);
      bestFit->SetLineColor(kBlue+1);
      bestFit->SetLineWidth(2);
      bestFit->Draw("SAME");
    }

    TH1D* hCorRaw = nullptr;
    if (hCor3D){
      hCor3D->GetZaxis()->SetRange(zLo,zHi);
      hCor3D->GetYaxis()->SetRange(1,hCor3D->GetNbinsY());
      std::unique_ptr<TH1D> hCor( static_cast<TH1D*>(hCor3D->Project3D("x")) );
      if(hCor){
        hCor->SetDirectory(nullptr);
        hCor->SetName(Form("hCorEta_E%d",i));
        hCor->SetMarkerStyle(21); hCor->SetMarkerColor(kRed);
        hCor->SetLineColor(kRed); hCor->Draw("SAME E");
        hCorRaw=hCor.release();
      }
    }

    TLegend lg(0.54,0.75,0.88,0.84); lg.SetBorderSize(0);
    lg.SetFillStyle(0); lg.SetTextSize(0.035);
    lg.AddEntry(hUnc.get(), "Uncorrected", "lp");
    if(hCorRaw) lg.AddEntry(hCorRaw, "Corrected", "lp");
    if(bestFit) lg.AddEntry(bestFit, Form("asinh fit (b=%.3g)", best.val), "l");
    lg.Draw();
  }

  cFitsEta.SaveAs(Form("%s/LocalEtaFits_4by2.png",outDir));
}

// ======================= Write per-slice b-values to file ===================

static void WriteBValuesTxt(TH3F* hUnc3D,
                            const std::vector<std::pair<double,double>>& eEdges,
                            std::ofstream& bOut,
                            const char* etaVariantKey,
                            const char* zVariantKey)
{
  if (!hUnc3D) return;
  if (!bOut.is_open()) return;

  const int N = static_cast<int>(eEdges.size());
  for (int i=0;i<N;++i)
  {
    const double eLo = eEdges[i].first;
    const double eHi = eEdges[i].second;

    const int zLo = std::max(1, hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
    const int zHi = std::min(hUnc3D->GetNbinsZ(), hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));

    // φ
    hUnc3D->GetZaxis()->SetRange(zLo, zHi);
    hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX());
    std::unique_ptr<TH1D> hPhi( static_cast<TH1D*>(hUnc3D->Project3D("y")) );
    if (hPhi) {
      BRes bPhi = FitAsinh1D(hPhi.get());
      if (std::isfinite(bPhi.val))
        bOut << Form("PHI [%.1f,%.1f)  %.6f  %s  %s\n", eLo, eHi, bPhi.val, etaVariantKey, zVariantKey);
    }

    // η
    hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());
    std::unique_ptr<TH1D> hEta( static_cast<TH1D*>(hUnc3D->Project3D("x")) );
    if (hEta) {
      BRes bEta = FitAsinh1D(hEta.get());
      if (std::isfinite(bEta.val))
        bOut << Form("ETA [%.1f,%.1f)  %.6f  %s  %s\n", eLo, eHi, bEta.val, etaVariantKey, zVariantKey);
    }
  }
}


static void SaveOverlayAcrossVariants(
    const std::string& outBaseDir,
    const std::vector<double>& ecenters,
    const std::map<std::string, std::vector<double>>& byVar,
    const std::map<std::string, std::vector<double>>& byVarErr,
    bool isPhi,
    const char* outSuffix)
{
  if (byVar.empty()) return;

  // ----------------------------- configuration -----------------------------
  const double E0   = 3.0;                       // reference energy (GeV)
  const double xLo  = E_edges[0] - 0.5;          // assume these globals exist
  const double xHi  = E_edges[N_E] - 0.5;

  // pretty label: try η keys, then fine-z, then coarse-z; else echo key
  auto prettyName = [&](const std::string& key)->const char*
  {
    auto itEta = kEtaPretty.find(key);
    if (itEta != kEtaPretty.end()) return itEta->second.c_str();
    auto itZf  = kZFinePretty.find(key);
    if (itZf  != kZFinePretty.end()) return itZf->second.c_str();
    auto itZc  = kZPretty.find(key);
    return (itZc == kZPretty.end()) ? key.c_str() : itZc->second.c_str();
  };

    // ----------------------------- y-range (±σ if available) ------------------
    double ymin = +1e9, ymax = -1e9;
    for (const auto& kv : byVar)
    {
      const auto& key  = kv.first;
      const auto& vals = kv.second;
      const auto itErr = byVarErr.find(key);
      const std::vector<double>* errs = (itErr == byVarErr.end()) ? nullptr : &itErr->second;

      for (size_t i = 0; i < vals.size(); ++i)
      {
        if (!std::isfinite(vals[i])) continue;
        const double e = (errs && i < errs->size() && std::isfinite((*errs)[i])) ? (*errs)[i] : 0.0;
        ymin = std::min(ymin, vals[i] - e);
        ymax = std::max(ymax, vals[i] + e);
      }
    }
    if (!std::isfinite(ymin) || !std::isfinite(ymax)) { ymin = 0.0; ymax = 1.0; }
    double yLo = ymin - 0.05 * std::fabs(ymin);
    double yHi = ymax + 0.10 * std::fabs(ymax);

    // Force axis range for η overlays with "_zOnly" suffix
    // Top fixed at 0.42; add extra bottom buffer down to 0.14
    if (!isPhi && outSuffix && std::string(outSuffix).find("_zOnly") != std::string::npos) {
      yHi = 0.42;
      yLo = 0.14;
      // safety: if yLo accidentally crosses yHi, keep at least a small span
      if (yLo >= yHi - 1e-6) yLo = std::max(0.0, yHi - 0.10);
    }

    // ----------------------------- canvas & frame -----------------------------
    TCanvas c("cbOverlay",
              isPhi ? "b_{#varphi}(E) overlay" : "b_{#eta}(E) overlay",
              900, 700);
    TH2F frame("frame", "", 100, xLo, xHi, 100, yLo, yHi);
  frame.SetTitle(isPhi ? "best-fit  b_{#varphi}  vs  E;E  [GeV];b_{#varphi}"
                       : "best-fit  b_{#eta}     vs  E;E  [GeV];b_{#eta}");
  frame.SetStats(0);
  frame.Draw();

  const int colors[8] = { kBlack, kBlue+1, kRed+1, kGreen+2,
                          kMagenta+1, kOrange+7, kCyan+2, kGray+2 };
  const int    kMarkerStyle   = 20;
  const double kMarkerSize    = 1.2;

  TLegend leg(0.16, 0.65, 0.42, 0.89);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.040);

  // TLatex overlay disabled: no on-canvas fit text
  std::vector<std::unique_ptr<TGraphErrors>> keepPoints;
  std::vector<std::unique_ptr<TF1>>          keepFits;

    // ----------------------------- terminal banner ----------------------------
    {
      const char* B = "\033[1m";   // ANSI bold
      const char* R = "\033[0m";   // ANSI reset

      std::ostringstream t;
      t << (isPhi ? "[PHI] " : "[ETA] ")
        << "b(E) = b0 + m * ln(E/E0),  E0 = " << E0 << " GeV";
      std::cout << "\n" << B << t.str() << R << "\n";
      std::cout << "Weighted least squares per variant: w_i = 1/sigma_i^2 (if provided) or 1 otherwise.\n"
                << "Parameter errors: unscaled if errors provided; scaled by sqrt(χ²/ndf) when unit weights are used.\n\n";

      // Bold, clearly separated header
      std::cout << B
                << std::left  << std::setw(18) << "variant" << R
                << " | " << std::right << std::setw(6)  << "N"
                << " | " << std::setw(9)  << "E[min]"
                << " | " << std::setw(9)  << "E[max]"
                << " | " << std::setw(12) << "b0@E0"
                << " | " << std::setw(10) << "±b0"
                << " | " << std::setw(12) << "m (ln E)"
                << " | " << std::setw(10) << "±m"
                << " | " << std::setw(9)  << "z(m)"
                << " | " << std::setw(12) << "m/decade"
                << " | " << std::setw(12) << "Δb_fit"
                << " | " << std::setw(12) << "Δb_data"
                << " | " << std::setw(8)  << "R^2"
                << " | " << std::setw(10) << "RMSE_w"
                << " | " << std::setw(10) << "χ²/ndf"
                << "\n";

      // Horizontal rule
      std::cout << std::string(
                     18+2 + 6+3 + 9+3 + 9+3 + 12+3 + 10+3 + 12+3 + 10+3 +
                     9+3 + 12+3 + 12+3 + 8+3 + 10+3 + 10, '-') << "\n";

      // Numeric formatting for table rows
      std::cout.setf(std::ios::fixed);
      std::cout << std::setprecision(6);
    }
  // Desired η-order for legends (if η-keys are used)
  auto containsKey = [](const std::vector<std::string>& v, const std::string& s)->bool {
    for (const auto& x : v) if (x == s) return true; return false;
  };
  std::vector<std::string> desiredEta = {"etaCore","etaMid","etaEdge","fullEta","originalEta"};

    // Build ordered key list for PLOTTING, then extend a second list for TABULATION.
    //  - If all fine-z tags → use kZFineOrder
    //  - else if all coarse-z tags → use kZOrder
    //  - else (η overlays) → EXCLUDE "originalEta" from plotting order,
    //                        but APPEND it to the table order.
    std::vector<std::string> keysOrdered;     // for plotting
    std::vector<std::string> keysForTable;    // for printing/CSV (may include originalEta)
    bool allCoarseZ = true, allFineZ = true;
    for (const auto& kv : byVar) {
      if (kZPretty.find(kv.first)     == kZPretty.end())     allCoarseZ = false;
      if (kZFinePretty.find(kv.first) == kZFinePretty.end()) allFineZ   = false;
    }

    if (allFineZ) {
      for (const auto& z : kZFineOrder)
        if (byVar.find(z) != byVar.end()) keysOrdered.push_back(z);
      keysForTable = keysOrdered;  // z-only: table == plot
    } else if (allCoarseZ) {
      for (const auto& z : kZOrder)
        if (byVar.find(z) != byVar.end()) keysOrdered.push_back(z);
      keysForTable = keysOrdered;  // z-only: table == plot
    } else {
      // η overlays: plot excludes "originalEta"
      for (const auto& k : desiredEta)
        if (k != "originalEta" && byVar.find(k) != byVar.end())
          keysOrdered.push_back(k);
      for (const auto& kv : byVar)
        if (kv.first != "originalEta" && !containsKey(keysOrdered, kv.first))
          keysOrdered.push_back(kv.first);

      // Table should ALSO include "originalEta" if present
      keysForTable = keysOrdered;
      if (byVar.find("originalEta") != byVar.end() &&
          !containsKey(keysForTable, "originalEta")) {
        keysForTable.push_back("originalEta");
      }
    }

    // Collect CSV rows (one line per variant)
    std::vector<std::string> csvRows;
    {
      std::ostringstream hdr;
      hdr << "variant,N,E_min,E_max,b0@E0,db0,m_lnE,dm_lnE,z_m,m_per_decade,Delta_b_fit,Delta_b_data,R2w,RMSEw,chi2_ndf";
      csvRows.push_back(hdr.str());
    }

    int idx = 0;
    for (const auto& key : keysForTable)
    {
      const auto& vals = byVar.at(key);
      if (vals.size() != ecenters.size()) { ++idx; continue; }

      auto itErr = byVarErr.find(key);
      std::vector<double> errs = (itErr != byVarErr.end()) ? itErr->second
                                                           : std::vector<double>(vals.size(), 0.0);
      if (errs.size() != vals.size()) errs.assign(vals.size(), 0.0);

      std::unique_ptr<TGraphErrors> g(new TGraphErrors(static_cast<int>(vals.size())));
      double Emin = +1e9, Emax = -1e9;
      double yMin = +1e9, yMax = -1e9;

      // Load points to graph and track data range
      for (size_t i = 0; i < vals.size(); ++i)
      {
        const double Ei = ecenters[i];
        const double yi = vals[i];
        const double si = (std::isfinite(errs[i]) ? errs[i] : 0.0);
        if (!std::isfinite(Ei) || !std::isfinite(yi)) continue;

        g->SetPoint(static_cast<int>(i), Ei, yi);
        g->SetPointError(static_cast<int>(i), 0.0, si);
        Emin = std::min(Emin, Ei);
        Emax = std::max(Emax, Ei);
        yMin = std::min(yMin, yi);
        yMax = std::max(yMax, yi);
      }

        const int col = colors[idx % 8];
        bool drawThis = (allCoarseZ || allFineZ) || (key != "originalEta"); // do not draw originalEta on η overlays

        if (drawThis) {
          g->SetMarkerStyle(kMarkerStyle);
          g->SetMarkerSize(kMarkerSize);
          g->SetMarkerColor(col);
          g->SetLineColor(col);
          g->Draw("P SAME");
          leg.AddEntry(g.get(), prettyName(key), "p");
          keepPoints.emplace_back(std::move(g));
        } else {
          // keep the object alive for lifetime safety even if not drawn
          keepPoints.emplace_back(std::move(g));
        }

      // closed-form weighted least squares for line in log(E/E0)
      double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
      int     Nuse = 0;

      std::vector<double> xi; xi.reserve(vals.size());
      std::vector<double> yi; yi.reserve(vals.size());
      std::vector<double> wi; wi.reserve(vals.size());

      bool haveYerrsAll = true; // if any sigma == 0, we'll scale errors by chi2/ndf
      for (size_t i = 0; i < vals.size(); ++i)
      {
        const double E  = ecenters[i];
        const double y  = vals[i];
        const double sy = errs[i];

        if (!(std::isfinite(E) && E > 0.0 && std::isfinite(y))) continue;

        const double x = std::log(E / E0);
        double w = 1.0;
        if (std::isfinite(sy) && sy > 0.0) w = 1.0/(sy*sy); else haveYerrsAll = false;

        S   += w;
        Sx  += w * x;
        Sy  += w * y;
        Sxx += w * x * x;
        Sxy += w * x * y;

        xi.push_back(x);
        yi.push_back(y);
        wi.push_back(w);
        ++Nuse;
      }

      double b0 = std::numeric_limits<double>::quiet_NaN();
      double m  = std::numeric_limits<double>::quiet_NaN();
      double db0= std::numeric_limits<double>::quiet_NaN();
      double dm = std::numeric_limits<double>::quiet_NaN();
      double chi2=0.0, R2w=0.0, RMSEw=0.0, redChi2=0.0;
      double dBfit=0.0, dBdata=0.0, z_m=0.0, m_per_dec=0.0;

      const double Delta = S*Sxx - Sx*Sx;
      if (Nuse >= 2 && (Delta > 0.0))
      {
        // Fit parameters
        m  = (S * Sxy - Sx * Sy) / Delta;
        b0 = (Sxx * Sy - Sx * Sxy) / Delta;

        // Parameter uncertainties from WLS normal equations
        const double var_m  =  S / Delta;
        const double var_b0 = Sxx / Delta;

        // Goodness metrics
        double ybar_w = Sy / S;
        double TSSw = 0.0; // weighted total sum of squares
        for (int i = 0; i < Nuse; ++i)
        {
          const double yhat = b0 + m * xi[i];
          const double r    = yi[i] - yhat;
          chi2  += wi[i] * r * r;
          TSSw  += wi[i] * (yi[i] - ybar_w) * (yi[i] - ybar_w);
          if (i == 0) { dBfit = yhat; dBdata = yi[i]; }
          else {
            dBfit = std::max(dBfit, yhat);
          }
        }
        // Δb_data from raw values
        dBdata = yMax - yMin;

        const int ndf = std::max(0, Nuse - 2);
        redChi2 = (ndf > 0) ? (chi2 / ndf) : 0.0;
        R2w     = (TSSw > 0) ? (1.0 - chi2 / TSSw) : 0.0;
        RMSEw   = (S > 0.0) ? std::sqrt(chi2 / S) : 0.0;

        // If we did not have per-point uncertainties, scale parameter errors by sqrt(chi2/ndf)
        const double scale = (!haveYerrsAll && ndf > 0) ? std::sqrt(redChi2) : 1.0;
        dm  = std::sqrt(std::max(0.0, var_m )) * scale;
        db0 = std::sqrt(std::max(0.0, var_b0)) * scale;

        // Report extras
        z_m        = (dm > 0.0 ? m / dm : 0.0);
        m_per_dec  = m * std::log(10.0);

        // Recompute Δb_fit using predicted range over actually used energy range (min/max xi)
        if (!xi.empty())
        {
          const double x_min = *std::min_element(xi.begin(), xi.end());
          const double x_max = *std::max_element(xi.begin(), xi.end());
          const double yhat_min = b0 + m * x_min;
          const double yhat_max = b0 + m * x_max;
          dBfit = yhat_max - yhat_min;
        }

        // Draw the fitted line
        std::unique_ptr<TF1> f(new TF1(Form("f_%s_%d", key.c_str(), idx),
                                       "[0] + [1]*log(x/[2])", xLo, xHi));
        f->SetParameters(b0, m, E0);
        f->FixParameter(2, E0);
          f->SetLineColor(colors[idx % 8]);
          f->SetLineWidth(2);
          if (drawThis) f->Draw("SAME");
          keepFits.emplace_back(std::move(f));
      }

        // Print one clean, aligned row per variant
        {
          std::ostringstream ln;
          ln.setf(std::ios::fixed);
          ln << std::setprecision(6);

          ln << std::left << std::setw(18) << key
             << " | " << std::right << std::setw(6)  << Nuse
             << " | " << std::setw(9)  << (std::isfinite(Emin)?Emin:0.0)
             << " | " << std::setw(9)  << (std::isfinite(Emax)?Emax:0.0)
             << " | " << std::setw(12) << (std::isfinite(b0)?b0:0.0)
             << " | " << std::setw(10) << (std::isfinite(db0)?db0:0.0)
             << " | " << std::setw(12) << (std::isfinite(m)?m:0.0)
             << " | " << std::setw(10) << (std::isfinite(dm)?dm:0.0)
             << " | " << std::setw(9)  << (std::isfinite(z_m)?z_m:0.0)
             << " | " << std::setw(12) << (std::isfinite(m_per_dec)?m_per_dec:0.0)
             << " | " << std::setw(12) << (std::isfinite(dBfit)?dBfit:0.0)
             << " | " << std::setw(12) << (std::isfinite(dBdata)?dBdata:0.0)
             << " | " << std::setw(8)  << (std::isfinite(R2w)?R2w:0.0)
             << " | " << std::setw(10) << (std::isfinite(RMSEw)?RMSEw:0.0)
             << " | " << std::setw(10) << (std::isfinite(redChi2)?redChi2:0.0);

          std::cout << ln.str() << "\n";
        }
      // CSV line (same fields)
      {
        std::ostringstream ln;
        ln << key << ","
           << Nuse << ","
           << (std::isfinite(Emin)?Emin:0.0) << ","
           << (std::isfinite(Emax)?Emax:0.0) << ","
           << (std::isfinite(b0)?b0:0.0) << ","
           << (std::isfinite(db0)?db0:0.0) << ","
           << (std::isfinite(m)?m:0.0) << ","
           << (std::isfinite(dm)?dm:0.0) << ","
           << (std::isfinite(z_m)?z_m:0.0) << ","
           << (std::isfinite(m_per_dec)?m_per_dec:0.0) << ","
           << (std::isfinite(dBfit)?dBfit:0.0) << ","
           << (std::isfinite(dBdata)?dBdata:0.0) << ","
           << (std::isfinite(R2w)?R2w:0.0) << ","
           << (std::isfinite(RMSEw)?RMSEw:0.0) << ","
           << (std::isfinite(redChi2)?redChi2:0.0);
        csvRows.push_back(ln.str());
      }

      ++idx;
    }

    // Write CSV once per overlay call
    {
      const char* which = isPhi ? "phi" : "eta";
      std::string suff  = (outSuffix ? std::string(outSuffix) : std::string());
      std::string csvPath = outBaseDir + "/bFitParams_" + which + (suff.empty()? "" : suff) + ".csv";
      std::ofstream fout(csvPath);
      if (fout)
      {
        for (const auto& s : csvRows) fout << s << "\n";
        fout.close();
        std::cout << "[b-fit] Wrote CSV: " << csvPath << "\n";
      }
      else
      {
        std::cerr << "[b-fit] WARNING: Could not write CSV to " << csvPath << "\n";
      }
    }

    leg.Draw();

  // Single overlay image (NO redundant second save)
  TString outMain = Form("%s/%s%s.png", outBaseDir.c_str(),
                         isPhi ? "bValuesPhiOverlay" : "bValuesEtaOverlay",
                         outSuffix ? outSuffix : "");
  c.SaveAs(outMain);
}


/*
 CODE FOR RESIDUALS
 */


// Make this visible to any function that draws residual panels.
struct ResidualSpec {
    TH1F*       h;         // histogram to draw (original; caller controls lifetime)
    const char* label;     // legend label
    Color_t     col;       // marker/line color
    Style_t     mk;        // marker style
    int         ls;        // (unused in current panel, kept for completeness)
};

// Lowercase utility
static inline std::string toLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

// Return the name of the first non-null histogram in a vector
static inline std::string firstHistName(const std::vector<TH1F*>& V) {
  for (TH1F* h : V) if (h) return std::string(h->GetName());
  return std::string();
}

// Decide whether a set of hists is φ or η from names/titles
static inline std::string detectResidual(const std::vector<TH1F*>& V) {
  const std::string n = toLower(firstHistName(V));
  if (n.find("phi") != std::string::npos) return "phi";
  if (n.find("eta") != std::string::npos) return "eta";
  for (TH1F* h : V) if (h) {
    std::string t = toLower(h->GetTitle());
    if (t.find("phi") != std::string::npos) return "phi";
    if (t.find("eta") != std::string::npos) return "eta";
  }
  return "resid";
}

// Map a histogram name to the variant label used elsewhere in the file
static inline std::string variantFromHistName(const std::string& hn) {
  const std::string n = toLower(hn);
  // Clusterizer family
  if (n.find("cpcorrea") != std::string::npos) return "CorrectPosition(EA), cluster";
  if (n.find("cpbcorr")  != std::string::npos) return "b(E) corr, cluster";
  if (n.find("cpcorr")   != std::string::npos) return "CorrectPosition, cluster";
  if (n.find("cpraw")    != std::string::npos) return "no corr, cluster";
  // From-scratch (PDC) family
  if (n.find("_corr_")   != std::string::npos) return "b(E) corr, scratch";
  if (n.find("_raw_")    != std::string::npos) return "no corr, scratch";
  return "unknown";
}

// Save multiple μ/σ(E) series to a CSV next to the plot file.
// CSV columns: series,variant,residual,E,mu,dmu,sigma,dsigma
static inline void saveMuSigmaCSV(
  const std::string& outArtifact,
  const std::vector<std::string>& variants,
  const std::vector<std::string>& residuals,
  const std::vector<std::vector<double>>& eCtrs,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& dmus,
  const std::vector<std::vector<double>>& sigmas,
  const std::vector<std::vector<double>>& dsigmas)
{
  // Derive foo.csv from foo.png/foo.pdf/etc.
  std::string csv = outArtifact;
  const auto dot = csv.find_last_of('.');
  if (dot != std::string::npos) csv.replace(dot, std::string::npos, ".csv");
  else                          csv.append(".csv");

  std::ofstream f(csv);
  if (!f.is_open()) {
    std::cerr << "[Playground][ERROR] Cannot write CSV: " << csv << "\n";
    return;
  }

  f << "series,variant,residual,E,mu,dmu,sigma,dsigma\n";

  const size_t S = std::max({
    variants.size(), residuals.size(),
    eCtrs.size(), mus.size(), dmus.size(), sigmas.size(), dsigmas.size()
  });

  for (size_t s = 0; s < S; ++s) {
    const auto &E   = (s < eCtrs.size() ? eCtrs[s] : std::vector<double>{});
    const auto &Mu  = (s < mus.size()   ? mus[s]   : std::vector<double>{});
    const auto &dMu = (s < dmus.size()  ? dmus[s]  : std::vector<double>{});
    const auto &Sg  = (s < sigmas.size()? sigmas[s]: std::vector<double>{});
    const auto &dSg = (s < dsigmas.size()? dsigmas[s]: std::vector<double>{});

    const size_t N = std::min({E.size(), Mu.size(), dMu.size(), Sg.size(), dSg.size()});
    for (size_t i = 0; i < N; ++i) {
      f << s << ','
        << (s < variants.size()  ? variants[s]  : "") << ','
        << (s < residuals.size() ? residuals[s] : "") << ','
        << E[i]   << ','
        << Mu[i]  << ','
        << dMu[i] << ','
        << Sg[i]  << ','
        << dSg[i] << '\n';
    }
  }

  f.close();
  std::cout << "[Playground] Wrote " << csv << "\n";
}

// Template parameters let us accept your existing lambdas by perfect forwarding.
// No std::function, no indirection, no behavior change.
template <class DrawResidualPanelT,
          class MakeEqualPadsT,
          class StatsFromHistRangeT>
static void MakeThreeWaySummaries(
    const char* outDir,
    const std::vector<std::pair<double,double>>& eEdges,
    double xMin, double xMax,

    // CLUS (Constant-b, "CP family") histogram sets
    //   φ : h_phi_diff_cpRaw_%s, h_phi_diff_cpCorr_%s, h_phi_diff_cpCorrEA_%s
    //   η : h_eta_diff_cpRaw_%s, h_eta_diff_cpCorr_%s, h_eta_diff_cpCorrEA_%s
    const std::vector<TH1F*>& Hphi_cpRaw,
    const std::vector<TH1F*>& Hphi_cpCorr,
    const std::vector<TH1F*>& Hphi_cpCorrEA,
    const std::vector<TH1F*>& Heta_cpRaw,
    const std::vector<TH1F*>& Heta_cpCorr,
    const std::vector<TH1F*>& Heta_cpCorrEA,

    // helpers you already defined above inside your playground
    DrawResidualPanelT&&      drawResidualPanel,
    MakeEqualPadsT&&          makeEqualHeightTwinPads,
    StatsFromHistRangeT&&     statsFromHistRange)
{
    // =====================================================================
    // NEW: 3WAY summary overlays + μ/σ vs E (RAW‑CP, CP‑corr, CP‑corr(EA))
    //       Outputs into: <outDir>/3WAYneededSUMMARIES
    // =====================================================================
    const std::string dir3Way = std::string(outDir) + "/3WAYneededSUMMARIES";
    ensureDir(dir3Way.c_str());

    auto nonEmpty = [](TH1F* h){ return h && h->Integral() > 0; };

    // Find the N-th energy slice where BOTH hists in a pair are non-empty
    auto nthNonEmptyPair = [&](const std::vector<TH1F*>& A,
                               const std::vector<TH1F*>& B,
                               int nth)->int
    {
        int count = 0;
        for (std::size_t i=0; i<eEdges.size(); ++i) {
            if (i<A.size() && i<B.size() && nonEmpty(A[i]) && nonEmpty(B[i])) {
                ++count;
                if (count == nth) return static_cast<int>(i);
            }
        }
        return -1;
    };

    // Find the N-th energy slice where ALL THREE are non-empty
    auto nthNonEmptyTriplet = [&](const std::vector<TH1F*>& A,
                                  const std::vector<TH1F*>& B,
                                  const std::vector<TH1F*>& C,
                                  int nth)->int
    {
        int count = 0;
        for (std::size_t i=0; i<eEdges.size(); ++i) {
            if (i<A.size() && i<B.size() && i<C.size() &&
                nonEmpty(A[i]) && nonEmpty(B[i]) && nonEmpty(C[i]))
            {
                ++count;
                if (count == nth) return static_cast<int>(i);
            }
        }
        return -1;
    };

    // 2-series (cpRaw vs cpCorr or cpRaw vs cpCorrEA) for the N-th non-empty energy bin
    auto overlayTwoSeriesNthBin = [&](bool isPhi,
                                      const std::vector<TH1F*>& Vraw,
                                      const std::vector<TH1F*>& Valt,
                                      const char* altShort,          // "CP-corr" or "CP-corr(EA)"
                                      Color_t color,                  // kRed+1 for φ, kBlue+1 for η
                                      const std::string& outPngStemBase,
                                      int nth)
    {
        const int iE = nthNonEmptyPair(Vraw, Valt, nth);
        if (iE < 0) {
            std::cerr << "[Playground][WARN] 3WAY overlay skipped for "
                      << (isPhi? "phi" : "eta")
                      << " - no non-empty bin #" << nth << " for " << altShort << " vs RAW-CP.\n";
            return;
        }

        const char* suffix =
            (nth == 1) ? "_FirstBin" :
            (nth == 2) ? "_SecondBin" : Form("_Bin%d", nth);

        TCanvas c(Form("c3Way_%s_%s_bin%d", isPhi? "phi" : "eta", altShort, nth),
                  "3WAY Nth-bin overlay", 1000, 750);

        const char* delta = isPhi ? "#Delta#phi" : "#Delta#eta";
        std::vector<ResidualSpec> specs = {
            { Vraw[iE], Form("%s RAW-CP",      delta), color, (Style_t)20, 1 },  // closed circle
            { Valt[iE], Form("%s %s",          delta, altShort), color, (Style_t)24, 1 }   // open circle
        };

        std::string header = Form("%s  -  First Bin Overlay: RAW-CP vs %s",
                                  delta, altShort);

        drawResidualPanel(&c,
                          header,
                          eEdges[iE].first, eEdges[iE].second,
                          specs, xMin, xMax, 1.0);

        const std::string outPng = dir3Way + "/" + outPngStemBase + std::string(suffix) + ".png";
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // 3-series (cpRaw, cpCorr, cpCorrEA) for the N-th non-empty energy bin
    auto overlayThreeSeriesNthBin = [&](bool isPhi,
                                        const std::vector<TH1F*>& Vraw,
                                        const std::vector<TH1F*>& Vcor,
                                        const std::vector<TH1F*>& VcorEA,
                                        Color_t color,
                                        const std::string& outPngStemBase,
                                        int nth)
    {
        const int iE = nthNonEmptyTriplet(Vraw, Vcor, VcorEA, nth);
        if (iE < 0) {
            std::cerr << "[Playground][WARN] 3WAY triplet overlay skipped for "
                      << (isPhi? "phi" : "eta")
                      << " - no non-empty bin #" << nth << " with RAW-CP, CP-corr, CP-corr(EA).\n";
            return;
        }

        const char* suffix =
            (nth == 1) ? "_FirstBin" :
            (nth == 2) ? "_SecondBin" : Form("_Bin%d", nth);

        TCanvas c(Form("c3Way_trip_%s_bin%d", isPhi? "phi" : "eta", nth),
                  "3WAY Nth-bin triplet overlay", 1000, 750);

        const char* delta = isPhi ? "#Delta#phi" : "#Delta#eta";
        std::vector<ResidualSpec> specs = {
            { Vraw  [iE], Form("%s RAW-CP",        delta), color, (Style_t)20, 1 }, // closed circle
            { Vcor  [iE], Form("%s CP-corr",       delta), color, (Style_t)24, 1 }, // open circle
            { VcorEA[iE], Form("%s CP-corr(EA)",   delta), color, (Style_t)25, 1 }  // open square
        };

        std::string header = Form("%s  -  First Bin Overlay: RAW-CP / CP-corr / CP-corr(EA)",
                                  delta);

        drawResidualPanel(&c,
                          header,
                          eEdges[iE].first, eEdges[iE].second,
                          specs, xMin, xMax, 1.0);

        const std::string outPng = dir3Way + "/" + outPngStemBase + std::string(suffix) + ".png";
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // ---- Three-variant μ/σ vs E summary (φ and η) to the same folder ----
    struct MuSigSeries { std::vector<double> E, mu, dmu, sg, dsg; };

    auto makeMuSigSeries = [&](const std::vector<TH1F*>& V)->MuSigSeries {
        MuSigSeries S;
        for (std::size_t i=0; i<eEdges.size(); ++i){
            if (i>=V.size() || !nonEmpty(V[i])) continue;
            double m, me, r, re;
            std::tie(m, me, r, re) = statsFromHistRange(V[i], xMin, xMax);
            S.E .push_back(0.5*(eEdges[i].first + eEdges[i].second));
            S.mu.push_back(m);  S.dmu.push_back(me);
            S.sg.push_back(r);  S.dsg.push_back(re);
        }
        return S;
    };

    // Supports forced y-ranges so φ and η plots can share identical axes.
    auto drawMuSigmaThree = [&](bool isPhi,
                                const MuSigSeries& RAW,
                                const MuSigSeries& COR,
                                const MuSigSeries& CEA,
                                const std::string& outPng,
                                bool forceMuRange, double muYmin, double muYmax,
                                bool forceSgRange, double sgYmin, double sgYmax)
    {
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

        TCanvas c(isPhi ? "c3Way_MuSig_Phi" : "c3Way_MuSig_Eta",
                  "3WAY μ/σ vs E", 1000, 850);

        auto pads = makeEqualHeightTwinPads(
            c,
            std::string("ThreeVarMuSig_") + (isPhi ? "phi" : "eta"),
            0.15, 0.05, 0.14, 0.02, 0.04, 0.16
        );
        TPad* pTop = pads.first;
        TPad* pBot = pads.second;

        auto makeGE = [](const MuSigSeries& S, Color_t col, Style_t mk, bool useMu){
            const size_t N = S.E.size();
            std::vector<double> ex(N, 0.0);
            auto* g = new TGraphErrors(
                (int)N,
                const_cast<double*>(S.E.data()),
                const_cast<double*>((useMu ? S.mu : S.sg).data()),
                ex.data(),
                const_cast<double*>((useMu ? S.dmu : S.dsg).data())
            );
            g->SetMarkerStyle(mk);
            g->SetMarkerSize(1.05);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            return g;
        };

        const Color_t col = isPhi ? (kRed+1) : (kBlue+1);

        // ---------- μ(E)
        pTop->cd();
        double muMin = 0.0, muMax = 0.0;
        if (!forceMuRange)
        {
            double muLo = +1e30, muHi = -1e30;
            auto updMu = [&](const MuSigSeries& S){
                for (size_t i=0; i<S.mu.size(); ++i) {
                    const double lo = S.mu[i] - (S.dmu.empty()?0.0:S.dmu[i]);
                    const double hi = S.mu[i] + (S.dmu.empty()?0.0:S.dmu[i]);
                    muLo = std::min(muLo, lo);
                    muHi = std::max(muHi, hi);
                }
            };
            updMu(RAW); updMu(COR); updMu(CEA);
            if (muLo > muHi) { muLo = -1e-3; muHi = 1e-3; }
            const double absMax = std::max(std::fabs(muLo), std::fabs(muHi));
            const double minHalfRange = 3e-4;
            const double half = std::max(absMax * 1.10, minHalfRange);
            muMin = -half; muMax = +half;
        }
        else { muMin = muYmin; muMax = muYmax; }

        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(muMin);
        frU.SetMaximum(muMax);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

        // RAW=20 (closed), CORR=24 (open circle), CORR(EA)=25 (open square)
        TGraphErrors* gR = makeGE(RAW, col, 20, true);
        TGraphErrors* gC = makeGE(COR, col, 24, true);
        TGraphErrors* gA = makeGE(CEA, col, 25, true);
        gR->Draw("P SAME"); gC->Draw("P SAME"); gA->Draw("P SAME");

        // ---------- σ(E)
        pBot->cd();
        double sgMinFinal = 0.0, sgMaxFinal = 0.0;
        if (!forceSgRange)
        {
            double sgHi = 0.0;
            auto updSg = [&](const MuSigSeries& S){ for (double v : S.sg) sgHi = std::max(sgHi, v); };
            updSg(RAW); updSg(COR); updSg(CEA);
            if (sgHi <= 0) sgHi = 1.0;
            sgMinFinal = 0.0;
            sgMaxFinal = 1.15*sgHi;
        }
        else { sgMinFinal = sgYmin; sgMaxFinal = sgYmax; }

        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
        frL.SetMinimum(sgMinFinal);
        frL.SetMaximum(sgMaxFinal);
        frL.Draw("AXIS");

        TGraphErrors* hR = makeGE(RAW, col, 20, false);
        TGraphErrors* hC = makeGE(COR, col, 24, false);
        TGraphErrors* hA = makeGE(CEA, col, 25, false);
        hR->Draw("P SAME"); hC->Draw("P SAME"); hA->Draw("P SAME");

        // Legend in the μ(E) panel
        pTop->cd();
        TLegend lg(0.17, 0.73, 0.62, 0.86);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.040);
        lg.SetNColumns(3);
        lg.AddEntry(gR, "RAW-CP",     "p");
        lg.AddEntry(gC, "CP-corr",    "p");
        lg.AddEntry(gA, "CP-corr(EA)","p");
        lg.Draw();

        // Header
        c.cd();
        TLatex hh; hh.SetNDC(); hh.SetTextAlign(22); hh.SetTextSize(0.028);
        hh.DrawLatex(0.52,0.975,
                     Form("%s  -  Mean / RMS vs E  (Clusterizer CP family: RAW-CP, CP-corr, CP-corr(EA))",
                          isPhi ? "#Delta#phi" : "#Delta#eta"));

        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // --- Copy with μ(E) clamped to ±0.025 rad ---
    auto drawMuSigmaThree_muClamp = [&](bool isPhi,
                                        const MuSigSeries& RAW,
                                        const MuSigSeries& COR,
                                        const MuSigSeries& CEA,
                                        const std::string& outPng)
    {
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

        TCanvas c(isPhi ? "c3Way_MuSig_Phi_muClamp" : "c3Way_MuSig_Eta_muClamp",
                  "3WAY μ/σ vs E (μ clamp ±0.025)", 1000, 850);

        auto pads = makeEqualHeightTwinPads(
            c,
            std::string("ThreeVarMuSigClamp_") + (isPhi ? "phi" : "eta"),
            0.15, 0.05, 0.14, 0.02, 0.04, 0.16
        );
        TPad* pTop = pads.first;
        TPad* pBot = pads.second;

        auto makeGE = [](const MuSigSeries& S, Color_t col, Style_t mk, bool useMu){
            const size_t N = S.E.size();
            std::vector<double> ex(N, 0.0);
            auto* g = new TGraphErrors(
                (int)N,
                const_cast<double*>(S.E.data()),
                const_cast<double*>((useMu ? S.mu : S.sg).data()),
                ex.data(),
                const_cast<double*>((useMu ? S.dmu : S.dsg).data())
            );
            g->SetMarkerStyle(mk);
            g->SetMarkerSize(1.05);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            return g;
        };

        const Color_t col = isPhi ? (kRed+1) : (kBlue+1);

        // ---------- μ(E), clamped
        pTop->cd();
        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(-0.025);
        frU.SetMaximum(+0.025);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

        TGraphErrors* gR = makeGE(RAW, col, 20, true);
        TGraphErrors* gC = makeGE(COR, col, 24, true);
        TGraphErrors* gA = makeGE(CEA, col, 25, true);
        gR->Draw("P SAME"); gC->Draw("P SAME"); gA->Draw("P SAME");

        // Legend (top)
        pTop->cd();
        TLegend lg(0.37, 0.70, 0.92, 0.86);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.041);
        lg.SetNColumns(3);
        lg.AddEntry(gR, "RAW-CP",     "p");
        lg.AddEntry(gC, "CP-corr",    "p");
        lg.AddEntry(gA, "CP-corr(EA)","p");
        lg.Draw();

        // ---------- σ(E)
        pBot->cd();
        double sgHi = 0.0;
        auto updSg = [&](const MuSigSeries& S){ for (double v : S.sg) sgHi = std::max(sgHi, v); };
        updSg(RAW); updSg(COR); updSg(CEA);
        if (sgHi <= 0) sgHi = 1.0;

        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
        frL.SetMinimum(0.0);
        frL.SetMaximum(1.15*sgHi);
        frL.Draw("AXIS");

        TGraphErrors* hR = makeGE(RAW, col, 20, false);
        TGraphErrors* hC = makeGE(COR, col, 24, false);
        TGraphErrors* hA = makeGE(CEA, col, 25, false);
        hR->Draw("P SAME"); hC->Draw("P SAME"); hA->Draw("P SAME");

        // Header
        c.cd();
        TLatex hh; hh.SetNDC(); hh.SetTextAlign(22); hh.SetTextSize(0.028);
        hh.DrawLatex(0.52,0.975,
                     Form("%s  -  Mean / RMS vs E  (Clusterizer CP family: RAW-CP, CP-corr, CP-corr(EA))",
                          isPhi ? "#Delta#phi" : "#Delta#eta"));

        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // Build series for φ and η
    MuSigSeries RAW_phi = makeMuSigSeries(Hphi_cpRaw);
    MuSigSeries COR_phi = makeMuSigSeries(Hphi_cpCorr);
    MuSigSeries CEA_phi = makeMuSigSeries(Hphi_cpCorrEA);

    MuSigSeries RAW_eta = makeMuSigSeries(Heta_cpRaw);
    MuSigSeries COR_eta = makeMuSigSeries(Heta_cpCorr);
    MuSigSeries CEA_eta = makeMuSigSeries(Heta_cpCorrEA);

    // ---- compute shared μ(E) range across φ and η (symmetric, floor ±3e-4 rad, +10% headroom)
    auto muAbsMax3 = [&](const MuSigSeries& A, const MuSigSeries& B, const MuSigSeries& C)->double {
        double m = 0.0;
        auto upd = [&](const MuSigSeries& S){
            for (size_t i=0; i<S.mu.size(); ++i) {
                const double lo = S.mu[i] - (S.dmu.empty()?0.0:S.dmu[i]);
                const double hi = S.mu[i] + (S.dmu.empty()?0.0:S.dmu[i]);
                m = std::max(m, std::max(std::fabs(lo), std::fabs(hi)));
            }
        };
        upd(A); upd(B); upd(C);
        return m;
    };

    const double muAbsPhi = muAbsMax3(RAW_phi, COR_phi, CEA_phi);
    const double muAbsEta = muAbsMax3(RAW_eta, COR_eta, CEA_eta);
    const double muAbsAll = std::max(muAbsPhi, muAbsEta);
    const double muHalf   = std::max(muAbsAll * 1.10, 3e-4);
    const double muYmin   = -muHalf;
    const double muYmax   = +muHalf;

    // ---- compute shared σ(E) range: 0 → 1.15 × max(σ) across φ and η
    auto sgMax3 = [&](const MuSigSeries& A, const MuSigSeries& B, const MuSigSeries& C)->double {
        double s = 0.0;
        auto upd = [&](const MuSigSeries& S){ for (double v : S.sg) s = std::max(s, v); };
        upd(A); upd(B); upd(C);
        return s;
    };

    const double sgHiPhi = sgMax3(RAW_phi, COR_phi, CEA_phi);
    const double sgHiEta = sgMax3(RAW_eta, COR_eta, CEA_eta);
    const double sgYmin  = 0.0;
    const double sgYmax  = 1.15 * std::max(sgHiPhi, sgHiEta);

    // ---- draw with forced shared ranges
    drawMuSigmaThree(true,
                     RAW_phi, COR_phi, CEA_phi,
                     dir3Way + "/ThreeVariants_MeanSigmaVsE_Phi.png",
                     /*forceMuRange=*/true, muYmin, muYmax,
                     /*forceSgRange=*/true, sgYmin, sgYmax);

    drawMuSigmaThree(false,
                     RAW_eta, COR_eta, CEA_eta,
                     dir3Way + "/ThreeVariants_MeanSigmaVsE_Eta.png",
                     /*forceMuRange=*/true, muYmin, muYmax,
                     /*forceSgRange=*/true, sgYmin, sgYmax);

    // Parallel PNGs with μ(E) clamped to ±0.025 rad (top panel only)
    drawMuSigmaThree_muClamp(true,
                             RAW_phi, COR_phi, CEA_phi,
                             dir3Way + "/ThreeVariants_MeanSigmaVsE_Phi_muClamp_pm0p025.png");

    drawMuSigmaThree_muClamp(false,
                             RAW_eta, COR_eta, CEA_eta,
                             dir3Way + "/ThreeVariants_MeanSigmaVsE_Eta_muClamp_pm0p025.png");

    // ---- First non-empty energy bin overlays (exactly as requested)
    // φ: CP-corr vs CP-raw
    overlayTwoSeriesNthBin(true,  Hphi_cpRaw, Hphi_cpCorr,   "CP-corr",     kRed+1,
                           "Phi_cpCorr_vs_cpRaw", 1);
    // φ: CP-corr(EA) vs CP-raw
    overlayTwoSeriesNthBin(true,  Hphi_cpRaw, Hphi_cpCorrEA, "CP-corr(EA)", kRed+1,
                           "Phi_cpCorrEA_vs_cpRaw", 1);
    // η: CP-corr vs CP-raw
    overlayTwoSeriesNthBin(false, Heta_cpRaw, Heta_cpCorr,   "CP-corr",     kBlue+1,
                           "Eta_cpCorr_vs_cpRaw", 1);
    // η: CP-corr(EA) vs CP-raw
    overlayTwoSeriesNthBin(false, Heta_cpRaw, Heta_cpCorrEA, "CP-corr(EA)", kBlue+1,
                           "Eta_cpCorrEA_vs_cpRaw", 1);

    // ---- First-bin triplet overlays (RAW-CP / CP-corr / CP-corr(EA))
    overlayThreeSeriesNthBin(true,  Hphi_cpRaw, Hphi_cpCorr, Hphi_cpCorrEA, kRed+1,
                             "Phi_cpRaw_cpCorr_cpCorrEA", 1);
    overlayThreeSeriesNthBin(false, Heta_cpRaw, Heta_cpCorr, Heta_cpCorrEA, kBlue+1,
                             "Eta_cpRaw_cpCorr_cpCorrEA", 1);
}


// Template parameters let us accept your existing lambdas by perfect forwarding.
// No std::function, no indirection, no behavior change.
template <class DrawResidualPanelT,
          class MakeEqualPadsT,
          class StatsFromHistRangeT>
static void MakeFourWaySummaries(
    const char* outDir,
    const std::vector<std::pair<double,double>>& eEdges,
    double xMin, double xMax,

    // histogram sets
    const std::vector<TH1F*>& Hphi,
    const std::vector<TH1F*>& HphiCorr,
    const std::vector<TH1F*>& HphiScratchRaw,
    const std::vector<TH1F*>& HphiScratchCorr,
    const std::vector<TH1F*>& Heta,
    const std::vector<TH1F*>& HetaCorr,
    const std::vector<TH1F*>& HetaScratchRaw,
    const std::vector<TH1F*>& HetaScratchCorr,

    // helpers you already defined above inside your playground
    DrawResidualPanelT&&      drawResidualPanel,
    MakeEqualPadsT&&          makeEqualHeightTwinPads,
    StatsFromHistRangeT&&     statsFromHistRange)
{
    // =====================================================================
    // NEW: 4WAY summary overlays (first/second non-empty E bin) + singlets
    //       Outputs into: <outDir>/4WAYneededSUMMARIES
    // =====================================================================
    const std::string dir4Way = std::string(outDir) + "/4WAYneededSUMMARIES";
    ensureDir(dir4Way.c_str());

    // Find the N-th energy slice where BOTH hists in a pair are non-empty
    auto nthNonEmptyPair = [&](const std::vector<TH1F*>& A,
                               const std::vector<TH1F*>& B,
                               int nth)->int
    {
        int count = 0;
        for (std::size_t i=0; i<eEdges.size(); ++i) {
            if (i<A.size() && i<B.size() && A[i] && B[i] &&
                A[i]->Integral()>0 && B[i]->Integral()>0)
            {
                ++count;
                if (count == nth) return static_cast<int>(i);
            }
        }
        return -1;
    };

    // 2-series overlay for the N-th non-empty energy bin
    auto overlayTwoSeriesNthBin = [&](const char* baseLabel,   // "Variable-b (Benchmark)" or "Constant-b (Production)"
                                      Color_t color,           // kBlue+1 or kRed+1
                                      bool isPhi,              // true: phi, false: eta
                                      const std::vector<TH1F*>& Vraw,
                                      const std::vector<TH1F*>& Vcor,
                                      const std::string& outPngStemBase,
                                      int nth)
    {
        const int iE = nthNonEmptyPair(Vraw, Vcor, nth);
        if (iE < 0) {
            std::cerr << "[Playground][WARN] 4WAY overlay skipped for "
                      << baseLabel << (isPhi? " (phi)" : " (eta)")
                      << " - no non-empty bin #" << nth << ".\n";
            return;
        }

        const char* suffix =
            (nth == 1) ? "_FirstBin" :
            (nth == 2) ? "_SecondBin" : Form("_Bin%d", nth);

        TCanvas c(Form("c4Way_%s_%s_bin%d", baseLabel, isPhi? "phi" : "eta", nth),
                  "4WAY Nth-bin overlay", 1000, 750);

        const char* symShort = isPhi ? "#phi" : "#eta";
        std::vector<ResidualSpec> specs = {
            { Vraw[iE], Form("%s Uncorr, %s", baseLabel, symShort),  color, (Style_t)20, 1 },
            { Vcor[iE], Form("%s Corr, %s",   baseLabel, symShort),  color, (Style_t)24, 1 }
        };

        const char* residTitle = isPhi ? "#Delta#phi" : "#Delta#eta";
        const bool  isConstB   = (std::string(baseLabel).find("Constant-b") != std::string::npos);
        const char* familyTitle= isConstB ? "Constant b-Correction" : "Variable b-Correction";
        std::string header     = Form("%s Uncorrected/Corrected Overlay for %s", residTitle, familyTitle);

        drawResidualPanel(&c,
                          header,
                          eEdges[iE].first, eEdges[iE].second,
                          specs, xMin, xMax, 1.0);

        const std::string outPng = dir4Way + "/" + outPngStemBase + std::string(suffix) + ".png";
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // 1-series (stand-alone) plot(s) for the same N-th non-empty energy bin
    auto singleSeriesNthBin = [&](const char* baseLabel,    // same labels as overlay
                                  Color_t color,
                                  bool isPhi,
                                  const std::vector<TH1F*>& Vraw,
                                  const std::vector<TH1F*>& Vcor,
                                  const std::string& outPngStemBase,
                                  int nth)
    {
        const int iE = nthNonEmptyPair(Vraw, Vcor, nth);
        if (iE < 0) {
            std::cerr << "[Playground][WARN] 4WAY singlet skipped for "
                      << baseLabel << (isPhi? " (phi)" : " (eta)")
                      << " - no non-empty bin #" << nth << ".\n";
            return;
        }

        const char* suffix =
            (nth == 1) ? "_FirstBin" :
            (nth == 2) ? "_SecondBin" : Form("_Bin%d", nth);

        const char* symShort   = isPhi ? "#phi" : "#eta";
        const char* residTitle = isPhi ? "#Delta#phi" : "#Delta#eta";
        const bool  isConstB   = (std::string(baseLabel).find("Constant-b") != std::string::npos);
        const char* familyTitle= isConstB ? "Constant b-Correction" : "Variable b-Correction";

        auto drawOne = [&](TH1F* H, const char* tag, const char* stemTag, Style_t mstyle)
        {
            if (!H || H->Integral() <= 0) return;

            TCanvas c(Form("c4Way_%s_%s_%s_bin%d", baseLabel, isPhi? "phi":"eta", stemTag, nth),
                      "4WAY Nth-bin single", 1000, 750);

            std::vector<ResidualSpec> specs = {
                { H, Form("%s %s, %s", baseLabel, tag, symShort), color, mstyle, 1 }
            };

            std::string header = Form("%s %s Only for %s", residTitle, tag, familyTitle);

            drawResidualPanel(&c,
                              header,
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 1.0);

            const std::string outPng = dir4Way + "/" + outPngStemBase + "_" + stemTag + std::string(suffix) + ".png";
            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // Produce both stand-alone views
        drawOne(Vraw[iE], "Uncorr", "Uncorr", (Style_t)20);
        drawOne(Vcor[iE], "Corr",   "Corr",   (Style_t)24);
    };

    // ---- First and Second non-empty bin overlays + singlets ----

    // PDC (2×2 Block-Local), φ and η
    overlayTwoSeriesNthBin("Variable-b (Benchmark)", kBlue+1,  true,
                           HphiScratchRaw,  HphiScratchCorr,
                           "PDC_phi_RawVsCorr", 1);
    singleSeriesNthBin     ("Variable-b (Benchmark)", kBlue+1,  true,
                           HphiScratchRaw,  HphiScratchCorr,
                           "PDC_phi_RawVsCorr", 1);

    overlayTwoSeriesNthBin("Variable-b (Benchmark)", kBlue+1,  true,
                           HphiScratchRaw,  HphiScratchCorr,
                           "PDC_phi_RawVsCorr", 2);
    singleSeriesNthBin     ("Variable-b (Benchmark)", kBlue+1,  true,
                           HphiScratchRaw,  HphiScratchCorr,
                           "PDC_phi_RawVsCorr", 2);

    overlayTwoSeriesNthBin("Variable-b (Benchmark)", kBlue+1,  false,
                           HetaScratchRaw,  HetaScratchCorr,
                           "PDC_eta_RawVsCorr", 1);
    singleSeriesNthBin     ("Variable-b (Benchmark)", kBlue+1,  false,
                           HetaScratchRaw,  HetaScratchCorr,
                           "PDC_eta_RawVsCorr", 1);

    overlayTwoSeriesNthBin("Variable-b (Benchmark)", kBlue+1,  false,
                           HetaScratchRaw,  HetaScratchCorr,
                           "PDC_eta_RawVsCorr", 2);
    singleSeriesNthBin     ("Variable-b (Benchmark)", kBlue+1,  false,
                           HetaScratchRaw,  HetaScratchCorr,
                           "PDC_eta_RawVsCorr", 2);

    // CLUS (Constant-b (Production)), φ and η
    overlayTwoSeriesNthBin("Constant-b (Production)", kRed+1,  true,
                           Hphi, HphiCorr,
                           "CLUS_phi_RawVsCorr", 1);
    singleSeriesNthBin     ("Constant-b (Production)", kRed+1,  true,
                           Hphi, HphiCorr,
                           "CLUS_phi_RawVsCorr", 1);

    overlayTwoSeriesNthBin("Constant-b (Production)", kRed+1,  true,
                           Hphi, HphiCorr,
                           "CLUS_phi_RawVsCorr", 2);
    singleSeriesNthBin     ("Constant-b (Production)", kRed+1,  true,
                           Hphi, HphiCorr,
                           "CLUS_phi_RawVsCorr", 2);

    overlayTwoSeriesNthBin("Constant-b (Production)", kRed+1,  false,
                           Heta, HetaCorr,
                           "CLUS_eta_RawVsCorr", 1);
    singleSeriesNthBin     ("Constant-b (Production)", kRed+1,  false,
                           Heta, HetaCorr,
                           "CLUS_eta_RawVsCorr", 1);

    overlayTwoSeriesNthBin("Constant-b (Production)", kRed+1,  false,
                           Heta, HetaCorr,
                           "CLUS_eta_RawVsCorr", 2);
    singleSeriesNthBin     ("Constant-b (Production)", kRed+1,  false,
                           Heta, HetaCorr,
                           "CLUS_eta_RawVsCorr", 2);



    // ---- Four-variant μ/σ vs E summary (φ and η) to the same folder ----
    struct MuSigSeries { std::vector<double> E, mu, dmu, sg, dsg; };

    auto makeMuSigSeries = [&](const std::vector<TH1F*>& V)->MuSigSeries {
        MuSigSeries S;
        for (std::size_t i=0; i<eEdges.size(); ++i){
            if (i>=V.size() || !V[i] || V[i]->Integral()<=0) continue;
            double m, me, r, re;
            std::tie(m, me, r, re) = statsFromHistRange(V[i], xMin, xMax);
            S.E .push_back(0.5*(eEdges[i].first + eEdges[i].second));
            S.mu.push_back(m);  S.dmu.push_back(me);
            S.sg.push_back(r);  S.dsg.push_back(re);
        }
        return S;
    };

    // Supports forced y-ranges so φ and η plots can share identical axes.
    auto drawMuSigmaFour = [&](bool isPhi,
                               const MuSigSeries& CLUSraw,
                               const MuSigSeries& CLUScor,
                               const MuSigSeries& PDCraw,
                               const MuSigSeries& PDCcor,
                               const std::string& outPng,
                               bool forceMuRange, double muYmin, double muYmax,
                               bool forceSgRange, double sgYmin, double sgYmax)
    {
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

        TCanvas c(isPhi ? "c4Way_MuSig_Phi" : "c4Way_MuSig_Eta",
                  "4WAY μ/σ vs E", 1000, 850);

        auto pads = makeEqualHeightTwinPads(
            c,
            std::string("FourVarMuSig_") + (isPhi ? "phi" : "eta"),
            0.15, 0.05, 0.14, 0.02, 0.04, 0.16
        );
        TPad* pTop = pads.first;
        TPad* pBot = pads.second;

        auto makeGE = [](const MuSigSeries& S, Color_t col, Style_t mk, bool useMu){
            const size_t N = S.E.size();
            std::vector<double> ex(N, 0.0);
            auto* g = new TGraphErrors(
                (int)N,
                const_cast<double*>(S.E.data()),
                const_cast<double*>((useMu ? S.mu : S.sg).data()),
                ex.data(),
                const_cast<double*>((useMu ? S.dmu : S.dsg).data())
            );
            g->SetMarkerStyle(mk);
            g->SetMarkerSize(1.05);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            return g;
        };

        // ---------- μ(E)
        pTop->cd();
        double muMin = 0.0, muMax = 0.0;
        if (!forceMuRange)
        {
            double muLo = +1e30, muHi = -1e30;
            auto updMu = [&](const MuSigSeries& S){
                for (size_t i=0; i<S.mu.size(); ++i) {
                    const double lo = S.mu[i] - (S.dmu.empty()?0.0:S.dmu[i]);
                    const double hi = S.mu[i] + (S.dmu.empty()?0.0:S.dmu[i]);
                    muLo = std::min(muLo, lo);
                    muHi = std::max(muHi, hi);
                }
            };
            updMu(CLUSraw); updMu(CLUScor); updMu(PDCraw); updMu(PDCcor);
            if (muLo > muHi) { muLo = -1e-3; muHi = 1e-3; }

            // Symmetric around zero; enforce at least ±3e-4 rad
            const double absMax = std::max(std::fabs(muLo), std::fabs(muHi));
            const double minHalfRange = 3e-4; // ±0.3×10⁻³ rad
            const double half = std::max(absMax * 1.10, minHalfRange); // +10% headroom
            muMin = -half; muMax = +half;
        }
        else
        {
            muMin = muYmin; muMax = muYmax;
        }

        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(muMin);
        frU.SetMaximum(muMax);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

        // Styles: CLUS=red, PDC=blue; Uncorr=filled (20), Corr=open (24)
        TGraphErrors* gTCu = makeGE(CLUSraw, kRed+1,  20, true);
        TGraphErrors* gTCc = makeGE(CLUScor, kRed+1,  24, true);
        TGraphErrors* gBLu = makeGE(PDCraw,  kBlue+1, 20, true);
        TGraphErrors* gBLc = makeGE(PDCcor,  kBlue+1, 24, true);
        gBLu->Draw("P SAME"); gBLc->Draw("P SAME"); gTCu->Draw("P SAME"); gTCc->Draw("P SAME");

        // ---------- σ(E)
        pBot->cd();
        double sgMinFinal = 0.0, sgMaxFinal = 0.0;
        if (!forceSgRange)
        {
            double sgHi = 0.0;
            auto updSg = [&](const MuSigSeries& S){ for (double v : S.sg) sgHi = std::max(sgHi, v); };
            updSg(CLUSraw); updSg(CLUScor); updSg(PDCraw); updSg(PDCcor);
            if (sgHi <= 0) sgHi = 1.0;
            sgMinFinal = 0.0;
            sgMaxFinal = 1.15*sgHi;
        }
        else
        {
            sgMinFinal = sgYmin;
            sgMaxFinal = sgYmax;
        }

        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
        frL.SetMinimum(sgMinFinal);
        frL.SetMaximum(sgMaxFinal);
        frL.Draw("AXIS");

        TGraphErrors* hTCu = makeGE(CLUSraw, kRed+1,  20, false);
        TGraphErrors* hTCc = makeGE(CLUScor, kRed+1,  24, false);
        TGraphErrors* hBLu = makeGE(PDCraw,  kBlue+1, 20, false);
        TGraphErrors* hBLc = makeGE(PDCcor,  kBlue+1, 24, false);
        hBLu->Draw("P SAME"); hBLc->Draw("P SAME"); hTCu->Draw("P SAME"); hTCc->Draw("P SAME");

        // Legend in the μ(E) panel (upper-right), 2 columns
        pTop->cd();
        TLegend lg(0.17, 0.73, 0.65, 0.86);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.038);
        lg.SetNColumns(2);
        lg.AddEntry(gTCu, "Constant-b (Production) Uncorr", "p");
        lg.AddEntry(gBLu, "Variable-b (Benchmark) Uncorr",  "p");
        lg.AddEntry(gTCc, "Constant-b (Production) Corr",   "p");
        lg.AddEntry(gBLc, "Variable-b (Benchmark) Corr",    "p");
        lg.Draw();

        // Header
        c.cd();
        TLatex hh; hh.SetNDC(); hh.SetTextAlign(22); hh.SetTextSize(0.028);
        hh.DrawLatex(0.52,0.975,
                     Form("%s  -  Mean / RMS vs E  (Constant-b (Production) & Variable-b (Benchmark): Uncorr vs Corr)",
                          isPhi ? "#Delta#phi" : "#Delta#eta"));

        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    
    // --- Copy of drawMuSigmaFour, but μ(E) y-axis is clamped to ±0.025 rad ---
    auto drawMuSigmaFour_muClamp = [&](bool isPhi,
                                       const MuSigSeries& CLUSraw,
                                       const MuSigSeries& CLUScor,
                                       const MuSigSeries& PDCraw,
                                       const MuSigSeries& PDCcor,
                                       const std::string& outPng)
    {
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

        TCanvas c(isPhi ? "c4Way_MuSig_Phi_muClamp" : "c4Way_MuSig_Eta_muClamp",
                  "4WAY μ/σ vs E (μ clamp ±0.025)", 1000, 850);

        auto pads = makeEqualHeightTwinPads(
            c,
            std::string("FourVarMuSigClamp_") + (isPhi ? "phi" : "eta"),
            0.15, 0.05, 0.14, 0.02, 0.04, 0.16
        );
        TPad* pTop = pads.first;
        TPad* pBot = pads.second;

        auto makeGE = [](const MuSigSeries& S, Color_t col, Style_t mk, bool useMu){
            const size_t N = S.E.size();
            std::vector<double> ex(N, 0.0);
            auto* g = new TGraphErrors(
                (int)N,
                const_cast<double*>(S.E.data()),
                const_cast<double*>((useMu ? S.mu : S.sg).data()),
                ex.data(),
                const_cast<double*>((useMu ? S.dmu : S.dsg).data())
            );
            g->SetMarkerStyle(mk);
            g->SetMarkerSize(1.05);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            return g;
        };

        // ---------- μ(E)  (CLAMPED to ±0.025 rad) ----------
        pTop->cd();
        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(-0.025);
        frU.SetMaximum(+0.025);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

        // Styles: CLUS=red, PDC=blue; Uncorr=filled (20), Corr=open (24)
        TGraphErrors* gTCu = makeGE(CLUSraw, kRed+1,  20, true);
        TGraphErrors* gTCc = makeGE(CLUScor, kRed+1,  24, true);
        TGraphErrors* gBLu = makeGE(PDCraw,  kBlue+1, 20, true);
        TGraphErrors* gBLc = makeGE(PDCcor,  kBlue+1, 24, true);
        // draw order (blue first to keep red on top if overlapping)
        gBLu->Draw("P SAME"); gBLc->Draw("P SAME"); gTCu->Draw("P SAME"); gTCc->Draw("P SAME");

        // Put legend on the TOP pad (same as your revised version)
        pTop->cd();
        TLegend lg(0.37, 0.7, 0.92, 0.86);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.041);
        lg.SetNColumns(2);
        lg.AddEntry(gTCu, "Constant-b (Production) Uncorr", "p");
        lg.AddEntry(gBLu, "Variable-b (Benchmark) Uncorr",  "p");
        lg.AddEntry(gTCc, "Constant-b (Production) Corr",   "p");
        lg.AddEntry(gBLc, "Variable-b (Benchmark) Corr",    "p");
        lg.Draw();

        // ---------- σ(E) ----------
        pBot->cd();
        double sgHi = 0.0;
        auto updSg = [&](const MuSigSeries& S){ for (double v : S.sg) sgHi = std::max(sgHi, v); };
        updSg(CLUSraw); updSg(CLUScor); updSg(PDCraw); updSg(PDCcor);
        if (sgHi <= 0) sgHi = 1.0;

        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
        frL.SetMinimum(0.0);
        frL.SetMaximum(1.15*sgHi);
        frL.Draw("AXIS");

        TGraphErrors* hTCu = makeGE(CLUSraw, kRed+1,  20, false);
        TGraphErrors* hTCc = makeGE(CLUScor, kRed+1,  24, false);
        TGraphErrors* hBLu = makeGE(PDCraw,  kBlue+1, 20, false);
        TGraphErrors* hBLc = makeGE(PDCcor,  kBlue+1, 24, false);
        hBLu->Draw("P SAME"); hBLc->Draw("P SAME"); hTCu->Draw("P SAME"); hTCc->Draw("P SAME");

        // Header
        c.cd();
        TLatex hh; hh.SetNDC(); hh.SetTextAlign(22); hh.SetTextSize(0.028);
        hh.DrawLatex(0.52,0.975,
                     Form("%s  -  Mean / RMS vs E  (Constant-b (Production) & Variable-b (Benchmark): Uncorr vs Corr)",
                          isPhi ? "#Delta#phi" : "#Delta#eta"));

        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // Build series for φ and η (both first so we can harmonize y-ranges)
    MuSigSeries CLUSraw_phi = makeMuSigSeries(Hphi);
    MuSigSeries CLUScor_phi = makeMuSigSeries(HphiCorr);
    MuSigSeries PDCraw_phi  = makeMuSigSeries(HphiScratchRaw);
    MuSigSeries PDCcor_phi  = makeMuSigSeries(HphiScratchCorr);

    MuSigSeries CLUSraw_eta = makeMuSigSeries(Heta);
    MuSigSeries CLUScor_eta = makeMuSigSeries(HetaCorr);
    MuSigSeries PDCraw_eta  = makeMuSigSeries(HetaScratchRaw);
    MuSigSeries PDCcor_eta  = makeMuSigSeries(HetaScratchCorr);

    // ---- compute shared μ(E) range: symmetric, floor ±3e-4 rad, +10% headroom
    auto muAbsMax = [&](const MuSigSeries& A, const MuSigSeries& B,
                        const MuSigSeries& C, const MuSigSeries& D)->double {
      double m = 0.0;
      auto upd = [&](const MuSigSeries& S){
        for (size_t i=0; i<S.mu.size(); ++i) {
          const double lo = S.mu[i] - (S.dmu.empty()?0.0:S.dmu[i]);
          const double hi = S.mu[i] + (S.dmu.empty()?0.0:S.dmu[i]);
          m = std::max(m, std::max(std::fabs(lo), std::fabs(hi)));
        }
      };
      upd(A); upd(B); upd(C); upd(D);
      return m;
    };

    const double muAbsPhi = muAbsMax(CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi);
    const double muHalf   = std::max(muAbsPhi * 1.10, 3e-4);
    const double muYmin   = -muHalf;
    const double muYmax   = +muHalf;

    // ---- compute shared σ(E) range: 0 → 1.15 × max(σ) across φ and η
    auto sgMax = [&](const MuSigSeries& A, const MuSigSeries& B,
                     const MuSigSeries& C, const MuSigSeries& D)->double {
      double s = 0.0;
      auto upd = [&](const MuSigSeries& S){ for (double v : S.sg) s = std::max(s, v); };
      upd(A); upd(B); upd(C); upd(D);
      return s;
    };

    const double sgHiPhi = sgMax(CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi);
    const double sgHiEta = sgMax(CLUSraw_eta, CLUScor_eta, PDCraw_eta, PDCcor_eta);
    const double sgYmin  = 0.0;
    const double sgYmax  = 1.15 * std::max(sgHiPhi, sgHiEta);

    // ---- draw with forced shared ranges (matches the 12-arg lambda signature)
    drawMuSigmaFour(true,
                    CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi,
                    dir4Way + "/FourVariants_MeanSigmaVsE_Phi.png",
                    /*forceMuRange=*/true, muYmin, muYmax,
                    /*forceSgRange=*/true, sgYmin, sgYmax);

    drawMuSigmaFour(false,
                    CLUSraw_eta, CLUScor_eta, PDCraw_eta, PDCcor_eta,
                    dir4Way + "/FourVariants_MeanSigmaVsE_Eta.png",
                    /*forceMuRange=*/true, muYmin, muYmax,
                    /*forceSgRange=*/true, sgYmin, sgYmax);

    // Parallel PNGs with μ(E) clamped to ±0.025 rad (top panel only)
    drawMuSigmaFour_muClamp(true,
                            CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi,
                            dir4Way + "/FourVariants_MeanSigmaVsE_Phi_muClamp_pm0p025.png");

    drawMuSigmaFour_muClamp(false,
                            CLUSraw_eta, CLUScor_eta, PDCraw_eta, PDCcor_eta,
                            dir4Way + "/FourVariants_MeanSigmaVsE_Eta_muClamp_pm0p025.png");
}



void MakeDeltaPhiEtaPlayground(
    const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/PositionDep_sim_ALL.root",
    const char* outDir = "/Users/patsfan753/Desktop/scratchPDC",
    double      xMin   = -0.04,
    double      xMax   =  0.04)
{
    // -------------------- global plot style & folders --------------------
    gStyle->SetOptStat(0);
    ensureDir(outDir);
    
    const std::string dirOverlayNoCorr        = std::string(outDir) + "/etaPhiCompareNoCorrection";
    const std::string dirSeparateNoCorr       = std::string(outDir) + "/seperatePhiEtaNoCorrection";
    const std::string dirRawVsCorr            = std::string(outDir) + "/Clusterizer_RawVsCorr_EtaPhi";
    const std::string dirScratchRawVsCorr     = std::string(outDir) + "/fromScratch_RawVsCorr_EtaPhi";
    const std::string dirFSvsClusNoCorr       = std::string(outDir) + "/fromScratchVsClusterizerNoCorrection";
    const std::string dirFSvsClusWithCorr     = std::string(outDir) + "/fromScratchVsClusterizerWithCorrection";
    const std::string dirIndivPerVar          = std::string(outDir) + "/invididualPerVariantSummary";
    const std::string dirEtaPhi2D             = std::string(outDir) + "/2DetaDphiDep";
    const std::string dirFirstBinOverlaysOnly = std::string(outDir) + "/firstBinOverlaysOnly";
    ensureDir(dirOverlayNoCorr.c_str());
    ensureDir(dirSeparateNoCorr.c_str());
    ensureDir(dirRawVsCorr.c_str());
    ensureDir(dirScratchRawVsCorr.c_str());
    ensureDir(dirFSvsClusNoCorr.c_str());
    ensureDir(dirFSvsClusWithCorr.c_str());
    ensureDir(dirIndivPerVar.c_str());
    ensureDir(dirEtaPhi2D.c_str());
    ensureDir(dirFirstBinOverlaysOnly.c_str());

    // -------------------- input file --------------------
    std::cout << "[Playground] Opening: " << inFile << "\n";
    std::unique_ptr<TFile> fIn(TFile::Open(inFile,"READ"));
    if (!fIn || fIn->IsZombie()){
        std::cerr << "[Playground][ERROR] Cannot open file: " << inFile << "\n";
        return;
    }

    // -------------------- binning mode detection (as in PDCanalysis) --------------------
    EBinningMode binMode = EBinningMode::kRange;
    if (fIn->Get("h3_blockCoord_E_range") == nullptr &&
        fIn->Get("h3_blockCoord_E_disc")  != nullptr) {
        binMode = EBinningMode::kDiscrete;
        std::cout << "[Playground] Detected DISCRETE energy binning\n";
    } else {
        std::cout << "[Playground] Detected RANGE energy binning\n";
    }

    std::vector<std::pair<double,double>> eEdges;
    eEdges.reserve(N_E);
    for (int i = 0; i < N_E; ++i)
        eEdges.emplace_back(E_edges[i], E_edges[i+1]);
    
    // tag builder identical to PDCanalysis
    auto makeSliceTag = [binMode](std::size_t /*iE*/, const std::pair<double,double>& edge){
        return (binMode == EBinningMode::kRange)
               ? TString::Format("%.0f_%.0f", edge.first, edge.second)
               : TString::Format("E%.0f",     edge.first);
    };

    // =====================================================================
    //                      micro-helpers
    // =====================================================================

    auto setPadMargins = [](double L=0.12,double R=0.04,double T=0.12,double B=0.12){
        gPad->SetLeftMargin(L); gPad->SetRightMargin(R);
        gPad->SetTopMargin(T);  gPad->SetBottomMargin(B);
    };

    auto cloneCountsHist = [](TH1F* h)->TH1F*{
        if (!h) return nullptr;
        auto* c = static_cast<TH1F*>(h->Clone());
        c->SetDirectory(nullptr);
        return c;
    };

    auto styleCountsHist = [](TH1F* h, Color_t col, Style_t mk=20, double msz=0.7){
        if (!h) return;
        h->SetMarkerStyle(mk);
        h->SetMarkerColor(col);
        h->SetLineColor  (col);
        h->SetMarkerSize (msz);
    };

    auto computeYmaxCounts = [](const std::initializer_list<const TH1F*>& hh)->double{
        double yMax = 0.0;
        for (auto* h : hh){
            if (!h) continue;
            for (int b=1; b<=h->GetNbinsX(); ++b)
                yMax = std::max(yMax, h->GetBinContent(b) + h->GetBinError(b));
        }
        return yMax;
    };

    auto drawNDCEnergy = [](double x, double y, double eLo, double eHi, double ts=0.040){
        TLatex t; t.SetNDC(); t.SetTextSize(ts); t.SetTextAlign(33);
        t.DrawLatex(x,y,Form("[%.0f, %.0f) GeV", eLo, eHi));
    };

    auto makeNDCLegend = [](double x1,double y1,double x2,double y2)->TLegend*{
        auto* lg = new TLegend(x1,y1,x2,y2,"","NDC");
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
        return lg;
    };

    auto drawHeaderNDC = [](const char* txt, double y=0.975, double ts=0.045){
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(ts);
        h.DrawLatex(0.5, y, txt);
    };

    auto styleFit = [](TF1* f, Color_t col, int lw=2, int npx=400){
        if (!f) return;
        f->SetLineColor(col); f->SetLineWidth(lw); f->SetNpx(npx);
    };

    // Build two vertically stacked pads whose *inner drawable heights* are identical.
    // L/R are shared margins. TT/BT are top-pad top/bottom margins. TB/BB are bottom-pad top/bottom margins.
    // Returns {pTop, pBot}.
    auto makeEqualHeightTwinPads = [&](TCanvas& c,
                                       const std::string& stem,
                                       double L = 0.15, double R = 0.05,
                                       double TT = 0.14, double BT = 0.02,   // top pad margins
                                       double TB = 0.04, double BB = 0.16)   // bottom pad margins
        -> std::pair<TPad*,TPad*>
    {
        const double hTop   = 1.0 - TT - BT;            // inner height fraction (top)
        const double hBot   = 1.0 - TB - BB;            // inner height fraction (bottom)
        const double ySplit = hTop / (hTop + hBot);     // bottom pad height in NDC

        const std::string nmTop = stem + "_top";
        const std::string nmBot = stem + "_bot";

        if (auto* oldTop = dynamic_cast<TPad*>(gROOT->FindObject(nmTop.c_str()))) oldTop->Delete();
        if (auto* oldBot = dynamic_cast<TPad*>(gROOT->FindObject(nmBot.c_str()))) oldBot->Delete();

        TPad* pBot = new TPad(nmBot.c_str(), "", 0.0, 0.00, 1.0, ySplit);
        TPad* pTop = new TPad(nmTop.c_str(), "", 0.0, ySplit, 1.0, 1.0);

        pTop->SetLeftMargin(L); pTop->SetRightMargin(R);
        pTop->SetTopMargin(TT); pTop->SetBottomMargin(BT);
        pBot->SetLeftMargin(L); pBot->SetRightMargin(R);
        pBot->SetTopMargin(TB); pBot->SetBottomMargin(BB);

        pTop->SetFillColor(0); pBot->SetFillColor(0);
        pTop->SetBorderMode(0); pBot->SetBorderMode(0);

        pTop->Draw();
        pBot->Draw();
        c.Modified(); c.Update();

        return std::make_pair(pTop, pBot);
    };




    // ---- Range-aware (xMin,xMax) Mean/RMS using ROOT built-ins ----
    auto statsFromHistRange = [&](TH1F* h, double lo, double hi)
        -> std::tuple<double,double,double,double>
    {
        if (!h || h->Integral() <= 0) return std::make_tuple(0.0,0.0,0.0,0.0);
        TAxis* ax = h->GetXaxis();
        const int nbin    = ax->GetNbins();
        const int oldLo   = ax->GetFirst();
        const int oldHi   = ax->GetLast();
        const int iLo     = std::max(1,       ax->FindFixBin(lo + 1e-9));
        const int iHi     = std::min(nbin,    ax->FindFixBin(hi - 1e-9));
        if (iLo > iHi)    return std::make_tuple(0.0,0.0,0.0,0.0);
        ax->SetRange(iLo, iHi);
        const double m    = h->GetMean();
        const double me   = h->GetMeanError();
        const double r    = h->GetRMS();
        const double re   = h->GetRMSError();
        ax->SetRange(oldLo, oldHi);
        return std::make_tuple(m,me,r,re);
    };
    

    auto drawResidualPanel = [&](TVirtualPad* pad,
                                 const std::string& plotName,
                                 double eLo, double eHi,
                                 const std::vector<ResidualSpec>& specs,
                                 double xMin, double xMax,
                                 double textScale = 1.0)
    {
        if (!pad || specs.empty()) return;
        pad->cd();

        // Consistent margins for all panels (increase bottom margin so x-axis title doesn't clip)
        setPadMargins(0.14, 0.06, 0.06, 0.16);

        // Clone/style hists; compute yMax from points and from fitted curves
        std::vector<TH1F*> clones;
        clones.reserve(specs.size());
        double yMax = 0.0;

        for (const auto& s : specs) {
            if (!s.h || s.h->Integral() == 0) continue;
            TH1F* c = cloneCountsHist(s.h);
            styleCountsHist(c, s.col, s.mk, 1.0);
            clones.push_back(c);
            yMax = std::max(yMax, computeYmaxCounts({c}));
        }
        if (clones.empty()) return;

        // Decide x-axis title based on what we're overlaying:
        //   • Only φ  -> "#Delta#phi  [rad]"
        //   • Only η  -> "#Delta#eta  [rad]"
        //   • Both    -> "#Delta#phi/#eta  [rad]"
        bool sawPhi = false, sawEta = false;
        for (const auto& s : specs) {
            // check legend label if present
            if (s.label) {
                TString lab(s.label); lab.ToLower();
                if (lab.Contains("#phi") || lab.Contains("phi")) sawPhi = true;
                if (lab.Contains("#eta") || lab.Contains("eta")) sawEta = true;
            }
            // also inspect histogram name as a fallback
            if (s.h) {
                TString hn(s.h->GetName()); hn.ToLower();
                if (hn.Contains("phi")) sawPhi = true;
                if (hn.Contains("eta")) sawEta = true;
            }
        }
        const char* xTitle = "#Delta  [rad]";
        if (sawPhi && sawEta)      xTitle = "#Delta#phi/#eta  [rad]";
        else if (sawPhi)           xTitle = "#Delta#phi  [rad]";
        else if (sawEta)           xTitle = "#Delta#eta  [rad]";

        // Draw histograms with tuned axis cosmetics
        TH1F* base = clones.front();
        base->SetTitle("");
        base->GetXaxis()->SetTitle(xTitle);
        base->GetYaxis()->SetTitle("Counts");
        base->GetXaxis()->SetRangeUser(xMin, xMax);
        base->GetYaxis()->SetRangeUser(0.0, 1.20*yMax);

        // Make sure the x-axis title sits fully inside the pad
        base->GetXaxis()->SetTitleSize(0.045 * textScale);
        base->GetXaxis()->SetLabelSize(0.035 * textScale);
        base->GetXaxis()->SetTitleOffset(1.10);

        base->Draw("E");
        for (std::size_t i=1; i<clones.size(); ++i) clones[i]->Draw("E SAME");

        // Top-left: Plot name + energy
        TLatex head; head.SetNDC(); head.SetTextFont(42);
        head.SetTextAlign(13); head.SetTextSize(0.033 * textScale);
        head.DrawLatex(0.17, 0.89, Form("%s, [%.0f, %.0f) GeV", plotName.c_str(), eLo, eHi));

        // Top-right legend (points only, no error bars) - keep alive on the pad
        TLegend* lg = new TLegend(0.16, 0.75, 0.3, 0.83, "", "NDC");
        lg->SetBorderSize(0);
        lg->SetFillStyle(0);
        lg->SetTextSize(0.028 * textScale);
        lg->SetTextColor(kBlack);
        for (std::size_t i=0; i<clones.size(); ++i)
            lg->AddEntry(clones[i], specs[i].label, "p");
        lg->Draw();
        gPad->Update();

        pad->Modified(); pad->Update();
    };


    // =====================================================================
    //                         data loaders + fitter
    // =====================================================================

    auto loadSet = [&](const char* hPat) -> std::vector<TH1F*>
    {
        std::vector<TH1F*> H(eEdges.size(), nullptr);
        for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
        {
            const TString tag = makeSliceTag(iE, eEdges[iE]);
            TH1F* h = dynamic_cast<TH1F*>( fIn->Get(Form(hPat, tag.Data())) );
            if (!h)
                std::cerr << "[Playground][WARN] Missing hist for slice " << tag
                          << " (" << hPat << ")\n";
            H[iE] = h;
        }
        return H;
    };

    // =====================================================================
    //                    per‑residual producer (keeps outputs)
    // =====================================================================
    auto produceResidualSet = [&](const char* stem,          // "DeltaPhi" or "DeltaEta"
                                  const char* hPat,          // "h_phi_diff_cpRaw_%s" / "h_eta_diff_cpRaw_%s"
                                  const char* sym)           // "#Delta#phi" / "#Delta#eta"
    {
        const int v = 2;     // "no corr, cluster" (for colors/markers in μ/σ plots)
        auto H = loadSet(hPat);

        // ---------- PNG #1 : First-bin single plot (unified panel drawer)
        if (!H.empty() && H[0] && H[0]->Integral()>0)
        {
            TCanvas c1(Form("cPlay_%s_First",stem),
                       Form("%s first slice (Clusterizer Output, No Position Correction)",stem),1000,750);

            const bool isPhi = (TString(stem) == "DeltaPhi");
            const int  col   = isPhi ? (kRed+1) : (kBlue+1);

            std::vector<ResidualSpec> specs = {
                { H[0], sym, (Color_t)col, (Style_t)20, 1 }
            };

            drawResidualPanel(&c1, "Clusterizer Output, No Position Correction",
                              eEdges.front().first, eEdges.front().second,
                              specs, xMin, xMax, 1.0);

            std::string out1 = dirSeparateNoCorr + "/" + std::string(stem) + "_noCorrCluster_FirstBin.png";
            c1.SaveAs(out1.c_str());
            c1.Close();
            std::cout << "[Playground] Wrote " << out1 << "\n";
        }
        else {
            std::cerr << "[Playground][ERROR] First energy‑bin histogram missing/empty for " << stem
                      << " - skipped PNG #1.\n";
        }


        // ---------- PNG #2 : 4×2 table (this residual only) ----------
        {
            const int nCol=4, nRow=2, maxPads=nCol*nRow;
            TCanvas cOT(Form("cPlay_%s_Table", stem),
                        Form("%s - all slices (Clusterizer Output, No Position Correction)", sym),
                        1600,900);
            cOT.SetTopMargin(0.10);
            cOT.Divide(nCol,nRow,0,0);

            const bool isPhi = (TString(stem) == "DeltaPhi");
            const int  col   = isPhi ? (kRed+1) : (kBlue+1);

            int pad=1;
            for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE)
            {
                if (!H[iE] || H[iE]->Integral()==0) { ++pad; continue; }

                cOT.cd(pad++);
                setPadMargins();

                std::vector<ResidualSpec> specs = {
                    { H[iE], sym, (Color_t)col, (Style_t)20, 1 }
                };

                // Use a slightly smaller text scale on table pads
                drawResidualPanel(gPad, "Clusterizer Output, No Position Correction",
                                  eEdges[iE].first, eEdges[iE].second,
                                  specs, xMin, xMax, 0.85);
            }

            cOT.cd(0);
            drawHeaderNDC(Form("%s  -  Clusterizer Output, No Position Correction  (all energy bins)", sym), 0.975, 0.045);

            std::string outOT = dirSeparateNoCorr + "/" + std::string(stem) + "_noCorrCluster_AllBins_Table.png";
            cOT.SaveAs(outOT.c_str());
            cOT.Close();
            std::cout << "[Playground] Wrote " << outOT << "\n";
        }

        // ---------- PNG #3 : μ(E) & σ(E) summary ----------
        {
            std::vector<double> eCtr; eCtr.reserve(eEdges.size());
            std::vector<double> mu, dmu, sg, dsg;

            for (std::size_t iE=0; iE<eEdges.size(); ++iE)
            {
                eCtr.push_back( 0.5*(eEdges[iE].first + eEdges[iE].second) );

                if (!H[iE] || H[iE]->Integral()==0){
                    mu .push_back(0.0); dmu.push_back(0.0);
                    sg .push_back(0.0); dsg.push_back(0.0);
                    continue;
                }

                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(H[iE], xMin, xMax);
                mu .push_back(m );  dmu.push_back(me);
                sg .push_back(r );  dsg.push_back(re);

            }

            if (!eCtr.empty())
            {
                const double xAxisMin = 0.0;
                const double xAxisMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
                std::vector<double> ex(eCtr.size(), 0.0);

                TCanvas cS(Form("cPlay_%s_MuSig",stem),
                           Form("%s #mu/#sigma_{RMS} vs E (Clusterizer Output, No Position Correction)",sym),
                           900,800);
                auto pads = makeEqualHeightTwinPads(cS, "MuSigSingle", 0.15,0.05, 0.14,0.02, 0.04,0.16);
                TPad* pTop = pads.first;
                TPad* pBot = pads.second;

                // μ(E)
                pTop->cd();
                gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.04); // keep your style
                double muLo=1e30, muHi=-1e30;
                for (std::size_t i=0;i<mu.size();++i){
                    muLo = std::min(muLo, mu[i]-dmu[i]);
                    muHi = std::max(muHi, mu[i]+dmu[i]);
                }
                const double padMu = 0.2*(muHi-muLo);
                TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
                frU.SetMinimum(muLo - padMu);
                frU.SetMaximum(muHi + padMu);
                frU.Draw("AXIS");
                TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

                TGraphErrors gMu(eCtr.size(), eCtr.data(), mu.data(), ex.data(), dmu.data());
                gMu.SetMarkerStyle(kMk[v]); gMu.SetMarkerColor(kCol[v]);
                gMu.SetLineColor  (kCol[v]); gMu.SetMarkerSize(1.1);
                gMu.Draw("P SAME");

                // σ(E)
                pBot->cd();
                gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);     // keep your style
                double sgHi=-1e30; for (double s : sg) sgHi = std::max(sgHi,s);
                TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
                frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

                TGraphErrors gSg(eCtr.size(), eCtr.data(), sg.data(), ex.data(), dsg.data());
                gSg.SetMarkerStyle(kMk[v]); gSg.SetMarkerColor(kCol[v]);
                gSg.SetLineColor  (kCol[v]); gSg.SetMarkerSize(1.1);
                gSg.Draw("P SAME");

                // Title + save
                cS.cd(0);
                TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.038);
                h.DrawLatex(0.5,0.98,Form("%s  (Clusterizer Output, No Position Correction)  -  Mean / RMS vs E", sym));

                std::string out3 = dirSeparateNoCorr + "/" + std::string(stem) + "_noCorrCluster_MeanSigmaVsE.png";
                cS.SaveAs(out3.c_str());
                cS.Close();
                std::cout << "[Playground] Wrote " << out3 << "\n";

                // CSV: single series (“no corr, cluster”), residual from stem
                const std::string resid = (TString(stem)=="DeltaPhi") ? "phi" : "eta";
                saveMuSigmaCSV(out3,
                               { "no corr, cluster" },         // variants
                               { resid },                       // residuals
                               { eCtr },                        // E centers
                               { mu  }, { dmu },                // μ, dμ
                               { sg  }, { dsg });               // σ, dσ

            }
        }

        return H;    // for overlays later
    };

    // =====================================================================
    //                    Produce Δφ and Δη sets (unchanged outputs)
    // =====================================================================
    auto Hphi = produceResidualSet("DeltaPhi", "h_phi_diff_cpRaw_%s", "#Delta#phi");
    auto Heta = produceResidualSet("DeltaEta", "h_eta_diff_cpRaw_%s", "#Delta#eta");

    // Corrected (CorrectPosition, cluster) sets
    auto HphiCorr = loadSet("h_phi_diff_cpCorr_%s");
    auto HetaCorr = loadSet("h_eta_diff_cpCorr_%s");
    
    auto Hphi_cpCorrEA = loadSet("h_phi_diff_cpCorrEA_%s");
    auto Heta_cpCorrEA = loadSet("h_eta_diff_cpCorrEA_%s");

    // -------------------- from-scratch (PDC) sets --------------------
    // "no corr, scratch"  vs  "b(E) corr, scratch"
    //   φ : h_phi_diff_raw_%s   vs  h_phi_diff_corr_%s
    //   η : h_eta_diff_raw_%s   vs  h_eta_diff_corr_%s
    auto HphiScratchRaw  = loadSet("h_phi_diff_raw_%s");
    auto HphiScratchCorr = loadSet("h_phi_diff_corr_%s");
    auto HetaScratchRaw  = loadSet("h_eta_diff_raw_%s");
    auto HetaScratchCorr = loadSet("h_eta_diff_corr_%s");

    // =========================================================================================
    // NEW: First-bin overlays (four PNGs) into firstBinOverlaysOnly
    //       • PDC (Variable-b (Benchmark)): blue;   RAW = closed circle (20), CORR = open circle (24)
    //       • CLUS (Constant-b (Production)):   red;    RAW = closed circle (20), CORR = open circle (24)
    //       • Legend labels include residual tag (", #phi" or ", #eta")
    // =========================================================================================
    if (!eEdges.empty()) {
        const double eLo = eEdges.front().first;
        const double eHi = eEdges.front().second;

        // Small helper to draw a two-series first-bin overlay with the standardized styles/labels.
        auto drawFirstBinTwoSeries = [&](TH1F* hRaw, TH1F* hCorr,
                                         const char* variantBase,  // "Variable-b (Benchmark)" or "Constant-b (Production)"
                                         const char* residLatex,   // "#phi" or "#eta"
                                         Color_t color,
                                         const std::string& outPng)
        {
            if (!hRaw || !hCorr) return;
            if (hRaw->Integral() <= 0 || hCorr->Integral() <= 0) return;

            TCanvas c(Form("cFB_%s_%s", variantBase, residLatex), "first-bin overlay", 1000, 750);

            // Two series with required styles and legend text
            std::vector<ResidualSpec> specs = {
                { hRaw,  Form("%s Uncorr, %s", variantBase, residLatex), color, (Style_t)20, 1 }, // RAW: closed
                { hCorr, Form("%s Corr, %s",   variantBase, residLatex), color, (Style_t)24, 1 }  // CORR: open
            };

            drawResidualPanel(&c,
                              "First-bin Overlay (raw vs corr)",
                              eLo, eHi,
                              specs,
                              xMin, xMax,
                              1.0);

            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // Paths
        const std::string pngPDCphi = std::string(dirFirstBinOverlaysOnly) + "/PDC_RawVsCorr_Phi_FirstBin.png";
        const std::string pngPDCeta = std::string(dirFirstBinOverlaysOnly) + "/PDC_RawVsCorr_Eta_FirstBin.png";
        const std::string pngCLUSphi= std::string(dirFirstBinOverlaysOnly) + "/CLUS_RawVsCorr_Phi_FirstBin.png";
        const std::string pngCLUSeta= std::string(dirFirstBinOverlaysOnly) + "/CLUS_RawVsCorr_Eta_FirstBin.png";

        // Draw the four overlays (guarded by presence of first-bin hists)
        if (!HphiScratchRaw.empty() && !HphiScratchCorr.empty()
            && HphiScratchRaw[0] && HphiScratchCorr[0])
        {
            drawFirstBinTwoSeries(HphiScratchRaw[0], HphiScratchCorr[0],
                                  "Variable-b (Benchmark)", "#phi", kBlue+1, pngPDCphi);
        }
        if (!HetaScratchRaw.empty() && !HetaScratchCorr.empty()
            && HetaScratchRaw[0] && HetaScratchCorr[0])
        {
            drawFirstBinTwoSeries(HetaScratchRaw[0], HetaScratchCorr[0],
                                  "Variable-b (Benchmark)", "#eta", kBlue+1, pngPDCeta);
        }
        if (!Hphi.empty() && !HphiCorr.empty()
            && Hphi[0] && HphiCorr[0])
        {
            drawFirstBinTwoSeries(Hphi[0], HphiCorr[0],
                                  "Constant-b (Production)", "#phi", kRed+1, pngCLUSphi);
        }
        if (!Heta.empty() && !HetaCorr.empty()
            && Heta[0] && HetaCorr[0])
        {
            drawFirstBinTwoSeries(Heta[0], HetaCorr[0],
                                  "Constant-b (Production)", "#eta", kRed+1, pngCLUSeta);
        }
    }
    
    // =====================================================================
    // NEW: Headline σ(Δφ) vs E and σ(Δη) vs E with four curves
    //       • Constant-b (Production):    CLUS-RAW  = Hphi/Heta,        CLUS-CORR  = HphiCorr/HetaCorr
    //       • 2×2 block-local: PDC-RAW   = HphiScratchRaw/HetaScratchRaw
    //                           PDC-CORR  = HphiScratchCorr/HetaScratchCorr
    // Saves directly to outDir:
    //   - Headline_Sigma_DeltaPhi_vsE.png
    //   - Headline_Sigma_DeltaEta_vsE.png
    // =====================================================================
    struct SigSeries { std::vector<double> E, S, dS; };

    auto buildSigmaSeries = [&](const std::vector<TH1F*>& V)->SigSeries {
        SigSeries S;
        S.E.reserve(V.size()); S.S.reserve(V.size()); S.dS.reserve(V.size());
        for (std::size_t i=0; i<V.size(); ++i) {
            TH1F* h = V[i];
            if (!h || h->Integral()<=0) continue;
            double m, me, r, re;
            std::tie(m, me, r, re) = statsFromHistRange(h, xMin, xMax); // same range as rest of code
            S.E .push_back(0.5*(eEdges[i].first + eEdges[i].second));
            S.S .push_back(r);
            S.dS.push_back(re);
        }
        return S;
    };

    // Utility to draw the four-series σ(E) overlay for a given residual
    auto drawHeadlineSigma = [&](const SigSeries& sCLUSraw,
                                 const SigSeries& sCLUScor,
                                 const SigSeries& sPDCraw,
                                 const SigSeries& sPDCcor,
                                 const char* yTitle,
                                 const char* pngName)
    {
        // x-range from binning
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0
                          : eEdges.back().second; // upper edge of last slice

        // y-range (max across all four)
        auto ymax = [&](const SigSeries& SS)->double {
            double v=0.0;
            for (std::size_t i=0;i<SS.S.size();++i)
                v = std::max(v, SS.S[i] + (SS.dS.empty()?0.0:SS.dS[i]));
            return v;
        };
        double yMax = std::max( std::max(ymax(sCLUSraw), ymax(sCLUScor)),
                                std::max(ymax(sPDCraw),  ymax(sPDCcor)) );
        if (yMax <= 0.0) yMax = 1.0;

        // Build graphs
        auto makeGE = [](const SigSeries& SS, Color_t col, Style_t mk){
            auto* g = new TGraphErrors((int)SS.E.size(),
                                       const_cast<double*>(SS.E.data()),
                                       const_cast<double*>(SS.S.data()),
                                       nullptr,
                                       const_cast<double*>(SS.dS.data()));
            g->SetMarkerStyle(mk);
            g->SetMarkerSize(1.05);
            g->SetMarkerColor(col);
            g->SetLineColor(col);
            return g;
        };

        // Styling convention:
        //   Constant-b (Production) (CLUS) = red;    PDC (block-local) = blue
        //   Uncorr = closed circle (20);  Corr = open circle (24)
        TGraphErrors* gCLUSraw = makeGE(sCLUSraw, kRed+1,  20);
        TGraphErrors* gCLUScor = makeGE(sCLUScor, kRed+1,  24);
        TGraphErrors* gPDCraw  = makeGE(sPDCraw,  kBlue+1, 20);
        TGraphErrors* gPDCcor  = makeGE(sPDCcor,  kBlue+1, 24);

        // Canvas & frame  (increase bottom margin + explicit title/label sizes/offset)
        TCanvas c("cHeadlineSigma","Headline σ(E)",900,600);
        c.SetLeftMargin(0.15);
        c.SetRightMargin(0.06);
        c.SetTopMargin(0.10);
        c.SetBottomMargin(0.18); // was 0.13: prevents x-axis title clipping

        TH1F fr("fr","",1,xMinE,xMaxE);
        fr.SetTitle("");

        // Axis titles and labels
        fr.GetXaxis()->SetTitle("E  [GeV]");
        fr.GetYaxis()->SetTitle(yTitle); // "#sigma_{RMS}  [rad]"

        // Ensure the x-axis title fits inside the pad
        fr.GetXaxis()->SetTitleSize(0.050);
        fr.GetXaxis()->SetLabelSize(0.038);
        fr.GetXaxis()->SetTitleOffset(1.05); // mild offset keeps it clear but inside the pad

        // Keep y-axis sizing consistent
        fr.GetYaxis()->SetTitleSize(0.050);
        fr.GetYaxis()->SetLabelSize(0.038);

        fr.SetMinimum(0.0);
        fr.SetMaximum(1.15*yMax);
        fr.Draw("AXIS");
        // Draw
        gPDCraw ->Draw("P SAME");
        gPDCcor ->Draw("P SAME");
        gCLUSraw->Draw("P SAME");
        gCLUScor->Draw("P SAME");

        // Legend with requested labels
        TLegend lg(0.165,0.185,0.64,0.33);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.033);
        lg.SetNColumns(2);

        // Row 1
        lg.AddEntry(gCLUSraw, "Constant-b (Production) Uncorr",    "p"); // CLUS-RAW
        lg.AddEntry(gPDCraw,  "Variable-b (Benchmark) Uncorr", "p"); // PDC-RAW
        // Row 2
        lg.AddEntry(gCLUScor, "Constant-b (Production) Corr",      "p"); // CLUS-CORR
        lg.AddEntry(gPDCcor,  "Variable-b (Benchmark) Corr",   "p"); // PDC-CORR

        lg.Draw();

        // Subheader (kept subtle)
        TLatex h; h.SetNDC(); h.SetTextAlign(13); h.SetTextSize(0.045);
        h.DrawLatex(0.4,0.955,"Residual RMS vs energy");

        const std::string outPng = std::string(outDir) + "/" + pngName;
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // Build σ-series for the four variants, for φ and η
    SigSeries sPhi_CLUSraw  = buildSigmaSeries(Hphi);
    SigSeries sPhi_CLUScor  = buildSigmaSeries(HphiCorr);
    SigSeries sPhi_PDCraw   = buildSigmaSeries(HphiScratchRaw);
    SigSeries sPhi_PDCcor   = buildSigmaSeries(HphiScratchCorr);

    SigSeries sEta_CLUSraw  = buildSigmaSeries(Heta);
    SigSeries sEta_CLUScor  = buildSigmaSeries(HetaCorr);
    SigSeries sEta_PDCraw   = buildSigmaSeries(HetaScratchRaw);
    SigSeries sEta_PDCcor   = buildSigmaSeries(HetaScratchCorr);

    // Draw & save the two headline PNGs in the base output directory
    drawHeadlineSigma(sPhi_CLUSraw, sPhi_CLUScor, sPhi_PDCraw, sPhi_PDCcor,
                      "#sigma_{RMS}(#Delta#phi)  [rad]",
                      "Headline_Sigma_DeltaPhi_vsE.png");

    drawHeadlineSigma(sEta_CLUSraw, sEta_CLUScor, sEta_PDCraw, sEta_PDCcor,
                      "#sigma_{RMS}(#Delta#eta)  [rad]",
                      "Headline_Sigma_DeltaEta_vsE.png");

    // =====================================================================
    // NEW: Ratio plots (corrected Constant-b (Production) / corrected Variable-b (Benchmark))
    //        • One series per residual (φ, η), black open-circle markers (24)
    //        • Error propagation: r = A/B,  dr = r * sqrt( (dA/A)^2 + (dB/B)^2 )
    // Saves directly to outDir as:
    //   - Headline_Ratio_Sigma_DeltaPhi_vsE.png
    //   - Headline_Ratio_Sigma_DeltaEta_vsE.png
    // =====================================================================
    auto makeRatioSeries = [&](const SigSeries& A, const SigSeries& B)->SigSeries {
        // Match by E-center with a small tolerance (1e-6 GeV)
        const long long SCALE = 1000000;
        std::unordered_map<long long, std::pair<double,std::pair<double,double>>> mb; // key -> (E, (B, dB))
        for (std::size_t j=0; j<B.E.size(); ++j) {
            const long long key = llround(B.E[j]*SCALE);
            const double Bj  = B.S[j];
            const double dBj = (B.dS.empty()? 0.0 : B.dS[j]);
            mb[key] = {B.E[j], {Bj, dBj}};
        }

        SigSeries R;
        for (std::size_t i=0; i<A.E.size(); ++i) {
            const long long key = llround(A.E[i]*SCALE);
            auto it = mb.find(key);
            if (it == mb.end()) continue;

            const double E   = A.E[i];
            const double Ai  = A.S[i];
            const double dAi = (A.dS.empty()? 0.0 : A.dS[i]);

            const double Bj  = it->second.second.first;
            const double dBj = it->second.second.second;

            if (Ai<=0.0 || Bj<=0.0) continue;

            const double r   = Ai / Bj;
            const double dr  = r * std::sqrt( (dAi>0? (dAi/Ai)*(dAi/Ai) : 0.0) +
                                              (dBj>0? (dBj/Bj)*(dBj/Bj) : 0.0) );

            R.E .push_back(E);
            R.S .push_back(r);
            R.dS.push_back(dr);
        }
        return R;
    };

    auto drawHeadlineRatio = [&](const SigSeries& R,
                                 const char* yTitle,
                                 const char* pngName)
    {
        if (R.E.empty()) {
            std::cerr << "[Playground][WARN] Ratio series empty for " << pngName << "\n";
            return;
        }

        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? (R.E.back()+1.0) : eEdges.back().second;

        // y-range with padding around 1
        double yLo = +1e30, yHi = -1e30;
        for (std::size_t i=0;i<R.S.size();++i) {
            yLo = std::min(yLo, R.S[i] - (R.dS.empty()?0.0:R.dS[i]));
            yHi = std::max(yHi, R.S[i] + (R.dS.empty()?0.0:R.dS[i]));
        }
        // Guard rails; ratios are expected near unity
        const double pad = 0.12*(yHi - yLo + 1e-9);
        yLo = std::min(yLo, 0.85);
        yHi = std::max(yHi, 1.15);
        yLo -= pad; yHi += pad;

        TCanvas c("cHeadlineRatio","Headline ratio σ(E)",900,600);
        c.SetLeftMargin(0.15);
        c.SetRightMargin(0.06);
        c.SetTopMargin(0.10);
        c.SetBottomMargin(0.18);

        TH1F fr("fr","",1,xMinE,xMaxE);
        fr.SetTitle("");
        fr.GetXaxis()->SetTitle("E  [GeV]");
        fr.GetYaxis()->SetTitle(yTitle);

        fr.GetXaxis()->SetTitleSize(0.050);
        fr.GetXaxis()->SetLabelSize(0.038);
        fr.GetXaxis()->SetTitleOffset(1.05);

        fr.GetYaxis()->SetTitleSize(0.050);
        fr.GetYaxis()->SetLabelSize(0.038);

        fr.SetMinimum(yLo);
        fr.SetMaximum(yHi);
        fr.Draw("AXIS");

        // reference line at 1
        TLine l1(xMinE, 1.0, xMaxE, 1.0);
        l1.SetLineStyle(2);
        l1.SetLineColor(kGray+2);
        l1.Draw();

        std::vector<double> ex(R.E.size(), 0.0);
        TGraphErrors gR((int)R.E.size(),
                        const_cast<double*>(R.E.data()),
                        const_cast<double*>(R.S.data()),
                        ex.data(),
                        const_cast<double*>(R.dS.data()));
        gR.SetMarkerStyle(24);         // "corrected" style
        gR.SetMarkerSize(1.05);
        gR.SetMarkerColor(kBlack);     // black ratio points
        gR.SetLineColor(kBlack);
        gR.Draw("P SAME");

        // Small caption-like header
        TLatex h; h.SetNDC(); h.SetTextAlign(13); h.SetTextSize(0.040);
        h.DrawLatex(0.16,0.955,"RMS Ratio Constant-b (Production)/Variable-b (Benchmark)");

        const std::string outPng = std::string(outDir) + "/" + pngName;
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // Build ratios (corrected Constant-b (Production) / corrected Variable-b (Benchmark))
    SigSeries rPhi = makeRatioSeries(sPhi_CLUScor, sPhi_PDCcor);
    SigSeries rEta = makeRatioSeries(sEta_CLUScor, sEta_PDCcor);

    // Draw ratio PNGs (black open circles, corrected style)
    drawHeadlineRatio(rPhi,
                      "Ratio  #sigma_{RMS}(#Delta#phi)_{const-b} / #sigma_{RMS}(#Delta#phi)_{var-b}",
                      "Headline_Ratio_Sigma_DeltaPhi_vsE.png");

    drawHeadlineRatio(rEta,
                      "Ratio  #sigma_{RMS}(#Delta#eta)_{const-b} / #sigma_{RMS}(#Delta#eta)_{var-b}",
                      "Headline_Ratio_Sigma_DeltaEta_vsE.png");

    // =====================================================================
    //           helpers for overlays (first-bin, μ/σ two-series, four-series)
    // =====================================================================
    auto makeFirstBinVariantOverlay = [&](TH1F* hrefRaw, TH1F* hrefCorr,
                                          Color_t col, const char* deltaSym,
                                          const char* outPng,
                                          // legend box (defaults kept for backward compatibility)
                                          double legX1 = 0.80, double legY1 = 0.78,
                                          double legX2 = 0.93, double legY2 = 0.90,
                                          // optional override header (nullptr -> use generic title)
                                          const char* headerOverride = nullptr,
                                          double headerX = 0.94, double headerY = 0.92,
                                          // μ/σ/χ²/NDF text block coords (kept for compatibility)
                                          double statsX = 0.94,
                                          double statsYraw  = 0.18,
                                          double statsYcorr = 0.13)
    {
        if (!hrefRaw || !hrefCorr || hrefRaw->Integral()==0 || hrefCorr->Integral()==0) return;

        // silence unused-parameter warnings (we now standardize layout in drawResidualPanel)
        (void)legX1; (void)legY1; (void)legX2; (void)legY2;
        (void)headerX; (void)headerY; (void)statsX; (void)statsYraw; (void)statsYcorr;

        TCanvas c(Form("cVar_%s", outPng), "variant first slice overlay", 1000, 750);

        const std::string title = (headerOverride && headerOverride[0] != '\0')
            ? std::string(headerOverride)
            : std::string("Clusterizer with/without Position Correction");

        std::vector<ResidualSpec> specs = {
            { hrefRaw,  Form("%s without Correction", deltaSym), (Color_t)col, (Style_t)20, 1 },
            { hrefCorr, Form("%s with Correction",     deltaSym), (Color_t)col, (Style_t)24, 2 }
        };

        drawResidualPanel(&c, title,
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);

        c.SaveAs(outPng);
        c.Close();
    };


    // ---- μ/σ vs E overlay (two series: raw vs corr) ----
    auto makeMuSigmaVariantOverlay = [&](const std::vector<TH1F*>& Vraw,
                                         const std::vector<TH1F*>& Vcor,
                                         Color_t col, Style_t mkRaw, Style_t mkCor,
                                         const char* ymuTitle,
                                         const char* outPng)
    {
        std::vector<double> eCtr, muR, dmuR, sgR, dsgR, muC, dmuC, sgC, dsgC;
        for (std::size_t i=0;i<eEdges.size();++i){
            if (!Vraw[i] || !Vcor[i] || Vraw[i]->Integral()==0 || Vcor[i]->Integral()==0) continue;
            eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second));

            double mR, meR, rR, reR, mC, meC, rC, reC;
            std::tie(mR, meR, rR, reR) = statsFromHistRange(Vraw[i], xMin, xMax);
            std::tie(mC, meC, rC, reC) = statsFromHistRange(Vcor[i], xMin, xMax);

            muR.push_back(mR); dmuR.push_back(meR);
            sgR.push_back(rR); dsgR.push_back(reR);
            muC.push_back(mC); dmuC.push_back(meC);
            sgC.push_back(rC); dsgC.push_back(reC);
        }
        if (eCtr.empty()) return;

        const double xAxisMin = 0.0;
        const double xAxisMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
        std::vector<double> ex(eCtr.size(), 0.0);

        TCanvas c("cVar_MuSig", "variant mu/sigma", 900, 800);
        auto pads = makeEqualHeightTwinPads(c, "VarMuSig", 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first; TPad* pBot = pads.second;

        // μ(E)
        pTop->cd();
        double muLo=1e30, muHi=-1e30;
        for (std::size_t i=0;i<muR.size();++i){
            muLo = std::min(muLo, std::min(muR[i]-dmuR[i], muC[i]-dmuC[i]));
            muHi = std::max(muHi, std::max(muR[i]+dmuR[i], muC[i]+dmuC[i]));
        }
        const double padMu = 0.25*(muHi-muLo);
        TH1F frU("frU",Form(";E  [GeV];%s  [rad]", ymuTitle),1,xAxisMin,xAxisMax);
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

        TGraphErrors gMuR(eCtr.size(), eCtr.data(), muR.data(), ex.data(), dmuR.data());
        gMuR.SetMarkerStyle(mkRaw); gMuR.SetMarkerColor(col); gMuR.SetLineColor(col); gMuR.SetMarkerSize(1.0);
        gMuR.Draw("P SAME");
        TGraphErrors gMuC(eCtr.size(), eCtr.data(), muC.data(), ex.data(), dmuC.data());
        gMuC.SetMarkerStyle(mkCor); gMuC.SetMarkerColor(col); gMuC.SetLineColor(col); gMuC.SetMarkerSize(1.0);
        gMuC.Draw("P SAME");

        TLegend legU(0.72,0.75,0.93,0.92);
        legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.035);
        legU.AddEntry(&gMuR,"no corr, cluster","p");
        legU.AddEntry(&gMuC,"CorrectPosition, cluster","p");
        legU.Draw();

        // σ(E)
        pBot->cd();
        double sgHi=-1e30; for (double s : sgR) sgHi = std::max(sgHi,s); for (double s : sgC) sgHi = std::max(sgHi,s);
        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

        TGraphErrors gSgR(eCtr.size(), eCtr.data(), sgR.data(), ex.data(), dsgR.data());
        gSgR.SetMarkerStyle(mkRaw); gSgR.SetMarkerColor(col); gSgR.SetLineColor(col); gSgR.SetMarkerSize(1.0);
        gSgR.Draw("P SAME");
        TGraphErrors gSgC(eCtr.size(), eCtr.data(), sgC.data(), ex.data(), dsgC.data());
        gSgC.SetMarkerStyle(mkCor); gSgC.SetMarkerColor(col); gSgC.SetLineColor(col); gSgC.SetMarkerSize(1.0);
        gSgC.Draw("P SAME");

        // Header
        c.cd();
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.040);
        h.DrawLatex(0.50,0.985,"no corr, cluster  vs  CorrectPosition, cluster  -  Mean / RMS vs E");


        c.SaveAs(outPng);
        c.Close();

        // CSV: two series, residual auto-detected from input hist names
        const std::string resid = detectResidual(Vraw);
        saveMuSigmaCSV(outPng,
                       { "no corr, cluster", "CorrectPosition, cluster" },
                       { resid, resid },
                       { eCtr, eCtr },
                       { muR,  muC  }, { dmuR, dmuC },
                       { sgR,  sgC  }, { dsgR, dsgC });
    };

    // ---- labeled μ/σ vs E overlay (two series) for custom legends/headers ----
    auto makeMuSigmaVariantOverlayWithLabels = [&](const std::vector<TH1F*>& Vraw,
                                                   const std::vector<TH1F*>& Vcor,
                                                   Color_t col, Style_t mkRaw, Style_t mkCor,
                                                   const char* ymuTitle,
                                                   const char* outPng,
                                                   const char* labelRaw,
                                                   const char* labelCorr,
                                                   const char* headerText)
    {
        std::vector<double> eCtr, muR, dmuR, sgR, dsgR, muC, dmuC, sgC, dsgC;
        for (std::size_t i=0;i<eEdges.size();++i){
            if (!Vraw[i] || !Vcor[i] || Vraw[i]->Integral()==0 || Vcor[i]->Integral()==0) continue;
            eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second));

            double mR, meR, rR, reR, mC, meC, rC, reC;
            std::tie(mR, meR, rR, reR) = statsFromHistRange(Vraw[i], xMin, xMax);
            std::tie(mC, meC, rC, reC) = statsFromHistRange(Vcor[i], xMin, xMax);

            muR.push_back(mR); dmuR.push_back(meR);
            sgR.push_back(rR); dsgR.push_back(reR);
            muC.push_back(mC); dmuC.push_back(meC);
            sgC.push_back(rC); dsgC.push_back(reC);

        }
        if (eCtr.empty()) return;

        const double xAxisMin = 0.0;
        const double xAxisMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
        std::vector<double> ex(eCtr.size(), 0.0);

        TCanvas c("cVar_MuSig_lbl", "variant mu/sigma (labeled)", 900, 800);
        auto pads = makeEqualHeightTwinPads(c, "VarMuSig_lbl", 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first; TPad* pBot = pads.second;

        // μ(E)
        pTop->cd();
        double muLo=1e30, muHi=-1e30;
        for (std::size_t i=0;i<muR.size();++i){
            muLo = std::min(muLo, std::min(muR[i]-dmuR[i], muC[i]-dmuC[i]));
            muHi = std::max(muHi, std::max(muR[i]+dmuR[i], muC[i]+dmuC[i]));
        }
        const double padMu = 0.25*(muHi-muLo);
        TH1F frU("frU",Form(";E  [GeV];%s  [rad]", ymuTitle),1,xAxisMin,xAxisMax);
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

        TGraphErrors gMuR(eCtr.size(), eCtr.data(), muR.data(), ex.data(), dmuR.data());
        gMuR.SetMarkerStyle(mkRaw); gMuR.SetMarkerColor(col); gMuR.SetLineColor(col); gMuR.SetMarkerSize(1.0);
        gMuR.Draw("P SAME");
        TGraphErrors gMuC(eCtr.size(), eCtr.data(), muC.data(), ex.data(), dmuC.data());
        gMuC.SetMarkerStyle(mkCor); gMuC.SetMarkerColor(col); gMuC.SetLineColor(col); gMuC.SetMarkerSize(1.0);
        gMuC.Draw("P SAME");

        TLegend legU(0.68,0.75,0.93,0.92);
        legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.035);
        legU.AddEntry(&gMuR,labelRaw,"p");
        legU.AddEntry(&gMuC,labelCorr,"p");
        legU.Draw();

        // σ(E)
        pBot->cd();
        double sgHi=-1e30; for (double s : sgR) sgHi = std::max(sgHi,s); for (double s : sgC) sgHi = std::max(sgHi,s);
        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

        TGraphErrors gSgR(eCtr.size(), eCtr.data(), sgR.data(), ex.data(), dsgR.data());
        gSgR.SetMarkerStyle(mkRaw); gSgR.SetMarkerColor(col); gSgR.SetLineColor(col); gSgR.SetMarkerSize(1.0);
        gSgR.Draw("P SAME");
        TGraphErrors gSgC(eCtr.size(), eCtr.data(), sgC.data(), ex.data(), dsgC.data());
        gSgC.SetMarkerStyle(mkCor); gSgC.SetMarkerColor(col); gSgC.SetLineColor(col); gSgC.SetMarkerSize(1.0);
        gSgC.Draw("P SAME");

        // Header
        c.cd();
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.032);
        h.DrawLatex(0.50,0.985, headerText);

        c.SaveAs(outPng);
        c.Close();

        // CSV: honor custom labels
        const std::string resid = detectResidual(Vraw);
        saveMuSigmaCSV(outPng,
                       { std::string(labelRaw), std::string(labelCorr) },
                       { resid, resid },
                       { eCtr, eCtr },
                       { muR,  muC  }, { dmuR, dmuC },
                       { sgR,  sgC  }, { dsgR, dsgC });
    };


    // ---- μ/σ vs E overlay (four series: φ/η raw/corr) ----
    auto makeMuSigmaFourSeriesOverlay = [&](const std::vector<TH1F*>& Praw,
                                            const std::vector<TH1F*>& Pcorr,
                                            const std::vector<TH1F*>& Eraw,
                                            const std::vector<TH1F*>& Ecorr,
                                            const char* outPng)
    {
        std::vector<double> eCtrP, muPR, dmuPR, sgPR, dsgPR, muPC, dmuPC, sgPC, dsgPC;
        std::vector<double> eCtrE, muER, dmuER, sgER, dsgER, muEC, dmuEC, sgEC, dsgEC;

        for (std::size_t i=0;i<eEdges.size();++i){
            if (Praw[i] && Praw[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Praw[i], xMin, xMax);
                eCtrP.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muPR.push_back(m);  dmuPR.push_back(me);
                sgPR.push_back(r);  dsgPR.push_back(re);
            }
            if (Pcorr[i] && Pcorr[i]->Integral()>0 && Praw[i] && Praw[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Pcorr[i], xMin, xMax);
                muPC.push_back(m);  dmuPC.push_back(me);
                sgPC.push_back(r);  dsgPC.push_back(re);
            }

            if (Eraw[i] && Eraw[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Eraw[i], xMin, xMax);
                eCtrE.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muER.push_back(m);  dmuER.push_back(me);
                sgER.push_back(r);  dsgER.push_back(re);
            }
            if (Ecorr[i] && Ecorr[i]->Integral()>0 && Eraw[i] && Eraw[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Ecorr[i], xMin, xMax);
                muEC.push_back(m);  dmuEC.push_back(me);
                sgEC.push_back(r);  dsgEC.push_back(re);
            }
        }

        if (eCtrP.empty() && eCtrE.empty()) return;

        const double xAxisMin = 0.0;
        const double xAxisMax = (eEdges.back().second + eEdges.back().first)*0.5 + 0.5*(eEdges.back().second - eEdges.back().first);
        TCanvas c("cVar_MuSig4", "four series mu/sigma", 1000, 850);

        auto pads = makeEqualHeightTwinPads(c, "MuSig4", 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first; TPad* pBot = pads.second;

        // μ(E)
        pTop->cd();
        double muLo=1e30, muHi=-1e30;
        auto upd = [&](const std::vector<double>& m, const std::vector<double>& dm){
            for (std::size_t i=0;i<m.size();++i){ muLo = std::min(muLo, m[i]-dm[i]); muHi = std::max(muHi, m[i]+dm[i]); }
        };
        upd(muPR,dmuPR); upd(muPC,dmuPC); upd(muER,dmuER); upd(muEC,dmuEC);
        const double padMu = 0.25*(muHi-muLo);
        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
        frU.SetMinimum(muLo - padMu); frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0); frU.GetXaxis()->SetTitleSize(0.0); frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

        std::vector<double> exP(eCtrP.size(),0.0), exE(eCtrE.size(),0.0);
        TGraphErrors gMuPR(eCtrP.size(), eCtrP.data(), muPR.data(), exP.data(), dmuPR.data());
        gMuPR.SetMarkerStyle(20); gMuPR.SetMarkerColor(kRed+1);  gMuPR.SetLineColor(kRed+1);  gMuPR.Draw("P SAME");
        TGraphErrors gMuPC(eCtrP.size(), eCtrP.data(), muPC.data(), exP.data(), dmuPC.data());
        gMuPC.SetMarkerStyle(24); gMuPC.SetMarkerColor(kRed+1);  gMuPC.SetLineColor(kRed+1);  gMuPC.Draw("P SAME");
        TGraphErrors gMuER(eCtrE.size(), eCtrE.data(), muER.data(), exE.data(), dmuER.data());
        gMuER.SetMarkerStyle(20); gMuER.SetMarkerColor(kBlue+1); gMuER.SetLineColor(kBlue+1); gMuER.Draw("P SAME");
        TGraphErrors gMuEC(eCtrE.size(), eCtrE.data(), muEC.data(), exE.data(), dmuEC.data());
        gMuEC.SetMarkerStyle(24); gMuEC.SetMarkerColor(kBlue+1); gMuEC.SetLineColor(kBlue+1); gMuEC.Draw("P SAME");

        TLegend legU(0.65, 0.72, 0.93, 0.85);  // wider box for 2 columns
        legU.SetBorderSize(0);
        legU.SetFillStyle(0);
        legU.SetTextSize(0.04);
        legU.SetNColumns(2);  // <--- two-column layout

        legU.AddEntry(&gMuPR,"#Delta#phi no corr","p");
        legU.AddEntry(&gMuER,"#Delta#eta no corr","p");
        legU.AddEntry(&gMuPC,"#Delta#phi with corr","p");
        legU.AddEntry(&gMuEC,"#Delta#eta with corr","p");

        legU.Draw();

        // σ(E)
        pBot->cd();
        double sgHi=-1e30;
        auto upds = [&](const std::vector<double>& s){ for(double v: s) sgHi = std::max(sgHi, v); };
        upds(sgPR); upds(sgPC); upds(sgER); upds(sgEC);
        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

        TGraphErrors gSgPR(eCtrP.size(), eCtrP.data(), sgPR.data(), exP.data(), dsgPR.data());
        gSgPR.SetMarkerStyle(20); gSgPR.SetMarkerColor(kRed+1);  gSgPR.SetLineColor(kRed+1);  gSgPR.Draw("P SAME");
        TGraphErrors gSgPC(eCtrP.size(), eCtrP.data(), sgPC.data(), exP.data(), dsgPC.data());
        gSgPC.SetMarkerStyle(24); gSgPC.SetMarkerColor(kRed+1);  gSgPC.SetLineColor(kRed+1);  gSgPC.Draw("P SAME");
        TGraphErrors gSgER(eCtrE.size(), eCtrE.data(), sgER.data(), exE.data(), dsgER.data());
        gSgER.SetMarkerStyle(20); gSgER.SetMarkerColor(kBlue+1); gSgER.SetLineColor(kBlue+1); gSgER.Draw("P SAME");
        TGraphErrors gSgEC(eCtrE.size(), eCtrE.data(), sgEC.data(), exE.data(), dsgEC.data());
        gSgEC.SetMarkerStyle(24); gSgEC.SetMarkerColor(kBlue+1); gSgEC.SetLineColor(kBlue+1); gSgEC.Draw("P SAME");

        c.cd();
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.032);
        h.DrawLatex(0.50,0.985,"#Delta#phi / #Delta#eta  (not corr vs corr)  -  Mean / RMS vs E");


        c.SaveAs(outPng);
        c.Close();

        // CSV: 4 series (phi raw/corr, eta raw/corr) - cluster variants
        saveMuSigmaCSV(outPng,
                       { "no corr, cluster", "CorrectPosition, cluster",
                         "no corr, cluster", "CorrectPosition, cluster" },
                       { "phi", "phi", "eta", "eta" },
                       { eCtrP, eCtrP, eCtrE, eCtrE },
                       { muPR,  muPC,  muER,  muEC  }, { dmuPR, dmuPC, dmuER, dmuEC },
                       { sgPR,  sgPC,  sgER,  sgEC  }, { dsgPR, dsgPC, dsgER, dsgEC });
    };

    // =====================================================================
    //                               PNGs
    // =====================================================================

    // --------------------- Clusterizer (existing) ---------------------
    if (!Hphi.empty() && !HphiCorr.empty())
        makeFirstBinVariantOverlay(
            Hphi[0], HphiCorr[0], kRed+1, "#Delta#phi",
            (dirRawVsCorr + "/DeltaPhi_cpRaw_vs_cpCorr_FirstBin_Overlay.png").c_str(),
            // legend box (shifted left)
            0.62, 0.78, 0.85, 0.90,
            // centered header override at (0.50, 0.88)
            Form("Clusterizer with/without Position Correction, [%.0f, %.0f) GeV",
                 eEdges.front().first, eEdges.front().second),
            0.50, 0.88,
            // μ/σ/χ² text block near top-right (moved up)
            0.94, 0.84, 0.79
        );

    if (!Heta.empty() && !HetaCorr.empty())
        makeFirstBinVariantOverlay(
            Heta[0], HetaCorr[0], kBlue+1, "#Delta#eta",
            (dirRawVsCorr + "/DeltaEta_cpRaw_vs_cpCorr_FirstBin_Overlay.png").c_str(),
            // legend box (shifted left)
            0.62, 0.78, 0.85, 0.90,
            // centered header override at (0.50, 0.88)
            Form("Clusterizer with/without Position Correction, [%.0f, %.0f) GeV",
                 eEdges.front().first, eEdges.front().second),
            0.50, 0.88,
            // μ/σ/χ² text block near top-right (moved up)
            0.94, 0.84, 0.79
        );

    // μ/σ vs E overlays per residual (two series)
    if (!Hphi.empty() && !HphiCorr.empty())
        makeMuSigmaVariantOverlay(Hphi, HphiCorr, kRed+1, 20, 24, "#mu",
            (dirRawVsCorr + "/DeltaPhi_cpRaw_vs_cpCorr_MeanSigmaVsE_Overlay.png").c_str());

    if (!Heta.empty() && !HetaCorr.empty())
        makeMuSigmaVariantOverlay(Heta, HetaCorr, kBlue+1, 20, 24, "#mu",
            (dirRawVsCorr + "/DeltaEta_cpRaw_vs_cpCorr_MeanSigmaVsE_Overlay.png").c_str());

    // μ/σ vs E overlay for FOUR series (Δφ/Δη raw+corr)
    if (!Hphi.empty() && !HphiCorr.empty() && !Heta.empty() && !HetaCorr.empty())
        makeMuSigmaFourSeriesOverlay(Hphi, HphiCorr, Heta, HetaCorr,
            (dirRawVsCorr + "/DeltaEtaPhi_cpRaw_vs_cpCorr_MeanSigmaVsE_Overlay_4Series.png").c_str());

    // --------------------- From-scratch (NEW) ------------------------
    // First-bin overlays (COUNTS) for PDC raw vs PDC corr
    if (!HphiScratchRaw.empty() && !HphiScratchCorr.empty()
        && HphiScratchRaw[0] && HphiScratchCorr[0]
        && HphiScratchRaw[0]->Integral()>0 && HphiScratchCorr[0]->Integral()>0)
        makeFirstBinVariantOverlay(HphiScratchRaw[0], HphiScratchCorr[0], kRed+1,  "#Delta#phi",
            (dirScratchRawVsCorr + "/DeltaPhi_raw_vs_corr_FirstBin_Overlay.png").c_str());

    if (!HetaScratchRaw.empty() && !HetaScratchCorr.empty()
        && HetaScratchRaw[0] && HetaScratchCorr[0]
        && HetaScratchRaw[0]->Integral()>0 && HetaScratchCorr[0]->Integral()>0)
        makeFirstBinVariantOverlay(HetaScratchRaw[0], HetaScratchCorr[0], kBlue+1, "#Delta#eta",
            (dirScratchRawVsCorr + "/DeltaEta_raw_vs_corr_FirstBin_Overlay.png").c_str());

    // μ/σ vs E overlays per residual (two series), with custom labels:
    //   "no corr, scratch"   vs   "b(E) corr, scratch"
    if (!HphiScratchRaw.empty() && !HphiScratchCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HphiScratchRaw, HphiScratchCorr,
            kRed+1, 20, 24, "#mu",
            (dirScratchRawVsCorr + "/DeltaPhi_raw_vs_corr_MeanSigmaVsE_Overlay.png").c_str(),
            "no corr, scratch", "b(E) corr, scratch",
            "no corr, scratch  vs  b(E) corr, scratch  -  Mean / RMS vs E");

    if (!HetaScratchRaw.empty() && !HetaScratchCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HetaScratchRaw, HetaScratchCorr,
            kBlue+1, 20, 24, "#mu",
            (dirScratchRawVsCorr + "/DeltaEta_raw_vs_corr_MeanSigmaVsE_Overlay.png").c_str(),
            "no corr, scratch", "b(E) corr, scratch",
            "no corr, scratch  vs  b(E) corr, scratch  -  Mean / RMS vs E");

    // μ/σ vs E overlay for FOUR series (Δφ/Δη raw+corr, from-scratch)
    if (!HphiScratchRaw.empty() && !HphiScratchCorr.empty()
        && !HetaScratchRaw.empty() && !HetaScratchCorr.empty())
        makeMuSigmaFourSeriesOverlay(HphiScratchRaw, HphiScratchCorr,
                                     HetaScratchRaw, HetaScratchCorr,
            (dirScratchRawVsCorr + "/DeltaEtaPhi_raw_vs_corr_MeanSigmaVsE_Overlay_4Series.png").c_str());

    // Helper: four‑series μ/σ(E) overlay with custom labels for “from scratch vs clusterizer”
    auto makeMuSigmaFourSeriesOverlayFSvClus = [&](const std::vector<TH1F*>& Pfs,
                                                   const std::vector<TH1F*>& Pcl,
                                                   const std::vector<TH1F*>& Efs,
                                                   const std::vector<TH1F*>& Ecl,
                                                   const char* outPng,
                                                   const char* headerText)
    {
        std::vector<double> eCtrPfs, muPfs, dmuPfs, sgPfs, dsgPfs;
        std::vector<double> eCtrPcl, muPcl, dmuPcl, sgPcl, dsgPcl;
        std::vector<double> eCtrEfs, muEfs, dmuEfs, sgEfs, dsgEfs;
        std::vector<double> eCtrEcl, muEcl, dmuEcl, sgEcl, dsgEcl;

        for (std::size_t i=0;i<eEdges.size();++i){
            if (Pfs[i] && Pfs[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Pfs[i], xMin, xMax);
                eCtrPfs.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muPfs.push_back(m);  dmuPfs.push_back(me);
                sgPfs.push_back(r);  dsgPfs.push_back(re);
            }
            if (Pcl[i] && Pcl[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Pcl[i], xMin, xMax);
                eCtrPcl.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muPcl.push_back(m);  dmuPcl.push_back(me);
                sgPcl.push_back(r);  dsgPcl.push_back(re);
            }
            if (Efs[i] && Efs[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Efs[i], xMin, xMax);
                eCtrEfs.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muEfs.push_back(m);  dmuEfs.push_back(me);
                sgEfs.push_back(r);  dsgEfs.push_back(re);
            }
            if (Ecl[i] && Ecl[i]->Integral()>0){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(Ecl[i], xMin, xMax);
                eCtrEcl.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                muEcl.push_back(m);  dmuEcl.push_back(me);
                sgEcl.push_back(r);  dsgEcl.push_back(re);
            }

        }
        if (eCtrPfs.empty() && eCtrPcl.empty() && eCtrEfs.empty() && eCtrEcl.empty()) return;

        const double xAxisMin = 0.0;
        const double xAxisMax = eEdges.back().second - 0.5*(eEdges.back().second - eEdges.back().first)
                              + 0.5*(eEdges.back().second - eEdges.back().first);

        std::vector<double> exPfs(eCtrPfs.size(),0.0), exPcl(eCtrPcl.size(),0.0);
        std::vector<double> exEfs(eCtrEfs.size(),0.0), exEcl(eCtrEcl.size(),0.0);

        TCanvas c("cFSvClus_4Series","FS vs Clusterizer - four series μ/σ(E)",1000,850);
        auto pads = makeEqualHeightTwinPads(c, "FSvClus4", 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first; TPad* pBot = pads.second;
        // ---- μ(E)
        pTop->cd();
        double muLo=1e30, muHi=-1e30;
        auto updRange = [&](const std::vector<double>& m, const std::vector<double>& dm){
            for (std::size_t i=0;i<m.size();++i){ muLo = std::min(muLo, m[i]-dm[i]); muHi = std::max(muHi, m[i]+dm[i]); }
        };
        updRange(muPfs,dmuPfs); updRange(muPcl,dmuPcl); updRange(muEfs,dmuEfs); updRange(muEcl,dmuEcl);
        const double padMu = 0.25*(muHi-muLo);
        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

        TGraphErrors gMuPfs(eCtrPfs.size(), eCtrPfs.data(), muPfs.data(), exPfs.data(), dmuPfs.data());
        gMuPfs.SetMarkerStyle(20); gMuPfs.SetMarkerColor(kRed+1);  gMuPfs.SetLineColor(kRed+1);  gMuPfs.Draw("P SAME");
        TGraphErrors gMuEfs(eCtrEfs.size(), eCtrEfs.data(), muEfs.data(), exEfs.data(), dmuEfs.data());
        gMuEfs.SetMarkerStyle(20); gMuEfs.SetMarkerColor(kBlue+1); gMuEfs.SetLineColor(kBlue+1); gMuEfs.Draw("P SAME");
        TGraphErrors gMuPcl(eCtrPcl.size(), eCtrPcl.data(), muPcl.data(), exPcl.data(), dmuPcl.data());
        gMuPcl.SetMarkerStyle(24); gMuPcl.SetMarkerColor(kRed+1);  gMuPcl.SetLineColor(kRed+1);  gMuPcl.Draw("P SAME");
        TGraphErrors gMuEcl(eCtrEcl.size(), eCtrEcl.data(), muEcl.data(), exEcl.data(), dmuEcl.data());
        gMuEcl.SetMarkerStyle(24); gMuEcl.SetMarkerColor(kBlue+1); gMuEcl.SetLineColor(kBlue+1); gMuEcl.Draw("P SAME");

        TLegend legU(0.62, 0.77, 0.9, 0.86);
        legU.SetBorderSize(0);
        legU.SetFillStyle(0);
        legU.SetTextSize(0.033);
        legU.SetNColumns(2);  // 2-column legend as requested
        legU.AddEntry(&gMuPfs, "#Delta#phi from scratch", "p");
        legU.AddEntry(&gMuEfs, "#Delta#eta from scratch", "p");
        legU.AddEntry(&gMuPcl, "#Delta#phi clusterizer",  "p");
        legU.AddEntry(&gMuEcl, "#Delta#eta clusterizer",  "p");
        legU.Draw();

        // ---- σ(E)
        pBot->cd();
        double sgHi=-1e30;
        auto updSg = [&](const std::vector<double>& s){ for(double v: s) sgHi = std::max(sgHi, v); };
        updSg(sgPfs); updSg(sgPcl); updSg(sgEfs); updSg(sgEcl);
        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

        TGraphErrors gSgPfs(eCtrPfs.size(), eCtrPfs.data(), sgPfs.data(), exPfs.data(), dsgPfs.data());
        gSgPfs.SetMarkerStyle(20); gSgPfs.SetMarkerColor(kRed+1);  gSgPfs.SetLineColor(kRed+1);  gSgPfs.Draw("P SAME");
        TGraphErrors gSgEfs(eCtrEfs.size(), eCtrEfs.data(), sgEfs.data(), exEfs.data(), dsgEfs.data());
        gSgEfs.SetMarkerStyle(20); gSgEfs.SetMarkerColor(kBlue+1); gSgEfs.SetLineColor(kBlue+1); gSgEfs.Draw("P SAME");
        TGraphErrors gSgPcl(eCtrPcl.size(), eCtrPcl.data(), sgPcl.data(), exPcl.data(), dsgPcl.data());
        gSgPcl.SetMarkerStyle(24); gSgPcl.SetMarkerColor(kRed+1);  gSgPcl.SetLineColor(kRed+1);  gSgPcl.Draw("P SAME");
        TGraphErrors gSgEcl(eCtrEcl.size(), eCtrEcl.data(), sgEcl.data(), exEcl.data(), dsgEcl.data());
        gSgEcl.SetMarkerStyle(24); gSgEcl.SetMarkerColor(kBlue+1); gSgEcl.SetLineColor(kBlue+1); gSgEcl.Draw("P SAME");

        c.cd();
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.032);
        h.DrawLatex(0.50,0.985, headerText);

        c.SaveAs(outPng);
        c.Close();

        // CSV: infer variants from histogram names (scratch vs clusterizer; raw or corr)
        const std::string vPfs = variantFromHistName(firstHistName(Pfs));
        const std::string vPcl = variantFromHistName(firstHistName(Pcl));
        const std::string vEfs = variantFromHistName(firstHistName(Efs));
        const std::string vEcl = variantFromHistName(firstHistName(Ecl));

        saveMuSigmaCSV(outPng,
                       { vPfs, vPcl, vEfs, vEcl },
                       { "phi", "phi", "eta", "eta" },
                       { eCtrPfs, eCtrPcl, eCtrEfs, eCtrEcl },
                       { muPfs,   muPcl,   muEfs,   muEcl   }, { dmuPfs, dmuPcl, dmuEfs, dmuEcl },
                       { sgPfs,   sgPcl,   sgEfs,   sgEcl   }, { dsgPfs, dsgPcl, dsgEfs, dsgEcl });

    };


    // ---------------------------------------------------------------------
    // From‑scratch vs Clusterizer  (NO position correction)
    // ---------------------------------------------------------------------

    // Δφ - first‑bin overlay (counts)
    if (!HphiScratchRaw.empty() && !Hphi.empty()
        && HphiScratchRaw[0] && Hphi[0]
        && HphiScratchRaw[0]->Integral()>0 && Hphi[0]->Integral()>0)
    {
        TCanvas cFSvCl_Phi_First("cFSvCl_Phi_First_NoCorr",
                                 "#Delta #phi from-scratch vs clusterizer (no corr) - first bin",
                                 1000, 750);
        std::vector<ResidualSpec> specs = {
            { HphiScratchRaw[0], "#Delta#phi  from scratch", (Color_t)(kRed+1),  (Style_t)20, 1 },
            { Hphi[0],           "#Delta#phi  clusterizer",  (Color_t)(kRed+1),  (Style_t)24, 2 }
        };
        drawResidualPanel(&cFSvCl_Phi_First,
                          "From-scratch vs Clusterizer, No Position Correction",
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);
        cFSvCl_Phi_First.SaveAs((dirFSvsClusNoCorr + "/DeltaPhi_FSvsClus_FirstBin_Overlay.png").c_str());
        cFSvCl_Phi_First.Close();
    }

    // Δη - first‑bin overlay (counts)
    if (!HetaScratchRaw.empty() && !Heta.empty()
        && HetaScratchRaw[0] && Heta[0]
        && HetaScratchRaw[0]->Integral()>0 && Heta[0]->Integral()>0)
    {
        TCanvas cFSvCl_Eta_First("cFSvCl_Eta_First_NoCorr",
                                 "#Delta #eta from-scratch vs clusterizer (no corr) - first bin",
                                 1000, 750);
        std::vector<ResidualSpec> specs = {
            { HetaScratchRaw[0], "#Delta#eta  from scratch", (Color_t)(kBlue+1), (Style_t)20, 1 },
            { Heta[0],           "#Delta#eta  clusterizer",  (Color_t)(kBlue+1), (Style_t)24, 2 }
        };
        drawResidualPanel(&cFSvCl_Eta_First,
                          "From-scratch vs Clusterizer, No Position Correction",
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);
        cFSvCl_Eta_First.SaveAs((dirFSvsClusNoCorr + "/DeltaEta_FSvsClus_FirstBin_Overlay.png").c_str());
        cFSvCl_Eta_First.Close();
    }

    // Δφ - 4×2 table overlay (all energy bins)
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Phi_Table_NoCorr",
                   "Δφ from-scratch vs clusterizer (no corr) - table", 1600, 900);
        cT.SetTopMargin(0.10);
        cT.Divide(nCol, nRow, 0, 0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE) {
            if (!HphiScratchRaw[iE] || !Hphi[iE] ||
                HphiScratchRaw[iE]->Integral()==0 || Hphi[iE]->Integral()==0) { ++pad; continue; }

            cT.cd(pad++);
            setPadMargins();

            std::vector<ResidualSpec> specs = {
                { HphiScratchRaw[iE], "#Delta#phi  from scratch", (Color_t)(kRed+1),  (Style_t)20, 1 },
                { Hphi[iE],           "#Delta#phi  clusterizer",  (Color_t)(kRed+1),  (Style_t)24, 2 }
            };

            drawResidualPanel(gPad,
                              "From-scratch vs Clusterizer, No Position Correction",
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 0.85);
        }
        cT.cd(0);
        drawHeaderNDC("#Delta#phi  -  From-scratch vs Clusterizer  (No Position Correction)", 0.975, 0.045);
        cT.SaveAs((dirFSvsClusNoCorr + "/DeltaPhi_FSvsClus_AllBins_Table.png").c_str());
        cT.Close();
    }

    // Δη - 4×2 table overlay (all energy bins)
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Eta_Table_NoCorr",
                   "#Delta #eta from-scratch vs clusterizer (no corr) - table", 1600, 900);
        cT.SetTopMargin(0.10);
        cT.Divide(nCol, nRow, 0, 0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE) {
            if (!HetaScratchRaw[iE] || !Heta[iE] ||
                HetaScratchRaw[iE]->Integral()==0 || Heta[iE]->Integral()==0) { ++pad; continue; }

            cT.cd(pad++);
            setPadMargins();

            std::vector<ResidualSpec> specs = {
                { HetaScratchRaw[iE], "#Delta#eta  from scratch", (Color_t)(kBlue+1), (Style_t)20, 1 },
                { Heta[iE],           "#Delta#eta  clusterizer",  (Color_t)(kBlue+1), (Style_t)24, 2 }
            };

            drawResidualPanel(gPad,
                              "From-scratch vs Clusterizer, No Position Correction",
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 0.85);
        }
        cT.cd(0);
        drawHeaderNDC("#Delta#eta  -  From-scratch vs Clusterizer  (No Position Correction)", 0.975, 0.045);
        cT.SaveAs((dirFSvsClusNoCorr + "/DeltaEta_FSvsClus_AllBins_Table.png").c_str());
        cT.Close();
    }

    // Δφ - μ/σ(E) overlay (two series): from‑scratch vs clusterizer, no corr
    if (!HphiScratchRaw.empty() && !Hphi.empty())
        makeMuSigmaVariantOverlayWithLabels(HphiScratchRaw, Hphi,
            kRed+1, 20, 24, "#mu",
            (dirFSvsClusNoCorr + "/DeltaPhi_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch", "clusterizer",
            "#Delta#phi  -  from scratch vs clusterizer  (No Position Correction)  -  Mean / RMS vs E");

    // Δη - μ/σ(E) overlay (two series): from‑scratch vs clusterizer, no corr
    if (!HetaScratchRaw.empty() && !Heta.empty())
        makeMuSigmaVariantOverlayWithLabels(HetaScratchRaw, Heta,
            kBlue+1, 20, 24, "#mu",
            (dirFSvsClusNoCorr + "/DeltaEta_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch", "clusterizer",
            "#Delta#eta  -  from scratch vs clusterizer  (No Position Correction)  -  Mean / RMS vs E");

    // NEW: 4-series μ/σ(E) overlay (from‑scratch vs clusterizer), no corr
    if (!HphiScratchRaw.empty() && !Hphi.empty()
        && !HetaScratchRaw.empty() && !Heta.empty())
        makeMuSigmaFourSeriesOverlayFSvClus(HphiScratchRaw, Hphi,
                                            HetaScratchRaw, Heta,
            (dirFSvsClusNoCorr + "/DeltaEtaPhi_FSvsClus_MeanSigmaVsE_Overlay_4Series.png").c_str(),
            "#Delta#phi / #Delta#eta  -  from scratch vs clusterizer  (No Position Correction)  -  Mean / RMS vs E");


    // Δφ - first‑bin overlay (counts), corrected
    if (!HphiScratchCorr.empty() && !HphiCorr.empty()
        && HphiScratchCorr[0] && HphiCorr[0]
        && HphiScratchCorr[0]->Integral()>0 && HphiCorr[0]->Integral()>0)
    {
        TCanvas cFSvCl_Phi_FirstC("cFSvCl_Phi_First_WithCorr",
                                  "#Delta #phi from-scratch vs clusterizer (with corr) - first bin",
                                  1000, 750);
        std::vector<ResidualSpec> specs = {
            { HphiScratchCorr[0], "#Delta#phi  from scratch (corr)", (Color_t)(kRed+1),  (Style_t)20, 1 },
            { HphiCorr[0],        "#Delta#phi  clusterizer (corr)",  (Color_t)(kRed+1),  (Style_t)24, 2 }
        };
        drawResidualPanel(&cFSvCl_Phi_FirstC,
                          "From-scratch vs Clusterizer, With Position Correction",
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);
        cFSvCl_Phi_FirstC.SaveAs((dirFSvsClusWithCorr + "/DeltaPhi_FSvsClus_FirstBin_Overlay.png").c_str());
        cFSvCl_Phi_FirstC.Close();
    }

    // Δη - first‑bin overlay (counts), corrected
    if (!HetaScratchCorr.empty() && !HetaCorr.empty()
        && HetaScratchCorr[0] && HetaCorr[0]
        && HetaScratchCorr[0]->Integral()>0 && HetaCorr[0]->Integral()>0)
    {
        TCanvas cFSvCl_Eta_FirstC("cFSvCl_Eta_First_WithCorr",
                                  "#Delta #eta from-scratch vs clusterizer (with corr) - first bin",
                                  1000, 750);
        std::vector<ResidualSpec> specs = {
            { HetaScratchCorr[0], "#Delta#eta  from scratch (corr)", (Color_t)(kBlue+1), (Style_t)20, 1 },
            { HetaCorr[0],        "#Delta#eta  clusterizer (corr)",  (Color_t)(kBlue+1), (Style_t)24, 2 }
        };
        drawResidualPanel(&cFSvCl_Eta_FirstC,
                          "From-scratch vs Clusterizer, With Position Correction",
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);
        cFSvCl_Eta_FirstC.SaveAs((dirFSvsClusWithCorr + "/DeltaEta_FSvsClus_FirstBin_Overlay.png").c_str());
        cFSvCl_Eta_FirstC.Close();
    }

    // Δφ - 4×2 table overlay (all energy bins), corrected
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Phi_Table_WithCorr",
                   "#Delta #phi from-scratch vs clusterizer (with corr) - table", 1600, 900);
        cT.SetTopMargin(0.10);
        cT.Divide(nCol, nRow, 0, 0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE) {
            if (!HphiScratchCorr[iE] || !HphiCorr[iE] ||
                HphiScratchCorr[iE]->Integral()==0 || HphiCorr[iE]->Integral()==0) { ++pad; continue; }

            cT.cd(pad++);
            setPadMargins();

            std::vector<ResidualSpec> specs = {
                { HphiScratchCorr[iE], "#Delta#phi  from scratch (corr)", (Color_t)(kRed+1),  (Style_t)20, 1 },
                { HphiCorr[iE],        "#Delta#phi  clusterizer (corr)",  (Color_t)(kRed+1),  (Style_t)24, 2 }
            };

            drawResidualPanel(gPad,
                              "From-scratch vs Clusterizer, With Position Correction",
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 0.85);
        }
        cT.cd(0);
        drawHeaderNDC("#Delta#phi  - From-scratch vs Clusterizer  (With Position Correction)", 0.975, 0.045);
        cT.SaveAs((dirFSvsClusWithCorr + "/DeltaPhi_FSvsClus_AllBins_Table.png").c_str());
        cT.Close();
    }

    // Δη - 4×2 table overlay (all energy bins), corrected
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Eta_Table_WithCorr",
                   "Δη from-scratch vs clusterizer (with corr) - table", 1600, 900);
        cT.SetTopMargin(0.10);
        cT.Divide(nCol, nRow, 0, 0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE) {
            if (!HetaScratchCorr[iE] || !HetaCorr[iE] ||
                HetaScratchCorr[iE]->Integral()==0 || HetaCorr[iE]->Integral()==0) { ++pad; continue; }

            cT.cd(pad++);
            setPadMargins();

            std::vector<ResidualSpec> specs = {
                { HetaScratchCorr[iE], "#Delta#eta  from scratch (corr)", (Color_t)(kBlue+1), (Style_t)20, 1 },
                { HetaCorr[iE],        "#Delta#eta  clusterizer (corr)",  (Color_t)(kBlue+1), (Style_t)24, 2 }
            };

            drawResidualPanel(gPad,
                              "From-scratch vs Clusterizer, With Position Correction",
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 0.85);
        }
        cT.cd(0);
        drawHeaderNDC("#Delta#eta  -  From-scratch vs Clusterizer  (With Position Correction)", 0.975, 0.045);
        cT.SaveAs((dirFSvsClusWithCorr + "/DeltaEta_FSvsClus_AllBins_Table.png").c_str());
        cT.Close();
    }

    // Δφ - μ/σ(E) overlay (two series): from‑scratch vs clusterizer, with corr
    if (!HphiScratchCorr.empty() && !HphiCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HphiScratchCorr, HphiCorr,
            kRed+1, 20, 24, "#mu",
            (dirFSvsClusWithCorr + "/DeltaPhi_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch (corr)", "clusterizer (corr)",
            "#Delta#phi  -  from scratch vs clusterizer  (With Position Correction)  -  Mean / RMS vs E");


    // Δη - μ/σ(E) overlay (two series): from‑scratch vs clusterizer, with corr
    if (!HetaScratchCorr.empty() && !HetaCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HetaScratchCorr, HetaCorr,
            kBlue+1, 20, 24, "#mu",
            (dirFSvsClusWithCorr + "/DeltaEta_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch (corr)", "clusterizer (corr)",
            "#Delta#eta  -  from scratch vs clusterizer  (With Position Correction)  -  Mean / RMS vs E");

    // NEW: 4-series μ/σ(E) overlay (from‑scratch vs clusterizer), with corr
    if (!HphiScratchCorr.empty() && !HphiCorr.empty()
        && !HetaScratchCorr.empty() && !HetaCorr.empty())
        makeMuSigmaFourSeriesOverlayFSvClus(HphiScratchCorr, HphiCorr,
                                            HetaScratchCorr, HetaCorr,
            (dirFSvsClusWithCorr + "/DeltaEtaPhi_FSvsClus_MeanSigmaVsE_Overlay_4Series.png").c_str(),
            "#Delta#phi / #Delta#eta  -  from scratch vs clusterizer  (With Position Correction)  -  Mean / RMS vs E");



    // =====================================================================
    //                   OVERLAYS  (φ=red, η=blue)
    // =====================================================================
    // ---- Overlay #1 : first-bin histogram (both residuals) ----------
    if (!Hphi.empty() && !Heta.empty() && Hphi[0] && Heta[0] &&
        Hphi[0]->Integral()>0 && Heta[0]->Integral()>0)
    {
        TCanvas cO1("cPlay_Overlay_First","DeltaEta/DeltaPhi first slice overlay",1000,750);

        std::vector<ResidualSpec> specs = {
            { Hphi[0], "#Delta#phi", (Color_t)(kRed+1),  (Style_t)20, 1 },
            { Heta[0], "#Delta#eta", (Color_t)(kBlue+1), (Style_t)20, 1 }
        };

        drawResidualPanel(&cO1, "Clusterizer Output, No Position Correction",
                          eEdges.front().first, eEdges.front().second,
                          specs, xMin, xMax, 1.0);

        std::string outO1 = dirOverlayNoCorr + "/DeltaEtaPhi_noCorrCluster_FirstBin_Overlay.png";
        cO1.SaveAs(outO1.c_str());
        cO1.Close();
        std::cout << "[Playground] Wrote " << outO1 << "\n";

    } else {
        std::cerr << "[Playground][WARN] First-bin overlay skipped (missing/empty hist).\n";
    }

    // ---- Overlay #2 : 4×2 table (both residuals on each pad) ----------
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cOT("cPlay_Overlay_Table","DeltaEta/DeltaPhi overlay - all slices",1600,900);
        cOT.SetTopMargin(0.10);
        cOT.Divide(nCol,nRow,0,0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE)
        {
            if (!Hphi[iE] || !Heta[iE] ||
                Hphi[iE]->Integral()==0 || Heta[iE]->Integral()==0) { ++pad; continue; }

            cOT.cd(pad++);
            setPadMargins();

            std::vector<ResidualSpec> specs = {
                { Hphi[iE], "#Delta#phi", (Color_t)(kRed+1),  (Style_t)20, 1 },
                { Heta[iE], "#Delta#eta", (Color_t)(kBlue+1), (Style_t)20, 1 }
            };

            drawResidualPanel(gPad, "Clusterizer Output, No Position Correction",
                              eEdges[iE].first, eEdges[iE].second,
                              specs, xMin, xMax, 0.85);
        }
        cOT.cd(0);
        drawHeaderNDC("#Delta#phi & #Delta#eta  -  Clusterizer Output, No Position Correction  (all energy bins)", 0.975, 0.045);

        std::string outOT = dirOverlayNoCorr + "/DeltaEtaPhi_noCorrCluster_AllBins_Table_Overlay.png";
        cOT.SaveAs(outOT.c_str());
        cOT.Close();
        std::cout << "[Playground] Wrote " << outOT << "\n";
    }

    // ---- Overlay #2b : 2×4 ratio table (Δφ/Δη per slice) ----------
    {
        const int nColR=2, nRowR=4, maxPadsR=nColR*nRowR;
        TCanvas cRT("cPlay_Ratio_Table","#Delta#phi / #Delta#eta - ratio (no corr, cluster) - all slices",1600,900);
        cRT.SetTopMargin(0.10);
        cRT.Divide(nColR,nRowR,0,0);

        int padR=1;
        for (std::size_t iE=0; iE<eEdges.size() && padR<=maxPadsR; ++iE)
        {
            if (!Hphi[iE] || !Heta[iE] ||
                Hphi[iE]->Integral()==0 || Heta[iE]->Integral()==0) { ++padR; continue; }

            cRT.cd(padR++);
            setPadMargins();

            TH1F* hP = cloneCountsHist(Hphi[iE]);
            TH1F* hE = cloneCountsHist(Heta[iE]);
            if (!hP || !hE) continue;

            // ensure proper error propagation (guard Sumw2 to avoid ROOT warnings)
            if (hP->GetSumw2N() == 0) hP->Sumw2();
            if (hE->GetSumw2N() == 0) hE->Sumw2();

            // build ratio histogram
            TH1F* hR = static_cast<TH1F*>(hP->Clone(Form("h_ratio_%zu", iE)));
            hR->SetDirectory(nullptr);
            hR->SetTitle("");
            if (hR->GetSumw2N() == 0) hR->Sumw2();

            hR->Divide(hE);

            // style & axes
            hR->SetMarkerStyle(20);
            hR->SetMarkerSize(0.8);
            hR->SetMarkerColor(kBlack);
            hR->SetLineColor(kBlack);

            hR->GetXaxis()->SetTitle("#Delta  [rad]");
            hR->GetYaxis()->SetTitle("#Delta#phi / #Delta#eta");
            hR->GetXaxis()->SetRangeUser(xMin, xMax);

            double yMaxR = computeYmaxCounts({hR});
            if (yMaxR <= 0.) yMaxR = 1.0;
            hR->GetYaxis()->SetRangeUser(0., 1.20*yMaxR);
            hR->Draw("E");

            // reference line at ratio = 1
            TLine l1(xMin, 1.0, xMax, 1.0);
            l1.SetLineStyle(2);
            l1.Draw();

            // energy label (NDC)
            drawNDCEnergy(0.96, 0.92, eEdges[iE].first, eEdges[iE].second, 0.042);

            // legend (NDC)
            TLegend* lgR = makeNDCLegend(0.70,0.78,0.93,0.90);
            lgR->AddEntry(hR, "#Delta#phi / #Delta#eta", "p");
            lgR->Draw();

            gPad->Modified();
            gPad->Update();
        }

        cRT.cd(0);
        drawHeaderNDC("#Delta#phi / #Delta#eta  -  Clusterizer Output, No Position Correction  (ratio in each energy bin)", 0.975, 0.045);

        std::string outRT = dirOverlayNoCorr + "/DeltaEtaPhi_noCorrCluster_AllBins_Table_Ratio.png";
        cRT.SaveAs(outRT.c_str());
        cRT.Close();
        std::cout << "[Playground] Wrote " << outRT << "\n";
    }


    // ---- Overlay #3 : μ(E) & σ(E) (both residuals on shared pads) ----
    {
        std::vector<double> eCtr;
        std::vector<double> muP, dmuP, sgP, dsgP;  // phi
        std::vector<double> muE, dmuE, sgE, dsgE;  // eta

        for (std::size_t iE=0;iE<eEdges.size();++iE)
        {
            if (!Hphi[iE] || !Heta[iE] ||
                Hphi[iE]->Integral()==0 || Heta[iE]->Integral()==0)
                continue;

            eCtr.push_back( 0.5*(eEdges[iE].first + eEdges[iE].second) );

            double mP, meP, rP, reP, mE, meE, rE, reE;
            std::tie(mP, meP, rP, reP) = statsFromHistRange(Hphi[iE], xMin, xMax);
            std::tie(mE, meE, rE, reE) = statsFromHistRange(Heta[iE], xMin, xMax);

            muP.push_back(mP); dmuP.push_back(meP);
            sgP.push_back(rP); dsgP.push_back(reP);

            muE.push_back(mE); dmuE.push_back(meE);
            sgE.push_back(rE); dsgE.push_back(reE);
        }

        if (!eCtr.empty())
        {
            const double xAxisMin = 0.0;
            const double xAxisMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
            std::vector<double> ex(eCtr.size(), 0.0);

            TCanvas cS("cPlay_Overlay_MuSig","DeltaEta/DeltaPhi overlay - Mean / RMS vs E",900,800);

            auto pads = makeEqualHeightTwinPads(cS, "OverlayMuSig", 0.15,0.05, 0.14,0.02, 0.04,0.16);
            TPad* pTop = pads.first; TPad* pBot = pads.second;

            // μ(E)
            pTop->cd();
            double muLo=1e30, muHi=-1e30;
            for (std::size_t i=0;i<muP.size();++i){
                muLo = std::min(muLo, std::min(muP[i]-dmuP[i], muE[i]-dmuE[i]));
                muHi = std::max(muHi, std::max(muP[i]+dmuP[i], muE[i]+dmuE[i]));
            }
            const double padMu = 0.25*(muHi-muLo);

            TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
            frU.SetMinimum(muLo - padMu);
            frU.SetMaximum(muHi + padMu);
            frU.GetXaxis()->SetLabelSize(0.0);
            frU.GetXaxis()->SetTitleSize(0.0);
            frU.GetXaxis()->SetTickLength(0.0);
            frU.Draw("AXIS");
            TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

            TGraphErrors gMuP(eCtr.size(), eCtr.data(), muP.data(), ex.data(), dmuP.data());
            gMuP.SetMarkerStyle(21); gMuP.SetMarkerColor(kRed+1);
            gMuP.SetLineColor  (kRed+1); gMuP.SetMarkerSize(1.0);
            gMuP.Draw("P SAME");

            TGraphErrors gMuE(eCtr.size(), eCtr.data(), muE.data(), ex.data(), dmuE.data());
            gMuE.SetMarkerStyle(20); gMuE.SetMarkerColor(kBlue+1);
            gMuE.SetLineColor  (kBlue+1); gMuE.SetMarkerSize(1.0);
            gMuE.Draw("P SAME");

            TLegend legU(0.8,0.75,0.9,0.85);
            legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.035);
            legU.AddEntry(&gMuP,"#Delta#phi","p");
            legU.AddEntry(&gMuE,"#Delta#eta","p");
            legU.Draw();

            // σ(E)
            pBot->cd();

            double sgHi=-1e30;
            for (double s : sgP) sgHi = std::max(sgHi,s);
            for (double s : sgE) sgHi = std::max(sgHi,s);

            TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
            frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi);
            frL.Draw("AXIS");

            TGraphErrors gSgP(eCtr.size(), eCtr.data(), sgP.data(), ex.data(), dsgP.data());
            gSgP.SetMarkerStyle(21); gSgP.SetMarkerColor(kRed+1);
            gSgP.SetLineColor  (kRed+1); gSgP.SetMarkerSize(1.0);
            gSgP.Draw("P SAME");

            TGraphErrors gSgE(eCtr.size(), eCtr.data(), sgE.data(), ex.data(), dsgE.data());
            gSgE.SetMarkerStyle(20); gSgE.SetMarkerColor(kBlue+1);
            gSgE.SetLineColor  (kBlue+1); gSgE.SetMarkerSize(1.0);
            gSgE.Draw("P SAME");

            // Header
            cS.cd();
            TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.038);
            h.DrawLatex(0.50,0.97,"#Delta#phi & #Delta#eta - Mean / RMS vs E, Clusterizer with No Position Correction");


            std::string outOS = dirOverlayNoCorr + "/DeltaEtaPhi_noCorrCluster_MeanSigmaVsE_Overlay.png";
            cS.SaveAs(outOS.c_str());
            cS.Close();
            std::cout << "[Playground] Wrote " << outOS << "\n";
        } else {
            std::cerr << "[Playground][WARN] Overlay μ/σ summary skipped (no common non-empty slices).\n";
        }
    }

    // ================================================================
    //  Δφ variants - μ/σ vs E overlays + μ vs ln(E) fits (as before)
    // ================================================================
    {
        const std::string dirPhiLog = std::string(outDir) + "/PhiLogFitSummary";
        ensureDir(dirPhiLog.c_str());
        // ---- All Δφ variants: PDC (raw/corr), CLUS (raw/corr/b-corr), and 4× EA flavours
        constexpr int NPHI = 9;

        const std::array<const char*,NPHI> phiPat = {
            "h_phi_diff_raw_%s",                       // 0: PDC raw
            "h_phi_diff_corr_%s",                      // 1: PDC corrected
            "h_phi_diff_cpRaw_%s",                     // 2: CLUS raw (no corr)
            "h_phi_diff_cpCorr_%s",                    // 3: CLUS CorrectPosition
            "h_phi_diff_cpBcorr_%s",                   // 4: CLUS b(E) corr
            "h_phi_diff_cpCorrEA_geom_%s",             // 5: CLUS-CP(|vz| + E)
            "h_phi_diff_cpCorrEA_fitEtaDep_%s",        // 6: CLUS-CP(EA |eta|+E fits)
            "h_phi_diff_cpCorrEA_fitEnergyOnly_%s",    // 7: CLUS-CP(EA E-only fits)
            "h_phi_diff_cpCorrEA_fitPhiE_etaEtaDep_%s" // 8: CLUS-CP(EA φ:E-only, η:|η|+E)
        };
        const std::array<const char*,NPHI> phiLab = {
            "PDC raw",
            "PDC corrected",
            "no corr, cluster",
            "CorrectPosition, cluster",
            "b(E) corr, cluster",
            "CorrectPosition(EA geom), cluster",
            "CorrectPosition(EA |#eta|+E), cluster",
            "CorrectPosition(EA E-only), cluster",
            "CorrectPosition(EA #varphi:E-only, #eta:|#eta|+E), cluster"
        };
        const std::array<const char*,NPHI> phiEnum = {
            "ETiltVariant::PDC_RAW",                 // 0
            "ETiltVariant::PDC_CORR",                // 1
            "ETiltVariant::CLUS_RAW",                // 2
            "ETiltVariant::CLUS_CP",                 // 3
            "ETiltVariant::CLUS_BCORR",              // 4
            "ETiltVariant::CLUS_CP_EA_FIT_ZDEP",     // 5
            "ETiltVariant::CLUS_CP_EA_FIT_ETADEP",   // 6
            "ETiltVariant::CLUS_CP_EA_FIT_EONLY",    // 7
            "ETiltVariant::CLUS_CP_EA_MIX"           // 8
        };


        std::array<std::vector<TH1F*>,NPHI> PHI{};
        for (int v=0; v<NPHI; ++v) PHI[v] = loadSet(phiPat[v]);

        std::vector<double> eCtr; eCtr.reserve(eEdges.size());
        for (std::size_t iE=0; iE<eEdges.size(); ++iE)
            eCtr.push_back(0.5*(eEdges[iE].first + eEdges[iE].second));

        std::array<std::vector<double>,NPHI> MU, dMU, SG, dSG;
        std::vector<int> shown; shown.reserve(NPHI);

        for (int v=0; v<NPHI; ++v) {
            MU[v].reserve(eEdges.size());
            dMU[v].reserve(eEdges.size());
            SG[v].reserve(eEdges.size());
            dSG[v].reserve(eEdges.size());

            bool any = false;
            for (std::size_t iE=0; iE<eEdges.size(); ++iE) {
                if (!PHI[v][iE] || PHI[v][iE]->Integral()==0) {
                    MU[v].push_back(0.0);  dMU[v].push_back(0.0);
                    SG[v].push_back(0.0);  dSG[v].push_back(0.0);
                    continue;
                }
                any = true;
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(PHI[v][iE], xMin, xMax);
                MU[v].push_back(m);   dMU[v].push_back(me);
                SG[v].push_back(r);   dSG[v].push_back(re);
            }
            if (any) shown.push_back(v);
        }

        // Local style maps for up to 9 series
        const Color_t colByV[NPHI] = {
            kBlue+1,  // PDC raw
            kBlue+1,  // PDC corrected
            kRed+1,   // CLUS raw
            kRed+1,   // CLUS CP
            kRed+1,   // CLUS b(E)
            kRed+3,   // EA geom
            kRed+2,   // EA |eta|+E
            kRed+4,   // EA E-only
            kMagenta+1// EA mix
        };
        const Style_t mkByV [NPHI] = {
            20,       // PDC raw     (filled)
            24,       // PDC corr    (open)
            20,       // CLUS raw    (filled)
            24,       // CLUS CP     (open)
            21,       // CLUS b(E)   (filled square)
            24,       // EA geom     (open)
            25,       // EA |eta|+E  (open diamond)
            27,       // EA E-only   (open triangle)
            28        // EA mix      (open star)
        };



        // μ/σ vs E overlay (all available variants)
        if (!shown.empty()) {
            const double xMinE = 0.0;
            const double xMaxE = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
            std::vector<double> ex(eCtr.size(), 0.0);

            TCanvas c("cPhiVar_MuSig","DeltaPhi - variants μ/σ vs E",1000,850);
            auto pads = makeEqualHeightTwinPads(c, "PhiVarMuSig", 0.15,0.05, 0.14,0.02, 0.04,0.16);
            TPad* pTop = pads.first; TPad* pBot = pads.second;

            // Top pad: μ(E)
            pTop->cd();
            double muLo=1e30, muHi=-1e30;
            for (int v : shown)
                for (std::size_t i=0; i<MU[v].size(); ++i) {
                    muLo = std::min(muLo, MU[v][i] - dMU[v][i]);
                    muHi = std::max(muHi, MU[v][i] + dMU[v][i]);
                }
            const double padMu = 0.25*(muHi-muLo);
            TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
            frU.SetMinimum(muLo - padMu);
            frU.SetMaximum(muHi + padMu);
            frU.GetXaxis()->SetLabelSize(0.0);
            frU.GetXaxis()->SetTitleSize(0.0);
            frU.GetXaxis()->SetTickLength(0.0);
            frU.Draw("AXIS");
            TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

            TLegend legU(0.62,0.74,0.93,0.92);
            legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.032);

            std::vector<std::unique_ptr<TGraphErrors>> gMu;
            for (int v : shown) {
                auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                        eCtr.data(), MU[v].data(),
                                                        ex.data(),   dMU[v].data());
                g->SetMarkerStyle(mkByV[v]);
                g->SetMarkerColor(colByV[v]);
                g->SetLineColor  (colByV[v]);
                g->SetMarkerSize(1.0);
                g->Draw("P SAME");
                legU.AddEntry(g.get(), phiLab[v], "p");
                gMu.emplace_back(std::move(g));
            }
            legU.Draw();

            // Bottom pad: σ(E)
            pBot->cd();
            double sgHi = -1e30;
            for (int v : shown) for (double s : SG[v]) sgHi = std::max(sgHi, s);
            TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
            frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

            std::vector<std::unique_ptr<TGraphErrors>> gSg;
            for (int v : shown) {
                auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                        eCtr.data(), SG[v].data(),
                                                        ex.data(),   dSG[v].data());
                g->SetMarkerStyle(mkByV[v]);
                g->SetMarkerColor(colByV[v]);
                g->SetLineColor  (colByV[v]);
                g->SetMarkerSize(1.0);
                g->Draw("P SAME");
                gSg.emplace_back(std::move(g));
            }

            c.cd();
            TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.040);
            h.DrawLatex(0.50,0.985,"#Delta#phi  -  Mean / RMS vs E  (variants)");
            const std::string outP = dirPhiLog + "/DeltaPhi_MeanSigmaVsE_9Variants.png";

            c.SaveAs(outP.c_str());
            c.Close();

            // CSV: include only shown variants, labels from phiLab[v]
            std::vector<std::string> Vlab; Vlab.reserve(shown.size());
            std::vector<std::vector<double>> EctrS, MUv, dMUv, SGv, dSGv;
            for (int v : shown){
                Vlab.emplace_back(phiLab[v]);
                EctrS.push_back(eCtr);
                MUv  .push_back(MU[v]);   dMUv.push_back(dMU[v]);
                SGv  .push_back(SG[v]);   dSGv.push_back(dSG[v]);
            }
            saveMuSigmaCSV(outP,
                           Vlab,
                           std::vector<std::string>(Vlab.size(), "phi"),
                           EctrS, MUv, dMUv, SGv, dSGv);

            // ---------------- μ vs ln(E) overlay and fits: ALL variants (no CLUS_CP_EA alias) ----------------
            if (!eCtr.empty()) {
                // Build ln(E)
                std::vector<double> lnEAll(eCtr.size());
                std::transform(eCtr.begin(), eCtr.end(), lnEAll.begin(),
                               [](double e){ return std::log(e); });

                const double xLoAll = *std::min_element(lnEAll.begin(), lnEAll.end()) - 0.05;
                const double xHiAll = *std::max_element(lnEAll.begin(), lnEAll.end()) + 0.05;

                // y-range across all MU[v]
                double yLoAll = +1e30, yHiAll = -1e30;
                for (int v = 0; v < NPHI; ++v) {
                    if (MU[v].empty()) continue;
                    yLoAll = std::min(yLoAll, *std::min_element(MU[v].begin(), MU[v].end()));
                    yHiAll = std::max(yHiAll, *std::max_element(MU[v].begin(), MU[v].end()));
                }
                const double padAll = 0.15*(yHiAll - yLoAll);
                yLoAll -= padAll;
                yHiAll += 0.35*(yHiAll - yLoAll);

                TCanvas cAll("cPhi_AllVar_mu_vs_lnE","Δφ - μ vs ln(E) (ALL variants)",1100,780);
                cAll.SetLeftMargin(0.15); cAll.SetRightMargin(0.06);

                TH1F frAll("frAll",";ln E  [GeV];#mu  [rad]",1,xLoAll,xHiAll);
                frAll.SetMinimum(yLoAll); frAll.SetMaximum(yHiAll);
                frAll.Draw("AXIS");

                std::vector<double> ex0All(eCtr.size(), 0.0);
                TLegend lgAll(0.17,0.72,0.92,0.90);
                lgAll.SetBorderSize(0); lgAll.SetFillStyle(0); lgAll.SetTextSize(0.030);
                lgAll.SetNColumns(2);

                // -------- Fit all variants, cache (intercept, slope), then write in required order --------
                std::array<double,NPHI> fitB{};     // intercepts
                std::array<double,NPHI> fitM{};     // slopes
                std::array<bool,  NPHI> haveFit{};  // availability

                for (int v = 0; v < NPHI; ++v) {
                    bool has = false;
                    for (auto* h : PHI[v]) { if (h && h->Integral()>0) { has = true; break; } }
                    if (!has) { haveFit[v] = false; continue; }

                    TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                       lnEAll.data(), MU[v].data(),
                                                       ex0All.data(), dMU[v].data());
                    g->SetMarkerStyle(mkByV[v]);
                    g->SetMarkerColor(colByV[v]);
                    g->SetLineColor  (colByV[v]);
                    g->SetMarkerSize(1.05);
                    g->Draw("P SAME");

                    TF1 f(Form("fAll_%d",v), "pol1", xLoAll, xHiAll);
                    f.SetLineColor (colByV[v]);
                    f.SetLineStyle (2);
                    f.SetLineWidth (2);
                    g->Fit(&f, "Q");
                    f.Draw("SAME");

                    fitB[v]    = f.GetParameter(0);   // intercept
                    fitM[v]    = f.GetParameter(1);   // slope
                    haveFit[v] = true;

                    lgAll.AddEntry(g, Form("%s  (m=%.2e, b=%.2e)", phiLab[v], fitM[v], fitB[v]), "p");
                }

                // Open the file fresh and write in the exact order & format you asked for
                std::ofstream fAll(dirPhiLog + "/DeltaPhi_MuVsLogE_fit.txt", std::ios::trunc);
                fAll << "# ---- Δφ: μ vs ln(E) linear fits (intercept first, slope second) ----\n";

                // Pretty comment strings exactly as in your spec
                const std::array<const char*,NPHI> outComment = {
                    "no corr, scratch",                              // 0: PDC raw
                    "b(E) corr, scratch",                            // 1: PDC corrected
                    "no corr, cluster",                              // 2: CLUS raw
                    "CorrectPosition, cluster",                      // 3: CLUS CP
                    "b(E) corr, cluster",                            // 4: CLUS b(E)
                    "CorrectPosition(EA geom), cluster",             // 5
                    "CorrectPosition(EA |#eta|+E), cluster",         // 6
                    "CorrectPosition(EA E-only), cluster",           // 7
                    "CorrectPosition(EA #varphi:E-only, #eta:|#eta|+E), cluster" // 8
                };

                // Required emission order:
                //   CLUS_RAW, CLUS_CP, EA geom, EA |eta|+E, EA E-only, EA mix, CLUS_BCORR, PDC_RAW, PDC_CORR
                const int vCLUSraw  = 2;
                const int vCLUScp   = 3;
                const int vEAgeom   = 5;
                const int vEAetaE   = 6;
                const int vEAeOnly  = 7;
                const int vEAmix    = 8;
                const int vBCorr    = 4;
                const int vPDCraw   = 0;
                const int vPDCcorr  = 1;

                const std::vector<int> writeOrder = {
                    vCLUSraw, vCLUScp, vEAgeom, vEAetaE, vEAeOnly, vEAmix, vBCorr, vPDCraw, vPDCcorr
                };

                for (int v : writeOrder) {
                    if (!haveFit[v]) continue;  // skip if that variant has no data
                    // Write as { intercept, slope } so the RHS entry (slope) is non-negative in your case
                    fAll << "  {" << phiEnum[v] << ", { "
                         << Form("%.6e", fitB[v]) << ", " << Form("%.6e", fitM[v])
                         << " }},  // " << outComment[v] << "\n";
                }

                fAll.close();
                lgAll.Draw();
                cAll.SaveAs((dirPhiLog + "/DeltaPhi_MuVsLogE_AllVariants.png").c_str());
                cAll.Close();
            }
        }


        // ==================== CLUS-only overlays (RAW, CORR, CORR(EA_geom)) + EA fits ====================
        {
            // indices in our expanded variant table
            const int vCLUSraw    = 2; // no corr, cluster
            const int vCLUScorr   = 3; // CorrectPosition, cluster
            const int vEA_geom    = 5; // EA geometry
            const int vEA_etaDep  = 6; // EA |eta|+E fits
            const int vEA_eOnly   = 7; // EA E-only fits
            const int vEA_mix     = 8; // EA mixed (phi:E-only, eta:|eta|+E)

            auto hasAny = [&](int vidx)->bool {
                for (auto* h : PHI[vidx]) if (h && h->Integral()>0) return true;
                return false;
            };

            // 3-way plot uses: CLUS RAW / CLUS CP / EA_geom as representative
            std::vector<int> use3;
            for (int v : {vCLUSraw, vCLUScorr, vEA_geom}) if (hasAny(v)) use3.push_back(v);

            if (!use3.empty() && eCtr.size() >= 2) {
                // ---------------- μ/σ vs E (3-way) ----------------
                const double xMinE = 0.0;
                const double xMaxE = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
                std::vector<double> ex(eCtr.size(), 0.0);

                TCanvas c3("cPhiVar_MuSig_3","DeltaPhi - CLUS RAW/CORR/CORR(EA_geom) μ/σ vs E",1000,850);
                auto pads3 = makeEqualHeightTwinPads(c3, "PhiVarMuSig3", 0.15,0.05, 0.14,0.02, 0.04,0.16);
                TPad* pTop3 = pads3.first; TPad* pBot3 = pads3.second;

                // --- Top: μ(E)
                pTop3->cd();
                double muLo=+1e30, muHi=-1e30;
                for (int v : use3) {
                    for (std::size_t i=0; i<MU[v].size(); ++i) {
                        muLo = std::min(muLo, MU[v][i] - dMU[v][i]);
                        muHi = std::max(muHi, MU[v][i] + dMU[v][i]);
                    }
                }
                const double padMu = 0.25*(muHi - muLo);
                TH1F frU3("frU3",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
                frU3.SetMinimum(muLo - padMu);
                frU3.SetMaximum(muHi + padMu);
                frU3.GetXaxis()->SetLabelSize(0.0);
                frU3.GetXaxis()->SetTitleSize(0.0);
                frU3.GetXaxis()->SetTickLength(0.0);
                frU3.Draw("AXIS");
                TLine l0_3(xMinE,0.0,xMaxE,0.0); l0_3.SetLineStyle(2); l0_3.Draw();

                TLegend legU3(0.60,0.74,0.93,0.92);
                legU3.SetBorderSize(0); legU3.SetFillStyle(0); legU3.SetTextSize(0.033);

                std::vector<std::unique_ptr<TGraphErrors>> gMu3;
                gMu3.reserve(use3.size());
                for (int v : use3) {
                    auto g = std::make_unique<TGraphErrors>(
                        (int)eCtr.size(), eCtr.data(), MU[v].data(), ex.data(), dMU[v].data());
                    g->SetMarkerStyle(mkByV[v]);
                    g->SetMarkerColor(colByV[v]);
                    g->SetLineColor  (colByV[v]);
                    g->SetMarkerSize(1.05);
                    g->Draw("P SAME");
                    legU3.AddEntry(g.get(), phiLab[v], "p");
                    gMu3.emplace_back(std::move(g));
                }
                legU3.Draw();

                // --- Bottom: σ(E)
                pBot3->cd();
                double sgHi=-1e30;
                for (int v : use3) for (double s : SG[v]) sgHi = std::max(sgHi, s);
                TH1F frL3("frL3",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
                frL3.SetMinimum(0.0); frL3.SetMaximum(1.15*sgHi); frL3.Draw("AXIS");

                std::vector<std::unique_ptr<TGraphErrors>> gSg3;
                gSg3.reserve(use3.size());
                for (int v : use3) {
                    auto g = std::make_unique<TGraphErrors>(
                        (int)eCtr.size(), eCtr.data(), SG[v].data(), ex.data(), dSG[v].data());
                    g->SetMarkerStyle(mkByV[v]);
                    g->SetMarkerColor(colByV[v]);
                    g->SetLineColor  (colByV[v]);
                    g->SetMarkerSize(1.05);
                    g->Draw("P SAME");
                    gSg3.emplace_back(std::move(g));
                }

                c3.cd();
                TLatex h3; h3.SetNDC(); h3.SetTextAlign(22); h3.SetTextSize(0.038);
                h3.DrawLatex(0.50,0.985,"#Delta#phi  -  Clusterizer RAW / CORR / CORR(EA_{geom})  -  Mean / RMS vs E");

                const std::string outP3 = dirPhiLog + "/DeltaPhi_MeanSigmaVsE_Clusterizer3.png";
                c3.SaveAs(outP3.c_str());
                c3.Close();

                // ---------------- μ vs ln(E) with fits (3-way) ----------------
                std::vector<double> lnE(eCtr.size());
                std::transform(eCtr.begin(), eCtr.end(), lnE.begin(),
                               [](double e){ return std::log(e); });

                const double xLo = *std::min_element(lnE.begin(),lnE.end()) - 0.05;
                const double xHi = *std::max_element(lnE.begin(),lnE.end()) + 0.05;

                double yLo=+1e30, yHi=-1e30;
                for (int v : use3) {
                    yLo = std::min(yLo, *std::min_element(MU[v].begin(), MU[v].end()));
                    yHi = std::max(yHi, *std::max_element(MU[v].begin(), MU[v].end()));
                }
                const double pad=0.15*(yHi-yLo); yLo-=pad; yHi+=0.35*(yHi-yLo);

                TCanvas c3ln("cLn3","μ vs lnE (CLUS RAW/CORR/CORR(EA_geom))",900,640);
                c3ln.SetLeftMargin(0.15); c3ln.SetRightMargin(0.06);

                TH1F fr3("",";ln E  [GeV];#mu  [rad]",1,xLo,xHi);
                fr3.SetMinimum(yLo); fr3.SetMaximum(yHi); fr3.Draw("AXIS");

                std::vector<double> ex0(eCtr.size(), 0.0);
                TLegend lg3(0.17,0.74,0.92,0.90);
                lg3.SetBorderSize(0); lg3.SetFillStyle(0); lg3.SetTextSize(0.032);

                // Draw the three series with linear fits; no file writes here (consolidated block handles all).
                for (int v : use3) {
                    TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                       lnE.data(), MU[v].data(),
                                                       ex0.data(), dMU[v].data());
                    g->SetMarkerStyle(mkByV[v]);
                    g->SetMarkerColor(colByV[v]);
                    g->SetLineColor  (colByV[v]);
                    g->SetMarkerSize(1.1);
                    g->Draw("P SAME");

                    TF1 f(Form("f3_%d",v), "pol1", xLo, xHi);
                    f.SetLineColor (colByV[v]);
                    f.SetLineStyle (2);
                    f.SetLineWidth (2);
                    g->Fit(&f, "Q");
                    f.Draw("SAME");

                    // Legend as before
                    lg3.AddEntry(g, Form("%s  (m=%.2e, b=%.2e)",
                                         phiLab[v], f.GetParameter(1), f.GetParameter(0)), "p");
                }

                lg3.Draw();
                c3ln.SaveAs((dirPhiLog + "/DeltaPhi_MuVsLogE_Clusterizer3.png").c_str());
                c3ln.Close();
            }
        }

        // μ vs ln(E) - FOUR variants (CLUS-RAW, CLUS-CORR, PDC-RAW, PDC-CORR)
        // Also produce μ/σ vs E overlay with the same four series and required styles.
        {
            // Variant indices in PHI / phiLab ordering above
            const int vPDCraw  = 0; // "PDC raw"               -> Variable-b (Benchmark) Uncorr
            const int vPDCcorr = 1; // "PDC corrected"         -> Variable-b (Benchmark) Corr
            const int vCLUSraw = 2; // "no corr, cluster"      -> Constant-b (Production) Uncorr
            const int vCLUScor = 3; // "CorrectPosition, cluster" -> Constant-b (Production) Corr

            // Which series actually have content?
            auto hasAny = [&](int vidx)->bool {
                for (auto* h : PHI[vidx]) if (h && h->Integral()>0) return true;
                return false;
            };
            std::vector<int> useV;
            for (int v : {vCLUSraw, vCLUScor, vPDCraw, vPDCcorr})
                if (hasAny(v)) useV.push_back(v);

            if (!useV.empty() && eCtr.size() >= 2) {
                // ----------- Common helpers (styles, labels) -----------
                auto seriesColor = [&](int v)->Color_t {
                    // Clusterizer = red; PDC = blue
                    return (v==vCLUSraw || v==vCLUScor) ? (kRed+1) : (kBlue+1);
                };
                auto seriesMarker = [&](int v)->Style_t {
                    // RAW = closed circle (20); CORR = open circle (24)
                    return (v==vCLUSraw || v==vPDCraw) ? (Style_t)20 : (Style_t)24;
                };
                auto prettyLabel = [&](int v)->const char* {
                    if (v==vPDCraw)  return "Variable-b (Benchmark) Uncorr";
                    if (v==vPDCcorr) return "Variable-b (Benchmark) Corr";
                    if (v==vCLUSraw) return "Constant-b (Production) Uncorr";
                    if (v==vCLUScor) return "Constant-b (Production) Corr";
                    return phiLab[v];
                };

                // ==================== Plot 1: μ vs ln(E) (4 series) ====================
                std::vector<double> lnE(eCtr.size());
                std::transform(eCtr.begin(), eCtr.end(), lnE.begin(),
                               [](double e){ return std::log(e); });

                const double xLo = *std::min_element(lnE.begin(),lnE.end()) - 0.05;
                const double xHi = *std::max_element(lnE.begin(),lnE.end()) + 0.05;

                double yLo=+1e30, yHi=-1e30;
                for (int v : useV) {
                    yLo = std::min(yLo, *std::min_element(MU[v].begin(), MU[v].end()));
                    yHi = std::max(yHi, *std::max_element(MU[v].begin(), MU[v].end()));
                }
                const double pad=0.15*(yHi-yLo); yLo-=pad; yHi+=0.35*(yHi-yLo);

                TCanvas c4ln("cLn4","μ vs lnE (four variants)",900,640);
                c4ln.SetLeftMargin(0.15); c4ln.SetRightMargin(0.06);

                TH1F frLn("",";ln E  [GeV];#mu  [rad]",1,xLo,xHi);
                frLn.SetMinimum(yLo); frLn.SetMaximum(yHi); frLn.Draw("AXIS");

                std::vector<double> ex0(eCtr.size(), 0.0);

                TLegend lg4(0.175,0.72,0.45,0.91);
                lg4.SetBorderSize(0); lg4.SetFillStyle(0); lg4.SetTextSize(0.03);

                // Keep local TF1 objects alive by storing them
                std::vector<std::unique_ptr<TF1>> fits;
                fits.reserve(useV.size());

                // Plot-only fits (no file writes here; consolidated ALL-variants block handles the file).
                for (int v : useV) {
                    TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                       lnE.data(), MU[v].data(),
                                                       ex0.data(), dMU[v].data());
                    g->SetMarkerStyle(seriesMarker(v));
                    g->SetMarkerColor(seriesColor(v));
                    g->SetLineColor  (seriesColor(v));
                    g->SetMarkerSize(1.1);
                    g->Draw("P SAME");

                    TF1* f = new TF1(Form("f_ln_%d", v), "pol1", xLo, xHi);
                    f->SetLineColor (seriesColor(v));
                    f->SetLineStyle (2);
                    f->SetLineWidth (2);
                    g->Fit(f, "Q");
                    f->Draw("SAME");
                    fits.emplace_back(f);

                    const double b = f->GetParameter(0); // intercept
                    const double m = f->GetParameter(1); // slope

                    // Legend unchanged
                    lg4.AddEntry(g, Form("%s  (m=%.2e, b=%.2e)", prettyLabel(v), m, b), "p");
                }

                // TLatex summary
                TLatex t4; t4.SetNDC(); t4.SetTextSize(0.044); t4.SetTextAlign(13);
                t4.DrawLatex(0.67, 0.24, "#mu(E) = b + m #upoint ln E");

                lg4.Draw();

                c4ln.SaveAs((dirPhiLog + "/DeltaPhi_MuVsLogE_RawVsCorr.png").c_str());
                c4ln.Close();

                // ==================== Plot 2: μ(E) & σ(E) (4 series) ====================
                const double xMinE = 0.0;
                const double xMaxE = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
                std::vector<double> exE(eCtr.size(), 0.0);

                TCanvas c4("cPhiVar_MuSig_4","DeltaPhi - four variants μ/σ vs E",1000,850);
                auto pads4 = makeEqualHeightTwinPads(c4, "PhiVarMuSig4", 0.15,0.05, 0.14,0.02, 0.04,0.16);
                TPad* pTop4 = pads4.first; TPad* pBot4 = pads4.second;

                // Top: μ(E)
                pTop4->cd();
                double muLo=+1e30, muHi=-1e30;
                for (int v : useV) {
                    for (std::size_t i=0;i<eCtr.size();++i) {
                        muLo = std::min(muLo, MU[v][i] - dMU[v][i]);
                        muHi = std::max(muHi, MU[v][i] + dMU[v][i]);
                    }
                }
                const double padMu = 0.25*(muHi - muLo);
                TH1F frU4("frU4",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
                frU4.SetMinimum(muLo - padMu);
                frU4.SetMaximum(muHi + padMu);
                frU4.GetXaxis()->SetLabelSize(0.0);
                frU4.GetXaxis()->SetTitleSize(0.0);
                frU4.GetXaxis()->SetTickLength(0.0);
                frU4.Draw("AXIS");
                TLine l0U4(xMinE,0.0,xMaxE,0.0); l0U4.SetLineStyle(2); l0U4.Draw();

                // Keep graph objects alive; also capture specific variant pointers for legend ordering
                std::vector<std::unique_ptr<TGraphErrors>> gMu4;
                gMu4.reserve(useV.size());

                TGraphErrors* gTC_Uncorr = nullptr;  // "Constant-b (Production) Uncorr"
                TGraphErrors* gBL_Uncorr = nullptr;  // "Variable-b (Benchmark) Uncorr"
                TGraphErrors* gTC_Corr   = nullptr;  // "Constant-b (Production) Corr"
                TGraphErrors* gBL_Corr   = nullptr;  // "Variable-b (Benchmark) Corr"

                for (int v : useV) {
                    auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                            eCtr.data(), MU[v].data(),
                                                            exE.data(), dMU[v].data());
                    g->SetMarkerStyle(seriesMarker(v));     // RAW=20 (closed), CORR=24 (open)
                    g->SetMarkerColor(seriesColor(v));      // CLUS=red, PDC=blue
                    g->SetLineColor  (seriesColor(v));
                    g->SetMarkerSize(1.0);
                    g->Draw("P SAME");

                    const std::string lab = prettyLabel(v);
                    if (lab == "Constant-b (Production) Uncorr")         gTC_Uncorr = g.get();
                    else if (lab == "Variable-b (Benchmark) Uncorr") gBL_Uncorr = g.get();
                    else if (lab == "Constant-b (Production) Corr")      gTC_Corr   = g.get();
                    else if (lab == "Variable-b (Benchmark) Corr")   gBL_Corr   = g.get();

                    gMu4.emplace_back(std::move(g));
                }

                // Legend moved to the μ(E) (top) panel, upper-right; row-major order
                pTop4->cd();  // ensure legend is drawn on the top pad
                TLegend legS4(0.17, 0.72, 0.7, 0.87);   // upper-right in top pad
                legS4.SetBorderSize(0);
                legS4.SetFillStyle(0);
                legS4.SetTextSize(0.039);
                legS4.SetNColumns(2);

                // Row 1
                if (gTC_Uncorr) legS4.AddEntry(gTC_Uncorr, "Constant-b (Production) Uncorr", "p");
                if (gBL_Uncorr) legS4.AddEntry(gBL_Uncorr, "Variable-b (Benchmark) Uncorr", "p");
                // Row 2
                if (gTC_Corr)   legS4.AddEntry(gTC_Corr,   "Constant-b (Production) Corr", "p");
                if (gBL_Corr)   legS4.AddEntry(gBL_Corr,   "Variable-b (Benchmark) Corr", "p");

                legS4.Draw();

                // Bottom: σ(E)
                pBot4->cd();
                double sgHi=-1e30;
                for (int v : useV) for (double s : SG[v]) sgHi = std::max(sgHi, s);
                TH1F frL4("frL4",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
                frL4.SetMinimum(0.0); frL4.SetMaximum(1.15*sgHi); frL4.Draw("AXIS");

                std::vector<std::unique_ptr<TGraphErrors>> gSg4;
                gSg4.reserve(useV.size());
                for (int v : useV) {
                    auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                            eCtr.data(), SG[v].data(),
                                                            exE.data(), dSG[v].data());
                    g->SetMarkerStyle(seriesMarker(v));
                    g->SetMarkerColor(seriesColor(v));
                    g->SetLineColor  (seriesColor(v));
                    g->SetMarkerSize(1.0);
                    g->Draw("P SAME");
                    gSg4.emplace_back(std::move(g));
                }

                // Header
                c4.cd();
                TLatex h4; h4.SetNDC(); h4.SetTextAlign(22); h4.SetTextSize(0.036);
                h4.DrawLatex(0.52,0.975,"#Delta#phi  -  Mean / RMS vs E  (Constant-b (Production) & Variable-b (Benchmark): Uncorr vs Corr)");

                const std::string outP4 = dirPhiLog + "/DeltaPhi_MeanSigmaVsE_RawVsCorr.png";
                c4.SaveAs(outP4.c_str());
                c4.Close();

                // CSV for the four-series μ/σ(E) overlay
                std::vector<std::string> vlabCSV; vlabCSV.reserve(useV.size());
                std::vector<std::string> residCSV(useV.size(), "phi");
                std::vector<std::vector<double>> ECSV(useV.size(), eCtr);
                std::vector<std::vector<double>> MUC(useV.size());
                std::vector<std::vector<double>> DMUC(useV.size());
                std::vector<std::vector<double>> SGC(useV.size());
                std::vector<std::vector<double>> DSGC(useV.size());
                for (std::size_t i=0;i<useV.size();++i){
                    int v = useV[i];
                    vlabCSV.push_back(prettyLabel(v));
                    MUC[i]  = MU[v];   DMUC[i] = dMU[v];
                    SGC[i]  = SG[v];   DSGC[i] = dSG[v];
                }
                saveMuSigmaCSV(outP4, vlabCSV, residCSV, ECSV, MUC, DMUC, SGC, DSGC);
                
            }
        }
    }
    
    // =====================================================================================
    // NEW: CLUSEAvariantsVSRAW - CLUS-RAW vs CLUS-CP and vs each CLUS-CP(EA) variant
    //      • First-bin overlays, all-bins tables, μ/σ vs E summaries, and σ ratios
    //      • Plus: single μ/σ(E) overlays with ALL variants over RAW
    //      • Outputs go to:  outDir/CLUSEAvariantsVSRAW
    // =====================================================================================
    auto MakeClusEAvariantsVsRaw =
    [&](const std::vector<TH1F*>& Hphi_raw,
        const std::vector<TH1F*>& Hphi_cp,
        const std::vector<TH1F*>& Heta_raw,
        const std::vector<TH1F*>& Heta_cp)
    {
        const std::string dir = std::string(outDir) + "/CLUSEAvariantsVSRAW";
        const std::string dirFirstBin = dir + "/firstBinOverlays";
        const std::string dirMean     = dir + "/meanSummaries";
        ensureDir(dir.c_str());
        ensureDir(dirFirstBin.c_str());
        ensureDir(dirMean.c_str());

        const double eLo = eEdges.front().first;
        const double eHi = eEdges.front().second;

        // ---------------- Robust ratio maker: match by energy-bin index ----------------
        auto makeRatioSeriesByIndex =
        [&](const std::vector<TH1F*>& A, const std::vector<TH1F*>& B)->SigSeries
        {
            SigSeries R;
            if (A.empty() || B.empty()) return R;
            const std::size_t n = std::min(A.size(), B.size());
            R.E.reserve(n); R.S.reserve(n); R.dS.reserve(n);

            for (std::size_t i=0; i<n; ++i) {
                TH1F* hA = A[i];
                TH1F* hB = B[i];
                if (!hA || !hB) continue;
                if (hA->Integral()<=0 || hB->Integral()<=0) continue;

                double mA, meA, sA, seA;
                double mB, meB, sB, seB;
                std::tie(mA, meA, sA, seA) = statsFromHistRange(hA, xMin, xMax);
                std::tie(mB, meB, sB, seB) = statsFromHistRange(hB, xMin, xMax);
                if (sA<=0.0 || sB<=0.0) continue;

                const double Ectr = 0.5*(eEdges[i].first + eEdges[i].second);
                const double r    = sA / sB;
                const double dr   = r * std::sqrt( (seA>0 ? (seA/sA)*(seA/sA) : 0.0) +
                                                   (seB>0 ? (seB/sB)*(seB/sB) : 0.0) );
                R.E .push_back(Ectr);
                R.S .push_back(r);
                R.dS.push_back(dr);
            }
            return R;
        };

        // ---- small helpers (custom legend text) ----
        auto drawFirstBinOverlayCustom =
        [&](TH1F* hA, TH1F* hB,
            const char* title,
            const char* residLatex,
            const char* labelA,
            const char* labelB,
            Color_t col,
            const std::string& outPng)
        {
            if (!hA || !hB) return;
            if (hA->Integral()<=0 || hB->Integral()<=0) return;

            TCanvas c("cEA_FB","",1000,750);
            std::vector<ResidualSpec> specs = {
                { hA, Form("%s, %s", labelA, residLatex), col, (Style_t)20, 1 },
                { hB, Form("%s, %s", labelB, residLatex), col, (Style_t)24, 1 }
            };
            drawResidualPanel(&c, title, eLo, eHi, specs, xMin, xMax, 1.0);
            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        auto drawTableOverlayCustom =
        [&](const std::vector<TH1F*>& A,
            const std::vector<TH1F*>& B,
            const char* panelTitle,
            const char* residLatex,
            const char* labelA,
            const char* labelB,
            Color_t col,
            const std::string& outPng)
        {
            const int nCol=4, nRow=2, maxPads=nCol*nRow;
            TCanvas cT("cEA_Table","",1600,900);
            cT.SetTopMargin(0.10);
            cT.Divide(nCol,nRow,0,0);

            int pad=1;
            for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE)
            {
                if (!A[iE] || !B[iE] || A[iE]->Integral()==0 || B[iE]->Integral()==0) { ++pad; continue; }
                cT.cd(pad++);
                setPadMargins();

                std::vector<ResidualSpec> specs = {
                    { A[iE], Form("%s, %s", labelA, residLatex), col, (Style_t)20, 1 },
                    { B[iE], Form("%s, %s", labelB, residLatex), col, (Style_t)24, 1 }
                };
                drawResidualPanel(gPad, panelTitle,
                                  eEdges[iE].first, eEdges[iE].second,
                                  specs, xMin, xMax, 0.85);
            }
            cT.cd(0);
            drawHeaderNDC(Form("%s  -  %s", residLatex, panelTitle), 0.975, 0.045);
            cT.SaveAs(outPng.c_str());
            cT.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        auto doPair =
        [&](const std::vector<TH1F*>& RAW,
            const std::vector<TH1F*>& VAR,
            const char* residLatex,
            Color_t col,
            const char* labelRaw,
            const char* labelVar,
            const char* stem)
        {
            if (RAW.empty() || VAR.empty()) return;

            // first bin
            if (RAW[0] && VAR[0] && RAW[0]->Integral()>0 && VAR[0]->Integral()>0) {
                drawFirstBinOverlayCustom(RAW[0], VAR[0],
                    "Clusterizer: RAW vs variant",
                    residLatex, labelRaw, labelVar, col,
                    dirFirstBin + "/" + std::string(stem) + "_FirstBin.png");
            }

            // table (kept in the CLUSEAvariantsVSRAW root)
            drawTableOverlayCustom(RAW, VAR,
                "Clusterizer: RAW vs variant (all energy bins)",
                residLatex, labelRaw, labelVar, col,
                dir + "/" + std::string(stem) + "_AllBins_Table.png");

            // μ/σ vs E
            makeMuSigmaVariantOverlayWithLabels(RAW, VAR, col, 20, 24, "#mu",
                (dirMean + "/" + std::string(stem) + "_MeanSigmaVsE.png").c_str(),
                labelRaw, labelVar,
                Form("%s  -  %s vs %s  -  Mean / RMS vs E",
                     residLatex, labelRaw, labelVar));
        };

        // ---------------- CLUS-CP vs CLUS-RAW ----------------
        doPair(Hphi_raw, Hphi_cp,
               "#Delta#phi", kRed+1,
               "Constant-b (Production) Uncorr",
               "Constant-b (Production) Corr",
               "CLUSvsCP_DeltaPhi");
        doPair(Heta_raw, Heta_cp,
               "#Delta#eta", kBlue+1,
               "Constant-b (Production) Uncorr",
               "Constant-b (Production) Corr",
               "CLUSvsCP_DeltaEta");

        // ---------------- Each CLUS-CP(EA) variant vs CLUS-RAW ----------------
        struct EAdef { const char* key; const char* pretty; };
        const std::vector<EAdef> eaDefs = {
            {"geom",              "(Production) CP(|vz| + E)"},
            {"fitEtaDep",         "(Production) CP(EA |#eta|+E)"},
            {"fitEnergyOnly",     "(Production) CP(EA E-only)"},
            {"fitPhiE_etaEtaDep", "(Production) CP(EA #varphi:E-only, #eta:|#eta|+E)"}
        };

        std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>> phiRatioEAseries;
        std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>> etaRatioEAseries;

        // Also cache EA sets so we can make an "all variants over RAW" μ/σ overlay.
        struct SeriesSet { std::string label; std::vector<TH1F*> H; };
        std::vector<SeriesSet> phiEAsets, etaEAsets;

        auto pickColorForIdx = [](int i)->Color_t {
            switch(i){ case 0: return kRed+2; case 1: return kGreen+2; case 2: return kMagenta+1; default: return kOrange+7; }
        };
        auto pickMarkerForIdx = [](int i)->Style_t {
            switch(i){ case 0: return 24; case 1: return 25; case 2: return 27; default: return 28; }
        };

        for (std::size_t iv=0; iv<eaDefs.size(); ++iv)
        {
            const auto& d = eaDefs[iv];
            const std::string patP = std::string("h_phi_diff_cpCorrEA_") + d.key + "_%s";
            const std::string patE = std::string("h_eta_diff_cpCorrEA_") + d.key + "_%s";
            auto P = loadSet(patP.c_str());
            auto E = loadSet(patE.c_str());

            bool anyP=false; for (auto* h : P) if (h && h->Integral()>0) { anyP=true; break; }
            bool anyE=false; for (auto* h : E) if (h && h->Integral()>0) { anyE=true; break; }

            if (anyP) {
                doPair(Hphi_raw, P, "#Delta#phi", kRed+1,
                       "Constant-b (Production) Uncorr", d.pretty,
                       (std::string("CLUSEA_") + d.key + "_vs_RAW_DeltaPhi").c_str());

                // Try index-based ratio; if empty, fall back to σ(E) ratio
                SigSeries r = makeRatioSeriesByIndex(P, Hphi_raw);
                if (r.E.empty()) {
                    SigSeries sEA = buildSigmaSeries(P);
                    SigSeries sR  = buildSigmaSeries(Hphi_raw);
                    // classic centre-match ratio (already in your code base)
                    r = makeRatioSeries(sEA, sR);
                }
                if (!r.E.empty())
                    phiRatioEAseries.push_back({r, d.pretty, pickColorForIdx((int)iv), pickMarkerForIdx((int)iv)});
                phiEAsets.push_back({d.pretty, P});
            }

            if (anyE) {
                doPair(Heta_raw, E, "#Delta#eta", kBlue+1,
                       "Constant-b (Production) Uncorr", d.pretty,
                       (std::string("CLUSEA_") + d.key + "_vs_RAW_DeltaEta").c_str());

                SigSeries r = makeRatioSeriesByIndex(E, Heta_raw);
                if (r.E.empty()) {
                    SigSeries sEA = buildSigmaSeries(E);
                    SigSeries sR  = buildSigmaSeries(Heta_raw);
                    r = makeRatioSeries(sEA, sR);
                }
                if (!r.E.empty())
                    etaRatioEAseries.push_back({r, d.pretty, pickColorForIdx((int)iv), pickMarkerForIdx((int)iv)});
                etaEAsets.push_back({d.pretty, E});
            }
        }

        // Build CP/RAW ratios (to be added alongside all EA/RAW ratios) - with the same fallback
        SigSeries rCP_phi = makeRatioSeriesByIndex(Hphi_cp, Hphi_raw);
        if (rCP_phi.E.empty()) {
            SigSeries sCP = buildSigmaSeries(Hphi_cp);
            SigSeries sR  = buildSigmaSeries(Hphi_raw);
            rCP_phi       = makeRatioSeries(sCP, sR);
        }

        SigSeries rCP_eta = makeRatioSeriesByIndex(Heta_cp, Heta_raw);
        if (rCP_eta.E.empty()) {
            SigSeries sCP = buildSigmaSeries(Heta_cp);
            SigSeries sR  = buildSigmaSeries(Heta_raw);
            rCP_eta       = makeRatioSeries(sCP, sR);
        }
        // Helper: draw ALL ratios (CP and every EA) on ONE canvas (graphs kept alive)
        auto drawAllRatiosOneCanvas =
        [&](const std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>>& eaSeriesIn,
            const SigSeries& rCPIn,
            const char* lblCP,
            const char* yTitle,
            const std::string& outPngStem)
        {
            // Keep only non-empty series
            std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>> eaSeries;
            eaSeries.reserve(eaSeriesIn.size());
            for (const auto& t : eaSeriesIn) {
                const SigSeries& R = std::get<0>(t);
                if (!R.E.empty()) eaSeries.push_back(t);
            }
            const bool hasCP = !rCPIn.E.empty();
            if (!hasCP && eaSeries.empty()) {
                std::cerr << "[Playground][WARN] No ratio points to draw for " << outPngStem << "\n";
                return;
            }

            // X range from binning
            const double xMinE = 0.0;
            const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

            // Y range from available points; fallback → [0.85, 1.15]
            double yLo = +1e30, yHi = -1e30;
            auto scan = [&](const SigSeries& R){
                for (std::size_t i=0;i<R.S.size();++i){
                    const double err = (R.dS.empty()?0.0:R.dS[i]);
                    yLo = std::min(yLo, R.S[i] - err);
                    yHi = std::max(yHi, R.S[i] + err);
                }
            };
            if (hasCP) scan(rCPIn);
            for (const auto& t : eaSeries) scan(std::get<0>(t));

            if (!(yLo < yHi)) { yLo = 0.85; yHi = 1.15; }
            else {
                const double pad = 0.10*(yHi - yLo + 1e-9);
                yLo -= pad; yHi += pad;
            }

            TCanvas c("cEA_RatioAll","All ratios over RAW",900,600);
            c.SetLeftMargin(0.15);
            c.SetRightMargin(0.06);
            c.SetTopMargin(0.10);
            c.SetBottomMargin(0.18);

            TH1F fr("fr","",1,xMinE,xMaxE);
            fr.SetTitle("");
            fr.GetXaxis()->SetTitle("E  [GeV]");
            fr.GetYaxis()->SetTitle(yTitle);
            fr.GetXaxis()->SetTitleSize(0.050);
            fr.GetXaxis()->SetLabelSize(0.038);
            fr.GetXaxis()->SetTitleOffset(1.05);
            fr.GetYaxis()->SetTitleSize(0.050);
            fr.GetYaxis()->SetLabelSize(0.038);
            fr.SetMinimum(yLo); fr.SetMaximum(yHi);
            fr.Draw("AXIS");

            TLine l1(xMinE,1.0,xMaxE,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

            // Smaller legend text; two columns; slightly lower so it never overlaps
            TLegend lg(0.3,0.75,0.88,0.88);
            lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.025); lg.SetNColumns(2);

            std::vector<std::unique_ptr<TGraphErrors>> keep;

            // Draw one series as CLOSED circles (style=20), color-coded; points only (no lines)
            auto drawOne = [&](const SigSeries& R, Color_t col, const char* lab){
                if (R.E.empty()) return;
                std::vector<double> ex(R.E.size(), 0.0);
                auto g = std::make_unique<TGraphErrors>(
                             (int)R.E.size(),
                             const_cast<double*>(R.E.data()),
                             const_cast<double*>(R.S.data()),
                             ex.data(),
                             const_cast<double*>(R.dS.data()));
                g->SetMarkerStyle(20);        // closed circle for all
                g->SetMarkerSize(1.15);
                g->SetMarkerColor(col);
                g->SetLineColor(0);           // no line
                g->Draw("P SAME");            // points only
                lg.AddEntry(g.get(), lab, "p");
                keep.emplace_back(std::move(g));
            };

            if (hasCP) drawOne(rCPIn, kRed+1, lblCP);
            for (const auto& t : eaSeries) {
                const SigSeries& R = std::get<0>(t);
                const char*      L = std::get<1>(t);
                const Color_t    C = std::get<2>(t);
                drawOne(R, C, L);
            }

            lg.Draw();
            TLatex h; h.SetNDC(); h.SetTextAlign(13); h.SetTextSize(0.035);
            h.DrawLatex(0.16,0.955,"RMS Ratio: CP & EA variants over RAW");

            c.SaveAs((dir + "/" + outPngStem).c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << dir + "/" + outPngStem << "\n";
        };

        // One canvas for Δφ
        drawAllRatiosOneCanvas(
            phiRatioEAseries,
            rCP_phi,
            "(Production) Corr / RAW",
            "Ratio  #sigma_{RMS}(#Delta#phi) / #sigma_{RMS}(#Delta#phi)_{RAW}",
            "Ratio_AllVariants_over_RAW_DeltaPhi.png"
        );

        // One canvas for Δη
        drawAllRatiosOneCanvas(
            etaRatioEAseries,
            rCP_eta,
            "(Production) Corr / RAW",
            "Ratio  #sigma_{RMS}(#Delta#eta) / #sigma_{RMS}(#Delta#eta)_{RAW}",
            "Ratio_AllVariants_over_RAW_DeltaEta.png"
        );

        // ---- Extra: emit a version that EXCLUDES the CP(|vz| + E) variant -----------------
        auto filterOutGeom =
            [&](const std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>>& in)
                -> std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>>
        {
            std::vector<std::tuple<SigSeries,const char*,Color_t,Style_t>> out;
            out.reserve(in.size());
            for (const auto& t : in) {
                const char* lbl = std::get<1>(t);
                // skip any label that contains "geom" (case-insensitive)
                TString L(lbl); L.ToLower();
                if (L.Contains("geom")) continue;
                out.push_back(t);
            }
            return out;
        };

        auto phiNoGeom = filterOutGeom(phiRatioEAseries);
        auto etaNoGeom = filterOutGeom(etaRatioEAseries);

        drawAllRatiosOneCanvas(
            phiNoGeom,
            rCP_phi,
            "Constant-b (Production) Corr / RAW",
            "Ratio  #sigma_{RMS}(#Delta#phi) / #sigma_{RMS}(#Delta#phi)_{RAW}",
            "Ratio_AllVariants_over_RAW_DeltaPhi_noGeom.png"
        );

        drawAllRatiosOneCanvas(
            etaNoGeom,
            rCP_eta,
            "Constant-b (Production) Corr / RAW",
            "Ratio  #sigma_{RMS}(#Delta#eta) / #sigma_{RMS}(#Delta#eta)_{RAW}",
            "Ratio_AllVariants_over_RAW_DeltaEta_noGeom.png"
        );


        // ---------------- μ/σ(E) overlays with ALL variants over RAW (single canvas per residual) ---------
        auto drawAllVariantsMuSigmaOverRaw =
        [&](const std::vector<TH1F*>& RAW,
            const std::vector<TH1F*>& CP,
            const std::vector<SeriesSet>& EAsets,
            const char* residLatex,
            Color_t baseColor,
            const std::string& outPng)
        {
            // Defensive: ensure the parent directory of outPng exists
            auto ensureParentDir = [&](const std::string& p){
                std::string d = p;
                std::size_t pos = d.find_last_of("/\\");
                if (pos != std::string::npos) {
                    d.resize(pos);
                    gSystem->mkdir(d.c_str(), /*recursive*/true);
                }
            };
            ensureParentDir(outPng);

            // Build μ/σ arrays for RAW and CP
            std::vector<double> eCtr; eCtr.reserve(eEdges.size());
            std::vector<double> muR, dmuR, sgR, dsgR;
            std::vector<double> muC, dmuC, sgC, dsgC;

            auto pushIf = [&](TH1F* h, std::vector<double>& mu, std::vector<double>& dmu,
                              std::vector<double>& sg, std::vector<double>& dsg, std::size_t i){
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(h, xMin, xMax);
                mu.push_back(m);  dmu.push_back(me);
                sg.push_back(r);  dsg.push_back(re);
                if (eCtr.size() < i+1) eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second));
            };

            const std::size_t N = eEdges.size();
            for (std::size_t i=0;i<N;++i){
                if (RAW[i] && RAW[i]->Integral()>0){
                    pushIf(RAW[i], muR, dmuR, sgR, dsgR, i);
                } else {
                    // still advance eCtr to keep axis scaling sane if anything else exists in this bin
                    bool hasEA = false;
                    for (const auto& s : EAsets) { if (i<s.H.size() && s.H[i] && s.H[i]->Integral()>0){ hasEA=true; break; } }
                    if (hasEA || (CP[i] && CP[i]->Integral()>0)){
                        if (eCtr.size() < i+1) eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second));
                    }
                }
                if (CP[i] && CP[i]->Integral()>0){
                    pushIf(CP[i],  muC, dmuC, sgC, dsgC, i);
                } else {
                    // keep vector sizes aligned with RAW positions
                    if (!muR.empty()) { muC.push_back(0); dmuC.push_back(0); sgC.push_back(0); dsgC.push_back(0); }
                }
            }

            // Pre-compute EA series μ/σ
            struct SeriesMU { std::string label; std::vector<double> mu, dmu, sg, dsg; };
            std::vector<SeriesMU> eaMU; eaMU.reserve(EAsets.size());
            for (const auto& s : EAsets){
                SeriesMU S; S.label = s.label;
                S.mu.reserve(N); S.dmu.reserve(N); S.sg.reserve(N); S.dsg.reserve(N);
                for (std::size_t i=0;i<N;++i){
                    if (!s.H[i] || s.H[i]->Integral()<=0) { S.mu.push_back(0); S.dmu.push_back(0); S.sg.push_back(0); S.dsg.push_back(0); continue; }
                    double m, me, r, re; std::tie(m, me, r, re) = statsFromHistRange(s.H[i], xMin, xMax);
                    S.mu.push_back(m); S.dmu.push_back(me); S.sg.push_back(r); S.dsg.push_back(re);
                    if (eCtr.size() < i+1) eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second)); // make sure axis has this bin
                }
                eaMU.push_back(std::move(S));
            }

            // If we still have zero centers, synthesize a frame over the full energy range
            if (eCtr.empty()) {
                for (std::size_t i=0;i<N;++i) eCtr.push_back(0.5*(eEdges[i].first + eEdges[i].second));
            }

            const double xAxisMin = 0.0;
            const double xAxisMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
            std::vector<double> ex(eCtr.size(), 0.0);

            TCanvas c("cAllVar_MuSig","All variants over RAW - Mean/RMS",1000,850);
            auto pads = makeEqualHeightTwinPads(c, "AllVarMuSig", 0.15,0.05, 0.14,0.02, 0.04,0.16);
            TPad* pTop = pads.first; TPad* pBot = pads.second;

            // ---------------- μ(E)
            pTop->cd();
            double muLo=+1e30, muHi=-1e30;
            auto updMuRange = [&](double val, double err){ muLo = std::min(muLo, val-err); muHi = std::max(muHi, val+err); };

            for (std::size_t i=0;i<muR.size();++i){ updMuRange(muR[i], dmuR[i]); }
            for (std::size_t i=0;i<muC.size();++i){ updMuRange(muC[i], dmuC[i]); }
            for (const auto& S : eaMU)
                for (std::size_t i=0;i<S.mu.size();++i) updMuRange(S.mu[i], S.dmu[i]);

            if (!(muLo<muHi)) { muLo=-1e-4; muHi=+1e-4; } // robust fallback
            const double padMu = 0.25*(muHi-muLo);
            TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
            frU.SetMinimum(muLo - padMu);
            frU.SetMaximum(muHi + padMu);
            frU.GetXaxis()->SetLabelSize(0.0);
            frU.GetXaxis()->SetTitleSize(0.0);
            frU.GetXaxis()->SetTickLength(0.0);
            frU.Draw("AXIS");
            TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

            // Keep graphs alive
            std::vector<std::unique_ptr<TGraphErrors>> keepMu, keepSg;

            // RAW and CP (baseline color) - draw if we have any points
            if (!muR.empty()){
                auto gMuR = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), muR.data(), ex.data(), dmuR.data());
                gMuR->SetMarkerStyle(20); gMuR->SetMarkerColor(baseColor); gMuR->SetLineColor(baseColor); gMuR->SetMarkerSize(1.0); gMuR->Draw("P SAME");
                keepMu.emplace_back(std::move(gMuR));
            }
            if (!muC.empty()){
                auto gMuC = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), muC.data(), ex.data(), dmuC.data());
                gMuC->SetMarkerStyle(24); gMuC->SetMarkerColor(baseColor); gMuC->SetLineColor(baseColor); gMuC->SetMarkerSize(1.0); gMuC->Draw("P SAME");
                keepMu.emplace_back(std::move(gMuC));
            }

            // EA variants
            for (std::size_t iv=0; iv<eaMU.size(); ++iv){
                auto& S = eaMU[iv];
                auto g = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), S.mu.data(), ex.data(), S.dmu.data());
                g->SetMarkerStyle(pickMarkerForIdx((int)iv));
                g->SetMarkerColor(pickColorForIdx((int)iv));
                g->SetLineColor  (pickColorForIdx((int)iv));
                g->SetMarkerSize(1.0);
                g->Draw("P SAME");
                keepMu.emplace_back(std::move(g));
            }

            // Legend
            TLegend legU(0.54,0.70,0.93,0.92);
            legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.033); legU.SetNColumns(2);
            if (!keepMu.empty()) {
                if (keepMu.size()>=1) legU.AddEntry(keepMu[0].get(), "Constant-b (Production) Uncorr", "p");
                if (keepMu.size()>=2) legU.AddEntry(keepMu[1].get(), "Constant-b (Production) Corr",   "p");
            }
            for (std::size_t iv=0; iv<eaMU.size(); ++iv)
                legU.AddEntry(keepMu[ (keepMu.size()>2? 2+iv : iv) ].get(), eaMU[iv].label.c_str(), "p");
            legU.Draw();

            // ---------------- σ(E)
            pBot->cd();
            double sgLo=+1e30, sgHi=-1e30;
            auto updSgRange = [&](double val, double err){ sgLo = std::min(sgLo, val-err); sgHi = std::max(sgHi, val+err); };

            for (std::size_t i=0;i<sgR.size();++i){ updSgRange(sgR[i], dsgR[i]); }
            for (std::size_t i=0;i<sgC.size();++i){ updSgRange(sgC[i], dsgC[i]); }
            for (const auto& S : eaMU)
                for (std::size_t i=0;i<S.sg.size();++i) updSgRange(S.sg[i], S.dsg[i]);

            if (!(sgLo<sgHi)) { sgLo=0.0; sgHi=1.0; }
            TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
            frL.SetMinimum(std::max(0.0, sgLo - 0.15*(sgHi-sgLo)));
            frL.SetMaximum(sgHi + 0.15*(sgHi-sgLo));
            frL.Draw("AXIS");

            if (!sgR.empty()){
                auto gSgR = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), sgR.data(), ex.data(), dsgR.data());
                gSgR->SetMarkerStyle(20); gSgR->SetMarkerColor(baseColor); gSgR->SetLineColor(baseColor); gSgR->SetMarkerSize(1.0); gSgR->Draw("P SAME");
                keepSg.emplace_back(std::move(gSgR));
            }
            if (!sgC.empty()){
                auto gSgC = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), sgC.data(), ex.data(), dsgC.data());
                gSgC->SetMarkerStyle(24); gSgC->SetMarkerColor(baseColor); gSgC->SetLineColor(baseColor); gSgC->SetMarkerSize(1.0); gSgC->Draw("P SAME");
                keepSg.emplace_back(std::move(gSgC));
            }
            for (std::size_t iv=0; iv<eaMU.size(); ++iv){
                auto& S = eaMU[iv];
                auto g = std::make_unique<TGraphErrors>((int)eCtr.size(), eCtr.data(), S.sg.data(), ex.data(), S.dsg.data());
                g->SetMarkerStyle(pickMarkerForIdx((int)iv));
                g->SetMarkerColor(pickColorForIdx((int)iv));
                g->SetLineColor  (pickColorForIdx((int)iv));
                g->SetMarkerSize(1.0);
                g->Draw("P SAME");
                keepSg.emplace_back(std::move(g));
            }

            c.cd();
            TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.038);
            h.DrawLatex(0.50,0.985,Form("%s  -  RAW (filled) vs all variants (open/custom)  -  Mean / RMS vs E", residLatex));

            // Always write the file, even if it ends up empty - avoids “missing file” surprises
            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // Emit the “all variants over RAW” μ/σ overlays
        drawAllVariantsMuSigmaOverRaw(Hphi_raw, Hphi_cp, phiEAsets,
            "#Delta#phi", kRed+1,  dirMean + "/AllVariantsOverRAW_DeltaPhi_MeanSigmaVsE.png");
        drawAllVariantsMuSigmaOverRaw(Heta_raw, Heta_cp, etaEAsets,
            "#Delta#eta", kBlue+1, dirMean + "/AllVariantsOverRAW_DeltaEta_MeanSigmaVsE.png");
    };  // <<< CLOSES the MakeClusEAvariantsVsRaw lambda
    // Invoke the new block (uses the CLUS raw/corrected sets we already loaded)
    MakeClusEAvariantsVsRaw(Hphi, HphiCorr, Heta, HetaCorr);

    // =====================================================================================
    // NEW: Emit a curated bundle of plots under /Users/patsfan753/Desktop/FinalPDCPlots
    //       • For φ and η in parallel
    //       • Mean/σ vs E overlays for (EA|η|+E, EA E-only, CP, RAW)
    //       • Mean/σ vs E overlays for (EA|η|+E, CP, RAW)  and  (EA E-only, CP, RAW)
    //       • First-bin overlays (counts) and 4×2 table overlays for the same bundles
    //       • Ratio plots: (each variant / RAW); closed circles only; no connecting lines
    //       • NEW: σ(E) overlay (top) + σ/σ_RAW ratio (bottom) for each 3-variant bundle
    // =====================================================================================
    {
      const std::string finalBase = "/Users/patsfan753/Desktop/FinalPDCPlots";
      auto Ensure = [&](const std::string& d){ gSystem->mkdir(d.c_str(), /*recursive*/true); };

      const std::string dPhi = finalBase + "/phi";
      const std::string dEta = finalBase + "/eta";
      Ensure(finalBase); Ensure(dPhi); Ensure(dEta);
      Ensure(dPhi + "/meansigma"); Ensure(dEta + "/meansigma");
      Ensure(dPhi + "/firstbin");  Ensure(dEta + "/firstbin");
      Ensure(dPhi + "/tables");    Ensure(dEta + "/tables");
      Ensure(dPhi + "/ratios");    Ensure(dEta + "/ratios");
      Ensure(dPhi + "/sigmaRatio");Ensure(dEta + "/sigmaRatio"); // NEW

      // Helper: draw Mean/σ vs E for N variants (N ≥ 2) on twin pads
      auto drawMuSigmaMulti = [&](const std::string& outPng,
                                  const char* header,
                                  const std::vector<std::pair<std::string,std::vector<TH1F*>>>& V,
                                  Color_t /*baseColorIf2*/ = kRed+1)
      {
        if (V.size() < 2) return;

        // Build μ/σ arrays
        std::vector<double> eCtr; eCtr.reserve(eEdges.size());
        std::vector<std::vector<double>> MU(V.size()), dMUv(V.size()), SG(V.size()), dSGv(V.size());

        for (std::size_t i=0;i<eEdges.size();++i)
        {
          bool rowHasAny=false;
          for (std::size_t k=0;k<V.size();++k)
          {
            TH1F* h = (i < V[k].second.size()) ? V[k].second[i] : nullptr;
            if (h && h->Integral()>0) rowHasAny=true;
          }
          if (!rowHasAny) continue;

          const double eC = 0.5*(eEdges[i].first + eEdges[i].second);
          eCtr.push_back(eC);

          for (std::size_t k=0;k<V.size();++k)
          {
            TH1F* h = (i < V[k].second.size()) ? V[k].second[i] : nullptr;
              double m=0,me=0,r=0,re=0;
              if (h && h->Integral()>0) std::tie(m,me,r,re) = statsFromHistRange(h, xMin, xMax);

              // Fallbacks to ensure mean/stddev errors are populated if statsFromHistRange returns 0
              if (h && me <= 0) {
                double me_root = h->GetMeanError();          // ROOT’s mean error
                if (me_root > 0) me = me_root;
              }
              if (h && re <= 0) {
                double re_sd    = h->GetStdDevError();       // error of std dev (preferred)
                double re_rms   = h->GetRMSError();          // error of RMS (alt)
                if (re_sd  > 0) re = re_sd;
                else if (re_rms > 0) re = re_rms;
              }

              MU [k].push_back(m);
              dMUv[k].push_back(me);
              SG [k].push_back(r);
              dSGv[k].push_back(re);
          }
        }
        if (eCtr.empty()) return;

        const double xMinE = 0.0;
        const double xMaxE = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
        std::vector<double> ex(eCtr.size(), 0.0);

        // style palette for multiple variants
        const Color_t cols[] = { kGreen+2, kMagenta+1, kRed+1, kBlue+1, kOrange+7, kCyan+2, kGray+2 };
        const Style_t mks [] = { 20, 20, 20, 20, 20, 20, 20 }; // closed circles everywhere
        const int NC = sizeof(cols)/sizeof(cols[0]);

        TCanvas c("cMuSigMulti","μ/σ overlays",1000,850);
        auto pads = makeEqualHeightTwinPads(c, "MuSigMulti", 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first; TPad* pBot = pads.second;

        // ---- μ(E)
        pTop->cd();
        double muLo=+1e30, muHi=-1e30;
        for (std::size_t k=0;k<V.size();++k)
          for (std::size_t i=0;i<MU[k].size();++i)
          { muLo = std::min(muLo, MU[k][i]-dMUv[k][i]); muHi = std::max(muHi, MU[k][i]+dMUv[k][i]); }
        if (!(muLo<muHi)) { muLo=-1e-4; muHi=+1e-4; }
        const double padMu = 0.25*(muHi-muLo);

        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(muLo - padMu); frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xMinE,0.0,xMaxE,0.0); l0.SetLineStyle(2); l0.Draw();

          // Draw μ(E) graphs first (store pointers for legend entries)
          std::vector<std::unique_ptr<TGraphErrors>> keepMu, keepSg;
          for (std::size_t k=0;k<V.size();++k)
          {
            auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                    eCtr.data(), MU[k].data(),
                                                    ex.data(),   dMUv[k].data());
              g->SetMarkerStyle(mks[k%NC]);       // closed circles
              g->SetMarkerColor(cols[k%NC]);
              g->SetLineColor(cols[k%NC]);        // make error bar same color as marker (no white slit)
              g->SetLineWidth(1);                 // thin bar
              g->SetMarkerSize(1.0);
              g->Draw("P SAME");                  // points only (no connecting line)
            keepMu.emplace_back(std::move(g));
          }

          // Deterministic two-column legends for 3-variant overlays (heap-allocated legends).
          // Expected order of V: [0]=EA variant, [1]=CP ("Constant-b (legacy)"), [2]=RAW ("No Correction")
          pTop->cd();

          if (V.size() == 3)
          {
            const int idxEA  = 0;
            const int idxCP  = 1;
            const int idxRAW = 2;

              // Put both legends on the same Y span so first rows align,
              // and bring them closer horizontally.
              const double y1 = 0.73, y2 = 0.83;   // common vertical band for both columns
              const double xL1 = 0.18, xL2 = 0.24; // left column (RAW only)
              const double xR1 = 0.32, xR2 = 0.38; // right column (CP on top, EA below)

              // Left legend: RAW only
              TLegend* legLeft = new TLegend(xL1, y1, xL2, y2, "", "brNDC");
              legLeft->SetBorderSize(0);
              legLeft->SetFillStyle(0);
              legLeft->SetTextSize(0.04);
              legLeft->SetEntrySeparation(0.0);
              legLeft->SetMargin(0.20);
              if ((int)keepMu.size() > idxRAW && keepMu[idxRAW])
                legLeft->AddEntry(keepMu[idxRAW].get(), V[idxRAW].first.c_str(), "p");
              legLeft->Draw();

              // Right legend: CP on top, then EA
              TLegend* legRight = new TLegend(xR1, y1, xR2, y2, "", "brNDC");
              legRight->SetBorderSize(0);
              legRight->SetFillStyle(0);
              legRight->SetTextSize(0.04);
              legRight->SetEntrySeparation(0.0);
              legRight->SetMargin(0.20);
              if ((int)keepMu.size() > idxCP  && keepMu[idxCP])
                legRight->AddEntry(keepMu[idxCP].get(),  V[idxCP].first.c_str(),  "p");
              if ((int)keepMu.size() > idxEA  && keepMu[idxEA])
                legRight->AddEntry(keepMu[idxEA].get(),  V[idxEA].first.c_str(),  "p");
              legRight->Draw();

          }
          else
          {
            // Fallback: single legend with all entries (also heap-allocated)
            TLegend* legU = new TLegend(0.16, 0.73, 0.92, 0.92, "", "brNDC");
            legU->SetBorderSize(0);
            legU->SetFillStyle(0);
            legU->SetTextSize(0.033);
            legU->SetNColumns(1);
            for (std::size_t k = 0; k < V.size(); ++k)
              legU->AddEntry(keepMu[k].get(), V[k].first.c_str(), "p");
            legU->Draw();
          }

          // Force the top pad to register the newly drawn legends before switching pads
          pTop->Modified();
          pTop->Update();


        // ---- σ(E)
        pBot->cd();
        double sgHi=-1e30;
        for (std::size_t k=0;k<V.size();++k)
          for (double s : SG[k]) sgHi = std::max(sgHi, s);
        if (!(sgHi>0)) sgHi = 1.0;

        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

        for (std::size_t k=0;k<V.size();++k)
        {
          auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                  eCtr.data(), SG[k].data(),
                                                  ex.data(),   dSGv[k].data());
            g->SetMarkerStyle(mks[k%NC]);
            g->SetMarkerColor(cols[k%NC]);
            g->SetLineColor(cols[k%NC]);
            g->SetLineWidth(1);
            g->SetMarkerSize(1.0);
            g->Draw("P SAME");
          keepSg.emplace_back(std::move(g));
        }

          // header (smaller font only for SigmaRatio_EA_vs_CP_* outputs)
          c.cd();
          TLatex h; h.SetNDC(); h.SetTextAlign(22);
          {
            TString outName(outPng.c_str());
            const bool isEAvsCP = outName.Contains("SigmaRatio_EA_vs_CP_");
            h.SetTextSize(isEAvsCP ? 0.028 : 0.036);
          }
          h.DrawLatex(0.52,0.975, header);

          // Single save/close (avoid double SaveAs/Close that can segfault)
          c.SaveAs(outPng.c_str());
          c.Close();
      };

      // Helper: first-bin multi-overlay (counts) for N histograms
      auto drawFirstBinMulti = [&](const std::string& outPng,
                                   const char* panelTitle,
                                   const std::vector<std::pair<std::string,TH1F*>>& items,
                                   double eLo, double eHi,
                                   bool isPhi)
      {
        if (items.empty()) return;
        TCanvas c("cFBmulti","first-bin overlay",1000,750);
        std::vector<ResidualSpec> specs; specs.reserve(items.size());
        const Color_t cols[] = { kGreen+2, kMagenta+1, kRed+1, kBlue+1, kOrange+7, kCyan+2 };
        for (std::size_t i=0;i<items.size();++i)
          specs.push_back( { items[i].second,
                             (isPhi?"#Delta#phi ":"#Delta#eta ") + TString(items[i].first.c_str()),
                             cols[i%6], (Style_t)20, 1 } );
        drawResidualPanel(&c, panelTitle, eLo, eHi, specs, xMin, xMax, 1.0);
        c.SaveAs(outPng.c_str());
        c.Close();
      };

      // Helper: 4×2 table overlay for N variants
      auto drawTableMulti = [&](const std::string& outPng,
                                const char* panelTitle,
                                const std::vector<std::pair<std::string,std::vector<TH1F*>>>& sets,
                                bool isPhi)
      {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cTblMulti","table overlay",1600,900);
        cT.SetTopMargin(0.10);
        cT.Divide(nCol,nRow,0,0);

        int pad=1;
        for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE)
        {
          // Gather non-empty hists for this bin
          std::vector<ResidualSpec> specs;
          const Color_t cols[] = { kGreen+2, kMagenta+1, kRed+1, kBlue+1, kOrange+7, kCyan+2 };
          for (std::size_t k=0;k<sets.size();++k)
          {
            TH1F* h = (iE<sets[k].second.size()) ? sets[k].second[iE] : nullptr;
            if (!h || h->Integral()==0) continue;
            specs.push_back( { h,
                               (isPhi?"#Delta#phi ":"#Delta#eta ") + TString(sets[k].first.c_str()),
                               cols[k%6], (Style_t)20, 1 } );
          }
          if (specs.empty()) { ++pad; continue; }

          cT.cd(pad++);
          setPadMargins();
          drawResidualPanel(gPad, panelTitle, eEdges[iE].first, eEdges[iE].second, specs, xMin, xMax, 0.85);
        }

        cT.cd(0);
        drawHeaderNDC(panelTitle, 0.975, 0.045);
        cT.SaveAs(outPng.c_str());
        cT.Close();
      };

      // Helper: ratio overlay (each variant / RAW) - closed circles only, small legend, no lines
      auto ratioOverlay = [&](const std::string& outPng,
                              const char* yTitle,
                              const std::vector<std::pair<std::string,std::vector<TH1F*>>>& variants,
                              const std::vector<TH1F*>& RAW)
      {
        // Make by-index ratios
        auto makeRatioSeriesByIndex = [&](const std::vector<TH1F*>& A, const std::vector<TH1F*>& B)->SigSeries
        {
          SigSeries R;
          const std::size_t n = std::min(A.size(), B.size());
          for (std::size_t i=0;i<n;++i)
          {
            TH1F* hA=A[i], *hB=B[i];
            if (!hA || !hB) continue;
            if (hA->Integral()<=0 || hB->Integral()<=0) continue;
            double mA,meA,sA,seA, mB,meB,sB,seB;
            std::tie(mA,meA,sA,seA) = statsFromHistRange(hA, xMin, xMax);
            std::tie(mB,meB,sB,seB) = statsFromHistRange(hB, xMin, xMax);
            if (!(sA>0 && sB>0)) continue;
            const double Ectr = 0.5*(eEdges[i].first + eEdges[i].second);
            const double r = sA/sB;
            const double dr = r * std::sqrt( (seA>0? (seA/sA)*(seA/sA) : 0.0) +
                                             (seB>0? (seB/sB)*(seB/sB) : 0.0) );
            R.E.push_back(Ectr); R.S.push_back(r); R.dS.push_back(dr);
          }
          return R;
        };

        // Graph frame
        const double xMinE = 0.0;
        const double xMaxE = eEdges.empty() ? 1.0 : eEdges.back().second;

        // Collect series
        struct RS { SigSeries s; std::string lab; Color_t col; };
        const Color_t cols[] = { kGreen+2, kMagenta+1, kRed+1 };
        std::vector<RS> rows; rows.reserve(variants.size());
        for (std::size_t k=0;k<variants.size();++k)
        {
          RS rs; rs.lab = variants[k].first; rs.col = cols[k%3];
          rs.s = makeRatioSeriesByIndex(variants[k].second, RAW);
          rows.push_back(std::move(rs));
        }

        // y-range
        double yLo=+1e30, yHi=-1e30;
        for (const auto& r : rows)
          for (std::size_t i=0;i<r.s.S.size();++i)
          { const double e = (r.s.dS.empty()?0.0:r.s.dS[i]);
            yLo = std::min(yLo, r.s.S[i]-e);
            yHi = std::max(yHi, r.s.S[i]+e); }
        if (!(yLo<yHi)) { yLo=0.85; yHi=1.15; } else { const double pad=0.10*(yHi-yLo+1e-9); yLo-=pad; yHi+=pad; }

        TCanvas c("cRatio","ratios",900,600);
        c.SetLeftMargin(0.15); c.SetRightMargin(0.06); c.SetTopMargin(0.10); c.SetBottomMargin(0.18);

        TH1F fr("fr","",1,xMinE,xMaxE);
        fr.SetTitle("");
        fr.GetXaxis()->SetTitle("E  [GeV]");
        fr.GetYaxis()->SetTitle(yTitle);
        fr.GetXaxis()->SetTitleSize(0.050);
        fr.GetXaxis()->SetLabelSize(0.038);
        fr.GetXaxis()->SetTitleOffset(1.05);
        fr.GetYaxis()->SetTitleSize(0.050);
        fr.GetYaxis()->SetLabelSize(0.038);
        fr.SetMinimum(yLo); fr.SetMaximum(yHi);
        fr.Draw("AXIS");

        TLine l1(xMinE,1.0,xMaxE,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

        TLegend lg(0.30,0.75,0.88,0.88);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.025); lg.SetNColumns(2);

        std::vector<std::unique_ptr<TGraphErrors>> keep;
        for (const auto& r : rows)
        {
          const SigSeries& R = r.s;
          if (R.E.empty()) continue;
          std::vector<double> ex(R.E.size(), 0.0);
          auto g = std::make_unique<TGraphErrors>((int)R.E.size(),
                            const_cast<double*>(R.E.data()),
                            const_cast<double*>(R.S.data()),
                            ex.data(),
                            const_cast<double*>(R.dS.data()));
          g->SetMarkerStyle(20);   // closed circle only
          g->SetMarkerSize(1.15);
          g->SetMarkerColor(r.col);
          g->SetLineColor(0);      // no connecting line
          g->Draw("P SAME");
          lg.AddEntry(g.get(), r.lab.c_str(), "p");
          keep.emplace_back(std::move(g));
        }

        lg.Draw();
        c.SaveAs(outPng.c_str());
        c.Close();
      };

        // NEW: σ(E) top + ratio bottom for N variants (N ≥ 2).
        // Denominator is provided by `denomHists` and named by `denomLabelToken` ("RAW", "CP", ...).
        auto drawSigmaWithRatioMulti =
          [&](const std::string& outPng,
              const char* header,
              const std::vector<std::pair<std::string,std::vector<TH1F*>>>& V,
              const std::vector<TH1F*>& denomHists,
              const char* denomLabelToken,
              bool isPhi)
        {
          if (V.size() < 2) return;

          auto toLower = [](std::string s){ for (auto& c : s) c = std::tolower(c); return s; };

          // Build σ arrays for all variants and denominator
          std::vector<double> eCtr; eCtr.reserve(eEdges.size());
          std::vector<std::vector<double>> SG(V.size()), dSG(V.size());
          std::vector<double> Sden, dSden;

          for (std::size_t i=0;i<eEdges.size();++i)
          {
            bool rowHasAny=false;
            for (std::size_t k=0;k<V.size();++k)
            {
              TH1F* h = (i < V[k].second.size()) ? V[k].second[i] : nullptr;
              if (h && h->Integral()>0) { rowHasAny = true; break; }
            }
            if (!rowHasAny) continue;

            const double eC = 0.5*(eEdges[i].first + eEdges[i].second);
            eCtr.push_back(eC);

              for (std::size_t k=0;k<V.size();++k)
              {
                TH1F* h = (i < V[k].second.size()) ? V[k].second[i] : nullptr;
                double m=0,me=0,s=0,se=0;
                if (h && h->Integral()>0) std::tie(m,me,s,se) = statsFromHistRange(h, xMin, xMax);

                // Fallbacks to guarantee non-zero σ errors when available from ROOT
                if (h && se <= 0) {
                  const double se_sd  = h->GetStdDevError();  // preferred
                  const double se_rms = h->GetRMSError();     // alternative
                  if (se_sd  > 0) se = se_sd;
                  else if (se_rms > 0) se = se_rms;
                }

                SG .at(k).push_back(s);
                dSG.at(k).push_back(se);
              }

              {
                TH1F* hD = (i < denomHists.size()) ? denomHists[i] : nullptr;
                double m=0,me=0,s=0,se=0;
                if (hD && hD->Integral()>0) std::tie(m,me,s,se) = statsFromHistRange(hD, xMin, xMax);

                // Same fallback for denominator errors
                if (hD && se <= 0) {
                  const double se_sd  = hD->GetStdDevError();
                  const double se_rms = hD->GetRMSError();
                  if (se_sd  > 0) se = se_sd;
                  else if (se_rms > 0) se = se_rms;
                }

                Sden.push_back(s);
                dSden.push_back(se);
              }

          }
          if (eCtr.empty()) return;

          const double xMinE = 0.0;
          const double xMaxE = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
          std::vector<double> ex(eCtr.size(), 0.0);

          const Color_t cols[] = { kGreen+2, kMagenta+1, kRed+1, kBlue+1, kOrange+7, kCyan+2, kGray+2 };
          const Style_t mk   = 20; // closed circle
          const int NC = sizeof(cols)/sizeof(cols[0]);

          TCanvas c("cSigRat","σ & ratio",1000,850);
          auto pads = makeEqualHeightTwinPads(c, "SigRat", 0.15,0.05, 0.14,0.02, 0.04,0.16);
          TPad* pTop = pads.first; TPad* pBot = pads.second;

            pTop->cd();
            // Tight auto-zoom using data ± error
            // Tight auto-zoom using data ± error
            double sLo=+1e30, sHi=-1e30;
            for (std::size_t k=0;k<V.size();++k)
              for (std::size_t i=0;i<SG[k].size();++i) {
                const double s  = SG[k][i];
                const double ds = (i<dSG[k].size()? dSG[k][i] : 0.0);
                sLo = std::min(sLo, s - ds);
                sHi = std::max(sHi, s + ds);
              }
            if (!(sLo < sHi)) { sLo = 0.0; sHi = 1.0; }

            // Compute data range
            const double range = (sHi - sLo + 1e-12);

            // Base padding
            double padLo = 0.07 * range;
            double padHi = 0.07 * range;

            // For Δφ give more headroom by default
            if (isPhi) padHi = 0.18 * range;   // was 0.07; ~2.6× more space

            // File-specific override: SigmaRatio_EAeOnly_CP_RAW_phi needs extra room for the legend
            {
              TString outName(outPng.c_str());
              if (isPhi && outName.Contains("SigmaRatio_EAeOnly_CP_RAW_phi")) {
                // Ensure at least 8×10^-4 rad of absolute headroom
                padHi = std::max(padHi, 0.00080);
              }
            }

            // Freeze the axis with the enlarged headroom
            TH1F frU("frU",";E  [GeV];#sigma_{RMS}  [rad]",1,xMinE,xMaxE);
            frU.SetMinimum(sLo - padLo);
            frU.SetMaximum(sHi + padHi);
            frU.GetXaxis()->SetLabelSize(0.0);
            frU.GetXaxis()->SetTitleSize(0.0);
            frU.GetXaxis()->SetTickLength(0.0);
            frU.Draw("AXIS");


            std::vector<std::unique_ptr<TGraphErrors>> keepTop, keepBot;
            for (std::size_t k=0;k<V.size();++k)
            {
              auto g = std::make_unique<TGraphErrors>((int)eCtr.size(),
                                                      eCtr.data(), SG[k].data(),
                                                      ex.data(),   dSG[k].data());
              g->SetMarkerStyle(mk);
              g->SetMarkerColor(cols[k%NC]);
              g->SetLineColor(cols[k%NC]);  // color error bars same as marker
              g->SetLineWidth(1);
              g->SetMarkerSize(1.0);
              g->Draw("P SAME");
              keepTop.emplace_back(std::move(g));
            }

            // Split legend (left column = denominator only) if denom is also in V and V has 3 entries.
            int denIdxInV = -1;
            {
              const std::string token = toLower(denomLabelToken ? denomLabelToken : "");
              for (std::size_t k=0;k<V.size();++k)
              {
                const std::string lab = toLower(V[k].first);
                const bool isRAW = (token.find("raw")!=std::string::npos) &&
                                   (lab.find("raw")!=std::string::npos || lab.find("no correction")!=std::string::npos);
                const bool isCP  = (token.find("cp")!=std::string::npos) &&
                                   (lab.find("constant-b")!=std::string::npos || lab.find("constant b")!=std::string::npos ||
                                    lab.find("legacy")!=std::string::npos || lab=="cp");
                if (isRAW || isCP) { denIdxInV = (int)k; break; }
              }
            }

            // ----- TOP legend(s): enforce a global, consistent ordering everywhere -----
            auto rankOf = [&](const std::string& s)->int {
              std::string t = toLower(s);
              // 0 = Constant-b/CP (always first)
              if (t.find("constant-b")!=std::string::npos || t.find("constant b")!=std::string::npos ||
                  t.find("legacy")!=std::string::npos || t=="cp")
                return 0;
              // 1 = Energy-Dep b (EA energy-only)
              if (t.find("energy-dep")!=std::string::npos || t.find("energy dep")!=std::string::npos ||
                  t.find("e-only")!=std::string::npos || t.find("e only")!=std::string::npos)
                return 1;
              // 2 = Energy + |η|-dep b (EA eta+energy)
              if (t.find("|#eta|-dep")!=std::string::npos || t.find("eta-dep")!=std::string::npos ||
                  t.find("|η|-dep")!=std::string::npos)
                return 2;
              // 3 = RAW / No Correction
              if (t.find("no correction")!=std::string::npos || t=="raw")
                return 3;
              // 4 = anything else
              return 4;
            };

            auto buildOrdered = [&](bool excludeDen)->std::vector<int>{
              std::vector<int> idx;
              idx.reserve(V.size());
              for (int k=0;k<(int)V.size();++k){
                if (excludeDen && k==denIdxInV) continue;
                idx.push_back(k);
              }
              std::sort(idx.begin(), idx.end(), [&](int a,int b){
                int ra = rankOf(V[a].first), rb = rankOf(V[b].first);
                if (ra!=rb) return ra<rb;
                return a<b; // stable tiebreaker
              });
              return idx;
            };

            if ((int)V.size()==3 && denIdxInV>=0)
            {
              const double y1=0.75, y2=0.85;
              const double xL1=0.62, xL2=0.73;
              const double xR1=0.78, xR2=0.86;

              // LEFT: denominator only
              TLegend* legLeft  = new TLegend(xL1,y1,xL2,y2,"","brNDC");
              legLeft->SetBorderSize(0); legLeft->SetFillStyle(0); legLeft->SetTextSize(0.04);
              legLeft->SetEntrySeparation(0.0); legLeft->SetMargin(0.20);
              legLeft->AddEntry(keepTop[denIdxInV].get(), V[denIdxInV].first.c_str(), "p");
              legLeft->Draw();

              // RIGHT: ordered entries (exclude denominator from the list)
              TLegend* legRight = new TLegend(xR1,y1,xR2,y2,"","brNDC");
              legRight->SetBorderSize(0); legRight->SetFillStyle(0); legRight->SetTextSize(0.04);
              legRight->SetEntrySeparation(0.0); legRight->SetMargin(0.20);

              for (int idx : buildOrdered(/*excludeDen=*/true))
                legRight->AddEntry(keepTop[idx].get(), V[idx].first.c_str(), "p");

              legRight->Draw();
            }
            else
            {
              // Single legend case: order everything (denominator may be included here)
              TLegend* legU = new TLegend(0.16,0.73,0.92,0.92,"","brNDC");
              legU->SetBorderSize(0); legU->SetFillStyle(0); legU->SetTextSize(0.033);
              legU->SetNColumns(1);

              for (int idx : buildOrdered(/*excludeDen=*/false))
                legU->AddEntry(keepTop[idx].get(), V[idx].first.c_str(), "p");

              legU->Draw();
            }
            pTop->Modified(); pTop->Update();

          // ---------- BOTTOM: σ / σ_denominator ----------
          pBot->cd();
          struct RS{ std::vector<double> E,S,dS; Color_t col; std::string lab; };
          std::vector<RS> rows;

          for (std::size_t k=0;k<V.size();++k)
          {
            // do NOT include denominator if it is present in V
            if ((int)k==denIdxInV) continue;
            RS rs; rs.col = cols[k%NC]; rs.lab = V[k].first;
            rs.E = eCtr; rs.S.resize(eCtr.size()); rs.dS.resize(eCtr.size());
            for (std::size_t i=0;i<eCtr.size();++i)
            {
              const double s  = SG[k][i];
              const double ds = dSG[k][i];
              const double r  = Sden[i];
              const double dr = dSden[i];
              if (s>0 && r>0)
              {
                const double ratio = s/r;
                const double erel  = std::sqrt( (ds>0? (ds/s)*(ds/s):0.0) + (dr>0? (dr/r)*(dr/r):0.0) );
                rs.S[i]  = ratio;
                rs.dS[i] = ratio * erel;
              } else { rs.S[i]=0.0; rs.dS[i]=0.0; }
            }
            rows.push_back(std::move(rs));
          }

            // Smart zoom: tight range with small pad; clamp around unity if ratios are close to 1
            double yLo=+1e30, yHi=-1e30;
            for (const auto& r : rows)
              for (std::size_t i=0;i<r.S.size();++i) {
                yLo = std::min(yLo, r.S[i] - r.dS[i]);
                yHi = std::max(yHi, r.S[i] + r.dS[i]);
              }
            if (!(yLo < yHi)) { yLo = 0.98; yHi = 1.02; }
            double padBot = 0.03*(yHi - yLo + 1e-12);
            double ymin   = yLo - padBot;
            double ymax   = yHi + padBot;

            // Extra zoom for two-variant CP or RAW ratios (EA-only vs EA|η|): focus near 1.00
            {
              std::string tok = toLower(denomLabelToken ? denomLabelToken : "");
              bool tightCase = (V.size() <= 3) &&
                               (tok.find("cp")  != std::string::npos || tok.find("raw") != std::string::npos);
              if (tightCase && ymin > 0.90 && ymax < 1.05) {
                ymin = std::max(0.94, yLo - 0.01);
                ymax = std::min(1.02, yHi + 0.01);
              }
            }

            // EXTRA BUFFER only for SigmaRatio_EA_vs_CP_* outputs to avoid cutting error bars
            {
              TString outName(outPng.c_str());
              const bool isEAvsCP = outName.Contains("SigmaRatio_EA_vs_CP_");
              if (isEAvsCP) {
                const double extra = 0.015;  // ~1.5% headroom above/below
                ymin -= extra;
                ymax += extra;
              }
            }

            TH1F frL("frL","",1,xMinE,xMaxE);
            frL.SetTitle("");
            frL.GetXaxis()->SetTitle("E  [GeV]");


            // Denominator-aware Y title
            std::string tok = toLower(denomLabelToken ? denomLabelToken : "");
            std::string denomTok = (tok.find("raw")!=std::string::npos) ? "RAW"
                                  : (tok.find("cp") !=std::string::npos) ? "CP"
                                  : denomLabelToken ? denomLabelToken : "REF";
            std::string yTitle = isPhi ?
              std::string("Ratio  #sigma_{RMS}(#Delta#phi) / #sigma_{RMS}(#Delta#phi)_{") + denomTok + "}"
            : std::string("Ratio  #sigma_{RMS}(#Delta#eta) / #sigma_{RMS}(#Delta#eta)_{") + denomTok + "}";
            frL.GetYaxis()->SetTitle(yTitle.c_str());

            frL.SetMinimum(ymin); frL.SetMaximum(ymax);
            frL.Draw("AXIS");

            TLine l1(xMinE,1.0,xMaxE,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

            // Draw ratio points ONLY (no bottom legend)
            for (const auto& r : rows)
            {
              auto g = std::make_unique<TGraphErrors>((int)r.E.size(),
                                                      const_cast<double*>(r.E.data()),
                                                      const_cast<double*>(r.S.data()),
                                                      ex.data(),
                                                      const_cast<double*>(r.dS.data()));
              g->SetMarkerStyle(mk);
              g->SetMarkerSize(1.15);
              g->SetMarkerColor(r.col);
              g->SetLineColor(r.col);   // color error bars
              g->SetLineWidth(1);
              g->Draw("P SAME");
              keepBot.emplace_back(std::move(g));
            }
            pBot->Modified(); pBot->Update();


            // header (smaller font only for SigmaRatio_EA_vs_CP_* outputs)
            c.cd();
            TLatex h; h.SetNDC(); h.SetTextAlign(22);
            {
              TString outName(outPng.c_str());
              const bool isEAvsCP = outName.Contains("SigmaRatio_EA_vs_CP_");
              h.SetTextSize(isEAvsCP ? 0.024 : 0.036);  // tweak 0.024 to taste
            }
            h.DrawLatex(0.52,0.975, header);

            c.SaveAs(outPng.c_str());
            c.Close();

        };

      // ---------- Prepare the specific variant sets you requested ----------
      // φ variants to load (EA|η|+E, EA E-only) - use your established names:
      auto HphiEAetaE   = loadSet("h_phi_diff_cpCorrEA_fitEtaDep_%s");
      auto HphiEAeOnly  = loadSet("h_phi_diff_cpCorrEA_fitEnergyOnly_%s");
      // η counterparts:
      auto HetaEAetaE   = loadSet("h_eta_diff_cpCorrEA_fitEtaDep_%s");
      auto HetaEAeOnly  = loadSet("h_eta_diff_cpCorrEA_fitEnergyOnly_%s");

      // Shorthand labels
      const std::string L_EAetaE  = "Energy + |#eta|-dep b";
      const std::string L_EAeOnly = "Energy-Dep b";
      const std::string L_CP      = "Constant-b (legacy)";
      const std::string L_RAW     = "No Correction";

      // ==================================== φ ====================================
      {
        // 4-variant mean/σ overlays
        drawMuSigmaMulti(dPhi + "/meansigma/MeanSigma_4Variants_phi.png",
                         "#Delta#phi - Mean / RMS vs E (No Corr, Constant b, Energy-Dep b, Energy + |#eta|-dep b)",
                         {
                           {L_EAetaE,  HphiEAetaE},
                           {L_EAeOnly, HphiEAeOnly},
                           {L_CP,      HphiCorr},
                           {L_RAW,     Hphi}
                         });

        // 3-variant (EA|η|+E, CP, RAW)
        drawMuSigmaMulti(dPhi + "/meansigma/MeanSigma_EAetaE_CP_RAW_phi.png",
                         "#Delta#phi - Mean / RMS vs E (No Corr, Constant b, Energy + |#eta|-dep b)",
                         {
                           {L_EAetaE,  HphiEAetaE},
                           {L_CP,      HphiCorr},
                           {L_RAW,     Hphi}
                         });

        // 3-variant (EA E-only, CP, RAW)
        drawMuSigmaMulti(dPhi + "/meansigma/MeanSigma_EAeOnly_CP_RAW_phi.png",
                         "#Delta#phi - Mean / RMS vs E (No Corr, Constant b, Energy-Dep b)",
                         {
                           {L_EAeOnly, HphiEAeOnly},
                           {L_CP,      HphiCorr},
                           {L_RAW,     Hphi}
                         });

        // First-bin multi-overlay
        if (!eEdges.empty())
        {
          const double eLo0 = eEdges.front().first, eHi0 = eEdges.front().second;
          drawFirstBinMulti(dPhi + "/firstbin/FirstBin_phi_EAetaE_CP_RAW.png",
                            "First-bin Overlay (No Corr, Constant b, Energy + |#eta|-dep b)", {
                              {L_EAetaE,  (HphiEAetaE.empty()? nullptr : HphiEAetaE[0])},
                              {L_CP,      (HphiCorr.empty()?    nullptr : HphiCorr[0])},
                              {L_RAW,     (Hphi.empty()?        nullptr : Hphi[0])},
                            }, eLo0, eHi0, /*isPhi*/true);

          drawFirstBinMulti(dPhi + "/firstbin/FirstBin_phi_EAeOnly_CP_RAW.png",
                            "First-bin Overlay (No Corr, Constant b, Energy-Dep b)", {
                              {L_EAeOnly, (HphiEAeOnly.empty()? nullptr : HphiEAeOnly[0])},
                              {L_CP,      (HphiCorr.empty()?    nullptr : HphiCorr[0])},
                              {L_RAW,     (Hphi.empty()?        nullptr : Hphi[0])},
                            }, eLo0, eHi0, /*isPhi*/true);
        }

        // 4×2 table overlays
        drawTableMulti(dPhi + "/tables/Table_phi_EAetaE_CP_RAW.png",
                       "#Delta#phi - RAW/CP/EA|#eta|+E (all energy bins)",
                       { {L_EAetaE, HphiEAetaE}, {L_CP, HphiCorr}, {L_RAW, Hphi} }, /*isPhi*/true);

        drawTableMulti(dPhi + "/tables/Table_phi_EAeOnly_CP_RAW.png",
                       "#Delta#phi - RAW/CP/EA E-only (all energy bins)",
                       { {L_EAeOnly, HphiEAeOnly}, {L_CP, HphiCorr}, {L_RAW, Hphi} }, /*isPhi*/true);

          // Ratio plots vs RAW (keep if you still want RAW-normalized summaries)
          ratioOverlay(dPhi + "/ratios/Ratio_phi_EAetaE_CP_over_RAW.png",
                       "Ratio  #sigma_{RMS}(#Delta#phi) / #sigma_{RMS}(#Delta#phi)_{RAW}",
                       { {L_EAetaE, HphiEAetaE}, {L_CP, HphiCorr} },
                       Hphi);

          ratioOverlay(dPhi + "/ratios/Ratio_phi_EAeOnly_CP_over_RAW.png",
                       "Ratio  #sigma_{RMS}(#Delta#phi) / #sigma_{RMS}(#Delta#phi)_{RAW}",
                       { {L_EAeOnly, HphiEAeOnly}, {L_CP, HphiCorr} },
                       Hphi);

          // σ(E) overlay (top) = {CP, EA E-only, EA |η|+E}; bottom = each / CP
          drawSigmaWithRatioMulti(dPhi + "/sigmaRatio/SigmaRatio_EA_vs_CP_phi.png",
                                  "#Delta#phi  -  #sigma_{RMS}(E) overlay: CP vs EA variants  (top);   ratios / CP  (bottom)",
                                  {
                                    { L_EAeOnly, HphiEAeOnly },
                                    { L_EAetaE,  HphiEAetaE  },
                                    { L_CP,      HphiCorr    }
                                  },
                                  HphiCorr,   // Denominator = CLUS-CP
                                  "CP",
                                  /*isPhi*/true);

          // NEW: this is the file you’re trying to tweak; it uses RAW as denominator.
          // Your function already contains filename-specific extra headroom for this PNG.
          drawSigmaWithRatioMulti(dPhi + "/sigmaRatio/SigmaRatio_EAeOnly_CP_RAW_phi.png",
                                  "#Delta#phi  -  #sigma_{RMS}(E) (top)   and   #sigma/#sigma_{RAW} (bottom)",
                                  {
                                    { L_EAeOnly, HphiEAeOnly },
                                    { L_CP,      HphiCorr    },
                                    { L_RAW,     Hphi        }
                                  },
                                  Hphi,       // Denominator = RAW
                                  "RAW",
                                  /*isPhi*/true);

      }

      // ==================================== η ====================================
      {
        // 4-variant mean/σ overlays
        drawMuSigmaMulti(dEta + "/meansigma/MeanSigma_4Variants_eta.png",
                         "#Delta#eta - Mean / RMS vs E (EA|#eta|+E, EA E-only, CP, RAW)",
                         {
                           {L_EAetaE,  HetaEAetaE},
                           {L_EAeOnly, HetaEAeOnly},
                           {L_CP,      HetaCorr},
                           {L_RAW,     Heta}
                         });

        // 3-variant (EA|η|+E, CP, RAW)
        drawMuSigmaMulti(dEta + "/meansigma/MeanSigma_EAetaE_CP_RAW_eta.png",
                         "#Delta#eta - Mean / RMS vs E (EA|#eta|+E, CP, RAW)",
                         {
                           {L_EAetaE,  HetaEAetaE},
                           {L_CP,      HetaCorr},
                           {L_RAW,     Heta}
                         });

        // 3-variant (EA E-only, CP, RAW)
        drawMuSigmaMulti(dEta + "/meansigma/MeanSigma_EAeOnly_CP_RAW_eta.png",
                         "#Delta#eta - Mean / RMS vs E (EA E-only, CP, RAW)",
                         {
                           {L_EAeOnly, HetaEAeOnly},
                           {L_CP,      HetaCorr},
                           {L_RAW,     Heta}
                         });

        // First-bin multi-overlay
        if (!eEdges.empty())
        {
          const double eLo0 = eEdges.front().first, eHi0 = eEdges.front().second;
          drawFirstBinMulti(dEta + "/firstbin/FirstBin_eta_EAetaE_CP_RAW.png",
                            "First-bin Overlay (EA|#eta|+E, CP, RAW)", {
                              {L_EAetaE,  (HetaEAetaE.empty()? nullptr : HetaEAetaE[0])},
                              {L_CP,      (HetaCorr.empty()?   nullptr : HetaCorr[0])},
                              {L_RAW,     (Heta.empty()?       nullptr : Heta[0])},
                            }, eLo0, eHi0, /*isPhi*/false);

          drawFirstBinMulti(dEta + "/firstbin/FirstBin_eta_EAeOnly_CP_RAW.png",
                            "First-bin Overlay (EA E-only, CP, RAW)", {
                              {L_EAeOnly, (HetaEAeOnly.empty()? nullptr : HetaEAeOnly[0])},
                              {L_CP,      (HetaCorr.empty()?   nullptr : HetaCorr[0])},
                              {L_RAW,     (Heta.empty()?       nullptr : Heta[0])},
                            }, eLo0, eHi0, /*isPhi*/false);
        }

        // 4×2 table overlays
        drawTableMulti(dEta + "/tables/Table_eta_EAetaE_CP_RAW.png",
                       "#Delta#eta - RAW/CP/EA|#eta|+E (all energy bins)",
                       { {L_EAetaE, HetaEAetaE}, {L_CP, HetaCorr}, {L_RAW, Heta} }, /*isPhi*/false);

        drawTableMulti(dEta + "/tables/Table_eta_EAeOnly_CP_RAW.png",
                       "#Delta#eta - RAW/CP/EA E-only (all energy bins)",
                       { {L_EAeOnly, HetaEAeOnly}, {L_CP, HetaCorr}, {L_RAW, Heta} }, /*isPhi*/false);

        // Ratio plots vs RAW
        ratioOverlay(dEta + "/ratios/Ratio_eta_EAetaE_CP_over_RAW.png",
                     "Ratio  #sigma_{RMS}(#Delta#eta) / #sigma_{RMS}(#Delta#eta)_{RAW}",
                     { {L_EAetaE, HetaEAetaE}, {L_CP, HetaCorr} },
                     Heta);

        ratioOverlay(dEta + "/ratios/Ratio_eta_EAeOnly_CP_over_RAW.png",
                     "Ratio  #sigma_{RMS}(#Delta#eta) / #sigma_{RMS}(#Delta#eta)_{RAW}",
                     { {L_EAeOnly, HetaEAeOnly}, {L_CP, HetaCorr} },
                     Heta);

          // CP-normalized: top σ(E) overlay {EA E-only, EA |η|-dep, CP}; bottom: /CP
          drawSigmaWithRatioMulti(dEta + "/sigmaRatio/SigmaRatio_EA_vs_CP_eta.png",
                                  "#Delta#eta  -  #sigma_{RMS}(E) overlay: Legacy vs Energy Dep and With/Without |#eta|-dep  (top);   ratios / legacy  (bottom)",
                                  {
                                    { L_EAeOnly, HetaEAeOnly },
                                    { L_EAetaE,  HetaEAetaE  },
                                    { L_CP,      HetaCorr    }
                                  },
                                  HetaCorr,   // Denominator = CLUS-CP for eta
                                  "CP",
                                  /*isPhi*/false);

          // NEW: RAW-normalized: this produces SigmaRatio_EAeOnly_CP_RAW_eta.png
          // so the filename-based headroom logic can take effect.
          drawSigmaWithRatioMulti(dEta + "/sigmaRatio/SigmaRatio_EAeOnly_CP_RAW_eta.png",
                                  "#Delta#eta  -  #sigma_{RMS}(E) (top)   and   #sigma/#sigma_{RAW} (bottom)",
                                  {
                                    { L_EAeOnly, HetaEAeOnly },
                                    { L_CP,      HetaCorr    },
                                    { L_RAW,     Heta        }
                                  },
                                  Heta,        // Denominator = RAW
                                  "RAW",
                                  /*isPhi*/false);

      }

      std::cout << "[FinalPDCPlots] Wrote curated outputs under: " << finalBase << "\n";
    }



    // ---------------- existing summary calls remain unchanged ----------------
    MakeThreeWaySummaries(
        outDir,
        eEdges,
        xMin, xMax,
        /* φ (cpRaw, cpCorr, cpCorrEA) */
        Hphi, HphiCorr, Hphi_cpCorrEA,
        /* η (cpRaw, cpCorr, cpCorrEA) */
        Heta, HetaCorr, Heta_cpCorrEA,
        /* reuse your existing helpers */
        drawResidualPanel,
        makeEqualHeightTwinPads,
        statsFromHistRange
    );

    MakeFourWaySummaries(
        outDir,
        eEdges,
        xMin, xMax,
        /* φ sets */
        Hphi, HphiCorr, HphiScratchRaw, HphiScratchCorr,
        /* η sets */
        Heta, HetaCorr, HetaScratchRaw, HetaScratchCorr,
        /* reuse your existing lambdas for identical rendering */
        drawResidualPanel,
        makeEqualHeightTwinPads,
        statsFromHistRange
    );

    // =====================================================================
    // NEW: Two summary tables of "smallest RMS per energy bin"
    //      Table A: Across *all* variants that are loaded/read in this job
    //      Table B: Restricted to CLUS-RAW, CLUS-CP, CLUS-CP(EA*) variants
    //      Both are printed to stdout and written as CSV in outDir.
    // =====================================================================

    auto pushIfAny = [&](std::vector<std::pair<std::string, std::vector<TH1F*>>>& dst,
                         const std::string& label, const char* pat)
    {
        auto S = loadSet(pat);
        bool any=false;
        for (auto* h : S) { if (h && h->Integral()>0) { any = true; break; } }
        if (any) dst.emplace_back(label, std::move(S));
    };

    auto computeMinRMSPerBin = [&](const std::vector<std::pair<std::string, std::vector<TH1F*>>>& V,
                                   std::vector<std::string>& bestLabel,
                                   std::vector<double>&     bestSigma)
    {
        bestLabel.assign(eEdges.size(), "--");
        bestSigma.assign(eEdges.size(),  0.0);
        for (std::size_t iE=0; iE<eEdges.size(); ++iE) {
            double best = 1e30; std::string lab = "--";
            for (const auto& kv : V) {
                const auto& label = kv.first;
                const auto& vecH  = kv.second;
                if (iE >= vecH.size()) continue;
                TH1F* h = vecH[iE];
                if (!h || h->Integral()==0) continue;
                double m, me, r, re;
                std::tie(m, me, r, re) = statsFromHistRange(h, xMin, xMax);
                if (r>0.0 && r < best) { best = r; lab = label; }
            }
            if (best < 1e29) { bestSigma[iE] = best; bestLabel[iE] = lab; }
        }
    };

    // ---------- ANSI styles (kept subtle; works in most terminals) ----------
    constexpr const char* RST  = "\033[0m";
    constexpr const char* BOLD = "\033[1m";
    constexpr const char* DIM  = "\033[2m";
    constexpr const char* CYAN = "\033[36m";
    constexpr const char* YELL = "\033[33m";
    constexpr const char* GRN  = "\033[32m";
    constexpr const char* MAG  = "\033[35m";

    auto printMinRMSTable = [&](const char* title,
                                const std::vector<double>&     sPhi,
                                const std::vector<std::string>& lPhi,
                                const std::vector<double>&     sEta,
                                const std::vector<std::string>& lEta)
    {
        const std::string bar(118,'-');
        std::cout << "\n" << CYAN << BOLD << "┌" << bar << "┐" << RST << "\n";
        std::cout << CYAN << BOLD << "│ " << title << RST << "\n";
        std::cout << CYAN << BOLD << "├" << bar << "┤" << RST << "\n";

        // Header
        std::cout << BOLD
                  << "│ " << std::left  << std::setw(15) << "E bin [GeV]"
                  << "│ " << std::right << std::setw(18) << "σ_RMS(Δφ)"
                  << " │ " << std::left << std::setw(42) << "variant (phi)"
                  << "│ " << std::right<< std::setw(18) << "σ_RMS(Δη)"
                  << " │ " << std::left << std::setw(42) << "variant (eta)"
                  << RST << "\n";
        std::cout << DIM << "│ " << std::string(15,'-') << "│ "
                  << std::string(18,'-') << " │ " << std::string(42,'-')
                  << "│ " << std::string(18,'-') << " │ " << std::string(42,'-')
                  << RST << "\n";

        for (std::size_t iE=0; iE<eEdges.size(); ++iE) {
            const double eLo = eEdges[iE].first;
            const double eHi = eEdges[iE].second;
            const double p = (iE<sPhi.size()? sPhi[iE] : 0.0);
            const double e = (iE<sEta.size()? sEta[iE] : 0.0);
            const std::string lp = (iE<lPhi.size()? lPhi[iE] : "--");
            const std::string le = (iE<lEta.size()? lEta[iE] : "--");

            std::ostringstream eb;
            eb << "[" << std::fixed << std::setprecision(0) << eLo << "," << eHi << ")";

            std::cout << "│ " << std::left  << std::setw(15) << eb.str()
                      << "│ " << std::right << std::setw(18) << std::setprecision(6) << p
                      << " │ " << std::left << std::setw(42) << lp
                      << "│ " << std::right << std::setw(18) << std::setprecision(6) << e
                      << " │ " << std::left << std::setw(42) << le
                      << "\n";
        }
        std::cout << CYAN << BOLD << "└" << bar << "┘" << RST << "\n";
    };

    // Pretty, ranked table (smallest → largest) per energy bin for Δφ and Δη
    auto printTableWithSecond = [&](const char* title,
                                    const std::vector<std::vector<std::pair<std::string,double>>>& phiSorted,
                                    const std::vector<std::vector<std::pair<std::string,double>>>& etaSorted)
    {
        const std::string bar(118,'-');
        std::cout << "\n" << CYAN << BOLD << "╔" << bar << "╗" << RST << "\n";
        std::cout << CYAN << BOLD << "║ " << title << RST << "\n";
        std::cout << CYAN << BOLD << "╠" << bar << "╣" << RST << "\n";

        for (std::size_t iE=0;iE<eEdges.size();++iE)
        {
            // Side-by-side ranking lists
            const auto& P = phiSorted[iE];
            const auto& E = etaSorted[iE];
            const std::size_t nRows = std::max(P.size(), E.size());

            // Energy-bin banner
            std::ostringstream eb;
            eb << "[" << std::fixed << std::setprecision(0)
               << eEdges[iE].first << "," << eEdges[iE].second << ") GeV";

            std::cout << YELL << BOLD << "  " << eb.str() << RST << "\n";
            std::cout << BOLD
                      << "  " << std::setw(4)  << "idx"
                      << "  " << std::setw(42) << "Δφ variant"
                      << "  " << std::setw(12) << "σ_RMS(Δφ)"
                      << "   │  "
                      << std::setw(4)  << "idx"
                      << "  " << std::setw(42) << "Δη variant"
                      << "  " << std::setw(12) << "σ_RMS(Δη)"
                      << RST << "\n";
            std::cout << DIM
                      << "  " << std::string(4,'-')  << "  " << std::string(42,'-') << "  "
                      << std::string(12,'-') << "   │  "
                      << std::string(4,'-')  << "  " << std::string(42,'-') << "  "
                      << std::string(12,'-')
                      << RST << "\n";

            for (std::size_t r=0; r<nRows; ++r)
            {
                auto drawCell = [&](const std::vector<std::pair<std::string,double>>& v,
                                    std::size_t idx, bool isPhi)
                {
                    if (idx >= v.size()) {
                        std::cout << std::setw(4)  << ""
                                  << "  " << std::setw(42) << ""
                                  << "  " << std::setw(12) << "";
                        return;
                    }
                    const bool isBest = (idx==0);
                    const bool is2nd  = (idx==1);
                    const char* color = isBest ? GRN : (is2nd ? MAG : RST);

                    std::ostringstream val; val << std::setprecision(6) << v[idx].second;

                    std::cout << color
                              << std::setw(4)  << idx+1 << RST
                              << "  " << color << std::left  << std::setw(42) << v[idx].first << RST
                              << "  " << color << std::right << std::setw(12) << val.str() << RST;
                };

                std::cout << "  ";
                drawCell(P, r, true);
                std::cout << "   │  ";
                drawCell(E, r, false);
                std::cout << "\n";
            }

            // Gap line (2nd - best) if we have at least two entries
            auto gapIf = [&](const std::vector<std::pair<std::string,double>>& v)->std::string{
                if (v.size()<2) return "n/a";
                std::ostringstream oss; oss << std::setprecision(6) << (v[1].second - v[0].second);
                return oss.str();
            };
            std::cout << DIM
                      << "  gap(Δφ) = " << gapIf(P)
                      << "   |   gap(Δη) = " << gapIf(E)
                      << RST << "\n";

            std::cout << CYAN << DIM << std::string(118,'-') << RST << "\n";
        }

        std::cout << CYAN << BOLD << "╚" << bar << "╝" << RST << "\n";
    };

    auto writeCSVwithSecond = [&](const std::string& path,
                                  const std::vector<std::vector<std::pair<std::string,double>>>& phiSorted,
                                  const std::vector<std::vector<std::pair<std::string,double>>>& etaSorted)
    {
        std::ofstream ofs(path);
        ofs << "E_lo,E_hi,"
               "phi_best_label,phi_best_rms,phi_second_label,phi_second_rms,phi_gap,"
               "eta_best_label,eta_best_rms,eta_second_label,eta_second_rms,eta_gap\n";
        for (std::size_t iE=0;iE<eEdges.size();++iE){
            auto get = [&](const std::vector<std::pair<std::string,double>>& v)->std::tuple<std::string,double,std::string,double,double>{
                if (v.empty()) return {"",NAN,"",NAN,NAN};
                const auto& b = v[0];
                if (v.size()<2) return {b.first,b.second,"",NAN,NAN};
                const auto& s = v[1]; return {b.first,b.second,s.first,s.second,s.second-b.second};
            };
            auto [pl,ps,pl2,ps2,pd] = get(phiSorted[iE]);
            auto [el,es,el2,es2,ed] = get(etaSorted[iE]);

            ofs << std::fixed << std::setprecision(3) << eEdges[iE].first << ","
                << eEdges[iE].second << ","
                << "\"" << pl << "\"," << std::setprecision(6) << ps << ","
                << "\"" << pl2 << "\"," << ps2 << "," << (std::isfinite(pd)?pd:NAN) << ","
                << "\"" << el << "\"," << es << ","
                << "\"" << el2 << "\"," << es2 << "," << (std::isfinite(ed)?ed:NAN) << "\n";
        }
        ofs.close();
    };

    // ---------------- Collect ALL available variants (phi & eta) ----------------
    std::vector<std::pair<std::string, std::vector<TH1F*>>> phiAll, etaAll;

    // PDC (from scratch)
    pushIfAny(phiAll, "no corr, scratch",     "h_phi_diff_raw_%s");
    pushIfAny(phiAll, "b(E) corr, scratch",   "h_phi_diff_corr_%s");
    pushIfAny(etaAll, "no corr, scratch",     "h_eta_diff_raw_%s");
    pushIfAny(etaAll, "b(E) corr, scratch",   "h_eta_diff_corr_%s");

    // CLUS: RAW / CP / b(E)
    pushIfAny(phiAll, "no corr, cluster",     "h_phi_diff_cpRaw_%s");
    pushIfAny(phiAll, "CorrectPosition, cluster", "h_phi_diff_cpCorr_%s");
    pushIfAny(phiAll, "b(E) corr, cluster",   "h_phi_diff_cpBcorr_%s");
    pushIfAny(etaAll, "no corr, cluster",     "h_eta_diff_cpRaw_%s");
    pushIfAny(etaAll, "CorrectPosition, cluster", "h_eta_diff_cpCorr_%s");
    pushIfAny(etaAll, "b(E) corr, cluster",   "h_eta_diff_cpBcorr_%s");

    // CLUS: CP(EA) generic + explicit EA flavors
    pushIfAny(phiAll, "CP(EA), cluster",                          "h_phi_diff_cpCorrEA_%s");
    pushIfAny(phiAll, "CP(|vz| + E), cluster",                     "h_phi_diff_cpCorrEA_geom_%s");
    pushIfAny(phiAll, "CP(EA |#eta|+E), cluster",                 "h_phi_diff_cpCorrEA_fitEtaDep_%s");
    pushIfAny(phiAll, "CP(EA E-only), cluster",                   "h_phi_diff_cpCorrEA_fitEnergyOnly_%s");
    pushIfAny(phiAll, "CP(EA #varphi:E-only, #eta:|#eta|+E), cluster",
                       "h_phi_diff_cpCorrEA_fitPhiE_etaEtaDep_%s");

    pushIfAny(etaAll, "CP(EA), cluster",                          "h_eta_diff_cpCorrEA_%s");
    pushIfAny(etaAll, "CP(|vz| + E), cluster",                     "h_eta_diff_cpCorrEA_geom_%s");
    pushIfAny(etaAll, "CP(EA |#eta|+E), cluster",                 "h_eta_diff_cpCorrEA_fitEtaDep_%s");
    pushIfAny(etaAll, "CP(EA E-only), cluster",                   "h_eta_diff_cpCorrEA_fitEnergyOnly_%s");
    pushIfAny(etaAll, "CP(EA #varphi:E-only, #eta:|#eta|+E), cluster",
                       "h_eta_diff_cpCorrEA_fitPhiE_etaEtaDep_%s");

    // --------- Sort all candidates per E-bin (ALL variants) ----------
    // Sort every variant’s RMS per E-bin ascending; return per-bin vectors of {label, rms}
    auto sortCandidatesPerBin =
        [&](const std::vector<std::pair<std::string,std::vector<TH1F*>>>& V)
            -> std::vector<std::vector<std::pair<std::string,double>>>
    {
        const std::size_t NE = eEdges.size();
        std::vector<std::vector<std::pair<std::string,double>>> out(NE);

        for (std::size_t iE=0;iE<NE;++iE){
            for (const auto& kv : V){
                const std::string& lab = kv.first;
                const auto& vec = kv.second;
                if (iE>=vec.size()) continue;
                TH1F* h = vec[iE];
                if (!h || h->Integral()<=0) continue;
                double m,me,r,re;
                std::tie(m,me,r,re) = statsFromHistRange(h, xMin, xMax);
                if (r>0.0) out[iE].emplace_back(lab, r);
            }
            std::sort(out[iE].begin(), out[iE].end(),
                      [](const auto& a, const auto& b){ return a.second < b.second; });
        }
        return out;
    };

    auto phiSorted_All = sortCandidatesPerBin(phiAll);
    auto etaSorted_All = sortCandidatesPerBin(etaAll);

    // ---------------- Collect CLUS-only family (RAW, CP, CP(EA*)) ----------------
    std::vector<std::pair<std::string, std::vector<TH1F*>>> phiClus, etaClus;
    pushIfAny(phiClus, "no corr, cluster",         "h_phi_diff_cpRaw_%s");
    pushIfAny(phiClus, "CorrectPosition, cluster", "h_phi_diff_cpCorr_%s");
    pushIfAny(phiClus, "CP(EA), cluster",          "h_phi_diff_cpCorrEA_%s");
    pushIfAny(phiClus, "CP(|vz| + E), cluster",     "h_phi_diff_cpCorrEA_geom_%s");
    pushIfAny(phiClus, "CP(EA |#eta|+E), cluster", "h_phi_diff_cpCorrEA_fitEtaDep_%s");
    pushIfAny(phiClus, "CP(EA E-only), cluster",   "h_phi_diff_cpCorrEA_fitEnergyOnly_%s");
    pushIfAny(phiClus, "CP(EA #varphi:E-only, #eta:|#eta|+E), cluster",
                        "h_phi_diff_cpCorrEA_fitPhiE_etaEtaDep_%s");

    pushIfAny(etaClus, "no corr, cluster",         "h_eta_diff_cpRaw_%s");
    pushIfAny(etaClus, "CorrectPosition, cluster", "h_eta_diff_cpCorr_%s");
    pushIfAny(etaClus, "CP(EA), cluster",          "h_eta_diff_cpCorrEA_%s");
    pushIfAny(etaClus, "CP(|vz| + E), cluster",     "h_eta_diff_cpCorrEA_geom_%s");
    pushIfAny(etaClus, "CP(EA |#eta|+E), cluster", "h_eta_diff_cpCorrEA_fitEtaDep_%s");
    pushIfAny(etaClus, "CP(EA E-only), cluster",   "h_eta_diff_cpCorrEA_fitEnergyOnly_%s");
    pushIfAny(etaClus, "CP(EA #varphi:E-only, #eta:|#eta|+E), cluster",
                        "h_eta_diff_cpCorrEA_fitPhiE_etaEtaDep_%s");

    // --------- Sort per E-bin (CLUS-only family) ----------
    auto phiSorted_Clus = sortCandidatesPerBin(phiClus);
    auto etaSorted_Clus = sortCandidatesPerBin(etaClus);

    // --------- Pick best & second; compute gaps ----------
    auto pickBestAndSecond =
        [&](const std::vector<std::vector<std::pair<std::string,double>>>& sorted,
            std::vector<std::string>& bestLab,
            std::vector<double>&      bestSig,
            std::vector<std::string>& secLab,
            std::vector<double>&      secSig,
            std::vector<double>&      gap)
    {
        const std::size_t NE = sorted.size();
        bestLab.assign(NE, "-");
        bestSig.assign(NE, std::numeric_limits<double>::quiet_NaN());
        secLab .assign(NE, "-");
        secSig .assign(NE, std::numeric_limits<double>::quiet_NaN());
        gap    .assign(NE, std::numeric_limits<double>::quiet_NaN());

        for (std::size_t iE=0; iE<NE; ++iE) {
            const auto& vec = sorted[iE];
            if (vec.empty()) continue;
            bestLab[iE] = vec[0].first;
            bestSig[iE] = vec[0].second;
            if (vec.size() >= 2) {
                secLab[iE] = vec[1].first;
                secSig[iE] = vec[1].second;
                gap[iE]    = secSig[iE] - bestSig[iE];
            }
        }
    };

    std::vector<std::string> bestLabPhi_All, bestLabEta_All, secLabPhi_All, secLabEta_All;
    std::vector<double>      bestSigPhi_All, bestSigEta_All, secSigPhi_All, secSigEta_All, gapPhi_All, gapEta_All;
    pickBestAndSecond(phiSorted_All, bestLabPhi_All, bestSigPhi_All, secLabPhi_All, secSigPhi_All, gapPhi_All);
    pickBestAndSecond(etaSorted_All, bestLabEta_All, bestSigEta_All, secLabEta_All, secSigEta_All, gapEta_All);

    std::vector<std::string> bestLabPhi_Clus, bestLabEta_Clus, secLabPhi_Clus, secLabEta_Clus;
    std::vector<double>      bestSigPhi_Clus, bestSigEta_Clus, secSigPhi_Clus, secSigEta_Clus, gapPhi_Clus, gapEta_Clus;
    pickBestAndSecond(phiSorted_Clus, bestLabPhi_Clus, bestSigPhi_Clus, secLabPhi_Clus, secSigPhi_Clus, gapPhi_Clus);
    pickBestAndSecond(etaSorted_Clus, bestLabEta_Clus, bestSigEta_Clus, secLabEta_Clus, secSigEta_Clus, gapEta_Clus);


    // --------- Print tables (with 2nd & gap) + full rankings ----------
    printTableWithSecond("TABLE A - Minimum RMS per energy bin across ALL variants (phi & eta)",
                         phiSorted_All, etaSorted_All);

    printTableWithSecond("TABLE B - Minimum RMS per energy bin within CLUS-RAW / CLUS-CP / CLUS-CP(EA*) family (phi & eta)",
                         phiSorted_Clus, etaSorted_Clus);

    // --------- Write CSVs (with 2nd & gap) ----------
    const std::string csvAll  = std::string(outDir) + "/MinRMS_ALL_variants.csv";
    const std::string csvClus = std::string(outDir) + "/MinRMS_CLUS_only.csv";
    writeCSVwithSecond(csvAll,  phiSorted_All,  etaSorted_All);
    writeCSVwithSecond(csvClus, phiSorted_Clus, etaSorted_Clus);

    std::cout << "[Playground] Wrote CSV: " << csvAll  << "\n";
    std::cout << "[Playground] Wrote CSV: " << csvClus << "\n";
    std::cout << "[Playground] Completed Δφ & Δη outputs into: " << outDir << "\n";

}

// ============================================================================
//  MakeResidualsSuite
//  * Robust energy-bin resolution (X-axis title → Y-axis → title → name)
//  * Verbose scan + tabulated summary (console + file)
//  * Stronger error handling and diagnostics
//  * Modularized "PlaneResiduals" writer in a separate function
//  * Keeps your directory structure and overlay products
//  * NEW:
//    - Title/header now appends vz-window (", |v_{z}| < 60 cm" for 0_60vz, otherwise ", |v_{z}| < 10 cm")
//    - /RESIDUALS/phi/Overlays/phiTILTfits: (i) 2-panel μ/σ overlay (all variants), column legend
//      (ii) μ vs ln(E) per-variant fits, and (iii) enumerated fit-constants .txt, as requested
// ============================================================================

namespace {
  // ---------- small utility: string pad ----------
  inline std::string pad(const std::string& s, std::size_t w, bool left=false) {
    if (s.size() >= w) return s;
    if (left) return std::string(w - s.size(), ' ') + s;
    return s + std::string(w - s.size(), ' ');
  }

  // ---------- robust path creation ----------
  inline void ensure_dir(const std::string& p){ EnsureDir(p); }

  // ---------- vz note ----------
  inline const char* vz_note(bool isVz60){ return isVz60 ? ", |v_{z}| < 60 cm" : ", |v_{z}| < 10 cm"; }

  // ---------- color / styles ----------
  const int kMkStyle = 20;  // closed circle
  const int kLnWidth = 2;

  // Stats containers (moved to file scope so helpers can use them)
  struct ResStats { std::vector<double> mean, dmean, rms, drms; };
  using Var2Stats = std::map<std::string, ResStats>;
  // ---------- variants ----------
  struct VarDef {
    std::string key;        // canonical key used in folders & legends
    std::string pretty;     // legend label
    std::string phiPrefix;  // histogram name prefix for φ residuals
    std::string etaPrefix;  // histogram name prefix for η residuals
    Color_t     color;      // ROOT color
  };
//  static const Color_t kVarColor[] = { kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kGray+2 };

  // Palette consistent with your earlier request
  const Color_t C_CLUSRAW        = kBlack;
  const Color_t C_CLUSCP         = kRed+1;
  const Color_t C_CLUSCP_EA      = kMagenta+1;
  const Color_t C_CLUSCP_EA_vz   = kBlue+1;
  const Color_t C_CLUSCP_EA_eta  = kGreen+2;  // purple
  const Color_t C_PDCraw         = kBlack;
  const Color_t C_PDCcor         = kGray+2;
  const Color_t C_CP_BVALS       = kOrange+7;
  const Color_t C_CLUSCP_EA_mx   = kCyan+2;

  const std::vector<VarDef> kVariants = {
    // key                    pretty                                             phiPrefix                                     etaPrefix
    {"CLUSRAW",              "No Correction",                                   "h_phi_diff_cpRaw_",                           "h_eta_diff_cpRaw_",                          C_CLUSRAW},
    {"CLUSCP",               "Constant-b (legacy)",                             "h_phi_diff_cpCorr_",                          "h_eta_diff_cpCorr_",                         C_CLUSCP},
    {"CLUSCP_EA",            "Energy-Dep b",                                    "h_phi_diff_cpCorrEA_fitEnergyOnly_",          "h_eta_diff_cpCorrEA_fitEnergyOnly_",         C_CLUSCP_EA},
    {"CLUSCP_EA_vzDep",      "Energy + |v_{z}|-dep b",                          "h_phi_diff_cpCorrEA_geom_",                   "h_eta_diff_cpCorrEA_geom_",                  C_CLUSCP_EA_vz}, // (ZDEP)
    {"CLUSCP_EA_etaDep",     "Energy + |#eta|-dep b",                           "h_phi_diff_cpCorrEA_fitEtaDep_",              "h_eta_diff_cpCorrEA_fitEtaDep_",             C_CLUSCP_EA_eta},
    {"PDCraw",               "2x2 Block-Local No Correction",                   "h_phi_diff_raw_",                              "h_eta_diff_raw_",                             C_PDCraw},
    {"PDCcor",               "2x2 Block-Local With Corr",                       "h_phi_diff_corr_",                             "h_eta_diff_corr_",                            C_PDCcor},
    {"CLUScpDirectBvals",    "2x2 block b-values applied to No Corr Prod",      "h_phi_diff_cpBcorr_",                          "h_eta_diff_cpBcorr_",                         C_CP_BVALS},
    {"CLUSCP_EA_mixed",      "Energy-dep Only #phi and |#eta|-dep in #eta",     "h_phi_diff_cpCorrEA_fitPhiE_etaEtaDep_",       "h_eta_diff_cpCorrEA_fitPhiE_etaEtaDep_",      C_CLUSCP_EA_mx}
  };

  // ---------- φ / η axis descriptors ----------
  struct KindDef { std::string tag; std::string axisTitle; };
  const KindDef kPhi = {"phi", "#Delta#phi  [rad]"};
  const KindDef kEta = {"eta", "#Delta#eta"};

  // ---------- containers ----------
  using HVec      = std::vector<TH1*>;
  using Var2Hists = std::map<std::string, HVec>;

  // ---------- histogram styling ----------
  inline void style_h(TH1* h, Color_t col) {
    if (!h) return;
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(kMkStyle);
    h->SetLineWidth(kLnWidth);
  }

  // ---------- resolve energy bin from axis/title/name ----------
  // returns [-1] if not found. Records where it was resolved from in *src (optional)
  static int resolve_energy_bin(const TH1* h,
                                const std::string& objName,
                                const std::vector<std::pair<double,double>>& eSlices,
                                const std::vector<double>& eCenters,
                                std::string* src = nullptr)
  {
    if (!h) return -1;
    const int NB = static_cast<int>(eSlices.size());
    auto try_match_range = [&](const std::string& s)->int{
      static const std::regex rg(R"(([-+]?\d+(?:\.\d+)?)\s*<\s*E\s*<\s*([-+]?\d+(?:\.\d+)?))");
      std::smatch m;
      if (std::regex_search(s, m, rg) && m.size() == 3) {
        const double a = std::atof(m[1].str().c_str());
        const double b = std::atof(m[2].str().c_str());
        // exact slice (within tolerance)
        for (int i=0;i<NB;++i){
          if (std::fabs(a - eSlices[i].first)  < 1e-3 &&
              std::fabs(b - eSlices[i].second) < 1e-3) return i;
        }
        // fallback: nearest center to (a+b)/2
        const double c = 0.5*(a+b);
        int best=-1; double dmin=1e9;
        for (int i=0;i<NB;++i){
          const double d = std::fabs(c - eCenters[i]);
          if (d<dmin){dmin=d;best=i;}
        }
        return best;
      }
      return -1;
    };
    auto try_match_equal = [&](const std::string& s)->int{
      static const std::regex rg(R"(E\s*=\s*([-+]?\d+(?:\.\d+)?))");
      std::smatch m;
      if (std::regex_search(s, m, rg) && m.size() == 2) {
        const double c = std::atof(m[1].str().c_str());
        int best=-1; double dmin=1e9;
        for (int i=0;i<NB;++i){
          const double d = std::fabs(c - eCenters[i]);
          if (d<dmin){dmin=d;best=i;}
        }
        return best;
      }
      return -1;
    };
    auto try_match_name = [&](const std::string& s)->int{
      // Accept patterns like "..._2_4", "..._10_12", or "...2to4..."
      static const std::regex rg(R"((?:^|_)(\d+)\s*(?:to|_)\s*(\d+)(?:_|$))", std::regex::icase);
      std::smatch m;
      if (std::regex_search(s, m, rg) && m.size() == 3) {
        const double a = std::atof(m[1].str().c_str());
        const double b = std::atof(m[2].str().c_str());
        // exact slice first
        for (int i=0;i<NB;++i){
          if (std::fabs(a - eSlices[i].first)  < 1e-3 &&
              std::fabs(b - eSlices[i].second) < 1e-3) return i;
        }
        // fallback: nearest center
        const double c = 0.5*(a+b);
        int best=-1; double dmin=1e9;
        for (int i=0;i<NB;++i){
          const double d = std::fabs(c - eCenters[i]);
          if (d<dmin){ dmin=d; best=i; }
        }
        return best;
      }
      return -1;
    };

    // 1) X-axis title
    {
      const std::string s = h->GetXaxis() ? h->GetXaxis()->GetTitle() : "";
      int idx = try_match_range(s); if (idx>=0){ if(src)*src="XAXIS"; return idx; }
      idx = try_match_equal(s);     if (idx>=0){ if(src)*src="XAXIS"; return idx; }
    }
    // 2) Y-axis title (defensive)
    {
      const std::string s = h->GetYaxis() ? h->GetYaxis()->GetTitle() : "";
      int idx = try_match_range(s); if (idx>=0){ if(src)*src="YAXIS"; return idx; }
      idx = try_match_equal(s);     if (idx>=0){ if(src)*src="YAXIS"; return idx; }
    }
    // 3) Histogram title
    {
      const std::string s = h->GetTitle() ? h->GetTitle() : "";
      int idx = try_match_range(s); if (idx>=0){ if(src)*src="TITLE"; return idx; }
      idx = try_match_equal(s);     if (idx>=0){ if(src)*src="TITLE"; return idx; }
    }
    // 4) Object name (E2to4 / E2_4 variants)
    {
      int idx = try_match_name(objName);
      if (idx>=0){ if(src)*src="NAME"; return idx; }
    }
    if (src) *src = "UNRESOLVED";
    return -1;
  }

  // ---------- human-friendly color / pretty lookups ----------
  inline Color_t color_of(const std::string& key){
    for (const auto& v : kVariants) if (v.key==key) return v.color;
    return kBlack;
  }
  inline const char* pretty_of(const std::string& key){
    for (const auto& v : kVariants) if (v.key==key) return v.pretty.c_str();
    return key.c_str();
  }

  // ---------- enumerator name (for fit table) ----------
  inline const char* enum_of(const std::string& key){
    if (key=="CLUSRAW")             return "ETiltVariant::CLUS_RAW";
    if (key=="CLUSCP")              return "ETiltVariant::CLUS_CP";
    if (key=="CLUSCP_EA_vzDep")     return "ETiltVariant::CLUS_CP_EA_FIT_ZDEP";      // (geom → ZDEP)
    if (key=="CLUSCP_EA_etaDep")    return "ETiltVariant::CLUS_CP_EA_FIT_ETADEP";
    if (key=="CLUSCP_EA")           return "ETiltVariant::CLUS_CP_EA_FIT_EONLY";
    if (key=="CLUSCP_EA_mixed")     return "ETiltVariant::CLUS_CP_EA_MIX";
    if (key=="CLUScpDirectBvals")   return "ETiltVariant::CLUS_BCORR";
    if (key=="PDCraw")              return "ETiltVariant::PDC_RAW";
    if (key=="PDCcor")              return "ETiltVariant::PDC_CORR";
    return "ETiltVariant::UNKNOWN";
  }

  // ---------- plane-residuals writer (modularized) ----------
  static void WritePlaneResiduals(const KindDef& kind,
                                  const std::string& variantKey,
                                  const std::vector<std::pair<double,double>>& eSlices,
                                  const HVec& hv,
                                  const std::string& baseOutDir,
                                  bool isVz60) // NEW: control title suffix
  {
    const std::string planeDir = baseOutDir + "/PlaneResiduals";
    ensure_dir(planeDir);

    int savedSingles = 0;

    // -- per-bin PNGs
    for (int i=0;i<(int)hv.size();++i){
      TH1* h = hv[i];
      if (!h) continue;
      const double entries = h->GetEntries();
      if (entries<=0) {
        std::cerr << "[PlaneResiduals][" << kind.tag << "][" << variantKey
                  << "] Ebin " << i << " exists but has zero entries; skipping single plot.\n";
        continue;
      }

      TCanvas c("c_res","",900,700);
      gPad->SetRightMargin(0.04);
      gPad->SetLeftMargin (0.12);
      gPad->SetBottomMargin(0.12);
      style_h(h, color_of(variantKey));

      // axis titles (overwrite any accidental axis text in the source)
      h->GetXaxis()->SetTitle(kind.axisTitle.c_str());
      h->GetYaxis()->SetTitle("Counts");

      const double ymax = 1.25 * std::max(1.0, h->GetMaximum());
      h->GetYaxis()->SetRangeUser(0.0, ymax);
      h->Draw("E1");

      const auto& e = eSlices[i];
      TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.045); hdr.SetTextAlign(13);
      hdr.DrawLatex(0.12, 0.94, Form("%s   —   %.0f < E < %.0f GeV%s",
                                     pretty_of(variantKey), e.first, e.second, vz_note(isVz60)));

      // stats (mean / RMS)
      TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.040); st.SetTextAlign(33);
      st.DrawLatex(0.96, 0.86, Form("<res>=%.5f #pm %.5f", h->GetMean(),   h->GetMeanError()));
      st.DrawLatex(0.96, 0.82, Form("RMS =%.5f #pm %.5f",  h->GetRMS(),    h->GetRMSError()));

      const std::string png = Form("%s/%s_Ebin%02d.png", planeDir.c_str(), variantKey.c_str(), i);
      c.SaveAs(png.c_str());
      ++savedSingles;
    }

    // -- 2×4 mosaic
    bool any = false; for (auto* h : hv) if (h && h->GetEntries()>0) { any=true; break; }
    TCanvas cTab("c_res_table","",1600,1000);
    cTab.Divide(4,2);
    for (int p=1;p<=8;++p){
      TPad* pad = (TPad*)cTab.cd(p);
      pad->SetRightMargin(0.04);
      pad->SetLeftMargin (0.12);
      pad->SetBottomMargin(0.12);
    }
    for (int i=0;i<(int)hv.size();++i){
      cTab.cd(i+1);
      TH1* h = hv[i];
      if (h && h->GetEntries()>0){
        style_h(h, color_of(variantKey));
        h->GetXaxis()->SetTitle(kind.axisTitle.c_str());
        h->GetYaxis()->SetTitle("Counts");
        h->Draw("E1");
        const auto& e = eSlices[i];
        TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.040); hdr.SetTextAlign(13);
        hdr.DrawLatex(0.12, 0.94, Form("%.0f<E<%.0f GeV%s", e.first, e.second, vz_note(isVz60)));
      } else {
        TH1F fr("fr","",10,-1,1);
        fr.SetStats(0);
        fr.GetXaxis()->SetTitle(kind.axisTitle.c_str());
        fr.GetYaxis()->SetTitle("Counts");
        fr.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.040); tx.SetTextAlign(22);
        tx.DrawLatex(0.50, 0.55, "No data");
      }
    }
    const std::string png = Form("%s/%s_Ebins_2x4.png", planeDir.c_str(), variantKey.c_str());
    cTab.SaveAs(png.c_str());

    std::cout << "[PlaneResiduals][" << kind.tag << "][" << variantKey << "] "
              << "Saved " << savedSingles << " per-bin PNGs + mosaic -> " << planeDir << "\n";
  }

  // ---------- pretty printing of the scan results ----------
  struct ScanRow {
    std::string kind;       // phi / eta
    std::string variant;    // key
    std::string suffix;     // "" or "0_60vz"
    int         ebin = -1;  // 0..7 or -1
    std::string ewin;       // "a-b"
    std::string name;       // object name in file
    std::string src;        // XAXIS / YAXIS / TITLE / NAME / UNRESOLVED
    Long64_t    entries = 0;
    double      mean = std::numeric_limits<double>::quiet_NaN();
    double      rms  = std::numeric_limits<double>::quiet_NaN();
    bool        used = false; // got placed in final vector
  };

  static void print_and_write_summary(const std::vector<ScanRow>& rows,
                                      const std::string& outBaseDir)
  {
    if (rows.empty()) return;

    auto print_header = [](){
      std::cout << "\n--------------------------------------------------------------------------------------------\n";
      std::cout << "  KIND  VARIANT                 SUFFIX   BIN   E-Range    ENTRIES     MEAN       RMS   SRC  NAME\n";
      std::cout << "--------------------------------------------------------------------------------------------\n";
    };
    print_header();

    // Also write to a file
    const std::string outText = outBaseDir + "/RESIDUALS/_scan_summary.txt";
    std::ofstream fout(outText);
    if (fout) {
      fout << "KIND,VARIANT,SUFFIX,BIN,E_low,E_high,ENTRIES,MEAN,RMS,SRC,NAME\n";
    }

    auto pr = [&](const ScanRow& r){
      std::ostringstream ew; ew<<std::fixed<<std::setprecision(0);
      ew<<pad(std::to_string((int)std::round(r.ebin>=0 ? 0.0 : 0.0)),0); // dummy to force stream init
      ew.str(""); // clear
    };

    for (const auto& r : rows){
      const std::string eTxt = (r.ebin>=0 && !r.ewin.empty()) ? r.ewin : "--";
      std::cout
        << "  " << pad(r.kind, 4)
        << "  " << pad(r.variant, 22)
        << "  " << pad(r.suffix, 7)
        << "  " << pad( (r.ebin>=0 ? std::to_string(r.ebin) : std::string("--")), 3, true)
        << "  " << pad(eTxt, 9)
        << "  " << pad(std::to_string((long long)r.entries), 9, true)
        << "  " << pad(Form("%.5g", r.mean), 9, true)
        << "  " << pad(Form("%.5g", r.rms ), 8, true)
        << "  " << pad(r.src, 5)
        << "  " << r.name << "\n";

      if (fout){
        double eLo=-1, eHi=-1;
        if (!r.ewin.empty()){
          auto pos = r.ewin.find('-');
          if (pos != std::string::npos){
            eLo = atof(r.ewin.substr(0,pos).c_str());
            eHi = atof(r.ewin.substr(pos+1).c_str());
          }
        }
        fout << r.kind << "," << r.variant << "," << r.suffix << ","
             << (r.ebin>=0 ? std::to_string(r.ebin) : "")
             << "," << (eLo>=0?Form("%.0f",eLo):"")
             << "," << (eHi>=0?Form("%.0f",eHi):"")
             << "," << r.entries << ","
             << r.mean << "," << r.rms  << ","
             << r.src  << ","
             << r.name << "\n";
      }
    }
    std::cout << "--------------------------------------------------------------------------------------------\n";
    std::cout << "[MakeResidualsSuite] Wrote scan summary -> " << outText << "\n";
  }

    // ------------------------------------------------------------------
    //  EmitMinRmsSummary: pretty console tables + CSVs for min RMS per E-bin
    // ------------------------------------------------------------------
    static void EmitMinRmsSummary(const std::string& outDir,
                                  int NB,
                                  const std::vector<std::pair<double,double>>& eSlices,
                                  const std::map<std::string, ResStats>& S_PH,
                                  const std::map<std::string, ResStats>& S_ET)
    {
      ensure_dir(outDir);

      auto collectSortedPerBin =
        [&](const std::vector<std::string>& keys,
            const std::map<std::string,ResStats>& S)
            -> std::vector<std::vector<std::pair<std::string,double>>>
      {
        std::vector<std::vector<std::pair<std::string,double>>> out(NB);
        for (int i=0;i<NB;++i){
          for (const auto& k : keys){
            auto it = S.find(k);
            if (it==S.end()) continue;
            const double r = (i<(int)it->second.rms.size()? it->second.rms[i] : std::numeric_limits<double>::quiet_NaN());
            if (std::isfinite(r))
              out[i].emplace_back(std::string(pretty_of(k)), r);
          }
          std::sort(out[i].begin(), out[i].end(),
                    [](const auto& a, const auto& b){ return a.second < b.second; });
        }
        return out;
      };

      auto writeCSVwithSecond =
        [&](const std::string& path,
            const std::vector<std::vector<std::pair<std::string,double>>>& phiSorted,
            const std::vector<std::vector<std::pair<std::string,double>>>& etaSorted)
      {
        std::ofstream ofs(path);
        ofs << "E_lo,E_hi,phi_best_label,phi_best_rms,phi_second_label,phi_second_rms,phi_gap,"
               "eta_best_label,eta_best_rms,eta_second_label,eta_second_rms,eta_gap\n";
        for (int i=0;i<NB;++i){
          auto gapRow = [&](const std::vector<std::pair<std::string,double>>& v,
                            std::string& l1,double& s1,std::string& l2,double& s2,double& gap)
          {
            l1.clear(); l2.clear(); s1=NAN; s2=NAN; gap=NAN;
            if (v.empty()) return;
            l1 = v[0].first; s1 = v[0].second;
            if (v.size()>1){ l2 = v[1].first; s2 = v[1].second; gap = s2 - s1; }
          };

          std::string pl,pl2,el,el2; double ps,ps2,pg, es,es2,eg;
          gapRow(phiSorted[i], pl,ps, pl2,ps2, pg);
          gapRow(etaSorted[i], el,es, el2,es2, eg);

          ofs << std::fixed << std::setprecision(0) << eSlices[i].first << ","
              << eSlices[i].second << ","
              << "\"" << pl << "\"," << std::setprecision(6) << ps << ","
              << "\"" << pl2 << "\"," << ps2 << "," << (std::isfinite(pg)?pg:NAN) << ","
              << "\"" << el << "\"," << es << ","
              << "\"" << el2 << "\"," << es2 << "," << (std::isfinite(eg)?eg:NAN) << "\n";
        }
        ofs.close();
      };

      auto printRankTable =
        [&](const char* title,
            const std::vector<std::vector<std::pair<std::string,double>>>& phiSorted,
            const std::vector<std::vector<std::pair<std::string,double>>>& etaSorted)
      {
        std::cout << "\n========================================================================================\n";
        std::cout << title << "\n";
        std::cout << "----------------------------------------------------------------------------------------\n";
        for (int i=0;i<NB;++i)
        {
          std::cout << "E bin [" << std::fixed << std::setprecision(0)
                    << eSlices[i].first << "," << eSlices[i].second << ") GeV\n";

          auto printOne = [&](const char* tag,
                              const std::vector<std::pair<std::string,double>>& v)
          {
            std::cout << "  " << tag << ":";
            if (v.empty()) { std::cout << "  (no entries)\n"; return; }
            std::cout << "  best = {" << v[0].first << ", " << std::setprecision(6) << v[0].second << "}";
            if (v.size()>1) {
              const double gap = v[1].second - v[0].second;
              std::cout << "   second = {" << v[1].first << ", " << v[1].second << "}   gap = " << gap;
            }
            std::cout << "\n";
            if (v.size()>2){
              std::cout << "    ranking:";
              for (size_t r=0;r<v.size();++r)
                std::cout << (r?" → ":" ") << "[" << r+1 << "] " << v[r].first << " (" << v[r].second << ")";
              std::cout << "\n";
            }
          };
          printOne("Δφ", phiSorted[i]);
          printOne("Δη", etaSorted[i]);
          std::cout << "----------------------------------------------------------------------------------------\n";
        }
      };

      // ALL variants list (use keys from kVariants)
      std::vector<std::string> keysAll;
      for (const auto& v : kVariants) keysAll.push_back(v.key);

      // CLUS family list (RAW, CP, CP(EA*))
      const std::vector<std::string> keysClus = {
        "CLUSRAW","CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_mixed"
      };

      // Build sorted rankings from internal stats
      auto phiAllSorted  = collectSortedPerBin(keysAll , S_PH);
      auto etaAllSorted  = collectSortedPerBin(keysAll , S_ET);
      auto phiCluSorted  = collectSortedPerBin(keysClus, S_PH);
      auto etaCluSorted  = collectSortedPerBin(keysClus, S_ET);

      // Print to stdout
      printRankTable("TABLE A - Minimum RMS per energy bin across ALL variants (phi & eta)",
                     phiAllSorted, etaAllSorted);
      printRankTable("TABLE B - Minimum RMS per energy bin within CLUS-RAW / CLUS-CP / CLUS-CP(EA*) (phi & eta)",
                     phiCluSorted, etaCluSorted);

      // Write CSVs
      writeCSVwithSecond(outDir + "/MinRMS_ALL_variants.csv",  phiAllSorted, etaAllSorted);
      writeCSVwithSecond(outDir + "/MinRMS_CLUS_only.csv",     phiCluSorted, etaCluSorted);

      std::cout << "[MakeResidualsSuite] Wrote CSVs:\n"
                << "  • " << outDir + "/MinRMS_ALL_variants.csv\n"
                << "  • " << outDir + "/MinRMS_CLUS_only.csv\n";
    }
} // end anon namespace


// ============================================================================
//  MakeResidualsSuite (entry point)
// ============================================================================
static void MakeResidualsSuite(TFile* fin, const std::string& outBaseDir)
{
  if (!fin || fin->IsZombie()) {
    std::cerr << "[MakeResidualsSuite] FATAL: invalid input TFile.\n";
    return;
  }
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // ---------------- energy bins (8) ----------------
  const auto eSlices  = MakeEnergySlices();   // vector<pair<double,double>>
  const auto eCenters = MakeEnergyCenters();  // vector<double>
  const int  NB       = static_cast<int>(eSlices.size());
  if (NB != 8) {
    std::cerr << "[MakeResidualsSuite][WARN] Expected 8 energy bins; found " << NB << ". Proceeding.\n";
  }

  // ---------------- directory scaffold ----------------
  const std::string RESID = outBaseDir + "/RESIDUALS";
  ensure_dir(RESID);
  const std::string DIR_PH  = RESID + "/phi";
  const std::string DIR_ET  = RESID + "/eta";
  ensure_dir(DIR_PH); ensure_dir(DIR_ET);

  auto scaffold_kind = [&](const std::string& base){
    for (const auto& v : kVariants) {
      ensure_dir(base + "/" + v.key);
      ensure_dir(base + "/" + v.key + "/PlaneResiduals");
    }
    // z<=60
    const std::string z60 = base + "/zTo60";
    ensure_dir(z60);
    for (const auto& v : kVariants) {
      ensure_dir(z60 + "/" + v.key + "_vz0to60");
      ensure_dir(z60 + "/" + v.key + "_vz0to60/PlaneResiduals");
    }
    // overlays
    ensure_dir(base + "/Overlays/firstBin");
    ensure_dir(base + "/Overlays/meansigma");
    ensure_dir(base + "/Overlays/sigmaRatio");
    ensure_dir(base + "/Overlays_0_60vz/firstBin");
    ensure_dir(base + "/Overlays_0_60vz/meansigma");
    ensure_dir(base + "/Overlays_0_60vz/sigmaRatio");

    // NEW: φ tilt-fit products
    if (base == DIR_PH) {
      ensure_dir(base + "/Overlays/phiTILTfits");
    }
  };
  scaffold_kind(DIR_PH);
  scaffold_kind(DIR_ET);

  // ---------------- collect histograms (robust parse) ----------------
  using Var2Hists = std::map<std::string, HVec>;
  Var2Hists PH, PH_60, ET, ET_60;              // variant -> [NB] hists
  for (const auto& v : kVariants) {
    PH    [v.key] = HVec(NB,nullptr);
    PH_60 [v.key] = HVec(NB,nullptr);
    ET    [v.key] = HVec(NB,nullptr);
    ET_60 [v.key] = HVec(NB,nullptr);
  }

  std::vector<ScanRow> scanRows;
  scanRows.reserve(1024);

  // iterate top-level keys
  int nKeys = 0, nH1 = 0;
  TIter next(fin->GetListOfKeys());
  while (TKey* k = static_cast<TKey*>(next())) {
    ++nKeys;
    const TString nm = k->GetName();
    TObject* obj = k->ReadObj();
    std::unique_ptr<TObject> guard(obj);

    TH1* h = dynamic_cast<TH1*>(obj);
    if (!h) continue;
    ++nH1;

    const std::string name = nm.Data();
    const bool isVz60 = (name.size() >= 7 && name.rfind("_0_60vz") == name.size()-7);

    // Determine kind & variant by prefix match
    auto try_attach = [&](const KindDef& kind, const VarDef& vdef)->bool {
      const std::string& pref = (kind.tag=="phi" ? vdef.phiPrefix : vdef.etaPrefix);
      if (!TString(name.c_str()).BeginsWith(pref.c_str())) return false;

      // Skip per-vz debug histos like "..._vz0_5", "..._vz10_15", etc.
      // Keep only the special suite "..._0_60vz".
      if (name.find("_vz") != std::string::npos && name.find("_0_60vz") == std::string::npos)
          return false;

      std::string src;
      const int ib = resolve_energy_bin(h, name, eSlices, eCenters, &src);
      const auto& targetES = eSlices;
      const std::string ewin = (ib>=0 && ib<(int)targetES.size())
                               ? (Form("%.0f-%.0f", targetES[ib].first, targetES[ib].second))
                               : "";

      ScanRow row;
      row.kind    = kind.tag;
      row.variant = vdef.key;
      row.suffix  = (isVz60 ? "0_60vz" : "");
      row.ebin    = ib;
      row.ewin    = ewin;
      row.name    = name;
      row.src     = src.empty() ? "UNRES" : src;
      row.entries = (Long64_t)h->GetEntries();
      row.mean    = h->GetMean();
      row.rms     = h->GetRMS();

      bool placed = false;
      if (ib>=0 && ib<NB) {
        TH1* hc = dynamic_cast<TH1*>(h->Clone());
        if (hc) {
          hc->SetDirectory(nullptr);
          HVec& vec = (kind.tag=="phi"
                        ? (isVz60 ? PH_60[vdef.key] : PH[vdef.key])
                        : (isVz60 ? ET_60[vdef.key] : ET[vdef.key]));
          if (vec[ib]!=nullptr) {
            std::cerr << "[MakeResidualsSuite][WARN] Duplicate for ("<<kind.tag<<","<<vdef.key<<","<<ib
                      << (isVz60? ",0_60vz":"") << "). Keeping the first and ignoring: " << name << "\n";
            delete hc;
          } else {
            vec[ib] = hc;
            placed = true;
          }
        }
      } else {
        std::cerr << "[MakeResidualsSuite][WARN] Could not resolve energy bin for " << name
                  << " (src=" << row.src << "). Skipping placement.\n";
      }
      row.used = placed;
      scanRows.push_back(std::move(row));
      return true;
    };

    for (const auto& vdef : kVariants) {
      if (try_attach(kPhi, vdef)) break;
      if (try_attach(kEta, vdef)) break;
    }
  }

  std::cout << "[MakeResidualsSuite] Scanned " << nKeys << " keys; usable TH1: " << nH1 << "\n";

  // ------------- verbose summary table (console + file) -------------
  print_and_write_summary(scanRows, outBaseDir);

  // ------------- compute stats (mean/RMS) for overlays -------------
  auto init_stats = [&](Var2Stats& S){
    for (const auto& v : kVariants) {
      S[v.key].mean .assign(NB, std::numeric_limits<double>::quiet_NaN());
      S[v.key].dmean.assign(NB, std::numeric_limits<double>::quiet_NaN());
      S[v.key].rms  .assign(NB, std::numeric_limits<double>::quiet_NaN());
      S[v.key].drms .assign(NB, std::numeric_limits<double>::quiet_NaN());
    }
  };
  auto compute_stats = [&](const Var2Hists& M, Var2Stats& S){
    init_stats(S);
    for (const auto& kv : M){
      const std::string& key = kv.first;
      const HVec& vec = kv.second;
      for (int i=0;i<NB;++i){
        TH1* h = vec[i];
        if (!h || h->GetEntries()<=0) continue;
        S[key].mean[i] = h->GetMean();
        S[key].dmean[i]= h->GetMeanError();
        S[key].rms [i] = h->GetRMS();
        // Prefer StdDevError; fallback to RMSError
        double eSD = h->GetStdDevError();
        double eRM = h->GetRMSError();
        S[key].drms[i]= (eSD>0? eSD : (eRM>0? eRM : 0.0));
      }
    }
  };

  Var2Stats S_PH, S_PH_60, S_ET, S_ET_60;
  compute_stats(PH,    S_PH);
  compute_stats(PH_60, S_PH_60);
  compute_stats(ET,    S_ET);
  compute_stats(ET_60, S_ET_60);

  // ------------- per-variant plane residuals (defaults + 0-60) -------------
  for (const auto& v : kVariants) {
    // φ defaults
    if (PH.count(v.key)) {
      const std::string base = DIR_PH + "/" + v.key;
      WritePlaneResiduals(kPhi, v.key, eSlices, PH.at(v.key), base, /*isVz60=*/false);
    }
    // φ 0–60
    if (PH_60.count(v.key)) {
      const std::string base = DIR_PH + "/zTo60/" + v.key + "_vz0to60";
      WritePlaneResiduals(kPhi, v.key, eSlices, PH_60.at(v.key), base, /*isVz60=*/true);
    }
    // η defaults
    if (ET.count(v.key)) {
      const std::string base = DIR_ET + "/" + v.key;
      WritePlaneResiduals(kEta, v.key, eSlices, ET.at(v.key), base, /*isVz60=*/false);
    }
    // η 0–60
    if (ET_60.count(v.key)) {
      const std::string base = DIR_ET + "/zTo60/" + v.key + "_vz0to60";
      WritePlaneResiduals(kEta, v.key, eSlices, ET_60.at(v.key), base, /*isVz60=*/true);
    }
  }

  // ---------------- overlays (same groupings as before) ----------------
  auto mkGraph = [&](const std::vector<double>& x,
                     const std::vector<double>& y,
                     const std::vector<double>& ye,
                     Color_t c)->std::unique_ptr<TGraphErrors>
  {
    std::vector<double> xs, ys, yes, xerr;
    xs.reserve(x.size()); ys.reserve(x.size()); yes.reserve(x.size()); xerr.reserve(x.size());
    for (size_t i=0;i<x.size();++i){
      const double yy  = y[i];
      const double yee = ye[i];
      if (std::isfinite(yy)) { xs.push_back(x[i]); ys.push_back(yy); yes.push_back(std::isfinite(yee)?yee:0.0); xerr.push_back(0.0); }
    }
    auto g = std::make_unique<TGraphErrors>((int)xs.size(),
                                            xs.data(), ys.data(),
                                            xerr.data(), yes.data());
    g->SetMarkerStyle(kMkStyle);
    g->SetLineWidth(kLnWidth);
    g->SetMarkerColor(c);
    g->SetLineColor(c);
    return g;
  };

  auto meansigma_overlay = [&](const KindDef& kind,
                                 const Var2Stats& S,
                                 const std::vector<std::string>& which,
                                 const std::string& outPng,
                                 const char* vzNoteStr)
    {
      // Canvas + pads like the reference block
      TCanvas c("cMuSg","",1000,850);
      TPad top ("top" ,"",0,0.50,1,1);  top.SetBottomMargin(0.06); top.SetLeftMargin(0.15); top.SetRightMargin(0.06); top.Draw();
      TPad bot ("bot" ,"",0,0.00,1,0.50); bot.SetTopMargin(0.06);  bot.SetLeftMargin(0.15); bot.SetRightMargin(0.06); bot.SetBottomMargin(0.18); bot.Draw();


      const double xMin = eSlices.front().first - 0.5;
      const double xMax = eSlices.back().second + 0.5;

      // ---------------------- top: mean (with dashed y=0 and symmetric scaling) ----------------------
      top.cd();

      // Compute tight mean range from data ± errors
      double muLo = +1e30, muHi = -1e30;
      for (const auto& k : which){
        const auto it = S.find(k); if (it==S.end()) continue;
        for (size_t i=0;i<it->second.mean.size();++i){
          const double m  = it->second.mean[i];
          const double dm = (i<it->second.dmean.size()? it->second.dmean[i] : 0.0);
          if (!std::isfinite(m)) continue;
          muLo = std::min(muLo, m - (std::isfinite(dm)?dm:0.0));
          muHi = std::max(muHi, m + (std::isfinite(dm)?dm:0.0));
        }
      }
      if (!(muLo < muHi)) { muLo = -1e-4; muHi = +1e-4; }

      // Symmetric about zero with modest pad (matches your reference look)
      const double absMax = std::max(std::fabs(muLo), std::fabs(muHi));
      const double padMu  = 0.15 * (2.0*absMax);
      const double yMinU  = -(absMax + padMu);
      const double yMaxU  =  (absMax + padMu);

        TH1F frU("frU","",1,xMin,xMax);
        frU.SetTitle(Form("Residual means vs E (%s)%s;E  [GeV];#mu(%s) [rad]",
                          kind.tag.c_str(),
                          (vzNoteStr?vzNoteStr:""),
                          (kind.tag=="phi"?"#Delta#phi":"#Delta#eta")));
        frU.SetMinimum(yMinU);
        frU.SetMaximum(yMaxU);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);

        // >>> NEW: larger Y title for the TOP panel (tweak these two numbers as you like)
        frU.GetYaxis()->SetTitleSize(0.058);
        frU.GetYaxis()->SetTitleOffset(0.90);
        // (optional) slightly bigger Y labels so the scale matches
        frU.GetYaxis()->SetLabelSize(0.042);

        frU.Draw("AXIS");


      // y=0 dashed line (requested)
      TLine l0(xMin,0.0,xMax,0.0); l0.SetLineStyle(2); l0.SetLineColor(kGray+2); l0.Draw();

        // helper: override only the vz label to "z_{vtx}"
        auto legend_label = [&](const std::string& key)->const char* {
          if (key == "CLUSCP_EA_vzDep") return "Energy + |z_{vtx}|-dep b";
          return pretty_of(key);
        };

        // draw graphs and keep a handle by key so we can place them in custom columns
        std::vector<std::unique_ptr<TGraphErrors>> keepU;
        std::map<std::string,TGraphErrors*> gByKey;

        for (const auto& k : which){
          const auto it = S.find(k);
          if (it==S.end()) continue;
          auto g = mkGraph(eCenters, it->second.mean, it->second.dmean, color_of(k));
          if (g->GetN()>0) {
            g->SetLineColor(g->GetMarkerColor()); // show error bars
            g->SetLineWidth(1);
            g->Draw("P SAME");
            gByKey[k] = g.get();
            keepU.emplace_back(std::move(g));
          }
        }

        // Mirror the μ(E) legend layout so the columns match visually
        TLegend* legLeft  = new TLegend(0.46, 0.82, 0.62, 0.88, "", "brNDC");
        TLegend* legRight = new TLegend(0.65, 0.73, 0.92, 0.9, "", "brNDC");
        for (auto* L : {legLeft, legRight}) {
          L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.040);
          L->SetNColumns(1); L->SetEntrySeparation(0.0); L->SetMargin(0.18);
        }

        auto addIf = [&](TLegend* L, const char* key){
          auto it = gByKey.find(key);
          if (it != gByKey.end() && it->second) L->AddEntry(it->second, legend_label(key), "p");
        };

        // Base columns always start with RAW (left) and CP (right)
        std::vector<std::string> leftKeys  = {"CLUSRAW"};
        std::vector<std::string> rightKeys = {"CLUSCP"};

        // Detect which EA flavours are present
        const bool hasEta   = gByKey.count("CLUSCP_EA_etaDep");
        const bool hasVz    = gByKey.count("CLUSCP_EA_vzDep");
        const bool hasEonly = gByKey.count("CLUSCP_EA");
        const int  nEA      = (hasEta?1:0) + (hasVz?1:0) + (hasEonly?1:0);

        // If there is exactly ONE EA variant, put it on the RIGHT under CP
        // (this fixes the CLUSRAW_CLUSCP_CLUSCP_EA_etaDep.png case)
        if (nEA == 1) {
          if (hasEta)   rightKeys.push_back("CLUSCP_EA_etaDep");
          if (hasVz)    rightKeys.push_back("CLUSCP_EA_vzDep");
          if (hasEonly) rightKeys.push_back("CLUSCP_EA");
        } else {
          // Otherwise use the original two-column layout
          if (hasEta)   leftKeys.push_back("CLUSCP_EA_etaDep");
          if (hasEonly) leftKeys.push_back("CLUSCP_EA");
          if (hasVz)    rightKeys.push_back("CLUSCP_EA_vzDep");
        }

        // Emit in order
        for (const auto& k : leftKeys)  addIf(legLeft,  k.c_str());
        for (const auto& k : rightKeys) addIf(legRight, k.c_str());

        legLeft->Draw();
        legRight->Draw();

      // ---------------------- bottom: RMS (tight zoom, slight headroom) ----------------------
      bot.cd();
      double sLo = +1e30, sHi = -1e30;
      for (const auto& k : which){
        const auto it = S.find(k); if (it==S.end()) continue;
        for (size_t i=0;i<it->second.rms.size();++i){
          const double s  = it->second.rms[i];
          const double ds = (i<it->second.drms.size()? it->second.drms[i] : 0.0);
          if (!std::isfinite(s)) continue;
          sLo = std::min(sLo, s - (std::isfinite(ds)?ds:0.0));
          sHi = std::max(sHi, s + (std::isfinite(ds)?ds:0.0));
        }
      }
      if (!(sLo < sHi)) { sLo = 0.0; sHi = 0.06; }
      const double padS = 0.12*(sHi - sLo + 1e-12);

        TH1F frL("frL","",1,xMin,xMax);
        frL.SetTitle(Form(";E  [GeV];#sigma_{RMS}(%s) [rad]", (kind.tag=="phi"?"#Delta#phi":"#Delta#eta")));
        frL.SetMinimum(std::max(0.0, sLo - padS));
        frL.SetMaximum(sHi + padS);

        // >>> NEW: larger Y title for the BOTTOM panel (tweak these two numbers as you like)
        frL.GetYaxis()->SetTitleSize(0.06);
        frL.GetYaxis()->SetTitleOffset(0.83);
        // (optional) matching label size
        frL.GetYaxis()->SetLabelSize(0.04);

        frL.Draw("AXIS");


        std::vector<std::unique_ptr<TGraphErrors>> keepL;
        for (const auto& k : which){
          const auto it = S.find(k);
          if (it==S.end()) continue;
          auto g = mkGraph(eCenters, it->second.rms, it->second.drms, color_of(k));
          if (g->GetN()>0) {
            g->SetLineColor(g->GetMarkerColor()); // show error bars
            g->SetLineWidth(1);
            g->Draw("P SAME");
            keepL.emplace_back(std::move(g));
          }
        }

      c.SaveAs(outPng.c_str());
    };
    
    
    std::function<void(const KindDef&,
                       const Var2Stats&,
                       const std::vector<std::string>&,
                       const std::string&,
                       const std::string&,
                       const char*)>
    sigma_ratio_overlay = [&](const KindDef& kind,
                              const Var2Stats& S,
                              const std::vector<std::string>& which,
                              const std::string& baselineKey,
                              const std::string& outPng,
                              const char* vzNoteStr)
    {
      const auto itB = S.find(baselineKey);
      if (itB==S.end()) return;
      const auto& B = itB->second;

      TCanvas c("cSR","",1000,850);
      TPad top ("top" ,"",0,0.35,1,1);  top.SetBottomMargin(0.02); top.SetLeftMargin(0.15); top.SetRightMargin(0.06); top.Draw();
      TPad bot ("bot" ,"",0,0.00,1,0.35); bot.SetTopMargin(0.04);  bot.SetLeftMargin(0.15); bot.SetRightMargin(0.06); bot.SetBottomMargin(0.18); bot.Draw();

      const double xMin = eSlices.front().first - 0.5;
      const double xMax = eSlices.back().second + 0.5;

      // ---------------- top: RMS(E) ----------------
      top.cd();

      // Tight auto-zoom using values ± errors
      double sLo = +1e30, sHi = -1e30;
      for (const auto& k : which){
        const auto it = S.find(k); if (it==S.end()) continue;
        for (size_t i=0;i<it->second.rms.size();++i){
          const double s  = it->second.rms[i];
          const double ds = (i<it->second.drms.size()? it->second.drms[i] : 0.0);
          if (!std::isfinite(s)) continue;
          sLo = std::min(sLo, s - (std::isfinite(ds)?ds:0.0));
          sHi = std::max(sHi, s + (std::isfinite(ds)?ds:0.0));
        }
      }
      if (!(sLo < sHi)) { sLo = 0.0; sHi = 0.06; }
      const double padS = 0.4*(sHi - sLo + 1e-12);

      TH1F frU("frU","",1,xMin,xMax);
      frU.SetTitle(Form("RMS vs E (%s)%s;E  [GeV];#sigma_{RMS}(%s) [rad]",
                        kind.tag.c_str(),
                        (vzNoteStr?vzNoteStr:""),
                        (kind.tag=="phi"?"#Delta#phi":"#Delta#eta")));
      frU.SetMinimum(std::max(0.0, sLo - padS));
      frU.SetMaximum(sHi + padS);
      frU.GetXaxis()->SetLabelSize(0.0);
      frU.GetXaxis()->SetTitleSize(0.0);
      frU.GetXaxis()->SetTickLength(0.0);
      frU.Draw("AXIS");

        auto legend_label = [&](const std::string& key)->const char* {
          if (key == "CLUSCP_EA_vzDep") return "Energy + |z_{vtx}|-dep b";
          return pretty_of(key);
        };

        std::vector<std::unique_ptr<TGraphErrors>> keepU;
        std::map<std::string,TGraphErrors*> gByKey;

        // Define once near the start of sigma_ratio_overlay() so it's visible everywhere in this function:
        const TString outNameTL(outPng.c_str());

        // TOP: draw graphs (with optional left/right bin offsets)
        {
          const bool useBinOffsets =
              outNameTL.Contains("sigma_vsE_ratio_to_CLUSCP_EA_vz_eta_withEonly.png");

            // ---- Tunable, even-in-bin offsets (TOP & BOTTOM) -----------------------
            // Single knob: set the per-slot step as a fraction of the bin width.
            // Suggested range: 0.06–0.18 (keeps markers well within the bin).
            const double kBinOffset = 0.12;

            // Discrete slot per variant: 0=center (kept empty for readability), ±1, ±2 ...
            auto slotOf = [&](const std::string& key)->int {
              if (key == "CLUSCP")            return -2; // Constant-b (legacy)
              if (key == "CLUSCP_EA_etaDep")  return -1; // Energy + |η|-dep b
              if (key == "CLUSCP_EA_vzDep")   return +1; // Energy + |z_vtx|-dep b
              if (key == "CLUSCP_EA")         return +2; // Energy-Dep b (E-only)
              return 0;                                   // others (if any) on center (usually unused)
            };

            // Fractional X-offset inside the bin; multiply by bin width (eHi-eLo)
            auto xOffsetFrac = [&](const std::string& key)->double {
              return slotOf(key) * kBinOffset;
            };

          for (const auto& k : which){
            if (!S.count(k)) continue;

            // Build shifted X, Y, dY
            std::vector<double> X, Y, dX, dY;
            X.reserve(NB); Y.reserve(NB); dX.assign(NB,0.0); dY.reserve(NB);
            for (int i=0;i<NB;++i){
              const double y  = (i<(int)S.at(k).rms.size()? S.at(k).rms[i]  : std::numeric_limits<double>::quiet_NaN());
              const double dy = (i<(int)S.at(k).drms.size()? S.at(k).drms[i]: 0.0);
              if (!std::isfinite(y)) continue;
              const double eLo  = eSlices[i].first;
              const double eHi  = eSlices[i].second;
              const double eCtr = 0.5*(eLo + eHi);
              const double eWid = (eHi - eLo);
                const double xOff = (useBinOffsets ? xOffsetFrac(k) : 0.0) * eWid;
              X.push_back(eCtr + xOff);
              Y.push_back(y);
              dY.push_back(dy);
            }

            auto g = std::make_unique<TGraphErrors>(
              (int)X.size(),
              X.empty()?nullptr:&X[0],
              Y.empty()?nullptr:&Y[0],
              dX.empty()?nullptr:&dX[0],
              dY.empty()?nullptr:&dY[0]
            );
              if (g->GetN()>0){
                g->SetMarkerColor( color_of(k) );
                g->SetLineColor(   color_of(k) );
                g->SetLineWidth(1);
                g->SetMarkerStyle(kMkStyle);   // ensure visible filled circles
                g->SetMarkerSize(1.10);        // enlarge so they show up over error bars
                g->Draw("P SAME");
                gByKey[k] = g.get();
                keepU.emplace_back(std::move(g));
              }
          }
        }

        // outNameTL is already defined above.
        // Match both the base name and the "_withEonly" variant by testing the stem
        const bool useSingleRightLegend =
            outNameTL.Contains("sigma_vsE_ratio_to_CLUSCP_EA_vz_eta");


        auto addIf = [&](TLegend* L, const char* key){
          auto it = gByKey.find(key);
          if (it != gByKey.end() && it->second)
            L->AddEntry(it->second, legend_label(key), "p");
        };

        if (useSingleRightLegend) {
          // --- Single column legend in the top-right (requested) ---
          TLegend* leg1 = new TLegend(0.62, 0.68, 0.92, 0.88, "", "brNDC");
          leg1->SetBorderSize(0);
          leg1->SetFillStyle(0);
          leg1->SetTextSize(0.044);
          leg1->SetNColumns(1);
          leg1->SetEntrySeparation(0.0);
          leg1->SetMargin(0.18);

            // Deterministic order: CP → EA(|η|+E) → EA(|z|+E) → EA(E-only) → RAW
            addIf(leg1, "CLUSCP");
            addIf(leg1, "CLUSCP_EA_etaDep");
            addIf(leg1, "CLUSCP_EA_vzDep");
            addIf(leg1, "CLUSCP_EA");       // Energy-Dep b (E-only)
            addIf(leg1, "CLUSRAW");


          leg1->Draw();
        } else {
          // --- Default: two-column, row-aligned legends ---
          TLegend* legLeft  = new TLegend(0.4, 0.73, 0.52, 0.92, "", "brNDC");
          TLegend* legRight = new TLegend(0.56, 0.73, 0.92, 0.9, "", "brNDC");
          for (auto* L : {legLeft, legRight}) {
            L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.040);
            L->SetNColumns(1); L->SetEntrySeparation(0.0); L->SetMargin(0.18);
          }

          std::vector<std::string> leftKeys  = {"CLUSRAW"};
          std::vector<std::string> rightKeys = {"CLUSCP"};

          const bool hasEta   = gByKey.count("CLUSCP_EA_etaDep");
          const bool hasVz    = gByKey.count("CLUSCP_EA_vzDep");
          const bool hasEonly = gByKey.count("CLUSCP_EA");
          const int  nEA      = (hasEta?1:0) + (hasVz?1:0) + (hasEonly?1:0);

          if (nEA == 1) {
            if (hasEta)   rightKeys.push_back("CLUSCP_EA_etaDep");
            if (hasVz)    rightKeys.push_back("CLUSCP_EA_vzDep");
            if (hasEonly) rightKeys.push_back("CLUSCP_EA");
          } else {
            if (hasEta)   leftKeys.push_back("CLUSCP_EA_etaDep");
            if (hasEonly) leftKeys.push_back("CLUSCP_EA");
            if (hasVz)    rightKeys.push_back("CLUSCP_EA_vzDep");
          }

          for (const auto& k : leftKeys)  addIf(legLeft,  k.c_str());
          for (const auto& k : rightKeys) addIf(legRight, k.c_str());

          legLeft->Draw();
          legRight->Draw();
        }
      // ---------------- bottom: RMS / baseline ----------------
      bot.cd();

      // Build ratios to baseline
      std::vector<std::unique_ptr<TGraphErrors>> keepL;
      std::vector<double> yminCand, ymaxCand;

        // BOTTOM: build ratio graphs (use the same offset pattern as TOP)
        {
          const bool useBinOffsets =
              outNameTL.Contains("sigma_vsE_ratio_to_CLUSCP_EA_vz_eta_withEonly.png");

          // Bottom-panel offset scheme (mirrors TOP). Center slot intentionally left empty.
          const double kBinOffsetBottom = 0.12;
          auto slotOfBot = [&](const std::string& key)->int {
            if (key == "CLUSCP")            return -2; // Constant-b (legacy)
            if (key == "CLUSCP_EA_etaDep")  return -1; // Energy + |η|-dep b
            if (key == "CLUSCP_EA_vzDep")   return +1; // Energy + |z_vtx|-dep b
            if (key == "CLUSCP_EA")         return +2; // Energy-Dep b (E-only)
            return 0;
          };
          auto xOffsetFracBot = [&](const std::string& key)->double {
            return slotOfBot(key) * kBinOffsetBottom;
          };

          for (const auto& k : which){
            if (k==baselineKey || !S.count(k)) continue;

            std::vector<double> X, Y, dX, dY;
            X.reserve(NB); Y.reserve(NB); dX.assign(NB,0.0); dY.reserve(NB);


            for (int i=0;i<NB;++i){
              const double s  = S.at(k).rms[i];
              const double se = S.at(k).drms[i];
              const double r  = B.rms[i];
              const double re = B.drms[i];
              if (std::isfinite(s) && std::isfinite(r) && r>0.0){
                const double ratio = s/r;
                const double erel2 = ( (se>0&&s>0? (se/s)*(se/s) : 0.0) +
                                       (re>0&&r>0? (re/r)*(re/r) : 0.0) );
                const double eLo  = eSlices[i].first;
                const double eHi  = eSlices[i].second;
                const double eCtr = 0.5*(eLo + eHi);
                const double eWid = (eHi - eLo);
                const double xOff = (useBinOffsets ? xOffsetFracBot(k) : 0.0) * eWid;


                X.push_back(eCtr + xOff);
                Y.push_back(ratio);
                dY.push_back(ratio * std::sqrt(erel2));
                yminCand.push_back(ratio - dY.back());
                ymaxCand.push_back(ratio + dY.back());
              }
            }

            auto g = std::make_unique<TGraphErrors>(
              (int)X.size(),
              X.empty()?nullptr:&X[0],
              Y.empty()?nullptr:&Y[0],
              dX.empty()?nullptr:&dX[0],
              dY.empty()?nullptr:&dY[0]
            );
              if (g->GetN()>0){
                g->SetMarkerColor( color_of(k) );
                g->SetLineColor(   color_of(k) ); // keep error bars
                g->SetLineWidth(1);
                g->SetMarkerStyle(kMkStyle);   // same visible marker as top
                g->SetMarkerSize(1.10);
                keepL.emplace_back(std::move(g));
              }
          }
        }

      // Auto y-range centered around 1 with small padding; fallback if empty
      double yLo = 0.98, yHi = 1.02;
      if (!yminCand.empty() && !ymaxCand.empty()){
        yLo = *std::min_element(yminCand.begin(), yminCand.end());
        yHi = *std::max_element(ymaxCand.begin(), ymaxCand.end());
        const double pad = 0.07*(yHi - yLo + 1e-12);
        yLo -= pad; yHi += pad;
      }

        TH1F frL("frL","",1,xMin,xMax);
        frL.SetTitle(";E  [GeV];#sigma_{RMS} Corr / #sigma_{RMS} Raw");
        frL.SetMinimum(yLo);
        frL.SetMaximum(yHi);

        // Enlarge Y-axis ticks and title on the bottom subpanel
        frL.GetYaxis()->SetLabelSize(0.055);   // tick label size
        frL.GetYaxis()->SetTitleSize(0.065);   // title font size
        frL.GetYaxis()->SetTitleOffset(0.8);  // bring title closer

        frL.Draw("AXIS");

      TLine l1(xMin,1.0,xMax,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

      for (auto& g : keepL) g->Draw("P SAME");

        // -------- Console table: ratios to baseline for all overlaid variants --------
        {
          // Build the list of variants that are actually being drawn (exclude the baseline)
          std::vector<std::string> vars;
          vars.reserve(which.size());
          for (const auto& k : which) if (k != baselineKey && S.count(k)) vars.push_back(k);

          // Sort into a stable, readable order if present
          auto pos = [&](const std::string& k)->int {
            if (k=="CLUSCP")            return 0;  // Legacy constant-b
            if (k=="CLUSCP_EA_etaDep")  return 1;  // |η|-dep
            if (k=="CLUSCP_EA_vzDep")   return 2;  // |z_vtx|-dep
            if (k=="CLUSCP_EA")         return 3;  // Energy-only
            if (k=="CLUSRAW")           return 4;  // (shouldn't be here, but keep stable)
            return 5;
          };
          std::sort(vars.begin(), vars.end(), [&](const std::string& a, const std::string& b){
            if (pos(a)!=pos(b)) return pos(a)<pos(b);
            return a<b;
          });

          // Friendly column names (reuse legend label helper defined above)
          auto colName = [&](const std::string& k)->std::string {
            return std::string(legend_label(k));
          };

          // Header
          std::cout << "\n\033[1m[RMS ratio table]  baseline = " << baselineKey
                    << "   (" << kind.tag << (vzNoteStr?vzNoteStr:"") << ")\033[0m\n";

          const int Wbin = 12;
          const int Wcol = 18;

          // First header row
          std::cout << std::left << std::setw(Wbin) << "E-bin";
          for (const auto& k : vars) {
            std::ostringstream h; h << colName(k);
            std::cout << " | " << std::setw(Wcol) << h.str();
          }
          std::cout << "\n";

          // Second header row: what each column shows
          std::cout << std::left << std::setw(Wbin) << "";
          for (size_t j=0;j<vars.size();++j) {
            std::cout << " | " << std::setw(Wcol) << "ratio (Δ%)";
          }
          std::cout << "\n";

          // Rule
          std::cout << std::string(Wbin + (int)vars.size()*(3+Wcol), '-') << "\n";

          // Per energy bin rows
          for (int i=0; i<NB; ++i) {
            // skip empty baseline points
            if (!(i<(int)B.rms.size()) || !std::isfinite(B.rms[i]) || B.rms[i]<=0.0) continue;

            // E-bin label from the configured slices
            std::ostringstream ebin; ebin.setf(std::ios::fixed);
            ebin << std::setprecision(0) << eSlices[i].first << "-" << eSlices[i].second;
            std::cout << std::left << std::setw(Wbin) << ebin.str();

            for (const auto& k : vars) {
              const auto &V = S.at(k);
              double ratio = std::numeric_limits<double>::quiet_NaN();
              if (i<(int)V.rms.size() && std::isfinite(V.rms[i]) && V.rms[i]>0.0) {
                ratio = V.rms[i] / B.rms[i];
              }
              std::ostringstream cell; cell.setf(std::ios::fixed);
              if (std::isfinite(ratio)) {
                const double dPct = (ratio - 1.0) * 100.0;
                cell << std::setprecision(4) << ratio
                     << " (" << std::showpos << std::setprecision(2) << dPct << "%)";
              } else {
                cell << "  n/a";
              }
              std::string s = cell.str();
              if ((int)s.size() < Wcol) s += std::string(Wcol - (int)s.size(), ' ');
              std::cout << " | " << s;
            }
            std::cout << "\n";
          }
          std::cout << std::endl;
        }
        // -----------------------------------------------------------------------------


        // Save the main file
        c.SaveAs(outPng.c_str());


        // --- Also emit a second version that includes the Energy-only series (CLUSCP_EA) ---
        // Only trigger for ".../sigma_vsE_ratio_to_CLUSCP_EA_vz_eta.png" and only if CLUSCP_EA
        // is not already in the 'which' list.
        {
          TString outName(outPng.c_str());
          const bool isTarget = outName.Contains("sigma_vsE_ratio_to_CLUSCP_EA_vz_eta.png");

          bool hasEonlyAlready = false;
          for (const auto& w : which) { if (w == "CLUSCP_EA") { hasEonlyAlready = true; break; } }

          if (isTarget && !hasEonlyAlready) {
            std::vector<std::string> whichPlus = which;
            whichPlus.push_back("CLUSCP_EA");  // Energy-Dep b (E-only)

            TString outPlus = outName;
            outPlus.ReplaceAll(".png", "_withEonly.png");

            // Re-enter this routine with the expanded variant set & new filename
            sigma_ratio_overlay(kind, S, whichPlus, baselineKey, outPlus.Data(), vzNoteStr);
          }
        }
    };
    
    
  // ------- helper to emit overlay groups for one kind -------
  auto save_meansigma_groups = [&](const KindDef& kind,
                                   const Var2Stats& S,
                                   const std::string& baseOut,
                                   bool isVz60)
  {
    const char* vzn = vz_note(isVz60);
    // all variants
    { std::vector<std::string> grp; for (const auto& v : kVariants) grp.push_back(v.key);
      meansigma_overlay(kind, S, grp, baseOut + "/allVariants.png", vzn);
    }
    // PDCraw vs PDCcor
    meansigma_overlay(kind, S, {"PDCraw","PDCcor"}, baseOut + "/PDCraw_vs_PDCcor.png", vzn);
    // PDCraw, PDCcor, CLUSRAW, CLUSCP
    meansigma_overlay(kind, S, {"PDCraw","PDCcor","CLUSRAW","CLUSCP"}, baseOut + "/PDC_PDC_CLUSRAW_CLUSCP.png", vzn);
    // CLUSRAW, CLUSCP, CLUSCP_EA, CLUSCP_EA_vzDep, CLUSCP_EA_etaDep
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep"},
                      baseOut + "/CLUS_five_way.png", vzn);
    // triples
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP"},                         baseOut + "/CLUSRAW_vs_CLUSCP.png", vzn);
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA_vzDep"},       baseOut + "/CLUSRAW_CLUSCP_CLUSCP_EA_vzDep.png", vzn);
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA_etaDep"},      baseOut + "/CLUSRAW_CLUSCP_CLUSCP_EA_etaDep.png", vzn);
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA"},             baseOut + "/CLUSRAW_CLUSCP_CLUSCP_EA.png", vzn);
    meansigma_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA"},
                      baseOut + "/CLUSRAW_CLUSCP_EA_all3.png", vzn);
  };

  auto save_sigma_groups = [&](const KindDef& kind,
                               const Var2Stats& S,
                               const std::string& baseOut,
                               bool isVz60)
  {
    const char* vzn = vz_note(isVz60);
    sigma_ratio_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA"}, "CLUSRAW",
                        baseOut + "/sigma_vsE_ratio_to_CLUSRAW_CP_EA.png", vzn);
    sigma_ratio_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep"}, "CLUSRAW",
                        baseOut + "/sigma_vsE_ratio_to_CLUSRAW_CP_EA_vz_eta.png", vzn);
    sigma_ratio_overlay(kind, S, {"CLUSCP","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep"}, "CLUSCP",
                        baseOut + "/sigma_vsE_ratio_to_CLUSCP_EA_vz_eta.png", vzn);
    sigma_ratio_overlay(kind, S, {"CLUSRAW","PDCcor"}, "CLUSRAW",
                        baseOut + "/sigma_vsE_ratio_PDCcor_to_CLUSRAW.png", vzn);
  };

  // Emit overlays (defaults + 0–60) with vz notes in titles
  save_meansigma_groups(kPhi, S_PH,    DIR_PH + "/Overlays/meansigma", false);
  save_meansigma_groups(kEta, S_ET,    DIR_ET + "/Overlays/meansigma", false);
  save_meansigma_groups(kPhi, S_PH_60, DIR_PH + "/Overlays_0_60vz/meansigma", true);
  save_meansigma_groups(kEta, S_ET_60, DIR_ET + "/Overlays_0_60vz/meansigma", true);

  save_sigma_groups(kPhi, S_PH,    DIR_PH + "/Overlays/sigmaRatio", false);
  save_sigma_groups(kEta, S_ET,    DIR_ET + "/Overlays/sigmaRatio", false);
  save_sigma_groups(kPhi, S_PH_60, DIR_PH + "/Overlays_0_60vz/sigmaRatio", true);
  save_sigma_groups(kEta, S_ET_60, DIR_ET + "/Overlays_0_60vz/sigmaRatio", true);

  // ------------- first-bin pair overlays (normalized shapes) -------------
  auto first_bin_pair = [&](const KindDef& kind,
                            const Var2Hists& M,
                            const std::string& outDir,
                            const std::string& A,
                            const std::string& B,
                            const std::pair<double,double>& e0,
                            bool isVz60)
  {
    ensure_dir(outDir);
    TH1* hA = (M.count(A) && M.at(A)[0]) ? M.at(A)[0] : nullptr;
    TH1* hB = (M.count(B) && M.at(B)[0]) ? M.at(B)[0] : nullptr;
    if (!hA || !hB || hA->GetEntries()<=0 || hB->GetEntries()<=0) return;

    auto* a = (TH1*)hA->Clone(); a->SetDirectory(nullptr);
    auto* b = (TH1*)hB->Clone(); b->SetDirectory(nullptr);
    if (a->Integral()!=0) a->Scale(1.0 / a->Integral());
    if (b->Integral()!=0) b->Scale(1.0 / b->Integral());

    style_h(a, color_of(A)); style_h(b, color_of(B));

    TCanvas c("cFB","",900,700);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetRightMargin(0.04);

    const double ymax = 1.2 * std::max(a->GetMaximum(), b->GetMaximum());
    a->GetYaxis()->SetRangeUser(0.0, ymax);
    a->GetXaxis()->SetTitle(kind.axisTitle.c_str());
    a->GetYaxis()->SetTitle("Normalized counts");
    a->Draw("HIST");
    b->Draw("HIST SAME");

    TLegend lg(0.60,0.75,0.92,0.90); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.040);
    lg.AddEntry(a, pretty_of(A), "l");
    lg.AddEntry(b, pretty_of(B), "l");
    lg.Draw();

    TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.042); hdr.SetTextAlign(13);
    hdr.DrawLatex(0.12,0.94, Form("First energy bin (%.0f < E < %.0f GeV%s)",
                                  e0.first, e0.second, vz_note(isVz60)));

    std::string stem = outDir + "/" + A + "_vs_" + B + ".png";
    c.SaveAs(stem.c_str());
    delete a; delete b;
  };

  // φ
  first_bin_pair(kPhi, PH,    DIR_PH + "/Overlays/firstBin", "CLUSRAW","CLUSCP",                 eSlices[0], false);
  first_bin_pair(kPhi, PH,    DIR_PH + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA",              eSlices[0], false);
  first_bin_pair(kPhi, PH,    DIR_PH + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_vzDep",        eSlices[0], false);
  first_bin_pair(kPhi, PH,    DIR_PH + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_etaDep",       eSlices[0], false);
  first_bin_pair(kPhi, PH,    DIR_PH + "/Overlays/firstBin", "PDCraw","PDCcor",                  eSlices[0], false);
  first_bin_pair(kPhi, PH_60, DIR_PH + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP",          eSlices[0], true);
  first_bin_pair(kPhi, PH_60, DIR_PH + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA",       eSlices[0], true);
  first_bin_pair(kPhi, PH_60, DIR_PH + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA_vzDep", eSlices[0], true);
  first_bin_pair(kPhi, PH_60, DIR_PH + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA_etaDep",eSlices[0], true);
  first_bin_pair(kPhi, PH_60, DIR_PH + "/Overlays_0_60vz/firstBin", "PDCraw","PDCcor",           eSlices[0], true);

  // η
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP",                 eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA",              eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_vzDep",        eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_etaDep",       eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "PDCraw","PDCcor",                  eSlices[0], false);
  first_bin_pair(kEta, ET_60, DIR_ET + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP",          eSlices[0], true);
  first_bin_pair(kEta, ET_60, DIR_ET + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA",       eSlices[0], true);
  first_bin_pair(kEta, ET_60, DIR_ET + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA_vzDep", eSlices[0], true);
  first_bin_pair(kEta, ET_60, DIR_ET + "/Overlays_0_60vz/firstBin", "CLUSRAW","CLUSCP_EA_etaDep",eSlices[0], true);
  first_bin_pair(kEta, ET_60, DIR_ET + "/Overlays_0_60vz/firstBin", "PDCraw","PDCcor",           eSlices[0], true);

  // ---------------- CSV summaries (defaults only) ----------------
  auto write_csv = [&](const KindDef& kind, const Var2Stats& S, const std::string& outCsv){
    std::ofstream fout(outCsv);
    if (!fout) { std::cerr << "[MakeResidualsSuite][WARN] Cannot write " << outCsv << "\n"; return; }
    fout << "variant,Ecenter,mean,meanErr,RMS,RMSErr\n";
    for (const auto& v : kVariants) {
      const auto it = S.find(v.key);
      if (it==S.end()) continue;
      for (int i=0;i<NB;++i){
        const double m  = it->second.mean[i];
        const double dm = it->second.dmean[i];
        const double s  = it->second.rms[i];
        const double ds = it->second.drms[i];
        if (!std::isfinite(m) || !std::isfinite(s)) continue;
        fout << v.key << "," << eCenters[i] << "," << m << "," << dm << "," << s << "," << ds << "\n";
      }
    }
    std::cout << "[MakeResidualsSuite] Wrote " << outCsv << "\n";
  };
  write_csv(kPhi, S_PH, DIR_PH + "/Overlays/meansigma/summary_phi.csv");
  write_csv(kEta, S_ET, DIR_ET + "/Overlays/meansigma/summary_eta.csv");

  // ======================================================================
  //  NEW: φ tilt fits (all variants): μ/σ vs E overlay + μ vs ln(E) fits
  // ======================================================================
  {
    const std::string outDir = DIR_PH + "/Overlays/phiTILTfits";
    ensure_dir(outDir);

    // ----- 2-panel μ/σ vs E overlay (all variants) with a columnated legend
    {
      TCanvas c("cPhiTiltMuSig","",1000,850);
      TPad top ("top" ,"",0,0.35,1,1);  top.SetBottomMargin(0.02); top.SetGrid(0,0);  top.Draw();
      TPad bot ("bot" ,"",0,0.00,1,0.35); bot.SetTopMargin(0.04);  bot.SetGrid(0,0);  bot.Draw();

      const double xMin = eSlices.front().first - 0.5;
      const double xMax = eSlices.back().second + 0.5;

      // Upper (μ)
      top.cd();
      double muLo=+1e30, muHi=-1e30;
      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;
        for (int i=0;i<NB;++i){
          const double m=it->second.mean[i], dm=it->second.dmean[i];
          if (!std::isfinite(m)) continue;
          muLo = std::min(muLo, m - (std::isfinite(dm)?dm:0.0));
          muHi = std::max(muHi, m + (std::isfinite(dm)?dm:0.0));
        }
      }
      if (!(std::isfinite(muLo)&&std::isfinite(muHi))){ muLo=-0.03; muHi=0.03; }
      const double padMu = 0.25*(muHi - muLo);
      TH2F frU("frU","",100,xMin,xMax,100, muLo - padMu, muHi + padMu);
      frU.SetTitle(Form("Residual means vs E (phi)%s;E [GeV];mean(#Delta#phi)", vz_note(false)));
      frU.Draw();

      TLegend legU(0.14,0.73,0.94,0.92); legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.034);
      legU.SetNColumns(3);

      std::vector<std::unique_ptr<TGraphErrors>> gKeepMu;
      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;
        auto g = mkGraph(eCenters, it->second.mean, it->second.dmean, color_of(v.key));
        if (g->GetN()>0){ g->Draw("P SAME"); legU.AddEntry(g.get(), pretty_of(v.key), "p"); gKeepMu.emplace_back(std::move(g)); }
      }
      legU.Draw();

      // Lower (σ)
      bot.cd();
      double sgHi = -1e30;
      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;
        for (double s : it->second.rms){ if (std::isfinite(s)) sgHi = std::max(sgHi, s); }
      }
      if (!std::isfinite(sgHi)) sgHi = 0.06;
      TH2F frL("frL","",100,xMin,xMax,100, 0.0, 1.15*sgHi);
      frL.SetTitle(";E [GeV];RMS(#Delta#phi)");
      frL.Draw();

      std::vector<std::unique_ptr<TGraphErrors>> gKeepSg;
      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;
        auto g = mkGraph(eCenters, it->second.rms, it->second.drms, color_of(v.key));
        if (g->GetN()>0){ g->Draw("P SAME"); gKeepSg.emplace_back(std::move(g)); }
      }

      c.SaveAs((outDir + "/DeltaPhi_MeanSigmaVsE_AllVariants.png").c_str());
    }

    // ----- μ vs ln(E) overlay with per-variant linear fits; export enumerated .txt
    {
      std::vector<double> lnE; lnE.reserve(eCenters.size());
      for (double e : eCenters) if (e>0) lnE.push_back(std::log(e));

      double xl=+1e30, xh=-1e30;
      for (double e : eCenters) if (e>0){ const double x = std::log(e); xl=std::min(xl,x); xh=std::max(xh,x); }
      if (!(std::isfinite(xl)&&std::isfinite(xh))){ xl=0.0; xh=1.0; }

      double yl=+1e30, yh=-1e30;
      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;
        for (double m : it->second.mean){
          if (!std::isfinite(m)) continue;
          yl = std::min(yl, m);
          yh = std::max(yh, m);
        }
      }
      if (!(std::isfinite(yl)&&std::isfinite(yh))){ yl=-0.01; yh=0.01; }
      const double pad = 0.15*(yh-yl); yl -= pad; yh += 0.35*(yh-yl);

      TCanvas c("cPhiLn","Δφ - μ vs ln(E) (ALL variants)",1100,780);
      c.SetLeftMargin(0.15); c.SetRightMargin(0.06);
      TH1F fr("fr","",1,xl-0.05,xh+0.05);
      fr.SetMinimum(yl); fr.SetMaximum(yh);
      fr.SetTitle(";ln E  [GeV];#mu  [rad]");
      fr.Draw("AXIS");

      TLegend lg(0.17,0.72,0.92,0.90);
      lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.030);
      lg.SetNColumns(2);

      std::map<std::string, std::pair<double,double>> fitBM; // key -> {b,m}
      std::vector<std::unique_ptr<TGraphErrors>> keepG;
      std::vector<std::unique_ptr<TF1>>          keepF;

      auto make_ln_graph = [&](const std::vector<double>& Ec,
                               const std::vector<double>& MU,
                               const std::vector<double>& dMU)->std::unique_ptr<TGraphErrors>
      {
        std::vector<double> x,y,xe,ye;
        for (size_t i=0;i<Ec.size();++i){
          double e = Ec[i];
          double m = (i<MU.size()? MU[i] : std::numeric_limits<double>::quiet_NaN());
          double dm= (i<dMU.size()? dMU[i]: 0.0);
          if (e>0 && std::isfinite(m)){
            x.push_back(std::log(e));
            y.push_back(m);
            xe.push_back(0.0);
            ye.push_back(std::isfinite(dm)?dm:0.0);
          }
        }
        return std::make_unique<TGraphErrors>((int)x.size(), x.data(), y.data(), xe.data(), ye.data());
      };

      for (const auto& v : kVariants){
        const auto it = S_PH.find(v.key); if (it==S_PH.end()) continue;

        auto g = make_ln_graph(eCenters, it->second.mean, it->second.dmean);
        if (g->GetN() < 2) continue; // need at least 2 points to fit

        g->SetMarkerStyle(kMkStyle);
        g->SetMarkerColor(color_of(v.key));
        g->SetLineColor  (color_of(v.key));
        g->SetMarkerSize(1.05);
        g->Draw("P SAME");

        TGraphErrors* gptr = g.get();
        keepG.emplace_back(std::move(g));

        TF1* f = new TF1(Form("f_phi_%s", v.key.c_str()), "pol1", xl, xh);
        f->SetLineColor (color_of(v.key));
        f->SetLineStyle (2);
        f->SetLineWidth (2);
        gptr->Fit(f, "Q");
        f->Draw("SAME");
        keepF.emplace_back(f);

        const double b = f->GetParameter(0); // intercept
        const double m = f->GetParameter(1); // slope
        fitBM[v.key] = {b,m};

        lg.AddEntry(gptr, Form("%s  (m=%.2e, b=%.2e)", pretty_of(v.key), m, b), "p");
      }

      lg.Draw();
      c.SaveAs((outDir + "/DeltaPhi_MuVsLogE_AllVariants.png").c_str());

      struct OutLine { const char* key; const char* enumName; const char* comment; };
      const std::vector<OutLine> order = {
        {"CLUSRAW",            "ETiltVariant::CLUS_RAW",               "no corr, cluster"},
        {"CLUSCP",             "ETiltVariant::CLUS_CP",                "CorrectPosition, cluster"},
        {"CLUSCP_EA_vzDep",    "ETiltVariant::CLUS_CP_EA_FIT_ZDEP",    "CorrectPosition(EA geom), cluster"},
        {"CLUSCP_EA_etaDep",   "ETiltVariant::CLUS_CP_EA_FIT_ETADEP",  "CorrectPosition(EA |#eta|+E), cluster"},
        {"CLUSCP_EA",          "ETiltVariant::CLUS_CP_EA_FIT_EONLY",   "CorrectPosition(EA E-only), cluster"},
        {"CLUSCP_EA_mixed",    "ETiltVariant::CLUS_CP_EA_MIX",         "CorrectPosition(EA #varphi:E-only, #eta:|#eta|+E), cluster"},
        {"CLUScpDirectBvals",  "ETiltVariant::CLUS_BCORR",             "b(E) corr, cluster"},
        {"PDCraw",             "ETiltVariant::PDC_RAW",                "no corr, scratch"},
        {"PDCcor",             "ETiltVariant::PDC_CORR",               "b(E) corr, scratch"}
      };

      std::ofstream fout(outDir + "/DeltaPhi_MuVsLogE_fit.txt", std::ios::trunc);
      fout << "# ---- Δφ: μ vs ln(E) linear fits (intercept first, slope second) ----\n";
        for (const auto& L : order){
          auto it = fitBM.find(L.key);
          if (it == fitBM.end()) continue; // skip if variant not present / not fitted
          const double b = it->second.first;
          const double m = it->second.second;

          // force first entry (intercept) to be non-negative in the output
          const double bOut = std::fabs(b);

          fout << "  {" << L.enumName << ", { "
               << Form("%.6e", bOut) << ", " << Form("%.6e", m)
               << " }},  // " << L.comment << "\n";
        }
      fout.close();
    }
  } // end φ tilt fits wrapper

  // ---- Min-RMS summary tables (console + CSVs)
  {
    const std::string summaryDir = DIR_PH + std::string("/Overlays/phiTILTfits");
    EmitMinRmsSummary(summaryDir, NB, eSlices, S_PH, S_ET);
  }

  std::cout << "[MakeResidualsSuite] Residuals suite written under: " << RESID << "\n";
}






















static inline std::string SanitizeForPath(TString s)
{
  std::string out = s.Data();
  for (char& c : out)
  {
    if (c==' ' || c=='/' || c=='\\' || c=='(' || c==')' ||
        c=='<' || c=='>' || c=='=' || c==':' || c==','  ||
        c=='[' || c==']' || c=='|' || c=='\"' || c=='\'')
      c = '_';
  }
  return out;
}

// Try to extract a compact energy-slice label from the histogram name.
// Falls back to the full name if no simple pattern is found.
static inline std::string ExtractEnergyLabel(const TString& hname)
{
  const std::string s = hname.Data();
  // common suffix patterns people use:  _EbinNN  /  _E_NN  /  _E NN
  int pos = -1;
  if ((pos = s.rfind("_Ebin")) != (int)std::string::npos) return s.substr(pos+1);
  if ((pos = s.rfind("_E_"))   != (int)std::string::npos) return s.substr(pos+1);
  if ((pos = s.rfind("_E"))    != (int)std::string::npos) return s.substr(pos+1);
  return s;
}

void AnalyzeBvsResCLUSCPEA(const char* inFilePath, const char* outBaseDir)
{
  if (!inFilePath || !outBaseDir) { std::cerr << "[AnalyzeBvsRes] bad args\n"; return; }

  // ---- output dirs
  const std::string base   = std::string(outBaseDir) + "/analyzeBvsResCLUSCPEA";
  const std::string dPhi2D = base + "/phi/2D";
  const std::string dPhiPr = base + "/phi/profile";
  const std::string dPhi1D = base + "/phi/projections";
  const std::string dEta2D = base + "/eta/2D";
  const std::string dEtaPr = base + "/eta/profile";
  const std::string dEta1D = base + "/eta/projections";

  EnsureDir(base);  EnsureDir(base + "/phi"); EnsureDir(base + "/eta");
  EnsureDir(dPhi2D); EnsureDir(dPhiPr); EnsureDir(dPhi1D);
  EnsureDir(dEta2D); EnsureDir(dEtaPr); EnsureDir(dEta1D);

  // ---- open file
  std::unique_ptr<TFile> fin(TFile::Open(inFilePath,"READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[AnalyzeBvsRes] cannot open " << inFilePath << "\n";
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  // 2D heatmap
  auto draw2D = [&](TH2* h2, const std::string& outPng, bool isPhi)
  {
    if (!h2) return;
    // Style + titles
    std::unique_ptr<TH2> h( static_cast<TH2*>(h2->Clone(Form("c2_%s", h2->GetName()))) );
    h->SetDirectory(nullptr);
    h->SetStats(0);
    h->SetContour(255);
    h->GetXaxis()->SetTitle(isPhi ? "b_{#varphi}" : "b_{#eta}");
    h->GetYaxis()->SetTitle(isPhi ? "#Delta#varphi (folded) [rad]" : "#Delta#eta");

    TCanvas c("c2","",1000,800);
    c.cd();
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin (0.12);
    gPad->SetBottomMargin(0.12);
    h->Draw("COLZ");
    c.SaveAs(outPng.c_str());
  };

  // Profile + RMS band
  auto drawProfilePanel = [&](TH2* h2, const std::string& outPng, bool isPhi)
  {
    if (!h2) return;

    std::unique_ptr<TProfile> prof( h2->ProfileX(Form("prof_%s", h2->GetName()), 1, -1, "S") );
    if (!prof) return;
    prof->SetLineColor(kBlue+2);
    prof->SetMarkerColor(kBlue+2);
    prof->SetMarkerStyle(20);
    prof->SetMarkerSize(0.9);
    prof->SetTitle("");

    const int nx = h2->GetNbinsX();
    std::unique_ptr<TGraphAsymmErrors> band(new TGraphAsymmErrors(nx));
    band->SetFillColorAlpha(kBlue, 0.25);
    band->SetLineColor(kBlue);
    band->SetFillStyle(1001);

    for (int ix=1; ix<=nx; ++ix)
    {
      std::unique_ptr<TH1D> slice( h2->ProjectionY(Form("slice_%s_%d", h2->GetName(), ix), ix, ix) );
      const double x   = h2->GetXaxis()->GetBinCenter(ix);
      double ymean = 0.0, yrms = 0.0;

      if (slice && slice->GetEntries() > 0)
      {
        ymean = slice->GetMean();
        yrms  = slice->GetRMS();
      }
      band->SetPoint(ix-1, x, ymean);
      band->SetPointError(ix-1, 0.0, 0.0, yrms, yrms);
    }

    double xLo = h2->GetXaxis()->GetXmin();
    double xHi = h2->GetXaxis()->GetXmax();
    double yLo = h2->GetYaxis()->GetXmin();
    double yHi = h2->GetYaxis()->GetXmax();
    const double yPad = 0.08*(yHi - yLo);
    yLo -= yPad; yHi += yPad;

    TCanvas c("cProf","",1000,800);
    c.cd();
    gPad->SetRightMargin(0.06);
    gPad->SetLeftMargin (0.12);
    gPad->SetBottomMargin(0.12);

    TH2F frame("frame","", 100, xLo, xHi, 100, yLo, yHi);
    frame.SetStats(0);
    frame.GetXaxis()->SetTitle(isPhi ? "b_{#varphi}" : "b_{#eta}");
    frame.GetYaxis()->SetTitle(isPhi ? "#Delta#varphi (folded) [rad]" : "#Delta#eta");
    frame.Draw();

    TLine zline(xLo, 0.0, xHi, 0.0);
    zline.SetLineColor(kGray+2);
    zline.SetLineStyle(2);
    zline.Draw();

    band->Draw("E3 SAME");
    prof->Draw("P SAME");

    TLegend leg(0.16,0.80,0.50,0.90);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);
    leg.AddEntry(prof.get(),  "<#Delta> vs b  (profile)", "lp");
    leg.AddEntry(band.get(),  "RMS envelope",            "f");
    leg.Draw();

    c.SaveAs(outPng.c_str());
  };

  // 1D projections
  auto draw1DProjections = [&](TH2* h2, const std::string& outDir, bool isPhi)
  {
    if (!h2) return;

    std::unique_ptr<TH1D> hX( h2->ProjectionX(Form("prjX_%s",h2->GetName()), 1, -1) );
    std::unique_ptr<TH1D> hY( h2->ProjectionY(Form("prjY_%s",h2->GetName()), 1, -1) );

    if (hX)
    {
      TCanvas cX("cX","",900,650);
      cX.cd(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      hX->SetLineColor(kBlack);
      hX->SetMarkerStyle(20);
      hX->SetMarkerColor(kBlack);
      hX->SetTitle("");
      hX->GetXaxis()->SetTitle(isPhi ? "b_{#varphi}" : "b_{#eta}");
      hX->GetYaxis()->SetTitle("entries");
      hX->Draw("E1");
      const std::string nm = outDir + "/" + SanitizeForPath(h2->GetName()) + "_projX_b.png";
      cX.SaveAs(nm.c_str());
    }

    if (hY)
    {
      TCanvas cY("cY","",900,650);
      cY.cd(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      hY->SetLineColor(kBlack);
      hY->SetMarkerStyle(20);
      hY->SetMarkerColor(kBlack);
      hY->SetTitle("");
      hY->GetXaxis()->SetTitle(isPhi ? "#Delta#varphi (folded) [rad]" : "#Delta#eta");
      hY->GetYaxis()->SetTitle("entries");
      hY->Draw("E1");
      const std::string nm = outDir + "/" + SanitizeForPath(h2->GetName()) + "_projY_resid.png";
      cY.SaveAs(nm.c_str());
    }
  };

  // NEW: per-b dissection table (printed to terminal)
  auto printPerBinTable = [&](TH2* h2, bool isPhi)
  {
    if (!h2) return;

    const int nx = h2->GetNbinsX();
    const int ny = h2->GetNbinsY();

    // Totals
    Long64_t tot = static_cast<Long64_t>(h2->GetEntries());
    Long64_t nonEmptyXBins = 0;
    for (int ix=1; ix<=nx; ++ix)
    {
      double sumW = 0;
      for (int iy=1; iy<=ny; ++iy) sumW += h2->GetBinContent(ix, iy);
      if (sumW > 0) ++nonEmptyXBins;
    }

    const TString hname = h2->GetName();
    const std::string eLbl = ExtractEnergyLabel(hname);
    const std::string side = isPhi ? "Δφ vs bφ" : "Δη vs bη";

    std::cout << "\n";
    std::cout << "────────────────────────────────────────────────────────────────────────────\n";
    std::cout << "  " << side << "  |  hist: " << hname
              << "  |  energy-slice: " << eLbl << "\n";
    std::cout << "  X bins = " << nx << "  |  non-empty X bins = " << nonEmptyXBins
              << "  |  total entries = " << tot << "\n";
    std::cout << "────────────────────────────────────────────────────────────────────────────\n";

    // Header
    std::cout << std::left
              << std::setw(6)  << "ix"
              << std::setw(14) << "b_center"
              << std::setw(12) << "entries"
              << std::setw(16) << "<resid>"
              << std::setw(16) << "RMS(resid)"
              << std::setw(16) << "SE(<resid>)"
              << "\n";
    std::cout << std::string(6+14+12+16+16+16, '-') << "\n";

    // Rows
    for (int ix=1; ix<=nx; ++ix)
    {
      std::unique_ptr<TH1D> slice(
        h2->ProjectionY(Form("slice_%s_%d_tab", h2->GetName(), ix), ix, ix)
      );

      const double bCenter = h2->GetXaxis()->GetBinCenter(ix);
      Long64_t n = 0;
      double mean = 0.0, rms = 0.0, se = 0.0;

      if (slice)
      {
        n    = static_cast<Long64_t>(slice->GetEntries());
        mean = (n > 0) ? slice->GetMean() : 0.0;
        rms  = (n > 1) ? slice->GetRMS()  : 0.0;
        se   = (n > 0) ? (rms / std::sqrt(static_cast<double>(n))) : 0.0;
      }

      std::cout << std::left
                << std::setw(6)  << ix
                << std::setw(14) << Form("%.6f", bCenter)
                << std::setw(12) << n
                << std::setw(16) << Form("%+.6e", mean)
                << std::setw(16) << Form("%.6e",  rms)
                << std::setw(16) << Form("%.6e",  se)
                << "\n";
    }

    std::cout << std::string(6+14+12+16+16+16, '-') << "\n";
    std::cout << std::flush;
  };

  // ---- iterate keys and process all matching TH2
  TIter nextKey(fin->GetListOfKeys());
  while (TKey* k = static_cast<TKey*>(nextKey()))
  {
    TObject* obj = k->ReadObj();
    if (!obj) continue;

    TH2* h2 = dynamic_cast<TH2*>(obj);
    if (!h2) { delete obj; continue; }

    const TString hname = h2->GetName();
    const bool isPhi = hname.BeginsWith("h2_dphiEA_vs_bphi_");
    const bool isEta = hname.BeginsWith("h2_detaEA_vs_beta_");
    if (!isPhi && !isEta) { delete obj; continue; }

    // Detach from file so we can safely delete the original object
    h2->SetDirectory(nullptr);

    // Output paths
    const std::string out2D = (isPhi ? dPhi2D : dEta2D) + "/" + SanitizeForPath(hname) + ".png";
    const std::string outPR = (isPhi ? dPhiPr : dEtaPr) + "/" + SanitizeForPath(hname) + "_profile.png";
    const std::string out1D = (isPhi ? dPhi1D : dEta1D);

    // Plots
    draw2D(h2, out2D, isPhi);
    drawProfilePanel(h2, outPR, isPhi);
    draw1DProjections(h2, out1D, isPhi);

    // NEW: per-b dissection (printed to terminal)
    printPerBinTable(h2, isPhi);

    delete obj; // original read object (we cloned or projected as needed)
  }

  std::cout << "[AnalyzeBvsRes] PNGs written under: " << base << "\n";
}




// ==================== ASH/LOG RAW‑RMS scan reader & plots ====================
//
// This block reads the per‑parameter scan histograms produced by your macro:
//   • ASH:  h_dx_ash_b{bVal4dec}_{E‑label}
//   • LOG:  h_dx_log_w0{w0Val2dec}_{E‑label}
// and produces publication‑style overlays using **RAW RMS** only (no Gaussians).
//
// All outputs are written under:  <outBaseDir>/ashLogScans
// It is self‑contained; nothing else in your file is modified.
//

// Small palette for E‑bin overlays
static const Color_t kAshPaletteRMS[5] =
{ kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7 };

/* ---------- helpers to format labels and find histograms ---------- */
static inline TString fmtElabelRange(int iE, const double* edges)
{
  // Example: E_2to4, E_4to6, ...
  const int lo = static_cast<int>(std::lround(edges[iE]));
  const int hi = static_cast<int>(std::lround(edges[iE+1]));
  return Form("E_%dto%d", lo, hi);
}

static inline TString fmtElabelEbin(int iE)
{
  // Example: Ebin0, Ebin1, ...
  return Form("Ebin%d", iE);
}

static inline TString fmtB4(double b)  { return Form("%.4f", std::round(b*10000.)/10000.); }
static inline TString fmtW2(double w0) { return Form("%.2f",  std::round(w0*100.)/100.);   }

static TH1* tryGet1D(TFile* f, const std::vector<TString>& names)
{
  for (const auto& n : names)
  {
    if (TH1* h = dynamic_cast<TH1*>( f->Get(n) ))
    {
      if (h->GetEntries() > 0 && h->Integral() > 0) return h;
    }
  }
  return nullptr;
}

/* ----------------------- RAW‑RMS metric (no fits) ----------------------- */
double rawRMS(TH1* h, const TString& /*debugPNG*/)
{
  if (!h) return 0.0;
  if (h->GetEntries() < 2) return 0.0;
  return h->GetRMS();
}

/* -------------- containers for scan results (one TGraph per E bin) ------ */
struct ScanResultsRMS
{
  std::vector<TGraph*> tgVec;      // one TGraph per E bin
  std::vector<double>  bestParam;  // best b or best w0 per E bin
  std::vector<double>  bestSigma;  // minimal RMS per E bin
  double minY = DBL_MAX, maxY = -DBL_MAX;
  int totalHist=0, missingHist=0, zeroEntry=0, usedHist=0;
};

/* --------------------------- do ASH (RAW‑RMS) --------------------------- */
static ScanResultsRMS
doAshScan_RMS(TFile* f,
              int N_E,
              const double* E_edges,
              const std::vector<double>& bScan_vals,
              const TString& suffixLabel,
              const TString& histOutDir,
              bool skipLowest3)
{
  const int iE_start = skipLowest3 ? 3 : 0;

  ScanResultsRMS R;
  R.tgVec.resize(N_E);

  for (int iE = iE_start; iE < N_E; ++iE)
  {
    TGraph* g = new TGraph;
    g->SetName(Form("gAsh_%s_E%d", suffixLabel.Data(), iE));

    for (double bValRaw : bScan_vals)
    {
      R.totalHist++;

      const TString bStr = fmtB4(bValRaw);
      const TString el1  = fmtElabelRange(iE, E_edges);
      const TString el2  = fmtElabelEbin(iE);

      // Try a few common naming patterns used across macros
      std::vector<TString> cands = {
        Form("h_dx_ash_b%s_%s", bStr.Data(), el1.Data()),
        Form("h_dx_ash_b%s_%s", bStr.Data(), el2.Data()),
        Form("h_dx_ash_b%s_E%d", bStr.Data(), iE),
        Form("h_dx_ash_b%s",    bStr.Data()) // last‑ditch
      };

      TH1* h = tryGet1D(f, cands);
      if (!h) { R.missingHist++; continue; }
      if (h->GetEntries() <= 0 || h->Integral() <= 0) { R.zeroEntry++; continue; }
      R.usedHist++;

      const double sVal = rawRMS(h, "");
      if (sVal > 0.0)
      {
        g->SetPoint(g->GetN(), bValRaw, sVal);
        R.minY = std::min(R.minY, sVal);
        R.maxY = std::max(R.maxY, sVal);

        if ((int)R.bestSigma.size() < N_E) { R.bestSigma.resize(N_E, DBL_MAX); R.bestParam.resize(N_E, 0.0); }
        if (sVal < R.bestSigma[iE]) { R.bestSigma[iE] = sVal; R.bestParam[iE] = bValRaw; }
      }
    }

    const int idx = (iE - iE_start) % 5;
    const int col = kAshPaletteRMS[idx];
    g->SetMarkerStyle(20 + idx);
    g->SetMarkerColor(col);
    g->SetLineColor  (col);
    g->SetMarkerSize(1.2);

    R.tgVec[iE] = g;
  }

  if (R.minY == DBL_MAX) { R.minY = 0.0; R.maxY = 1.0; }
  return R;
}

/* --------------------------- do LOG (RAW‑RMS) --------------------------- */
static ScanResultsRMS
doLogScan_RMS(TFile* f,
              int N_E,
              const double* E_edges,
              const std::vector<double>& w0Scan_vals,
              const TString& suffixLabel,
              const TString& histOutDir,
              bool skipLowest3)
{
  const int iE_start = skipLowest3 ? 3 : 0;

  ScanResultsRMS R;
  R.tgVec.resize(N_E);

  for (int iE = iE_start; iE < N_E; ++iE)
  {
    TGraph* g = new TGraph;
    g->SetName(Form("gLog_%s_E%d", suffixLabel.Data(), iE));

    for (double w0Raw : w0Scan_vals)
    {
      R.totalHist++;

      const TString wStr = fmtW2(w0Raw);
      const TString el1  = fmtElabelRange(iE, E_edges);
      const TString el2  = fmtElabelEbin(iE);

      std::vector<TString> cands = {
        Form("h_dx_log_w0%s_%s", wStr.Data(), el1.Data()),
        Form("h_dx_log_w0%s_%s", wStr.Data(), el2.Data()),
        Form("h_dx_log_w0%s_E%d", wStr.Data(), iE),
        Form("h_dx_log_w0%s",    wStr.Data())
      };

      TH1* h = tryGet1D(f, cands);
      if (!h) { R.missingHist++; continue; }
      if (h->GetEntries() <= 0 || h->Integral() <= 0) { R.zeroEntry++; continue; }
      R.usedHist++;

      const double sVal = rawRMS(h, "");
      if (sVal > 0.0)
      {
        g->SetPoint(g->GetN(), w0Raw, sVal);
        R.minY = std::min(R.minY, sVal);
        R.maxY = std::max(R.maxY, sVal);

        if ((int)R.bestSigma.size() < N_E) { R.bestSigma.resize(N_E, DBL_MAX); R.bestParam.resize(N_E, 0.0); }
        if (sVal < R.bestSigma[iE]) { R.bestSigma[iE] = sVal; R.bestParam[iE] = w0Raw; }
      }
    }

    const int idx = (iE - iE_start) % 5;
    const int col = kAshPaletteRMS[idx];
    g->SetMarkerStyle(20 + idx);
    g->SetMarkerColor(col);
    g->SetLineColor  (col);
    g->SetMarkerSize(1.2);

    R.tgVec[iE] = g;
  }

  if (R.minY == DBL_MAX) { R.minY = 0.0; R.maxY = 1.0; }
  return R;
}

/* --------------------- side‑by‑side overlay (RAW‑RMS) --------------------- */
static void drawAshLogSideBySide_RMS( const ScanResultsRMS& ashRes,
                                      const ScanResultsRMS& logRes,
                                      const double* E_edges,
                                      int N_E,
                                      const TString& outBaseDir,
                                      bool skipLowest3 )
{
  const int iE_start = skipLowest3 ? 3 : 0;

  // Canvas with two pads
  TCanvas cSide("cSideBySide_RMS","ASH vs LOG - RAW RMS",1600,600);
  cSide.Divide(2,1);

  // ---- Left: ASH ----
  cSide.cd(1);
  gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetGrid();

  TGraph* gAshRef = nullptr;
  for (int iE=iE_start;iE<N_E && !gAshRef;++iE) if (ashRes.tgVec[iE] && ashRes.tgVec[iE]->GetN()) gAshRef = ashRes.tgVec[iE];
  if (gAshRef)
  {
    gAshRef->Draw("ALP");
    gAshRef->GetXaxis()->SetTitle("b  (local tower units)");
    gAshRef->GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (local tower units)");
    gAshRef->GetYaxis()->SetRangeUser(std::max(0.0, ashRes.minY*0.90), ashRes.maxY*1.10);
    for (int iE=iE_start;iE<N_E;++iE) if (ashRes.tgVec[iE] && ashRes.tgVec[iE]->GetN()) ashRes.tgVec[iE]->Draw("LP SAME");
  }

  TLegend LA(0.13,0.65,0.48,0.88); LA.SetBorderSize(0); LA.SetFillStyle(0);
  for (int iE=iE_start;iE<N_E;++iE)
  {
    TGraph* g = ashRes.tgVec[iE];
    const double bBest = (iE < (int)ashRes.bestParam.size()) ? ashRes.bestParam[iE] : 0.0;
    const double sBest = (iE < (int)ashRes.bestSigma.size()) ? ashRes.bestSigma[iE] : 0.0;
    if (!g || !g->GetN())
      LA.AddEntry((TObject*)nullptr, Form("%.0f #leq E < %.0f GeV - no data", E_edges[iE], E_edges[iE+1]), "");
    else
      LA.AddEntry(g, Form("%.0f #leq E < %.0f GeV  (b=%.2f,  #sigma_{RMS}=%.3f)",
                          E_edges[iE], E_edges[iE+1], bBest, sBest), "lp");
  }
  LA.Draw();
  TLatex tA; tA.SetNDC(); tA.SetTextSize(0.035); tA.DrawLatex(0.12,0.93,"Ash scan [RMS]");

  // ---- Right: LOG ----
  cSide.cd(2);
  gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetGrid();

  TGraph* gLogRef = nullptr;
  for (int iE=iE_start;iE<N_E && !gLogRef;++iE) if (logRes.tgVec[iE] && logRes.tgVec[iE]->GetN()) gLogRef = logRes.tgVec[iE];
  if (gLogRef)
  {
    gLogRef->Draw("ALP");
    gLogRef->GetXaxis()->SetTitle("w_{0}");
    gLogRef->GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (local tower units)");
    gLogRef->GetYaxis()->SetRangeUser(std::max(0.0, logRes.minY*0.90), logRes.maxY*1.10);
    for (int iE=iE_start;iE<N_E;++iE) if (logRes.tgVec[iE] && logRes.tgVec[iE]->GetN()) logRes.tgVec[iE]->Draw("LP SAME");
  }

  TLegend LB(0.60,0.65,0.85,0.88); LB.SetBorderSize(0); LB.SetFillStyle(0);
  for (int iE=iE_start;iE<N_E;++iE)
  {
    TGraph* g = logRes.tgVec[iE];
    const double wBest = (iE < (int)logRes.bestParam.size()) ? logRes.bestParam[iE] : 0.0;
    const double sBest = (iE < (int)logRes.bestSigma.size()) ? logRes.bestSigma[iE] : 0.0;
    if (!g || !g->GetN())
      LB.AddEntry((TObject*)nullptr, Form("%.0f #leq E < %.0f GeV - no data", E_edges[iE], E_edges[iE+1]), "");
    else
      LB.AddEntry(g, Form("%.0f #leq E < %.0f GeV  (w_{0}=%.2f,  #sigma_{RMS}=%.3f)",
                          E_edges[iE], E_edges[iE+1], wBest, sBest), "lp");
  }
  LB.Draw();
  TLatex tB; tB.SetNDC(); tB.SetTextSize(0.040); tB.DrawLatex(0.15,0.93,"Log scan [RMS]");

  TString out = outBaseDir + "/ashLogScans/SideBySide_RMS.png";
  gSystem->mkdir(outBaseDir + "/ashLogScans", true);
  cSide.SaveAs(out);
  std::cout << "[ASH/LOG‑RMS] wrote " << out << "\n";

  // Single‑panel publication style (ASH then LOG)
  auto makeSingle = [&](const ScanResultsRMS& R, const char* tag, const char* xTit, const char* yTit, double xMaxCut)
  {
    TGraph* gRef = nullptr;
    for (int iE=iE_start;iE<N_E && !gRef;++iE) if (R.tgVec[iE] && R.tgVec[iE]->GetN()) gRef = R.tgVec[iE];
    if (!gRef) return;

    TCanvas c(Form("c%sScan_RMS",tag), Form("%s scan (RMS)",tag), 900,650);
    c.SetLeftMargin(0.13); c.SetBottomMargin(0.18);
    gRef->Draw("ALP");
    gRef->GetXaxis()->SetTitle(xTit);
    gRef->GetYaxis()->SetTitle(yTit);

    double yMin=1e30,yMax=-1e30,x,y;
    for (int i=iE_start;i<N_E;++i)
      if (R.tgVec[i] && R.tgVec[i]->GetN())
        for (int j=0;j<R.tgVec[i]->GetN();++j)
        { R.tgVec[i]->GetPoint(j,x,y);
          if (x<=xMaxCut){ yMin=std::min(yMin,y); yMax=std::max(yMax,y);} }
    const double pad = 0.05*(yMax-yMin);
    gRef->GetYaxis()->SetRangeUser(std::max(0.0,yMin-pad), yMax+pad);

    for (int iE=iE_start;iE<N_E;++iE) if (R.tgVec[iE] && R.tgVec[iE]->GetN()) R.tgVec[iE]->Draw("LP SAME");


    TLegend L(0.63, 0.70, 0.90, 0.92);
    L.SetBorderSize(0);  L.SetFillStyle(0);  L.SetTextSize(0.022);
    for (int iE=iE_start;iE<N_E;++iE) if (R.tgVec[iE] && R.tgVec[iE]->GetN())
        L.AddEntry(R.tgVec[iE], Form("%.0f #leq E < %.0f GeV", E_edges[iE], E_edges[iE+1]), "lp");

    L.Draw();

    TString out1 = outBaseDir + "/ashLogScans/" + TString::Format("%sScan_RMS.png",tag);
    c.SaveAs(out1);
    std::cout << "[ASH/LOG‑RMS] wrote " << out1 << "\n";
  };

  makeSingle(ashRes, "Ash", "b  (local tower units)", "#sigma_{RMS}_{x}  (local tower units)", 1e9 /* no cut */);
  makeSingle(logRes, "Log", "w_{0}",                  "#sigma_{RMS}_{x}  (local tower units)",  5.5);
}

/* ---- Min‑RMS b(E) overlay vs “no |η| dep” best‑fit b(E) from bValues.txt ---- */
static void drawMinBOverlay_RMS_vs_OriginalEta( const ScanResultsRMS& ashRes,
                                                const double* E_edges, int N_E,
                                                const TString& outBaseDir,
                                                bool skipLowest3 )
{
  const int iE_start = skipLowest3 ? 3 : 0;

  // Read reference b(E) for originalEta from the existing text file
  std::vector<double> bPhiOrig(N_E, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> bEtaOrig(N_E, std::numeric_limits<double>::quiet_NaN());

  const std::string bfile = std::string(outBaseDir) + "/bValues.txt";
  std::ifstream fin(bfile);
  if (fin.is_open())
  {
    std::string line;
    while (std::getline(fin, line))
    {
      if (line.empty() || line[0]=='#') continue;
      char axis[8] = {0}, key[128] = {0};
      double eLo=0, eHi=0, bval=0;
      int n = std::sscanf(line.c_str(), " %7s [%lf , %lf ) %lf %127s", axis, &eLo, &eHi, &bval, key);
      if (n != 5) continue;
      if (std::string(key) != "originalEta") continue;

      int idx = -1;
      for (int i=0;i<N_E;++i)
        if (std::fabs(E_edges[i]-eLo)<1e-6 && std::fabs(E_edges[i+1]-eHi)<1e-6) { idx = i; break; }
      if (idx<0) continue;

      if      (std::string(axis) == "PHI") bPhiOrig[idx] = bval;
      else if (std::string(axis) == "ETA") bEtaOrig[idx] = bval;
    }
  }
  else
  {
    std::cerr << "[ASH/LOG‑RMS] WARNING: cannot open " << bfile << " - overlays will show RMS only.\n";
  }

  auto makeOverlay = [&](const char* tag, const std::vector<double>& bFitOrig)
  {
    std::vector<double> eC, bMin, bRef;
    for (int i=iE_start;i<N_E;++i)
    {
      const double eCtr = 0.5*(E_edges[i] + E_edges[i+1]);
      const bool haveMin = (i < (int)ashRes.bestParam.size()) && std::isfinite(ashRes.bestParam[i]);
      const bool haveRef = (i < (int)bFitOrig.size())        && std::isfinite(bFitOrig[i]);
      if (!haveMin && !haveRef) continue;
      eC.push_back(eCtr);
      bMin.push_back(haveMin ? ashRes.bestParam[i] : std::numeric_limits<double>::quiet_NaN());
      bRef.push_back(haveRef ? bFitOrig[i]         : std::numeric_limits<double>::quiet_NaN());
    }
    if (eC.empty()) return;

    // Build graphs
    TGraph gRMS((int)eC.size());
    TGraph gFit((int)eC.size());
    int ip=0, jf=0;
    for (std::size_t k=0;k<eC.size();++k)
    {
      if (std::isfinite(bMin[k])) { gRMS.SetPoint(ip++, eC[k], bMin[k]); }
      if (std::isfinite(bRef[k])) { gFit.SetPoint(jf++, eC[k], bRef[k]); }
    }
    gRMS.SetMarkerStyle(24); gRMS.SetMarkerSize(1.2); gRMS.SetMarkerColor(kBlue+1); gRMS.SetLineColor(kBlue+1);
    gFit.SetMarkerStyle(20); gFit.SetMarkerSize(1.2); gFit.SetMarkerColor(kBlack);  gFit.SetLineColor(kBlack);

    // Frame
    double xMin = E_edges[iE_start]-0.5, xMax = E_edges[N_E]-0.5;
    double yMin = +1e30, yMax = -1e30;
    for (int i=0;i<gRMS.GetN();++i){ double x,y; gRMS.GetPoint(i,x,y); yMin=std::min(yMin,y); yMax=std::max(yMax,y); }
    for (int i=0;i<gFit.GetN();++i){ double x,y; gFit.GetPoint(i,x,y); yMin=std::min(yMin,y); yMax=std::max(yMax,y); }
    if (!std::isfinite(yMin) || !std::isfinite(yMax)) { yMin=0.0; yMax=1.0; }
    const double pad = 0.10*(yMax-yMin);
    yMin = std::max(0.0, yMin - pad);
    yMax += pad;

    TCanvas c(Form("cMinB_%s",tag), Form("Min‑RMS b vs E - %s",tag), 900,650);
    c.SetLeftMargin(0.13); c.SetBottomMargin(0.15);
    TH2F fr("fr","",100,xMin,xMax,100,yMin,yMax);
    fr.SetStats(0);
    fr.GetXaxis()->SetTitle("E  [GeV]");
    fr.GetYaxis()->SetTitle("b  (local tower units)");
    fr.Draw();

    gRMS.Draw("P SAME");
    if (gFit.GetN()>0) gFit.Draw("P SAME");

    TLegend L(0.52,0.74,0.88,0.89);
    L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextSize(0.033);
    L.AddEntry(&gRMS, "Min RMS b (Ash, RAW RMS)", "p");
    if (gFit.GetN()>0)
      L.AddEntry(&gFit, "no |#eta| dep fit  b(E)", "p");
    L.Draw();

    gSystem->mkdir(outBaseDir + "/ashLogScans", true);
    TString outP = outBaseDir + "/ashLogScans/MinRMS_b_vs_E_OVERLAY_";
    outP += tag; outP += ".png";
    c.SaveAs(outP);
    std::cout << "[ASH/LOG‑RMS] wrote " << outP << "\n";
  };

  makeOverlay("phi", bPhiOrig);
  makeOverlay("eta", bEtaOrig);
}

/* ------------ PHENIX‑style σx(E) summary from ASH minima (RAW‑RMS) --------- */
static void makePhenixLikeSummary_RMS( const ScanResultsRMS& ashRes,
                                       const double* E_edges, int N_E,
                                       double cellSize_cm,
                                       const TString& outBaseDir,
                                       bool skipLowest3 )
{
  const int iE_start = skipLowest3 ? 3 : 0;

  TGraph gAsh; gAsh.SetName("gAshBest_RMS");
  int nGood = 0; double xMin=1e30, xMax=-1e30;

  for (int iE=iE_start;iE<N_E;++iE)
  {
    const bool has = (iE < (int)ashRes.bestSigma.size() && std::isfinite(ashRes.bestSigma[iE]));
    if (!has) continue;
    const double eCtr   = 0.5*(E_edges[iE]+E_edges[iE+1]);
    const double sigma_mm = ashRes.bestSigma[iE] * cellSize_cm * 10.0; // tower→cm→mm
    gAsh.SetPoint(gAsh.GetN(), eCtr, sigma_mm);
    ++nGood; xMin=std::min(xMin,eCtr); xMax=std::max(xMax,eCtr);
  }
  if (nGood < 3) { std::cerr << "[ASH/LOG‑RMS] PHENIX summary needs ≥3 points.\n"; return; }

  TCanvas c("cPhenixAshOnly_RMS","PHENIX σx(E) - Ash (RMS)",600,600);
  c.SetLeftMargin(0.17); c.SetBottomMargin(0.15);
  gAsh.SetMarkerStyle(24); gAsh.SetMarkerSize(1.4);
  gAsh.SetMarkerColor(kBlack); gAsh.SetLineColor(kBlack);

  gAsh.Draw("AP");
  gAsh.GetXaxis()->SetTitle("E  (GeV)");
  gAsh.GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (mm)");
  gAsh.GetXaxis()->SetLimits(xMin*0.9, xMax*1.1);

  double yMin=1e30,yMax=-1e30,x,y;
  for (int i=0;i<gAsh.GetN();++i){ gAsh.GetPoint(i,x,y); yMin=std::min(yMin,y); yMax=std::max(yMax,y); }
  const double pad = 0.15*(yMax-yMin);
  gAsh.GetYaxis()->SetRangeUser(std::max(0.0,yMin-pad), yMax+pad);

  // Optional illustrative fit p0 + p1/sqrt(E) to the **RMS** trend
  TF1 fFit("fFit_RMS","[0]+[1]/TMath::Sqrt(x)", xMin*0.9 , xMax*1.1);
  fFit.SetParameters(1.0,6.0);
  TFitResultPtr fr = gAsh.Fit(&fFit,"QRNS");
  if (fr && fr->IsValid()) { fFit.SetLineStyle(2); fFit.Draw("SAME"); }

  if (fr && fr->IsValid())
  {
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.037);
    tl.DrawLatex(0.18,0.25, Form("#sigma_{RMS}_{x}=%.2f + %.2f E^{-1/2}  (mm)",
                                 fFit.GetParameter(0), fFit.GetParameter(1)));
  }

  gSystem->mkdir(outBaseDir + "/ashLogScans", true);
  TString out = outBaseDir + "/ashLogScans/PhenixSummary_AshOnly_RMS.png";
  c.SaveAs(out);
  std::cout << "[ASH/LOG‑RMS] wrote " << out << "\n";
}

/* ---------------- orchestrator (RMS‑only) called from MAIN ---------------- */
static void runAshLogRMSOnly(const char* inFilePath, const char* outBaseDir)
{
  if (!inFilePath || !outBaseDir) return;

  // Build scans (match booking granularity)
  std::vector<double> bScan;
  for (double b = 0.01; b <= 0.50 + 1e-9; b += 0.01)
    bScan.push_back( std::round(b*10000.)/10000. );

  std::vector<double> w0Scan;
  for (double w = 1.5; w <= 7.0 + 1e-9; w += 0.1)
    w0Scan.push_back( std::round(w*100.)/100. );

  std::unique_ptr<TFile> fIn( TFile::Open(inFilePath,"READ") );
  if (!fIn || fIn->IsZombie())
  {
    std::cerr << "[ASH/LOG‑RMS] Cannot open '" << inFilePath << "'. Skip RMS scans.\n";
    return;
  }

  const TString baseDir = TString(outBaseDir) + "/ashLogScans";
  gSystem->mkdir(baseDir, true);
  const TString debugDir = baseDir + "/DEBUG"; gSystem->mkdir(debugDir, true);

  // Scans (skip the lowest 3 E bins as in your example; set false if you want all)
  const bool skipLowest3 = true;
  auto ashRMS = doAshScan_RMS(fIn.get(), N_E, E_edges, bScan, "rms", debugDir, skipLowest3);
  auto logRMS = doLogScan_RMS(fIn.get(), N_E, E_edges, w0Scan, "rms", debugDir, skipLowest3);

  drawAshLogSideBySide_RMS(ashRMS, logRMS, E_edges, N_E, outBaseDir, skipLowest3);
  makePhenixLikeSummary_RMS(ashRMS, E_edges, N_E, 5.55 /*cell cm*/, outBaseDir, skipLowest3);
  drawMinBOverlay_RMS_vs_OriginalEta(ashRMS, E_edges, N_E, outBaseDir, skipLowest3);

  // Also emit a tiny look‑up table: E_center → b_opt (RMS)
  std::ofstream lut( std::string(outBaseDir) + "/ashLogScans/bRMS_opt_vs_E.txt" );
  if (lut.is_open())
  {
    lut << "# E_center_GeV   b_opt_towerUnits   sigma_RMS_towerUnits\n";
    for (int i=0;i<N_E;++i)
    {
      const double eCtr = 0.5*(E_edges[i] + E_edges[i+1]);
      const bool ok = (i < (int)ashRMS.bestParam.size() && i < (int)ashRMS.bestSigma.size()
                       && std::isfinite(ashRMS.bestParam[i]) && std::isfinite(ashRMS.bestSigma[i]));
      lut << std::fixed << std::setprecision(3) << eCtr << "  "
          << std::setprecision(4) << (ok?ashRMS.bestParam[i]:NAN) << "  "
          << std::setprecision(6) << (ok?ashRMS.bestSigma[i]:NAN) << "\n";
    }
    lut.close();
    std::cout << "[ASH/LOG‑RMS] wrote LUT → " << (std::string(outBaseDir) + "/ashLogScans/bRMS_opt_vs_E.txt") << "\n";
  }
}

// ================================ MAIN ======================================

void PDCAnalysisPrime()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // --- paths
  const std::string inFilePath = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/PositionDep_sim_ALL.root";
  const std::string outBaseDir = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/SimOutputPrime";

  EnsureDir(outBaseDir);

  // Open file
  std::unique_ptr<TFile> fin(TFile::Open(inFilePath.c_str(),"READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open input file: " << inFilePath << "\n";
    return;
  }

  // Energy slices
  const auto eEdges = MakeEnergySlices();
  const auto eCent  = MakeEnergyCenters();

    RunPi0MassAnalysis(inFilePath.c_str(), outBaseDir.c_str());
    
    // Prepare global b-values text file
    const std::string bValFile = outBaseDir + "/bValues.txt";
    std::ofstream bOut(bValFile);
    if (!bOut.is_open())
    {
      std::cerr << "[WARN] Cannot open " << bValFile << " => won't save b-values.\n";
    }
    else
    {
      bOut << "# E range   best-b   (PHI or ETA)   etaVariant   zVariant\n";
    }
    
    MakeResidualsSuite(fin.get(), outBaseDir);

  // Variant definitions: list of possible names (range / disc) for each histogram
  struct VSpec {
    std::string key;
    std::vector<TString> uncNames;
    std::vector<TString> corNames;
  };

  const std::vector<VSpec> variants = {
    {"originalEta",
      {"h3_blockCoord_E_range", "h3_blockCoord_E_disc"},
      {"h3_blockCoord_Ecorr_range", "h3_blockCoord_Ecorr_disc"}
    },
    {"fullEta",
      {"h3_blockCoord_E_full_range", "h3_blockCoord_E_full_disc"},
      {"h3_blockCoord_Ecorr_full_range", "h3_blockCoord_Ecorr_full_disc"}
    },
    {"etaCore",
      {"h3_blockCoord_E_etaCore_le0p20_range", "h3_blockCoord_E_etaCore_le0p20_disc"},
      {"h3_blockCoord_Ecorr_etaCore_le0p20_range", "h3_blockCoord_Ecorr_etaCore_le0p20_disc"}
    },
    {"etaMid",
      {"h3_blockCoord_E_etaMid_0p20to0p70_range", "h3_blockCoord_E_etaMid_0p20to0p70_disc"},
      {"h3_blockCoord_Ecorr_etaMid_0p20to0p70_range", "h3_blockCoord_Ecorr_etaMid_0p20to0p70_disc"}
    },
    {"etaEdge",
      {"h3_blockCoord_E_etaEdge_0p70to1p10_range", "h3_blockCoord_E_etaEdge_0p70to1p10_disc"},
      {"h3_blockCoord_Ecorr_etaEdge_0p70to1p10_range", "h3_blockCoord_Ecorr_etaEdge_0p70to1p10_disc"}
    }
  };

    // To build the global overlays across variants:
    std::map<std::string, std::vector<double>> phiByVariant;     // key -> bφ(E)
    std::map<std::string, std::vector<double>> etaByVariant;     // key -> bη(E)
    std::map<std::string, std::vector<double>> phiErrByVariant;  // key -> σ_bφ(E)
    std::map<std::string, std::vector<double>> etaErrByVariant;  // key -> σ_bη(E)

    // NEW: store per-variant counts of clusters per energy bin (UNCORRECTED)
    std::map<std::string, std::vector<long long>> entriesByVariant;

    for (const auto& v : variants)
      {
        const std::string outDir = outBaseDir + "/" + v.key;
        EnsureDir(outDir);

        // --- run your normal per-variant products (lego, 2D, overlays, b vs E) ---
        TH3F* hUnc = GetTH3FByNames(fin.get(), v.uncNames);
        TH3F* hCor = nullptr;
        if (!isFirstPass) hCor = GetTH3FByNames(fin.get(), v.corNames);
        if (!hUnc) { std::cerr << "[WARN] missing " << v.key << " - skipped\n"; continue; }
        hUnc->SetDirectory(nullptr); if (hCor) hCor->SetDirectory(nullptr);

        SaveTH3Lego(hUnc, (outDir + "/lego_unc.png").c_str(), kEtaPretty.at(v.key).c_str());
        Plot2DBlockEtaPhi(hUnc, hCor, isFirstPass, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
        OverlayUncorrPhiEta (hUnc, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
          // 1) Make the standard (filled) b(E) plot for this view and cache values
          BVecs bv = MakeBvaluesVsEnergyPlot(hUnc, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
          phiByVariant[v.key]    = bv.bphi;  phiErrByVariant[v.key] = bv.bphiErr;
          etaByVariant[v.key]    = bv.beta;  etaErrByVariant[v.key] = bv.betaErr;
          if (bOut.is_open()) WriteBValuesTxt(hUnc, eEdges, bOut, v.key.c_str(), "originalZRange");

//          // 2) For originalEta ONLY: read sim file, compute b(E), and save an overlay PNG
//          if (v.key == "originalEta")
//          {
//            const char* simPath = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/PositionDep_sim_ALL.root";
//            std::unique_ptr<TFile> fSim(TFile::Open(simPath, "READ"));
//            if (fSim && !fSim->IsZombie())
//            {
//              // try both booking name styles for the uncorrected 3D (same as your production)
//              TH3F* hSimUnc = GetTH3FByNames(fSim.get(), {
//                "h3_blockCoord_E_range", "h3_blockCoord_E_disc"
//              });
//              if (hSimUnc)
//              {
//                hSimUnc->SetDirectory(nullptr);
//
//                // --- compute sim b(E) without making its own plot
//                auto computeBvecsNoPlot = [&](TH3F* h3, const std::vector<std::pair<double,double>>& edges)->BVecs
//                {
//                  BVecs R;
//                  if (!h3) return R;
//                  const int Nbin = static_cast<int>(edges.size());
//                  R.ecenters.reserve(Nbin); R.bphi.reserve(Nbin); R.beta.reserve(Nbin);
//                  R.bphiErr.reserve(Nbin); R.betaErr.reserve(Nbin);
//
//                  for (int i=0; i<Nbin; ++i)
//                  {
//                    const double eLo = edges[i].first, eHi = edges[i].second;
//                    const int zLo = std::max(1, h3->GetZaxis()->FindBin(eLo + 1e-9));
//                    const int zHi = std::min(h3->GetNbinsZ(), h3->GetZaxis()->FindBin(eHi - 1e-9));
//
//                    // φ (Y projection)
//                    h3->GetZaxis()->SetRange(zLo, zHi);
//                    h3->GetXaxis()->SetRange(1, h3->GetNbinsX());
//                    TH1D* hPhi = static_cast<TH1D*>(h3->Project3D("y"));
//                    if (!hPhi) continue;
//                    hPhi->SetDirectory(nullptr);
//
//                    // η (X projection)
//                    h3->GetYaxis()->SetRange(1, h3->GetNbinsY());
//                    TH1D* hEta = static_cast<TH1D*>(h3->Project3D("x"));
//                    if (!hEta) { delete hPhi; continue; }
//                    hEta->SetDirectory(nullptr);
//
//                    const BRes bPhi = FitAsinh1D(hPhi);
//                    const BRes bEta = FitAsinh1D(hEta);
//
//                    R.ecenters.push_back(0.5*(eLo+eHi));
//                    R.bphi.push_back(bPhi.val);
//                    R.beta.push_back(bEta.val);
//                    R.bphiErr.push_back(bPhi.err);
//                    R.betaErr.push_back(bEta.err);
//
//                    delete hPhi;
//                    delete hEta;
//                  }
//                  return R;
//                };
//
//                BVecs bSim = computeBvecsNoPlot(hSimUnc, eEdges);
//
//                // --- build an overlay canvas with BOTH: data (filled) & sim (open circle)
//                TCanvas cOvB("bValuesOverlay","b values overlay (data vs sim)", 900, 700);
//
//                // Find sensible Y range across both sets (go lower than before)
//                double ymin=+1e9, ymax=-1e9;
//                auto upd = [&](const std::vector<double>& v){
//                  for (double q : v) if (std::isfinite(q)) { ymin = std::min(ymin,q); ymax = std::max(ymax,q); }
//                };
//                upd(bv.bphi); upd(bv.beta); upd(bSim.bphi); upd(bSim.beta);
//                if (!(std::isfinite(ymin)&&std::isfinite(ymax))) { ymin=0.0; ymax=1.0; }
//                const double span = std::max(1e-6, ymax - ymin);
//                const double yLo  = ymin - 0.50*span;             // ↓ more headroom below than before
//                const double yHi  = ymax + 0.25*span;
//
//                const double xLo = E_edges[0]-0.5, xHi = E_edges[eEdges.size()]-0.5;
//
//                TH2F* fr = new TH2F("bFrameOverlay",
//                                    Form("best-fit  b  vs  E  (#scale[0.8]{%s});E  [GeV];b", kEtaPretty.at(v.key).c_str()),
//                                    100, xLo, xHi, 100, yLo, yHi);
//                fr->SetStats(0);
//                fr->Draw();
//
//                // Data (filled) — use CLOSED CIRCLES for b_phi (red) and keep b_eta (blue) as closed circles
//                TGraphErrors gPhi_data(bv.ecenters.size(), &bv.ecenters[0], &bv.bphi[0],  nullptr, &bv.bphiErr[0]);
//                TGraphErrors gEta_data(bv.ecenters.size(), &bv.ecenters[0], &bv.beta[0],  nullptr, &bv.betaErr[0]);
//                gPhi_data.SetMarkerStyle(20); gPhi_data.SetMarkerColor(kRed+1);   gPhi_data.SetLineColor(kRed+1);   // closed red circle
//                gEta_data.SetMarkerStyle(20); gEta_data.SetMarkerColor(kBlue+1);  gEta_data.SetLineColor(kBlue+1);  // closed blue circle
//
//                // Simulation (open circle markers), same colors
//                TGraphErrors gPhi_sim(bSim.ecenters.size(), &bSim.ecenters[0], &bSim.bphi[0],  nullptr, &bSim.bphiErr[0]);
//                TGraphErrors gEta_sim(bSim.ecenters.size(), &bSim.ecenters[0], &bSim.beta[0],  nullptr, &bSim.betaErr[0]);
//                gPhi_sim.SetMarkerStyle(24); gPhi_sim.SetMarkerColor(kRed+1);   gPhi_sim.SetLineColor(kRed+1);     // open red circle
//                gEta_sim.SetMarkerStyle(24); gEta_sim.SetMarkerColor(kBlue+1);  gEta_sim.SetLineColor(kBlue+1);    // open blue circle
//
//                // Draw order: data first, then sim on top
//                gPhi_data.Draw("P SAME");
//                gEta_data.Draw("P SAME");
//                gPhi_sim.Draw("P SAME");
//                gEta_sim.Draw("P SAME");
//
//                // Legend (top-right) — TWO COLUMNS: left=DATA, right=SIM; top row=b_phi, bottom row=b_eta
//                TLegend legB(0.58,0.78,0.9,0.90);
//                legB.SetBorderSize(0); legB.SetFillStyle(0); legB.SetTextSize(0.040);
//                legB.SetNColumns(2);
//                legB.SetColumnSeparation(0.10);
//                // row 1: b_phi (data, sim)
//                legB.AddEntry(&gPhi_data, "b_{#varphi} data", "lp");
//                legB.AddEntry(&gPhi_sim , "b_{#varphi}  sim", "lp");
//                // row 2: b_eta (data, sim)
//                legB.AddEntry(&gEta_data, "b_{#eta} data", "lp");
//                legB.AddEntry(&gEta_sim , "b_{#eta}  sim", "lp");
//                legB.Draw();
//
//                // Save alongside the standard product in originalEta/
//                TString outOverlay = Form("%s/bValuesOverlay.png", outDir.c_str());
//                cOvB.SaveAs(outOverlay);
//
//                delete fr;
//                delete hSimUnc;
//              } // if hSimUnc
//            }   // if fSim ok
//          }     // if originalEta


        if (!isFirstPass && hCor) { SaveTH3Lego(hCor, (outDir + "/lego_cor.png").c_str(), kEtaPretty.at(v.key).c_str());
                                    FitLocalPhiEta(hUnc, hCor, false, eEdges, outDir.c_str()); }

        // --- NEW: run the Playground only for originalEta into .../originalEta/scratchPDC ---
        if (v.key == "originalEta")
        {
          const std::string scratch = outDir + "/scratchPDC";
          EnsureDir(scratch);
          MakeDeltaPhiEtaPlayground(inFilePath.c_str(), scratch.c_str(),
                                    /*xMin*/ -0.035, /*xMax*/ 0.035);
        }

        // --- NEW: compute entries per energy bin (UNCORRECTED 3D → sum over X,Y in that Z slice) ---
        {
          const int nx = hUnc->GetNbinsX();
          const int ny = hUnc->GetNbinsY();
          std::vector<long long> cntVec;
          cntVec.reserve(eEdges.size());
          for (size_t i = 0; i < eEdges.size(); ++i) {
            const double eLo = eEdges[i].first;
            const double eHi = eEdges[i].second;
            const int zLo = std::max(1, hUnc->GetZaxis()->FindBin(eLo + 1e-9));
            const int zHi = std::min(hUnc->GetNbinsZ(), hUnc->GetZaxis()->FindBin(eHi - 1e-9));
            // Include full XY range, and only the Z slice corresponding to the energy bin:
            const double sum = hUnc->Integral(1, nx, 1, ny, zLo, zHi);
            const long long cnt = static_cast<long long>(std::llround(sum));
            cntVec.push_back(cnt);
          }
          entriesByVariant[v.key] = std::move(cntVec);
        }

        delete hUnc; if (hCor) delete hCor;
    }

    // NEW: Print entries summary tables (UNCORRECTED)
    //   1) Across η-views (you already computed entriesByVariant above)
    //   2) Across |z_vtx| slices by reading the h3_blockCoord_E_zXXtoYY_* histos from file
    {
      // ------------------ 1) η-views table ------------------
      {
        // Preferred column order; print only those we actually filled:
        std::vector<std::string> colOrder = {"etaCore","etaMid","etaEdge","fullEta","originalEta"};
        std::vector<std::string> cols; cols.reserve(colOrder.size());
        for (const auto& k : colOrder) {
          if (entriesByVariant.count(k)) cols.push_back(k);
        }

        if (!cols.empty()) {
          std::cout << "\n";
          std::cout << "╔═════════════════════════════════════════════════════════════════════════════════════╗\n";
          std::cout << "║  Entries summary: number of clusters per E-bin and |η| bin (uncorrected Z-slices)  ║\n";
          std::cout << "╚═════════════════════════════════════════════════════════════════════════════════════╝\n";
          // Header
          std::cout << std::left << std::setw(20) << "E-bin";
          for (const auto& key : cols) {
            std::cout << " | " << std::setw(14) << kEtaPretty.at(key);
          }
          std::cout << "\n";
          std::cout << std::string(20 + (int)cols.size() * (3 + 14), '-') << "\n";

          // Rows
          for (size_t i = 0; i < eEdges.size(); ++i) {
            std::cout << std::left << std::setw(20)
                      << Form("[%.1f,%.1f) GeV", eEdges[i].first, eEdges[i].second)
                      << std::right;
            for (const auto& key : cols) {
              const auto& vec = entriesByVariant.at(key);
              const long long val = (i < vec.size()) ? vec[i] : 0LL;
              std::cout << " | " << std::setw(14) << val;
            }
            std::cout << "\n";
          }

          // Totals row
          std::cout << std::string(20 + (int)cols.size() * (3 + 14), '-') << "\n";
          std::cout << std::left << std::setw(20) << "Total" << std::right;
          for (const auto& key : cols) {
            long long tot = 0;
            for (const auto& v : entriesByVariant.at(key)) tot += v;
            std::cout << " | " << std::setw(14) << tot;
          }
          std::cout << "\n\n";
        }
      }

        // ------------------ 2) |z_vtx|-slice table (COARSE) ------------------
        {
          // Try to read each z-slice TH3 (either *_range or *_disc) and compute per-E-bin counts
          // The order/labels come from your kZOrder / kZPretty.
          std::map<std::string, std::vector<long long>> zEntries;   // zTag -> counts per E-bin
          std::vector<std::string> zCols; zCols.reserve(kZOrder.size());

          for (const auto& ztag : kZOrder)
          {
            std::vector<TString> zNames = {
              Form("h3_blockCoord_E_%s_range", ztag.c_str()),
              Form("h3_blockCoord_E_%s_disc" , ztag.c_str())
            };
            TH3F* hZ = GetTH3FByNames(fin.get(), zNames);
            if (!hZ) continue;

            hZ->SetDirectory(nullptr);
            const int nx = hZ->GetNbinsX();
            const int ny = hZ->GetNbinsY();

            std::vector<long long> cntVec; cntVec.reserve(eEdges.size());
            for (size_t i = 0; i < eEdges.size(); ++i) {
              const double eLo = eEdges[i].first;
              const double eHi = eEdges[i].second;
              const int zLo = std::max(1, hZ->GetZaxis()->FindBin(eLo + 1e-9));
              const int zHi = std::min(hZ->GetNbinsZ(), hZ->GetZaxis()->FindBin(eHi - 1e-9));
              const double sum = hZ->Integral(1, nx, 1, ny, zLo, zHi);
              const long long cnt = static_cast<long long>(std::llround(sum));
              cntVec.push_back(cnt);
            }

            zEntries[ztag] = std::move(cntVec);
            zCols.push_back(ztag);
            delete hZ;
          }

          if (!zCols.empty()) {
            std::cout << "\n";
            std::cout << "╔══════════════════════════════════════════════════════════════════════════════════════╗\n";
            std::cout << "║  Entries summary: number of clusters per E-bin and |z_vtx| bin (uncorrected z-slices)║\n";
            std::cout << "╚══════════════════════════════════════════════════════════════════════════════════════╝\n";
            // Header
            std::cout << std::left << std::setw(20) << "E-bin";
            for (const auto& ztag : zCols) {
              std::cout << " | " << std::setw(14) << kZPretty.at(ztag);
            }
            std::cout << "\n";
            std::cout << std::string(20 + (int)zCols.size() * (3 + 14), '-') << "\n";

            // Rows
            for (size_t i = 0; i < eEdges.size(); ++i) {
              std::cout << std::left << std::setw(20)
                        << Form("[%.1f,%.1f) GeV", eEdges[i].first, eEdges[i].second)
                        << std::right;
              for (const auto& ztag : zCols) {
                const auto& vec = zEntries.at(ztag);
                const long long val = (i < vec.size()) ? vec[i] : 0LL;
                std::cout << " | " << std::setw(14) << val;
              }
              std::cout << "\n";
            }

            // Totals row
            std::cout << std::string(20 + (int)zCols.size() * (3 + 14), '-') << "\n";
            std::cout << std::left << std::setw(20) << "Total" << std::right;
            for (const auto& ztag : zCols) {
              long long tot = 0;
              for (const auto& v : zEntries.at(ztag)) tot += v;
              std::cout << " | " << std::setw(14) << tot;
            }
            std::cout << "\n\n";
          }
        }

        // ------------------ 2b) |z_vtx|-slice table (FINE 0–10 cm) ------------------
        {
          // Read each fine z-slice TH3 (z00to02 ... z08to10) and compute per-E-bin counts
          // The order/labels come from kZFineOrder / kZFinePretty.
          std::map<std::string, std::vector<long long>> zEntriesFine;   // fine zTag -> counts
          std::vector<std::string> zColsFine; zColsFine.reserve(kZFineOrder.size());

          for (const auto& ztag : kZFineOrder)
          {
            std::vector<TString> zNamesFine = {
              Form("h3_blockCoord_E_%s_range", ztag.c_str()),
              Form("h3_blockCoord_E_%s_disc" , ztag.c_str())
            };
            TH3F* hZf = GetTH3FByNames(fin.get(), zNamesFine);
            if (!hZf) continue;

            hZf->SetDirectory(nullptr);
            const int nx = hZf->GetNbinsX();
            const int ny = hZf->GetNbinsY();

            std::vector<long long> cntVec; cntVec.reserve(eEdges.size());
            for (size_t i = 0; i < eEdges.size(); ++i) {
              const double eLo = eEdges[i].first;
              const double eHi = eEdges[i].second;
              const int zLo = std::max(1, hZf->GetZaxis()->FindBin(eLo + 1e-9));
              const int zHi = std::min(hZf->GetNbinsZ(), hZf->GetZaxis()->FindBin(eHi - 1e-9));
              const double sum = hZf->Integral(1, nx, 1, ny, zLo, zHi);
              const long long cnt = static_cast<long long>(std::llround(sum));
              cntVec.push_back(cnt);
            }

            zEntriesFine[ztag] = std::move(cntVec);
            zColsFine.push_back(ztag);
            delete hZf;
          }

          if (!zColsFine.empty()) {
            std::cout << "\n";
            std::cout << "╔════════════════════════════════════════════════════════════════════════════════════════════╗\n";
            std::cout << "║  Entries summary: number of clusters per E-bin and |z_vtx| bin (uncorrected fine z-slices 0–10 cm) ║\n";
            std::cout << "╚════════════════════════════════════════════════════════════════════════════════════════════╝\n";
            // Header
            std::cout << std::left << std::setw(20) << "E-bin";
            for (const auto& ztag : zColsFine) {
              std::cout << " | " << std::setw(14) << kZFinePretty.at(ztag);
            }
            std::cout << "\n";
            std::cout << std::string(20 + (int)zColsFine.size() * (3 + 14), '-') << "\n";

            // Rows
            for (size_t i = 0; i < eEdges.size(); ++i) {
              std::cout << std::left << std::setw(20)
                        << Form("[%.1f,%.1f) GeV", eEdges[i].first, eEdges[i].second)
                        << std::right;
              for (const auto& ztag : zColsFine) {
                const auto& vec = zEntriesFine.at(ztag);
                const long long val = (i < vec.size()) ? vec[i] : 0LL;
                std::cout << " | " << std::setw(14) << val;
              }
              std::cout << "\n";
            }

            // Totals row
            std::cout << std::string(20 + (int)zColsFine.size() * (3 + 14), '-') << "\n";
            std::cout << std::left << std::setw(20) << "Total" << std::right;
            for (const auto& ztag : zColsFine) {
              long long tot = 0;
              for (const auto& v : zEntriesFine.at(ztag)) tot += v;
              std::cout << " | " << std::setw(14) << tot;
            }
            std::cout << "\n\n";
          }
        }
      }
    // ----------------------- Per-|z_vtx| slices (uncorr always; corrected if available) -----------------------
    {
      // ========== COARSE |z| bins: z00to10, z10to20, z20to30, z30to45, z45to60 ==========
      std::map<std::string, std::vector<double>> zPhiByVar, zEtaByVar;
      std::map<std::string, std::vector<double>> zPhiErrByVar, zEtaErrByVar;

      for (const auto& ztag : kZOrder)
      {
        // Uncorrected and corrected names
        std::vector<TString> zNamesUnc = {
          Form("h3_blockCoord_E_%s_range", ztag.c_str()),
          Form("h3_blockCoord_E_%s_disc" , ztag.c_str())
        };
        std::vector<TString> zNamesCor = {
          Form("h3_blockCoord_Ecorr_%s_range", ztag.c_str()),
          Form("h3_blockCoord_Ecorr_%s_disc" , ztag.c_str())
        };

        TH3F* hZ_unc = GetTH3FByNames(fin.get(), zNamesUnc);
        if (!hZ_unc) {
          std::cerr << "[WARN] missing z-slice TH3F (uncorr): " << ztag << " - skipped\n";
          continue;
        }
        hZ_unc->SetDirectory(nullptr);

        TH3F* hZ_cor = (!isFirstPass ? GetTH3FByNames(fin.get(), zNamesCor) : nullptr);
        if (hZ_cor) hZ_cor->SetDirectory(nullptr);

        const std::string outDirZ = outBaseDir + "/" + ztag;
        EnsureDir(outDirZ);

        // Uncorrected products (always)
        SaveTH3Lego(hZ_unc, (outDirZ + "/lego_unc.png").c_str(), kZPretty.at(ztag).c_str());
        Plot2DBlockEtaPhi(hZ_unc, hZ_cor, /*firstPass*/isFirstPass, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());
        OverlayUncorrPhiEta(hZ_unc, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());
        BVecs zbv = MakeBvaluesVsEnergyPlot(hZ_unc, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());

        // Persist per-z b(E) for overlays
        zPhiByVar[ztag]     = zbv.bphi;  zPhiErrByVar[ztag] = zbv.bphiErr;
        zEtaByVar[ztag]     = zbv.beta;  zEtaErrByVar[ztag] = zbv.betaErr;

        // Write b-values tagged by z-slice
        if (bOut.is_open()) WriteBValuesTxt(hZ_unc, eEdges, bOut, "originalEta", ztag.c_str());

        // Corrected products (only if second pass and histogram exists)
        if (hZ_cor && !isFirstPass)
        {
          SaveTH3Lego(hZ_cor, (outDirZ + "/lego_cor.png").c_str(), kZPretty.at(ztag).c_str());
          FitLocalPhiEta(hZ_unc, hZ_cor, /*isFirstPassFlag*/false, eEdges, outDirZ.c_str());
        }

        delete hZ_unc;
        if (hZ_cor) delete hZ_cor;
      }

      if (!zPhiByVar.empty())
      {
        SaveOverlayAcrossVariants(outBaseDir, eCent, zPhiByVar, zPhiErrByVar, /*isPhi*/true , "_zOnly");
        SaveOverlayAcrossVariants(outBaseDir, eCent, zEtaByVar, zEtaErrByVar, /*isPhi*/false, "_zOnly");
      }

      // ========== FINE |z| bins: 0–2, 2–4, 4–6, 6–8, 8–10 cm ==========
      std::map<std::string, std::vector<double>> zPhiFineByVar, zEtaFineByVar;
      std::map<std::string, std::vector<double>> zPhiFineErrByVar, zEtaFineErrByVar;

      for (const auto& ztag : kZFineOrder)
      {
        std::vector<TString> zNamesFineUnc = {
          Form("h3_blockCoord_E_%s_range", ztag.c_str()),
          Form("h3_blockCoord_E_%s_disc" , ztag.c_str())
        };
        std::vector<TString> zNamesFineCor = {
          Form("h3_blockCoord_Ecorr_%s_range", ztag.c_str()),
          Form("h3_blockCoord_Ecorr_%s_disc" , ztag.c_str())
        };

        TH3F* hZf_unc = GetTH3FByNames(fin.get(), zNamesFineUnc);
        if (!hZf_unc) continue; // fine z is optional
        hZf_unc->SetDirectory(nullptr);

        TH3F* hZf_cor = (!isFirstPass ? GetTH3FByNames(fin.get(), zNamesFineCor) : nullptr);
        if (hZf_cor) hZf_cor->SetDirectory(nullptr);

        const std::string outDirZf = outBaseDir + "/" + ztag;
        EnsureDir(outDirZf);

        // Uncorrected fine products
        SaveTH3Lego(hZf_unc, (outDirZf + "/lego_unc.png").c_str(), kZFinePretty.at(ztag).c_str());
        Plot2DBlockEtaPhi(hZf_unc, hZf_cor, /*firstPass*/isFirstPass, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());
        OverlayUncorrPhiEta(hZf_unc, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());
        BVecs zbvFine = MakeBvaluesVsEnergyPlot(hZf_unc, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());

        zPhiFineByVar[ztag]     = zbvFine.bphi;  zPhiFineErrByVar[ztag] = zbvFine.bphiErr;
        zEtaFineByVar[ztag]     = zbvFine.beta;  zEtaFineErrByVar[ztag] = zbvFine.betaErr;

        // Log fine-z b-values
        if (bOut.is_open()) WriteBValuesTxt(hZf_unc, eEdges, bOut, "originalEta", ztag.c_str());

        // Corrected fine products (second pass only)
        if (hZf_cor && !isFirstPass)
        {
          SaveTH3Lego(hZf_cor, (outDirZf + "/lego_cor.png").c_str(), kZFinePretty.at(ztag).c_str());
          FitLocalPhiEta(hZf_unc, hZf_cor, /*isFirstPassFlag*/false, eEdges, outDirZf.c_str());
        }

        delete hZf_unc;
        if (hZf_cor) delete hZf_cor;
      }

      if (!zPhiFineByVar.empty())
      {
        SaveOverlayAcrossVariants(outBaseDir, eCent, zPhiFineByVar, zPhiFineErrByVar, /*isPhi*/true ,  "_zFine");
        SaveOverlayAcrossVariants(outBaseDir, eCent, zEtaFineByVar, zEtaFineErrByVar, /*isPhi*/false, "_zFine");
      }
    }

    // Close b-values after all z-slice writes (coarse + fine)
    if (bOut.is_open()) bOut.close();

    // ----------------------- Global overlays across η-variants (base path) -----------------------
    SaveOverlayAcrossVariants(outBaseDir, eCent, phiByVariant, phiErrByVariant, /*isPhi*/true , ""  );
    SaveOverlayAcrossVariants(outBaseDir, eCent, etaByVariant, etaErrByVariant, /*isPhi*/false, ""  );

    AnalyzeBvsResCLUSCPEA(inFilePath.c_str(), outBaseDir.c_str());


    
  std::cout << "[DONE] Outputs written under: " << outBaseDir << "\n";
}

// Auto-run when invoked as a ROOT script:  root -b -q -l PDCAnalysisPrime.cpp
PDCAnalysisPrime();
