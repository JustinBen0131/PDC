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
  {"z00to10",   "|z_{vtx}| 0-10 cm"},
  {"z10to20",   "|z_{vtx}| 10-20 cm"},
  {"z20to30",   "|z_{vtx}| 20-30 cm"},
  {"z30to45",   "|z_{vtx}| 30-45 cm"},
  {"z45to60",   "|z_{vtx}| 45-60 cm"}
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
  trial.SetParLimits(1, 1e-5, 5.0);

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

    // Build scatter plot + draw WLS fits and annotate m, b0 on the canvas
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

    // --- Weighted least squares on x = ln(E/E0) to match SaveOverlayAcrossVariants ---
    const double E0 = 3.0;

    auto wlsFit = [&](const std::vector<double>& E,
                      const std::vector<double>& y,
                      const std::vector<double>& sy,
                      double& m, double& b0, double& dm, double& db0) -> bool
    {
      double S=0, Sx=0, Sy_=0, Sxx=0, Sxy=0;
      int Nuse=0; bool haveYerrsAll=true;
      std::vector<double> xi, yi, wi; xi.reserve(y.size()); yi.reserve(y.size()); wi.reserve(y.size());
      for (size_t i=0;i<y.size() && i<E.size(); ++i)
      {
        const double Ei = E[i], yi_ = y[i];
        if (!(std::isfinite(Ei) && Ei>0.0 && std::isfinite(yi_))) continue;
        const double xi_ = std::log(Ei/E0);
        double w = 1.0;
        if (i<sy.size() && std::isfinite(sy[i]) && sy[i]>0.0) w = 1.0/(sy[i]*sy[i]); else haveYerrsAll=false;
        S+=w; Sx+=w*xi_; Sy_+=w*yi_; Sxx+=w*xi_*xi_; Sxy+=w*xi_*yi_;
        xi.push_back(xi_); yi.push_back(yi_); wi.push_back(w); ++Nuse;
      }
      const double Delta = S*Sxx - Sx*Sx;
      if (Nuse<2 || !(Delta>0.0)) return false;

      m  = (S*Sxy - Sx*Sy_) / Delta;
      b0 = (Sxx*Sy_ - Sx*Sxy) / Delta;

      const double var_m  =  S / Delta;
      const double var_b0 = Sxx / Delta;

      // Goodness-of-fit and error scaling if we had to use unit weights
      double chi2=0.0, ybar_w = Sy_/S; const int ndf = std::max(0, Nuse-2);
      for (int i=0;i<Nuse; ++i) {
        const double yhat = b0 + m*xi[i];
        const double r    = yi[i] - yhat;
        chi2 += wi[i]*r*r;
      }
      const double scale = (!haveYerrsAll && ndf>0) ? std::sqrt(chi2/std::max(1,ndf)) : 1.0;

      dm  = std::sqrt(std::max(0.0, var_m )) * scale;
      db0 = std::sqrt(std::max(0.0, var_b0)) * scale;
      return true;
    };

    double m_phi=std::numeric_limits<double>::quiet_NaN(), b0_phi=std::numeric_limits<double>::quiet_NaN();
    double dm_phi=std::numeric_limits<double>::quiet_NaN(), db0_phi=std::numeric_limits<double>::quiet_NaN();
    double m_eta=std::numeric_limits<double>::quiet_NaN(), b0_eta=std::numeric_limits<double>::quiet_NaN();
    double dm_eta=std::numeric_limits<double>::quiet_NaN(), db0_eta=std::numeric_limits<double>::quiet_NaN();

    const bool okPhi = wlsFit(R.ecenters, R.bphi, R.bphiErr, m_phi, b0_phi, dm_phi, db0_phi);
    const bool okEta = wlsFit(R.ecenters, R.beta, R.betaErr, m_eta, b0_eta, dm_eta, db0_eta);

    // Draw the fit curves and keep them alive until after SaveAs
    std::vector<std::unique_ptr<TF1>> keepFitLines;
    keepFitLines.reserve(2);

    if (okPhi) {
      const auto fPhiName = Form("fPhi_line_%p", static_cast<void*>(gPhi));
      auto fPhi = std::make_unique<TF1>(fPhiName, "[0] + [1]*log(x/[2])", xLo, xHi);
      fPhi->SetParameters(b0_phi, m_phi, E0);
      fPhi->FixParameter(2, E0);
      fPhi->SetLineColor(kRed+1);
      fPhi->SetLineWidth(2);
      fPhi->SetNpx(800);
      fPhi->Draw("SAME");
      keepFitLines.emplace_back(std::move(fPhi));
    }
    if (okEta) {
      const auto fEtaName = Form("fEta_line_%p", static_cast<void*>(gEta));
      auto fEta = std::make_unique<TF1>(fEtaName, "[0] + [1]*log(x/[2])", xLo, xHi);
      fEta->SetParameters(b0_eta, m_eta, E0);
      fEta->FixParameter(2, E0);
      fEta->SetLineColor(kBlue+1);
      fEta->SetLineWidth(2);
      fEta->SetNpx(800);
      fEta->Draw("SAME");
      keepFitLines.emplace_back(std::move(fEta));
    }

    // Legend for the points (as before)
    TLegend leg(0.20,0.78,0.52,0.90);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
    leg.AddEntry(gPhi,"b_{#varphi}","lp");
    leg.AddEntry(gEta,"b_{#eta}","lp");
    leg.Draw();

    // Stamp fit parameters in the top-right (matching colors)
    {
      TLatex tl; tl.SetNDC(); tl.SetTextAlign(33); tl.SetTextSize(0.040);
      if (okPhi) { tl.SetTextColor(kRed+1);
        tl.DrawLatex(0.88, 0.88, Form("#phi: b_{0}=%.3f,  m=%.3f",
                                      b0_phi, m_phi));
      }
      if (okEta) { tl.SetTextColor(kBlue+1);
        tl.DrawLatex(0.88, 0.82, Form("#eta: b_{0}=%.3f,  m=%.3f",
                                      b0_eta, m_eta));
      }
      tl.SetTextColor(kBlack); // restore
    }

    cB.Modified();
    cB.Update();
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

  const double E0   = 3.0;
  const double xLo  = E_edges[0] - 0.5;
  const double xHi  = E_edges[N_E] - 0.5;
  const std::string baseSfx = outSuffix ? std::string(outSuffix) : std::string();

    // ---------- master TXT (6-column) path & header (always written at top-level) ----------
    auto topBaseDir = outBaseDir;
    const std::string token = "/etaAndzDepBlocks/";
    const auto p = topBaseDir.find(token);
    if (p != std::string::npos) topBaseDir = topBaseDir.substr(0, p); // strip /etaAndzDepBlocks/<etaKey>
    const std::string masterTXT = topBaseDir + "/bFit_master.txt";

    auto ensureMasterHeader = [&](const std::string& path)
    {
      bool need = true;
      {
        std::ifstream fin(path);
        if (fin.good()) { fin.seekg(0, std::ios::end); need = (fin.tellg() == 0); }
      }
      if (need) {
        std::ofstream out(path, std::ios::app);
        // space-delimited header
        out << "variantName etaOrPhi etaRange zRange m b0\n";
      }
    };
    ensureMasterHeader(masterTXT);

  // Detect if this call is under .../etaAndzDepBlocks/<etaKey> and fetch <etaKey>
  std::string etaKeyFromPath;
  if (p != std::string::npos) {
    size_t start = p + token.size();
    size_t slash = outBaseDir.find('/', start);
    etaKeyFromPath = (slash == std::string::npos)
                     ? outBaseDir.substr(start)
                     : outBaseDir.substr(start, slash - start);
  }

  // ---------- helpers ----------
  auto prettyName = [&](const std::string& key)->const char*
  {
    auto itEta = kEtaPretty.find(key);
    if (itEta != kEtaPretty.end()) return itEta->second.c_str();
    auto itZf  = kZFinePretty.find(key);
    if (itZf  != kZFinePretty.end()) return itZf->second.c_str();
    auto itZc  = kZPretty.find(key);
    return (itZc == kZPretty.end()) ? key.c_str() : itZc->second.c_str();
  };

  auto containsKey = [](const std::vector<std::string>& v, const std::string& s)->bool {
    for (const auto& x : v) if (x == s) return true; return false;
  };

  // ---------- classification of the provided map keys ----------
  bool hasOnlyCoarse = true, hasOnlyFine = true, hasAnyZ = false, hasAnyEta = false;
  for (const auto& kv : byVar) {
    if (kZPretty.find(kv.first)     == kZPretty.end()) hasOnlyCoarse = false; else hasAnyZ = true;
    if (kZFinePretty.find(kv.first) == kZFinePretty.end()) hasOnlyFine = false; else hasAnyZ = true;
    if (kEtaPretty.find(kv.first)   != kEtaPretty.end()) hasAnyEta = true;
  }

  // Decide where overlay images/CSVs should go for THIS call
  //  - z-only (coarse or fine) at top level --> outBaseDir + "/zDepOnlybvalues"
  //  - eta-only (across eta variants)       --> outBaseDir + "/etaDepOnlyBlockAna"
  //  - inside etaAndzDepBlocks/<etaKey>     --> keep in outBaseDir (and label filenames with _<etaKey>)
  std::string overlayOutDir;
  std::string suffixLabelTag;  // e.g. "_etaCore" when called under etaAndzDepBlocks/etaCore
  const bool inEtaAndZContext = (!etaKeyFromPath.empty() && hasAnyZ);

  if (inEtaAndZContext) {
    overlayOutDir   = outBaseDir;                // already .../etaAndzDepBlocks/<etaKey>
    suffixLabelTag  = "_" + etaKeyFromPath;      // tag filenames with the eta-range
  } else if (hasAnyEta && !hasAnyZ) {
    overlayOutDir   = outBaseDir + "/etaDepOnlyBlockAna";
    suffixLabelTag  = "";                        // no per-eta tag for the global eta-only overlay
  } else {
    overlayOutDir   = outBaseDir + "/zDepOnlybvalues";
    suffixLabelTag  = "";                        // z-only overlays at top level
  }
  EnsureDir(overlayOutDir);

  // ---------- local WLS helper: slope m, intercept b0 at E0 ----------
  auto wlsFit = [&](const std::vector<double>& E,
                    const std::vector<double>& y,
                    const std::vector<double>* sy) -> std::pair<double,double>
  {
    double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
    int Nuse=0;
    for (size_t i=0;i<y.size() && i<E.size();++i)
    {
      const double Ei = E[i];
      const double yi = y[i];
      if (!(std::isfinite(Ei) && Ei>0.0 && std::isfinite(yi))) continue;
      const double xi = std::log(Ei/E0);
      double wi = 1.0;
      if (sy && i<sy->size() && std::isfinite((*sy)[i]) && (*sy)[i] > 0.0) wi = 1.0/((*sy)[i]*(*sy)[i]);
      S   += wi; Sx += wi*xi; Sy += wi*yi; Sxx += wi*xi*xi; Sxy += wi*xi*yi;
      ++Nuse;
    }
    const double Delta = S*Sxx - Sx*Sx;
    if (Nuse < 2 || !(Delta > 0.0)) return { std::numeric_limits<double>::quiet_NaN(),
                                             std::numeric_limits<double>::quiet_NaN() };
    const double m  = (S*Sxy - Sx*Sy) / Delta;
    const double b0 = (Sxx*Sy - Sx*Sxy) / Delta;
    return {m, b0};
  };

    // ---------- write master 6-column TXT rows for ALL keys in this call (DE-DUPED) ----------
    {
      // Use a run-local set keyed by the tuple (variantName, axis, etaRange, zRange).
      // std::map is used here (already included via <map>) to avoid needing <set>/<unordered_set>.
      static std::map<std::string, char> g_master_seen;  // value unused; presence is enough

      std::ofstream mout(masterTXT, std::ios::app);
      if (mout)
      {
        for (const auto& kv : byVar)
        {
          const std::string& key  = kv.first;
          const auto&         vec = kv.second;
          auto itE = byVarErr.find(key);
          const std::vector<double>* err = (itE == byVarErr.end()) ? nullptr : &itE->second;

          const auto fit = wlsFit(ecenters, vec, err);
          if (!std::isfinite(fit.first) || !std::isfinite(fit.second)) continue;

          std::string variantName, etaRange, zRange;

          if (!etaKeyFromPath.empty() && hasAnyZ) {
            // called from .../etaAndzDepBlocks/<etaKey> with z keys
            variantName = "zAndEtaAndEnergyDep";
            etaRange    = etaKeyFromPath;
            zRange      = key;
          } else if (hasOnlyCoarse || hasOnlyFine || hasAnyZ) {
            // z-only overlays at top level (coarse or fine)
            variantName = "zAndEnergyDep";
            etaRange    = "originalEta";
            zRange      = key;
          } else {
            // η-only overlays
            if (key == "originalEta") variantName = "energyDepOnly";
            else                      variantName = "etaAndEnergyDep";
            etaRange = key;
            zRange   = "originalZRange";
          }

          // Build a canonical signature for de-duplication (ignore numeric fit values on purpose).
          const std::string axis = (isPhi ? "phi" : "eta");
          const std::string sig  = variantName + "|" + axis + "|" + etaRange + "|" + zRange;

          // Only write if this (variant,axis,etaRange,zRange) hasn't been written before in this run.
          if (!g_master_seen.count(sig)) {
            g_master_seen[sig] = 1;
            // space-delimited row
            mout << variantName << " "
                 << axis       << " "
                 << etaRange   << " "
                 << zRange     << " "
                 << fit.first  << " "
                 << fit.second << "\n";
          }
        }
      }
    }

    
  // ---------- Reusable renderer for ONE overlay (subset) ----------
  auto renderOverlay = [&](const std::string& sfx,
                           const std::map<std::string, std::vector<double>>& M,
                           const std::map<std::string, std::vector<double>>& Merr)
  {
    if (M.empty()) return;

    // y-range across the subset (±σ if available)
    double ymin = +1e9, ymax = -1e9;
    for (const auto& kv : M)
    {
      const auto& vals = kv.second;
      auto itErr = Merr.find(kv.first);
      const std::vector<double>* errs = (itErr == Merr.end()) ? nullptr : &itErr->second;
      for (size_t i=0;i<vals.size();++i)
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

    // For η overlays with z-only suffixes, relax autoscale slightly
    if (!isPhi && sfx.find("_zOnly") != std::string::npos) {
      const double padTop = 0.06 * std::fabs(ymax);
      const double padBot = 0.03 * std::fabs(ymin);
      yLo = std::max(0.0, ymin - padBot);
      yHi = ymax + padTop;
      if (yHi - yLo < 1e-4) { yLo -= 0.01; yHi += 0.01; }
    }

      // Auto-scale y-axis for η overlay saved as ..._zOnly_upto60.png
      // Add +40% extra top buffer relative to the existing pad.
      if (!isPhi && sfx.find("_zOnly_upto60") != std::string::npos) {
          const double baseTop = std::max(0.12 * std::fabs(ymax), 0.005);  // ≥0.005 headroom
          const double padTop  = baseTop * 1.95;                           // +40% extra
          const double padBot  = 0.04 * std::fabs(ymin);
          yLo = std::max(0.0, ymin - padBot);
          yHi = ymax + padTop;
          if (yHi - yLo < 1e-4) {                      // guard for near-flat cases
              yLo = std::max(0.0, ymin - 0.01);
              yHi = ymax + 0.1;
          }
      }

    // Keys classification for ordering
    bool allCoarseZ = true, allFineZ = true;
    for (const auto& kv : M) {
      if (kZPretty.find(kv.first)     == kZPretty.end())     allCoarseZ = false;
      if (kZFinePretty.find(kv.first) == kZFinePretty.end()) allFineZ   = false;
    }

    std::vector<std::string> desiredEta = {"etaCore","etaMid","etaEdge","fullEta","originalEta"};
    std::vector<std::string> keysOrdered, keysForTable;

    if (allFineZ) {
      for (const auto& z : kZFineOrder)
        if (M.find(z) != M.end()) keysOrdered.push_back(z);
      keysForTable = keysOrdered;
    } else if (allCoarseZ) {
      for (const auto& z : kZOrder)
        if (M.find(z) != M.end()) keysOrdered.push_back(z);
      keysForTable = keysOrdered;
    } else {
      for (const auto& k : desiredEta)
        if (k != "originalEta" && M.find(k) != M.end())
          keysOrdered.push_back(k);
      for (const auto& kv : M)
        if (kv.first != "originalEta" && !containsKey(keysOrdered, kv.first))
          keysOrdered.push_back(kv.first);

      keysForTable = keysOrdered;
      if (M.find("originalEta") != M.end() &&
          !containsKey(keysForTable, "originalEta")) {
        keysForTable.push_back("originalEta");
      }
    }

    // Canvas & frame
    gStyle->SetOptStat(0);
    TCanvas c("cbOverlay",
              isPhi ? "b_{#varphi}(E) overlay" : "b_{#eta}(E) overlay",
              900, 700);
    TH2F frame("frame","",100,xLo,xHi,100,yLo,yHi);
    frame.SetTitle(isPhi ? "best-fit  b_{#varphi}  vs  E;E  [GeV];b_{#varphi}"
                         : "best-fit  b_{#eta}     vs  E;E  [GeV];b_{#eta}");
    frame.SetStats(0);
    frame.Draw();

    const int colors[8] = { kBlack, kBlue+1, kRed+1, kGreen+2,
                            kMagenta+1, kOrange+7, kCyan+2, kGray+2 };
    const int    kMarkerStyle   = 20;
    const double kMarkerSize    = 1.2;

    TLegend leg(0.16,0.65,0.42,0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);

    // Terminal banner
    {
      const char* B = "\033[1m";
      const char* R = "\033[0m";
      std::ostringstream t;
      t << (isPhi ? "[PHI] " : "[ETA] ")
        << "b(E) = b0 + m * ln(E/E0),  E0 = " << E0 << " GeV";
      std::cout << "\n" << B << t.str() << R << "\n";
      std::cout << "Weighted least squares per variant: w_i = 1/sigma_i^2 (if provided) or 1 otherwise.\n"
                << "Parameter errors: unscaled if errors provided; scaled by sqrt(χ²/ndf) when unit weights are used.\n\n";
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
                << " | " << std::setw(10) << "χ²/ndf" << "\n";
      std::cout << std::string(
                     18+2 + 6+3 + 9+3 + 9+3 + 12+3 + 10+3 + 12+3 + 10+3 +
                     9+3 + 12+3 + 12+3 + 8+3 + 10+3 + 10, '-') << "\n";
      std::cout.setf(std::ios::fixed);
      std::cout << std::setprecision(6);
    }

    std::vector<std::unique_ptr<TGraphErrors>> keepPoints;
    std::vector<std::unique_ptr<TF1>>          keepFits;

    // CSV rows (legacy per-call CSV still produced)
    std::vector<std::string> csvRows;
    {
      std::ostringstream hdr;
      hdr << "variant,N,E_min,E_max,b0@E0,db0,m_lnE,dm_lnE,z_m,m_per_decade,Delta_b_fit,Delta_b_data,R2w,RMSEw,chi2_ndf";
      csvRows.push_back(hdr.str());
    }

    int idx = 0;
    for (const auto& kv : M)
    {
      const std::string& key  = kv.first;
      const auto&        vals = kv.second;
      if (vals.size() != ecenters.size()) { ++idx; continue; }

      auto itErr = Merr.find(key);
      std::vector<double> errs = (itErr != Merr.end()) ? itErr->second
                                                       : std::vector<double>(vals.size(), 0.0);
      if (errs.size() != vals.size()) errs.assign(vals.size(), 0.0);

      std::unique_ptr<TGraphErrors> g(new TGraphErrors(static_cast<int>(vals.size())));
      double Emin = +1e9, Emax = -1e9;
      double yMin = +1e9, yMax = -1e9;

      for (size_t i=0;i<vals.size();++i)
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

      // Determine if this is an η-only overlay
      bool hasAnyEtaLocal = false;
      for (const auto& kv2 : M) {
        if (kEtaPretty.find(kv2.first) != kEtaPretty.end()) { hasAnyEtaLocal = true; break; }
      }
      const bool etaOnlyOverlay = !(allCoarseZ || allFineZ) && hasAnyEtaLocal;

      auto colorForEtaKey = [&](const std::string& k)->Color_t {
        if (k == "etaCore") return kBlack;
        if (k == "etaMid")  return kBlue+1;
        if (k == "etaEdge") return kRed+1;
        if (k == "fullEta") return kGreen+2;
        return kGray+2;
      };
      auto legendTextForEtaKey = [&](const std::string& k)->const char* {
        if (k == "etaCore") return "|#eta| #leq 0.20";
        if (k == "etaMid")  return "0.20 < |#eta| #leq 0.70";
        if (k == "etaEdge") return "0.70 < |#eta| #leq 1.10";
        if (k == "fullEta") return "|#eta| #leq 1.10";
        return prettyName(k);
      };

      Color_t col = etaOnlyOverlay ? colorForEtaKey(key) : colors[idx % 8];
      bool drawThis = (allCoarseZ || allFineZ) || (key != "originalEta");

      if (drawThis) {
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetMarkerColor(col);
        g->SetLineColor(col);
        g->Draw("P SAME");

        const char* legTxt = etaOnlyOverlay ? legendTextForEtaKey(key) : prettyName(key);
        leg.AddEntry(g.get(), legTxt, "p");

        keepPoints.emplace_back(std::move(g));
      } else {
        keepPoints.emplace_back(std::move(g));
      }

      // Weighted least squares on x = ln(E/E0)
      double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
      int     Nuse = 0;
      std::vector<double> xi; xi.reserve(vals.size());
      std::vector<double> yi; yi.reserve(vals.size());
      std::vector<double> wi; wi.reserve(vals.size());
      bool haveYerrsAll = true;

      for (size_t i=0;i<vals.size();++i)
      {
        const double E  = ecenters[i];
        const double y  = vals[i];
        const double sy = errs[i];
        if (!(std::isfinite(E) && E>0.0 && std::isfinite(y))) continue;

        const double x = std::log(E / E0);
        double w = 1.0;
        if (std::isfinite(sy) && sy > 0.0) w = 1.0/(sy*sy); else haveYerrsAll = false;

        S   += w;  Sx  += w*x;  Sy  += w*y;
        Sxx += w*x*x;  Sxy += w*x*y;
        xi.push_back(x); yi.push_back(y); wi.push_back(w);
        ++Nuse;
      }

      double b0=std::numeric_limits<double>::quiet_NaN();
      double m =std::numeric_limits<double>::quiet_NaN();
      double db0=std::numeric_limits<double>::quiet_NaN();
      double dm =std::numeric_limits<double>::quiet_NaN();
      double chi2=0.0, R2w=0.0, RMSEw=0.0, redChi2=0.0;
      double dBfit=0.0, dBdata=0.0, z_m=0.0, m_per_dec=0.0;

      const double Delta = S*Sxx - Sx*Sx;
      if (Nuse >= 2 && Delta > 0.0)
      {
        m  = (S * Sxy - Sx * Sy) / Delta;
        b0 = (Sxx * Sy - Sx * Sxy) / Delta;

        const double var_m  =  S / Delta;
        const double var_b0 = Sxx / Delta;

        double ybar_w = Sy / S;
        double TSSw = 0.0;
        for (int i = 0; i < Nuse; ++i)
        {
          const double yhat = b0 + m * xi[i];
          const double r    = yi[i] - yhat;
          chi2  += wi[i] * r * r;
          TSSw  += wi[i] * (yi[i] - ybar_w) * (yi[i] - ybar_w);
        }
        dBdata = yMax - yMin;

        const int ndf = std::max(0, Nuse - 2);
        redChi2 = (ndf > 0) ? (chi2 / ndf) : 0.0;
        R2w     = (TSSw > 0) ? (1.0 - chi2 / TSSw) : 0.0;
        RMSEw   = (S > 0.0) ? std::sqrt(chi2 / S) : 0.0;

        const double scale = (!haveYerrsAll && ndf > 0) ? std::sqrt(redChi2) : 1.0;
        dm  = std::sqrt(std::max(0.0, var_m )) * scale;
        db0 = std::sqrt(std::max(0.0, var_b0)) * scale;

        z_m       = (dm > 0.0 ? m / dm : 0.0);
        m_per_dec = m * std::log(10.0);

        if (!xi.empty())
        {
          const double x_min = *std::min_element(xi.begin(), xi.end());
          const double x_max = *std::max_element(xi.begin(), xi.end());
          const double yhat_min = b0 + m * x_min;
          const double yhat_max = b0 + m * x_max;
          dBfit = yhat_max - yhat_min;
        }

        std::unique_ptr<TF1> f(new TF1(Form("f_%s_%d", key.c_str(), idx),
                                       "[0] + [1]*log(x/[2])", xLo, xHi));
        f->SetParameters(b0, m, E0);
        f->FixParameter(2, E0);
        f->SetLineColor(col);
        f->SetLineWidth(2);
        if (drawThis) f->Draw("SAME");
        keepFits.emplace_back(std::move(f));
      }

      // Print & legacy CSV line
      {
        std::ostringstream ln; ln.setf(std::ios::fixed); ln << std::setprecision(6);
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

        std::ostringstream csv;
        csv << key << ","
            << Nuse << ","
            << (std::isfinite(Emin)?Emin:0.0) << ","
            << (std::isfinite(Emax)?Emax:0.0) << ","
            << (std::isfinite(b0)?b0:0.0) << ","
            << (std::isfinite(dm)?dm:0.0) << ","
            << (std::isfinite(m)?m:0.0) << ","
            << (std::isfinite(dm)?dm:0.0) << ","
            << (std::isfinite(z_m)?z_m:0.0) << ","
            << (std::isfinite(m_per_dec)?m_per_dec:0.0) << ","
            << (std::isfinite(dBfit)?dBfit:0.0) << ","
            << (std::isfinite(dBdata)?dBdata:0.0) << ","
            << (std::isfinite(R2w)?R2w:0.0) << ","
            << (std::isfinite(RMSEw)?RMSEw:0.0) << ","
            << (std::isfinite(redChi2)?redChi2:0.0);
        csvRows.push_back(csv.str());
      }

      ++idx;
    }

    // Write per-call CSV next to the image (overlayOutDir)
    {
      const char* which = isPhi ? "phi" : "eta";
      std::string csvPath = overlayOutDir + "/bFitParams_" + std::string(which) + sfx + ".csv";
      std::ofstream fout(csvPath);
      if (fout) {
        for (const auto& s : csvRows) fout << s << "\n";
        fout.close();
        std::cout << "[b-fit] Wrote CSV: " << csvPath << "\n";
      } else {
        std::cerr << "[b-fit] WARNING: Could not write CSV to " << csvPath << "\n";
      }
    }

      // Save the overlay image to the decided folder
      {
        EnsureDir(overlayOutDir);

        // Draw legend first
        leg.Draw();  // ensure the legend (with z ranges) is painted

        // If saving inside .../etaAndzDepBlocks/<etaKey>, stamp the η-range on mid right
        if (!etaKeyFromPath.empty()) {
          std::string etaText = etaKeyFromPath;
          auto itEtaLbl = kEtaPretty.find(etaKeyFromPath);
          if (itEtaLbl != kEtaPretty.end()) etaText = itEtaLbl->second;

          TLatex latex;
          latex.SetNDC();          // x,y in normalized device coordinates
          latex.SetTextAlign(32);  // right-aligned, vertically centered
          latex.SetTextSize(0.045);
          // mid-right of canvas (tweak if you want it higher/lower)
          latex.DrawLatex(0.86, 0.55, etaText.c_str());
        }

        c.RedrawAxis();  // keep axes on top of drawn objects

        TString outMain = Form("%s/%s%s.png",
                               overlayOutDir.c_str(),
                               isPhi ? "bValuesPhiOverlay" : "bValuesEtaOverlay",
                               sfx.c_str());
        c.SaveAs(outMain);
      }

  }; // renderOverlay

  // ---------- image rendering logic ----------
  if (hasOnlyCoarse)
  {
    // Keep ONLY the ≤60 cm overlay for coarse z (remove any >60 overlays)
    std::map<std::string, std::vector<double>> M60, M60err;
    for (const auto& tag : std::vector<std::string>{"z00to10","z10to20","z20to30","z30to45","z45to60"})
    {
      auto it = byVar.find(tag); if (it != byVar.end()) {
        M60[tag] = it->second;
        auto itE = byVarErr.find(tag); if (itE != byVarErr.end()) M60err[tag] = itE->second;
      }
    }

    if (!M60.empty()) {
      // Suffix: ensure both z-only + upto60, and include per-eta tag when inside etaAndzDepBlocks/<etaKey>
      const std::string sfx = suffixLabelTag + "_zOnly_upto60";
      renderOverlay(sfx, M60, M60err);
    }

    return; // coarse case handled (no high-z overlay any more)
  }

  // Fine‑z or η overlays → single overlay image (with proper folder & labeling)
  {
    // For fine z-only, keep the "_zFine" suffix; add per-eta tag if in etaAndzDepBlocks/<etaKey>
    // For eta-only (global), keep baseSfx (usually ""), no per-eta tag, and write to etaDepOnlyBlockAna/
    const std::string sfx = (hasOnlyFine || inEtaAndZContext) ? (suffixLabelTag + baseSfx)
                                                              : (baseSfx);
    renderOverlay(sfx, byVar, byVarErr);
  }
}

// -----------------------------------------------------------------------------
// SaveIncidenceQA
//   1) from TH3: ⟨α_φ^{signed}⟩(z) & ⟨α_η^{signed}⟩(z) by |η| band and sign
//        FOLDERS under  <outBaseDir>/incidenceQA/ :
//          etaCorePositive/   etaCoreNegative/   etaCoreAbsolute/
//          etaMidPositive/    etaMidNegative/    etaMidAbsolute/
//          etaEdgePositive/   etaEdgeNegative/   etaEdgeAbsolute/
//          etaFullPositive/   etaFullNegative/   etaFullAbsolute/
//            └── z_0_10|20|30|45|60/  → alphaPhiEtaOverlay.png   (⟨αφ⟩ & ⟨αη⟩ overlaid)
//        Also in base folder:
//          alphaEtaPositiveEtaOverlay.png   (η bands + only)
//          alphaPhiPositiveEtaOverlay.png   (φ bands + only)
//   2) α spectra per energy (φ / η)             → incidence_alpha_spectra_phi|eta.png
//   3) ⟨sec α⟩(α) per energy + 1/cos overlay    → incidence_secAlpha_profiles_phi|eta.png
// -----------------------------------------------------------------------------
static void SaveIncidenceQA(TFile* fin, const std::string& outBaseDir)
{
  if (!fin || fin->IsZombie()) { std::cerr << "[incidenceQA] bad TFile\n"; return; }

  // ---------- base outdir -----------------------------------------------------
  const std::string outDir = outBaseDir + "/incidenceQA";
  gSystem->mkdir(outDir.c_str(), true);

  // ---------- small readers ---------------------------------------------------
  auto getH3 = [&](const char* n)->TH3F* {
    TH3F* h=nullptr; fin->GetObject(n, h);
    if (h) { h->SetDirectory(nullptr); std::cout << "[incidenceQA] found " << n << "\n"; }
    else   { std::cout << "[incidenceQA] MISSING " << n << "\n"; }
    return h;
  };
  auto scanByEnergy = [&](const char* base, double eLo, double eHi, const char* klass)->TObject*
  {
    const char* suff[3] = {"", "_range", "_disc"};
    for (int s=0; s<3; ++s) {
      TObject* obj=nullptr;
      obj = fin->Get(Form("%s_%.1f_%.1f%s", base, eLo, eHi, suff[s])); if (obj) return obj;
      obj = fin->Get(Form("%s_%g_%g%s",     base, eLo, eHi, suff[s])); if (obj) return obj;
      obj = fin->Get(Form("%s_%d_%d%s",     base, (int)std::lround(eLo), (int)std::lround(eHi), suff[s])); if (obj) return obj;
    }
    // full scan
    TIter it(fin->GetListOfKeys());
    std::regex re_with_suff(Form("^%s(?:_.*)?_([-0-9.]+)_([-0-9.]+)_(?:range|disc)$", base));
    std::regex re_no_suff  (Form("^%s(?:_.*)?_([-0-9.]+)_([-0-9.]+)$",                base));
    while (auto* o = it()) {
      auto* k = dynamic_cast<TKey*>(o); if (!k) continue;
      if (klass && std::string(k->GetClassName()) != klass) continue;
      const std::string nm = k->GetName();
      std::smatch m;
      bool ok = (std::regex_match(nm, m, re_with_suff) && m.size()==3)
             || (std::regex_match(nm, m, re_no_suff)   && m.size()==3);
      if (!ok) continue;
      const double lo = atof(m[1].str().c_str());
      const double hi = atof(m[2].str().c_str());
      if (std::fabs(lo-eLo) < 1e-3 && std::fabs(hi-eHi) < 1e-3) {
        TObject* obj = fin->Get(nm.c_str());
        if (obj) return obj;
      }
    }
    std::cout << Form("[incidenceQA] %s : no match for %.3f–%.3f\n", base, eLo, eHi) << "\n";
    return nullptr;
  };
  auto findH1 = [&](const char* base, double eLo, double eHi)->TH1F* {
    auto* o = scanByEnergy(base, eLo, eHi, "TH1F");
    if (!o) return nullptr;
    auto* h = dynamic_cast<TH1F*>(o);
    if (h) h->SetDirectory(nullptr);
    return h;
  };
  auto findProf = [&](const char* base, double eLo, double eHi)->TProfile* {
    auto* o = scanByEnergy(base, eLo, eHi, "TProfile");
    if (!o) return nullptr;
    auto* p = dynamic_cast<TProfile*>(o);
    if (p) p->SetDirectory(nullptr);
    return p;
  };

  // ---------- TH3 inputs (signed α on Z-axis) --------------------------------
  gStyle->SetOptStat(0);
  std::unique_ptr<TH3F> h3Phi(getH3("h3_alphaPhi_vsVz_vsEta"));
  std::unique_ptr<TH3F> h3Eta(getH3("h3_alphaEta_vsVz_vsEta"));
  if (!h3Phi && !h3Eta) { std::cerr << "[incidenceQA] no α(z,η) TH3 → abort\n"; return; }

  const double zRangeMin = -60.0, zRangeMax = +60.0;

  // eta band meta
  struct BandSpec { const char* tag; double lo; double hi; Color_t col; };
  // absolute variants (both signs) use |η| ranges; positive/negative pick a side
  const BandSpec bandCore { "etaCore", 0.00, 0.20, kGreen+2 };
  const BandSpec bandMid  { "etaMid" , 0.20, 0.70, kBlue+1  };
  const BandSpec bandEdge { "etaEdge", 0.70, 1.10, kRed+1   };
  const BandSpec bandFull { "etaFull", 0.00, 1.10, kBlack   };

  // folder recipes
  struct FolderRecipe { std::string name; const BandSpec* band; enum class Sign {Pos,Neg,Abs} sign; };
  std::vector<FolderRecipe> folders = {
    {"etaCorePositive", &bandCore, FolderRecipe::Sign::Pos},
    {"etaCoreNegative", &bandCore, FolderRecipe::Sign::Neg},
    {"etaCoreAbsolute", &bandCore, FolderRecipe::Sign::Abs},

    {"etaMidPositive",  &bandMid , FolderRecipe::Sign::Pos},
    {"etaMidNegative",  &bandMid , FolderRecipe::Sign::Neg},
    {"etaMidAbsolute",  &bandMid , FolderRecipe::Sign::Abs},

    {"etaEdgePositive", &bandEdge, FolderRecipe::Sign::Pos},
    {"etaEdgeNegative", &bandEdge, FolderRecipe::Sign::Neg},
    {"etaEdgeAbsolute", &bandEdge, FolderRecipe::Sign::Abs},

    {"etaFullPositive", &bandFull, FolderRecipe::Sign::Pos},
    {"etaFullNegative", &bandFull, FolderRecipe::Sign::Neg},
    {"etaFullAbsolute", &bandFull, FolderRecipe::Sign::Abs},
  };

  // z caps
  const std::vector<int> zCaps = {10,20,30,45,60};

  // helper: build profile ⟨α⟩(z) in a given |z|<cap, eta selection, signed α
  auto makeProfileVsZ = [&](TH3F* h3, const char* name,
                            double etaLoAbs, double etaHiAbs,
                            FolderRecipe::Sign sign, int zCap)->std::unique_ptr<TProfile>
  {
    if (!h3) return nullptr;

    auto* axZ  = h3->GetXaxis();
    auto* axEta= h3->GetYaxis();
    auto* axA  = h3->GetZaxis();

    const int xlo = std::max(1, axZ->FindFixBin(-zCap));
    const int xhi = std::min(axZ->GetNbins(), axZ->FindFixBin(+zCap));

    auto prof = std::make_unique<TProfile>(name, "", xhi-xlo+1,
                                           axZ->GetBinLowEdge(xlo), axZ->GetBinUpEdge(xhi));
    prof->GetXaxis()->SetTitle("z_{vtx}  (cm)");
    prof->GetYaxis()->SetTitle("#LT#alpha^{signed}#GT  [rad]");

    auto accum_eta = [&](int ix, int ylo, int yhi, double& sw, double& sm)
    {
      for (int iy = ylo; iy <= yhi; ++iy)
        for (int ia = 1; ia <= axA->GetNbins(); ++ia)
        {
          const double w = h3->GetBinContent(ix, iy, ia);
          if (w <= 0) continue;
          const double a = axA->GetBinCenter(ia); // signed
          sw += w; sm += w * a;
        }
    };

    // indices for η ranges
    const int yPosLo = std::max(1, axEta->FindFixBin(+etaLoAbs));
    const int yPosHi = std::min(axEta->GetNbins(), axEta->FindFixBin(+etaHiAbs));
    const int yNegLo = std::max(1, axEta->FindFixBin(-etaHiAbs));
    const int yNegHi = std::min(axEta->GetNbins(), axEta->FindFixBin(-etaLoAbs));

    for (int ix = xlo; ix <= xhi; ++ix)
    {
      double sumw=0.0, sum=0.0;
      if (sign == FolderRecipe::Sign::Pos)      accum_eta(ix, yPosLo, yPosHi, sumw, sum);
      else if (sign == FolderRecipe::Sign::Neg) accum_eta(ix, yNegLo, yNegHi, sumw, sum);
      else { // Abs → both sides
        accum_eta(ix, yPosLo, yPosHi, sumw, sum);
        accum_eta(ix, yNegLo, yNegHi, sumw, sum);
      }
      if (sumw > 0.0) prof->Fill(axZ->GetBinCenter(ix), sum/sumw, sumw);
    }

    return prof;
  };

    // draw overlay with a TITLE that encodes the η-selection and |z| cap
    auto drawOverlay = [&](TProfile* pPhi, TProfile* pEta,
                           const std::string& pngPath, int zCap,
                           const char* etaTitle)  // <- NEW: title fragment with η info
    {
      if (!pPhi && !pEta) return;

      // y autoscale with asymmetric padding
      double yMin = +1e9, yMax = -1e9;
      auto upd = [&](TProfile* p){
        if (!p) return;
        for (int i=1;i<=p->GetNbinsX();++i) {
          if (p->GetBinEntries(i) <= 0) continue;
          const double y = p->GetBinContent(i);
          yMin = std::min(yMin, y);
          yMax = std::max(yMax, y);
        }
      };
      upd(pPhi); upd(pEta);
      if (!(yMin < yMax)) { yMin = -0.02; yMax = +0.02; }
      const double span = (yMax - yMin);
      const double padLow  = 0.06 * span;  // small padding below
      const double padHigh = 0.22 * span;  // extra padding above
      yMin -= padLow;
      yMax += padHigh;

      TCanvas c("c_alpha_overlay","alpha overlay",900,600);
      c.SetLeftMargin(0.12); c.SetRightMargin(0.06);
      c.SetBottomMargin(0.12); c.SetTopMargin(0.08);

      // TITLE now includes η band/sign and |z| cap
      TH2F frame("frame",
                 Form("%s  ;z_{vtx} (cm) (|z|<%d);#LT#alpha^{signed}#GT [rad]",
                      etaTitle, zCap),
                 100, -zCap, +zCap, 100, yMin, yMax);
      frame.SetStats(0); frame.Draw();

      if (pPhi) { pPhi->SetLineColor(kBlue+1); pPhi->SetLineWidth(3); pPhi->SetMarkerStyle(0); pPhi->Draw("hist same"); }
      if (pEta) { pEta->SetLineColor(kRed+1 ); pEta->SetLineWidth(3); pEta->SetMarkerStyle(0); pEta->Draw("hist same"); }

      // compact legend
      TLegend lg(0.66,0.80,0.90,0.89);
      lg.SetBorderSize(0);
      lg.SetFillStyle(0);
      lg.SetTextSize(0.038);
      if (pPhi) lg.AddEntry(pPhi, "#LT#alpha_{#varphi}#GT(z)", "l");
      if (pEta) lg.AddEntry(pEta, "#LT#alpha_{#eta}#GT(z)"   , "l");
      lg.Draw();

      c.SaveAs(pngPath.c_str());
    };

  // produce per-folder / per-zcap overlays
  for (const auto& f : folders)
  {
    const std::string folderPath = outDir + "/" + f.name;
    gSystem->mkdir(folderPath.c_str(), true);

    for (int cap : zCaps)
    {
      const std::string capDir = folderPath + Form("/z_0_%d", cap);
      gSystem->mkdir(capDir.c_str(), true);

      std::unique_ptr<TProfile> pPhi, pEta;
      if (h3Phi) pPhi = makeProfileVsZ(h3Phi.get(),
                                       Form("pPhi_%s_cap%d", f.name.c_str(), cap),
                                       f.band->lo, f.band->hi, f.sign, cap);
      if (h3Eta) pEta = makeProfileVsZ(h3Eta.get(),
                                       Form("pEta_%s_cap%d", f.name.c_str(), cap),
                                       f.band->lo, f.band->hi, f.sign, cap);

        // Build the η-title fragment exactly as requested: +η, −η, or |η|, and include the |z| cap
        std::string etaTitle;
        if (f.sign == FolderRecipe::Sign::Pos) {
          // positive η:  η ∈ [lo, hi]
          etaTitle = Form("+#eta#in[%.2f, %.2f],  |z_{vtx}|<%d cm",
                          f.band->lo, f.band->hi, cap);
        } else if (f.sign == FolderRecipe::Sign::Neg) {
          // negative η:  η ∈ [−hi, −lo]
          etaTitle = Form("-#eta#in[%.2f, %.2f],  |z_{vtx}|<%d cm",
                          -f.band->hi, -f.band->lo, cap);
        } else {
          // absolute:  |η| ∈ [lo, hi]
          etaTitle = Form("|#eta|#in[%.2f, %.2f],  |z_{vtx}|<%d cm",
                          f.band->lo, f.band->hi, cap);
        }

        const std::string outPng = capDir + "/alphaPhiEtaOverlay.png";
        drawOverlay(pPhi.get(), pEta.get(), outPng, cap, etaTitle.c_str());
        std::cout << "[incidenceQA] wrote " << outPng << "\n";
    }
  }

    // base overlays: ABSOLUTE-|η| bands (Core/Mid/Edge) for each view; match eta*Absolute/
    auto makeAbsBandProfile = [&](TH3F* h3, const BandSpec& bs)->std::unique_ptr<TProfile>
    {
      if (!h3) return nullptr;
      return makeProfileVsZ(h3, Form("p_%s_Absolute_fullZ", bs.tag),
                            bs.lo, bs.hi, FolderRecipe::Sign::Abs, /*zCap=*/60);
    };
    // alphaEta_EtaOverlay.png — absolute-|η| bands to match eta*Absolute/
    // Show cumulative labels and add extra headroom at the top
    {
      std::unique_ptr<TProfile> pCore, pMid, pEdge;
      if (h3Eta) {
        pCore = makeAbsBandProfile(h3Eta.get(), bandCore);  // |η| ∈ [0.00,0.20]
        pMid  = makeAbsBandProfile(h3Eta.get(), bandMid );  // |η| ∈ (0.20,0.70]
        pEdge = makeAbsBandProfile(h3Eta.get(), bandEdge);  // |η| ∈ (0.70,1.10]
      }

      // y autoscale with bigger top padding
      double yMin=+1e9, yMax=-1e9;
      auto upd = [&](TProfile* p){
        if (!p) return;
        for (int i=1;i<=p->GetNbinsX();++i) {
          if (p->GetBinEntries(i)>0) {
            const double y=p->GetBinContent(i);
            yMin=std::min(yMin,y); yMax=std::max(yMax,y);
          }
        }
      };
      upd(pCore.get()); upd(pMid.get()); upd(pEdge.get());
      if (!(yMin<yMax)) { yMin=-0.02; yMax=+0.02; }
      const double span = (yMax-yMin>0? yMax-yMin : 0.20);
      const double padLow  = 0.04*span;
      const double padHigh = 0.30*span;   // more room on the upper y-axis
      yMin -= padLow;
      yMax += padHigh;

      TCanvas c("c_etaBands_abs","alphaEta |eta|-bands",900,600);
      c.SetLeftMargin(0.12); c.SetRightMargin(0.06);
      c.SetBottomMargin(0.12); c.SetTopMargin(0.08);

      TH2F frame("frameEta","|#eta| bands;z_{vtx} (cm)  (|z|<60);#LT#alpha_{#eta}^{signed}#GT [rad]",
                 100, zRangeMin, zRangeMax, 100, yMin, yMax);
      frame.SetStats(0); frame.Draw();

      if (pCore) { pCore->SetLineColor(bandCore.col); pCore->SetLineWidth(3); pCore->Draw("hist same"); }
      if (pMid ) { pMid ->SetLineColor(bandMid .col); pMid ->SetLineWidth(3); pMid ->Draw("hist same"); }
      if (pEdge) { pEdge->SetLineColor(bandEdge.col); pEdge->SetLineWidth(3); pEdge->Draw("hist same"); }

      // Legend: cumulative thresholds
      TLegend lg(0.60,0.70,0.92,0.90);
      lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.038);
      if (pCore) lg.AddEntry(pCore.get(), "|#eta| #leq 0.20", "l");
      if (pMid ) lg.AddEntry(pMid .get(), "|#eta| #leq 0.70", "l");
      if (pEdge) lg.AddEntry(pEdge.get(), "|#eta| #leq 1.10", "l");
      lg.Draw();

      const std::string outPng = outDir + "/alphaEta_EtaOverlay.png";
      c.SaveAs(outPng.c_str());
      std::cout << "[incidenceQA] wrote " << outPng << "\n";
    }
    
    // alphaPhi_EtaOverlay.png — absolute-|η| bands to match eta*Absolute/
    // Cumulative labels and extra top headroom
    {
      std::unique_ptr<TProfile> pCore, pMid, pEdge;
      if (h3Phi) {
        pCore = makeAbsBandProfile(h3Phi.get(), bandCore);
        pMid  = makeAbsBandProfile(h3Phi.get(), bandMid );
        pEdge = makeAbsBandProfile(h3Phi.get(), bandEdge);
      }

      // y autoscale with bigger top padding
      double yMin=+1e9, yMax=-1e9;
      auto upd = [&](TProfile* p){
        if (!p) return;
        for (int i=1;i<=p->GetNbinsX();++i) {
          if (p->GetBinEntries(i)>0) {
            const double y=p->GetBinContent(i);
            yMin=std::min(yMin,y); yMax=std::max(yMax,y);
          }
        }
      };
      upd(pCore.get()); upd(pMid.get()); upd(pEdge.get());
      if (!(yMin<yMax)) { yMin=-0.02; yMax=+0.02; }
      const double span = (yMax-yMin>0? yMax-yMin : 0.20);
      const double padLow  = 0.04*span;
      const double padHigh = 0.30*span;
      yMin -= padLow;
      yMax += padHigh;

      TCanvas c("c_phiBands_abs","alphaPhi |eta|-bands",900,600);
      c.SetLeftMargin(0.12); c.SetRightMargin(0.06);
      c.SetBottomMargin(0.12); c.SetTopMargin(0.08);

      TH2F frame("framePhi","|#eta| bands;z_{vtx} (cm)  (|z|<60);#LT#alpha_{#varphi}^{signed}#GT [rad]",
                 100, zRangeMin, zRangeMax, 100, yMin, yMax);
      frame.SetStats(0); frame.Draw();

      if (pCore) { pCore->SetLineColor(bandCore.col); pCore->SetLineWidth(3); pCore->Draw("hist same"); }
      if (pMid ) { pMid ->SetLineColor(bandMid .col); pMid ->SetLineWidth(3); pMid ->Draw("hist same"); }
      if (pEdge) { pEdge->SetLineColor(bandEdge.col); pEdge->SetLineWidth(3); pEdge->Draw("hist same"); }

      TLegend lg(0.60,0.70,0.92,0.90);
      lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.038);
      if (pCore) lg.AddEntry(pCore.get(), "|#eta| #leq 0.20", "l");
      if (pMid ) lg.AddEntry(pMid .get(), "|#eta| #leq 0.70", "l");
      if (pEdge) lg.AddEntry(pEdge.get(), "|#eta| #leq 1.10", "l");
      lg.Draw();

      const std::string outPng = outDir + "/alphaPhi_EtaOverlay.png";
      c.SaveAs(outPng.c_str());
      std::cout << "[incidenceQA] wrote " << outPng << "\n";
    }
  // =========================
  // (2) α spectra per energy
  // =========================
  {
    // NOTE: expects your scalar edges/labels in scope: E_edges[], N_E
    TCanvas cPhi("c_alpha_spec_phi","alpha phi spectra",1600,1200); cPhi.Divide(4,2);
    TCanvas cEta("c_alpha_spec_eta","alpha eta spectra",1600,1200); cEta.Divide(4,2);

    for (int i=0;i<N_E; ++i)
    {
      const double eLo = E_edges[i], eHi = E_edges[i+1];
      TH1F* h1 = findH1("h_alphaPhi", eLo, eHi);
      std::cout << Form("[incidenceQA] αφ spectrum  E=[%.1f,%.1f): %s\n", eLo,eHi, h1?"FOUND":"missing");
      if (h1) {
        cPhi.cd(i+1)->SetGrid();
        h1->SetLineColor(kBlue+1); h1->SetLineWidth(2);
        h1->GetXaxis()->SetTitle("#alpha_{#varphi} [rad]");
        h1->GetYaxis()->SetTitle("Counts");
        h1->SetTitle(Form("%.0f–%.0f GeV", eLo, eHi));
        h1->Draw("hist");
      }
      TH1F* h2 = findH1("h_alphaEta", eLo, eHi);
      std::cout << Form("[incidenceQA] αη spectrum  E=[%.1f,%.1f): %s\n", eLo,eHi, h2?"FOUND":"missing");
      if (h2) {
        cEta.cd(i+1)->SetGrid();
        h2->SetLineColor(kRed+1); h2->SetLineWidth(2);
        h2->GetXaxis()->SetTitle("#alpha_{#eta} [rad]");
        h2->GetYaxis()->SetTitle("Counts");
        h2->SetTitle(Form("%.0f–%.0f GeV", eLo, eHi));
        h2->Draw("hist");
      }
    }
    cPhi.SaveAs( (outDir + "/incidence_alpha_spectra_phi.png").c_str() );
    cEta.SaveAs( (outDir + "/incidence_alpha_spectra_eta.png").c_str() );
    std::cout << "[incidenceQA] wrote α spectra PNGs in " << outDir << "\n";
  }

  // ===========================================
  // (3) ⟨sec α⟩(α) + 1/cos overlays, per energy
  // ===========================================
  {
    const double aMin = 0.0, aMax = 0.30;
    TCanvas cpPhi("c_sec_prof_phi","<sec alpha_phi> vs alpha",1600,1200); cpPhi.Divide(4,2);
    TCanvas cpEta("c_sec_prof_eta","<sec alpha_eta> vs alpha",1600,1200); cpEta.Divide(4,2);

    for (int i=0;i<N_E; ++i) {
      const double eLo = E_edges[i], eHi = E_edges[i+1];

      if (TProfile* p1 = findProf("p_secAlpha_phi", eLo, eHi)) {
        cpPhi.cd(i+1)->SetGrid();
        p1->SetTitle(Form("%.0f–%.0f GeV", eLo, eHi));
        p1->GetXaxis()->SetTitle("#alpha_{#varphi} [rad]");
        p1->GetYaxis()->SetTitle("#LTsec #alpha_{#varphi}#GT");
        p1->GetXaxis()->SetRangeUser(aMin, aMax);
        p1->SetLineColor(kBlue+1); p1->SetMarkerStyle(20); p1->Draw("E1");
        TF1 f(Form("f_cos_phi_%d",i),"1.0/cos(x)", aMin, aMax);
        f.SetLineColor(kGray+2); f.SetLineStyle(2); f.Draw("same");
      }

      if (TProfile* p2 = findProf("p_secAlpha_eta", eLo, eHi)) {
        cpEta.cd(i+1)->SetGrid();
        p2->SetTitle(Form("%.0f–%.0f GeV", eLo, eHi));
        p2->GetXaxis()->SetTitle("#alpha_{#eta} [rad]");
        p2->GetYaxis()->SetTitle("#LTsec #alpha_{#eta}#GT");
        p2->GetXaxis()->SetRangeUser(aMin, aMax);
        p2->SetLineColor(kRed+1); p2->SetMarkerStyle(20); p2->Draw("E1");
        TF1 f(Form("f_cos_eta_%d",i),"1.0/cos(x)", aMin, aMax);
        f.SetLineColor(kGray+2); f.SetLineStyle(2); f.Draw("same");
      }
    }
    cpPhi.SaveAs( (outDir + "/incidence_secAlpha_profiles_phi.png").c_str() );
    cpEta.SaveAs( (outDir + "/incidence_secAlpha_profiles_eta.png").c_str() );
    std::cout << "[incidenceQA] wrote <sec α> profile PNGs in " << outDir << "\n";
  }

    // ================================
    // (4) 〈αφ^fold〉 vs E per |η| band  with  φ(E)=a−b·lnE fit
    // ================================
    {
      struct Band { const char* base; const char* tag; Color_t col; };
      // base prefixes must match your booking names:
      //   h_alphaPhi_sgn_mean_fold_core_%s_%s, _mid_, _edge_
      std::vector<Band> bands = {
        {"h_alphaPhi_sgn_mean_fold_core", "core (|#eta|#leq0.20)",  kGreen+2},
        {"h_alphaPhi_sgn_mean_fold_mid" , "mid  (0.20<|#eta|#leq0.70)", kBlue+1},
        {"h_alphaPhi_sgn_mean_fold_edge", "edge (0.70<|#eta|#leq1.10)", kRed+1}
      };

      const double Emin = E_edges[0];
      const double Emax = E_edges[N_E];

      for (const auto& B : bands)
      {
        // Gather points (E_center , <alpha_phi^fold>) and mean errors
        std::vector<double> vx, vy, vey;
        vx.reserve(N_E); vy.reserve(N_E); vey.reserve(N_E);

        for (int i = 0; i < N_E; ++i)
        {
          const double eLo = E_edges[i], eHi = E_edges[i+1];
          TH1F* h = findH1(B.base, eLo, eHi);            // folded, banded
          if (!h || h->GetEntries() <= 0) {
            std::cout << Form("[incidenceQA] %s: missing/empty for %.1f–%.1f GeV\n",
                              B.base, eLo, eHi) << "\n";
            continue;
          }
          const double eCtr = 0.5*(eLo + eHi);
          vx.push_back(eCtr);
          vy.push_back(h->GetMean());                    // 〈αφ^fold〉
          vey.push_back(h->GetMeanError());              // weight by mean error if available
        }

        if (vx.empty()) continue;

        // Canvas & graph (use errors if present)
        TCanvas c(Form("c_alphaPhi_mean_vs_E_%s", B.tag),
                  Form("<alpha_phi^{fold}> vs E — %s", B.tag), 900, 600);
        c.SetLeftMargin(0.12); c.SetRightMargin(0.06);
        c.SetBottomMargin(0.12); c.SetTopMargin(0.08);
        c.SetGrid();

        // choose weighted or unweighted graph depending on non-zero errors
        bool useErrors = false;
        for (auto ey : vey) { if (ey > 0) { useErrors = true; break; } }

        std::unique_ptr<TGraph> gr;
        if (useErrors) {
          auto* ge = new TGraphErrors((int)vx.size());
          for (int i=0;i<(int)vx.size();++i) {
            ge->SetPoint(i, vx[i], vy[i]);
            ge->SetPointError(i, 0.0, vey[i]);          // x errors = 0
          }
          gr.reset(ge);
        } else {
          gr.reset(new TGraph((int)vx.size(), vx.data(), vy.data()));
        }

        gr->SetTitle(Form("<#alpha_{#varphi}^{fold}> vs E — %s;E_{#gamma} [GeV];<#alpha_{#varphi}^{fold}> [rad]", B.tag));
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.2);
        gr->SetLineColor(B.col);
        gr->SetMarkerColor(B.col);
        gr->Draw("AP");  // A: axes, P: points  (ROOT drawing basics)  // see manual
        // (ROOT graphics manual)  // citation in text

        // Fit model: a - b*log(E); fit over [Emin,Emax]
        TF1 fFit("f_alpha_fold", "[0] - [1]*log(x)", Emin, Emax);   // TF1 expression
        fFit.SetParameters(0.0, 0.0);                                // initial guesses are fine
        fFit.SetLineColor(kGray+2); fFit.SetLineStyle(2); fFit.SetLineWidth(3);

        // use quiet fit; for TGraphErrors, y-errors act as weights if present
        if (useErrors) gr->Fit(&fFit, "Q"); else gr->Fit(&fFit, "Q");

        const double a = fFit.GetParameter(0);
        const double b = fFit.GetParameter(1);
        const double ea = fFit.GetParError(0);
        const double eb = fFit.GetParError(1);

        fFit.Draw("same");

        // Annotate a,b on the plot  (TLatex usage; math text)
        TLatex tx;
        tx.SetNDC(true); tx.SetTextSize(0.040);
        tx.DrawLatex(0.15, 0.88, Form("#font[62]{Fit: } #alpha_{#varphi}^{fold}(E)=a-b\\,\\ln E"));
        tx.DrawLatex(0.15, 0.82, Form("a = %.4g #pm %.2g", a, ea));
        tx.DrawLatex(0.15, 0.76, Form("b = %.4g #pm %.2g", b, eb));

        const std::string outPng = outDir + "/" + std::string("alphaPhi_mean_vs_E_")
                                 + std::string(B.base).substr(std::string("h_alphaPhi_sgn_mean_fold_").size())
                                 + ".png";
        c.SaveAs(outPng.c_str());
        std::cout << "[incidenceQA] wrote " << outPng << "\n";
      }
    }
}




// ======================== Step-4: build δ(E,α) and (c0,m) ========================
// ---- kernel fit with fixed beff → returns delta ± edelta ----
// Wrapped, shifted Jacobian kernel (beff fixed)
static double shiftedWrapped(double *x, double *p) {
  const double Norm = p[0], d = p[1], b = p[2];
  if (b <= 0) return 0.0;
  const double S = std::sinh(1.0/(2.0*b));
  auto K = [&](double u)->double {
    if (u < -0.5 || u > 0.5) return 0.0;
    return (2.0*b)/std::sqrt(1.0 + 4.0*u*u*S*S);
  };
  const double X = x[0];
  // measured X wraps by ±1 tower
  return Norm * ( K(X - d) + K(X - d + 1.0) + K(X - d - 1.0) );
}

static bool FitDeltaFixedBeff(TH1D* hX, double beff, double& delta, double& edelta,
                              double& chi2ndf, double& pval, int rebin=2)
{
  delta = edelta = chi2ndf = pval = std::numeric_limits<double>::quiet_NaN();
  if (!hX || hX->GetEntries() < 50 || beff <= 0) return false;

  // Gentle rebin for stable likelihood fits (QA-only knob)
  if (rebin>1) hX->Rebin(rebin);

  // Robust seed for delta: median of X (clamped)
  double xq[3] = {0.25, 0.50, 0.75}, q[3];
  hX->GetQuantiles(3, q, xq);
  double d0 = std::clamp(q[1], -0.20, +0.20);

  TF1 f("fShiftWrap", shiftedWrapped, -0.5, 0.5, 3);
  f.SetParNames("Norm","delta","beff");
  f.SetParameters(std::max(1.0, hX->GetMaximum()), d0, beff);
  f.SetParLimits(0, 1e-6, 1e12);
  f.SetParLimits(1, -0.25, +0.25);
  f.FixParameter(2, beff);
  f.SetNpx(1200);

  // Use likelihood for sparse slices, quiet, store result
  TFitResultPtr r = hX->Fit(&f, "QNRLS");  // Q=quiet, N=no draw, R=range, L=LLH, S=store
  if (!r.Get()) return false;
  const bool ok = r->IsValid() && r->Status()==0;
  if (!ok) return false;

  delta   = f.GetParameter(1);
  edelta  = f.GetParError(1);
  chi2ndf = (r->Ndf()>0 ? r->Chi2()/r->Ndf() : std::numeric_limits<double>::quiet_NaN());
  pval    = r->Prob();
  return std::isfinite(delta) && std::isfinite(edelta);
}

// ---- robust name finder for 2D hists: h2_Xmeas{Phi,Eta}_vsAlpha_* ----
// Accepts both "..._<lo>_<hi>_range" and "..._<lo>_<hi>_disc" (and integer/float forms)
static TH2F* FindXmeasVsAlpha(TFile* fin, const char* base, double eLo, double eHi)
{
  // 1) Try explicit candidates first (fast path)
  const char* suff[2] = {"range","disc"};
  for (int s=0; s<2; ++s)
  {
    TH2F* h=nullptr;
    // float with one decimal
    fin->GetObject(Form("%s_%.1f_%.1f_%s", base, eLo, eHi, suff[s]), h); if (h) return h;
    // generic %g
    fin->GetObject(Form("%s_%g_%g_%s",     base, eLo, eHi, suff[s]), h); if (h) return h;
    // integer fallback (2_4 etc)
    fin->GetObject(Form("%s_%d_%d_%s",     base, (int)std::lround(eLo), (int)std::lround(eHi), suff[s]), h); if (h) return h;
  }

  // 2) Finally, scan the file: accept any suffix AFTER the hi bound
  //    Pattern:  ^base(_anything)?_<lo>_<hi>_(range|disc)$
  TIter next(fin->GetListOfKeys());
  std::regex re(Form("^%s.*_([-0-9.]+)_([-0-9.]+)_(range|disc)$", base));
  while (TObject* obj = next()) {
    auto* k = dynamic_cast<TKey*>(obj); if (!k) continue;
    if (std::string(k->GetClassName()) != "TH2F") continue;
    std::smatch m; std::string nm = k->GetName();
    if (!std::regex_match(nm, m, re) || m.size() != 4) continue;
    const double lo = atof(m[1].str().c_str());
    const double hi = atof(m[2].str().c_str());
    if (std::fabs(lo - eLo) < 1e-3 && std::fabs(hi - eHi) < 1e-3) {
      TH2F* h=nullptr; fin->GetObject(nm.c_str(), h);
      if (h) return h;
    }
  }
  return nullptr;
}

// ---- load energy-only b(E) for originalEta × originalZRange out of your bValues.txt ----
static bool Load_bE_fromTxt(const std::string& path,
                            const char* etaSel, const char* zSel,
                            std::vector<double>& eLo,
                            std::vector<double>& eHi,
                            std::vector<double>& bphi,
                            std::vector<double>& beta)
{
  eLo.clear(); eHi.clear(); bphi.clear(); beta.clear();
  std::ifstream fin(path);
  if (!fin) { std::cerr << "[deltaFit] cannot open bValues file: " << path << "\n"; return false; }

  // Columns: PHI/ETA  [eLo,eHi)  bval  (opt etaLabel) (opt zLabel)
  std::regex rx(R"(^\s*(PHI|ETA)\s*\[\s*([0-9]*\.?[0-9]+)\s*,\s*([0-9]*\.?[0-9]+)\)\s*([0-9]*\.?[0-9]+)(?:\s+(\S+))?(?:\s+(\S+))?)",
                std::regex::icase);
  std::map<std::pair<double,double>, double> mapPhi, mapEta;
  std::string line;
  auto lower = [](std::string s){ for(char& c:s) c=std::tolower(c); return s; };

  while (std::getline(fin, line)) {
    std::smatch m; if (!std::regex_search(line, m, rx)) continue;
    const std::string what = lower(m[1]);
    const double lo = atof(m[2].str().c_str());
    const double hi = atof(m[3].str().c_str());
    const double bv = atof(m[4].str().c_str());
    const std::string labA = m[5].matched ? lower(m[5]) : "originaleta";
    const std::string labB = m[6].matched ? lower(m[6]) : "originalzrange";
    // accept only our target view & z tags
    const bool okEta = (labA==lower(etaSel) || labB==lower(etaSel));
    const bool okZ   = (labA==lower(zSel)  || labB==lower(zSel));
    if (!(okEta && okZ)) continue;

    if (what=="phi") mapPhi[{lo,hi}] = bv;
    else             mapEta[{lo,hi}] = bv;
  }
  if (mapPhi.empty() || mapEta.empty()) {
    std::cerr << "[deltaFit] ERROR: did not find energy-only rows for "
              << etaSel << " × " << zSel << " in " << path << "\n";
    return false;
  }
  for (auto& kv : mapPhi) {
    eLo.push_back(kv.first.first);
    eHi.push_back(kv.first.second);
    bphi.push_back(kv.second);
    auto it = mapEta.find(kv.first);
    beta.push_back( (it==mapEta.end() ? 0.15 : it->second) );
  }
  return true;
}

// ---- one-parameter fit δ vs sinα -> c1(E) with intercept constrained to 0 ----
static bool FitSlopeThroughOrigin(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& ey,
                                  double& c1, double& ec1)
{
  if (x.size() < 3) return false;

  TGraphErrors g(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    g.SetPoint(i, x[i], y[i]);
    g.SetPointError(i, 0.0, (ey[i] > 0 ? ey[i] : 1e-3));
  }

  TF1 f("f1", "[0]*x", 0, 1);
  // Return a TFitResult (needs 'S'); quiet & no drawing
  TFitResultPtr r = g.Fit(&f, "QNR S");
  if (int(r) != 0) return false;

  c1  = f.GetParameter(0);
  ec1 = f.GetParError(0);
  return std::isfinite(c1) && std::isfinite(ec1);
}



static void BuildDeltaAndWriteFits(TFile* fin,
                                   const std::string& outBaseDir,
                                   const std::vector<std::pair<double,double>>& eEdges,
                                   const std::string& bValuesTxt,
                                   const char* etaSel   = "originalEta",
                                   const char* zSel     = "originalZRange",
                                   const double  E0     = 3.0)
{
  // 0) read b(E) for the requested view (legacy: originalEta × originalZRange)
  std::vector<double> eLo, eHi, bphi, beta;
  if (!Load_bE_fromTxt(bValuesTxt, etaSel, zSel, eLo, eHi, bphi, beta)) return;

  // output directory for all Step-4 artifacts
  const std::string outDir = outBaseDir + "/deltaFit";
  gSystem->mkdir(outDir.c_str(), true);

  // plain tables (no plots except per-α slice QA)
  std::ofstream fPhi(outDir + "/delta_table_phi.txt");
  std::ofstream fEta(outDir + "/delta_table_eta.txt");

  // also build c1(E) numbers (no PNGs)
  std::vector<double> eCtr, c1phi, c1phiErr, c1eta, c1etaErr;

  gStyle->SetOptStat(0);

  // helper: find secα profile with either ..._range or ..._disc endings
  auto findSecProf = [&](const char* base, double lo, double hi)->TProfile*
  {
    const char* suff[2] = {"range","disc"};
    for (int s=0; s<2; ++s) {
      TProfile* p=nullptr;
      fin->GetObject(Form("%s_%.1f_%.1f_%s", base, lo, hi, suff[s]), p); if (p) return p;
      fin->GetObject(Form("%s_%g_%g_%s",     base, lo, hi, suff[s]), p); if (p) return p;
      fin->GetObject(Form("%s_%d_%d_%s",     base, (int)std::lround(lo), (int)std::lround(hi), suff[s]), p); if (p) return p;
    }
    // brute-force fallback
    TIter it(fin->GetListOfKeys());
    std::regex re(Form("^%s.*_([-0-9.]+)_([-0-9.]+)_(range|disc)$", base));
    while (TObject* obj = it()) {
      auto* k = dynamic_cast<TKey*>(obj); if (!k) continue;
      if (std::string(k->GetClassName()) != "TProfile") continue;
      std::smatch m; std::string nm = k->GetName();
      if (!std::regex_match(nm, m, re) || m.size()!=4) continue;
      const double L = atof(m[1].str().c_str());
      const double H = atof(m[2].str().c_str());
      if (std::fabs(L-lo)<1e-3 && std::fabs(H-hi)<1e-3) {
        TProfile* p=nullptr; fin->GetObject(nm.c_str(), p);
        if (p) return p;
      }
    }
    return nullptr;
  };

    // helper: accumulate δ points per (E, α) slice (no per-slice PNGs)
    // xsOut/ysOut/eysOut: vectors indexed by energy slice; each holds that E’s {sinα, δ, σδ} arrays
    auto accumulate_slices = [&](TH2F* h2, TProfile* pSec, size_t iE, bool isPhi,
                                 double lo, double hi, const std::vector<double>& bEvec,
                                 std::vector<std::vector<double>>& xsOut,
                                 std::vector<std::vector<double>>& ysOut,
                                 std::vector<std::vector<double>>& eysOut)
    {
      if (!h2) return;

      const int nA     = h2->GetNbinsX();
      const double Ei  = 0.5*(lo+hi);
      auto& xs  = xsOut[iE];
      auto& ys  = ysOut[iE];
      auto& eys = eysOut[iE];

      for (int j = 1; j <= nA; ++j)
      {
        const double aCtr = h2->GetXaxis()->GetBinCenter(j);
        std::unique_ptr<TH1D> hX( h2->ProjectionY(Form("%s_e%02zu_a%03d", isPhi?"py_phi":"py_eta", iE, j), j, j, "e") );
        if (!hX || hX->GetEntries() < 50) continue;

        // beff(E, α) = b(E) * <sec α>_bin (or analytic sec)
        const double bE   = (iE < bEvec.size() ? bEvec[iE] : 0.15);
        const double secA = (pSec && pSec->GetBinEntries(j)>0) ? pSec->GetBinContent(j)
                           : 1.0/std::max(1e-6, std::cos(aCtr));
        const double beff = std::max(1e-6, bE*secA);

        // wrapped, shifted kernel (beff fixed) with robust seeding
        double delta=std::numeric_limits<double>::quiet_NaN();
        double edelta=delta, chi2ndf=delta, pval=delta;
        if (!FitDeltaFixedBeff(hX.get(), beff, delta, edelta, chi2ndf, pval, /*rebin=*/2)) continue;

        // write the table line (phi/eta share same format)
        if (isPhi) fPhi << Form("%7.3f %7.3f  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                                lo, hi, aCtr, secA, beff, delta, edelta);
        else       fEta << Form("%7.3f %7.3f  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                                lo, hi, aCtr, secA, beff, delta, edelta);

        // accumulate point
        const double s = std::sin(aCtr);
        xs .push_back(s);
        ys .push_back(delta);
        eys.push_back(edelta > 0 ? edelta : 1e-3);
      }

      // also stash the energy center for later write_master
      eCtr.push_back(Ei);
    };

    // compact page (4×2) of δ vs sinα for all energies, with through-origin fit per panel
    auto plot_delta_vs_sinalpha_grid =
      [&](const char* outPng,
          const std::vector<std::pair<double,double>>& edges,
          const std::vector<std::vector<double>>& xsByE,
          const std::vector<std::vector<double>>& ysByE,
          const std::vector<std::vector<double>>& eysByE,
          bool isPhi,
          std::vector<double>& c1Out,
          std::vector<double>& c1ErrOut)
    {
      const int NE = static_cast<int>(edges.size());
      TCanvas c("c_delta_grid", isPhi ? "#delta_{#varphi} vs sin#alpha" : "#delta_{#eta} vs sin#alpha", 1600, 1200);
      c.Divide(4,2);
      c.SetFillColor(0);

      c1Out.clear(); c1ErrOut.clear(); c1Out.reserve(NE); c1ErrOut.reserve(NE);

      for (int i=0; i<NE; ++i) {
        c.cd(i+1)->SetGrid();
        const auto& x  = xsByE[i];
        const auto& y  = ysByE[i];
        const auto& ey = eysByE[i];
        if (x.size()<3) { TLatex t; t.SetNDC(); t.DrawLatex(0.2,0.5,"insufficient points"); continue; }

          // Keep objects alive until after SaveAs
          static std::vector<std::unique_ptr<TGraphErrors>> s_keepGraphs;
          static std::vector<std::unique_ptr<TF1>>          s_keepFits;

          auto gptr_up = std::make_unique<TGraphErrors>(static_cast<int>(x.size()));
          TGraphErrors* gptr = gptr_up.get();
          for (size_t k=0;k<x.size();++k) {
            gptr->SetPoint(static_cast<int>(k), x[k], y[k]);
            gptr->SetPointError(static_cast<int>(k), 0.0, (ey[k]>0?ey[k]:1e-3));
          }
          gptr->SetMarkerStyle(20);
          gptr->SetMarkerSize(1.0);
          gptr->SetMarkerColor(kBlack);
          gptr->SetLineColor(kBlack);
          gptr->SetTitle(Form("%.0f–%.0f GeV;sin#alpha;#delta",
                              edges[i].first, edges[i].second));
          gptr->Draw("AP");

          auto f_up = std::make_unique<TF1>(Form("f_lin_%d", i), "[0]*x", 0, 1);
          TF1* fptr = f_up.get();
          TFitResultPtr r = gptr->Fit(fptr, "QNR S");
          double c1  = fptr->GetParameter(0);
          double ec1 = fptr->GetParError(0);

          // retain until canvas is saved
          s_keepGraphs.emplace_back(std::move(gptr_up));
          s_keepFits.emplace_back(std::move(f_up));


        // quick R^2 (unweighted)
        double ybar=0; for (double yi: y) ybar+=yi; ybar/=y.size();
        double ss_tot=0, ss_res=0;
        for (size_t k=0;k<y.size();++k) {
          const double yhat = c1*x[k];
          ss_res += (y[k]-yhat)*(y[k]-yhat);
          ss_tot += (y[k]-ybar)*(y[k]-ybar);
        }
        const double R2 = (ss_tot>0? 1.0-ss_res/ss_tot : 0.0);

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.045);
        tl.DrawLatex(0.15,0.85,Form("c_{1}=%.4f #pm %.4f", c1, ec1));
        tl.DrawLatex(0.15,0.78,Form("R^{2}=%.3f", R2));

        c1Out.push_back(c1);
        c1ErrOut.push_back(ec1);
      }
        c.Modified();
        c.Update();
        c.SaveAs(outPng);

        
    };


    // 1) loop energy slices — accumulate δ vs sinα points (no per-slice PNGs)
    const size_t NE = eEdges.size();
    std::vector<std::vector<double>> xs_phi(NE), ys_phi(NE), eys_phi(NE);
    std::vector<std::vector<double>> xs_eta(NE), ys_eta(NE), eys_eta(NE);

    for (size_t i = 0; i < NE; ++i)
    {
      const double lo = eEdges[i].first;
      const double hi = eEdges[i].second;

      // pull the 2D maps
      TH2F* h2phi = FindXmeasVsAlpha(fin, "h2_XmeasPhi_vsAlpha", lo, hi);
      TH2F* h2eta = FindXmeasVsAlpha(fin, "h2_XmeasEta_vsAlpha", lo, hi);
      if (!h2phi && !h2eta) continue;

      // matching <sec α> profiles (or analytic sec if missing)
        TProfile* pSecPhi = findSecProf("p_secAlpha_phi", lo, hi);
        TProfile* pSecEta = findSecProf("p_secAlpha_eta", lo, hi);

      // accumulate δ points
      if (h2phi) accumulate_slices(h2phi, pSecPhi, i, /*isPhi=*/true,  lo, hi, bphi, xs_phi, ys_phi, eys_phi);
      if (h2eta) accumulate_slices(h2eta, pSecEta, i, /*isPhi=*/false, lo, hi, beta, xs_eta, ys_eta, eys_eta);
    }

    // close the per-slice tables
    fPhi.close();
    fEta.close();

    // 2) compact δ vs sinα pages (one per axis) + fill c1(E) arrays
    plot_delta_vs_sinalpha_grid( (outDir + "/delta_vs_sinAlpha_phi.png").c_str(),
                                 eEdges, xs_phi, ys_phi, eys_phi, /*isPhi=*/true,
                                 c1phi, c1phiErr);

    plot_delta_vs_sinalpha_grid( (outDir + "/delta_vs_sinAlpha_eta.png").c_str(),
                                 eEdges, xs_eta, ys_eta, eys_eta, /*isPhi=*/false,
                                 c1eta, c1etaErr);

    // eCtr already filled by accumulate_slices (Ei per energy)
    std::cout << "[deltaFit] wrote δ vs sinα grids and δ tables to " << outDir << "\n";
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
  const Color_t C_CLUSCP_EA_eta  = kGreen+2;
  const Color_t C_PDCraw         = kBlack;
  const Color_t C_PDCcor         = kGray+2;
  const Color_t C_CP_BVALS       = kOrange+7;
   // NEW colors for the new EA variants
  const Color_t C_CLUSCP_EA_vzEta = kViolet+1;  // |z| + |η| + E
  const Color_t C_CLUSCP_EA_Einc  = kCyan+2;    // E-only + incident-angle

  const std::vector<VarDef> kVariants = {
    // key                      pretty                                             phiPrefix                                                etaPrefix                                                 color
    {"CLUSRAW",                "No Correction",                                   "h_phi_diff_cpRaw_",                                      "h_eta_diff_cpRaw_",                                       C_CLUSRAW},
    {"CLUSCP",                 "Constant-b (legacy)",                             "h_phi_diff_cpCorr_",                                     "h_eta_diff_cpCorr_",                                      C_CLUSCP},
    {"CLUSCP_EA",              "Energy-Dep b",                                    "h_phi_diff_cpCorrEA_fitEnergyOnly_",                     "h_eta_diff_cpCorrEA_fitEnergyOnly_",                      C_CLUSCP_EA},
    {"CLUSCP_EA_vzDep",        "Energy + |v_{z}|-dep b",                          "h_phi_diff_cpCorrEA_fitZVTXDep_",                        "h_eta_diff_cpCorrEA_fitZVTXDep_",                         C_CLUSCP_EA_vz}, // (ZDEP)

    {"CLUSCP_EA_etaDep",       "Energy + |#eta|-dep b",                           "h_phi_diff_cpCorrEA_fitEtaDep_",                         "h_eta_diff_cpCorrEA_fitEtaDep_",                          C_CLUSCP_EA_eta},
    {"CLUSCP_EA_vzEtaDep",     "Energy + |z_{vtx}|+|#eta|-dep b",                 "h_phi_diff_cpCorrEA_fitZVTXEtaDep_",                     "h_eta_diff_cpCorrEA_fitZVTXEtaDep_",                      C_CLUSCP_EA_vzEta},
    {"CLUSCP_EA_EandIncident", "Energy-only + incident-angle",                    "h_phi_diff_cpCorrEA_fitEnergyOnly_AndIncidentAngle_",    "h_eta_diff_cpCorrEA_fitEnergyOnly_AndIncidentAngle_",     C_CLUSCP_EA_Einc},
    {"PDCraw",                 "2x2 Block-Local No Correction",                   "h_phi_diff_raw_",                                        "h_eta_diff_raw_",                                         C_PDCraw},
    {"PDCcor",                 "2x2 Block-Local With Corr",                       "h_phi_diff_corr_",                                       "h_eta_diff_corr_",                                        C_PDCcor},
    {"CLUScpDirectBvals",      "2x2 block b-values applied to No Corr Prod",      "h_phi_diff_cpBcorr_",                                    "h_eta_diff_cpBcorr_",                                     C_CP_BVALS}
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

    // ---------- auto x-range tightening (20% buffer around content) ----------
    inline void tighten_x_range(TH1* h, double frac = 0.20) {
      if (!h) return;
      TAxis* ax = h->GetXaxis();
      if (!ax) return;
      const int n = ax->GetNbins();
      int lo = 1, hi = n;
      while (lo <= hi && h->GetBinContent(lo) <= 0.0) ++lo;
      while (hi >= lo && h->GetBinContent(hi) <= 0.0) --hi;
      if (lo > hi) return; // no visible content
      const double xminC = ax->GetBinLowEdge(lo);
      const double xmaxC = ax->GetBinUpEdge(hi);
      const double w     = std::max(1e-12, xmaxC - xminC);
      const double pad   = frac * w;
      // clamp to the histogram’s native axis limits
      const double newMin = std::max(ax->GetXmin(), xminC - pad);
      const double newMax = std::min(ax->GetXmax(), xmaxC + pad);
      if (newMax > newMin) ax->SetRangeUser(newMin, newMax);
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
    if (key=="CLUSRAW")                 return "ETiltVariant::CLUS_RAW";
    if (key=="CLUSCP")                  return "ETiltVariant::CLUS_CP";
    if (key=="CLUSCP_EA_vzEtaDep")      return "ETiltVariant::CLUS_CP_EA_FIT_ZDEP_ETADEP";
    if (key=="CLUSCP_EA_vzDep")         return "ETiltVariant::CLUS_CP_EA_FIT_ZDEP";
    if (key=="CLUSCP_EA_etaDep")        return "ETiltVariant::CLUS_CP_EA_FIT_ETADEP";
    if (key=="CLUSCP_EA_EandIncident")  return "ETiltVariant::CLUS_CP_EA_FIT_EandIncident";
    if (key=="CLUSCP_EA")               return "ETiltVariant::CLUS_CP_EA_FIT_EONLY";
    if (key=="CLUScpDirectBvals")       return "ETiltVariant::CLUS_BCORR";
    if (key=="PDCraw")                  return "ETiltVariant::PDC_RAW";
    if (key=="PDCcor")                  return "ETiltVariant::PDC_CORR";
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

      // Match the 2×4 mosaic behavior: tighten X to content (with 20% pad)
      tighten_x_range(h, 0.20);

      // Compute Y max only over the currently visible X range
      int first = h->GetXaxis()->GetFirst();
      int last  = h->GetXaxis()->GetLast();
      double m  = 0.0;
      for (int b = first; b <= last; ++b) m = std::max(m, h->GetBinContent(b));
      const double ymax = 1.25 * std::max(1.0, m);
      h->GetYaxis()->SetRangeUser(0.0, ymax);

      h->Draw("E1");

      const auto& e = eSlices[i];
      TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.045); hdr.SetTextAlign(13);
      // Avoid Unicode em-dash (TLatex can't render it reliably)
      hdr.DrawLatex(0.12, 0.94, Form("%s   -   %.0f < E < %.0f GeV%s",
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

        // tighten X-range around content with a 20% buffer
        tighten_x_range(h, 0.20);

        // compute Y max within the visible X-range for tighter headroom (per pad)
        {
            int first = h->GetXaxis()->GetFirst();
            int last  = h->GetXaxis()->GetLast();
            double m  = 0.0;
            for (int b = first; b <= last; ++b) m = std::max(m, h->GetBinContent(b));
            h->GetYaxis()->SetRangeUser(0.0, 1.25 * std::max(1.0, m));
        }

        h->Draw("E1");

        // tighten X-range around content with a 20% buffer
        tighten_x_range(h, 0.20);

        // compute Y max within the currently visible X-range for tighter headroom
        double ymax = 1.0;
        {
            int first = h->GetXaxis()->GetFirst();
            int last  = h->GetXaxis()->GetLast();
            double m  = 0.0;
            for (int b = first; b <= last; ++b) m = std::max(m, h->GetBinContent(b));
            ymax = 1.25 * std::max(1.0, m);
        }
        h->GetYaxis()->SetRangeUser(0.0, ymax);
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
    //  EmitMinRmsSummary (QUIET):
    //   • Writes CSVs for min-RMS across bins (ALL / CLUS-only)
    //   • Also writes Δ-vs-CLUSRAW CSVs for φ and η
    //   • No console output (the final four ANSI tables are printed elsewhere)
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

      // Variant key lists
      std::vector<std::string> keysAll;
      for (const auto& v : kVariants) keysAll.push_back(v.key);
      const std::vector<std::string> keysClus = {
        "CLUSRAW","CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
      };

      // Sorted min-RMS rankings (for CSVs only)
      auto phiAllSorted  = collectSortedPerBin(keysAll , S_PH);
      auto etaAllSorted  = collectSortedPerBin(keysAll , S_ET);
      auto phiCluSorted  = collectSortedPerBin(keysClus, S_PH);
      auto etaCluSorted  = collectSortedPerBin(keysClus, S_ET);

      // CSVs for min-RMS
      writeCSVwithSecond(outDir + "/MinRMS_ALL_variants.csv",  phiAllSorted, etaAllSorted);
      writeCSVwithSecond(outDir + "/MinRMS_CLUS_only.csv",     phiCluSorted, etaCluSorted);

      // Δ-vs-CLUSRAW CSVs (quiet)
      auto writeDeltaCSV = [&](const std::string& csvPath,
                               const std::map<std::string,ResStats>& S,
                               const std::vector<std::string>& keys)
      {
        auto itB = S.find("CLUSRAW");
        if (itB==S.end()) return;
        const ResStats& B = itB->second;

        std::ofstream ofs(csvPath);
        ofs << "E_lo,E_hi";
        for (const auto& k : keys) ofs << "," << k << "_delta_pct," << k << "_delta_pct_err";
        ofs << "\n";

        for (int i=0;i<NB;++i){
          const double r  = (i<(int)B.rms.size()?  B.rms[i]  : NAN);
          const double re = (i<(int)B.drms.size()? B.drms[i] : 0.0);
          ofs << std::fixed << std::setprecision(0) << eSlices[i].first << "," << eSlices[i].second;

          for (const auto& k : keys){
            double s=NAN, se=0.0;
            auto it = S.find(k);
            if (it != S.end()){
              s  = (i<(int)it->second.rms.size()?  it->second.rms[i]  : NAN);
              se = (i<(int)it->second.drms.size()? it->second.drms[i] : 0.0);
            }
            if (std::isfinite(r) && r>0.0 && std::isfinite(s)){
              const double ratio   = s/r;
              const double varRat  = (se*se)/(r*r) + (s*s*re*re)/(r*r*r*r);
              const double dPct    = 100.0*(1.0 - ratio);
              const double dPctErr = 100.0*std::sqrt(std::max(0.0, varRat));
              ofs << "," << std::setprecision(3) << dPct << "," << dPctErr;
            } else {
              ofs << ",,";
            }
          }
          ofs << "\n";
        }
        ofs.close();
      };

      const std::vector<std::string> keysClusNoRaw = {
        "CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
      };
      writeDeltaCSV(outDir + "/Delta_vs_CLUSRAW_phi.csv", S_PH, keysClusNoRaw);
      writeDeltaCSV(outDir + "/Delta_vs_CLUSRAW_eta.csv", S_ET, keysClusNoRaw);
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
      // overlays (only |v_z| < 10 cm set)
      ensure_dir(base + "/Overlays/firstBin");
      ensure_dir(base + "/Overlays/meansigma");
      ensure_dir(base + "/Overlays/sigmaRatio");
      // NEW: percent-gain overlays vs CLUSCP baseline
      ensure_dir(base + "/Overlays/percentGain");

      // φ tilt-fit products
      if (base == DIR_PH) {
        ensure_dir(base + "/Overlays/phiTILTfits");
      }
    };
    
  scaffold_kind(DIR_PH);
  scaffold_kind(DIR_ET);

    using Var2Hists = std::map<std::string, HVec>;
    Var2Hists PH, ET;                            // variant -> [NB] hists
    for (const auto& v : kVariants) {
      PH[v.key] = HVec(NB,nullptr);
      ET[v.key] = HVec(NB,nullptr);
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
      const bool isVz60 = false; // no 0_60 variants in this production

      // Determine kind & variant by prefix match
      auto try_attach = [&](const KindDef& kind, const VarDef& vdef)->bool {
        const std::string& pref = (kind.tag=="phi" ? vdef.phiPrefix : vdef.etaPrefix);
        if (!TString(name.c_str()).BeginsWith(pref.c_str())) return false;

        std::string src;
        const int ib = resolve_energy_bin(h, name, eSlices, eCenters, &src);
        const auto& targetES = eSlices;
        const std::string ewin = (ib>=0 && ib<(int)targetES.size())
                                 ? (Form("%.0f-%.0f", targetES[ib].first, targetES[ib].second))
                                 : "";

        ScanRow row;
        row.kind    = kind.tag;
        row.variant = vdef.key;
        row.suffix  = ""; // no 0_60 suffix
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
            HVec& vec = (kind.tag=="phi" ? PH[vdef.key] : ET[vdef.key]);
            if (vec[ib]!=nullptr) {
              std::cerr << "[MakeResidualsSuite][WARN] Duplicate for ("<<kind.tag<<","<<vdef.key<<","<<ib<<"). Keeping the first and ignoring: " << name << "\n";
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


      // Longest-prefix match so generic prefixes don't swallow specific ones
      const VarDef* bestPhi = nullptr;  size_t bestPhiLen = 0;
      const VarDef* bestEta = nullptr;  size_t bestEtaLen = 0;

      for (const auto& vdef : kVariants) {
        const std::string& pphi = vdef.phiPrefix;
        if (!pphi.empty() && TString(name.c_str()).BeginsWith(pphi.c_str())) {
          if (pphi.size() > bestPhiLen) { bestPhi = &vdef; bestPhiLen = pphi.size(); }
        }
        const std::string& peta = vdef.etaPrefix;
        if (!peta.empty() && TString(name.c_str()).BeginsWith(peta.c_str())) {
          if (peta.size() > bestEtaLen) { bestEta = &vdef; bestEtaLen = peta.size(); }
        }
      }

      if (bestPhi) {
        try_attach(kPhi, *bestPhi);
      } else if (bestEta) {
        try_attach(kEta, *bestEta);
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

    Var2Stats S_PH, S_ET;
    compute_stats(PH, S_PH);
    compute_stats(ET, S_ET);


    for (const auto& v : kVariants) {
      if (PH.count(v.key)) {
        const std::string base = DIR_PH + "/" + v.key;
        WritePlaneResiduals(kPhi, v.key, eSlices, PH.at(v.key), base, /*isVz60=*/false);
      }
      if (ET.count(v.key)) {
        const std::string base = DIR_ET + "/" + v.key;
        WritePlaneResiduals(kEta, v.key, eSlices, ET.at(v.key), base, /*isVz60=*/false);
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

        // helper: provide human-friendly labels for EA flavours
        auto legend_label = [&](const std::string& key)->const char* {
          if (key == "CLUSCP_EA_vzDep")        return "Energy + |z_{vtx}|-dep b";
          if (key == "CLUSCP_EA_etaDep")       return "Energy + |#eta|-dep b";
          if (key == "CLUSCP_EA_vzEtaDep")     return "Energy + |z_{vtx}|+|#eta|-dep b";
          if (key == "CLUSCP_EA_EandIncident") return "Energy-only + incident-angle";
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
    
    
    // sigma_ratio_overlay with per-plot LEGEND POSITION and TEXT SIZE control
    // Usage: pass legend box (lx1,ly1,lx2,ly2) and text size (ltxt) for each output.
    std::function<void(const KindDef&,
                       const Var2Stats&,
                       const std::vector<std::string>&,
                       const std::string&,
                       const std::string&,
                       const char*,
                       double,double,double,double,double)>
    sigma_ratio_overlay = [&](const KindDef& kind,
                              const Var2Stats& S,
                              const std::vector<std::string>& which,
                              const std::string& baselineKey,
                              const std::string& outPng,
                              const char* vzNoteStr,
                              double lx1, double ly1, double lx2, double ly2, double ltxt)
    {
      const auto itB = S.find(baselineKey);
      if (itB==S.end()) return;
      const auto& B = itB->second;

      TCanvas c("cSR","",1000,850);
      TPad top ("top" ,"",0,0.35,1,1);  top.SetBottomMargin(0.02); top.SetLeftMargin(0.15); top.SetRightMargin(0.06); top.Draw();
      TPad bot ("bot" ,"",0,0.00,1,0.35); bot.SetTopMargin(0.04);  bot.SetLeftMargin(0.15); bot.SetRightMargin(0.06); bot.SetBottomMargin(0.18); bot.Draw();

      const double xMin = eSlices.front().first - 0.5;
      const double xMax = eSlices.back().second + 0.5;

      auto legend_label = [&](const std::string& key)->const char* {
        if (key == "CLUSCP_EA_vzDep")        return "Energy + |z_{vtx}|-dep b";
        if (key == "CLUSCP_EA_etaDep")       return "Energy + |#eta|-dep b";
        if (key == "CLUSCP_EA_vzEtaDep")     return "Energy + |z_{vtx}|+|#eta|-dep b";
        if (key == "CLUSCP_EA_EandIncident") return "Energy-only + incident-angle";
        return pretty_of(key);
      };

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

      // --- In-bin marker offsets (consistent across TOP/BOTTOM)
      const double kBinOffset = 0.12; // tune if you want more/less horizontal separation
      auto slotOf = [&](const std::string& key)->int {
        if (key == "CLUSCP")                 return -3; // Constant-b (legacy)
        if (key == "CLUSCP_EA_etaDep")       return -1; // |η|-dep
        if (key == "CLUSCP_EA_vzDep")        return +1; // |z|-dep
        if (key == "CLUSCP_EA_vzEtaDep")     return +3; // |z|+|η|-dep
        if (key == "CLUSCP_EA_EandIncident") return  0; // E-only + incident-angle (center)
        if (key == "CLUSCP_EA")              return +2; // E-only
        if (key == "CLUSRAW")                return -4; // raw (put slightly left of CP)
        return 0;
      };
      auto xOffsetFrac = [&](const std::string& key)->double {
        return slotOf(key) * kBinOffset;
      };

      // Draw top graphs
      std::vector<std::unique_ptr<TGraphErrors>> keepTop;
      std::map<std::string,TGraphErrors*> gByKey;
      for (const auto& k : which){
        auto it = S.find(k);
        if (it==S.end()) continue;

        std::vector<double> X,Y,EX,EY;
        for (int i=0;i<NB;++i){
          const double y  = it->second.rms[i];
          const double dy = it->second.drms[i];
          if (!std::isfinite(y)) continue;
          const double eLo  = eSlices[i].first;
          const double eHi  = eSlices[i].second;
          const double eCtr = 0.5*(eLo + eHi);
          const double eWid = (eHi - eLo);
          const double xOff = xOffsetFrac(k) * eWid;
          X.push_back(eCtr + xOff);
          Y.push_back(y);
          EX.push_back(0.0);
          EY.push_back(std::isfinite(dy)?dy:0.0);
        }

        auto g = std::make_unique<TGraphErrors>((int)X.size(),
                                                X.empty()?nullptr:&X[0],
                                                Y.empty()?nullptr:&Y[0],
                                                EX.empty()?nullptr:&EX[0],
                                                EY.empty()?nullptr:&EY[0]);
        if (g->GetN()>0){
          g->SetMarkerColor(color_of(k));
          g->SetLineColor  (color_of(k));
          g->SetMarkerStyle(kMkStyle);
          g->SetMarkerSize(1.10);
          g->SetLineWidth(1);
          g->Draw("P SAME");
          gByKey[k] = g.get();
          keepTop.emplace_back(std::move(g));
        }
      }

      // --- LEGEND (PER-PLOT TUNABLE) ---
      // Tune here per plot: lx1, ly1, lx2, ly2, ltxt (passed by caller)
      TLegend leg(lx1, ly1, lx2, ly2, "", "brNDC");
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(ltxt);
      leg.SetNColumns(1);
      leg.SetEntrySeparation(0.0);
      leg.SetMargin(0.18);
      for (const auto& k : which){
        auto it = gByKey.find(k);
        if (it!=gByKey.end() && it->second)
          leg.AddEntry(it->second, legend_label(k), "p");
      }
      leg.Draw();

      // ---------------- bottom: RMS / baseline ----------------
      bot.cd();

      std::vector<std::unique_ptr<TGraphErrors>> keepBot;
      std::vector<double> yminCand, ymaxCand;

      for (const auto& k : which){
        if (k==baselineKey) continue;
        auto it = S.find(k); if (it==S.end()) continue;

        std::vector<double> X,Y,EX,EY;
        for (int i=0;i<NB;++i){
          const double s  = it->second.rms[i];
          const double se = it->second.drms[i];
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
            const double xOff = xOffsetFrac(k) * eWid;

            X.push_back(eCtr + xOff);
            Y.push_back(ratio);
            EX.push_back(0.0);
            const double dy = ratio * std::sqrt(erel2);
            EY.push_back(dy);
            yminCand.push_back(ratio - dy);
            ymaxCand.push_back(ratio + dy);
          }
        }

        auto g = std::make_unique<TGraphErrors>((int)X.size(),
                                                X.empty()?nullptr:&X[0],
                                                Y.empty()?nullptr:&Y[0],
                                                EX.empty()?nullptr:&EX[0],
                                                EY.empty()?nullptr:&EY[0]);
        if (g->GetN()>0){
          g->SetMarkerColor(color_of(k));
          g->SetLineColor  (color_of(k));
          g->SetMarkerStyle(kMkStyle);
          g->SetMarkerSize(1.10);
          g->SetLineWidth(1);
          keepBot.emplace_back(std::move(g));
        }
      }

      double yLo = 0.98, yHi = 1.02;
      if (!yminCand.empty() && !ymaxCand.empty()){
        yLo = *std::min_element(yminCand.begin(), yminCand.end());
        yHi = *std::max_element(ymaxCand.begin(), ymaxCand.end());
        const double pad = 0.07*(yHi - yLo + 1e-12);
        yLo -= pad; yHi += pad;
      }

      TH1F frL("frL","",1,xMin,xMax);
      frL.SetTitle(Form(";E  [GeV];#sigma_{RMS}(%s) / #sigma_{RMS}(%s)",
                        (kind.tag=="phi"?"#Delta#phi":"#Delta#eta"),
                        pretty_of(baselineKey)));
      frL.SetMinimum(yLo);
      frL.SetMaximum(yHi);
      frL.GetYaxis()->SetLabelSize(0.055);
      frL.GetYaxis()->SetTitleSize(0.065);
      frL.GetYaxis()->SetTitleOffset(0.8);
      frL.Draw("AXIS");

      TLine l1(xMin,1.0,xMax,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

      for (auto& g : keepBot) g->Draw("P SAME");

      c.SaveAs(outPng.c_str());
    };

    
    
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

      // ------------------------------------------------------------------
      // NEW: two-panel plot in meansigma:
      //   Top: mean vs E for RAW, CP, and all EA variants
      //   Bottom: RMS ratio to CLUSRAW for CP and all EA variants
      //   (equal pad heights; proper error propagation; gentle x-offsets)
      // ------------------------------------------------------------------
      auto mean_plus_rmsratio_to_raw = [&](const std::vector<std::string>& keysMean,
                                           const std::vector<std::string>& keysRatio,
                                           const std::string& outPng)
      {
        // Baseline (RAW) for ratios
        auto itB = S.find("CLUSRAW");
        if (itB == S.end()) {
          std::cerr << "[mean_plus_rmsratio_to_raw] CLUSRAW baseline missing for " << kind.tag << "\n";
          return;
        }
        const ResStats& B = itB->second;

        TCanvas c("cMeanR","",1000,850);
        TPad top ("top" ,"",0,0.50,1,1);  top.SetBottomMargin(0.06); top.SetLeftMargin(0.15); top.SetRightMargin(0.06); top.Draw();
        TPad bot ("bot" ,"",0,0.00,1,0.50); bot.SetTopMargin(0.06);  bot.SetLeftMargin(0.15); bot.SetRightMargin(0.06); bot.SetBottomMargin(0.18); bot.Draw();

        const double xMin = eSlices.front().first - 0.5;
        const double xMax = eSlices.back().second + 0.5;

        // ---------- TOP: mean vs E ----------
        top.cd();
        double muLo = +1e30, muHi = -1e30;
        for (const auto& k : keysMean){
          auto it = S.find(k); if (it==S.end()) continue;
          for (size_t i=0;i<it->second.mean.size();++i){
            const double m  = it->second.mean[i];
            const double dm = (i<it->second.dmean.size()? it->second.dmean[i] : 0.0);
            if (!std::isfinite(m)) continue;
            muLo = std::min(muLo, m - (std::isfinite(dm)?dm:0.0));
            muHi = std::max(muHi, m + (std::isfinite(dm)?dm:0.0));
          }
        }
        if (!(muLo < muHi)) { muLo = -1e-4; muHi = +1e-4; }
        const double padMu = 0.15*(muHi - muLo + 1e-12);

        TH1F frU("frU","",1,xMin,xMax);
        frU.SetTitle(Form("Mean vs E (%s)%s;E  [GeV];#mu(%s) [rad]",
                          kind.tag.c_str(), vzn, (kind.tag=="phi"?"#Delta#phi":"#Delta#eta")));
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.GetYaxis()->SetTitleSize(0.058);
        frU.GetYaxis()->SetTitleOffset(0.90);
        frU.GetYaxis()->SetLabelSize(0.042);
        frU.Draw("AXIS");

        TLine l0(xMin,0.0,xMax,0.0); l0.SetLineStyle(2); l0.SetLineColor(kGray+2); l0.Draw();

          std::vector<std::unique_ptr<TGraphErrors>> keepU;
          TLegend legU(0.58, 0.72, 0.92, 0.90, "", "brNDC");
          legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.040);

          // Symmetric in-bin x-offsets around the bin center (7 series -> slots -3..0..+3)
          const int nMean = static_cast<int>(keysMean.size());
          const int halfMean = nMean / 2;                  // 3 for 7
          const bool oddMean = (nMean % 2 == 1);
          const double kOffsetFracMean = 0.12;             // fraction of bin width per slot
          auto slotByIndexMean = [&](int j)->int {
            return oddMean ? (j - halfMean)                // -half..0..+half (centered)
                           : ((j < halfMean) ? (j - halfMean) : (j - halfMean + 1));
          };

          for (int j=0; j<nMean; ++j){
            const std::string& k = keysMean[j];
            auto it = S.find(k); if (it==S.end()) continue;

            std::vector<double> X,Y,EX,EY;
            const int slot = slotByIndexMean(j);

            for (int i=0;i<NB;++i){
              const double m  = it->second.mean[i];
              const double dm = (i<(int)it->second.dmean.size()? it->second.dmean[i] : 0.0);
              if (!std::isfinite(m)) continue;

              const double eLo  = eSlices[i].first;
              const double eHi  = eSlices[i].second;
              const double eCtr = 0.5*(eLo + eHi);
              const double eWid = (eHi - eLo);
              const double xOff = slot * kOffsetFracMean * eWid;

              X.push_back(eCtr + xOff);
              Y.push_back(m);
              EX.push_back(0.0);
              EY.push_back(std::isfinite(dm) ? dm : 0.0);
            }

            auto g = std::make_unique<TGraphErrors>((int)X.size(),
                                                    X.empty()?nullptr:&X[0],
                                                    Y.empty()?nullptr:&Y[0],
                                                    EX.empty()?nullptr:&EX[0],
                                                    EY.empty()?nullptr:&EY[0]);
            if (g->GetN()>0){
              g->SetMarkerColor(color_of(k));
              g->SetLineColor  (color_of(k));
              g->SetLineWidth(1);
              g->SetMarkerSize(0.95);
              g->Draw("P SAME");
              legU.AddEntry(g.get(), pretty_of(k), "p");
              keepU.emplace_back(std::move(g));
            }
          }
          legU.Draw();

        // ---------- BOTTOM: RMS ratio to CLUSRAW ----------
        bot.cd();

          // Symmetric in-bin x-offsets with an empty center (even #series => no 0 slot)
          const int nRatio = static_cast<int>(keysRatio.size());   // expected 6
          const int halfR  = nRatio / 2;                           // 3
          const double kOffsetFracRatio = 0.12;
          auto slotByIndexRatio = [&](int j)->int {
            // Slots: [-halfR, ..., -1, +1, ..., +halfR]  (no 0)
            return (j < halfR) ? -(halfR - j) : (j - halfR + 1);
          };

          std::vector<std::unique_ptr<TGraphErrors>> keepL;
          double yLo = +1e30, yHi = -1e30;

          for (int j=0; j<nRatio; ++j){
            const std::string& k = keysRatio[j];
            auto it = S.find(k); if (it==S.end()) continue;

            std::vector<double> X,Y,EX,EY;
            const int slot = slotByIndexRatio(j);

            for (int i=0;i<NB;++i){
              const double s  = it->second.rms[i];
              const double se = (i<(int)it->second.drms.size()? it->second.drms[i] : 0.0);
              const double r  = B.rms[i];
              const double re = (i<(int)B.drms.size()? B.drms[i] : 0.0);

              if (std::isfinite(s) && std::isfinite(r) && r>0.0){
                const double ratio = s/r;
                const double erel2 = ((se>0&&s>0? (se/s)*(se/s) : 0.0) +
                                      (re>0&&r>0? (re/r)*(re/r) : 0.0));
                const double dy = ratio * std::sqrt(erel2);

                const double eLo  = eSlices[i].first;
                const double eHi  = eSlices[i].second;
                const double eCtr = 0.5*(eLo + eHi);
                const double eWid = (eHi - eLo);
                const double xOff = slot * kOffsetFracRatio * eWid;

                X.push_back(eCtr + xOff);
                Y.push_back(ratio);
                EX.push_back(0.0);
                EY.push_back(dy);

                yLo = std::min(yLo, ratio - dy);
                yHi = std::max(yHi, ratio + dy);
              }
            }

            auto g = std::make_unique<TGraphErrors>((int)X.size(),
                                                    X.empty()?nullptr:&X[0],
                                                    Y.empty()?nullptr:&Y[0],
                                                    EX.empty()?nullptr:&EX[0],
                                                    EY.empty()?nullptr:&EY[0]);
            if (g->GetN()>0){
              g->SetMarkerColor(color_of(k));
              g->SetLineColor  (color_of(k));
              g->SetMarkerStyle(kMkStyle);
              g->SetMarkerSize(0.95);
              g->SetLineWidth(1);
              keepL.emplace_back(std::move(g));
            }
          }

          if (!(yLo < yHi)) { yLo = 0.98; yHi = 1.02; }
          const double padS = 0.07*(yHi - yLo + 1e-12);

          TH1F frL("frL","",1,xMin,xMax);
          frL.SetTitle(Form(";E  [GeV];#sigma_{RMS}(%s) / #sigma_{RMS}(%s)",
                            (kind.tag=="phi"?"#Delta#phi":"#Delta#eta"), pretty_of("CLUSRAW")));
          frL.SetMinimum(yLo - padS);
          frL.SetMaximum(yHi + padS);
          frL.GetYaxis()->SetLabelSize(0.055);
          frL.GetYaxis()->SetTitleSize(0.065);
          frL.GetYaxis()->SetTitleOffset(0.80);
          frL.Draw("AXIS");

          TLine l1(xMin,1.0,xMax,1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();

          for (auto& g : keepL) g->Draw("P SAME");

          // NOTE: Legend intentionally omitted on the bottom subpanel.

        c.SaveAs(outPng.c_str());
      };

      // Keys for the new plot
      const std::vector<std::string> keysMean = {
        "CLUSRAW","CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
      };
      const std::vector<std::string> keysRatio = {
        "CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
      };

      // Emit the new two-panel plot (saved in the same meansigma folder)
      mean_plus_rmsratio_to_raw(keysMean, keysRatio,
                                baseOut + "/" + kind.tag + "__meanVsE__RMSratio_to_CLUSRAW.png");
    };
    
    

    // Saves the required sigmaRatio overlays with per-plot LEGEND tuning and
    // filenames built from the keys and the kind tag (phi/eta).
    auto save_sigma_groups = [&](const KindDef& kind,
                                 const Var2Stats& S,
                                 const std::string& baseOut,
                                 bool isVz60)
    {
      ensure_dir(baseOut);
      const char* vzn = vz_note(isVz60);

      // -------- helper: join keys for filename stem --------
      auto join_keys = [](const std::vector<std::string>& v)->std::string {
        std::string s;
        for (size_t i=0;i<v.size();++i){ s += v[i]; if (i+1<v.size()) s += "__"; }
        return s;
      };

      // -------- helper: emit one named overlay with legend tuning --------
      auto emit_named = [&](const std::vector<std::string>& keys,
                            const std::string& baselineKey,
                            double lx1, double ly1, double lx2, double ly2, double ltxt)
      {
        const std::string pngName = baseOut + "/" + kind.tag + "__" + join_keys(keys) + ".png";
        sigma_ratio_overlay(kind, S, keys, baselineKey, pngName, vzn,
                            lx1, ly1, lx2, ly2, ltxt);
      };

      // ====== NEW REQUESTED OUTPUTS (NAMES BY KEYS + KIND) ======

      // 1) Top: RMS overlay of { EandIncident, E-only, CP }, Bottom: ratios of {EandIncident, E-only} to CP
      //    Keys include baseline CLUSCP so the filename reflects all plotted series.
      {
        std::vector<std::string> keys = {"CLUSCP_EA_EandIncident","CLUSCP_EA","CLUSCP"};
        // LEGEND CONFIG for this PNG (tune here): x1,y1,x2,y2,textSize
        emit_named(keys, /*baseline=*/"CLUSCP", 0.62, 0.68, 0.92, 0.88, 0.044);
      }

      // 2) Top: RMS overlay of { E-only, |z|, |eta|, |z|+|eta|, EandIncident, CP }, Bottom: ratios of all (except CP) to CP
      {
        std::vector<std::string> keys = {
          "CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident","CLUSCP"
        };
        // LEGEND CONFIG for this PNG
        emit_named(keys, /*baseline=*/"CLUSCP", 0.60, 0.66, 0.92, 0.90, 0.040);
      }

      // 3) For each EA variant: Top overlays {RAW, CP, VARIANT}; Bottom: ratios of {CP, VARIANT} to RAW
      {
        const std::vector<std::string> variants = {
          "CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
        };
        for (const auto& vkey : variants){
          std::vector<std::string> keys = {"CLUSRAW","CLUSCP", vkey};
          // LEGEND CONFIG for each per-variant PNG
          emit_named(keys, /*baseline=*/"CLUSRAW", 0.58, 0.68, 0.92, 0.88, 0.042);
        }
      }

      // 4) Top: overlays {CP, all EA variants, RAW}; Bottom: ratios of {CP, all EA variants} to RAW
      {
        std::vector<std::string> keys = {
          "CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident","CLUSRAW"
        };
        // LEGEND CONFIG for this PNG
        emit_named(keys, /*baseline=*/"CLUSRAW", 0.55, 0.66, 0.92, 0.90, 0.038);
      }

      // ====== KEEP ORIGINAL (LEGACY-NAMED) OUTPUTS TO PRESERVE EXISTING PRODUCTS ======
      // You can tune their legend positions/sizes here as well.
      {
        sigma_ratio_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA"}, "CLUSRAW",
                            baseOut + "/sigma_vsE_ratio_to_CLUSRAW_CP_EA.png", vzn,
                            0.62, 0.68, 0.92, 0.88, 0.044);

        sigma_ratio_overlay(kind, S, {"CLUSRAW","CLUSCP","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep"}, "CLUSRAW",
                            baseOut + "/sigma_vsE_ratio_to_CLUSRAW_CP_EA_vz_eta.png", vzn,
                            0.60, 0.66, 0.92, 0.90, 0.040);

        sigma_ratio_overlay(kind, S,
                            {"CLUSCP","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident","CLUSCP_EA"},
                            "CLUSCP",
                            baseOut + "/sigma_vsE_ratio_to_CLUSCP_EA_vz_eta.png",
                            vzn,
                            0.55, 0.66, 0.92, 0.90, 0.038);

        sigma_ratio_overlay(kind, S, {"CLUSRAW","PDCcor"}, "CLUSRAW",
                            baseOut + "/sigma_vsE_ratio_PDCcor_to_CLUSRAW.png", vzn,
                            0.65, 0.70, 0.92, 0.88, 0.042);
      }
    };

    // Emit overlays (only |v_z| < 10 cm set)
    save_meansigma_groups(kPhi, S_PH, DIR_PH + "/Overlays/meansigma", false);
    save_meansigma_groups(kEta, S_ET, DIR_ET + "/Overlays/meansigma", false);

    save_sigma_groups(kPhi, S_PH, DIR_PH + "/Overlays/sigmaRatio", false);
    save_sigma_groups(kEta, S_ET, DIR_ET + "/Overlays/sigmaRatio", false);

    // ------------------------------------------------------------------
    // NEW: Percent-gain overlays wrt CLUSCP (RMS), for φ and η
    //   y = 100 * (1 - σ_variant / σ_CLUSCP)
    //   dy = 100 * sqrt( (σe_variant^2 / σ_CLUSCP^2) + (σ_variant^2 * σe_CLUSCP^2 / σ_CLUSCP^4) )
    //   (same horizontal in-bin offsetting scheme as sigma_ratio_overlay)
    // ------------------------------------------------------------------
    auto percent_gain_overlay = [&](const KindDef& kind,
                                    const Var2Stats& S,
                                    const std::vector<std::string>& variants,
                                    const std::string& baselineKey,
                                    const std::string& outPng,
                                    bool isVz60)
    {
      auto itB = S.find(baselineKey);
      if (itB == S.end()) {
        std::cerr << "[percentGain] Baseline " << baselineKey << " not found for " << kind.tag << "\n";
        return;
      }
      const ResStats& B = itB->second;

      TCanvas c("cPG","",1000,700);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.06);
      gPad->SetBottomMargin(0.14);

      const double xMin = eSlices.front().first - 0.5;
      const double xMax = eSlices.back().second + 0.5;

      // Legend labels consistent with other overlays
      auto legend_label = [&](const std::string& key)->const char* {
        if (key == "CLUSCP_EA_vzDep")        return "Energy + |z_{vtx}|-dep b";
        if (key == "CLUSCP_EA_etaDep")       return "Energy + |#eta|-dep b";
        if (key == "CLUSCP_EA_vzEtaDep")     return "Energy + |z_{vtx}|+|#eta|-dep b";
        if (key == "CLUSCP_EA_EandIncident") return "Energy-only + incident-angle";
        if (key == "CLUSCP_EA")              return "Energy-Dep b";
        return pretty_of(key);
      };

      // In-bin offset slots (reuse scheme for visual separation)
      const double kBinOffset = 0.12;
      auto slotOf = [&](const std::string& key)->int {
        if (key == "CLUSCP_EA_etaDep")       return -1; // |η|-dep
        if (key == "CLUSCP_EA_vzDep")        return +1; // |z|-dep
        if (key == "CLUSCP_EA_vzEtaDep")     return +3; // |z|+|η|-dep
        if (key == "CLUSCP_EA_EandIncident") return  0; // E-only + incident-angle
        if (key == "CLUSCP_EA")              return +2; // E-only
        return 0;
      };

      double yMin = +1e30, yMax = -1e30;
      std::vector<std::unique_ptr<TGraphErrors>> keep;
      std::map<std::string,TGraphErrors*> gByKey;

      for (const auto& key : variants) {
        auto it = S.find(key);
        if (it == S.end()) continue;

        std::vector<double> X, Y, EX, EY;
        X.reserve(NB); Y.reserve(NB); EX.reserve(NB); EY.reserve(NB);

        for (int i=0;i<NB;++i){
          const double s  = (i<(int)it->second.rms.size()?  it->second.rms[i]  : std::numeric_limits<double>::quiet_NaN());
          const double se = (i<(int)it->second.drms.size()? it->second.drms[i] : 0.0);
          const double r  = (i<(int)B.rms.size()?           B.rms[i]           : std::numeric_limits<double>::quiet_NaN());
          const double re = (i<(int)B.drms.size()?          B.drms[i]          : 0.0);
          if (!(std::isfinite(s) && std::isfinite(r) && r>0.0)) continue;

          const double imp    = 100.0 * (1.0 - s/r);
          const double varImp = (se*se)/(r*r) + (s*s*re*re)/(r*r*r*r);
          const double err    = 100.0 * std::sqrt(std::max(0.0, varImp));

          const double eLo  = eSlices[i].first;
          const double eHi  = eSlices[i].second;
          const double eCtr = 0.5*(eLo + eHi);
          const double eWid = (eHi - eLo);
          const double xOff = slotOf(key) * kBinOffset * eWid;

          X.push_back(eCtr + xOff);
          Y.push_back(imp);
          EX.push_back(0.0);
          EY.push_back(err);

          yMin = std::min(yMin, imp - err);
          yMax = std::max(yMax, imp + err);
        }

        auto g = std::make_unique<TGraphErrors>((int)X.size(),
                                                X.empty()?nullptr:&X[0],
                                                Y.empty()?nullptr:&Y[0],
                                                EX.empty()?nullptr:&EX[0],
                                                EY.empty()?nullptr:&EY[0]);
        if (g->GetN() > 0){
          g->SetMarkerColor(color_of(key));
          g->SetLineColor  (color_of(key));
          g->SetMarkerStyle(kMkStyle);
          g->SetMarkerSize(0.8);
          g->SetLineWidth(1);
          gByKey[key] = g.get();
          keep.emplace_back(std::move(g));
        }
      }

      if (!(yMin < yMax)) { yMin = -5.0; yMax = +5.0; }
      const double padY = 0.12*(yMax - yMin + 1e-12);

      TH1F fr("fr","",1,xMin,xMax);
      fr.SetTitle(Form("Percent RMS improvement vs E (%s)%s;E  [GeV];#Delta #sigma_{RMS}(%s) wrt Legacy CP  [%%]",
                        kind.tag.c_str(), vz_note(isVz60),
                        (kind.tag == "phi" ? "#phi" : "#eta")));
        
      fr.SetMinimum(yMin - padY);
      fr.SetMaximum(yMax + padY);
      fr.Draw("AXIS");

      TLine l0(xMin, 0.0, xMax, 0.0); l0.SetLineStyle(2); l0.SetLineColor(kGray+2); l0.Draw();

      for (auto& g : keep) g->Draw("P SAME");

      TLegend leg(0.60, 0.66, 0.92, 0.90, "", "brNDC");
      leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
      leg.SetNColumns(1); leg.SetEntrySeparation(0.0); leg.SetMargin(0.18);
      for (const auto& key : variants){
        auto itg = gByKey.find(key);
        if (itg != gByKey.end() && itg->second)
          leg.AddEntry(itg->second, legend_label(key), "p");
      }
      leg.Draw();

      c.SaveAs(outPng.c_str());
    };

    // Requested EA variants to compare against CLUSCP
    const std::vector<std::string> kEAkeys = {
      "CLUSCP_EA", "CLUSCP_EA_vzDep", "CLUSCP_EA_etaDep", "CLUSCP_EA_vzEtaDep", "CLUSCP_EA_EandIncident"
    };

    // One PNG each for φ and η; points at each energy bin with proper error propagation
    percent_gain_overlay(kPhi, S_PH, kEAkeys, "CLUSCP",
                         DIR_PH + "/Overlays/percentGain/percentGain_vs_CLUSCP_phi.png",
                         /*isVz60=*/false);

    percent_gain_overlay(kEta, S_ET, kEAkeys, "CLUSCP",
                         DIR_ET + "/Overlays/percentGain/percentGain_vs_CLUSCP_eta.png",
                         /*isVz60=*/false);

    
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

  // η
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP",                 eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA",              eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_vzDep",        eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "CLUSRAW","CLUSCP_EA_etaDep",       eSlices[0], false);
  first_bin_pair(kEta, ET,    DIR_ET + "/Overlays/firstBin", "PDCraw","PDCcor",                  eSlices[0], false);

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

      TLegend legU(0.14,0.73,0.94,0.92); legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.02);
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
          {"CLUSRAW",                 "ETiltVariant::CLUS_RAW",                     "no corr, cluster"},
          {"CLUSCP",                  "ETiltVariant::CLUS_CP",                      "CorrectPosition, cluster"},
          {"CLUSCP_EA_vzEtaDep",      "ETiltVariant::CLUS_CP_EA_FIT_ZDEP_ETADEP",   "CorrectPosition(EA |z|+|#eta|+E), cluster"},
          {"CLUSCP_EA_vzDep",         "ETiltVariant::CLUS_CP_EA_FIT_ZDEP",          "CorrectPosition(EA |z|+E), cluster"},
          {"CLUSCP_EA_etaDep",        "ETiltVariant::CLUS_CP_EA_FIT_ETADEP",        "CorrectPosition(EA |#eta|+E), cluster"},
          {"CLUSCP_EA_EandIncident",  "ETiltVariant::CLUS_CP_EA_FIT_EandIncident",  "CorrectPosition(EA E-only + incident-angle), cluster"},
          {"CLUSCP_EA",               "ETiltVariant::CLUS_CP_EA_FIT_EONLY",         "CorrectPosition(EA E-only), cluster"},
          {"CLUScpDirectBvals",       "ETiltVariant::CLUS_BCORR",                   "b(E) corr, cluster"},
          {"PDCraw",                  "ETiltVariant::PDC_RAW",                      "no corr, scratch"},
          {"PDCcor",                  "ETiltVariant::PDC_CORR",                     "b(E) corr, scratch"}
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

    // ---- Min-RMS CSVs (quiet) — φ only in phiTILTfits ----
    {
      const std::string summaryDir = DIR_PH + std::string("/Overlays/phiTILTfits");
      Var2Stats S_ET_dummy; // do not write any eta artifacts under phiTILTfits
      EmitMinRmsSummary(summaryDir, NB, eSlices, S_PH, S_ET_dummy);
    }

    auto printFinalDeltaTable = [&](const char* titleTag,
                                    const std::map<std::string, ResStats>& S,
                                    const std::vector<std::string>& keys)
    {
      auto itB = S.find("CLUSRAW");
      if (itB==S.end()) { std::cout << "[WARN] CLUSRAW baseline missing for " << titleTag << "\n"; return; }
      const ResStats& B = itB->second;

      const int Wbin = 9;    // E-bin field width
      const int Wcol = 24;   // per-variant column width (pad BEFORE coloring)
      const int Wwin = 18;   // Winner column width
        const char* BOLD  = "\033[1m";
        const char* GREEN = "\033[1;32m";
        const char* RED   = "\033[1;31m";
        const char* BLUE  = "\033[1;34m";
        const char* RESET = "\033[0m";

        // Locate Constant-b, EA(E-only), and EA(E)+θ_inc columns
        int idxCP = -1, idxEAonly = -1, idxEAinc = -1;
        for (size_t j=0; j<keys.size(); ++j) {
          if (keys[j] == "CLUSCP")                 idxCP     = (int)j;
          if (keys[j] == "CLUSCP_EA")              idxEAonly = (int)j;
          if (keys[j] == "CLUSCP_EA_EandIncident") idxEAinc  = (int)j;
        }


      // Compact, consistent header labels (keeps columns narrow & aligned)
      auto headerLabel = [&](const std::string& k)->std::string {
        if (k=="CLUSCP")                 return "CP (const-b)";
        if (k=="CLUSCP_EA")              return "EA (E-only)";
        if (k=="CLUSCP_EA_vzDep")        return "EA (|z|)";
        if (k=="CLUSCP_EA_etaDep")       return "EA (|#eta|)";
        if (k=="CLUSCP_EA_vzEtaDep")     return "EA (|z|+|#eta|)";
        if (k=="CLUSCP_EA_EandIncident") return "EA (E)+θ_inc";
        if (k=="PDCcor")                 return "PDC corr";
        if (k=="PDCraw")                 return "PDC raw";
        if (k=="CLUSRAW")                return "No corr";
        if (k=="CLUScpDirectBvals")      return "b(E) apply";
        return pretty_of(k);
      };
      auto truncateTo = [&](const std::string& s, int w)->std::string {
        if ((int)s.size() <= w) return s;
        if (w <= 1) return s.substr(0, std::max(0,w));
        return s.substr(0, w-1) + "…";
      };

      // ===== Header (two rows: labels, then "Δ% ± err") =====
      std::cout << "\n" << BOLD << "[FINAL ΔRMS vs CLUSRAW]  " << titleTag << RESET << "\n";

      // Row 1: variant labels
      std::cout << std::left << std::setw(Wbin) << "E-bin";
      for (const auto& k : keys) {
        const std::string lab = truncateTo(headerLabel(k), Wcol);
        std::cout << " | " << std::setw(Wcol) << lab;
      }
      std::cout << " | " << std::setw(Wwin) << "Winner" << "\n";

      // Row 2: unit/format row directly above values
      std::cout << std::left << std::setw(Wbin) << "";
      for (size_t j=0; j<keys.size(); ++j) {
        std::cout << " | " << std::setw(Wcol) << "Δ% ± err";
      }
      std::cout << " | " << std::setw(Wwin) << "" << "\n";

      // Rule
      std::cout << std::string(Wbin + (int)keys.size()*(3+Wcol) + 3 + Wwin, '-') << "\n";

      // ===== Body =====
      for (int i=0; i<NB; ++i)
      {
        const double r  = (i<(int)B.rms.size()?  B.rms[i]  : std::numeric_limits<double>::quiet_NaN());
        const double re = (i<(int)B.drms.size()? B.drms[i] : 0.0);

        std::vector<double> deltaPct(keys.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> deltaErr(keys.size(), std::numeric_limits<double>::quiet_NaN());

        // Per-variant Δ% vs CLUSRAW
        for (size_t j=0; j<keys.size(); ++j){
          auto it = S.find(keys[j]);
          if (it==S.end()) continue;
          const ResStats& V = it->second;
          if (!(i<(int)V.rms.size())) continue;
          const double s  = V.rms[i];
          const double se = (i<(int)V.drms.size()? V.drms[i] : 0.0);
          if (std::isfinite(r) && r>0.0 && std::isfinite(s)){
            const double ratio   = s/r;
            const double varRat  = (se*se)/(r*r) + (s*s*re*re)/(r*r*r*r);
            deltaPct[j] = 100.0*(1.0 - ratio);
            deltaErr[j] = 100.0*std::sqrt(std::max(0.0, varRat));
          }
        }

        // Winner (max Δ%)
        int winIdx = -1; double winVal = -1e300;
        for (size_t j=0; j<keys.size(); ++j){
          if (std::isfinite(deltaPct[j]) && deltaPct[j] > winVal) { winVal = deltaPct[j]; winIdx = (int)j; }
        }

          // Constant-b reference for this row (if available)
          const double cpRef     = (idxCP>=0     ? deltaPct[idxCP]     : std::numeric_limits<double>::quiet_NaN());
          const double eaOnlyRef = (idxEAonly>=0 ? deltaPct[idxEAonly] : std::numeric_limits<double>::quiet_NaN());

          // E-bin label
          std::ostringstream ebin; ebin.setf(std::ios::fixed);

        ebin << std::setprecision(0) << eSlices[i].first << "-" << eSlices[i].second;
        std::cout << std::left << std::setw(Wbin) << ebin.str();

          // Values row (pad → then apply color)
          for (size_t j=0; j<keys.size(); ++j){
            std::string cell;
            if (std::isfinite(deltaPct[j])){
              std::ostringstream ss; ss.setf(std::ios::fixed);
              ss << std::showpos << std::setprecision(2) << deltaPct[j] << "%" << std::noshowpos
                 << " ± " << std::setprecision(2) << deltaErr[j] << "%";
              cell = ss.str();
            } else {
              cell = "n/a";
            }
            if ((int)cell.size() < Wcol) cell += std::string(Wcol - (int)cell.size(), ' ');

            // Red highlight: underperforms Constant-b (strictly smaller Δ%)
            const bool underCP = (idxCP>=0 && (int)j!=idxCP &&
                                  std::isfinite(cpRef) && std::isfinite(deltaPct[j]) &&
                                  deltaPct[j] < cpRef);

            // Blue highlight: EA(E)+θ_inc worse than EA(E-only) in the same bin
            const bool worseThanEAonly = (idxEAinc>=0 && (int)j==idxEAinc &&
                                          idxEAonly>=0 &&
                                          std::isfinite(deltaPct[j]) &&
                                          std::isfinite(eaOnlyRef) &&
                                          deltaPct[j] < eaOnlyRef);

            std::string colored = cell;
            if (worseThanEAonly) {
              colored = std::string(BLUE) + cell + RESET;
            } else if (underCP) {
              colored = std::string(RED) + cell + RESET;
            }

            std::cout << " | " << colored;
          }


        // Winner column (green)
        if (winIdx >= 0){
          std::ostringstream w; w << GREEN << headerLabel(keys[winIdx]) << RESET;
          std::string wt = w.str();
          if ((int)wt.size() < Wwin) wt += std::string(Wwin - (int)wt.size(), ' ');
          std::cout << " | " << wt << "\n";
        } else {
          std::cout << " | " << std::setw(Wwin) << "-" << "\n";
        }
      }
    };


    // 1) CLUS family vs CLUSRAW (φ then η)
    const std::vector<std::string> keysClusNoRaw = {
      "CLUSCP","CLUSCP_EA","CLUSCP_EA_vzDep","CLUSCP_EA_etaDep","CLUSCP_EA_vzEtaDep","CLUSCP_EA_EandIncident"
    };
    printFinalDeltaTable("phi (Δφ) — CLUS family", S_PH, keysClusNoRaw);
    printFinalDeltaTable("eta (Δη) — CLUS family", S_ET, keysClusNoRaw);

    // 2) ALL variants vs CLUSRAW (φ then η) — exclude baseline itself
    std::vector<std::string> keysAllNoRaw;
    for (const auto& v : kVariants) if (v.key != "CLUSRAW") keysAllNoRaw.push_back(v.key);
    printFinalDeltaTable("phi (Δφ) — ALL variants", S_PH, keysAllNoRaw);
    printFinalDeltaTable("eta (Δη) — ALL variants", S_ET, keysAllNoRaw);

    std::cout << "\n[MakeResidualsSuite] Residuals suite written under: " << RESID << "\n";
}


void PDCAnalysisPRIMEPRIME()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // --- paths
  const std::string inFilePath = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/z_lessThen_10/PositionDep_sim_ALL_noPhiTilt.root";
  const std::string outBaseDir = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/z_lessThen_10/SimOutputPrimeNoPhiTilt";

  EnsureDir(outBaseDir);

  // Open file
  std::unique_ptr<TFile> fin(TFile::Open(inFilePath.c_str(),"READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open input file: " << inFilePath << "\n";
    return;
  }

    // Optional: dump all histograms in the input file
    const bool dumpHistos = true;  // <-- set true to print the full histogram inventory
    if (dumpHistos)
    {
      std::cout << "\n[HISTO DUMP] Listing histograms in: " << inFilePath << "\n";
      Long64_t nHist = 0;

      // recursive directory walker
      std::function<void(TDirectory*, const std::string&)> walk;
      walk = [&](TDirectory* dir, const std::string& prefix)
      {
        if (!dir) return;
        TIter it(dir->GetListOfKeys());
        while (TObject* o = it())
        {
          auto* key = dynamic_cast<TKey*>(o);
          if (!key) continue;

          // Read object; we delete it after use to avoid leaks
          std::unique_ptr<TObject> obj(key->ReadObj());
          if (!obj) continue;

          const char* kname = key->GetName();
          const char* cname = obj->ClassName();

          if (obj->InheritsFrom(TH1::Class()))
          {
            auto* h = static_cast<TH1*>(obj.get());
            std::cout << "  - " << prefix << kname
                      << "  [" << cname << "]  entries=" << static_cast<long long>(h->GetEntries())
                      << "\n";
            ++nHist;
          }
          else if (obj->InheritsFrom(TDirectory::Class()))
          {
            auto* sub = static_cast<TDirectory*>(obj.get());
            walk(sub, prefix + std::string(kname) + "/");
          }
        }
      };

      walk(fin.get(), "/");
      std::cout << "[HISTO DUMP] total histograms found: " << nHist << "\n\n";
    }

    // Energy slices
    const auto eEdges = MakeEnergySlices();
    const auto eCent  = MakeEnergyCenters();

    SaveIncidenceQA(fin.get(), outBaseDir);

    
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
        const std::string outDir = outBaseDir + "/etaDepOnlyBlockAna/" + v.key;
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
//            const char* simPath = "/Users/patsfan753/Desktop/PositionDependentCorrection/SINGLE_PHOTON_MC/z_lessThen_10/PositionDep_sim_ALL.root";
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

          // ALSO save the baseline (no extra η or z binning) into a dedicated folder
          if (v.key == "originalEta")
          {
            const std::string outDirBase = outBaseDir + "/originalEtaAndOriginalZ";
            EnsureDir(outDirBase);

            // Uncorrected products
            SaveTH3Lego(hUnc, (outDirBase + "/lego_unc.png").c_str(), kEtaPretty.at(v.key).c_str());
            // Plot2DBlockEtaPhi writes BOTH BlockCoord2D_E_unc.png and (if !isFirstPass && hCor) BlockCoord2D_E_cor.png
            Plot2DBlockEtaPhi(hUnc, hCor, isFirstPass, eEdges, outDirBase.c_str(), kEtaPretty.at(v.key).c_str());
            OverlayUncorrPhiEta(hUnc, eEdges, outDirBase.c_str(), kEtaPretty.at(v.key).c_str());
            MakeBvaluesVsEnergyPlot(hUnc, eEdges, outDirBase.c_str(), kEtaPretty.at(v.key).c_str());

            // b-values text labeling: originalEta × originalZRange
            if (bOut.is_open()) WriteBValuesTxt(hUnc, eEdges, bOut, "originalEta", "originalZRange");

            // Corrected (second pass only): save corrected lego AND 1D corrected-vs-uncorrected overlays
            if (!isFirstPass && hCor)
            {
              SaveTH3Lego(hCor, (outDirBase + "/lego_cor.png").c_str(), kEtaPretty.at(v.key).c_str());
              // 1D overlays (corrected vs uncorrected) in the same folder:
              //   LocalPhiFits_4by2.png   and   LocalEtaFits_4by2.png
              FitLocalPhiEta(hUnc, hCor, /*isFirstPassFlag*/false, eEdges, outDirBase.c_str());
            }
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

        // NEW output path for z-only: zDepOnlybvalues/<zTag>
        const std::string outDirZ = outBaseDir + "/zDepOnlybvalues/" + ztag;
        EnsureDir(outDirZ);

        // Uncorrected products (always)
        SaveTH3Lego(hZ_unc, (outDirZ + "/lego_unc.png").c_str(), kZPretty.at(ztag).c_str());
        Plot2DBlockEtaPhi(hZ_unc, hZ_cor, /*firstPass*/isFirstPass, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());
        OverlayUncorrPhiEta(hZ_unc, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());
        BVecs zbv = MakeBvaluesVsEnergyPlot(hZ_unc, eEdges, outDirZ.c_str(), kZPretty.at(ztag).c_str());

        // Persist per-z b(E) for overlays
        zPhiByVar[ztag]     = zbv.bphi;  zPhiErrByVar[ztag] = zbv.bphiErr;
        zEtaByVar[ztag]     = zbv.beta;  zEtaErrByVar[ztag] = zbv.betaErr;

        // b-values text labeling: originalEta × (this z-slice)
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

        // NEW output path for fine z-only: zDepOnlybvalues/<zTag>
        const std::string outDirZf = outBaseDir + "/zDepOnlybvalues/" + ztag;
        EnsureDir(outDirZf);

        // Uncorrected fine products
        SaveTH3Lego(hZf_unc, (outDirZf + "/lego_unc.png").c_str(), kZFinePretty.at(ztag).c_str());
        Plot2DBlockEtaPhi(hZf_unc, hZf_cor, /*firstPass*/isFirstPass, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());
        OverlayUncorrPhiEta(hZf_unc, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());
        BVecs zbvFine = MakeBvaluesVsEnergyPlot(hZf_unc, eEdges, outDirZf.c_str(), kZFinePretty.at(ztag).c_str());

        zPhiFineByVar[ztag]     = zbvFine.bphi;  zPhiFineErrByVar[ztag] = zbvFine.bphiErr;
        zEtaFineByVar[ztag]     = zbvFine.beta;  zEtaFineErrByVar[ztag] = zbvFine.betaErr;

        // b-values text labeling: originalEta × (this fine z-slice)
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

        // ========== NEW: η × |z| cross-slice processing (with per-η z-overlays) ==========
        {
          // η keys to consider (pretty labels come from kEtaPretty)
          const std::vector<std::string> etaKeys = {"fullEta","etaCore","etaMid","etaEdge","originalEta"};

          // Accumulators to build per-η overlays across z
          std::map<std::string, std::map<std::string, std::vector<double>>> etaCoarsePhi, etaCoarsePhiErr;
          std::map<std::string, std::map<std::string, std::vector<double>>> etaCoarseEta, etaCoarseEtaErr;
          std::map<std::string, std::map<std::string, std::vector<double>>> etaFinePhi,   etaFinePhiErr;
          std::map<std::string, std::map<std::string, std::vector<double>>> etaFineEta,   etaFineEtaErr;

          // Helper: find the best-available uncorr/corr histogram by names
          auto getXZ = [&](const std::string& etaKey, const std::string& ztag, bool corrected)->TH3F*
          {
            std::vector<TString> names;
            if (etaKey == "originalEta")
            {
              // fall back to z-only names for originalEta
              if (!corrected) {
                names = { Form("h3_blockCoord_E_%s_range", ztag.c_str()),
                          Form("h3_blockCoord_E_%s_disc" , ztag.c_str()) };
              } else {
                names = { Form("h3_blockCoord_Ecorr_%s_range", ztag.c_str()),
                          Form("h3_blockCoord_Ecorr_%s_disc" , ztag.c_str()) };
              }
            }
            else
            {
              // true η×z cross-slice names
              if (!corrected) {
                names = { Form("h3_blockCoord_E_%s_%s_range", etaKey.c_str(), ztag.c_str()),
                          Form("h3_blockCoord_E_%s_%s_disc" , etaKey.c_str(), ztag.c_str()) };
              } else {
                names = { Form("h3_blockCoord_Ecorr_%s_%s_range", etaKey.c_str(), ztag.c_str()),
                          Form("h3_blockCoord_Ecorr_%s_%s_disc" , etaKey.c_str(), ztag.c_str()) };
              }
            }
            return GetTH3FByNames(fin.get(), names);
          };

          // -------- COARSE z inside each etaKey --------
          for (const auto& etaKey : etaKeys)
          {
            const std::string baseEtaDir = outBaseDir + "/etaAndzDepBlocks/" + etaKey;
            EnsureDir(baseEtaDir);

            for (const auto& ztag : kZOrder)
            {
              TH3F* h_unc = getXZ(etaKey, ztag, /*corrected*/false);
              if (!h_unc) continue;
              h_unc->SetDirectory(nullptr);

              TH3F* h_cor = (!isFirstPass ? getXZ(etaKey, ztag, /*corrected*/true) : nullptr);
              if (h_cor) h_cor->SetDirectory(nullptr);

              const std::string outDirXZ = baseEtaDir + "/" + ztag;
              EnsureDir(outDirXZ);

              const char* prettyEta = kEtaPretty.count(etaKey) ? kEtaPretty.at(etaKey).c_str()
                                                               : etaKey.c_str();
              const char* prettyZ   = kZPretty.count(ztag) ? kZPretty.at(ztag).c_str()
                                                           : ztag.c_str();

              // Uncorrected products for each (η, z_coarse)
              SaveTH3Lego(h_unc, (outDirXZ + "/lego_unc.png").c_str(), prettyEta);
              Plot2DBlockEtaPhi(h_unc, h_cor, isFirstPass, eEdges, outDirXZ.c_str(), prettyZ);
              OverlayUncorrPhiEta(h_unc, eEdges, outDirXZ.c_str(), prettyEta);
              BVecs bXZ = MakeBvaluesVsEnergyPlot(h_unc, eEdges, outDirXZ.c_str(), prettyEta);

              // Accumulate b(E) for per-η z overlays
              etaCoarsePhi[etaKey][ztag]    = bXZ.bphi;
              etaCoarsePhiErr[etaKey][ztag] = bXZ.bphiErr;
              etaCoarseEta[etaKey][ztag]    = bXZ.beta;
              etaCoarseEtaErr[etaKey][ztag] = bXZ.betaErr;

              // b-values text labeling: (etaKey) × (ztag)
              if (bOut.is_open()) WriteBValuesTxt(h_unc, eEdges, bOut, etaKey.c_str(), ztag.c_str());

              // Corrected (second pass)
              if (h_cor && !isFirstPass)
              {
                SaveTH3Lego(h_cor, (outDirXZ + "/lego_cor.png").c_str(), prettyEta);
                FitLocalPhiEta(h_unc, h_cor, /*isFirstPassFlag*/false, eEdges, outDirXZ.c_str());
              }

              delete h_unc;
              if (h_cor) delete h_cor;
            }

            // After filling all coarse z for this η, write the two coarse overlays (upto60 & high-only)
            if (!etaCoarsePhi[etaKey].empty())
            {
              SaveOverlayAcrossVariants(baseEtaDir, eCent, etaCoarsePhi[etaKey], etaCoarsePhiErr[etaKey], /*isPhi*/true ,  "_zOnly");
              SaveOverlayAcrossVariants(baseEtaDir, eCent, etaCoarseEta[etaKey], etaCoarseEtaErr[etaKey], /*isPhi*/false, "_zOnly");
            }
          }

          // -------- FINE z inside each etaKey --------
          for (const auto& etaKey : etaKeys)
          {
            const std::string baseEtaDir = outBaseDir + "/etaAndzDepBlocks/" + etaKey;
            EnsureDir(baseEtaDir);

            for (const auto& ztag : kZFineOrder)
            {
              TH3F* h_unc = getXZ(etaKey, ztag, /*corrected*/false);
              if (!h_unc) continue;
              h_unc->SetDirectory(nullptr);

              TH3F* h_cor = (!isFirstPass ? getXZ(etaKey, ztag, /*corrected*/true) : nullptr);
              if (h_cor) h_cor->SetDirectory(nullptr);

              const std::string outDirXZf = baseEtaDir + "/" + ztag;
              EnsureDir(outDirXZf);

              const char* prettyEta = kEtaPretty.count(etaKey) ? kEtaPretty.at(etaKey).c_str()
                                                               : etaKey.c_str();
              const char* prettyZ   = kZFinePretty.count(ztag) ? kZFinePretty.at(ztag).c_str()
                                                               : ztag.c_str();

              // Uncorrected products for each (η, z_fine)
              SaveTH3Lego(h_unc, (outDirXZf + "/lego_unc.png").c_str(), prettyEta);
              Plot2DBlockEtaPhi(h_unc, h_cor, isFirstPass, eEdges, outDirXZf.c_str(), prettyZ);
              OverlayUncorrPhiEta(h_unc, eEdges, outDirXZf.c_str(), prettyEta);
              BVecs bXZf = MakeBvaluesVsEnergyPlot(h_unc, eEdges, outDirXZf.c_str(), prettyEta);

              // Accumulate b(E) for per-η fine-z overlays
              etaFinePhi[etaKey][ztag]    = bXZf.bphi;
              etaFinePhiErr[etaKey][ztag] = bXZf.bphiErr;
              etaFineEta[etaKey][ztag]    = bXZf.beta;
              etaFineEtaErr[etaKey][ztag] = bXZf.betaErr;

              // b-values text labeling: (etaKey) × (ztag fine)
              if (bOut.is_open()) WriteBValuesTxt(h_unc, eEdges, bOut, etaKey.c_str(), ztag.c_str());

              // Corrected (second pass)
              if (h_cor && !isFirstPass)
              {
                SaveTH3Lego(h_cor, (outDirXZf + "/lego_cor.png").c_str(), prettyEta);
                FitLocalPhiEta(h_unc, h_cor, /*isFirstPassFlag*/false, eEdges, outDirXZf.c_str());
              }

              delete h_unc;
              if (h_cor) delete h_cor;
            }

            // After filling all fine z for this η, write the fine overlay (third image in the trio)
            if (!etaFinePhi[etaKey].empty())
            {
              SaveOverlayAcrossVariants(baseEtaDir, eCent, etaFinePhi[etaKey], etaFinePhiErr[etaKey], /*isPhi*/true ,  "_zFine");
              SaveOverlayAcrossVariants(baseEtaDir, eCent, etaFineEta[etaKey], etaFineEtaErr[etaKey], /*isPhi*/false, "_zFine");
            }
          }
        }
    }

    // Close b-values after all z-slice writes (coarse + fine)
    if (bOut.is_open()) bOut.close();

    // ----------------------- Global overlays across η-variants (base path) -----------------------
    SaveOverlayAcrossVariants(outBaseDir, eCent, phiByVariant, phiErrByVariant, /*isPhi*/true , ""  );
    SaveOverlayAcrossVariants(outBaseDir, eCent, etaByVariant, etaErrByVariant, /*isPhi*/false, ""  );

    // ==================== NEW: Step-4 — build δ(E,α) & (c0,m) ====================
    // Use the b(E) lines we just wrote in outBaseDir + "/bValues.txt"
    BuildDeltaAndWriteFits(fin.get(),
                           outBaseDir,
                           eEdges,
                           outBaseDir + "/bValues.txt",
                           "originalEta", "originalZRange",
                           /*E0=*/3.0);

    std::cout << "[DONE] Outputs written under: " << outBaseDir << "\n";
}

// Auto-run when invoked as a ROOT script:  root -b -q -l PDCAnalysisPrime.cpp
PDCAnalysisPRIMEPRIME();
