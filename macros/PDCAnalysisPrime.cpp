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
    h2->SetTitle(Form("Uncorr: E=[%.1f,%.1f)  (%s)", eLo, eHi, etaPretty));
    h2->GetZaxis()->SetTitle("Entries");
    h2->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    h2->GetYaxis()->SetTitle("block #varphi_{local, 2#times 2}");
    tuneAxes(h2);

    c2D_unc.cd(i+1);
    h2->Draw("COLZ");

    TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
    tl.SetTextSize(0.045); tl.SetTextAlign(22);
    tl.DrawLatex(0.50, 0.97, hdrTxt);
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
      h2->SetTitle(Form("Corr: E=[%.1f,%.1f)  (%s)", eLo, eHi, etaPretty));
      h2->GetZaxis()->SetTitle("Entries");
      h2->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
      h2->GetYaxis()->SetTitle("block #varphi_{local, 2#times 2}");
      tuneAxes(h2);

      c2D_cor.cd(i+1);
      h2->Draw("COLZ");

      TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
      tl.SetTextSize(0.045); tl.SetTextAlign(22);
      tl.DrawLatex(0.50, 0.97, hdrTxt);
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

      TString sliceTitle = Form("Uncorrected: %.1f < E < %.1f  GeV   (#scale[0.8]{%s})",
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
                            const char* variantKey)
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
        bOut << Form("PHI [%.1f,%.1f)  %.6f  %s\n", eLo, eHi, bPhi.val, variantKey);
    }

    // η
    hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());
    std::unique_ptr<TH1D> hEta( static_cast<TH1D*>(hUnc3D->Project3D("x")) );
    if (hEta) {
      BRes bEta = FitAsinh1D(hEta.get());
      if (std::isfinite(bEta.val))
        bOut << Form("ETA [%.1f,%.1f)  %.6f  %s\n", eLo, eHi, bEta.val, variantKey);
    }
  }
}


static void SaveOverlayAcrossVariants(
    const std::string& outBaseDir,
    const std::vector<double>& ecenters,
    const std::map<std::string, std::vector<double>>& byVar,
    const std::map<std::string, std::vector<double>>& byVarErr,
    bool isPhi)
{
  if (byVar.empty()) return;

  // ----------------------------- configuration -----------------------------
  const double E0   = 3.0;                       // reference energy (GeV)
  const double xLo  = E_edges[0] - 0.5;          // assume these globals exist
  const double xHi  = E_edges[N_E] - 0.5;

  auto prettyName = [&](const std::string& key)->const char*
  {
    auto it = kEtaPretty.find(key);
    return (it == kEtaPretty.end()) ? key.c_str() : it->second.c_str();
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
  const double yLo = ymin - 0.05 * std::fabs(ymin);
  const double yHi = ymax + 0.10 * std::fabs(ymax);

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

  TLatex lt; lt.SetNDC(true); lt.SetTextSize(0.024);
  const double xTxt = 0.45, yTop = 0.87, yStep = 0.050;

  std::vector<std::unique_ptr<TGraphErrors>> keepPoints;
  std::vector<std::unique_ptr<TF1>>          keepFits;

  // ----------------------------- terminal banner ----------------------------
  std::cout << "\n"
            << (isPhi ? "[PHI] " : "[ETA] ")
            << "Log-fit model:  b(E) = b0 + m * ln(E/E0),   E0 = " << E0 << " GeV\n"
            << "Weighted least squares per variant: w_i = 1/sigma_i^2 (if sigma provided) or 1 otherwise.\n"
            << "Parameter errors: unscaled if errors provided; scaled by χ²/ndf when unit weights are used.\n\n";

  std::cout << std::left
            << std::setw(18) << "variant"
            << std::right
            << std::setw(4)  << "N"
            << std::setw(12) << "E[min]"
            << std::setw(12) << "E[max]"
            << std::setw(12) << "b0@E0"
            << std::setw(10) << "±b0"
            << std::setw(12) << "m (ln E)"
            << std::setw(10) << "±m"
            << std::setw(8)  << "z(m)"
            << std::setw(12) << "m/decade"
            << std::setw(14) << "Δb_fit"
            << std::setw(14) << "Δb_data"
            << std::setw(8)  << "R^2"
            << std::setw(12) << "RMSE_w"
            << std::setw(10) << "χ²/ndf"
            << "\n"
            << "──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n";

  // ----------------------------- per-variant loop ---------------------------
  int idx = 0;
  for (const auto& kv : byVar)
  {
    const std::string& key   = kv.first;
    const auto&        vals  = kv.second;
    if (vals.size() != ecenters.size()) continue;

    // errors (if any)
    auto itErr = byVarErr.find(key);
    std::vector<double> errs = (itErr != byVarErr.end()) ? itErr->second
                                                         : std::vector<double>(vals.size(), 0.0);
    if (errs.size() != vals.size()) errs.assign(vals.size(), 0.0);

    // points
    std::unique_ptr<TGraphErrors> g(new TGraphErrors(static_cast<int>(vals.size())));
    double Emin = +1e9, Emax = -1e9;
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
    }
    const int col = colors[idx % 8];
    g->SetMarkerStyle(kMarkerStyle);
    g->SetMarkerSize(kMarkerSize);
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->Draw("P SAME");
    leg.AddEntry(g.get(), prettyName(key), "p");
    keepPoints.emplace_back(std::move(g));

    // --------------------- closed-form weighted least squares ----------------
    // model: y = b0 + m * x, where x = ln(E/E0)
    // weights: w = 1/sigma^2 when available, else w = 1
    double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
    int     Nuse = 0;
    bool    haveYerrs = false;

    std::vector<double> xi; xi.reserve(vals.size());
    std::vector<double> yi; yi.reserve(vals.size());
    std::vector<double> wi; wi.reserve(vals.size());

    for (size_t i = 0; i < vals.size(); ++i)
    {
      const double E  = ecenters[i];
      const double y  = vals[i];
      const double sy = errs[i];

      if (!(std::isfinite(E) && E > 0.0 && std::isfinite(y))) continue;

      const double x = std::log(E / E0);
      const double w = (std::isfinite(sy) && sy > 0.0) ? (haveYerrs = true, 1.0/(sy*sy)) : 1.0;

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

    const double Delta = S*Sxx - Sx*Sx;
    if (Nuse < 2 || !(Delta > 0.0)) { ++idx; continue; }

    const double m  = (S * Sxy - Sx * Sy) / Delta;
    const double b0 = (Sxx * Sy - Sx * Sxy) / Delta;

    // residual sums
    double chi2 = 0.0, TSS_w = 0.0, Wsum = 0.0;
    const double ybar_w = Sy / S;
    for (int i = 0; i < Nuse; ++i)
    {
      const double yhat = b0 + m * xi[i];
      const double r    = yi[i] - yhat;
      chi2  += wi[i] * r * r;             // weighted SSE (χ²)
      TSS_w += wi[i] * (yi[i] - ybar_w) * (yi[i] - ybar_w);
      Wsum  += wi[i];
    }
    const int    ndf     = std::max(0, Nuse - 2);
    const double redChi2 = (ndf > 0) ? (chi2 / ndf) : 0.0;
    const double R2_w    = (TSS_w > 0) ? (1.0 - chi2 / TSS_w) : 0.0;
    const double RMSE_w  = (Wsum  > 0) ? std::sqrt(chi2 / Wsum) : 0.0;

    // covariance: (X^T W X)^{-1}; scale by redChi2 only if unit weights were used
    const double invFac  = 1.0 / Delta;
    const double baseVar = haveYerrs ? 1.0 : redChi2;    // key fix vs earlier version
    const double var_b0  = baseVar * (Sxx * invFac);
    const double var_m   = baseVar * (S   * invFac);
    const double sb0     = (var_b0 > 0.0) ? std::sqrt(var_b0) : 0.0;
    const double sm      = (var_m  > 0.0) ? std::sqrt(var_m ) : 0.0;
    const double z_m     = (sm > 0.0) ? (m / sm) : 0.0;

    // slopes per decade and Δb across observed span
    const double m_per_dec = m * std::log(10.0);
    const double db_fit    = (b0 + m * std::log(Emax / E0))
                           - (b0 + m * std::log(Emin / E0));
    // empirical Δb from first/last points in this vector (no model)
    double db_data = 0.0;
    if (!vals.empty()) {
      const double y1 = vals.front();
      const double y2 = vals.back();
      if (std::isfinite(y1) && std::isfinite(y2)) db_data = y2 - y1;
    }

    // draw fitted curve
    const double fitLo = std::max(xLo, 0.95 * Emin);
    const double fitHi = std::min(xHi, 1.05 * Emax);
    std::unique_ptr<TF1> f(new TF1(Form("f_%s_%d", key.c_str(), idx),
                                   "[0] + [1]*log(x/[2])", fitLo, fitHi));
    f->SetParameters(b0, m, E0);
    f->FixParameter(2, E0);
    f->SetLineColor(col);
    f->SetLineWidth(2);
    f->Draw("SAME");
    keepFits.emplace_back(std::move(f));

      lt.SetTextColor(col);
      lt.DrawLatex(xTxt, yTop - idx*yStep,
                   Form("%s: b_{0}=%.4f,  m=%.4f,  #chi^{2}/ndf=%.2f",
                        prettyName(key), b0, m, redChi2));


    // terminal row
    std::cout << std::left  << std::setw(18) << prettyName(key)
              << std::right << std::setw(4)  << Nuse
              << std::setw(12) << std::fixed << std::setprecision(3) << Emin
              << std::setw(12) << std::fixed << std::setprecision(3) << Emax
              << std::setw(12) << std::fixed << std::setprecision(6) << b0
              << std::setw(10) << std::fixed << std::setprecision(6) << sb0
              << std::setw(12) << std::fixed << std::setprecision(6) << m
              << std::setw(10) << std::fixed << std::setprecision(6) << sm
              << std::setw(8)  << std::fixed << std::setprecision(2) << z_m
              << std::setw(12) << std::fixed << std::setprecision(6) << m_per_dec
              << std::setw(14) << std::fixed << std::setprecision(6) << db_fit
              << std::setw(14) << std::fixed << std::setprecision(6) << db_data
              << std::setw(8)  << std::fixed << std::setprecision(3) << R2_w
              << std::setw(12) << std::fixed << std::setprecision(6) << RMSE_w
              << std::setw(10) << std::fixed << std::setprecision(3) << redChi2
              << "\n";

    ++idx;
  }

  std::cout << "──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n"
            << "(E0 fixed at " << E0 << " GeV; m/decade = m * ln 10. "
            << "Errors unscaled when σ_y provided; otherwise scaled by χ²/ndf.)\n\n";

  leg.Draw();

  TString out = Form("%s/%s", outBaseDir.c_str(),
                     isPhi ? "bValuesPhiOverlay.png" : "bValuesEtaOverlay.png");
  c.SaveAs(out);
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
                      << " – no non-empty bin #" << nth << " for " << altShort << " vs RAW-CP.\n";
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
                      << " – no non-empty bin #" << nth << " with RAW-CP, CP-corr, CP-corr(EA).\n";
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
                      << " – no non-empty bin #" << nth << ".\n";
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
                      << " – no non-empty bin #" << nth << ".\n";
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

    // ---------------------------------------------------------------------
}



void MakeDeltaPhiEtaPlayground(
    const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root",
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

        // Consistent margins for all panels (increase bottom margin so x–axis title doesn't clip)
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

        // Decide x–axis title based on what we're overlaying:
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

        // Make sure the x–axis title sits fully inside the pad
        base->GetXaxis()->SetTitleSize(0.045 * textScale);
        base->GetXaxis()->SetLabelSize(0.035 * textScale);
        base->GetXaxis()->SetTitleOffset(1.10);

        base->Draw("E");
        for (std::size_t i=1; i<clones.size(); ++i) clones[i]->Draw("E SAME");

        // Top-left: Plot name + energy
        TLatex head; head.SetNDC(); head.SetTextFont(42);
        head.SetTextAlign(13); head.SetTextSize(0.033 * textScale);
        head.DrawLatex(0.17, 0.89, Form("%s, [%.0f, %.0f) GeV", plotName.c_str(), eLo, eHi));

        // Top-right legend (points only, no error bars) — keep alive on the pad
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
                      << " – skipped PNG #1.\n";
        }


        // ---------- PNG #2 : 4×2 table (this residual only) ----------
        {
            const int nCol=4, nRow=2, maxPads=nCol*nRow;
            TCanvas cOT(Form("cPlay_%s_Table", stem),
                        Form("%s – all slices (Clusterizer Output, No Position Correction)", sym),
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
        c.SetBottomMargin(0.18); // was 0.13: prevents x–axis title clipping

        TH1F fr("fr","",1,xMinE,xMaxE);
        fr.SetTitle("");

        // Axis titles and labels
        fr.GetXaxis()->SetTitle("E  [GeV]");
        fr.GetYaxis()->SetTitle(yTitle); // "#sigma_{RMS}  [rad]"

        // Ensure the x–axis title fits inside the pad
        fr.GetXaxis()->SetTitleSize(0.050);
        fr.GetXaxis()->SetLabelSize(0.038);
        fr.GetXaxis()->SetTitleOffset(1.05); // mild offset keeps it clear but inside the pad

        // Keep y–axis sizing consistent
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
        h.DrawLatex(0.50,0.985,"no corr, cluster  vs  CorrectPosition, cluster  –  Mean / RMS vs E");


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

        // CSV: 4 series (phi raw/corr, eta raw/corr) – cluster variants
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
            "no corr, scratch  vs  b(E) corr, scratch  –  Mean / RMS vs E");

    if (!HetaScratchRaw.empty() && !HetaScratchCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HetaScratchRaw, HetaScratchCorr,
            kBlue+1, 20, 24, "#mu",
            (dirScratchRawVsCorr + "/DeltaEta_raw_vs_corr_MeanSigmaVsE_Overlay.png").c_str(),
            "no corr, scratch", "b(E) corr, scratch",
            "no corr, scratch  vs  b(E) corr, scratch  –  Mean / RMS vs E");

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

        TCanvas c("cFSvClus_4Series","FS vs Clusterizer – four series μ/σ(E)",1000,850);
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

    // Δφ — first‑bin overlay (counts)
    if (!HphiScratchRaw.empty() && !Hphi.empty()
        && HphiScratchRaw[0] && Hphi[0]
        && HphiScratchRaw[0]->Integral()>0 && Hphi[0]->Integral()>0)
    {
        TCanvas cFSvCl_Phi_First("cFSvCl_Phi_First_NoCorr",
                                 "#Delta #phi from-scratch vs clusterizer (no corr) – first bin",
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

    // Δη — first‑bin overlay (counts)
    if (!HetaScratchRaw.empty() && !Heta.empty()
        && HetaScratchRaw[0] && Heta[0]
        && HetaScratchRaw[0]->Integral()>0 && Heta[0]->Integral()>0)
    {
        TCanvas cFSvCl_Eta_First("cFSvCl_Eta_First_NoCorr",
                                 "#Delta #eta from-scratch vs clusterizer (no corr) – first bin",
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

    // Δφ — 4×2 table overlay (all energy bins)
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Phi_Table_NoCorr",
                   "Δφ from-scratch vs clusterizer (no corr) – table", 1600, 900);
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

    // Δη — 4×2 table overlay (all energy bins)
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Eta_Table_NoCorr",
                   "#Delta #eta from-scratch vs clusterizer (no corr) – table", 1600, 900);
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
        drawHeaderNDC("#Delta#eta  –  From-scratch vs Clusterizer  (No Position Correction)", 0.975, 0.045);
        cT.SaveAs((dirFSvsClusNoCorr + "/DeltaEta_FSvsClus_AllBins_Table.png").c_str());
        cT.Close();
    }

    // Δφ — μ/σ(E) overlay (two series): from‑scratch vs clusterizer, no corr
    if (!HphiScratchRaw.empty() && !Hphi.empty())
        makeMuSigmaVariantOverlayWithLabels(HphiScratchRaw, Hphi,
            kRed+1, 20, 24, "#mu",
            (dirFSvsClusNoCorr + "/DeltaPhi_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch", "clusterizer",
            "#Delta#phi  -  from scratch vs clusterizer  (No Position Correction)  –  Mean / RMS vs E");

    // Δη — μ/σ(E) overlay (two series): from‑scratch vs clusterizer, no corr
    if (!HetaScratchRaw.empty() && !Heta.empty())
        makeMuSigmaVariantOverlayWithLabels(HetaScratchRaw, Heta,
            kBlue+1, 20, 24, "#mu",
            (dirFSvsClusNoCorr + "/DeltaEta_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch", "clusterizer",
            "#Delta#eta  -  from scratch vs clusterizer  (No Position Correction)  –  Mean / RMS vs E");

    // NEW: 4-series μ/σ(E) overlay (from‑scratch vs clusterizer), no corr
    if (!HphiScratchRaw.empty() && !Hphi.empty()
        && !HetaScratchRaw.empty() && !Heta.empty())
        makeMuSigmaFourSeriesOverlayFSvClus(HphiScratchRaw, Hphi,
                                            HetaScratchRaw, Heta,
            (dirFSvsClusNoCorr + "/DeltaEtaPhi_FSvsClus_MeanSigmaVsE_Overlay_4Series.png").c_str(),
            "#Delta#phi / #Delta#eta  -  from scratch vs clusterizer  (No Position Correction)  -  Mean / RMS vs E");


    // Δφ — first‑bin overlay (counts), corrected
    if (!HphiScratchCorr.empty() && !HphiCorr.empty()
        && HphiScratchCorr[0] && HphiCorr[0]
        && HphiScratchCorr[0]->Integral()>0 && HphiCorr[0]->Integral()>0)
    {
        TCanvas cFSvCl_Phi_FirstC("cFSvCl_Phi_First_WithCorr",
                                  "#Delta #phi from-scratch vs clusterizer (with corr) – first bin",
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

    // Δη — first‑bin overlay (counts), corrected
    if (!HetaScratchCorr.empty() && !HetaCorr.empty()
        && HetaScratchCorr[0] && HetaCorr[0]
        && HetaScratchCorr[0]->Integral()>0 && HetaCorr[0]->Integral()>0)
    {
        TCanvas cFSvCl_Eta_FirstC("cFSvCl_Eta_First_WithCorr",
                                  "#Delta #eta from-scratch vs clusterizer (with corr) – first bin",
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

    // Δφ — 4×2 table overlay (all energy bins), corrected
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Phi_Table_WithCorr",
                   "#Delta #phi from-scratch vs clusterizer (with corr) – table", 1600, 900);
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

    // Δη — 4×2 table overlay (all energy bins), corrected
    {
        const int nCol=4, nRow=2, maxPads=nCol*nRow;
        TCanvas cT("cFSvCl_Eta_Table_WithCorr",
                   "Δη from-scratch vs clusterizer (with corr) – table", 1600, 900);
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

    // Δφ — μ/σ(E) overlay (two series): from‑scratch vs clusterizer, with corr
    if (!HphiScratchCorr.empty() && !HphiCorr.empty())
        makeMuSigmaVariantOverlayWithLabels(HphiScratchCorr, HphiCorr,
            kRed+1, 20, 24, "#mu",
            (dirFSvsClusWithCorr + "/DeltaPhi_FSvsClus_MeanSigmaVsE_Overlay.png").c_str(),
            "from scratch (corr)", "clusterizer (corr)",
            "#Delta#phi  -  from scratch vs clusterizer  (With Position Correction)  –  Mean / RMS vs E");


    // Δη — μ/σ(E) overlay (two series): from‑scratch vs clusterizer, with corr
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
        TCanvas cOT("cPlay_Overlay_Table","DeltaEta/DeltaPhi overlay – all slices",1600,900);
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
        TCanvas cRT("cPlay_Ratio_Table","#Delta#phi / #Delta#eta – ratio (no corr, cluster) – all slices",1600,900);
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

            TCanvas cS("cPlay_Overlay_MuSig","DeltaEta/DeltaPhi overlay – Mean / RMS vs E",900,800);

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
    //  Δφ variants – μ/σ vs E overlays + μ vs ln(E) fits (as before)
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
            "h_phi_diff_cpCorrEA_geom_%s",             // 5: CLUS-CP(EA geom)
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

            TCanvas c("cPhiVar_MuSig","DeltaPhi – variants μ/σ vs E",1000,850);
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
            h.DrawLatex(0.50,0.985,"#Delta#phi  –  Mean / RMS vs E  (variants)");
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

                TCanvas c3("cPhiVar_MuSig_3","DeltaPhi – CLUS RAW/CORR/CORR(EA_geom) μ/σ vs E",1000,850);
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
                h3.DrawLatex(0.50,0.985,"#Delta#phi  –  Clusterizer RAW / CORR / CORR(EA_{geom})  –  Mean / RMS vs E");

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

                // Append fits for 3-way (RAW/CP/EA_geom)
                std::ofstream fapp(dirPhiLog + "/DeltaPhi_MuVsLogE_fit.txt", std::ios::app);
                if (fapp.tellp() > 0) fapp << "\n# ---- Clusterizer 3-way (cpRaw / cpCorr / cpCorrEA_geom) ----\n";

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

                    fapp << Form("%-44s % .6e   % .6e\n",
                                 phiLab[v], f.GetParameter(0), f.GetParameter(1));

                    lg3.AddEntry(g, Form("%s  (m=%.2e, b=%.2e)",
                                         phiLab[v], f.GetParameter(1), f.GetParameter(0)), "p");
                }
                lg3.Draw();
                c3ln.SaveAs((dirPhiLog + "/DeltaPhi_MuVsLogE_Clusterizer3.png").c_str());
                c3ln.Close();

                // ---------------- EXTRA: μ vs ln(E) fits for all EA flavours ----------------
                std::vector<int> eaList;
                for (int v : {vEA_geom, vEA_etaDep, vEA_eOnly, vEA_mix}) if (hasAny(v)) eaList.push_back(v);
                if (!eaList.empty()) {
                    std::ofstream fea(dirPhiLog + "/DeltaPhi_MuVsLogE_fit.txt", std::ios::app);
                    fea << "\n# ---- CLUS-CP(EA) 4-way (geom / |eta|+E / E-only / mix) ----\n";
                    for (int v : eaList) {
                        TF1 f(Form("fEA_%d",v), "pol1", xLo, xHi);
                        // make a throwaway graph just to run the fit consistently (no draw)
                        TGraphErrors gEA((int)eCtr.size(),
                                         lnE.data(), MU[v].data(),
                                         ex0.data(), dMU[v].data());
                        gEA.Fit(&f, "Q");
                        fea << Form("%-44s % .6e   % .6e\n",
                                    phiLab[v], f.GetParameter(0), f.GetParameter(1));
                    }
                }
            }
        }

        // μ vs ln(E) – FOUR variants (CLUS-RAW, CLUS-CORR, PDC-RAW, PDC-CORR)
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

                TCanvas c4("cPhiVar_MuSig_4","DeltaPhi – four variants μ/σ vs E",1000,850);
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
    std::cout << "[Playground] Completed Δφ & Δη outputs into: " << outDir << "\n";
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



// ================================ MAIN ======================================

void PDCAnalysisPrime()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // ——— paths
  const std::string inFilePath = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root";
  const std::string outBaseDir = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutputPrime";

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

  // Prepare global b-values text file
  const std::string bValFile = outBaseDir + "/bValues.txt";
  std::ofstream bOut(bValFile);
  if (!bOut.is_open())
  {
    std::cerr << "[WARN] Cannot open " << bValFile << " => won't save b-values.\n";
  }
  else
  {
    bOut << "# E range   best-b   (PHI or ETA)   etaVariant\n";
  }

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

  for (const auto& v : variants)
    {
      const std::string outDir = outBaseDir + "/" + v.key;
      EnsureDir(outDir);

      // --- run your normal per-variant products (lego, 2D, overlays, b vs E) ---
      TH3F* hUnc = GetTH3FByNames(fin.get(), v.uncNames);
      TH3F* hCor = nullptr;
      if (!isFirstPass) hCor = GetTH3FByNames(fin.get(), v.corNames);
      if (!hUnc) { std::cerr << "[WARN] missing " << v.key << " – skipped\n"; continue; }
      hUnc->SetDirectory(nullptr); if (hCor) hCor->SetDirectory(nullptr);

      SaveTH3Lego(hUnc, (outDir + "/lego_unc.png").c_str(), kEtaPretty.at(v.key).c_str());
      Plot2DBlockEtaPhi(hUnc, hCor, isFirstPass, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
      OverlayUncorrPhiEta (hUnc, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
      BVecs bv = MakeBvaluesVsEnergyPlot(hUnc, eEdges, outDir.c_str(), kEtaPretty.at(v.key).c_str());
      phiByVariant[v.key]    = bv.bphi;  phiErrByVariant[v.key] = bv.bphiErr;
      etaByVariant[v.key]    = bv.beta;  etaErrByVariant[v.key] = bv.betaErr;
      if (bOut.is_open()) WriteBValuesTxt(hUnc, eEdges, bOut, v.key.c_str());
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

      delete hUnc; if (hCor) delete hCor;
  }

  if (bOut.is_open()) bOut.close();

  // Global overlays across variants (in base path)
  SaveOverlayAcrossVariants(outBaseDir, eCent, phiByVariant, phiErrByVariant, /*isPhi*/true);
  SaveOverlayAcrossVariants(outBaseDir, eCent, etaByVariant, etaErrByVariant, /*isPhi*/false);


  AnalyzeBvsResCLUSCPEA(inFilePath.c_str(), outBaseDir.c_str());
  std::cout << "[DONE] Outputs written under: " << outBaseDir << "\n";
}

// Auto-run when invoked as a ROOT script:  root -b -q -l PDCAnalysisPrime.cpp
PDCAnalysisPrime();
