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

static inline void EnsureDir(const std::string& path)
{
  gSystem->mkdir(path.c_str(), /*recursive*/true);
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


// ===================== Overlay b(E) across η-variants (with errors) ========

static void SaveOverlayAcrossVariants(
    const std::string& outBaseDir,
    const std::vector<double>& ecenters,
    const std::map<std::string, std::vector<double>>& byVar,
    const std::map<std::string, std::vector<double>>& byVarErr,
    bool isPhi)
{
  if (byVar.empty()) return;

  // y-range across all variants, include ±σ
  double ymin = +1e9, ymax = -1e9;
  for (const auto& kv : byVar)
  {
    const std::string& key = kv.first;
    const auto& vals = kv.second;
    const auto itErr = byVarErr.find(key);
    const std::vector<double>* errs = (itErr == byVarErr.end()) ? nullptr : &itErr->second;

    for (size_t i = 0; i < vals.size(); ++i)
    {
      if (!std::isfinite(vals[i])) continue;
      const double e = (errs && i < errs->size()) ? (*errs)[i] : 0.0;
      ymin = std::min(ymin, vals[i] - e);
      ymax = std::max(ymax, vals[i] + e);
    }
  }
  if (!std::isfinite(ymin) || !std::isfinite(ymax)) { ymin = 0.0; ymax = 1.0; }

  const double xLo = E_edges[0]-0.5, xHi = E_edges[N_E]-0.5;
  const double yLo = ymin - 0.05*std::fabs(ymin);
  const double yHi = ymax + 0.1*std::fabs(ymax);

  TCanvas c("cbOverlay", isPhi ? "b_{#varphi}(E) overlay" : "b_{#eta}(E) overlay", 900, 700);
  TH2F frame("frame", "", 100, xLo, xHi, 100, yLo, yHi);
  frame.SetTitle(isPhi ? "best-fit  b_{#varphi}  vs  E;E  [GeV];b_{#varphi}"
                       : "best-fit  b_{#eta}     vs  E;E  [GeV];b_{#eta}");
  frame.SetStats(0);
  frame.Draw();

  const int colors[5] = { kBlack, kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
  const int kMarkerStyle = 20;
  const double kMarkerSize = 1.2;

  TLegend leg(0.16, 0.65, 0.40, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.040);

  int idx = 0;
  std::vector<std::unique_ptr<TGraphErrors>> keep; // keep graphs alive

  for (const auto& kv : byVar)
  {
    const std::string& key = kv.first;
    const std::vector<double>& vals = kv.second;

    auto itPretty = kEtaPretty.find(key);
    const char* pretty = (itPretty == kEtaPretty.end() ? key.c_str() : itPretty->second.c_str());

    auto itErr = byVarErr.find(key);
    if (itErr == byVarErr.end()) continue;  // require errors present
    const std::vector<double>& errs = itErr->second;

    if (vals.size() != ecenters.size() || errs.size() != vals.size()) continue;

    std::unique_ptr<TGraphErrors> g(new TGraphErrors(static_cast<int>(vals.size())));
    for (size_t i = 0; i < vals.size(); ++i)
    {
      g->SetPoint(static_cast<int>(i), ecenters[i], vals[i]);
      g->SetPointError(static_cast<int>(i), 0.0, errs[i]); // no x-errors
    }

    const int col = colors[idx % 5];
    g->SetMarkerStyle(kMarkerStyle);
    g->SetMarkerSize(kMarkerSize);
    g->SetMarkerColor(col);
    g->SetLineColor(col);

    g->Draw("P SAME");            // points + symmetric y-errors by default
    leg.AddEntry(g.get(), pretty, "p");

    keep.emplace_back(std::move(g));
    ++idx;
  }

  leg.Draw();

  TString out = Form("%s/%s", outBaseDir.c_str(),
                     isPhi ? "bValuesPhiOverlay.png" : "bValuesEtaOverlay.png");
  c.SaveAs(out);
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

  // Loop over variants / folders
  for (const auto& v : variants)
  {
    // Output folder and pretty label
    const std::string outDir = outBaseDir + "/" + v.key;
    EnsureDir(outDir);
    const char* pretty = kEtaPretty.count(v.key) ? kEtaPretty.at(v.key).c_str() : v.key.c_str();

    // Retrieve uncorrected histogram
    TH3F* hUnc = GetTH3FByNames(fin.get(), v.uncNames);
    if (!hUnc)
    {
      std::cerr << "[WARN] Missing uncorrected TH3 for variant '" << v.key << "'. Skipping.\n";
      continue;
    }
    hUnc->SetDirectory(nullptr); // detach

    // Retrieve corrected histogram only for second pass
    TH3F* hCor = nullptr;
    if (!isFirstPass)
    {
      hCor = GetTH3FByNames(fin.get(), v.corNames);
      if (hCor) hCor->SetDirectory(nullptr);
    }

    // 1) lego_unc
    SaveTH3Lego(hUnc, (outDir + "/lego_unc.png").c_str(), pretty);

    // 2) 2D XY tables per E-bin (and corrected if second pass)
    Plot2DBlockEtaPhi(hUnc, hCor, isFirstPass, eEdges, outDir.c_str(), pretty);

    // 3) Overlay of uncorrected φ & η (per E-bin)
    OverlayUncorrPhiEta(hUnc, eEdges, outDir.c_str(), pretty);

    // 4) Per-variant bφ(E) & bη(E) + figure bValues_vs_E.png
    BVecs bv = MakeBvaluesVsEnergyPlot(hUnc, eEdges, outDir.c_str(), pretty);
    // store for global overlays (values + errors)
    phiByVariant    [v.key] = bv.bphi;
    etaByVariant    [v.key] = bv.beta;
    phiErrByVariant [v.key] = bv.bphiErr;
    etaErrByVariant [v.key] = bv.betaErr;

    // 5) Append to global text file
    if (bOut.is_open())
      WriteBValuesTxt(hUnc, eEdges, bOut, v.key.c_str());

    // 6) Second-pass only extras
    if (!isFirstPass)
    {
      if (hCor) SaveTH3Lego(hCor, (outDir + "/lego_cor.png").c_str(), pretty);
      FitLocalPhiEta(hUnc, hCor, /*isFirstPassFlag*/false, eEdges, outDir.c_str());
    }

    // Clean up pointers that were detached
    delete hUnc;
    if (hCor) delete hCor;
  }

  if (bOut.is_open()) bOut.close();

  // Global overlays across variants (in base path)
  SaveOverlayAcrossVariants(outBaseDir, eCent, phiByVariant, phiErrByVariant, /*isPhi*/true);
  SaveOverlayAcrossVariants(outBaseDir, eCent, etaByVariant, etaErrByVariant, /*isPhi*/false);


  std::cout << "[DONE] Outputs written under: " << outBaseDir << "\n";
}

// Auto-run when invoked as a ROOT script:  root -b -q -l PDCAnalysisPrime.cpp
PDCAnalysisPrime();
