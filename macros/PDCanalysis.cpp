#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TIterator.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TFitResult.h>
#include <TLatex.h>
#include <vector>
#include <utility>
#include <cmath>
#include <memory>
#include <cfloat>  // for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// 1) New energy‐bin array => 8 intervals
////////////////////////////////////////////////////////////////////////////////
constexpr double E_edges[] = {2, 4, 6, 8, 10, 12, 15, 20, 30};
constexpr int    N_E       = (sizeof(E_edges)/sizeof(E_edges[0])) - 1;
// => N_E=8

////////////////////////////////////////////////////////////////////////////////
// (A) "core Gaussian fit" method
////////////////////////////////////////////////////////////////////////////////
double coreGaussianSigma(TH1* h, const TString& pngSavePath)
{
  if(!h || h->GetEntries() < 20) return 0.;

  // quartiles [25%,50%,75%]
  double q[3], probs[3] = {0.25, 0.50, 0.75};
  h->GetQuantiles(3, q, probs);

  double xMin = q[0];
  double xMax = q[2];
  if(xMin >= xMax) return 0.;

  TF1 fitG("fitG","gaus", xMin, xMax);
  fitG.SetParameters(h->GetMaximum(), q[1], 0.3*(xMax-xMin));

  int fitStatus = h->Fit(&fitG, "QN0R");
  double sigma  = (fitStatus==0) ? fitG.GetParameter(2) : h->GetRMS();

  // If we want the visual diagnostic:
  if(!pngSavePath.IsNull())
  {
    TCanvas ctmp("ctmp","ctmp",600,500);
    h->Draw("E");

    fitG.SetLineColor(kRed);
    fitG.SetLineWidth(2);
    fitG.Draw("SAME");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.04);

    double fitMean = fitG.GetParameter(1);
    lat.SetTextColor(kRed);
    lat.DrawLatex(0.58, 0.82, Form("#mu=%.3f", fitMean));
    lat.DrawLatex(0.58, 0.74, Form("#sigma=%.3f", sigma));

    ctmp.SaveAs(pngSavePath);
    ctmp.Close();
  }

  return sigma;
}


////////////////////////////////////////////////////////////////////////////////
// (B) "raw RMS" method
////////////////////////////////////////////////////////////////////////////////
double rawRMS(TH1* h, const TString& /*unusedPNG*/)
{
  // We won't produce a diagnostic figure for raw RMS, or you can if you want.
  // Just return h->GetRMS().
  if(!h || h->GetEntries()<2) return 0.;
  return h->GetRMS();
}


////////////////////////////////////////////////////////////////////////////////
// Helper structure to store the results of scanning (Ash or Log)
////////////////////////////////////////////////////////////////////////////////
struct ScanResults
{
  // For each E‐bin, we store a TGraph of (parameter vs sigma)
  std::vector<TGraph*> tgVec;

  // In each E‐bin, we track best parameter + best sigma
  std::vector<double> bestParam;  // best b (cm) or best w0
  std::vector<double> bestSigma;

  // global min/max for the Y-axis across all E-bins
  double minY=DBL_MAX, maxY=-DBL_MAX;

  // counters
  int totalHist=0, missingHist=0, zeroEntry=0, usedHist=0;
};

////////////////////////////////////////////////////////////////////////////////
// doAshScan(...): scans b in [0.05..0.60], reads histograms, calls sigmaFunc
////////////////////////////////////////////////////////////////////////////////
ScanResults doAshScan(
    TFile* f,
    int N_E,
    const double* E_edges,
    const std::vector<double>& bScan,
    double cellSize,
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir
)
{
  ScanResults results;
  results.tgVec.resize(N_E);

  for(int iE=0; iE<N_E; iE++)
  {
    auto g = new TGraph;
    g->SetName( Form("gAsh_%s_E%d", suffix.Data(), iE) );

    for(double bValRaw : bScan)
    {
      // replicate 4-dec rounding
      double bValCell = std::round(bValRaw*10000.)/10000.;

      results.totalHist++;

      // histogram name => "h_dx_ash_b%.4f_E{iE}"
      TString hName = Form("h_dx_ash_b%.4f_E%d", bValCell, iE);
      TH1* h = dynamic_cast<TH1*>( f->Get(hName) );
      if(!h)
      {
        results.missingHist++;
        continue;
      }

      double ent = h->GetEntries();
      if(ent <= 0) {
        results.zeroEntry++;
      } else {
        results.usedHist++;
      }

      // create a little name for the optional PNG
      // e.g. "h_dx_ash_b0.2000_E2_fit.png"
      TString diagName = Form("%s/%s_%s.png", histOutDir.Data(), hName.Data(), suffix.Data());
      double sVal = sigmaFunc(h, diagName);

      double bCm = bValCell * cellSize;

      g->SetPoint(g->GetN(), bCm, sVal);
      if(sVal>0)
      {
        results.minY = std::min(results.minY, sVal);
        results.maxY = std::max(results.maxY, sVal);

        if(iE >= (int)results.bestSigma.size()) {
          results.bestSigma.resize(N_E, DBL_MAX);
          results.bestParam.resize(N_E, 0.);
        }
        if(sVal < results.bestSigma[iE]) {
          results.bestSigma[iE]  = sVal;
          results.bestParam[iE]  = bCm;
        }
      }
    } // end for bVal

    g->SetMarkerStyle(20 + (iE % 10));
    g->SetMarkerColor(1 + (iE % 10));
    g->SetLineColor  (1 + (iE % 10));
    g->SetMarkerSize(1.2);

    results.tgVec[iE] = g;
  } // end for iE

  if(results.minY == DBL_MAX) { results.minY=0.; results.maxY=1.; }

  return results;
}


////////////////////////////////////////////////////////////////////////////////
// doLogScan(...): scans w0 in [2..6], reads histos, calls sigmaFunc
////////////////////////////////////////////////////////////////////////////////
ScanResults doLogScan(
    TFile* f,
    int N_E,
    const double* E_edges,
    const std::vector<double>& w0Scan,
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir
)
{
  ScanResults results;
  results.tgVec.resize(N_E);

  for(int iE=0; iE<N_E; iE++)
  {
    auto g = new TGraph;
    g->SetName( Form("gLog_%s_E%d", suffix.Data(), iE) );

    for(double w0Raw : w0Scan)
    {
      double wVal = std::round(w0Raw*100.)/100.;

      results.totalHist++;

      // name => "h_dx_log_w0%.2f_E{iE}"
      TString hName = Form("h_dx_log_w0%.2f_E%d", wVal, iE);
      TH1* h = dynamic_cast<TH1*>( f->Get(hName) );
      if(!h)
      {
        results.missingHist++;
        continue;
      }

      double ent = h->GetEntries();
      if(ent<=0) {
        results.zeroEntry++;
      } else {
        results.usedHist++;
      }

      TString diagName = Form("%s/%s_%s.png", histOutDir.Data(), hName.Data(), suffix.Data());
      double sVal = sigmaFunc(h, diagName);

      g->SetPoint(g->GetN(), wVal, sVal);
      if(sVal>0)
      {
        results.minY = std::min(results.minY, sVal);
        results.maxY = std::max(results.maxY, sVal);

        if(iE >= (int)results.bestSigma.size()) {
          results.bestSigma.resize(N_E, DBL_MAX);
          results.bestParam.resize(N_E, 0.);
        }
        if(sVal < results.bestSigma[iE]) {
          results.bestSigma[iE] = sVal;
          results.bestParam[iE] = wVal;
        }
      }
    } // end w0Raw

    g->SetMarkerStyle(21 + (iE % 10));
    g->SetMarkerColor(2 + (iE % 10));
    g->SetLineColor  (2 + (iE % 10));
    g->SetMarkerSize(1.2);

    results.tgVec[iE] = g;
  } // end iE

  if(results.minY==DBL_MAX){ results.minY=0.; results.maxY=1.; }
  return results;
}


////////////////////////////////////////////////////////////////////////////////
// drawAshLogSideBySide(...):
//   Takes 2 ScanResults: one for Ash, one for Log, plus “method name” (e.g. "fit" or "rms")
//   => draws them on a single 1×2 canvas (left=Ash, right=Log) and saves a PNG
////////////////////////////////////////////////////////////////////////////////
void drawAshLogSideBySide(
    const ScanResults& ashRes,
    const ScanResults& logRes,
    const char* methodName,  // e.g. "FIT" or "RMS"
    const double* E_edges,
    int N_E,
    const TString& baseDir
)
{
  // 1) Make the 1×2 canvas
  TString canName = Form("cSideBySide_%s", methodName);
  TCanvas cSide(canName, canName, 1600, 600);
  cSide.Divide(2,1);

  // 2) Left pad => Ash
  cSide.cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  // if we have at least 1 E bin
  if(!ashRes.tgVec.empty())
  {
    ashRes.tgVec[0]->Draw("ALP");
    ashRes.tgVec[0]->GetXaxis()->SetTitle("b (cm)");
    ashRes.tgVec[0]->GetYaxis()->SetTitle("#sigma_{x} (cm)");
    ashRes.tgVec[0]->GetYaxis()->SetRangeUser(ashRes.minY - 0.1*fabs(ashRes.minY),
                                              ashRes.maxY + 0.1*fabs(ashRes.maxY));

    for(size_t i=1; i<ashRes.tgVec.size(); i++){
      ashRes.tgVec[i]->Draw("LP SAME");
    }
  }

  TLegend legA(0.15,0.65,0.48,0.88);
  legA.SetBorderSize(0);
  legA.SetFillStyle(0);
  for(int iE=0;iE<N_E;iE++){
    double bBest = (iE<(int)ashRes.bestParam.size()) ? ashRes.bestParam[iE] : 0.;
    double sBest = (iE<(int)ashRes.bestSigma.size()) ? ashRes.bestSigma[iE] : 0.;
    legA.AddEntry(
      (iE<(int)ashRes.tgVec.size()? ashRes.tgVec[iE] : nullptr),
      Form("%.1f< E<%.1f (b=%.3f, #sigma=%.3f)", E_edges[iE], E_edges[iE+1], bBest, sBest),
      "lp"
    );
  }
  legA.Draw();

  // top label => e.g. "Ash [FIT]" or "Ash [RMS]"
  TLatex latA; latA.SetNDC(true);
  latA.DrawLatex(0.2, 0.92, Form("Ash [%s]", methodName));

  // 3) Right pad => Log
  cSide.cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  if(!logRes.tgVec.empty())
  {
    logRes.tgVec[0]->Draw("ALP");
    logRes.tgVec[0]->GetXaxis()->SetTitle("w_{0}");
    logRes.tgVec[0]->GetYaxis()->SetTitle("#sigma_{x} (cm)");
    logRes.tgVec[0]->GetYaxis()->SetRangeUser(logRes.minY - 0.1*fabs(logRes.minY),
                                              logRes.maxY + 0.1*fabs(logRes.maxY));
    for(size_t i=1;i<logRes.tgVec.size(); i++){
      logRes.tgVec[i]->Draw("LP SAME");
    }
  }

  TLegend legB(0.15,0.65,0.48,0.88);
  legB.SetBorderSize(0);
  legB.SetFillStyle(0);
  for(int iE=0;iE<N_E;iE++){
    double wBest = (iE<(int)logRes.bestParam.size()) ? logRes.bestParam[iE] : 0.;
    double sBest = (iE<(int)logRes.bestSigma.size()) ? logRes.bestSigma[iE] : 0.;
    legB.AddEntry(
      (iE<(int)logRes.tgVec.size()? logRes.tgVec[iE] : nullptr),
      Form("%.1f< E<%.1f (w0=%.2f, #sigma=%.3f)", E_edges[iE], E_edges[iE+1], wBest, sBest),
      "lp"
    );
  }
  legB.Draw();

  TLatex latB; latB.SetNDC(true);
  latB.DrawLatex(0.2, 0.92, Form("Log [%s]", methodName));

  // 4) Save the 1×2 canvas
  TString outName = Form("%s/SideBySide_%s.png", baseDir.Data(), methodName);
  cSide.SaveAs(outName);
  std::cout << "[INFO] Wrote combined (Ash & Log) canvas => " << outName << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// The main function that does 2 methods => “fit” & “rms”
//   Each produces 1 row × 2 columns => left=Ash, right=Log
////////////////////////////////////////////////////////////////////////////////
void plotAshLogRMS_sideBySide(const char* infile="PositionDep_sim_ALL.root")
{
  // b-scan
  const double bMin=0.05, bMax=0.60, bStep=0.01;
  std::vector<double> bScan;
  for(double b=bMin; b<=bMax+1e-9; b+=bStep) bScan.push_back(b);

  // log-scan
  const double w0Min=2.0, w0Max=6.0, w0Step=0.1;
  std::vector<double> w0Scan;
  for(double w=w0Min; w<=w0Max+1e-9; w+=w0Step) w0Scan.push_back(w);

  // output
  TString baseDir    = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
  TString histOutDir = baseDir + "/ASH_LOG_PLOTS";
  gSystem->mkdir(baseDir,true);
  gSystem->mkdir(histOutDir,true);

  // 2) Open file
  std::unique_ptr<TFile> f(TFile::Open(infile,"READ"));
  if(!f || f->IsZombie()){
    std::cerr << "[ERROR] Could not open file => abort.\n";
    return;
  }
  std::cout << "[INFO] Successfully opened '" << infile << "'\n";

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);

  constexpr double cellSize=5.55; // cm

  // We'll define two "sigma" functions:
  //   (A) coreGaussianSigma
  //   (B) rawRMS  (we can create a small lambda that just calls rawRMS)
  auto rawRMSlambda = [&](TH1* h, const TString&) {
    if(!h || h->GetEntries()<2) return 0.;
    return h->GetRMS();
  };

  // 3) (FIT) doAshScan & doLogScan => produce (AshFit, LogFit)
  std::cout<<"\n=== Doing FIT-based scanning (coreGaussianSigma) ===\n";
  auto ashFit = doAshScan(f.get(), N_E, E_edges, bScan, cellSize, coreGaussianSigma, "fit", histOutDir);
  auto logFit = doLogScan(f.get(), N_E, E_edges, w0Scan, coreGaussianSigma, "fit", histOutDir);

  // Then draw them side-by-side in a single 1×2 canvas
  drawAshLogSideBySide(ashFit, logFit, "FIT", E_edges, N_E, baseDir);

  // 4) (RMS) doAshScan & doLogScan => produce (AshRMS, LogRMS)
  std::cout<<"\n=== Doing RMS-based scanning (h->GetRMS) ===\n";
  auto ashRMS = doAshScan(f.get(), N_E, E_edges, bScan, cellSize, rawRMSlambda, "rms", histOutDir);
  auto logRMS = doLogScan(f.get(), N_E, E_edges, w0Scan, rawRMSlambda, "rms", histOutDir);

  // Then draw them side-by-side in a single 1×2 canvas
  drawAshLogSideBySide(ashRMS, logRMS, "RMS", E_edges, N_E, baseDir);

  std::cout<<"\n[DONE] plotAshLogRMS_sideBySide completed.\n";
}


////////////////////////////////////////////////////////////////////////////////
// Asinh-based function for local phi fits:
//
//   Y(X) = Normalization * [ 2*b / sqrt(1 + 4*X^2 * sinh^2(1/(2*b))) ]
////////////////////////////////////////////////////////////////////////////////
double asinhModel(double *x, double *par)
{
  double Norm = par[0];  // overall scale
  double b    = par[1];  // the "b" parameter
  if (b <= 1e-9) return 1e-12; // avoid division by zero if b is small

  double Xcg = x[0];
  double arg = 1.0 + 4.0*Xcg*Xcg*TMath::SinH(1.0/(2.0*b))*TMath::SinH(1.0/(2.0*b));
  double denom = TMath::Sqrt(arg);

  double numer = 2.0 * b;
  double result = Norm * (numer / denom);
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// doLocalPhiEtaFits(...)
//
//  1) Projects local-φ or local-η from h3 in slices of pT (here E-slices).
//  2) Overlays all E-bins on one canvas => LocalPhiFits.png or LocalEtaFits.png
//  3) If coord=phi, also does single‐bin overlays of raw vs corrected Δφ
//  4) Finally, we produce a canvas with each bin's local‐φ distribution & fit.
//
// **Changed**: now it's a 4×2 canvas for 8 bins (instead of 2×2).
////////////////////////////////////////////////////////////////////////////////
void doLocalPhiEtaFits(TH3F* h3,
                       const double pT_bins_low[/*N_E*/],
                       const double pT_bins_high[/*N_E*/],
                       std::ofstream& bResults,
                       const TString& outDirPhi)
{
  // We have N_E=8 bins
  // Ensure output directory exists
  gSystem->mkdir(outDirPhi, /*recursive=*/true);

  /////////////////////////////////////////////////////////////////////////////
  // Equation strings for φ vs η
  /////////////////////////////////////////////////////////////////////////////
  TString eqPhiString =
    "Y(X) = Norm #times #frac{2b_{#phi}}{#sqrt{1 + 4X^{2} sinh^{2}(1/(2b_{#phi}))}}";
  TString eqEtaString =
    "Y(X) = Norm #times #frac{2b_{#eta}}{#sqrt{1 + 4X^{2} sinh^{2}(1/(2b_{#eta}))}}";

  // For drawing the equation text
  TLatex eqLatex;
  eqLatex.SetNDC(true);
  eqLatex.SetTextFont(42);
  eqLatex.SetTextSize(0.028);
  eqLatex.SetTextAlign(13);

  // Style arrays for the bins
  int colors[N_E]   = { kBlack, kRed+1, kBlue+1, kGreen+2,
                        kMagenta+1, kOrange+2, kAzure+2, kSpring+9 };
  int markers[N_E]  = {20, 20, 20, 20, 21, 21, 21, 21};

  // -------------------------------------------------------------------------
  // Helper lambda for the quadruple overlay
  // -------------------------------------------------------------------------
  auto fitCoordAndMakePlot =
    [&](const char* coordName,    // "phi" or "eta"
        const char* projectOpt,   // "y" for φ, "x" for η
        const TString& outPNGname,
        const TString& xTitle)
  {
    bool isPhi = (strcmp(coordName, "phi") == 0);
    TString eqString = isPhi ? eqPhiString : eqEtaString;

    // If writing to bResults, label which coordinate we’re fitting
    if (bResults.is_open())
    {
      if (isPhi)
      {
        bResults << "\n# -- Now processing local phi fits --\n"
                 << "# E_low  E_high  b_phi\n";
      }
      else
      {
        bResults << "\n# -- Now processing local eta fits --\n"
                 << "# E_low  E_high  b_eta\n";
      }
    }

    // 1) Create a canvas for the “overlay” of all N_E=8 bins
    TCanvas cCoord(Form("cCoord_%s", coordName),
                   Form("Local %s distributions", coordName),
                   900, 700);
    cCoord.SetLeftMargin(0.12);
    cCoord.SetRightMargin(0.05);
    cCoord.SetBottomMargin(0.12);
    cCoord.cd();

    TLegend leg(0.18, 0.68, 0.55, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.028);

    // We'll store the 1D projections and best-fit b-values
    std::vector<TH1D*> hLocalVec(N_E, nullptr);
    std::vector<double> bestBvals(N_E, 0.0);

    // 2) Loop over the 8 E-slices => single overlay
    for(int i=0; i<N_E; i++)
    {
      // Determine Z range for that E slice
      int zLo = h3->GetZaxis()->FindFixBin(pT_bins_low[i] + 1e-9);
      int zHi = h3->GetZaxis()->FindFixBin(pT_bins_high[i] - 1e-9);
      if(zLo < 1) zLo=1;
      if(zHi> h3->GetNbinsZ()) zHi= h3->GetNbinsZ();

      // Set the range & do the 1D projection
      h3->GetXaxis()->SetRange(1, h3->GetNbinsX());
      h3->GetZaxis()->SetRange(zLo, zHi);

      TH1D* hLocal = (TH1D*) h3->Project3D(projectOpt);
      if(!hLocal)
      {
        std::cerr << "[WARNING] Could not project " << coordName
                  << " for E bin i=" << i << std::endl;
        continue;
      }
      hLocalVec[i] = hLocal;

      hLocal->SetName(Form("hLocal%s_%.1fto%.1f",
                           coordName, pT_bins_low[i], pT_bins_high[i]));
      hLocal->SetTitle(Form(";%s; scaled counts", xTitle.Data()));

      // Normalize
      double integral = hLocal->Integral();
      if(integral > 1e-9) hLocal->Scale(1. / integral);

      // Style
      hLocal->SetMarkerStyle(markers[i]);
      hLocal->SetMarkerColor(colors[i]);
      hLocal->SetLineColor  (colors[i]);
      hLocal->SetMarkerSize(1.2);

      // Draw: one big overlay with all bins
      if (i == 0)
      {
        hLocal->Draw("E");
        hLocal->GetYaxis()->SetRangeUser(0, 0.4);
        hLocal->GetXaxis()->SetTitle(xTitle);
        hLocal->GetYaxis()->SetTitle("counts (normalized)");
      }
      else
      {
        hLocal->Draw("SAME E");
      }

      // ---------------------------------------------------------------------
      // Multi‐start fit approach to find the best b-parameter
      // ---------------------------------------------------------------------
      TF1* fAsinh = new TF1(Form("fAsinh_%s_%d", coordName, i),
                            asinhModel, -0.5, 0.5, 2);
      fAsinh->SetParNames("Norm", "bVal");
      // Reasonable parameter limits
      fAsinh->SetParLimits(0, 1e-6, 1e6);
      fAsinh->SetParLimits(1, 1e-5, 1.0);

      // Some initial guesses
      std::vector<std::pair<double,double>> initialGuesses = {
        {0.2, 0.14},
        {0.1, 0.10},
        {0.3, 0.20},
        {0.5, 0.08}
      };

      double bestChi2 = 1e15;
      TF1* bestFunc   = nullptr;

      double fitMinX = -0.5;
      double fitMaxX =  0.5;

      // Try each guess
      for(auto& guess : initialGuesses)
      {
        fAsinh->SetParameter(0, guess.first);
        fAsinh->SetParameter(1, guess.second);

        TFitResultPtr fitRes = hLocal->Fit(fAsinh, "RQE0S", "",
                                           fitMinX, fitMaxX);
        if(fitRes.Get() && fitRes->IsValid())
        {
          double thisChi2 = fitRes->Chi2();
          if(thisChi2 < bestChi2)
          {
            bestChi2 = thisChi2;
            if(bestFunc) delete bestFunc;
            bestFunc = (TF1*) fAsinh->Clone(
              Form("bestFunc_%s_%d", coordName, i)
            );
          }
        }
      }

      double bVal = 0.0;
      if(bestFunc)
      {
        bestFunc->SetLineColor(colors[i]);
        bestFunc->SetLineWidth(2);
        bestFunc->Draw("SAME");

        bVal = bestFunc->GetParameter(1);
        // Attach the final best-func to the histogram so we can retrieve it later
        hLocal->GetListOfFunctions()->Add(bestFunc);
      }
      else
      {
        std::cerr << "[WARNING] No valid multi‐start fit for "
                  << coordName << ", E bin i=" << i
                  << " => bVal=0.\n";
      }
      bestBvals[i] = bVal;

      // Legend entry
      TString legEntry = Form("E=[%.1f,%.1f] GeV : b=%.3g",
                              pT_bins_low[i], pT_bins_high[i], bVal);
      leg.AddEntry(hLocal, legEntry, "lp");

      // Possibly record b-values
      if (isPhi && bResults.is_open() && bestFunc
          && (pT_bins_high[i] > pT_bins_low[i])
          && (bVal > 1e-9))
      {
        bResults << pT_bins_low[i]  << "  "
                 << pT_bins_high[i] << "   "
                 << bVal << "\n";
      }
    } // end loop over i

    // Add eq label and legend on the big overlay
    leg.Draw();
    eqLatex.DrawLatex(0.57, 0.82, eqString);

    // Save that “overlay”
    TString mainOut = outDirPhi + "/" + outPNGname;
    cCoord.SaveAs(mainOut);
    std::cout << "[INFO] Wrote overlay " << coordName
              << " => " << mainOut << std::endl;

    // -----------------------------------------------------------------------
    // If it's local φ, also produce:
    //   (a) single‐bin overlay of raw vs. corrected Δφ
    //   (b) a **4×2** layout for 8 bin local‐φ distributions + asinh fits
    // -----------------------------------------------------------------------
    if (isPhi)
    {
        //----------------------------------------------------------------------
        // (a) Single‐bin overlay: raw vs corrected, if those histograms exist
        //----------------------------------------------------------------------
        for(int i=0; i<N_E; i++)
        {
            double ptLo = pT_bins_low[i];
            double ptHi = pT_bins_high[i];
            
            // 1) Retrieve the RAW Δφ histogram
            TString rawName = Form("h_phi_diff_raw_%.0f_%.0f", ptLo, ptHi);
            TH1F* hRaw = dynamic_cast<TH1F*>(gROOT->FindObject(rawName));
            if(!hRaw)
            {
                std::cerr << "[WARNING] No raw Δphi hist '" << rawName
                          << "' => skip overlay for E bin i=" << i << std::endl;
                continue;
            }
            
            // 2) Retrieve the CORRECTED Δφ histogram
            TString corrName = Form("h_phi_diff_corr_%.1f_%.1f", ptLo, ptHi);
            TH1F* hCorr = dynamic_cast<TH1F*>(gROOT->FindObject(corrName));
            if(!hCorr)
            {
                std::cerr << "[WARNING] No corrected Δphi hist '" << corrName
                          << "' => skip overlay for E bin i=" << i << std::endl;
                continue;
            }
            
            // 3) Make a canvas & style each histogram
            TString cName  = Form("cDeltaPhiCompare_%.1fto%.1f", ptLo, ptHi);
            TString cTitle = Form("#Delta#phi overlay [%.1f < E < %.1f]",
                                  ptLo, ptHi);
            
            TCanvas cSingle(cName, cTitle, 800, 600);
            cSingle.SetLeftMargin(0.12);
            cSingle.SetRightMargin(0.05);
            cSingle.SetBottomMargin(0.12);
            cSingle.cd();
            
            // Normalize each to area=1
            double rInt = hRaw->Integral();
            double cInt = hCorr->Integral();
            if(rInt  > 1e-9) hRaw ->Scale(1./rInt);
            if(cInt  > 1e-9) hCorr->Scale(1./cInt);
            
            // Style: raw in black
            hRaw->SetLineColor(kBlack);
            hRaw->SetMarkerColor(kBlack);
            hRaw->SetMarkerStyle(20);
            hRaw->SetMarkerSize(1);
            hRaw->SetTitle(cTitle);
            hRaw->GetXaxis()->SetTitle("#Delta#phi (reco - truth) [radians]");
            hRaw->GetYaxis()->SetTitle("counts (normalized)");
            hRaw->GetYaxis()->SetTitleOffset(1.2);
            
            // corrected in red
            hCorr->SetLineColor(kRed);
            hCorr->SetMarkerColor(kRed);
            hCorr->SetMarkerStyle(21);
            hCorr->SetMarkerSize(1);
            
            // 4) Draw them
            hRaw->Draw("E");
            double maxRaw  = hRaw->GetMaximum();
            double maxCorr = hCorr->GetMaximum();
            double newMax  = TMath::Max(maxRaw, maxCorr)*1.25;
            hRaw->GetYaxis()->SetRangeUser(0, newMax);
            
            hCorr->Draw("SAME E");
            
            // 5) Add a legend
            TLegend legSingle(0.55, 0.70, 0.85, 0.85);
            legSingle.SetBorderSize(0);
            legSingle.SetFillStyle(0);
            legSingle.SetTextSize(0.035);
            
            legSingle.AddEntry(hRaw,
                               Form("Raw #Delta#phi (%.1f< E<%.1f)", ptLo, ptHi),
                               "lp");
            legSingle.AddEntry(hCorr,
                               Form("Corrected #Delta#phi (%.1f< E<%.1f)", ptLo, ptHi),
                               "lp");
            legSingle.Draw();
            
            // 6) Save
            TString singleOut = Form("%s/DeltaPhiCompare_%.1fto%.1f.png",
                                     outDirPhi.Data(), ptLo, ptHi);
            cSingle.SaveAs(singleOut);
            
            std::cout << "[INFO] Single-bin Δphi overlay E=["
            << ptLo << "," << ptHi
            << "] => " << singleOut << std::endl;
        } // end loop over i
        
        //----------------------------------------------------------------------
        // (b) A **4×2** canvas showing each bin's local‐φ distribution & best-fit,
        //     possibly overlaying the corrected local‐φ projection
        //----------------------------------------------------------------------
        TCanvas c2by2("cLocalPhi4by2", "Local φ fits per E bin", 1600, 1000);
        c2by2.Divide(4, 2);  // 4 columns, 2 rows = 8 pads

        // Attempt to retrieve the corrected TH3
        TH3F* h3corr = dynamic_cast<TH3F*>(gROOT->FindObject("h2_cluster_block_cord_E_corrected"));
        if(!h3corr)
        {
            std::cout << "[WARNING] 'h2_cluster_block_cord_E_corrected' not found; no overlay.\n";
        }
        
        for (int i=0; i<N_E; i++)
        {
            c2by2.cd(i+1);
            // from the big overlay
            TH1D* histLocal = hLocalVec[i];
            if (!histLocal) continue;  // skip if empty
            
            // Put E range + b param in the title
            double bVal = 0.0;
            TF1* bf = dynamic_cast<TF1*>(
                       histLocal->GetListOfFunctions()->FindObject(
                         Form("bestFunc_phi_%d", i)
                       )
                     );
            if(bf) bVal = bf->GetParameter(1);

            TString binTitle = Form("Local #phi: E=[%.1f,%.1f], b=%.3g",
                                    pT_bins_low[i], pT_bins_high[i], bVal);
            histLocal->SetTitle(binTitle);
            histLocal->GetXaxis()->SetTitle("local #phi_{CG}");
            histLocal->GetYaxis()->SetTitle("counts (normalized)");
            histLocal->Draw("E");
            
            // Overlay the best function
            if(bf) bf->Draw("SAME");
            
            //-------------------------------------------------------------
            // Overlap the corrected local-φ from h3corr, if present
            //-------------------------------------------------------------
            if (h3corr)
            {
                // Define the pT bin range
                int zLoCorr = h3corr->GetZaxis()->FindFixBin(pT_bins_low[i] + 1e-9);
                int zHiCorr = h3corr->GetZaxis()->FindFixBin(pT_bins_high[i] - 1e-9);
                if(zLoCorr < 1) zLoCorr=1;
                if(zHiCorr> h3corr->GetNbinsZ()) zHiCorr= h3corr->GetNbinsZ();
                
                // Project in the φ dimension => "y"
                h3corr->GetXaxis()->SetRange(1, h3corr->GetNbinsX());
                h3corr->GetZaxis()->SetRange(zLoCorr, zHiCorr);
                
                TH1D* hCorrLocal = (TH1D*) h3corr->Project3D("y");
                if (hCorrLocal)
                {
                    double integC = hCorrLocal->Integral();
                    if (integC > 1e-9) hCorrLocal->Scale(1./integC);
                    
                    hCorrLocal->SetLineColor(kMagenta+1);
                    hCorrLocal->SetMarkerColor(kMagenta+1);
                    hCorrLocal->SetMarkerStyle(25);
                    hCorrLocal->SetMarkerSize(1.0);
                    hCorrLocal->SetTitle("");
                    
                    hCorrLocal->Draw("SAME E");
                    
                    // Possibly add a small local legend
                    TLegend legCorr(0.50, 0.65, 0.88, 0.82);
                    legCorr.SetBorderSize(0);
                    legCorr.SetFillStyle(0);
                    legCorr.SetTextSize(0.035);
                    legCorr.AddEntry(histLocal, "Uncorrected", "lp");
                    legCorr.AddEntry(hCorrLocal, "Corrected", "lp");
                    legCorr.Draw("SAME");
                }
                else
                {
                    std::cout << "[WARNING] Could not project corrected local-φ for bin " << i << std::endl;
                }
            } // end if(h3corr)
        }
        
        // Save the 4×2 layout
        TString out2by2 = outDirPhi + "/LocalPhiFits_4by2.png";
        c2by2.SaveAs(out2by2);
        std::cout << "[INFO] Wrote 4×2 local φ bin-by-bin fits => "
        << out2by2 << std::endl;
    }

  }; // end fitCoordAndMakePlot


  // -----------------------------------------------------------------------
  // (1) call for φ => “LocalPhiFits.png”
  // -----------------------------------------------------------------------
  fitCoordAndMakePlot(
    "phi",
    "y",
    "LocalPhiFits.png",
    "local #phi_{CG} in block"
  );

  // -----------------------------------------------------------------------
  // (2) call for η => “LocalEtaFits.png”
  // -----------------------------------------------------------------------
  fitCoordAndMakePlot(
    "eta",
    "x",
    "LocalEtaFits.png",
    "local #eta_{CG} in block"
  );

}


////////////////////////////////////////////////////////////////////////////////
// PDCanalysis()
//
//  1) Lists all histograms (like the original code).
//  2) Makes a 2D block-eta vs block-phi plot from h2_cluster_block_cord_E
//  3) Makes local-phi distributions for E slices [2,4], [4,6], [6,8], [8,10],
//     [10,12], [12,15], [15,20], [20,30], fits them, and records b-values.
//  4) Finally, saves *every* histogram as a PNG in:
//       /Users/.../SimOutput/allHistOutput
////////////////////////////////////////////////////////////////////////////////
void PDCanalysis()
{
  ///////////////////////////////////////////////////////////////////////
  // 0) Basic setup
  ///////////////////////////////////////////////////////////////////////
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
    
  // Root file with merged histograms
  const char* filename = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root";
  TFile* f = TFile::Open(filename,"READ");
  if(!f || f->IsZombie())
  {
    std::cerr << "[ERROR] Could not open file: " << filename << std::endl;
    return;
  }
  std::cout << "[INFO] Successfully opened: " << filename << std::endl;

  ///////////////////////////////////////////////////////////////////////
  // 1) Print a table of all histograms (the original function)
  ///////////////////////////////////////////////////////////////////////
  std::cout << std::left
            << std::setw(30) << "HISTOGRAM NAME"
            << std::setw(20) << "CLASS TYPE"
            << std::setw(15) << "ENTRIES"
            << std::endl;
  std::cout << std::string(65, '=') << std::endl;

  // We'll store the list of keys for reuse
  TList* keyList = f->GetListOfKeys();
  if(!keyList || keyList->IsEmpty())
  {
    std::cerr << "[WARNING] The file has no keys/histograms.\n";
  }
  else
  {
    TIter nextkey(keyList);
    TKey* key;
    while((key = (TKey*) nextkey()))
    {
      TClass* cl = gROOT->GetClass(key->GetClassName());
      if(!cl) continue;
      if(cl->InheritsFrom("TH1"))
      {
        TH1* hist = (TH1*) key->ReadObj();
        if(!hist) continue;
        std::cout << std::setw(30) << hist->GetName()
                  << std::setw(20) << hist->ClassName()
                  << std::setw(15) << (long long)hist->GetEntries()
                  << std::endl;
      }
    }
  }
  std::cout << "[INFO] Finished listing histograms.\n" << std::endl;

  ///////////////////////////////////////////////////////////////////////
  // 2) Grab the 3D histogram for block coords vs E
  ///////////////////////////////////////////////////////////////////////
  TH3F* h3 = (TH3F*) f->Get("h2_cluster_block_cord_E");
    
  if(!h3)
  {
    std::cerr << "[ERROR] 'h2_cluster_block_cord_E' not found => can't do block-eta vs. block-phi.\n";
    f->Close(); delete f;
    return;
  }

  // Output directories
  const TString baseDir   = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
  const TString outDir2D  = baseDir + "/2DPlots";
  const TString outDirPhi = baseDir + "/LocalPhiFits";
  const TString outDirAll = baseDir + "/allHistOutput";

  // Make sure directories exist
  gSystem->mkdir(baseDir,   true);
  gSystem->mkdir(outDir2D,  true);
  gSystem->mkdir(outDirPhi, true);
  gSystem->mkdir(outDirAll, true);

  ///////////////////////////////////////////////////////////////////////
  // 3) Create 2D plot from h3 by integrating over entire E
  ///////////////////////////////////////////////////////////////////////
    {
      // First, restore the Z-axis range to include all E
      h3->GetZaxis()->SetRange(1, h3->GetNbinsZ());
      TH2F* h2_block2D = (TH2F*) h3->Project3D("xy");
      h2_block2D->SetName("h2_blockCoord2D");
      h2_block2D->SetTitle("Block #eta vs Block #phi;block #eta;block #phi");

      TCanvas c2d("c2d","Block Eta-Phi (integrated)",900,700);
      c2d.SetLeftMargin(0.12);
      c2d.SetRightMargin(0.15);
      c2d.SetBottomMargin(0.12);

      h2_block2D->SetStats(false);
      h2_block2D->Draw("COLZ");

      c2d.SaveAs(outDir2D + "/BlockEtaPhi_2D.png");
      std::cout << "[INFO] Wrote 2D block-eta vs block-phi (all E) => "
                << (outDir2D + "/BlockEtaPhi_2D.png") << std::endl;
    }

    ///////////////////////////////////////////////////////////////////
    // 3b) Also produce 2D plots for each E slice
    ///////////////////////////////////////////////////////////////////
    {
      for(int i=0; i<N_E; i++)
      {
        // Determine Z-range in the 3D histogram (E axis)
        int zLo = h3->GetZaxis()->FindFixBin( E_edges[i] + 1e-9 );
        int zHi = h3->GetZaxis()->FindFixBin( E_edges[i+1] - 1e-9 );
        if(zLo < 1) zLo = 1;
        if(zHi > h3->GetNbinsZ()) zHi = h3->GetNbinsZ();

        // Apply that Z-range, then project onto X-Y => 2D
        h3->GetZaxis()->SetRange(zLo, zHi);
        TH2F* h2_slice = (TH2F*) h3->Project3D("xy");

        // Name & title
        h2_slice->SetName( Form("h2_blockCoord2D_E%.0fto%.0f", E_edges[i], E_edges[i+1]) );
        h2_slice->SetTitle( Form("Block #eta vs Block #phi (%.1f < E < %.1f);block #eta;block #phi",
                                   E_edges[i], E_edges[i+1]) );

        // Draw
        TCanvas c2dSlice( Form("c2d_E%.0fto%.0f", E_edges[i], E_edges[i+1]),
                          "Block Eta-Phi (E slice)", 900, 700);
        c2dSlice.SetLeftMargin(0.12);
        c2dSlice.SetRightMargin(0.15);
        c2dSlice.SetBottomMargin(0.12);

        h2_slice->SetStats(false);
        h2_slice->Draw("COLZ");

        // Save a PNG for each slice
        TString outName = Form("/BlockEtaPhi_2D_E%.0fto%.0f.png",
                               E_edges[i], E_edges[i+1]);
        c2dSlice.SaveAs(outDir2D + outName);

        std::cout << "[INFO] Wrote 2D block-eta vs block-phi => "
                  << (outDir2D + outName)
                  << " for E=[" << E_edges[i] << ", " << E_edges[i+1] << "]" << std::endl;
      }
    }

  ///////////////////////////////////////////////////////////////////////
  // 4) Local φ slices in E bins => do 1D projections, fits
  ///////////////////////////////////////////////////////////////////////

  // We'll record b-values to a text file
  std::ofstream bResults(Form("%s/bValues.txt", baseDir.Data()));
  if(!bResults.is_open()) {
    std::cerr << "[WARNING] Could not open bValues.txt for writing.\n";
  } else {
    bResults << "# E_low  E_high   bValue\n";
  }

  double pT_bins_low[N_E], pT_bins_high[N_E];
  for(int j=0; j<N_E; j++){
      pT_bins_low[j]  = E_edges[j];
      pT_bins_high[j] = E_edges[j+1];
  }
  doLocalPhiEtaFits(h3, pT_bins_low, pT_bins_high, bResults, outDirPhi);
    
  std::cout << "[INFO] Wrote local phi slice overlay => "
            << (outDirPhi + "/LocalPhiFits.png") << std::endl;

  if(bResults.is_open()) {
    bResults.close();
    std::cout << "[INFO] Wrote b-values to " << (baseDir + "/bValues.txt") << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////
  // 5) Save *every* histogram in the file as a PNG in /allHistOutput/
  ///////////////////////////////////////////////////////////////////////
  std::cout << "\n[INFO] Now saving *every* histogram as a PNG => " << outDirAll << std::endl;

  TIter nextAll(f->GetListOfKeys());
  TKey* keyAll;
  while( (keyAll = (TKey*) nextAll()) )
  {
    TClass* cl = gROOT->GetClass(keyAll->GetClassName());
    if(!cl) continue;

    TObject* obj = keyAll->ReadObj();
    if(!obj) continue;

    // We'll skip non-histogram objects
    if(!cl->InheritsFrom("TH1") && !cl->InheritsFrom("TProfile")) continue;

    TH1* htmp = dynamic_cast<TH1*>(obj);
    if(!htmp) continue;

    // Build output PNG name
    TString histName = htmp->GetName();
    TString outPNG   = outDirAll + "/" + histName + ".png";

    TCanvas ctemp("ctemp","",800,600);
    ctemp.SetLeftMargin(0.12);
    ctemp.SetRightMargin(0.15);
    ctemp.SetBottomMargin(0.12);

    // Decide how to draw
    if( htmp->InheritsFrom("TH2") ) {
      htmp->Draw("COLZ");
    }
    else if(htmp->InheritsFrom("TH3")) {
      // 3D is tricky to visualize. We'll do "BOX".
      htmp->Draw("BOX");
    }
    else if(htmp->InheritsFrom("TProfile")) {
      htmp->Draw("E1");
    }
    else {
      // presumably TH1
      htmp->Draw("E");
    }

    ctemp.SaveAs(outPNG);
    delete htmp;
  }

  // Finally call the side‐by‐side RMS vs b/w0
  plotAshLogRMS_sideBySide(filename);

  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
