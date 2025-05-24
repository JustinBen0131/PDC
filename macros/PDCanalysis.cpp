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

constexpr double pTedges[] = {0.5, 0.75, 1.25, 2, 3};
// Let the compiler figure out how many bins that is:
constexpr int N_PT = (sizeof(pTedges)/sizeof(pTedges[0])) - 1;

void plotAshLogRMS_sideBySide(const char* infile = "PositionDep_sim_ALL.root")
{
  // --------------------------------------------------------------------------------
  // 1) Parameter setup
  // --------------------------------------------------------------------------------
  // pT bins come from pTedges[] above.

  // Ash (b) scan
  const double bMin   = 0.50;
  const double bMax   = 2.00;
  const double bStep  = 0.05;

  std::vector<double> bScan;
  for(double b = bMin; b <= bMax + 1e-9; b += bStep) bScan.push_back(b);
  const int N_B = bScan.size();

  // Log (w0) scan
  const double w0Min  = 2.00;
  const double w0Max  = 6.00;
  const double w0Step = 0.10;

  std::vector<double> w0Scan;
  for(double w = w0Min; w <= w0Max + 1e-9; w += w0Step) w0Scan.push_back(w);
  const int N_W = w0Scan.size();

  // Output location
  const TString baseDir    = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
  const TString histOutDir = baseDir + "/ASH_LOG_PLOTS";

  gSystem->mkdir(baseDir,    /*parents=*/true);
  gSystem->mkdir(histOutDir, /*parents=*/true);

  // --------------------------------------------------------------------------------
  // 2) Input file
  // --------------------------------------------------------------------------------
  std::cout << "[INFO] Opening file '" << infile << "'..." << std::endl;
  std::unique_ptr<TFile> f(TFile::Open(infile,"READ"));
  if(!f || f->IsZombie()){
    std::cerr << "[ERROR] Cannot open file or file is Zombie.\n";
    return;
  }
  std::cout << "[INFO] Successfully opened '" << infile << "'." << std::endl;

  // --------------------------------------------------------------------------------
  // 3) ROOT style
  // --------------------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);

  // Setup color + marker for each pT bin
  std::vector<int> colors;
  std::vector<int> markers;
  {
    // Example color set (at least 11 distinct choices)
    const int baseColorList[] = {
      kBlack, kRed, kBlue, kMagenta+2, kOrange,
      kGreen+2, kCyan+2, kGray+1, kViolet+1, kAzure+2, kSpring+9
    };
    // Marker styles
    const int baseMarkerList[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

    for (int i = 0; i < N_PT; ++i) {
      colors.push_back( baseColorList[i % 11] );
      markers.push_back(baseMarkerList[i % 11]);
    }
  }

  // --------------------------------------------------------------------------------
  // 4) Ash (b) scan:  σx vs b
  // --------------------------------------------------------------------------------
  std::cout << "\n[INFO] ASH scan: We'll loop over " << N_B << " b-values "
            << "(" << bMin << " to " << bMax << " step " << bStep << ") "
            << " for " << N_PT << " pT bins."
            << "\n       => Expect " << (N_B * N_PT) << " histograms.\n";

  TCanvas cAsh("cAsh","Ash RMS vs b",800,600);
  cAsh.SetLeftMargin(0.12);
  cAsh.SetBottomMargin(0.12);
  cAsh.SetGrid();

  double gMinAsh =  DBL_MAX;
  double gMaxAsh = -DBL_MAX;

  std::vector<TGraph*> gAshVec;
  gAshVec.reserve(N_PT);

  // For each pT bin, track best (lowest) RMS
  std::vector<double> bestB(N_PT, 0.);
  std::vector<double> bestRMS_A(N_PT, DBL_MAX);

  // Counters for how many histograms are missing/found
  int totalHistAsh    = 0;
  int missingHistAsh  = 0;
  int zeroEntryAsh    = 0;  // found hist but 0 entries
  int usedHistAsh     = 0;  // found and non-empty

  for(int ipt = 0; ipt < N_PT; ++ipt)
  {
    auto g = new TGraph;
    g->SetName(Form("gAsh_pt%d", ipt));

    std::cout << "\n  [ASH-scan: pT bin " << ipt << "] Searching for histograms...\n";

    for(size_t ib=0; ib < bScan.size(); ++ib)
    {
      double bVal = bScan[ib];
      ++totalHistAsh; // We expect one histogram for each (ipt, bVal)

      TString hN  = Form("h_dx_ash_b%.3f_pt%d", bVal, ipt);
      TH1* h = dynamic_cast<TH1*>(f->Get(hN));

      if(!h) {
        std::cerr << "   [WARN] Missing histogram: " << hN << std::endl;
        ++missingHistAsh;
        continue;
      }
      // Found the histogram
      double entries = h->GetEntries();
      double integral= h->Integral();
      double rms     = h->GetRMS();

      if(entries == 0) {
        std::cout << "   [WARN] Found " << hN << " but it has 0 entries.\n";
        ++zeroEntryAsh;
      } else {
        ++usedHistAsh;
        std::cout << "   [INFO] Found " << hN << " => entries=" << entries
                  << ", integral=" << integral << ", RMS=" << rms << std::endl;
      }

      // Store point even if zero entries => might be 0 RMS
      g->SetPoint(g->GetN(), bVal, rms);

      if(rms > 0){
        gMinAsh = std::min(gMinAsh, rms);
        gMaxAsh = std::max(gMaxAsh, rms);
        if(rms < bestRMS_A[ipt]) {
          bestRMS_A[ipt] = rms;
          bestB[ipt]     = bVal;
        }
      }

      // Quick PNG snapshot
      // (Optional: disable if you don't need each histogram in a separate PNG)
      TCanvas ctmp;
      h->Draw("E");
      ctmp.SaveAs(histOutDir + "/" + hN + ".png");
    }

    g->SetMarkerStyle(markers[ipt]);
    g->SetMarkerColor(colors [ipt]);
    g->SetLineColor  (colors [ipt]);
    g->SetMarkerSize(1.3);
    gAshVec.push_back(g);
  }

  // If we never updated gMinAsh, it means no hist found at all
  if(gMinAsh == DBL_MAX) { gMinAsh=0; gMaxAsh=1; }

  cAsh.cd();
  if(!gAshVec.empty()) {
    gAshVec[0]->Draw("ALP");
    gAshVec[0]->GetXaxis()->SetTitle("b");
    gAshVec[0]->GetYaxis()->SetTitle("#sigma_{x}");
    gAshVec[0]->GetYaxis()->SetRangeUser(gMinAsh-0.1*fabs(gMinAsh),
                                         gMaxAsh+0.1*fabs(gMaxAsh));
    for(size_t i=1; i < gAshVec.size(); ++i) gAshVec[i]->Draw("LP SAME");
  }

  TLegend legA(0.15,0.70,0.45,0.85);
  legA.SetBorderSize(0);
  legA.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt) {
    legA.AddEntry(gAshVec[ipt],
      Form("%.1f < p_{T} < %.1f GeV  (best b=%.2f, RMS=%.4f)",
           pTedges[ipt], pTedges[ipt+1], bestB[ipt], bestRMS_A[ipt]),
      "lp");
  }
  legA.Draw();
  TLatex().DrawLatexNDC(0.2,0.92,"Ash");

  cAsh.SaveAs(baseDir+"/Ash_RMS_vs_b.png");
  std::cout << "\n[INFO] Saved " << (baseDir + "/Ash_RMS_vs_b.png") << std::endl;

  // Print summary for Ash
  std::cout << "\n[ASH SUMMARY]-----------------------------------------\n"
            << "  # pT bins             = " << N_PT << "\n"
            << "  # b-values scanned    = " << N_B << "\n"
            << "  => total histograms   = " << totalHistAsh << "\n"
            << "  => missing histograms = " << missingHistAsh << "\n"
            << "  => zero-entry hists   = " << zeroEntryAsh << "\n"
            << "  => used (non-empty)   = " << usedHistAsh << "\n";
  for(int ipt=0; ipt<N_PT; ++ipt){
    std::cout << "   pT bin " << ipt << " => best b=" << bestB[ipt]
              << ", best RMS=" << bestRMS_A[ipt] << "\n";
  }
  std::cout << "------------------------------------------------------\n";

  // --------------------------------------------------------------------------------
  // 5) Log (w0) scan:  σx vs w0
  // --------------------------------------------------------------------------------
  std::cout << "\n[INFO] LOG scan: We'll loop over " << N_W << " w0-values "
            << "(" << w0Min << " to " << w0Max << " step " << w0Step << ") "
            << " for " << N_PT << " pT bins."
            << "\n       => Expect " << (N_W * N_PT) << " histograms.\n";

  TCanvas cLog("cLog","Log RMS vs w0",800,600);
  cLog.SetLeftMargin(0.12);
  cLog.SetBottomMargin(0.12);
  cLog.SetGrid();

  double gMinLog =  DBL_MAX;
  double gMaxLog = -DBL_MAX;

  std::vector<TGraph*> gLogVec;
  gLogVec.reserve(N_PT);

  // For each pT bin, track best (lowest) RMS for the log scan
  std::vector<double> bestW(N_PT, 0.);
  std::vector<double> bestRMS_L(N_PT, DBL_MAX);

  // Counters for log hist usage
  int totalHistLog    = 0;
  int missingHistLog  = 0;
  int zeroEntryLog    = 0;
  int usedHistLog     = 0;

  for(int ipt=0; ipt<N_PT; ++ipt)
  {
    auto g = new TGraph;
    g->SetName(Form("gLog_pt%d", ipt));

    std::cout << "\n  [LOG-scan: pT bin " << ipt << "] Searching for histograms...\n";

    for(size_t iw=0; iw < w0Scan.size(); ++iw)
    {
      double wVal = w0Scan[iw];
      ++totalHistLog;

      TString hN  = Form("h_dx_log_w0%.2f_pt%d", wVal, ipt);
      TH1* h = dynamic_cast<TH1*>(f->Get(hN));
      if(!h) {
        std::cerr << "   [WARN] Missing histogram: " << hN << std::endl;
        ++missingHistLog;
        continue;
      }

      double entries = h->GetEntries();
      double integral= h->Integral();
      double rms     = h->GetRMS();

      if(entries == 0) {
        std::cout << "   [WARN] Found " << hN << " but 0 entries.\n";
        ++zeroEntryLog;
      } else {
        ++usedHistLog;
        std::cout << "   [INFO] Found " << hN << " => entries=" << entries
                  << ", integral=" << integral << ", RMS=" << rms << std::endl;
      }

      g->SetPoint(g->GetN(), wVal, rms);

      if(rms > 0){
        gMinLog = std::min(gMinLog, rms);
        gMaxLog = std::max(gMaxLog, rms);
        if(rms < bestRMS_L[ipt]) {
          bestRMS_L[ipt] = rms;
          bestW[ipt]     = wVal;
        }
      }

      TCanvas ctmp;
      h->Draw("E");
      ctmp.SaveAs(histOutDir + "/" + hN + ".png");
    }

    g->SetMarkerStyle(markers[ipt]);
    g->SetMarkerColor(colors[ipt]);
    g->SetLineColor  (colors[ipt]);
    g->SetMarkerSize(1.3);
    gLogVec.push_back(g);
  }

  if(gMinLog == DBL_MAX) { gMinLog=0; gMaxLog=1; }

  cLog.cd();
  if(!gLogVec.empty()) {
    gLogVec[0]->Draw("ALP");
    gLogVec[0]->GetXaxis()->SetTitle("w_{0}");
    gLogVec[0]->GetYaxis()->SetTitle("#sigma_{x}");
    gLogVec[0]->GetYaxis()->SetRangeUser(gMinLog-0.1*fabs(gMinLog),
                                         gMaxLog+0.1*fabs(gMaxLog));
    for(size_t i=1; i < gLogVec.size(); ++i) gLogVec[i]->Draw("LP SAME");
  }

  TLegend legB(0.15,0.65,0.45,0.85);
  legB.SetBorderSize(0);
  legB.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt) {
    legB.AddEntry(gLogVec[ipt],
      Form("%.1f < p_{T} < %.1f GeV (best w_{0}=%.2f, RMS=%.4f)",
           pTedges[ipt], pTedges[ipt+1], bestW[ipt], bestRMS_L[ipt]),
      "lp");
  }
  legB.Draw();
  TLatex().DrawLatexNDC(0.2,0.92,"Log");

  cLog.SaveAs(baseDir+"/Log_RMS_vs_w0.png");
  std::cout << "\n[INFO] Saved " << (baseDir + "/Log_RMS_vs_w0.png") << std::endl;

  // Print summary for Log
  std::cout << "\n[LOG SUMMARY]------------------------------------------\n"
            << "  # pT bins             = " << N_PT << "\n"
            << "  # w0-values scanned   = " << N_W << "\n"
            << "  => total histograms   = " << totalHistLog << "\n"
            << "  => missing histograms = " << missingHistLog << "\n"
            << "  => zero-entry hists   = " << zeroEntryLog << "\n"
            << "  => used (non-empty)   = " << usedHistLog << "\n";
  for(int ipt=0; ipt<N_PT; ++ipt){
    std::cout << "   pT bin " << ipt << " => best w0=" << bestW[ipt]
              << ", best RMS=" << bestRMS_L[ipt] << "\n";
  }
  std::cout << "--------------------------------------------------------\n";

  std::cout << "\n[DONE] plotAshLogRMS_sideBySide finished.\n" << std::endl;
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
//  1) Projects local-φ or local-η from h3 in slices of pT
//  2) Overlays all 4 bins on one canvas => LocalPhiFits.png or LocalEtaFits.png
//  3) If coord=phi, also does single‐bin overlays of raw vs corrected Δφ
//  4) **Now** also produces a 2×2 table of each per‐bin local‐φ distribution
//     with its best-fit, including pT range and b-parameter in the title.
////////////////////////////////////////////////////////////////////////////////
void doLocalPhiEtaFits(TH3F* h3,
                       const double pT_bins_low[4],
                       const double pT_bins_high[4],
                       std::ofstream& bResults,
                       const TString& outDirPhi)
{
  // We have exactly 4 pT bins:
  const int N_PT = 4;

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

  // Style arrays for the 4 bins
  int colors[N_PT]   = { kBlack, kRed+1, kBlue+1, kGreen+2 };
  int markers[N_PT]  = {20, 20, 20, 20}; // all circles

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
                 << "# pTlow  pThigh  b_phi\n";
      }
      else
      {
        bResults << "\n# -- Now processing local eta fits --\n"
                 << "# pTlow  pThigh  b_eta\n";
      }
    }

    // 1) Create a canvas for the “quadruple overlay” of the 4 pT bins
    TCanvas cCoord(Form("cCoord_%s", coordName),
                   Form("Local %s distributions", coordName),
                   900, 700);
    cCoord.SetLeftMargin(0.12);
    cCoord.SetRightMargin(0.05);
    cCoord.SetBottomMargin(0.12);
    cCoord.cd();

    TLegend leg(0.18, 0.68, 0.4, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.028);

    // We'll store the 1D projections and best-fit b-values
    std::vector<TH1D*> hLocalVec(N_PT, nullptr);
    std::vector<double> bestBvals(N_PT, 0.0);

    // 2) Loop over the 4 pT slices => single overlay
    for(int i=0; i<N_PT; i++)
    {
      // Determine Z range for that pT slice
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
                  << " for pT bin i=" << i << std::endl;
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
                  << coordName << ", pT bin i=" << i
                  << " => bVal=0.\n";
      }
      bestBvals[i] = bVal;

      // Legend entry
      TString legEntry = Form("p_{T}=[%.1f,%.1f] GeV : b=%.3g",
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

    // Save that “quadruple overlay”
    TString mainOut = outDirPhi + "/" + outPNGname;
    cCoord.SaveAs(mainOut);
    std::cout << "[INFO] Wrote quadruple-overlay " << coordName
              << " => " << mainOut << std::endl;

    // -----------------------------------------------------------------------
    // If it's local φ, also produce:
    //   (a) single‐bin overlay of raw vs. corrected Δφ
    //   (b) a 2×2 table of each bin’s local‐φ distribution & asinh fit
    // -----------------------------------------------------------------------
    if (isPhi)
    {
      //----------------------------------------------------------------------
      // (a) Single‐bin overlay: raw vs corrected, if those histograms exist
      //----------------------------------------------------------------------
      for(int i=0; i<N_PT; i++)
      {
        double ptLo = pT_bins_low[i];
        double ptHi = pT_bins_high[i];

        // 1) Retrieve the RAW Δφ histogram
        TString rawName = Form("h_phi_diff_raw_%.0f_%.0f", ptLo, ptHi);
        TH1F* hRaw = dynamic_cast<TH1F*>(gROOT->FindObject(rawName));
        if(!hRaw)
        {
          std::cerr << "[WARNING] No raw Δφ hist '" << rawName
                    << "' => skip overlay for pT bin i=" << i << std::endl;
          continue;
        }

        // 2) Retrieve the CORRECTED Δφ histogram
        TString corrName = Form("h_phi_diff_corr_%.1f_%.1f", ptLo, ptHi);
        TH1F* hCorr = dynamic_cast<TH1F*>(gROOT->FindObject(corrName));
        if(!hCorr)
        {
          std::cerr << "[WARNING] No corrected Δφ hist '" << corrName
                    << "' => skip overlay for pT bin i=" << i << std::endl;
          continue;
        }

        // 3) Make a canvas & style each histogram
        TString cName  = Form("cDeltaPhiCompare_%.1fto%.1f", ptLo, ptHi);
        TString cTitle = Form("#Delta#phi overlay [%.1f < p_{T} < %.1f]",
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
            Form("Raw #Delta#phi (%.1f< p_{T}<%.1f)", ptLo, ptHi),
            "lp");
        legSingle.AddEntry(hCorr,
            Form("Corrected #Delta#phi (%.1f< p_{T}<%.1f)", ptLo, ptHi),
            "lp");
        legSingle.Draw();

        // 6) Save
        TString singleOut = Form("%s/DeltaPhiCompare_%.1fto%.1f.png",
                                 outDirPhi.Data(), ptLo, ptHi);
        cSingle.SaveAs(singleOut);

        std::cout << "[INFO] Single-bin Δphi overlay pT=["
                  << ptLo << "," << ptHi
                  << "] => " << singleOut << std::endl;
      } // end loop over i

      //----------------------------------------------------------------------
      // (b) A 2×2 canvas showing each bin's local‐φ distribution & best-fit
      //----------------------------------------------------------------------
      TCanvas c2by2("cLocalPhi2by2", "Local φ fits per pT bin", 1000, 800);
      c2by2.Divide(2, 2);

      for (int i=0; i<N_PT; i++)
      {
        c2by2.cd(i+1);
        TH1D* histLocal = hLocalVec[i];
        if (!histLocal) continue;  // skip if empty

        // Put pT range + b param in the title
        double bVal = bestBvals[i];
        TString binTitle = Form("Local #phi: p_{T}=[%.1f,%.1f], b=%.3g",
                                pT_bins_low[i], pT_bins_high[i], bVal);
        histLocal->SetTitle(binTitle);
        histLocal->GetXaxis()->SetTitle(xTitle);
        histLocal->GetYaxis()->SetTitle("counts (normalized)");
        histLocal->Draw("E");

        // Overlay the best fit if it exists
        TF1* bestF = dynamic_cast<TF1*>(
          histLocal->GetListOfFunctions()->FindObject(
            Form("bestFunc_%s_%d", coordName, i)
          )
        );
        if(bestF) bestF->Draw("SAME");
      }

      // Save the 2×2 layout
      TString out2by2 = outDirPhi + "/LocalPhiFits_2by2.png";
      c2by2.SaveAs(out2by2);
      std::cout << "[INFO] Wrote 2×2 local φ bin-by-bin fits => "
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
// overlayUncorrCorrLocalPhi_noFits(...)
////////////////////////////////////////////////////////////////////////////////
//
// Overlays the uncorrected vs corrected local-φ distributions from two TH3F
// histograms, each shaped like (local-η binning, local-φ binning, pT binning).
//
// For each of the 4 pT bins, we:
//   1) Project the uncorrected histogram in pT range [pT_lo, pT_hi] onto "y"
//      (the local‐φ axis), producing TH1D hUnc.
//   2) Project the corrected histogram in the same pT range, producing TH1D hCorr.
//   3) Normalize each, overlay them on the same pad.
//   4) Add a legend labeling “Uncorrected” vs “Corrected” + pT range in the title.
//
// We then produce a single TCanvas with a 2×2 layout (one pad per pT bin) and
// save it to the specified “outPNG” file.
//
// NOTE: This version does *no fits*, so no b-parameters are computed or displayed.
//
////////////////////////////////////////////////////////////////////////////////
void overlayUncorrCorrLocalPhi_noFits(
    TH3F* h3_uncorr,
    TH3F* h3_corr,
    const double pT_low[4],
    const double pT_high[4],
    const TString& outPNG
)
{
  // Basic checks
  if(!h3_uncorr || !h3_corr)
  {
    std::cerr << "[ERROR] overlayUncorrCorrLocalPhi_noFits: null histogram pointer.\n";
    return;
  }

  // We'll produce a 2×2 TCanvas for four pT bins
  TCanvas c2by2("cLocalPhiUncorrCorr","Local #phi Overlay: Uncorr vs Corr",1200,1000);
  c2by2.Divide(2,2);

  // We assume exactly 4 pT bins
  const int N_BINS = 4;

  // For each pT bin => fill one pad
  for(int i=0; i<N_BINS; i++)
  {
    c2by2.cd(i+1);

    // pT range
    double pTlo = pT_low[i];
    double pThi = pT_high[i];

    //----------------------------------------------------------------------
    // (A) Project UNCORRECTED in that pT range
    //----------------------------------------------------------------------
    // 1) set the Z axis range for pT in [pTlo, pThi]
    int zLo = h3_uncorr->GetZaxis()->FindFixBin(pTlo + 1e-9);
    int zHi = h3_uncorr->GetZaxis()->FindFixBin(pThi - 1e-9);
    if(zLo < 1) zLo = 1;
    if(zHi > h3_uncorr->GetNbinsZ()) zHi = h3_uncorr->GetNbinsZ();
    h3_uncorr->GetZaxis()->SetRange(zLo, zHi);

    // 2) Project onto "y" => local φ
    TH1D* hUnc = (TH1D*) h3_uncorr->Project3D("y");
    if(!hUnc)
    {
      std::cerr << "[WARN] No uncorrected projection for pT bin " << i
                << " => skip.\n";
      continue;
    }

    // Normalize
    double intUnc = hUnc->Integral();
    if(intUnc > 1e-9) hUnc->Scale(1./intUnc);

    // Style
    hUnc->SetName( Form("hUnc_phi_ptBin%d", i) );
    hUnc->SetLineColor(kBlue+1);
    hUnc->SetMarkerColor(kBlue+1);
    hUnc->SetMarkerStyle(20);
    hUnc->SetMarkerSize(1);
    hUnc->SetTitle( Form("Local #phi: p_{T}=[%.1f,%.1f]", pTlo, pThi) );
    hUnc->GetXaxis()->SetTitle("local #phi");
    hUnc->GetYaxis()->SetTitle("counts (normalized)");

    // Draw first
    hUnc->Draw("E");
    double maxUnc = hUnc->GetMaximum();
    hUnc->GetYaxis()->SetRangeUser(0, 1.2 * maxUnc);

    //----------------------------------------------------------------------
    // (B) Project CORRECTED in the same pT range
    //----------------------------------------------------------------------
    zLo = h3_corr->GetZaxis()->FindFixBin(pTlo + 1e-9);
    zHi = h3_corr->GetZaxis()->FindFixBin(pThi - 1e-9);
    if(zLo < 1) zLo=1;
    if(zHi > h3_corr->GetNbinsZ()) zHi= h3_corr->GetNbinsZ();
    h3_corr->GetZaxis()->SetRange(zLo, zHi);

    TH1D* hCorr = (TH1D*) h3_corr->Project3D("y");
    if(!hCorr)
    {
      std::cerr << "[WARN] No corrected projection for pT bin " << i
                << " => skipping.\n";
      continue;
    }

    // Normalize
    double intCorr = hCorr->Integral();
    if(intCorr > 1e-9) hCorr->Scale(1./intCorr);

    // Style
    hCorr->SetName( Form("hCorr_phi_ptBin%d", i) );
    hCorr->SetLineColor(kRed);
    hCorr->SetMarkerColor(kRed);
    hCorr->SetMarkerStyle(21);
    hCorr->SetMarkerSize(1);

    // Overdraw
    hCorr->Draw("SAME E");

    //----------------------------------------------------------------------
    // (C) Legend
    //----------------------------------------------------------------------
    TLegend leg(0.50, 0.65, 0.88, 0.85);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hUnc, "Uncorrected", "lp");
    leg.AddEntry(hCorr,"Corrected",   "lp");
    leg.Draw();
  }

  //----------------------------------------------------------------------
  // Save the 2×2 layout
  //----------------------------------------------------------------------
  c2by2.SaveAs(outPNG);
  std::cout << "[INFO] overlayUncorrCorrLocalPhi_noFits => wrote 2×2 overlay to "
            << outPNG << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// PDCanalysis()
//
//  1) Lists all histograms (like the original code).
//  2) Makes a 2D block-eta vs block-phi plot from h2_cluster_block_cord_pt
//  3) Makes local-phi distributions for pT slices [2,4], [4,6], [6,8],
//     fits them, and records b-values.
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
  // 2) Grab the 3D histogram for block coords vs pT
  ///////////////////////////////////////////////////////////////////////
  TH3F* h3 = (TH3F*) f->Get("h2_cluster_block_cord_pt");
    
  if(!h3)
  {
    std::cerr << "[ERROR] 'h2_cluster_block_cord_pt' not found => can't do block-eta vs. block-phi.\n";
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
  // 3) Create 2D plot from h3 by integrating over entire pT
  ///////////////////////////////////////////////////////////////////////
    {
      // First, restore the Z-axis range to include all pT
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
      std::cout << "[INFO] Wrote 2D block-eta vs block-phi (all pT) => "
                << (outDir2D + "/BlockEtaPhi_2D.png") << std::endl;
    }

    ///////////////////////////////////////////////////////////////////
    // 3b) Also produce 2D plots for each pT slice [2,4], [4,6], [6,8]
    ///////////////////////////////////////////////////////////////////
    {
      for(int i=0; i<N_PT; i++)
      {
        // Determine Z-range in the 3D histogram (pT axis)
        int zLo = h3->GetZaxis()->FindFixBin( pTedges[i] + 1e-9 );
        int zHi = h3->GetZaxis()->FindFixBin( pTedges[i+1] - 1e-9 );
        if(zLo < 1) zLo = 1;
        if(zHi > h3->GetNbinsZ()) zHi = h3->GetNbinsZ();

        // Apply that Z-range, then project onto X-Y => 2D
        h3->GetZaxis()->SetRange(zLo, zHi);
        TH2F* h2_slice = (TH2F*) h3->Project3D("xy");

        // Name & title
        h2_slice->SetName( Form("h2_blockCoord2D_pt%.0fto%.0f", pTedges[i], pTedges[i+1]) );
        h2_slice->SetTitle( Form("Block #eta vs Block #phi (%.1f < p_{T} < %.1f);block #eta;block #phi",
                                   pTedges[i], pTedges[i+1]) );

        // Draw
        TCanvas c2dSlice( Form("c2d_pt%.0fto%.0f", pTedges[i], pTedges[i+1]),
                            "Block Eta-Phi (pT slice)", 900, 700);
        c2dSlice.SetLeftMargin(0.12);
        c2dSlice.SetRightMargin(0.15);
        c2dSlice.SetBottomMargin(0.12);

        h2_slice->SetStats(false);
        h2_slice->Draw("COLZ");

        // Save a PNG for each slice
        TString outName = Form("/BlockEtaPhi_2D_pt%.0fto%.0f.png",
                                 pTedges[i], pTedges[i+1]);
        c2dSlice.SaveAs(outDir2D + outName);

        std::cout << "[INFO] Wrote 2D block-eta vs block-phi => "
                  << (outDir2D + outName)
                  << " for pT=[" << pTedges[i] << ", " << pTedges[i+1] << "]" << std::endl;
      }
    }

  ///////////////////////////////////////////////////////////////////////
  // 4) Local φ slices in pT bins => do 1D projections, fits
  ///////////////////////////////////////////////////////////////////////

  // We'll record b-values to a text file
  std::ofstream bResults(Form("%s/bValues.txt", baseDir.Data()));
  if(!bResults.is_open()) {
    std::cerr << "[WARNING] Could not open bValues.txt for writing.\n";
  } else {
    bResults << "# pTlow  pThigh   bValue\n";
  }

  double pT_bins_low[4], pT_bins_high[4];
  for(int j=0; j<4; j++){
      pT_bins_low[j]  = pTedges[j];
      pT_bins_high[j] = pTedges[j+1];
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
  //
  //    We'll do a second iteration over keys, create a TCanvas,
  //    draw the histogram, and SaveAs(...). We'll do some minimal
  //    logic to handle TH1 vs TH2 vs TH3 vs TProfile.
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
    //   TH2 => "COLZ"
    //   TH3 => "BOX" or skip
    //   TProfile => normal "E1"
    //   TH1 => normal "E"
    //   We'll keep it simple
    if( htmp->InheritsFrom("TH2") ) {
      htmp->Draw("COLZ");
    }
    else if(htmp->InheritsFrom("TH3")) {
      // 3D is tricky to visualize. We'll do "BOX".
      // Or skip if you prefer. We'll just do BOX.
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
    // cleanup
    delete htmp;
  }
  plotAshLogRMS_sideBySide(filename);

    ////////////////////////////////////////////
    // (X) Overlays: uncorrected vs. corrected
    ////////////////////////////////////////////
    // 1) Retrieve both histograms from the file.
    //    Suppose the uncorrected is "h2_cluster_block_cord_pt"
    //    and the corrected is "h3_cluster_block_cord_corr_pt".
    //    Make sure these names match the ones in your code!
    TH3F* h3_uncorr = dynamic_cast<TH3F*>( f->Get("h2_cluster_block_cord_pt") );
    TH3F* h3_corr   = dynamic_cast<TH3F*>( f->Get("h3_cluster_block_cord_corr_pt") );
    if(!h3_uncorr || !h3_corr)
    {
      std::cerr << "[WARNING] Could not retrieve uncorrected or corrected TH3F. Skipping overlay." << std::endl;
    }
    else
    {
      // 2) Define four pT bins. For example:
      double pT_low[4]  = {2.0,  4.0,  6.0,  8.0};
      double pT_high[4] = {4.0,  6.0,  8.0, 12.0};

      // 3) Decide the output PNG name
      TString overlayPNG = baseDir + "/LocalPhi_Overlay_UncorrVsCorr_noFits.png";

      // 4) Call the function (make sure you've #included or have it in scope)
      overlayUncorrCorrLocalPhi_noFits(
        h3_uncorr,    // uncorrected TH3F
        h3_corr,      // corrected TH3F
        pT_low,       // array of 4 low edges
        pT_high,      // array of 4 high edges
        overlayPNG    // output PNG
      );

      std::cout << "[INFO] Done overlaying uncorrected vs corrected local-phi distributions.\n"
                << "       See: " << overlayPNG << std::endl;
    }

    
  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
