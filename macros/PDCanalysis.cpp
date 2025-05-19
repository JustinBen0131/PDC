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

void plotAshLogRMS_sideBySide(const char* infile = "PositionDep_sim_ALL.root")
{
  // -------------------------------------------------------------------------
  // 1) Configuration
  // -------------------------------------------------------------------------
  constexpr int    N_PT   = 4;
  const     double ptEdge[N_PT+1] = {2.0, 3.0, 5.0, 8.0, 12.0};

  // Ash (b) scan
  const double bMin  = 0.05, bMax = 0.30, bStep = 0.01;
  std::vector<double> bScan;
  for(double b=bMin; b<=bMax+1e-9; b+=bStep) bScan.push_back(b);
  const int N_B = bScan.size();

  // Log (w0) scan
  const double w0Min = 2.8,  w0Max = 5.5, w0Step = 0.10;
  std::vector<double> w0Scan;
  for(double w=w0Min; w<=w0Max+1e-9; w+=w0Step) w0Scan.push_back(w);
  const int N_W = w0Scan.size();

  // Output location
  const TString baseDir    = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
  const TString histOutDir = baseDir + "/ASH_LOG_PLOTS";
  gSystem->mkdir(baseDir,    /*parents=*/true);
  gSystem->mkdir(histOutDir, /*parents=*/true);

  // -------------------------------------------------------------------------
  // 2) Input file
  // -------------------------------------------------------------------------
  std::cout << "[INFO] Opening '" << infile << "'...\n";
  std::unique_ptr<TFile> f(TFile::Open(infile,"READ"));
  if(!f || f->IsZombie()){
    std::cerr << "[ERROR] Cannot open file.\n"; return;
  }

  // -------------------------------------------------------------------------
  // 3) ROOT style
  // -------------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);

  int  colors [4] = {kBlack, kRed, kBlue, kMagenta+2};
  int  markers[4] = {20,      21,   22,    23       };

  // ======================================================================
  // 4)  ASH  : σx  vs  b
  // ======================================================================
  TCanvas cAsh("cAsh","Ash RMS vs b",800,600);
  cAsh.SetLeftMargin(0.12);  cAsh.SetBottomMargin(0.12);  cAsh.SetGrid();

  double gMinAsh =  DBL_MAX, gMaxAsh = -DBL_MAX;
  std::vector<TGraph*> gAshVec;
  std::vector<double>  bestB (N_PT,0.), bestRMS_A(N_PT,DBL_MAX);

  for(int ipt=0; ipt<N_PT; ++ipt)
  {
    auto g = new TGraph; g->SetName(Form("gAsh_pt%d",ipt));
    for(size_t ib=0; ib<bScan.size(); ++ib)
    {
      double bVal = bScan[ib];
      TString hN  = Form("h_dx_ash_b%04.2f_pt%d",bVal,ipt);
      if(TH1* h = dynamic_cast<TH1*>(f->Get(hN)))
      {
        double r = h->GetRMS();
        g->SetPoint(g->GetN(), bVal, r);

        if(r>0){
          gMinAsh = std::min(gMinAsh,r);
          gMaxAsh = std::max(gMaxAsh,r);
          if(r<bestRMS_A[ipt]){ bestRMS_A[ipt]=r; bestB[ipt]=bVal; }
        }

        // quick PNG snapshot
        TCanvas ctmp; h->Draw("E"); ctmp.SaveAs(histOutDir+"/"+hN+".png");
      }
    }
    g->SetMarkerStyle(markers[ipt]);
    g->SetMarkerColor(colors [ipt]);
    g->SetLineColor  (colors [ipt]);
    g->SetMarkerSize(1.3);
    gAshVec.push_back(g);
  }

  // guard: empty?
  if(gMinAsh==DBL_MAX){ gMinAsh=0; gMaxAsh=1; }

  // draw onto the *correct* pad
  cAsh.cd();
  gAshVec[0]->Draw("ALP");
  gAshVec[0]->GetYaxis()->SetRangeUser(gMinAsh-0.1*fabs(gMinAsh),
                                       gMaxAsh+0.1*fabs(gMaxAsh));
  for(size_t i=1;i<gAshVec.size();++i) gAshVec[i]->Draw("LP SAME");

  TLegend legA(0.15,0.70,0.45,0.85); legA.SetBorderSize(0); legA.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt){
    legA.AddEntry(gAshVec[ipt],
      Form("%.1f<p_{T}<%.1f GeV  (best b=%.2f)",ptEdge[ipt],ptEdge[ipt+1],bestB[ipt]),
      "lp");
  }
  legA.Draw();
  TLatex().DrawLatexNDC(0.2,0.92,"Ash");

  cAsh.SaveAs(baseDir+"/Ash_RMS_vs_b.png");
  std::cout << "[INFO] Saved Ash_RMS_vs_b.png\n";

  // ======================================================================
  // 5)  LOG  : σx  vs  w0
  // ======================================================================
  TCanvas cLog("cLog","Log RMS vs w0",800,600);
  cLog.SetLeftMargin(0.12);  cLog.SetBottomMargin(0.12);  cLog.SetGrid();

  double gMinLog = DBL_MAX, gMaxLog = -DBL_MAX;
  std::vector<TGraph*> gLogVec;
  std::vector<double>  bestW (N_PT,0.), bestRMS_L(N_PT,DBL_MAX);

  for(int ipt=0; ipt<N_PT; ++ipt)
  {
    auto g = new TGraph; g->SetName(Form("gLog_pt%d",ipt));
    for(size_t iw=0; iw<w0Scan.size(); ++iw)
    {
      double wVal = w0Scan[iw];
      TString hN  = Form("h_dx_log_w0%04.2f_pt%d",wVal,ipt);
      if(TH1* h = dynamic_cast<TH1*>(f->Get(hN)))
      {
        double r = h->GetRMS();
        g->SetPoint(g->GetN(), wVal, r);

        if(r>0){
          gMinLog = std::min(gMinLog,r);
          gMaxLog = std::max(gMaxLog,r);
          if(r<bestRMS_L[ipt]){ bestRMS_L[ipt]=r; bestW[ipt]=wVal; }
        }

        TCanvas ctmp; h->Draw("E"); ctmp.SaveAs(histOutDir+"/"+hN+".png");
      }
    }
    g->SetMarkerStyle(markers[ipt]);
    g->SetMarkerColor(colors [ipt]);
    g->SetLineColor  (colors [ipt]);
    g->SetMarkerSize(1.3);
    gLogVec.push_back(g);
  }

  if(gMinLog==DBL_MAX){ gMinLog=0; gMaxLog=1; }

  // *** critical: make sure we're on the main canvas, not a stale ctmp ***
  cLog.cd();
  gLogVec[0]->Draw("ALP");
  gLogVec[0]->GetYaxis()->SetRangeUser(gMinLog-0.1*fabs(gMinLog),
                                       gMaxLog+0.1*fabs(gMaxLog));
  for(size_t i=1;i<gLogVec.size();++i) gLogVec[i]->Draw("LP SAME");

  TLegend legB(0.15,0.65,0.45,0.85); legB.SetBorderSize(0); legB.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt){
    legB.AddEntry(gLogVec[ipt],
      Form("%.1f<p_{T}<%.1f GeV  (best w_{0}=%.2f)",ptEdge[ipt],ptEdge[ipt+1],bestW[ipt]),
      "lp");
  }
  legB.Draw();
  TLatex().DrawLatexNDC(0.2,0.92,"Log");

  cLog.SaveAs(baseDir+"/Log_RMS_vs_w0.png");
  std::cout << "[INFO] Saved Log_RMS_vs_w0.png\n";

  std::cout << "[DONE] plotAshLogRMS_sideBySide finished.\n";
}

////////////////////////////////////////////////////////////////////////////////
// Asinh-based function for local phi fits:
//
//   Y(X) = Normalization * [ 2*b / sqrt(1 + 4*X^2 * sinh^2(1/(2*b))) ]
//
// We'll define a TF1 that has two parameters: [0]=Normalization, [1]=b
////////////////////////////////////////////////////////////////////////////////
double asinhModel(double *x, double *par)
{
  double Norm = par[0];  // overall scale
  double b    = par[1];  // the "b" parameter we want to extract

  // If b is extremely small or zero, avoid division by zero
  if (b <= 1e-9) return 1e-12;

  double Xcg = x[0];
  double arg = 1.0 + 4.0 * Xcg * Xcg * TMath::SinH(1.0/(2.0*b)) * TMath::SinH(1.0/(2.0*b));
  double denom = TMath::Sqrt(arg);

  double numer = 2.0 * b;
  double result = Norm * ( numer / denom );
  return result;
}
////////////////////////////////////////////////////////////////////////////////
// doLocalPhiEtaFits(...)
//
// Same as your original function for local φ/η slices. Now includes enhanced
// fit convergence by doing multiple initial guesses, setting parameter limits,
// and using "RQE" for the fit option to ensure more robust minimization.
//
//  - If the best attempt’s χ² is still poor, you could optionally do more
//    advanced logic (like further initial guesses or refine parLimits).
//  - Everything else (legend, single‐bin overlays for local‐φ, b‐value text file
//    labeling, etc.) is kept the same.
////////////////////////////////////////////////////////////////////////////////
void doLocalPhiEtaFits(TH3F* h3,
                       const double pT_bins_low[3],
                       const double pT_bins_high[3],
                       std::ofstream& bResults,
                       const TString& outDirPhi)
{
  // Ensure directory exists
  gSystem->mkdir(outDirPhi, /*recursive=*/true);

  /////////////////////////////////////////////////////////////////////////////
  // Two equation strings: "b_{#phi}" vs. "b_{#eta}"
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

  // Colors/markers for the three pT slices
  int colors[3]  = {kBlack, kRed, kBlue};
  int markers[3] = {20, 21, 22};

  // -------------------------------------------------------------------------
  // Helper lambda for the triple overlay, plus single‐bin overlays for φ
  // -------------------------------------------------------------------------
  auto fitCoordAndMakePlot =
    [&](const char* coordName,    // "phi" or "eta"
        const char* projectOpt,   // "y" for φ, "x" for η
        const TString& outPNGname,
        const TString& xTitle)
  {
    // Decide which equation string to use
    bool isPhi = (strcmp(coordName, "phi") == 0);
    TString eqString = (isPhi ? eqPhiString : eqEtaString);

    // Label in the text file
    if(bResults.is_open())
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

    // 1) Create a canvas for the triple overlay
    TCanvas cCoord(Form("cCoord_%s",coordName),
                   Form("Local %s distributions",coordName),
                   900,700);
    cCoord.SetLeftMargin(0.12);
    cCoord.SetRightMargin(0.05);
    cCoord.SetBottomMargin(0.12);
    cCoord.cd();

    TLegend leg(0.18,0.68,0.4,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.028);

    // We'll store the 3 uncorrected histograms to do single‐bin overlays
    std::vector<TH1D*> hLocalVec(3, nullptr);

    // 2) Loop over the 3 pT slices
    for(int i=0; i<3; i++)
    {
      int zLo = h3->GetZaxis()->FindFixBin(pT_bins_low[i] + 1e-9);
      int zHi = h3->GetZaxis()->FindFixBin(pT_bins_high[i] - 1e-9);
      if(zLo<1)   zLo=1;
      if(zHi>h3->GetNbinsZ()) zHi=h3->GetNbinsZ();

      // Restrict h3 to that pT slice
      h3->GetXaxis()->SetRange(1, h3->GetNbinsX());
      h3->GetZaxis()->SetRange(zLo, zHi);

      // Project onto local φ or η
      TH1D* hLocal = (TH1D*) h3->Project3D(projectOpt);
      if(!hLocal)
      {
        std::cerr << "[WARNING] Could not project " << coordName
                  << " pT bin " << i << std::endl;
        continue;
      }
      hLocalVec[i] = hLocal;

      hLocal->SetName(Form("hLocal%s_%.1fto%.1f",
                           coordName, pT_bins_low[i], pT_bins_high[i]));
      hLocal->SetTitle(Form(";%s; scaled counts", xTitle.Data()));

      // Normalize
      double integral = hLocal->Integral();
      if(integral > 1e-9) hLocal->Scale(1./integral);

      // Style
      hLocal->SetMarkerStyle(markers[i]);
      hLocal->SetMarkerColor(colors[i]);
      hLocal->SetLineColor(colors[i]);
      hLocal->SetMarkerSize(1.2);

      // Draw
      if(i == 0)
      {
        hLocal->Draw("E");
        hLocal->GetYaxis()->SetRangeUser(0, 0.4);
        hLocal->GetXaxis()->SetTitle(xTitle);
        hLocal->GetYaxis()->SetTitle("counts (normalized)");
      }
      else
      {
        hLocal->Draw("same E");
      }

      // ---------------------------------------------------------------------
      // Enhanced multi‐start fit approach to ensure robust convergence
      // ---------------------------------------------------------------------
      // We'll define a TF1 with param limits for bVal > 0.0, etc.
      TF1* fAsinh = new TF1(Form("fAsinh_%s_%d",coordName,i),
                            asinhModel, -0.5, 0.5, 2);
      fAsinh->SetParNames("Norm", "bVal");
      // Some parameter limits:
      // Norm can be from near zero up to large; bVal from e.g. 0.01 to 1.0
      // (If your bVal can be bigger than 1, increase the upper limit.)
      fAsinh->SetParLimits(0, 1e-6, 1e6);
      fAsinh->SetParLimits(1, 1e-5,  1.0);

      // We'll do a small multi‐start approach for the initial guesses:
      // (Norm, bVal)
      std::vector<std::pair<double,double>> initialGuesses = {
        {0.2, 0.14},
        {0.1, 0.10},
        {0.3, 0.20},
        {0.5, 0.08}
      };

      double bestChi2   = 1e15;
      double bestNorm   = 0.2;
      double bestBVal   = 0.14;
      double fitMinX    = -0.5;
      double fitMaxX    =  0.5;

      // We'll store the best TF1 from the tries
      TF1* bestFunc = nullptr;

      // Try each guess
      for(auto& guess : initialGuesses)
      {
        double guessNorm = guess.first;
        double guessB    = guess.second;

        // Set the guess
        fAsinh->SetParameter(0, guessNorm);
        fAsinh->SetParameter(1, guessB);

        // Perform the fit quietly + error estimation
        // "R" = use fit range, "Q"=quiet, "E"=error estimates
        TFitResultPtr fitRes = hLocal->Fit(fAsinh, "RQE0S", "",
                                           fitMinX, fitMaxX);

        // Check the resulting χ²
        if( fitRes.Get() && fitRes->IsValid() )
        {
          double thisChi2 = fitRes->Chi2();
          if(thisChi2 < bestChi2)
          {
            bestChi2 = thisChi2;
            bestNorm = fAsinh->GetParameter(0);
            bestBVal = fAsinh->GetParameter(1);

            // Keep a copy of the function
            if(bestFunc) delete bestFunc;
            bestFunc = (TF1*) fAsinh->Clone(Form("bestFunc_%s_%d",
                                                 coordName, i));
          }
        }
      } // end multi‐start

      // If we found a bestFunc, draw it & add to legend
      double bVal = 0.0;
      if(bestFunc)
      {
        // Update style
        bestFunc->SetLineColor(colors[i]);
        bestFunc->SetLineWidth(2);
        bestFunc->Draw("same");

        bVal = bestFunc->GetParameter(1);

        // We can optionally store it, or delete it:
        // For brevity, let's store it in the histogram's list so
        // it doesn't vanish when the function scope ends:
        hLocal->GetListOfFunctions()->Add(bestFunc);
      }
      else
      {
        std::cerr << "[WARNING] No valid multi‐start fit for pT bin i=" << i
                  << " => bVal=0.\n";
      }

      // Legend entry with final best bVal
      TString legEntry = Form("p_{T}=[%.0f,%.0f] GeV : b=%.3g",
                              pT_bins_low[i], pT_bins_high[i], bVal);
      leg.AddEntry(hLocal, legEntry, "lp");

      // Write b-value to file
      if(bResults.is_open())
      {
        bResults << pT_bins_low[i] << "  "
                 << pT_bins_high[i] << "   "
                 << bVal << "\n";
      }
    } // end loop over i

    // Draw eq label
    leg.Draw();
    eqLatex.DrawLatex(0.57, 0.82, eqString);

    // Save triple overlay
    TString mainOut = outDirPhi + "/" + outPNGname;
    cCoord.SaveAs(mainOut);
    std::cout << "[INFO] Wrote triple-overlay " << coordName
              << " => " << mainOut << std::endl;

    // -----------------------------------------------------------------------
    // 3) If local-φ, produce single‐bin overlays with corrected hists
    // -----------------------------------------------------------------------
    if(isPhi)
    {
      for(int i=0; i<3; i++)
      {
        TH1D* hUnc = hLocalVec[i];
        if(!hUnc) continue;

        // The corrected histogram name
        TString corrName = Form("h_localPhi_corrected_%.0f_%.0f",
                                pT_bins_low[i], pT_bins_high[i]);
        TH1F* hCorr = dynamic_cast<TH1F*>(gROOT->FindObject(corrName));
        if(!hCorr)
        {
          std::cerr << "[WARNING] No '"
                    << corrName << "' => skipping single-bin overlay for pT bin "
                    << i << "\n";
          continue;
        }

        // Build a new canvas
        TString cName  = Form("cLocalPhiCompare_%.1fto%.1f",
                              pT_bins_low[i], pT_bins_high[i]);
        TString cTitle = Form("Local-#phi overlay [%.1f< pT<%.1f]",
                              pT_bins_low[i], pT_bins_high[i]);
        TCanvas cSingle(cName, cTitle, 800,600);
        cSingle.SetLeftMargin(0.12);
        cSingle.SetRightMargin(0.05);
        cSingle.SetBottomMargin(0.12);
        cSingle.cd();

        // Clone the uncorrected so we can re‐normalize
        TH1D* hUncClone = (TH1D*)hUnc->Clone(Form("%s_clone", hUnc->GetName()));

        double uncInt  = hUncClone->Integral();
        double corrInt = hCorr->Integral();
        if(uncInt  > 1e-9) hUncClone->Scale(1./uncInt);
        if(corrInt > 1e-9) hCorr->Scale  (1./corrInt);

        // Style uncorrected
        hUncClone->SetLineColor(kBlack);
        hUncClone->SetMarkerColor(kBlack);
        hUncClone->SetMarkerStyle(20);
        hUncClone->SetMarkerSize(1);
        hUncClone->SetTitle(cTitle);
        hUncClone->GetXaxis()->SetTitle("local #phi_{CG} in block");
        hUncClone->GetYaxis()->SetTitle("counts (normalized)");
        hUncClone->GetYaxis()->SetTitleOffset(1.2);

        // Style corrected
        hCorr->SetLineColor(kRed);
        hCorr->SetMarkerColor(kRed);
        hCorr->SetMarkerStyle(21);
        hCorr->SetMarkerSize(1);

        // Draw uncorrected
        hUncClone->Draw("E");
        double maxUnc  = hUncClone->GetMaximum();
        double maxCorr = hCorr->GetMaximum();
        double newMax  = TMath::Max(maxUnc, maxCorr)*1.25;
        hUncClone->GetYaxis()->SetRangeUser(0, newMax);

        // Overlay corrected
        hCorr->Draw("SAME E");

        // Legend
        TLegend legSingle(0.55, 0.70, 0.85, 0.85);
        legSingle.SetBorderSize(0);
        legSingle.SetFillStyle(0);
        legSingle.SetTextSize(0.035);
        legSingle.AddEntry(hUncClone,
                Form("Raw local #phi (%.1f< p_{T}<%.1f)",
                     pT_bins_low[i], pT_bins_high[i]),
                "lp");
        legSingle.AddEntry(hCorr,
                Form("Corrected local #phi (%.1f< p_{T}<%.1f)",
                     pT_bins_low[i], pT_bins_high[i]),
                "lp");
        legSingle.Draw();

        // Save
        TString singleOut = Form("%s/LocalPhiCompare_%.0fto%.0f.png",
                                 outDirPhi.Data(), pT_bins_low[i], pT_bins_high[i]);
        cSingle.SaveAs(singleOut);
        std::cout << "[INFO] Single-bin local-#phi overlay pT=["
                  << pT_bins_low[i] << "," << pT_bins_high[i]
                  << "] => " << singleOut << std::endl;
      }
    } // end if(isPhi)
  }; // end lambda

  // -----------------------------------------------------------------------
  // 4) Call the helper for φ => “LocalPhiFits.png”
  // -----------------------------------------------------------------------
  fitCoordAndMakePlot(
    "phi",  // coordName
    "y",    // local φ is y-axis in TH3
    "LocalPhiFits.png",
    "local #phi_{CG} in block"
  );

  // -----------------------------------------------------------------------
  // 5) Then for η => “LocalEtaFits.png”
  // -----------------------------------------------------------------------
  fitCoordAndMakePlot(
    "eta",  // coordName
    "x",    // local η is x-axis in TH3
    "LocalEtaFits.png",
    "local #eta_{CG} in block"
  );
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
      double pT_bins_low[3]  = {2.0, 4.0, 6.0};
      double pT_bins_high[3] = {4.0, 6.0, 8.0};

      for(int i=0; i<3; i++)
      {
        // Determine Z-range in the 3D histogram (pT axis)
        int zLo = h3->GetZaxis()->FindFixBin(pT_bins_low[i] + 1e-9);
        int zHi = h3->GetZaxis()->FindFixBin(pT_bins_high[i] - 1e-9);
        if(zLo < 1) zLo = 1;
        if(zHi > h3->GetNbinsZ()) zHi = h3->GetNbinsZ();

        // Apply that Z-range, then project onto X-Y => 2D
        h3->GetZaxis()->SetRange(zLo, zHi);
        TH2F* h2_slice = (TH2F*) h3->Project3D("xy");

        // Name & title
        h2_slice->SetName( Form("h2_blockCoord2D_pt%.0fto%.0f",
                                pT_bins_low[i], pT_bins_high[i]) );
        h2_slice->SetTitle( Form("Block #eta vs Block #phi (%.1f < p_{T} < %.1f);block #eta;block #phi",
                                 pT_bins_low[i], pT_bins_high[i]) );

        // Draw
        TCanvas c2dSlice( Form("c2d_pt%.0fto%.0f", pT_bins_low[i], pT_bins_high[i]),
                          "Block Eta-Phi (pT slice)", 900, 700);
        c2dSlice.SetLeftMargin(0.12);
        c2dSlice.SetRightMargin(0.15);
        c2dSlice.SetBottomMargin(0.12);

        h2_slice->SetStats(false);
        h2_slice->Draw("COLZ");

        // Save a PNG for each slice
        TString outName = Form("/BlockEtaPhi_2D_pt%.0fto%.0f.png",
                               pT_bins_low[i], pT_bins_high[i]);
        c2dSlice.SaveAs(outDir2D + outName);

        std::cout << "[INFO] Wrote 2D block-eta vs block-phi => "
                  << (outDir2D + outName)
                  << " for pT=[" << pT_bins_low[i] << ", "
                  << pT_bins_high[i] << "]" << std::endl;
      }
    }

  ///////////////////////////////////////////////////////////////////////
  // 4) Local φ slices in pT bins => do 1D projections, fits
  ///////////////////////////////////////////////////////////////////////
  double pT_bins_low[3]  = {2.0, 4.0, 6.0};
  double pT_bins_high[3] = {4.0, 6.0, 8.0};

  // We'll record b-values to a text file
  std::ofstream bResults(Form("%s/bValues.txt", baseDir.Data()));
  if(!bResults.is_open()) {
    std::cerr << "[WARNING] Could not open bValues.txt for writing.\n";
  } else {
    bResults << "# pTlow  pThigh   bValue\n";
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

    
  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
