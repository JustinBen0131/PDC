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

constexpr double pTedges[] = {2.0, 3.0, 5.0, 8.0, 12.0};
// Let the compiler figure out how many bins that is:
constexpr int N_PT = (sizeof(pTedges)/sizeof(pTedges[0])) - 1;

void plotAshLogRMS_sideBySide(const char* infile = "PositionDep_sim_ALL.root")
{
  // Ash (b) scan

  const double bMin  = 0.18;
  const double bMax  = 0.32;
  const double bStep = 0.005;
  std::vector<double> bScan;
  for(double b = bMin; b <= bMax + 1e-9; b += bStep) bScan.push_back(b);
  const int N_B = bScan.size();


  // Log (w0) scan
  const double w0Min  = 3.0;
  const double w0Max  = 5.0;
  const double w0Step = 0.05;
  std::vector<double> w0Scan;
  for(double w = w0Min; w <= w0Max + 1e-9; w += w0Step) w0Scan.push_back(w);
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

    std::vector<int> colors;
    std::vector<int> markers;

    // For example, fill them in at runtime:
    {
      // This is just an example color set of length = N_PT
      // If you have fewer or more bins, it adjusts automatically
      const int baseColorList[] = {
        kBlack, kRed, kBlue, kMagenta+2, kOrange,
        kGreen+2, kCyan+2, kGray+1, kViolet+1, kAzure+2, kSpring+9
      };
      // We'll cycle through them if N_PT > 11, or just use the first N_PT if fewer.
      for (int i=0; i<N_PT; ++i) {
        colors.push_back( baseColorList[i % 11] );
      }

      // Similarly for marker styles:
      const int baseMarkerList[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
      for (int i=0; i<N_PT; ++i) {
        markers.push_back( baseMarkerList[i % 11] );
      }
    }
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
      TString hN  = Form("h_dx_ash_b%.3f_pt%d",   bVal, ipt);
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
  gAshVec[0]->GetXaxis()->SetTitle("b");
  gAshVec[0]->GetYaxis()->SetTitle("#sigma_{x}");
  gAshVec[0]->GetYaxis()->SetRangeUser(gMinAsh-0.1*fabs(gMinAsh),
                                       gMaxAsh+0.1*fabs(gMaxAsh));
  for(size_t i=1;i<gAshVec.size();++i) gAshVec[i]->Draw("LP SAME");

  TLegend legA(0.15,0.70,0.45,0.85); legA.SetBorderSize(0); legA.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt){
      legA.AddEntry(gAshVec[ipt],
                    Form("%.1f<p_{T}<%.1f GeV  (best b=%.2f)",
                         pTedges[ipt], pTedges[ipt+1], bestB[ipt]),
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
      TString hN  = Form("h_dx_log_w0%.2f_pt%d",wVal,ipt);
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
  gLogVec[0]->GetXaxis()->SetTitle("w_{0}");
  gLogVec[0]->GetYaxis()->SetTitle("#sigma_{x}");
  gLogVec[0]->GetYaxis()->SetRangeUser(gMinLog-0.1*fabs(gMinLog),
                                       gMaxLog+0.1*fabs(gMaxLog));
  for(size_t i=1;i<gLogVec.size();++i) gLogVec[i]->Draw("LP SAME");

  TLegend legB(0.15,0.65,0.45,0.85); legB.SetBorderSize(0); legB.SetFillStyle(0);
  for(int ipt=0; ipt<N_PT; ++ipt){
      legB.AddEntry(gLogVec[ipt],
                    Form("%.1f<p_{T}<%.1f GeV  (best w_{0}=%.2f)",
                         pTedges[ipt], pTedges[ipt+1], bestW[ipt]),
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
  double b    = par[1];  // the "b" parameter
  if (b <= 1e-9) return 1e-12; // avoid division by zero if b is small

  double Xcg = x[0];
  double arg = 1.0 + 4.0*Xcg*Xcg*TMath::SinH(1.0/(2.0*b))*TMath::SinH(1.0/(2.0*b));
  double denom = TMath::Sqrt(arg);

  double numer = 2.0 * b;
  double result = Norm * (numer/denom);
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// doLocalPhiEtaFits(...)
//
//  1) Projects local-φ or local-η from h3 in slices of pT
//  2) Overlays all 4 bins on one canvas => LocalPhiFits.png or LocalEtaFits.png
//  3) If coordinate == phi, also do single‐bin overlays raw vs corrected
//  4) **Now** also produce a 2×2 table of the final per‐bin fits
////////////////////////////////////////////////////////////////////////////////
void doLocalPhiEtaFits(TH3F* h3,
                       const double pT_bins_low[4],
                       const double pT_bins_high[4],
                       std::ofstream& bResults,
                       const TString& outDirPhi)
{
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

  const int NPT = 11;          // ← was 6


  int  colors [NPT] = { kBlack, kRed+1, kBlue+1, kGreen+2, kOrange+7,
                          kMagenta+2, kCyan+1, kViolet+3, kSpring+5, kGray+2, kAzure+2 };

    // use the SAME marker everywhere: filled circle (ROOT style 20)
  int  markers[NPT];
  std::fill_n(markers, NPT, 20);
    
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

    // 1) Create a canvas for the quadruple overlay
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

    std::vector<TH1D*> hLocalVec(NPT, nullptr);

    // 2) Loop over the 4 pT slices => single overlay
    std::vector<double> bestBvals(NPT, 0.0); // store b-values if we want them

    for(int i=0; i<NPT; i++)
    {
      // Determine Z range for that pT slice
      int zLo = h3->GetZaxis()->FindFixBin(pTedges[i] + 1e-9);
      int zHi = h3->GetZaxis()->FindFixBin(pTedges[i+1] - 1e-9);
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
      // Multi‐start fit approach
      // ---------------------------------------------------------------------
      TF1* fAsinh = new TF1(Form("fAsinh_%s_%d", coordName, i),
                            asinhModel, -0.5, 0.5, 2);
      fAsinh->SetParNames("Norm", "bVal");
      // Parameter limits
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
        hLocal->GetListOfFunctions()->Add(bestFunc); // store it so it doesn't vanish
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

    // Save that overlay
    TString mainOut = outDirPhi + "/" + outPNGname;
    cCoord.SaveAs(mainOut);
    std::cout << "[INFO] Wrote quadruple-overlay " << coordName
              << " => " << mainOut << std::endl;

      // -----------------------------------------------------------------------
      // If it's local φ, also produce the single‐bin overlay of raw vs. corrected
      // -----------------------------------------------------------------------
      if (isPhi)
      {
        // For each of the 4 pT bins
        for(int i=0; i<NPT; i++)
        {
          double ptLo = pT_bins_low[i];
          double ptHi = pT_bins_high[i];

          //----------------------------------------------------------------------
          // 1) Retrieve the RAW Δφ histogram
          //----------------------------------------------------------------------
          TString rawName = Form("h_phi_diff_raw_%.0f_%.0f", ptLo, ptHi);
          TH1F* hRaw = dynamic_cast<TH1F*>(gROOT->FindObject(rawName));
          if(!hRaw)
          {
            std::cerr << "[WARNING] No raw Δφ hist '" << rawName
                      << "' => skip overlay for pT bin i=" << i << std::endl;
            continue;
          }

          //----------------------------------------------------------------------
          // 2) Retrieve the CORRECTED Δφ histogram
          //----------------------------------------------------------------------
          TString corrName = Form("h_phi_diff_corr_%.1f_%.1f", ptLo, ptHi);
          TH1F* hCorr = dynamic_cast<TH1F*>(gROOT->FindObject(corrName));
          if(!hCorr)
          {
            std::cerr << "[WARNING] No corrected Δφ hist '" << corrName
                      << "' => skip overlay for pT bin i=" << i << std::endl;
            continue;
          }

          //----------------------------------------------------------------------
          // 3) Make a canvas & style each histogram
          //----------------------------------------------------------------------
          TString cName  = Form("cDeltaPhiCompare_%.1fto%.1f", ptLo, ptHi);
          TString cTitle = Form("Δ#phi overlay [%.1f < p_{T} < %.1f]", ptLo, ptHi);

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

          //----------------------------------------------------------------------
          // 4) Draw them
          //----------------------------------------------------------------------
          hRaw->Draw("E");
          double maxRaw  = hRaw->GetMaximum();
          double maxCorr = hCorr->GetMaximum();
          double newMax  = TMath::Max(maxRaw, maxCorr)*1.25;
          hRaw->GetYaxis()->SetRangeUser(0, newMax);

          hCorr->Draw("SAME E");

          //----------------------------------------------------------------------
          // 5) Add a legend
          //----------------------------------------------------------------------
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

          //----------------------------------------------------------------------
          // 6) Save
          //----------------------------------------------------------------------
          TString singleOut = Form("%s/DeltaPhiCompare_%.1fto%.1f.png",
                                   outDirPhi.Data(), ptLo, ptHi);
          cSingle.SaveAs(singleOut);

          std::cout << "[INFO] Single-bin Δφ overlay pT=["
                    << ptLo << "," << ptHi
                    << "] => " << singleOut << std::endl;
        } // end loop over i
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
  
  // Everything else remains the same
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

    
  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
