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

// ---------------------------------------------------------------------------
// doLocalPhiFits(...)
//
// This function encapsulates the "Local φ slices in pT bins" code snippet.
// It takes as inputs:
//   - h3 ..................... the 3D histogram pointer
//   - pT_bins_low, pT_bins_high arrays specifying the slice ranges
//   - bResults ............... reference to an open output file for b-values
//   - outDirPhi .............. path to the directory where the output PNG is saved
//
// Everything else (color arrays, legend, eqLatex, etc.) is created locally here.
// ---------------------------------------------------------------------------
void doLocalPhiFits(TH3F* h3,
                    const double pT_bins_low[3],
                    const double pT_bins_high[3],
                    std::ofstream& bResults,
                    const TString& outDirPhi)
{
  // Create canvas
  TCanvas cphi("cphi","Local phi distributions", 900,700);
  cphi.SetTitle("Local #phi_{CG} Distributions and asinh Fits");
  cphi.SetLeftMargin(0.12);
  cphi.SetRightMargin(0.05);
  cphi.SetBottomMargin(0.12);
  cphi.cd();

  int colors[3]  = {kBlack, kRed, kBlue};
  int markers[3] = {20, 21, 22};

  // For labeling the fit equation
  TLatex eqLatex;
  eqLatex.SetNDC(true);
  eqLatex.SetTextFont(42);
  eqLatex.SetTextSize(0.028);
  eqLatex.SetTextAlign(13);

  // The equation string shown on the plot
  TString eqString = "Y(X) = Norm #times #frac{2#it{b}}{#sqrt{1 + 4X^{2} sinh^{2}(1/(2b))}}";

  TLegend leg(0.18,0.68,0.4,0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);

  // Loop over the 3 pT slices
  for(int i=0; i<3; i++)
  {
    int zLo = h3->GetZaxis()->FindFixBin(pT_bins_low[i]+1e-9);
    int zHi = h3->GetZaxis()->FindFixBin(pT_bins_high[i]-1e-9);
    if(zLo<1) zLo=1;
    if(zHi>h3->GetNbinsZ()) zHi=h3->GetNbinsZ();

    // Restrict the 3D histogram to the desired pT slice
    h3->GetXaxis()->SetRange(1, h3->GetNbinsX());
    h3->GetZaxis()->SetRange(zLo, zHi);

    // Project onto the "y" axis => local φ
    TH1D* hLocalPhi = (TH1D*) h3->Project3D("y");
    if(!hLocalPhi)
    {
      std::cerr << "[WARNING] could not project pT bin " << i << std::endl;
      continue;
    }
    hLocalPhi->SetName(Form("hLocalPhi_%.1fto%.1f", pT_bins_low[i], pT_bins_high[i]));
    hLocalPhi->SetTitle(";local #phi_{CG}; scaled counts");

    double integral = hLocalPhi->Integral();
    if(integral > 1e-9) hLocalPhi->Scale(1./integral);

    // Style
    hLocalPhi->SetMarkerStyle(markers[i]);
    hLocalPhi->SetMarkerColor(colors[i]);
    hLocalPhi->SetLineColor(colors[i]);
    hLocalPhi->SetMarkerSize(1.2);

    // Draw on first iteration vs subsequent
    if(i == 0) {
      hLocalPhi->Draw("E");
      hLocalPhi->GetYaxis()->SetRangeUser(0, 0.4);
      hLocalPhi->GetXaxis()->SetTitle("local #phi_{CG} in block");
      hLocalPhi->GetYaxis()->SetTitle("counts (normalized)");
    }
    else {
      hLocalPhi->Draw("same E");
    }

    TF1* fAsinh = new TF1(Form("fAsinh_%d", i), asinhModel, -0.5, 0.5, 2);
    fAsinh->SetParameters(0.2, 0.14); // initial guesses
    fAsinh->SetParNames("Norm", "bVal");

    // Make the fit line match the points:
    fAsinh->SetLineColor(colors[i]);
    fAsinh->SetLineWidth(2);

    hLocalPhi->Fit(fAsinh, "Q", "", -0.5, 0.5);

    double bVal  = fAsinh->GetParameter(1);
    double bValE = fAsinh->GetParError(1);

    // Legend entry
    TString legEntry = Form("p_{T}=[%.0f,%.0f] GeV : b=%.3f #pm %.3f",
                            pT_bins_low[i], pT_bins_high[i], bVal, bValE);
    leg.AddEntry(hLocalPhi, legEntry, "lp");

    // Write b-values to file if open
    if(bResults.is_open()) {
      bResults << pT_bins_low[i] << "  "
               << pT_bins_high[i] << "   "
               << bVal << "\n";
    }
  }

  // Draw legend and eq label
  leg.Draw();
  eqLatex.DrawLatex(0.57, 0.82, eqString);

  // Save to file
  cphi.SaveAs(outDirPhi + "/LocalPhiFits.png");
  std::cout << "[INFO] Wrote local phi slice overlay => "
            << (outDirPhi + "/LocalPhiFits.png") << std::endl;
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
  h3->GetZaxis()->SetRange(1, h3->GetNbinsZ());
  TH2F* h2_block2D = (TH2F*) h3->Project3D("xy");
  h2_block2D->SetName("h2_blockCoord2D");
  h2_block2D->SetTitle("Block #eta vs Block #phi;block #eta;block #phi");

  TCanvas c2d("c2d","Block Eta-Phi",900,700);
  c2d.SetLeftMargin(0.12); c2d.SetRightMargin(0.15);
  c2d.SetBottomMargin(0.12);
  h2_block2D->SetStats(false);
  h2_block2D->Draw("COLZ");
  c2d.SaveAs(outDir2D + "/BlockEtaPhi_2D.png");
  std::cout << "[INFO] Wrote 2D block-eta vs block-phi => " << (outDir2D + "/BlockEtaPhi_2D.png") << std::endl;

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

  doLocalPhiFits(h3, pT_bins_low, pT_bins_high, bResults, outDirPhi);
    
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

  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
