#include <TFile.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

// ---------------------------------------------------------------------
// Toggle whether we only do uncorrected fits (isFirstPass=true)
//   or also overlay corrected (isFirstPass=false).
// ---------------------------------------------------------------------
static bool isFirstPass = true;

// ---------------------------------------------------------------------
// asinhModel(...)
//   Y(X) = Norm * [2*b / sqrt(1 + 4*X^2 * sinh^2(1/(2*b)))]
// ---------------------------------------------------------------------
double asinhModel(double* x, double* par)
{
  double Norm = par[0];
  double b    = par[1];
  if (b < 1e-9) return 1e-12;  // avoid pathological zero-division

  double X   = x[0];
  double sh  = TMath::SinH(1.0/(2.0*b));
  double arg = 1.0 + 4.0 * X*X * sh*sh;

  return Norm * ( (2.0 * b) / TMath::Sqrt(arg) );
}

// ---------------------------------------------------------------------
// LocalPhiOnly2by2() => now 4×2 => 8 bins
//
//   Reads two TH3Fs:
//       "h2_cluster_block_cord_E"           (UNCORRECTED)
//       "h2_cluster_block_cord_E_corrected" (CORRECTED, if isFirstPass==false).
//   Each TH3F has axes:
//       X => blockEta [-0.5..+1.5],
//       Y => blockPhi [-0.5..+1.5],
//       Z => cluster E up to e.g. 30, with 8 slices now:
//         [2..4), [4..6), [6..8), [8..10), [10..12), [12..15), [15..20), [20..30)
//
//   We produce:
//   (A) 2D block‐eta vs block‐phi (uncorr vs corr) in a 4×2 canvas
//   (B) 1D local‐phi distributions => asinh fits for each of the 8 energy slices
//   (C) 1D local‐eta distributions => asinh fits for each of the 8 energy slices
//
//   We write best‐fit b-values for both Phi and Eta to bValues.txt, e.g.:
//       PHI [2.0,4.0)  0.15891
//       ETA [2.0,4.0)  0.21735
//   ... now extended to 8 bins.
//
void LocalPhiOnly2by2()
{
  // 1) Turn off the default stats box
  gStyle->SetOptStat(0);

  // 2) Open the input file
  const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root";
  std::cout << "[DEBUG] Attempting to open file: " << inFile << std::endl;
  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie())
  {
    std::cerr << "[ERROR] Could not open file: " << inFile << std::endl;
    return;
  }
  std::cout << "[INFO] Successfully opened file: " << inFile << std::endl;

  // 3) Retrieve the UNCORRECTED TH3F
  TH3F* hUnc3D = dynamic_cast<TH3F*>(fIn->Get("h2_cluster_block_cord_E"));
  if (!hUnc3D)
  {
    std::cerr << "[ERROR] 'h2_cluster_block_cord_E' not found => cannot proceed.\n";
    fIn->Close();
    return;
  }
  hUnc3D->Sumw2();

  std::cout << "\n[DEBUG] === Uncorrected TH3F Info ===\n";
  hUnc3D->Print();
  std::cout << "=====================================\n" << std::endl;

  // 4) If isFirstPass==false => retrieve the CORRECTED TH3F
  TH3F* hCor3D = nullptr;
  if (!isFirstPass)
  {
    hCor3D = dynamic_cast<TH3F*>(fIn->Get("h2_cluster_block_cord_E_corrected"));
    if (!hCor3D)
    {
      std::cerr << "[WARN] isFirstPass=false but 'h2_cluster_block_cord_E_corrected' not found.\n"
                << "      => no corrected 2D hist available.\n";
    }
    else
    {
      hCor3D->Sumw2();
      std::cout << "[DEBUG] === Corrected TH3F Info ===\n";
      hCor3D->Print();
      std::cout << "===================================\n" << std::endl;
    }
  }

  // 5) The eight energy slices now:
  //    [2..4), [4..6), [6..8), [8..10), [10..12), [12..15), [15..20), [20..30)
  std::vector<std::pair<double,double>> eEdges = {
    {2.0, 4.0},
    {4.0, 6.0},
    {6.0, 8.0},
    {8.0, 10.0},
    {10.0, 12.0},
    {12.0, 15.0},
    {15.0, 20.0},
    {20.0, 30.0}
  };
  const int N_E = eEdges.size(); // =8

  // 6) Prepare a text file for best‐fit b-values
  const char* outDir = "/Users/patsfan753/Desktop";
  gSystem->mkdir(outDir, /*recursive=*/true);

  std::string bValFile = std::string(outDir) + "/bValues.txt";
  std::ofstream bOut(bValFile);
  if (!bOut.is_open())
  {
    std::cerr << "[WARN] Could not open " << bValFile << " => won't save b-values.\n";
  }
  else
  {
    // Label the columns more explicitly now that we have PHI and ETA
    bOut << "# E range    best-b  (PHI or ETA)\n";
  }

  // ===================================================================
  // ========== (A) 2D BlockCoord Histograms (now 4×2) =================
  // ===================================================================
  // If isFirstPass==true => only uncorrected
  // else => also produce a separate corrected canvas

  // Make them a bit larger to fit 8 subpads
  TCanvas c2D_unc("BlockCoord2D_E_unc","Uncorrected 2D block coords vs E",1600,1200);
  c2D_unc.Divide(4,2);  // 8 subpads

  TCanvas* c2D_cor = nullptr;
  if (!isFirstPass)
  {
    c2D_cor = new TCanvas("BlockCoord2D_E_cor","Corrected 2D block coords vs E",1600,1200);
    c2D_cor->Divide(4,2);
  }

  for (int i = 0; i < N_E; i++)
  {
    double eLo = eEdges[i].first;
    double eHi = eEdges[i].second;

    // find Z bins
    int zLoBin = hUnc3D->GetZaxis()->FindBin(eLo + 1e-9);
    int zHiBin = hUnc3D->GetZaxis()->FindBin(eHi - 1e-9);
    if (zLoBin < 1) zLoBin = 1;
    if (zHiBin > hUnc3D->GetNbinsZ()) zHiBin = hUnc3D->GetNbinsZ();

    // ---------- Uncorrected 2D slice ----------
    hUnc3D->GetZaxis()->SetRange(zLoBin, zHiBin);
    TH2D* h2dUnc = (TH2D*) hUnc3D->Project3D("xy");

    c2D_unc.cd(i+1);
    if (h2dUnc)
    {
      h2dUnc->SetName(Form("h2_unc_xy_Ebin%d", i));
      h2dUnc->SetTitle(Form("Uncorr: E=[%.1f,%.1f)  blockEta vs blockPhi", eLo, eHi));
      h2dUnc->GetXaxis()->SetTitle("blockEta_{CG}");
      h2dUnc->GetYaxis()->SetTitle("blockPhi_{CG}");
      h2dUnc->Draw("COLZ");
    }
    else
    {
      std::cerr << "[WARN] Could not project uncorrected 2D => bin i=" << i << std::endl;
    }

    // ---------- Corrected 2D slice if isFirstPass==false ----------
    if (!isFirstPass && hCor3D)
    {
      hCor3D->GetZaxis()->SetRange(zLoBin, zHiBin);
      TH2D* h2dCor = (TH2D*) hCor3D->Project3D("xy");

      c2D_cor->cd(i+1);
      if (h2dCor)
      {
        h2dCor->SetName(Form("h2_cor_xy_Ebin%d", i));
        h2dCor->SetTitle(Form("Corr: E=[%.1f,%.1f)  blockEta vs blockPhi", eLo, eHi));
        h2dCor->GetXaxis()->SetTitle("blockEta_{CG}");
        h2dCor->GetYaxis()->SetTitle("blockPhi_{CG}");
        h2dCor->Draw("COLZ");
      }
      else
      {
        std::cerr << "[WARN] Could not project corrected 2D => bin i=" << i << std::endl;
      }
    }
  }

  // Save uncorrected 2D
  TString out2D_unc = Form("%s/BlockCoord2D_E_unc.png", outDir);
  c2D_unc.SaveAs(out2D_unc);
  std::cout << "[INFO] Wrote uncorrected 2D => " << out2D_unc << std::endl;

  // Save corrected 2D if needed
  if (!isFirstPass && c2D_cor)
  {
    TString out2D_cor = Form("%s/BlockCoord2D_E_cor.png", outDir);
    c2D_cor->SaveAs(out2D_cor);
    std::cout << "[INFO] Wrote corrected 2D => " << out2D_cor << std::endl;

    delete c2D_cor;
    c2D_cor = nullptr;
  }

  // ===================================================================
  // ========== (B) 1D Local‐Phi Projections => asinh fits ==============
  // ===================================================================
  // We'll produce a single TCanvas => "LocalPhiFits_4by2.png" for 8 bins
  TCanvas cFitsPhi("LocalPhiFits_4by2","Local phi fits in 8 energy slices",1600,1200);
  cFitsPhi.Divide(4,2);

  for (int i = 0; i < N_E; i++)
  {
    double eLo = eEdges[i].first;
    double eHi = eEdges[i].second;
    std::cout << "\n[DEBUG] Projecting local phi for E slice [" << eLo << "," << eHi << ")\n";

    // range in uncorrected
    int zLoBin = hUnc3D->GetZaxis()->FindBin(eLo + 1e-9);
    int zHiBin = hUnc3D->GetZaxis()->FindBin(eHi - 1e-9);
    if (zLoBin < 1) zLoBin = 1;
    if (zHiBin > hUnc3D->GetNbinsZ()) zHiBin = hUnc3D->GetNbinsZ();

    cFitsPhi.cd(i+1);

    // uncorrected slice => project "y" => local phi
    hUnc3D->GetZaxis()->SetRange(zLoBin, zHiBin);
    // restore full X so we keep all blockEta
    hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX());

    TH1D* hUncPhi = (TH1D*) hUnc3D->Project3D("y");
    if (!hUncPhi)
    {
      std::cerr << "[WARN] Could not project uncorrected local-phi => i=" << i << std::endl;
      continue;
    }
    hUncPhi->SetName(Form("h1_unc_phi_E%d", i));
    hUncPhi->SetTitle(Form("Local #phi (Uncorr): E=[%.1f,%.1f)", eLo, eHi));
    hUncPhi->GetXaxis()->SetTitle("block #phi_{CG}");
    hUncPhi->GetYaxis()->SetTitle("counts");

    hUncPhi->SetLineColor(kBlack);
    hUncPhi->SetMarkerColor(kBlack);
    hUncPhi->SetMarkerStyle(20);
    hUncPhi->SetMarkerSize(1.2);
    hUncPhi->Draw("E");

    double maxU = hUncPhi->GetMaximum();
    if (maxU < 1e-12) maxU = 1.0;
    hUncPhi->GetYaxis()->SetRangeUser(0, maxU * 1.3);

    // Fit with asinh
    TF1  fAs("fAs", asinhModel, -0.5, 0.5, 2);
    fAs.SetParNames("Norm","bVal");
    fAs.SetParLimits(0, 1e-9, 1e9);
    fAs.SetParLimits(1, 1e-5, 2.0);

    double bestChi2 = 1e15;
    double bestB_phi= 0.0;
    TF1*   bestF    = nullptr;

    // multi-guess
    std::vector<std::pair<double,double>> guesses = {
      {0.1,0.10}, {0.2,0.14}, {0.3,0.20}, {0.5,0.08}
    };
    for (auto& g : guesses)
    {
      fAs.SetParameter(0, g.first);
      fAs.SetParameter(1, g.second);

      TFitResultPtr rr = hUncPhi->Fit(&fAs, "RQN0S", "", -0.5, 0.5);
      if (!rr.Get() || !rr->IsValid()) continue;

      double chi2 = rr->Chi2();
      if (chi2 < bestChi2)
      {
        bestChi2 = chi2;
        if (bestF) delete bestF;
        bestF = (TF1*) fAs.Clone(Form("bestF_E%d_phi", i));
        bestB_phi = bestF->GetParameter(1);
      }
    }

    // draw best fit
    if (bestF)
    {
      bestF->SetLineColor(kBlack);
      bestF->SetLineWidth(2);
      bestF->Draw("SAME");
    }

    // Possibly corrected overlay
    TH1D* hCorPhi = nullptr;
    if (!isFirstPass && hCor3D)
    {
      hCor3D->GetZaxis()->SetRange(zLoBin, zHiBin);
      hCor3D->GetXaxis()->SetRange(1, hCor3D->GetNbinsX());

      hCorPhi = (TH1D*) hCor3D->Project3D("y");
      if (hCorPhi)
      {
        hCorPhi->SetLineColor(kRed);
        hCorPhi->SetMarkerColor(kRed);
        hCorPhi->SetMarkerStyle(21);
        hCorPhi->SetMarkerSize(1.2);
        hCorPhi->SetTitle("");

        hCorPhi->Draw("SAME E");

        double maxC = hCorPhi->GetMaximum();
        double newMax = TMath::Max(maxU, maxC)*1.3;
        hUncPhi->GetYaxis()->SetRangeUser(0, newMax);
      }
    }

    // event counts
    double uncEntries = hUncPhi->GetEntries();
    double corEntries = (hCorPhi ? hCorPhi->GetEntries() : 0.0);

    // ----------------------------------------------------------------
    // OUTPUT ANNOTATIONS
    // ----------------------------------------------------------------
    if (isFirstPass)
    {
      // Just store the # events & b param in TLatex
      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.037);
      lat.SetTextColor(kBlue+3);

      lat.DrawLatex(0.59, 0.85, Form("Uncorr=%.0f evts", uncEntries));
      if (bestB_phi > 1e-9)
      {
        lat.DrawLatex(0.59, 0.78, Form("b=%.3g", bestB_phi));
      }
    }
    else
    {
      // TLegend plus bigger TLatex
      TLegend leg(0.55, 0.65, 0.89, 0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.035);

      leg.AddEntry(hUncPhi, "Uncorrected", "lp");
      if (hCorPhi)
      {
        leg.AddEntry(hCorPhi, "Corrected", "lp");
      }
      if (bestF)
      {
        TString lb = Form("Asinh fit(b=%.3g)", bestB_phi);
        leg.AddEntry(bestF, lb, "l");
      }
      leg.Draw("SAME");

      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.045);
      lat.SetTextColor(kBlue+3);

      double textX = 0.55;
      double textY = 0.83;
      lat.DrawLatex(textX, textY,        Form("Uncorr = %.0f", uncEntries));
      lat.DrawLatex(textX, textY - 0.08, Form("Corr = %.0f",   corEntries));
    }

    // record best b for local phi
    if (bOut.is_open() && bestB_phi > 1e-9)
    {
      bOut << Form("PHI [%.1f,%.1f)  %.5f\n", eLo, eHi, bestB_phi);
    }
  } // end loop over 8 E slices for phi

  // save the 1D PHI fits
  TString outFitsPhi = Form("%s/LocalPhiFits_4by2.png", outDir);
  cFitsPhi.SaveAs(outFitsPhi);
  std::cout << "[INFO] Wrote local phi fits => " << outFitsPhi << std::endl;


  // ===================================================================
  // ========== (C) 1D Local‐Eta Projections => asinh fits =============
  // ===================================================================
  // Another 4×2 canvas for the 8 bins
  TCanvas cFitsEta("LocalEtaFits_4by2","Local eta fits in 8 energy slices",1600,1200);
  cFitsEta.Divide(4,2);

  for (int i = 0; i < N_E; i++)
  {
    double eLo = eEdges[i].first;
    double eHi = eEdges[i].second;
    std::cout << "\n[DEBUG] Projecting local eta for E slice [" << eLo << "," << eHi << ")\n";

    int zLoBin = hUnc3D->GetZaxis()->FindBin(eLo + 1e-9);
    int zHiBin = hUnc3D->GetZaxis()->FindBin(eHi - 1e-9);
    if (zLoBin < 1) zLoBin = 1;
    if (zHiBin > hUnc3D->GetNbinsZ()) zHiBin = hUnc3D->GetNbinsZ();

    cFitsEta.cd(i+1);

    // uncorrected => local-eta => project "x"
    hUnc3D->GetZaxis()->SetRange(zLoBin, zHiBin);
    hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());

    TH1D* hUncEta = (TH1D*) hUnc3D->Project3D("x");
    if (!hUncEta)
    {
      std::cerr << "[WARN] Could not project uncorrected local-eta => i=" << i << std::endl;
      continue;
    }
    hUncEta->SetName(Form("h1_unc_eta_E%d", i));
    hUncEta->SetTitle(Form("Local #eta (Uncorr): E=[%.1f,%.1f)", eLo, eHi));
    hUncEta->GetXaxis()->SetTitle("block #eta_{CG}");
    hUncEta->GetYaxis()->SetTitle("counts");

    hUncEta->SetLineColor(kBlack);
    hUncEta->SetMarkerColor(kBlack);
    hUncEta->SetMarkerStyle(20);
    hUncEta->SetMarkerSize(1.2);
    hUncEta->Draw("E");

    double maxU = hUncEta->GetMaximum();
    if (maxU < 1e-12) maxU = 1.0;
    hUncEta->GetYaxis()->SetRangeUser(0, maxU * 1.3);

    // asinh fit
    TF1  fAs_eta("fAs_eta", asinhModel, -0.5, 0.5, 2);
    fAs_eta.SetParNames("Norm","bVal");
    fAs_eta.SetParLimits(0, 1e-9, 1e9);
    fAs_eta.SetParLimits(1, 1e-5, 2.0);

    double bestChi2 = 1e15;
    double bestB_eta= 0.0;
    TF1*   bestFeta = nullptr;

    std::vector<std::pair<double,double>> guesses_eta = {
      {0.1,0.10}, {0.2,0.14}, {0.3,0.20}, {0.5,0.08}
    };
    for (auto& g : guesses_eta)
    {
      fAs_eta.SetParameter(0, g.first);
      fAs_eta.SetParameter(1, g.second);

      TFitResultPtr rr = hUncEta->Fit(&fAs_eta, "RQN0S", "", -0.5, 0.5);
      if (!rr.Get() || !rr->IsValid()) continue;

      double chi2 = rr->Chi2();
      if (chi2 < bestChi2)
      {
        bestChi2 = chi2;
        if (bestFeta) delete bestFeta;
        bestFeta = (TF1*) fAs_eta.Clone(Form("bestF_E%d_eta", i));
        bestB_eta = bestFeta->GetParameter(1);
      }
    }

    if (bestFeta)
    {
      bestFeta->SetLineColor(kBlack);
      bestFeta->SetLineWidth(2);
      bestFeta->Draw("SAME");
    }

    // corrected overlay
    TH1D* hCorEta = nullptr;
    if (!isFirstPass && hCor3D)
    {
      hCor3D->GetZaxis()->SetRange(zLoBin, zHiBin);
      hCor3D->GetYaxis()->SetRange(1, hCor3D->GetNbinsY());

      hCorEta = (TH1D*) hCor3D->Project3D("x");
      if (hCorEta)
      {
        hCorEta->SetLineColor(kRed);
        hCorEta->SetMarkerColor(kRed);
        hCorEta->SetMarkerStyle(21);
        hCorEta->SetMarkerSize(1.2);
        hCorEta->SetTitle("");

        hCorEta->Draw("SAME E");

        double maxC = hCorEta->GetMaximum();
        double newMax = TMath::Max(maxU, maxC)*1.3;
        hUncEta->GetYaxis()->SetRangeUser(0, newMax);
      }
    }

    double uncEntries_eta = hUncEta->GetEntries();
    double corEntries_eta = (hCorEta ? hCorEta->GetEntries() : 0.0);

    if (isFirstPass)
    {
      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.037);
      lat.SetTextColor(kBlue+3);

      lat.DrawLatex(0.59, 0.85, Form("Uncorr=%.0f evts", uncEntries_eta));
      if (bestB_eta > 1e-9)
      {
        lat.DrawLatex(0.59, 0.78, Form("b=%.3g", bestB_eta));
      }
    }
    else
    {
      TLegend leg(0.55, 0.65, 0.89, 0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.035);

      leg.AddEntry(hUncEta, "Uncorrected", "lp");
      if (hCorEta)
      {
        leg.AddEntry(hCorEta, "Corrected", "lp");
      }
      if (bestFeta)
      {
        TString lb = Form("Asinh fit(b=%.3g)", bestB_eta);
        leg.AddEntry(bestFeta, lb, "l");
      }
      leg.Draw("SAME");

      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.045);
      lat.SetTextColor(kBlue+3);

      double textX = 0.55;
      double textY = 0.83;
      lat.DrawLatex(textX, textY,        Form("Uncorr = %.0f", uncEntries_eta));
      lat.DrawLatex(textX, textY - 0.08, Form("Corr = %.0f",   corEntries_eta));
    }

    if (bOut.is_open() && bestB_eta > 1e-9)
    {
      bOut << Form("ETA [%.1f,%.1f)  %.5f\n", eLo, eHi, bestB_eta);
    }
  } // end loop over 8 E slices for eta

  // save the 1D ETA fits
  TString outFitsEta = Form("%s/LocalEtaFits_4by2.png", outDir);
  cFitsEta.SaveAs(outFitsEta);
  std::cout << "[INFO] Wrote local eta fits => " << outFitsEta << std::endl;

  // finalize text file
  if (bOut.is_open())
  {
    bOut.close();
    std::cout << "[INFO] Wrote b-values => " << bValFile << "\n";
    std::cout << "[INFO] (Check file for lines labeled PHI vs ETA.)\n";
  }

  // close file
  fIn->Close();
  delete fIn;

  std::cout << "[INFO] LocalPhiOnly2by2() completed successfully.\n";
}
