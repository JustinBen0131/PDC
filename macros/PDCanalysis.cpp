#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include <TFile.h>
#include <algorithm>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TProfile2D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <functional>
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
#include <cmath>
#include <iomanip>
#include <map>

/**
 * \brief Toggle: If isFirstPass=false, we overlay corrected histograms.
 *                If true, we do uncorrected only.
 */
static bool isFirstPass = false;

struct BRes {
  double val{std::nan("")};   // best-fit b
  double err{0.0};            // 1 σ statistical error
};

/**
 * \brief 8 energy bins: E_edges[] = {2,4,6,8,10,12,15,20,30}
 *        => N_E=8
 */
constexpr double E_edges[] = {2,4,6,8,10,12,15,20,30};
constexpr int    N_E = (sizeof(E_edges)/sizeof(E_edges[0])) - 1;

// Global |z_vtx| bin edges used everywhere (absolute-|z| and ±z views)
static constexpr std::array<float,19> vzEdge = {
     0.f,  5.f, 10.f, 15.f, 20.f, 25.f, 30.f,
    40.f, 50.f, 60.f, 70.f, 80.f, 90.f,
   100.f,110.f,120.f,130.f,140.f,150.f
};
static constexpr int N_VzBins = static_cast<int>(vzEdge.size()) - 1;

/* ────────────────────────────────────────────────────────────────── */
/*  Global, compact palette – order matches the five energy slices   */
/*  ( BLUE , GREEN ,  RED ,  PURPLE , BLACK )                        */
/* ────────────────────────────────────────────────────────────────── */
static const int kAshPalette[5] =
{
    kBlue  + 1,      // deep blue
    kGreen + 2,      // vivid green
    kRed   + 1,      // bright red
    kMagenta + 2,    // purple
    kBlack           // black
};

//--------------------------------------------------------------------
//  Global enum + label helper
//--------------------------------------------------------------------
enum class EBinningMode { kRange, kDiscrete };
static EBinningMode binMode = EBinningMode::kRange;

/** Return the exact label that the producer used for E-slice *iE* */
inline std::string
makeLabel(int iE, EBinningMode mode, const double* E_edges /*N_E+1*/)
{
    if (mode == EBinningMode::kRange)
        return Form("%.0f_%.0f", E_edges[iE], E_edges[iE+1]); // e.g. "2_4"
    else
        return Form("E%.0f",      E_edges[iE]);                // e.g. "E6"
}

/**
 * \brief Enhanced Gaussian fit:
 *   - use quartiles for initial range
 *   - fallback to ±2 RMS if quartiles invalid
 *   - multiple guesses for amplitude, mean, sigma
 *   - debug prints
 */
double coreGaussianSigma(TH1* h, const TString& pngSavePath)
{
  std::cout << "\n[VERBOSE] coreGaussianSigma() called for hist='"
            << (h ? h->GetName() : "NULL") << "', pngSave='" << pngSavePath << "'\n";

  if(!h)
  {
    std::cout << "  [DEBUG] => Hist pointer is null => returning 0.\n";
    return 0.;
  }
  double nEntries = h->GetEntries();
  if(nEntries < 20)
  {
    std::cout << "  [DEBUG] => Hist '" << h->GetName()
              << "' has <20 entries => returning 0.\n";
    return 0.;
  }

  double inRange = h->Integral();
  if(inRange <= 0.)
  {
    std::cout << "  [DEBUG] => Hist '"<<h->GetName()
              <<"' has 0 total integral => returning 0.\n";
    return 0.;
  }

  // 1) Quartiles
  double q[3], probs[3] = {0.25, 0.50, 0.75};
  h->GetQuantiles(3, q, probs);
  double xMin = q[0];
  double xMax = q[2];

  std::cout << "  [DEBUG] => quartiles for '"<<h->GetName()<<"' => "
            << "Q25="<<q[0]<<", Q50="<<q[1]<<", Q75="<<q[2]
            <<". => Proposed fit range=["<<xMin<<","<<xMax<<"]\n";

  // fallback if quartiles are invalid
  if(xMin >= xMax)
  {
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    xMin = mean - 2.*rms;
    xMax = mean + 2.*rms;
    std::cout << "  [WARN] => quartiles invalid => using fallback ±2 RMS => mean="
              << mean <<", rms="<<rms<< " => range=["<<xMin<<","<<xMax<<"]\n";
  }

  // 2) Build TF1 in [xMin,xMax]
  TF1 fitG("fitG","gaus", xMin, xMax);

  // multiple guesses
  std::vector<std::tuple<double,double,double>> initGuesses = {
    { h->GetMaximum(),           h->GetMean(),        h->GetRMS()/2. },
    { h->GetMaximum()*1.2,       h->GetMean(),        h->GetRMS()    },
    { h->GetMaximum()*0.8,       h->GetMean()+0.3*h->GetRMS(), 0.5*h->GetRMS() },
    { h->GetMaximum()*1.5,       h->GetMean()-0.3*h->GetRMS(), 1.5*h->GetRMS()}
  };

  double bestSigma = 0.;
  bool fitOK       = false;
  double bestChi2  = 1e15;

  // we'll do a "S" to store fit result properly
  const char* fitOpts = "RQNS"; // R => use range, Q => quiet, N => no draw, S => store result

  for(size_t i=0; i<initGuesses.size(); i++)
  {
    double A0 = std::get<0>(initGuesses[i]);
    double M0 = std::get<1>(initGuesses[i]);
    double S0 = std::fabs(std::get<2>(initGuesses[i]));
    if(S0<1e-6) S0=0.5; // avoid zero

    fitG.SetParameter(0,A0);
    fitG.SetParameter(1,M0);
    fitG.SetParameter(2,S0);

    std::cout <<"    [DEBUG] Attempt #"<<i<<" => A0="<<A0<<", M0="<<M0<<", S0="<<S0
              <<" => Fit range=["<<xMin<<","<<xMax<<"]\n";

    TFitResultPtr rp = h->Fit(&fitG, fitOpts, "", xMin, xMax);
    int fitStatus = rp; // cast to int
    bool valid    = (rp.Get()!=nullptr && rp->IsValid());
    double chi2   = (valid ? rp->Chi2() : 1e15);

    std::cout <<"        => Fit result: fitStatus="<<fitStatus<<", chi2="<<chi2
              <<", valid="<<(valid?"Yes":"No")<<"\n";

    if(fitStatus==0 && valid && chi2<bestChi2)
    {
      bestChi2   = chi2;
      bestSigma  = fitG.GetParameter(2);
      fitOK      = true;
      std::cout <<"        => new bestSigma="<<bestSigma<<"\n";
    }
  }

  if(!fitOK)
  {
    bestSigma = h->GetRMS();
    std::cout <<"  [WARN] => all fits failed => returning RMS="<<bestSigma<<"\n";
  }
  else
  {
    std::cout <<"  [INFO] => final bestSigma="<<bestSigma
              <<", bestChi2="<<bestChi2<<"\n";
  }

  // 3) optional debug plot
  if(!pngSavePath.IsNull())
  {
    std::cout <<"  [DEBUG] => making debug canvas => '"<<pngSavePath<<"'\n";
    TCanvas ctmp("ctmp","coreGaussianSigma debug canvas",600,500);
    ctmp.cd();

    h->Draw("E");
    if(fitOK) fitG.Draw("SAME");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.04);
    lat.SetTextColor(kRed);

    double fitMean = fitG.GetParameter(1);
    lat.DrawLatex(0.58, 0.82, Form("#mu=%.3f", fitMean));
    lat.DrawLatex(0.58, 0.74, Form("#sigma_{RMS}=%.3f", bestSigma));

    ctmp.SaveAs(pngSavePath);
    std::cout <<"  [INFO] => Wrote debug PNG => "<<pngSavePath<<"\n";
  }

  return bestSigma;
}

/**
 * \brief rawRMS
 * Basic RMS approach, with debug prints
 */
double rawRMS(TH1* h, const TString& debugPNG)
{
  std::cout << "\n[VERBOSE] rawRMS() called for hist='"
            << (h ? h->GetName() : "NULL")
            <<"', debugPNG="<<debugPNG<<"\n";

  if(!h)
  {
    std::cout<<"  => hist pointer null => return 0.\n";
    return 0.;
  }
  double nE = h->GetEntries();
  if(nE<2)
  {
    std::cout<<"  => hist '"<<h->GetName()<<"' has <2 entries => return 0.\n";
    return 0.;
  }
  double val = h->GetRMS();
  std::cout<<"  => hist '"<<h->GetName()<<"' => RMS="<<val<<"\n";

  // (Optional) if you want a debug canvas for RMS only:
  /*
  if(!debugPNG.IsNull())
  {
    TCanvas ctmp("ctmpRMS","RMS debug",600,500);
    ctmp.cd();
    h->Draw("E");
    ctmp.SaveAs(debugPNG);
  }
  */

  return val;
}

/** \brief Container for scanning results across multiple E-slices */
struct ScanResults
{
  std::vector<TGraph*> tgVec;      ///< One TGraph per E bin
  std::vector<double> bestParam;   ///< best b or w0
  std::vector<double> bestSigma;   ///< best sigma
  double minY=DBL_MAX, maxY=-DBL_MAX;
  int totalHist=0, missingHist=0, zeroEntry=0, usedHist=0;
};

/**
 * \brief doAshScan
 *  - loops over bScan in [0..1.60] cm, steps of 0.05 cm (user-provided)
 *  - for each E-slice => tries to retrieve h_dx_ash_bXYZ_Ei
 *  - calls sigmaFunc on that histogram
 *  - saves TGraph of sigma vs bVal (in cm)
 */
ScanResults doAshScan(
    TFile* f,
    int N_E,
    const double* E_edges,
    EBinningMode mode,
    const std::vector<double>& bScan_cm,    // in actual cm steps (like 0.00..1.60)
    double cellSize,                        // 5.55 cm
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir,
    bool  skipLowest3 = false
)
{
  const int iE_start = skipLowest3 ? 3 : 0;
  ScanResults results;
  results.tgVec.resize(N_E);

  std::cout << "\n[DEBUG] doAshScan => ENTER. suffix='"<<suffix<<"', #Ebins="<<N_E
            <<", #bScan="<<bScan_cm.size()<<", cellSize="<<cellSize<<"\n";

  for (int iE = iE_start; iE < N_E; ++iE)
  {
    double eLo = E_edges[iE];
    double eHi = E_edges[iE+1];
    TGraph* g  = new TGraph;
    g->SetName( Form("gAsh_%s_E%d", suffix.Data(), iE) );

    std::cout <<"\n[DEBUG]  => E-bin="<<iE<<" => E=["<<eLo<<","<<eHi<<")\n"
              <<"           #bScan="<<bScan_cm.size()<<" steps.\n";

    int nPointsUsed = 0;
    for(double bVal : bScan_cm)             // bVal *already* is in tower units
    {
      results.totalHist++;

      // build histogram name
      const std::string lab = makeLabel(iE, binMode, E_edges);
      TString hName = Form("h_dx_ash_b%.4f_%s", bVal, lab.c_str());
        
      TH1* hptr = dynamic_cast<TH1*>( f->Get(hName) );

      std::cout<<"  [DEBUG] => bVal(tower)="<<bVal<<", => bVal="<<bVal
               <<" => histName='"<<hName<<"'\n";

      if(!hptr)
      {
        results.missingHist++;
        std::cout<<"     [WARN] => missing histogram => skip.\n";
        continue;
      }

      double ent    = hptr->GetEntries();
      double inrange= hptr->Integral();
      std::cout<<"     [INFO] => found '"<<hName<<"' => entries="<<ent
               <<", integral="<<inrange<<"\n";

      if(inrange<=0 || ent<=0)
      {
        results.zeroEntry++;
        std::cout<<"       => no in-range counts => skip.\n";
        continue;
      }
      results.usedHist++;

      // build diagName if you want to save a PNG
      TString diagName = Form("%s/%s_%s.png", histOutDir.Data(), hName.Data(), suffix.Data());

      // compute sigma
      double sVal = sigmaFunc(hptr, diagName);
      std::cout<<"       => sigmaFunc => sVal="<<sVal<<"\n";

      // fill TGraph if sVal>0
      if(sVal>0)
      {
        g->SetPoint( g->GetN(), bVal, sVal );
        nPointsUsed++;
        if(sVal<results.minY) results.minY=sVal;
        if(sVal>results.maxY) results.maxY=sVal;

        // update best param
        if( (int)results.bestSigma.size() < N_E )
        {
          results.bestSigma.resize(N_E, DBL_MAX);
          results.bestParam.resize(N_E, 0.);
        }
        if(sVal < results.bestSigma[iE])
        {
          results.bestSigma[iE] = sVal;
          results.bestParam[iE] = bVal;
          std::cout<<"         => new best sigma="<<sVal<<" at bVal="<<bVal<<"\n";
        }
      }
      else
      {
        std::cout<<"       => sVal<=0 => skip.\n";
      }
    } // end bScan

    // style
    int idx = (iE - iE_start) % 5;          // wrap every 5 slices
    int col = kAshPalette[idx];

    g->SetMarkerStyle(20 + idx);            // circles, squares, triangles, …
    g->SetMarkerColor(col);
    g->SetLineColor  (col);
    g->SetMarkerSize(1.2);

    results.tgVec[iE] = g;
    std::cout<<"   => TGraph for Ebin="<<iE<<" has n="<<nPointsUsed<<" points used.\n";
  }

  if(results.minY==DBL_MAX)
  {
    results.minY = 0.;
    results.maxY = 1.;
    std::cout<<"[WARN] doAshScan => no valid sigma => forcing [0..1] range.\n";
  }

  std::cout<<"\n[DEBUG] doAshScan => done. Summary:\n"
           <<"   totalHist="<<results.totalHist<<"\n"
           <<"   missing  ="<<results.missingHist<<"\n"
           <<"   zeroEntry="<<results.zeroEntry<<"\n"
           <<"   usedHist ="<<results.usedHist<<"\n"
           <<"   minY="<<results.minY<<", maxY="<<results.maxY<<"\n";

  return results;
}

/**
 * \brief doLogScan
 *   - w0 in [1.5..7.0], step=0.1
 *   - retrieve h_dx_log_w0XX_Ei
 *   - call sigmaFunc
 *   - fill TGraph => sigma vs w0
 */
ScanResults doLogScan(
    TFile* f,
    int N_E,
    const double* E_edges,
    EBinningMode mode,
    const std::vector<double>& w0Scan,
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir,
    bool  skipLowest3 = false
)
{
  const int iE_start = skipLowest3 ? 3 : 0;
  ScanResults results;
  results.tgVec.resize(N_E);
  results.minY=DBL_MAX;
  results.maxY=-DBL_MAX;

  std::cout << "\n[DEBUG] doAshScan => ENTER. suffix='"<<suffix<<"', #Ebins="
              <<N_E<<", skipping lowest 3: "<<std::boolalpha<<skipLowest3<<"\n";


  for (int iE = iE_start; iE < N_E; ++iE)
  {
    double eLo = E_edges[iE];
    double eHi = E_edges[iE+1];

    TGraph* g = new TGraph;
    g->SetName(Form("gLog_%s_E%d", suffix.Data(), iE));

    std::cout<<"\n[DEBUG] => Ebin="<<iE<<" => E=["<<eLo<<","<<eHi<<") => scanning w0...\n";

    for(double wraw : w0Scan)
    {
      results.totalHist++;

      double wVal = std::round(wraw*100.)/100.;
      const std::string lab = makeLabel(iE, binMode, E_edges);
      TString hN  = Form("h_dx_log_w0%.2f_%s", wVal, lab.c_str());

      TH1* hptr = dynamic_cast<TH1*>( f->Get(hN) );
      if(!hptr)
      {
        results.missingHist++;
        std::cout<<"   [WARN] => missing '"<<hN<<"'\n";
        continue;
      }
      double ent = hptr->GetEntries();
      double ing = hptr->Integral();
      std::cout<<"   [INFO] => histName='"<<hN<<"' => entries="<<ent<<", integral="<<ing<<"\n";

      if(ent<=0 || ing<=0)
      {
        results.zeroEntry++;
        std::cout<<"     => skip => no data.\n";
        continue;
      }
      results.usedHist++;

      TString diagName = Form("%s/%s_%s.png", histOutDir.Data(), hN.Data(), suffix.Data());

      double sVal = sigmaFunc(hptr, diagName);
      std::cout<<"     => sVal="<<sVal<<" => wVal="<<wVal<<"\n";

      if(sVal>0)
      {
        g->SetPoint(g->GetN(), wVal, sVal);
        if(sVal<results.minY) results.minY=sVal;
        if(sVal>results.maxY) results.maxY=sVal;

        if( (int)results.bestSigma.size()<N_E )
        {
          results.bestSigma.resize(N_E, DBL_MAX);
          results.bestParam.resize(N_E, 0.);
        }
        if(sVal<results.bestSigma[iE])
        {
          results.bestSigma[iE]=sVal;
          results.bestParam[iE]=wVal;
          std::cout<<"       => new bestSigma="<<sVal<<" at w0="<<wVal<<"\n";
        }
      }
      else
      {
        std::cout<<"     => sVal<=0 => skip.\n";
      }
    } // end w0Scan

    // style
    int idx = (iE - iE_start) % 5;          // wrap every 5 slices
    int col = kAshPalette[idx];

    g->SetMarkerStyle(20 + idx);            // circles, squares, triangles, …
    g->SetMarkerColor(col);
    g->SetLineColor  (col);
    g->SetMarkerSize(1.2);
    results.tgVec[iE] = g;
  }

  if(results.minY==DBL_MAX)
  {
    results.minY=0.;
    results.maxY=1.;
  }
  return results;
}


/* ==================================================================== */
/*  Side‑by‑side Ash‑vs‑Log overlay                                     */
/* ==================================================================== */
void drawAshLogSideBySide(
    const ScanResults& ashRes,
    const ScanResults& logRes,
    const char*        methodName,
    const double*      E_edges,
    int                N_E,
    const TString&     baseDir,
    bool               skipLowest3 /* = false */
)
{
  /* --------------------------------------------------------------- */
  const int iE_start = skipLowest3 ? 3 : 0;          // << core switch
  std::cout << "\n[DEBUG] drawAshLogSideBySide => method='" << methodName
            << "', skipLowest3=" << std::boolalpha << skipLowest3 << '\n';

  /* --------------------------------------------------------------- */
  TString cName = Form("cSideBySide_%s", methodName);
  TCanvas cSide(cName, cName, 1600, 600);
  cSide.Divide(2, 1);

  /* ============= LEFT PAD : ASH ============= */
  cSide.cd(1);
  gPad->SetLeftMargin  (0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  std::cout << "[DEBUG]  Ash  => minY=" << ashRes.minY
            << "  maxY=" << ashRes.maxY << '\n';

  /* ---- pick first kept slice as reference frame ---- */
  TGraph* gAshRef = nullptr;
  for (int iE = iE_start; iE < N_E && !gAshRef; ++iE)
    if (ashRes.tgVec[iE] && ashRes.tgVec[iE]->GetN())
      gAshRef = ashRes.tgVec[iE];

  if (gAshRef)
  {
    gAshRef->Draw("ALP");
    gAshRef->GetXaxis()->SetTitle("b  (local tower units)");
    gAshRef->GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (local tower units)");

    const double yLo = std::max(0.0, ashRes.minY * 0.90);
    const double yHi =               ashRes.maxY * 1.10;
    gAshRef->GetYaxis()->SetRangeUser(yLo, yHi);

    for (int iE = iE_start; iE < N_E; ++iE)
    {
      TGraph* g = ashRes.tgVec[iE];
      if (g && g->GetN()) g->Draw("LP SAME");
    }
  }
  else
  {
    std::cout << "[WARN]  no Ash graphs to draw – left pad left blank.\n";
  }

  /* ---- legend (Ash) ---- */
  TLegend legA(0.13, 0.65, 0.48, 0.88);
  legA.SetBorderSize(0);
  legA.SetFillStyle(0);

  for (int iE = iE_start; iE < N_E; ++iE)
  {
    TGraph* g = ashRes.tgVec[iE];
    const double bBest = (iE < (int)ashRes.bestParam .size()) ? ashRes.bestParam [iE] : 0.;
    const double sBest = (iE < (int)ashRes.bestSigma.size()) ? ashRes.bestSigma[iE] : 0.;

    if (!g || g->GetN() == 0)
    {
        legA.AddEntry((TObject*)nullptr,
          Form("%.0f #leq E < %.0f GeV  — no data",
               E_edges[iE], E_edges[iE+1]),
          "");
      continue;
    }
      legA.AddEntry(g,
        Form("%.0f #leq E < %.0f GeV  (b=%.2f,  #sigma_{RMS}=%.3f)",
             E_edges[iE], E_edges[iE+1], bBest, sBest),
        "lp");

  }
  legA.Draw();

  TLatex tA;  tA.SetNDC();  tA.SetTextSize(0.035);
  tA.DrawLatex(0.12, 0.93, Form("Ash scan [%s]", methodName));

  /* ============= RIGHT PAD : LOG ============= */
  cSide.cd(2);
  gPad->SetLeftMargin  (0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  std::cout << "[DEBUG]  Log  => minY=" << logRes.minY
            << "  maxY=" << logRes.maxY << '\n';

  TGraph* gLogRef = nullptr;
  for (int iE = iE_start; iE < N_E && !gLogRef; ++iE)
    if (logRes.tgVec[iE] && logRes.tgVec[iE]->GetN())
      gLogRef = logRes.tgVec[iE];

  if (gLogRef)
  {
    gLogRef->Draw("ALP");
    gLogRef->GetXaxis()->SetRangeUser( 3.5, 5.5 );
    gLogRef->GetXaxis()->SetTitle("w_{0}");
    gLogRef->GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (local tower units)");

    const double yLo = std::max(0.0, logRes.minY * 0.90);
    const double yHi =               logRes.maxY * 1.10;
    gLogRef->GetYaxis()->SetRangeUser(yLo, yHi);

    for (int iE = iE_start; iE < N_E; ++iE)
    {
      TGraph* g = logRes.tgVec[iE];
      if (g && g->GetN()) g->Draw("LP SAME");
    }
  }
  else
  {
    std::cout << "[WARN]  no Log graphs to draw – right pad left blank.\n";
  }

  /* ---- legend (Log) ---- */
  TLegend legB(0.60, 0.65, 0.85, 0.88);
  legB.SetBorderSize(0);
  legB.SetFillStyle(0);

  for (int iE = iE_start; iE < N_E; ++iE)
  {
    TGraph* g = logRes.tgVec[iE];
    const double wBest = (iE < (int)logRes.bestParam .size()) ? logRes.bestParam [iE] : 0.;
    const double sBest = (iE < (int)logRes.bestSigma.size()) ? logRes.bestSigma[iE] : 0.;

    if (!g || g->GetN() == 0)
    {
        legB.AddEntry((TObject*)nullptr,
          Form("%.0f #leq E < %.0f GeV  — no data",
               E_edges[iE], E_edges[iE+1]),
          "");
      continue;
    }
      legB.AddEntry(g,
        Form("%.0f #leq E < %.0f GeV  (w_{0}=%.2f,  #sigma_{RMS}=%.3f)",
             E_edges[iE], E_edges[iE+1], wBest, sBest),
        "lp");
  }
  legB.Draw();

  TLatex tB;  tB.SetNDC();  tB.SetTextSize(0.040);
  tB.DrawLatex(0.15, 0.93, Form("Log scan [%s]", methodName));

  /* --------------------------------------------------------------- */
  TString out = Form("%s/SideBySide_%s.png", baseDir.Data(), methodName);
  cSide.SaveAs(out);
  std::cout << "[INFO]  wrote side‑by‑side figure → " << out << '\n';

  /* ----------------------------------------------------------------
   *  ---  Publication‑style single‑panel plots  -------------------
   *  (exactly as before, but start loop from iE_start)
   * ---------------------------------------------------------------- */
  auto makeSinglePanel =
      [&](const ScanResults& R, const char* tag,
          const char* xTit, const char* yTit)
  {
    /* pick reference */
    TGraph* gRef = nullptr;
    for (int iE = iE_start; iE < N_E && !gRef; ++iE)
      if (R.tgVec[iE] && R.tgVec[iE]->GetN()) gRef = R.tgVec[iE];

    if (!gRef) return;                     // nothing to draw

    TString cname = Form("c%sScan_%s", tag, methodName);
    TCanvas c(cname, cname, 900, 650);
    c.SetLeftMargin(0.13);
    c.SetBottomMargin(0.18);

    gRef->Draw("ALP");
    gRef->GetXaxis()->SetTitle(xTit);
    gRef->GetYaxis()->SetTitle(yTit);
    /* --- dynamic Y‑range from visible points only -------------------- */
    const bool isLog = (std::strcmp(tag,"Log")==0);
    double yMin = 1e30 , yMax = -1e30 , xCut = (isLog ? 5.5 : 1e99);
    for (int i=iE_start;i<N_E;++i) if (R.tgVec[i]&&R.tgVec[i]->GetN())
    for (int j=0;j<R.tgVec[i]->GetN();++j){ double x,y; R.tgVec[i]->GetPoint(j,x,y);
        if (x<=xCut){ yMin = std::min(yMin,y); yMax = std::max(yMax,y);} }
    const double pad = 0.05*(yMax-yMin);          // ±5 % padding
    gRef->GetYaxis()->SetRangeUser(std::max(0.0,yMin-pad), yMax+pad);

    for (int iE = iE_start; iE < N_E; ++iE)
      if (R.tgVec[iE] && R.tgVec[iE]->GetN())
        R.tgVec[iE]->Draw("LP SAME");

    /* legend */
    TLegend L(0.63, 0.70, 0.90, 0.92);
    L.SetBorderSize(0);  L.SetFillStyle(0);  L.SetTextSize(0.022);
    for (int iE = iE_start; iE < N_E; ++iE)
      if (R.tgVec[iE] && R.tgVec[iE]->GetN())
          L.AddEntry(R.tgVec[iE],
            Form("%.0f #leq E < %.0f GeV", E_edges[iE], E_edges[iE+1]),
            "lp");
    L.Draw();

    TString outSingle = Form("%s/%sScan_%s.png",
                             baseDir.Data(), tag, methodName);
    c.SaveAs(outSingle);
    std::cout << "        wrote " << outSingle << '\n';
  };

  makeSinglePanel(ashRes, "Ash", "b  (local tower units)",
                                "#sigma_{RMS}_{x}  (local tower units)");
  makeSinglePanel(logRes, "Log", "w_{0}",
                                "#sigma_{RMS}_{x}  (local tower units)");
}


/* ==================================================================== */
/*  PHENIX‑style σₓ(E) summary – ASH‑only version                       */
/* ==================================================================== */
void makePhenixLikeSummary( const ScanResults& ashRes,
                            [[maybe_unused]] const ScanResults& logRes,
                            const double*      E_edges,
                            int                N_E,
                            double             cellSize_cm,   // 5.55
                            const TString&     outDir,
                            bool               skipLowest3 = false )
{
  std::cout << "\n[PHX‑Σ] ENTER  (cell=" << cellSize_cm
            << " cm;  skipLowest3=" << std::boolalpha << skipLowest3 << ")\n";

  const int iE_start = skipLowest3 ? 3 : 0;

  /* ------------------------------------------------------------------ */
  /* 1) Fill TGraph with Ash minima – convert ⇒ mm                      */
  /* ------------------------------------------------------------------ */
  TGraph gAsh;  gAsh.SetName("gAshBest");

  int    nGood = 0;
  double xMin  =  1e30, xMax = -1e30;

  for (int iE = iE_start; iE < N_E; ++iE)
  {
      const double eCtr = 0.5*(E_edges[iE] + E_edges[iE+1]);      // GeV

      const bool   has = (iE < (int)ashRes.bestSigma.size() &&
                          std::isfinite(ashRes.bestSigma[iE]));
      const double sig_mm = has ? ashRes.bestSigma[iE]*cellSize_cm*10.0 : NAN;

      std::cout << "  slice " << std::setw(2) << iE
                << "  centre=" << std::setw(6) << eCtr << " GeV"
                << "   σ=" << (has?std::to_string(sig_mm):"n/a") << "  (mm)\n";

      if (!has) continue;

      gAsh.SetPoint(gAsh.GetN(), eCtr, sig_mm);
      ++nGood;  xMin = std::min(xMin,eCtr);  xMax = std::max(xMax,eCtr);
  }

  if (nGood < 3)
  {
      std::cerr << "[PHX‑Σ][FATAL] need ≥3 Ash points, have "
                << nGood << " – abort.\n";
      return;
  }

  /* ------------------------------------------------------------------ */
  /* 2) Aesthetics                                                      */
  /* ------------------------------------------------------------------ */
  gAsh.SetMarkerStyle(24);  gAsh.SetMarkerSize(1.4);
  gAsh.SetMarkerColor(kBlack);  gAsh.SetLineColor(kBlack);

  /* ------------------------------------------------------------------ */
  /* 3)  Fit  p0 + p1 / √E                                              */
  /* ------------------------------------------------------------------ */
  TF1 fFit("fFit","[0]+[1]/TMath::Sqrt(x)", xMin*0.9 , xMax*1.1);
  fFit.SetParameters(1.0, 6.0);

  std::cout << "[PHX‑Σ] fit range  "
            << fFit.GetXmin() << " – " << fFit.GetXmax() << " GeV\n";
  TFitResultPtr fr = gAsh.Fit(&fFit,"QRNS");
  const bool fitOK = (fr && fr->IsValid());

  if (!fitOK)
      std::cerr << "[PHX‑Σ][WARN] fit failed – curve will still be drawn.\n";

  /* ------------------------------------------------------------------ */
  /* 4) Canvas & axes                                                   */
  /* ------------------------------------------------------------------ */
  TCanvas c("cPhenixAshOnly","PHENIX σx(E) – Ash only",600,600);
  c.SetLeftMargin(0.17);  c.SetBottomMargin(0.15);

  gAsh.Draw("AP");
  gAsh.GetXaxis()->SetTitle("E  (GeV)");
  gAsh.GetYaxis()->SetTitle("#sigma_{RMS}_{x}  (mm)");
  gAsh.GetXaxis()->SetLimits(xMin*0.9 , xMax*1.1);

  /* Y–range with padding                                               */
  double yMin=1e30,yMax=-1e30, x, y;
  for (int i=0;i<gAsh.GetN();++i){ gAsh.GetPoint(i,x,y);
                                   yMin=std::min(yMin,y);
                                   yMax=std::max(yMax,y);}
  const double pad = 0.15*(yMax-yMin);
  gAsh.GetYaxis()->SetRangeUser(std::max(0.0,yMin-pad), yMax+pad);

  if (fitOK) { fFit.SetLineStyle(2);  fFit.Draw("SAME"); }

  /* ------------------------------------------------------------------ */
  /* 5) Legend & fit text                                               */
  /* ------------------------------------------------------------------ */
  {
      TLegend L(0.55,0.75,0.88,0.88);
      L.SetBorderSize(0);  L.SetFillStyle(0);
      L.AddEntry(&gAsh,"Ash","p");
      if (fitOK) L.AddEntry(&fFit,"p_{0}+p_{1}E^{-1/2}","l");
      L.Draw();
  }
  TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.037);
  if (fitOK)
      tl.DrawLatex(0.18,0.25,
               Form("#sigma_{RMS}_{x}=%.2f + %.2f E^{-1/2}  (mm)",
                    fFit.GetParameter(0), fFit.GetParameter(1)));
  else
      tl.DrawLatex(0.18,0.25,"fit failed");

  /* ------------------------------------------------------------------ */
  /* 6) Output                                                          */
  /* ------------------------------------------------------------------ */
  if ( gSystem->AccessPathName(outDir,kWritePermission)!=kFALSE )
      gSystem->mkdir(outDir,kTRUE);

  TString outPng = Form("%s/PhenixSummary_AshOnly.png", outDir.Data());
  c.SaveAs(outPng);
  std::cout << "[PHX‑Σ] wrote  " << outPng << "\n[PHX‑Σ] DONE\n";
}


/**
 * \brief Main function that orchestrates the scanning and plotting
 *        of Ash-b and Log-w0 histograms.
 *
 * Usage:   root -l plotAshLogRMS_sideBySide.cpp\(\"PositionDep_sim_ALL.root\"\)
 */
std::map<double,double>
plotAshLogRMS_sideBySide(const char* infile = "PositionDep_sim_ALL.root")
{

  // 1) Build the b-scan in *tower-width units* 0.00 … 0.50, step 0.01
  //    (histograms were booked with exactly these values)
  std::vector<double> bScan_cm;               // keep the same variable name
  for (double b = 0.01; b <= 0.50 + 1e-9; b += 0.01)   // skip 0.00
      bScan_cm.push_back( std::round(b * 10000.) / 10000. );   // 4-dec precision

  // 2) Build the w0Scan => [1.5..7.0], step=0.1
  std::vector<double> w0Scan;
  for(double w=1.5; w<=7.0+1e-9; w+=0.1)
  {
    double val = std::round(w*100.)/100.;
    w0Scan.push_back(val);
  }

  // 3) Output directory
  TString baseDir    = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
  TString histOutDir = baseDir + "/ASH_LOG_PLOTS";

  gSystem->mkdir(baseDir, true);
  gSystem->mkdir(histOutDir,true);

  // 4) Open input file
  std::unique_ptr<TFile> fIn( TFile::Open(infile,"READ") );
  if(!fIn || fIn->IsZombie())
  {
    std::cerr<<"[ERROR] => Could not open file='"<<infile<<"'. Aborting.\n";
    return std::map<double,double>();   // or simply:  return {};
  }
  std::cout<<"[INFO] => Successfully opened '"<<infile<<"'.\n";

  // Global style
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);

  // 5) We'll do two approaches:
  //    (A) "FIT" => uses coreGaussianSigma
  //    (B) "RMS" => uses rawRMS
  // and plot side-by-side.

  // 5A) "FIT"
  std::cout<<"\n=== [STEP] Doing *FIT* approach with coreGaussianSigma ===\n";
  auto ashFit = doAshScan(fIn.get(), N_E, E_edges, binMode,
                            bScan_cm, 5.55,
                            coreGaussianSigma,
                            "fit", histOutDir,
                            /*skipLowest3=*/true);   // << add

  auto logFit = doLogScan(fIn.get(), N_E, E_edges, binMode,
                            w0Scan, coreGaussianSigma,
                            "fit", histOutDir,
                            /*skipLowest3=*/true);   // << add

  drawAshLogSideBySide(ashFit, logFit, "FIT",
                         E_edges, N_E, baseDir,
                         /*skipLowest3=*/true);      // << add

  // 5B) "RMS"
  std::cout<<"\n=== [STEP] Doing *RMS* approach with rawRMS ===\n";
  // Wrap rawRMS in a lambda if you like, or use directly
  auto ashRMS = doAshScan(fIn.get(),
                          N_E,
                          E_edges,
                          binMode,
                          bScan_cm,
                          5.55,
                          rawRMS,
                          "rms",
                          histOutDir,
                          true);
  auto logRMS = doLogScan(fIn.get(),
                          N_E,
                          E_edges,
                          binMode,
                          w0Scan,
                          rawRMS,
                          "rms",
                          histOutDir,
                          true);
  drawAshLogSideBySide(ashRMS, logRMS, "RMS", E_edges, N_E, baseDir, true);
    
  makePhenixLikeSummary(ashRMS, logRMS,
                          E_edges, N_E,
                          5.55,              // cell size (cm)
                          baseDir,           // output directory
                          /*skipLowest3=*/true);   // or false
    
  std::cout<<"\n[INFO] => plotAshLogRMS_sideBySide completed all tasks.\n\n";
    
  // ---------- NEW: build the look-up table for later use ---------- //
  std::map<double,double> bOpt;           // {E_slice_centre  →  b_RMS_opt}
  for (int i = 0; i < N_E; ++i)
  {
      const double eCtr = 0.5*(E_edges[i] + E_edges[i+1]);   // 2→4 → 3 GeV, …
      const double bVal = (i < (int)ashRMS.bestParam.size())
                           ? ashRMS.bestParam[i]               // optimal b from RMS scan
                           : std::numeric_limits<double>::quiet_NaN();
      bOpt[eCtr] = bVal;
  }
  return bOpt;                          // <── hand the table to the caller
}

/**
 * \brief asinhModel
 * Y(X) = Norm × [ 2·b / sqrt(1 + 4·X² · sinh²(1/(2·b)) ) ]
 */
double asinhModel(double* x, double* par)
{
  double Norm = par[0];
  double b    = par[1];

  if(b < 1e-9) return 1e-12; // prevent division-by-zero anomalies

  double X   = x[0];
  double sh  = TMath::SinH(1.0/(2.0*b));
  double arg = 1.0 + 4.0*X*X*sh*sh;

  return Norm * ( (2.0*b) / TMath::Sqrt(arg) );
}


// ===========================================================================
void FitLocalPhiEta(TH3F*  hUnc3D,
                    TH3F*  hCor3D,
                    bool   isFirstPass,
                    const std::vector<std::pair<double,double>>& eEdges,
                    const char* outDir,
                    std::ofstream& bOut)
// ===========================================================================
{
    // ----------( B )  Local-φ ------------------------------------------------
    TCanvas cFitsPhi("LocalPhiFits_4by2","Local #phi fits",1600,1200);
    cFitsPhi.Divide(4,2);

    const int N_E = static_cast<int>(eEdges.size());
    std::vector<std::unique_ptr<TH1D>> vPhiUnc , vPhiCor ;
    std::vector<std::unique_ptr<TF1>>  vPhiFit ;

    for (int i = 0; i < N_E; ++i)
    {
        double eLo = eEdges[i].first , eHi = eEdges[i].second;
        std::cout << "\n[Φ] slice " << i << "  E=[" << eLo << ',' << eHi << ")\n";
        cFitsPhi.cd(i+1);

        int zLo = std::max(1 , hUnc3D->GetZaxis()->FindBin(eLo+1e-9));
        int zHi = std::min(hUnc3D->GetNbinsZ(),
                           hUnc3D->GetZaxis()->FindBin(eHi-1e-9));

        hUnc3D->GetZaxis()->SetRange(zLo,zHi);
        hUnc3D->GetXaxis()->SetRange(1,hUnc3D->GetNbinsX());

        auto hUnc = std::unique_ptr<TH1D>(
                     static_cast<TH1D*>(hUnc3D->Project3D("y")));
        if(!hUnc){ std::cerr<<"[WARN] no uncor φ slice\n"; continue; }

        hUnc->SetDirectory(nullptr);
        hUnc->SetName (Form("hUncPhi_E%d",i));
        hUnc->SetTitle(Form("Local #phi : E=[%.1f,%.1f)",eLo,eHi));
        hUnc->GetXaxis()->SetTitle("block #phi_{local, 2#times 2}");
        hUnc->GetYaxis()->SetTitle("counts");
        hUnc->GetYaxis()->SetTitleOffset(1.81);
        hUnc->SetMarkerStyle(20);
        hUnc->Draw("E");
        
        // ---------- ❶ INSERT HEADER --------------------------------------
        { TLatex hd;  hd.SetNDC();
          hd.SetTextFont(44);             // CMS style
          hd.SetTextSize(0.045);
          hd.SetTextAlign(22);            // centred
          hd.DrawLatex(0.54,0.97,         // (x,y) in pad coords
                       Form("#bf{%.1f #minus %.1f  GeV}", eLo, eHi));
        }
        // -----------------------------------------------------------------

        vPhiUnc.emplace_back(std::move(hUnc));

        // ---------- 2)  asinh fit ------------------------------------------
        TF1 trial("asinhφ",asinhModel,-0.5,0.5,2);
        trial.SetParNames("Norm","bVal");
        trial.SetParLimits(0,1e-9,1e9);
        trial.SetParLimits(1,1e-5,2.0);

        double   bestChi2 = 1e50;       int bestNdf = 0;
        BRes     best;                      // central value + error
        std::unique_ptr<TF1> bestFit;

        for(auto g : {std::pair{0.1,0.10}, {0.2,0.14},
                      std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(g.first,g.second);
            TFitResultPtr r = vPhiUnc.back()->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2() < bestChi2){
                bestChi2 = r->Chi2();
                bestNdf  = r->Ndf();
                best.val = trial.GetParameter(1);
                best.err = trial.GetParError(1);
                
                bestFit = std::make_unique<TF1>(trial);
                bestFit->SetName(Form("bestF_phi_E%d", i));
            }
        }

        double maxU = vPhiUnc.back()->GetMaximum();
        vPhiUnc.back()->GetYaxis()->SetRangeUser(0,1.3*std::max(maxU,1.));

        if(bestFit){
            bestFit->SetLineColor(kBlue+1); bestFit->SetLineWidth(2);
            bestFit->Draw("SAME"); vPhiFit.emplace_back(std::move(bestFit));
        }

        // ---------- 4)  Corrected overlay ----------------------------------
        TH1D* hCorRaw=nullptr;
        if(!isFirstPass && hCor3D){
            hCor3D->GetZaxis()->SetRange(zLo,zHi);
            hCor3D->GetXaxis()->SetRange(1,hCor3D->GetNbinsX());
            auto hCor = std::unique_ptr<TH1D>(
                         static_cast<TH1D*>(hCor3D->Project3D("y")));
            if(hCor){
                hCor->SetDirectory(nullptr);
                hCor->SetName(Form("hCorPhi_E%d",i));
                hCor->SetMarkerStyle(21); hCor->SetMarkerColor(kRed);
                hCor->SetLineColor(kRed); hCor->Draw("SAME E");
                double maxC=hCor->GetMaximum();
                if(maxC>maxU) vPhiUnc.back()->GetYaxis()
                                      ->SetRangeUser(0,1.3*maxC);
                hCorRaw=hCor.get(); vPhiCor.emplace_back(std::move(hCor));
            }
        }

        // ---------- 5)  Legend & TLatex ------------------------------------
        double nUnc=vPhiUnc.back()->GetEntries();
        double nCor=hCorRaw? hCorRaw->GetEntries():0.;
        TLegend lg(0.54,0.75,0.88,0.8); lg.SetBorderSize(0);
        lg.SetFillStyle(0); lg.SetTextSize(0.035);
        lg.AddEntry(vPhiUnc.back().get(),"Uncorrected","lp");
        if(hCorRaw) lg.AddEntry(hCorRaw,"Corrected","lp");
        if(bestFit) lg.AddEntry(bestFit.get(),
                                Form("asinh fit (b=%.3g)", best.val),"l");
        lg.Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042);
        tl.DrawLatex(0.25,0.88,Form("b = %.3g #pm %.3g",best.val,best.err));
        tl.DrawLatex(0.25,0.85,Form("#chi^{2}_{LL}/NDF = %.2f",
                        bestNdf>0? bestChi2/bestNdf : 0.));
        tl.SetTextAlign(32);
        tl.DrawLatex(0.88,0.85,Form("Uncorr: %.0f",nUnc));
        if(hCorRaw) tl.DrawLatex(0.88,0.82,Form("Corr: %.0f",nCor));

        if (bOut.is_open() && best.val>0)
            bOut << Form("PHI [%.1f,%.1f) %.6f\n",eLo,eHi,best.val);
    } // end φ loop

    cFitsPhi.SaveAs(Form("%s/LocalPhiFits_4by2.png",outDir));
    std::cout << "[INFO] φ canvas written\n";

    // ----------( C )  Local-η  --------------------------------------------
    TCanvas cFitsEta("LocalEtaFits_4by2","Local #eta fits",1600,1200);
    cFitsEta.Divide(4,2);

    std::vector<std::unique_ptr<TH1D>> vEtaUnc , vEtaCor ;
    std::vector<std::unique_ptr<TF1>>  vEtaFit ;

    for (int i = 0; i < N_E; ++i)
    {
        double eLo = eEdges[i].first , eHi = eEdges[i].second;
        std::cout << "\n[η] slice " << i << "  E=[" << eLo << ',' << eHi << ")\n";
        cFitsEta.cd(i+1);

        int zLo = std::max(1 , hUnc3D->GetZaxis()->FindBin(eLo+1e-9));
        int zHi = std::min(hUnc3D->GetNbinsZ(),
                           hUnc3D->GetZaxis()->FindBin(eHi-1e-9));

        hUnc3D->GetZaxis()->SetRange(zLo,zHi);
        hUnc3D->GetYaxis()->SetRange(1,hUnc3D->GetNbinsY());

        auto hUnc = std::unique_ptr<TH1D>(
                     static_cast<TH1D*>(hUnc3D->Project3D("x")));
        if(!hUnc){ std::cerr<<"[WARN] no uncor η slice\n"; continue; }

        hUnc->SetDirectory(nullptr);
        hUnc->SetName(Form("hUncEta_E%d",i));
        hUnc->SetTitle(Form("Local #eta : E=[%.1f,%.1f)",eLo,eHi));
        hUnc->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
        hUnc->GetYaxis()->SetTitle("counts");
        hUnc->SetMarkerStyle(20);
        hUnc->Draw("E");
        vEtaUnc.emplace_back(std::move(hUnc));

        TF1 trial("asinhη",asinhModel,-0.5,0.5,2);
        trial.SetParNames("Norm","bVal");
        trial.SetParLimits(0,1e-9,1e9);
        trial.SetParLimits(1,1e-5,2.0);

        double bestChi2 = 1e50;   int bestNdf = 0;
        BRes   best;                    // <-- declare THE SAME STRUCT here
        std::unique_ptr<TF1> bestFit;

        for(auto g : {std::pair{0.1,0.10}, {0.2,0.14},
                      std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(g.first,g.second);
            TFitResultPtr r = vEtaUnc.back()->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2()<bestChi2){
                bestChi2 = r->Chi2();
                bestNdf  = r->Ndf();
                best.val = trial.GetParameter(1);
                best.err = trial.GetParError(1);
                
                bestFit = std::make_unique<TF1>(trial);
                bestFit->SetName(Form("bestF_eta_E%d", i));
            }
        }

        double maxU=vEtaUnc.back()->GetMaximum();
        vEtaUnc.back()->GetYaxis()->SetRangeUser(0,1.3*std::max(maxU,1.));

        if(bestFit){
            bestFit->SetLineColor(kBlue+1); bestFit->SetLineWidth(2);
            bestFit->Draw("SAME"); vEtaFit.emplace_back(std::move(bestFit));
        }

        TH1D* hCorRaw=nullptr;
        if(!isFirstPass && hCor3D){
            hCor3D->GetZaxis()->SetRange(zLo,zHi);
            hCor3D->GetYaxis()->SetRange(1,hCor3D->GetNbinsY());
            auto hCor = std::unique_ptr<TH1D>(
                         static_cast<TH1D*>(hCor3D->Project3D("x")));
            if(hCor){
                hCor->SetDirectory(nullptr);
                hCor->SetName(Form("hCorEta_E%d",i));
                hCor->SetMarkerStyle(21); hCor->SetMarkerColor(kRed);
                hCor->SetLineColor(kRed); hCor->Draw("SAME E");
                double maxC=hCor->GetMaximum();
                if(maxC>maxU) vEtaUnc.back()->GetYaxis()
                                    ->SetRangeUser(0,1.3*maxC);
                hCorRaw=hCor.get(); vEtaCor.emplace_back(std::move(hCor));
            }
        }

        double nUnc=vEtaUnc.back()->GetEntries();
        double nCor=hCorRaw? hCorRaw->GetEntries():0.;
        TLegend lg(0.54,0.75,0.88,0.8); lg.SetBorderSize(0);
        lg.SetFillStyle(0); lg.SetTextSize(0.035);
        lg.AddEntry(vEtaUnc.back().get(),"Uncorrected","lp");
        if(hCorRaw) lg.AddEntry(hCorRaw,"Corrected","lp");
        if(bestFit) lg.AddEntry(bestFit.get(),
                                Form("asinh fit (b=%.3g)",best.val),"l");
        lg.Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042);
        tl.DrawLatex(0.14,0.84,Form("b = %.3g #pm %.3g",best.val,best.err));
        tl.DrawLatex(0.14,0.76,Form("#chi^{2}_{LL}/NDF = %.2f",
                        bestNdf>0? bestChi2/bestNdf : 0.));
        tl.SetTextAlign(32);
        tl.DrawLatex(0.88,0.85,Form("Uncorr: %.0f",nUnc));
        if(hCorRaw) tl.DrawLatex(0.88,0.82,Form("Corr: %.0f",nCor));

        if(bOut.is_open() && best.val>0)
            bOut << Form("ETA [%.1f,%.1f)  %.6f\n", eLo, eHi, best.val);
    } // end η loop

    cFitsEta.SaveAs(Form("%s/LocalEtaFits_4by2.png",outDir));
    std::cout << "[INFO] η canvas written\n";
}




void PlotPhiShiftAndWidth(TH3F*  hUnc3D,
                          TH3F*  hCor3D,
                          const std::vector<std::pair<double,double>>& eEdges,
                          const char* outDir)
{
  if(!hUnc3D || !hCor3D){
      std::cerr << "[PhiShift] need both uncorrected *and* corrected 3-D hists – abort\n";
      return;
  }

  /* ————————————————————————
     Style tuned for publications
     ———————————————————————— */
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont  (42,"XYZ");
  gStyle->SetLabelFont  (42,"XYZ");
  gStyle->SetTitleSize  (0.045,"XYZ");
  gStyle->SetLabelSize  (0.040,"XYZ");
  gStyle->SetPadTickX(1);   // ticks on all four sides
  gStyle->SetPadTickY(1);

  const int N_E = (int)eEdges.size();
  std::vector<double> vEcen, vShift, vShiftErr, vRmsRatio;

  /* ———————————————————————————————————————————
     loop: project each (η,φ,E) cube on local φ
     ——————————————————————————————————————————— */
  for(int i=0;i<N_E;++i)
  {
      const double eLo=eEdges[i].first , eHi=eEdges[i].second ;
      const double eC =0.5*(eLo+eHi);

      auto slice = [&](TH3F* h)->std::unique_ptr<TH1D>
      {
          const int zLo = std::max(1 , h->GetZaxis()->FindBin(eLo+1e-9));
          const int zHi = std::min(h->GetNbinsZ(),
                                   h->GetZaxis()->FindBin(eHi-1e-9));
          h->GetZaxis()->SetRange(zLo,zHi);
          h->GetXaxis()->SetRange(1,h->GetNbinsX());   // full η range
          auto hphi = std::unique_ptr<TH1D>(
                         (TH1D*)h->Project3D("y") );
          hphi->SetDirectory(nullptr);
          return hphi;
      };

      auto hU = slice(hUnc3D);
      auto hC = slice(hCor3D);
      if(!hU || !hC) continue;

      const double meanU  = hU->GetMean();
      const double meanC  = hC->GetMean();
      const double rmsU   = hU->GetRMS();
      const double rmsC   = hC->GetRMS();
      const double nU     = std::max(1.0 , hU->GetEntries());

      vEcen    .push_back(eC);
      vShift   .push_back(meanC - meanU);
      vShiftErr.push_back(rmsU/std::sqrt(nU));      // statistical error
      vRmsRatio.push_back(rmsU>0 ? rmsC/rmsU : 0.0);
  }

  if(vEcen.empty()){
      std::cerr<<"[PhiShift] no valid slices\n"; return;
  }

  /* ———————————————————————————————————————————
     build graphs
     ——————————————————————————————————————————— */
  const int N = (int)vEcen.size();
  TGraphErrors gShift (N, &vEcen[0], &vShift[0]   , nullptr, &vShiftErr[0]);
  TGraph       gRms   (N, &vEcen[0], &vRmsRatio[0]);

  gShift.SetMarkerStyle(20);
  gShift.SetMarkerColor(kRed+1);  gShift.SetLineColor(kRed+1);
  gShift.SetLineWidth(2);

  gRms.SetMarkerStyle(25);
  gRms.SetMarkerColor(kAzure+2);  gRms.SetLineColor(kAzure+2);
  gRms.SetLineWidth(2);

  /* ———————————————————————————————————————————
     lay-out: two pads sharing the x-axis
     ——————————————————————————————————————————— */
  const double xMin = eEdges.front().first  - 0.5;
  const double xMax = eEdges.back ().second + 0.5;

  TCanvas c("PhiShift_vs_E","",1000,900);
  c.cd();
  TPad *p1 = new TPad("p1","",0,0.42,1,1);   // top 58 %
  TPad *p2 = new TPad("p2","",0,0   ,1,0.42); // bottom 42 %
  for(auto p : {p1,p2}){
      p->SetLeftMargin (0.14);
      p->SetRightMargin(0.04);
      p->SetTickx(); p->SetTicky();
      p->Draw();
  }

  /* ————————————————————
     (a) mean shift
     ———————————————————— */
  p1->cd();
  TH2F f1("f1",";E_{slice}  [GeV];#Delta#LT#varphi#GT  (corr - raw)",
          100,xMin,xMax, 100, -0.05, 0.05);
  f1.GetYaxis()->SetTitleOffset(1.2);
  f1.Draw();
  gShift.Draw("PZ SAME");

  TLine l0(xMin,0,xMax,0); l0.SetLineStyle(2); l0.Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.16,0.85,"#bf{(a)}  Mean shift");

  /* ————————————————————
     (b) width ratio
     ———————————————————— */
  p2->cd();
  TH2F f2("f2",";E_{slice}  [GeV];RMS_{corr} / RMS_{raw}",
          100,xMin,xMax, 100, 0.7, 1.25);
  f2.GetYaxis()->SetTitleOffset(1.35);
  f2.GetXaxis()->SetTitleOffset(1.2);
  f2.Draw();
  gRms.Draw("P SAME");

  TLine l1(xMin,1,xMax,1); l1.SetLineStyle(2); l1.Draw();

  lat.DrawLatex(0.16,0.84,"#bf{(b)}  Width ratio");

  /* ————————————————————
     Legend – unobtrusive, transparent
     ———————————————————— */
  TPaveText leg(0.70,0.78,0.93,0.90,"NDC");
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextAlign(12); leg.SetTextFont(42); leg.SetTextSize(0.036);
  leg.AddText( Form("#color[%d]{#square}  #Delta#LT#varphi#GT", kRed+1) );
  leg.AddText( Form("#color[%d]{#square}  RMS ratio",          kAzure+2) );
  p1->cd();      // place it in the upper pad
  leg.Draw();

  /* ————————————————————
     write file
     ———————————————————— */
  TString out = Form("%s/PhiShift_vs_E.png",outDir);
  c.SaveAs(out);
  std::cout<<"[PhiShift] wrote "<<out<<"\n";
}

void Plot2DBlockEtaPhi(TH3F* hUnc3D,
                       TH3F* hCor3D,
                       bool  isFirstPass,
                       const std::vector<std::pair<double,double>>& eEdges,
                       const char* outDir)
{
  // Number of energy slices
  const int N_E = eEdges.size(); // typically 8

  const double tSize = 0.045;     // title font size
  const double lSize = 0.035;     // label font size

  const double tOffX = 1.10;      // push X-title away from axis
  const double tOffY = 1.35;      // push Y-title away from axis

  gStyle->SetPadTickX(1);         // ticks on all sides
  gStyle->SetPadTickY(1);

  auto tuneAxes = [&](TH2* h)     // NEW helper
  {
      h->GetXaxis()->SetTitleSize(tSize);
      h->GetYaxis()->SetTitleSize(tSize);
      h->GetZaxis()->SetTitleSize(tSize);

      h->GetXaxis()->SetLabelSize(lSize);
      h->GetYaxis()->SetLabelSize(lSize);
      h->GetZaxis()->SetLabelSize(lSize);

      h->GetXaxis()->SetTitleOffset(tOffX);
      h->GetYaxis()->SetTitleOffset(tOffY);
      h->GetZaxis()->SetTitleOffset(1.8);       // Z needs more room
  };
  /* ------------------------------------------------------------ */


  // ===================================================================
  // (A) 2D blockEta vs blockPhi (4×2)
  // ===================================================================
  TCanvas c2D_unc("c2D_unc","Uncorrected 2D block coords vs E",2400,1200);
  c2D_unc.Divide(4,2);
  for (int p = 1; p <= 8; ++p) {
        TPad *pad = (TPad*) c2D_unc.cd(p);
        pad->SetRightMargin(0.18);   // or whatever fraction you like
        pad->SetLeftMargin (0.12);
        pad->SetBottomMargin(0.12);
  }

  TCanvas* c2D_cor = nullptr;
  if(!isFirstPass)
  {
    c2D_cor = new TCanvas("c2D_cor","Corrected 2D block coords vs E",2400,1200);
    c2D_cor->Divide(4,2);
    for (int p = 1; p <= 8; ++p) {
          TPad *pad = (TPad*) c2D_cor->cd(p);   //  <-  use ->
          pad->SetRightMargin(0.18);
          pad->SetLeftMargin (0.18);
          pad->SetBottomMargin(0.12);
    }
  }

  for(int i=0; i<N_E; i++)
  {
    double eLo = eEdges[i].first;
    double eHi = eEdges[i].second;

    const TString hdrTxt = Form("E = %.1f - %.1f  GeV", eLo, eHi);   // ASCII "-"

    // locate Z bins in uncorrected
    int zLo = hUnc3D->GetZaxis()->FindBin(eLo+1e-9);
    int zHi = hUnc3D->GetZaxis()->FindBin(eHi-1e-9);
    if(zLo<1) zLo=1;
    if(zHi>hUnc3D->GetNbinsZ()) zHi=hUnc3D->GetNbinsZ();
    hUnc3D->GetZaxis()->SetRange(zLo,zHi);

    TH2D* h2_unc = (TH2D*) hUnc3D->Project3D("xy");
    if (h2_unc) {
          h2_unc->SetName(Form("h2_unc_xy_Ebin%d", i));   // ❶ UNIQUE NAME
          h2_unc->SetTitle(Form("Uncorr: E=[%.1f,%.1f)", eLo, eHi));
          c2D_unc.cd(i+1);
          h2_unc->GetZaxis()->SetTitle("Energy [GeV]");
          h2_unc->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
          h2_unc->GetYaxis()->SetTitle("block #phi_{local, 2#times 2}");
          h2_unc->GetZaxis()->CenterTitle();          // optional
          h2_unc->GetZaxis()->SetTitleOffset(2.1);    // move label away
          tuneAxes(h2_unc);
          h2_unc->GetZaxis()->SetTitleOffset(1.3);
          h2_unc->Draw("COLZ");
          // do NOT delete here – wait until after SaveAs()
          /* TLatex header ------------------------------------------------ */
          TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
          tl.SetTextSize(0.045); tl.SetTextAlign(22);          // centred
          tl.DrawLatex(0.50, 0.97,
                     Form("#bf{UNCORR:  %s}", hdrTxt.Data()));   // E-range after the colon
    }
    else
    {
      std::cerr << "[WARN] Could not get uncorrected 2D for slice i=" << i << "\n";
    }

    // corrected if available
    if(!isFirstPass && hCor3D)
    {
      int zLoC = hCor3D->GetZaxis()->FindBin(eLo+1e-9);
      int zHiC = hCor3D->GetZaxis()->FindBin(eHi-1e-9);
      if(zLoC<1) zLoC=1;
      if(zHiC>hCor3D->GetNbinsZ()) zHiC=hCor3D->GetNbinsZ();
      hCor3D->GetZaxis()->SetRange(zLoC,zHiC);

      TH2D* h2_cor = (TH2D*) hCor3D->Project3D("xy");
      if (h2_cor) {
            h2_cor->SetName(Form("h2_cor_xy_Ebin%d", i));   // ❶ UNIQUE NAME
            h2_cor->SetTitle(Form("Corr: E=[%.1f,%.1f)", eLo, eHi));
            c2D_cor->cd(i+1);
            // NB: The original code uses "Uncorr" in the next line (likely a typo),
            // but we keep it EXACTLY as requested:
            h2_cor->SetTitle(Form("Corr: E=[%.1f,%.1f)", eLo, eHi));
            h2_cor->GetXaxis()->SetTitle("blockEta_{local, 2#times 2}");
            h2_cor->GetYaxis()->SetTitle("blockPhi_{local, 2#times 2}");
            h2_cor->GetZaxis()->SetTitle("Energy [GeV]");
            h2_cor->GetZaxis()->SetTitleOffset(2.1);
            h2_cor->GetZaxis()->CenterTitle();
            h2_cor->GetZaxis()->SetTitleOffset(1.3);
            tuneAxes(h2_cor);
            h2_cor->Draw("COLZ");
          
            TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
            tl.SetTextSize(0.045); tl.SetTextAlign(22);
            
            tl.DrawLatex(0.50, 0.97,
                       Form("#bf{CORR:    %s}", hdrTxt.Data()));
      }
      else
      {
        std::cerr << "[WARN] Could not project corrected 2D => i="<< i <<"\n";
      }
    }
  } // end 2D loop

  // Save results
  TString out2D_unc = Form("%s/BlockCoord2D_E_unc.png", outDir);
  c2D_unc.SaveAs(out2D_unc);
  std::cout << "[INFO] Wrote uncorrected => " << out2D_unc << "\n";

  if(!isFirstPass && c2D_cor)
  {
    TString out2D_cor=Form("%s/BlockCoord2D_E_cor.png", outDir);
    c2D_cor->SaveAs(out2D_cor);
    std::cout << "[INFO] Wrote corrected => " << out2D_cor << "\n";
    delete c2D_cor;
    c2D_cor = nullptr;
  }
}


// ---------------------------------------------------------------------------
//  Overlay of *UNCORRECTED* η and φ block-CG coordinates, in one canvas
//  • φ   (Y-proj)  = red   + red asinh fit
//  • η   (X-proj)  = blue  + blue asinh fit
//  • b-values (fit parameters) printed upper-right on every pad
//  • eight energy bins  → 4×2 pads  →  LocalPhiEtaOverlay_4by2.png
// ---------------------------------------------------------------------------
void OverlayUncorrPhiEta(TH3F*  hUnc3D,
                         const std::vector<std::pair<double,double>>& eEdges,
                         const char* outDir)
/* ------------------------------------------------------------------------ */
{
    if(!hUnc3D){
        std::cerr << "[Overlay] null hUnc3D – abort\n";
        return;
    }

    gStyle->SetOptStat(0);

    const int N_E = static_cast<int>(eEdges.size());   // usually 8

    TCanvas cOv("LocalPhiEtaOverlay_4by2",
                "Uncorr #varphi (red) and #eta (blue) overlays", 1600, 1200);
    cOv.Divide(4, 2);

    // keep everything alive until after SaveAs()
    std::vector<TObject*> keepAlive;

    //----------------------------------------------------------------------
    // helper: do a binned-LL fit to the asinh model and return (b, TF1*)
    //----------------------------------------------------------------------
    auto fitAsinh = [&](TH1D* h, Color_t col) -> std::pair<double, TF1*>
    {
        TF1 trial("trial", asinhModel, -0.5, 0.5, 2);
        trial.SetParNames("Norm", "b");
        trial.SetParLimits(0, 1e-9, 1e9);
        trial.SetParLimits(1, 1e-5, 2.0);
        double bestChi2 = 1e50;
        BRes   best;                 // <-- declare a local result holder
        TF1*   bestFit = nullptr;

        for (const auto& seed : std::initializer_list<std::pair<double,double>>
                                {{0.1,0.10},{0.2,0.14},{0.3,0.20},{0.5,0.08}})
        {
            trial.SetParameters(seed.first, seed.second);
            TFitResultPtr r = h->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;

            if(r->Chi2() < bestChi2){
                bestChi2 = r->Chi2();
                best.val = trial.GetParameter(1);
                best.err = trial.GetParError(1);

                if(bestFit) delete bestFit;
                bestFit = new TF1(trial);
            }
        }
        if(bestFit){
            bestFit->SetLineColor(col);
            bestFit->SetLineWidth(2);
        }
        return {best.val, bestFit};
    };

    //----------------------------------------------------------------------
    // main loop over energy slices
    //----------------------------------------------------------------------
    for(int i = 0; i < N_E; ++i)
    {
        const double eLo = eEdges[i].first;
        const double eHi = eEdges[i].second;

        std::cout << "[Overlay] slice " << i
                  << "  E=[" << eLo << ',' << eHi << ")\n";

        cOv.cd(i + 1);

        // ---- set Z-range for this energy slice --------------------------
        const int zLo = std::max(1,
              hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
        const int zHi = std::min(hUnc3D->GetNbinsZ(),
              hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));

        // -----------------------------------------------------------------
        //  φ  distribution  (Y-projection)
        // -----------------------------------------------------------------
        hUnc3D->GetZaxis()->SetRange(zLo, zHi);
        hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX()); // full η

        TH1D* hPhi = static_cast<TH1D*>( hUnc3D->Project3D("y") );
        if(!hPhi){ std::cerr << "[Overlay] missing φ hist\n"; continue; }
        hPhi->SetDirectory(nullptr);      // detach from current file
        hPhi->SetName(Form("hPhi_E%d", i));
        hPhi->SetMarkerStyle(21);
        hPhi->SetMarkerColor(kRed);
        hPhi->SetLineColor  (kRed);

        // -----------------------------------------------------------------
        //  η  distribution  (X-projection)
        // -----------------------------------------------------------------
        hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY()); // full φ

        TH1D* hEta = static_cast<TH1D*>( hUnc3D->Project3D("x") );
        if(!hEta){ std::cerr << "[Overlay] missing η hist\n"; continue; }
        hEta->SetDirectory(nullptr);
        hEta->SetName(Form("hEta_E%d", i));
        hEta->SetMarkerStyle(20);
        hEta->SetMarkerColor(kBlue + 1);
        hEta->SetLineColor  (kBlue + 1);

        // -----------------------------------------------------------------
        //  fits
        // -----------------------------------------------------------------
        auto [bPhi, fPhi] = fitAsinh(hPhi, kRed);
        auto [bEta, fEta] = fitAsinh(hEta, kBlue + 1);

        // -----------------------------------------------------------------
        //  draw – set a common y-range   <-- keep your existing lines
        // -----------------------------------------------------------------
        const double ymax = std::max(hPhi->GetMaximum(), hEta->GetMaximum());
        hEta->GetYaxis()->SetRangeUser(0.0, 1.30 * std::max(1.0, ymax));

        /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
        // (1) unified x–axis label for BOTH histograms
        hEta->GetXaxis()->SetTitle("#phi_{CG}/#eta_{CG}");
        hPhi->GetXaxis()->SetTitle("#phi_{CG}/#eta_{CG}");

        // (2) pad title: energy window + “uncorrected”
        TString sliceTitle = Form("Uncorrected: %.1f < E < %.1f  GeV", eLo, eHi);
        hEta->SetTitle(sliceTitle);   // the first histogram drawn defines the title
        hPhi->SetTitle(sliceTitle);   // keep both in sync
        /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

        hEta->Draw("E");            // blue first → defines axes & title
        hPhi->Draw("E SAME");       // then red
        if(fEta) fEta->Draw("SAME");
        if(fPhi) fPhi->Draw("SAME");


        TPad* pad = (TPad*) gPad;                  // current pad
        TLegend* leg = new TLegend(0.2, 0.8,     // (x1,y1,x2,y2) in NDC
                                   0.45, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.040);
        leg->AddEntry(hPhi, "local #varphi", "lp");   // red squares
        leg->AddEntry(hEta, "local #eta",     "lp");  // blue circles
        leg->Draw();

        /* keep it alive until after SaveAs() */
        keepAlive.push_back(leg);
        
        // -----------------------------------------------------------------
        //  text
        // -----------------------------------------------------------------
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.040);
        tl.SetTextAlign(32);               // right-aligned
        tl.DrawLatex(0.88, 0.88, Form("b_{#varphi}=%.3g", bPhi));
        tl.DrawLatex(0.88, 0.82, Form("b_{#eta}=%.3g",     bEta));

        // keep objects alive
        keepAlive.insert(keepAlive.end(),
                         {hPhi, hEta, fPhi, fEta});
    } // end slice loop

    //----------------------------------------------------------------------
    //  save and clean-up
    //----------------------------------------------------------------------
    TString outName = Form("%s/LocalPhiEtaOverlay_4by2.png", outDir);
    cOv.SaveAs(outName);
    std::cout << "[Overlay] wrote " << outName << '\n';

//    for(auto obj : keepAlive) delete obj;
}

// ---------------------------------------------------------------------------
//  Scatter-plot of best-fit 𝗯 values vs. energy slice
//  • red squares = φ    (Y-projection)   b_{φ}
//  • blue circles = η   (X-projection)   b_{η}
//  • y-axis limits chosen automatically to include *all* points
// ---------------------------------------------------------------------------
std::map<double,double>
PlotBvaluesVsEnergy(TH3F* hUnc3D,
                    const std::vector<std::pair<double,double>>& eEdges,
                    const char* outDir)
{
    std::map<double,double> bMap;
    if(!hUnc3D){
        std::cerr << "[bPlot] null hUnc3D – abort\n";
        return bMap;     // ← FIX: return an empty map, not “nothing”
    }

    gStyle->SetOptStat(0);

    const int N_E = static_cast<int>(eEdges.size());

    std::vector<double> vX, vBphi, vBeta, vBphiErr, vBetaErr;
    vX.reserve(N_E);  vBphi.reserve(N_E);  vBeta.reserve(N_E);
    vBphiErr.reserve(N_E);  vBetaErr.reserve(N_E);

    // ───────────────────────────────────────────────────────────────────
    // helper: identical to the fitter used elsewhere
    auto fitAsinh = [&](TH1D* h)->BRes
    {
        TF1 trial("trial", asinhModel, -0.5, 0.5, 2);
        trial.SetParNames("Norm","b");
        trial.SetParLimits(0,1e-9,1e9);
        trial.SetParLimits(1,1e-5,2.0);

        double bestChi2 = 1e50;
        BRes   best;                   // {val = NaN, err = 0}

        for(const auto& seed : std::initializer_list<std::pair<double,double>>
                               {{0.1,0.10},{0.2,0.14},{0.3,0.20},{0.5,0.08}})
        {
            trial.SetParameters(seed.first, seed.second);
            TFitResultPtr r = h->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2() < bestChi2){
                bestChi2 = r->Chi2();
                best.val = trial.GetParameter(1);
                best.err = trial.GetParError(1);
            }
        }
        return best;
    };

    std::vector<TObject*> keepAlive;        // stop auto-deletion

    // ───────────────────────────────────────────────────────────────────
    //   loop over energy bins – compute bφ and bη
    // ───────────────────────────────────────────────────────────────────
    for(int i=0; i<N_E; ++i)
    {
        const double eLo = eEdges[i].first;
        const double eHi = eEdges[i].second;

        const int zLo = std::max(1,
              hUnc3D->GetZaxis()->FindBin(eLo + 1e-9));
        const int zHi = std::min(hUnc3D->GetNbinsZ(),
              hUnc3D->GetZaxis()->FindBin(eHi - 1e-9));

        // φ  (Y-projection)
        hUnc3D->GetZaxis()->SetRange(zLo, zHi);
        hUnc3D->GetXaxis()->SetRange(1, hUnc3D->GetNbinsX());

        TH1D* hPhi = static_cast<TH1D*>( hUnc3D->Project3D("y") );
        if(!hPhi){ std::cerr << "[bPlot] missing φ hist, slice "<<i<<"\n"; continue;}
        hPhi->SetDirectory(nullptr);

        // η  (X-projection)
        hUnc3D->GetYaxis()->SetRange(1, hUnc3D->GetNbinsY());

        TH1D* hEta = static_cast<TH1D*>( hUnc3D->Project3D("x") );
        if(!hEta){ std::cerr << "[bPlot] missing η hist, slice "<<i<<"\n"; delete hPhi; continue;}
        hEta->SetDirectory(nullptr);

        const BRes bPhi = fitAsinh(hPhi);
        const BRes bEta = fitAsinh(hEta);

        vX   .push_back( 0.5*(eLo+eHi) );   // use bin centre
        vBphi.push_back(bPhi.val);
        vBeta.push_back(bEta.val);
        
        vBphiErr.push_back(bPhi.err);
        vBetaErr.push_back(bEta.err);

        keepAlive.insert(keepAlive.end(), {hPhi, hEta});
    }

    // ───────────────────────────────────────────────────────────────────
    //  build graphs
    // ───────────────────────────────────────────────────────────────────
    TGraphErrors* gPhi = new TGraphErrors( vX.size(), &vX[0], &vBphi[0],
                                          nullptr,  &vBphiErr[0] );
    
    TGraphErrors* gEta = new TGraphErrors( vX.size(), &vX[0], &vBeta[0],
                                          nullptr,  &vBetaErr[0] );
    
    gPhi->SetMarkerStyle(21);  gPhi->SetMarkerColor(kRed);
    gPhi->SetLineColor  (kRed);

    gEta->SetMarkerStyle(20);  gEta->SetMarkerColor(kBlue+1);
    gEta->SetLineColor  (kBlue+1);

    // y-axis limits
    const double ymin = std::min( *std::min_element(vBphi.begin(),vBphi.end()),
                                  *std::min_element(vBeta.begin(),vBeta.end()) );
    const double ymax = std::max( *std::max_element(vBphi.begin(),vBphi.end()),
                                  *std::max_element(vBeta.begin(),vBeta.end()) );

    const double yLo = 0.85 * ymin;
    const double yHi = 1.15 * ymax;

    // frame for axes
    TH2F* hFrame = new TH2F("bFrame", ";E_{slice}  [GeV]; best–fit  b",
                            100, vX.front()-1, vX.back()+1,
                            100, yLo, yHi);
    hFrame->SetStats(0);

    TCanvas cB("bValues_vs_E","b values vs energy", 800, 600);
    hFrame->Draw();
    gPhi->Draw("P SAME");
    gEta->Draw("P SAME");

    TLegend leg(0.2,0.78,0.40,0.90);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
    leg.AddEntry(gPhi,"b_{#varphi}","lp");
    leg.AddEntry(gEta,"b_{#eta}","lp");
    leg.Draw();

    keepAlive.insert(keepAlive.end(), {gPhi,gEta,hFrame});

    TString outName = Form("%s/bValues_vs_E.png", outDir);
    cB.SaveAs(outName);
    std::cout << "[bPlot] wrote " << outName << '\n';

//    for(auto obj : keepAlive) delete obj;
    
    // ---------- NEW single line ---------- //
    for(size_t i=0;i<vX.size();++i) bMap[vX[i]] = vBphi[i];
    return bMap;
}


/***********************************************************************
 * PlotChi2QA ‑‑ *v2.1* (2025‑06‑10)

 **********************************************************************/
void PlotChi2QA(const char* inFile,
                const char* outDir = "./Chi2_QA",
                int nWorstToPrint  = 12)
{
  /* ================================================================
   * 0) Style
   * ================================================================
   */
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(255);

  gSystem->mkdir(outDir,/*recursive=*/true);

  /* ================================================================
   * 1) Input histograms
   * ================================================================
   */
  std::unique_ptr<TFile> fin(TFile::Open(inFile,"READ"));
  if(!fin || fin->IsZombie()){
     std::cerr<<"[Chi2QA]  cannot open "<<inFile<<"\n"; return;
  }
  auto* hTot  = dynamic_cast<TH2F*    >(fin->Get("h2_chi2_tot_etaPhi"));
  auto* hRej  = dynamic_cast<TH2F*    >(fin->Get("h2_chi2_rej_etaPhi"));
  auto* pPass = dynamic_cast<TProfile2D*>(fin->Get("p_chi2_pass_etaPhi"));
  if(!hTot||!hRej||!pPass){
     std::cerr<<"[Chi2QA]  missing QA objects in "<<inFile<<"\n"; return;
  }
  hTot ->SetDirectory(nullptr);
  hRej ->SetDirectory(nullptr);
  pPass->SetDirectory(nullptr);
  fin->Close();

  const int nX = hTot->GetNbinsX();
  const int nY = hTot->GetNbinsY();

  /* ================================================================
   * 2) Reject‑fraction 2‑D map  (with binomial errors)
   * ================================================================
   */
  TH2D hRejFrac("h_chi2_rejFrac_etaPhi",
                "Reject fraction after #chi^{2} cut;"
                "Tower #eta index;Tower #varphi index;Reject fraction",
                nX,hTot->GetXaxis()->GetXmin(),hTot->GetXaxis()->GetXmax(),
                nY,hTot->GetYaxis()->GetXmin(),hTot->GetYaxis()->GetXmax());

  for(int ix=1;ix<=nX;++ix)
    for(int iy=1;iy<=nY;++iy){
        const double tot = hTot->GetBinContent(ix,iy);
        const double rej = hRej->GetBinContent(ix,iy);
        const double frac = (tot>0)? rej/tot : 0.;
        const double err  = (tot>0)
            ? TEfficiency::ClopperPearson((int)tot,(int)rej,0.683,false)-frac
            : 0.;
        hRejFrac.SetBinContent(ix,iy,frac);
        hRejFrac.SetBinError  (ix,iy,err );
    }

  /* ================================================================
   * 3) 1‑D projections  (reject & pass fractions)
   * ================================================================
   */
  auto projFrac = [&](bool alongEta, bool wantPass)->TH1D*
  {
      const char* axis  = alongEta? "eta":"phi";
      const char* kind  = wantPass? "Pass":"Rej";
      const int   nBins = alongEta? nX:nY;
      TH1D* h=new TH1D(Form("hChi2_%sFrac_vs%s",kind,axis),
                       Form("%s fraction vs tower #%s index",
                            wantPass? "Pass":"Reject",axis),
                       nBins,0,nBins);

      for(int i=1;i<=nBins;++i){
          double tot=0,acc=0;
          for(int j=1;j<=(alongEta?nY:nX);++j){
              int ix = alongEta? i:j;
              int iy = alongEta? j:i;
              const double t = hTot->GetBinContent(ix,iy);
              const double r = hRej->GetBinContent(ix,iy);
              const double p = pPass->GetBinContent(ix,iy);
              tot+=t;
              acc+= wantPass? (t*p) : r;
          }
          double frac=(tot>0)? acc/tot:0.;
          double err =(tot>0)? std::sqrt(frac*(1-frac)/tot):0.;
          h->SetBinContent(i,frac);
          h->SetBinError  (i,err );
      }
      h->GetYaxis()->SetTitle(Form("%s fraction / tower",wantPass?"Pass":"Reject"));
      h->GetXaxis()->SetTitle(Form("Tower #%s index",axis));
      return h;
  };

  TH1D* hRejVsEta = projFrac(true ,false);
  TH1D* hRejVsPhi = projFrac(false,false);
  TH1D* hPassVsEta= projFrac(true ,true );
  TH1D* hPassVsPhi= projFrac(false,true);

  /* ================================================================
   * 4)  list of worst towers
   * ================================================================
   */
  struct TowerStat{int iEta,iPhi;double frac,tot;};
  std::vector<TowerStat> v;
  v.reserve(nX*nY);
  for(int ix=1;ix<=nX;++ix)
    for(int iy=1;iy<=nY;++iy){
        const double tot=hTot->GetBinContent(ix,iy);
        if(tot<1) continue;
        v.push_back({ix-1,iy-1,hRejFrac.GetBinContent(ix,iy),tot});
    }
  std::sort(v.begin(),v.end(),[](auto&a,auto&b){return a.frac>b.frac;});
  std::cout<<"\n[Chi2QA]  Worst "<<nWorstToPrint<<" towers by reject fraction\n"
           <<"  η  φ   frac   (total)\n"
           <<"  -- --- ------ -------\n";
  for(int i=0;i<std::min(nWorstToPrint,(int)v.size());++i){
     const auto&t=v[i];
     std::cout<<std::setw(3)<<t.iEta<<std::setw(4)<<t.iPhi
              <<std::setw(8)<<std::fixed<<std::setprecision(3)<<t.frac
              <<" ("<<t.tot<<")\n";
  }

  /* ================================================================
   * 5)  drawing helpers  (no grids, centred title)
   * ================================================================
   */
  auto drawSave2D=[&](TH2* h,const char* cname,const char* png,
                      double zMin=0,double zMax=0){
      if(zMin<zMax){h->SetMinimum(zMin);h->SetMaximum(zMax);}
      TCanvas c(cname,cname,1000,760);
      c.SetRightMargin(0.15);

      h->GetXaxis()->SetTitleOffset(1.1);
      h->GetYaxis()->SetTitleOffset(1.1);
      h->GetZaxis()->SetTitleFont(42);
      h->GetZaxis()->SetTitleSize(0.045);

      if(std::string(h->GetName()).find("Pass")!=std::string::npos)
          h->GetZaxis()->SetTitle("Pass fraction");
      else if(std::string(h->GetName()).find("Frac")!=std::string::npos)
          h->GetZaxis()->SetTitle("Reject fraction");
      else
          h->GetZaxis()->SetTitle("Counts");

      h->Draw("COLZ");

      c.SaveAs(Form("%s/%s",outDir,png));
  };

  auto drawSave1D=[&](TH1* h,const char* cname,const char* png){
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerSize(0.8);
      TCanvas c(cname,cname,1000,600);
      h->Draw("E1");

      c.SaveAs(Form("%s/%s",outDir,png));
  };

  /* ================================================================
   * 6)  plots
   * ================================================================
   */
  drawSave2D(hTot       ,"cChi2Tot" ,"Chi2_Total.png"          );
  drawSave2D(hRej       ,"cChi2Rej" ,"Chi2_Rejected.png"       );
  drawSave2D(pPass      ,"cChi2Pass","Chi2_PassFraction.png",0,1);
  drawSave2D(&hRejFrac  ,"cChi2Frac","Chi2_RejectFraction.png",0,1);

  drawSave1D(hRejVsEta ,"cRejVsEta" ,"Chi2_RejFrac_vsEta.png" );
  drawSave1D(hRejVsPhi ,"cRejVsPhi" ,"Chi2_RejFrac_vsPhi.png" );
  drawSave1D(hPassVsEta,"cPassVsEta","Chi2_PassFrac_vsEta.png");
  drawSave1D(hPassVsPhi,"cPassVsPhi","Chi2_PassFrac_vsPhi.png");

  /* ================================================================
   * 7)  global numbers
   * ================================================================
   */
  const double tot=hTot->GetEntries();
  const double rej=hRej->GetEntries();
  std::cout<<"[Chi2QA]  GLOBAL: total="<<tot
           <<"  rejected="<<rej
           <<"  frac="<<(tot?rej/tot:0)<<"\n";

  /* ================================================================
   * 8)  save helper ROOT file
   * ================================================================
   */
  std::unique_ptr<TFile> fqa(TFile::Open(
      (std::string(inFile)+"_Chi2QA.root").c_str(),"RECREATE"));
  hTot->Write(); hRej->Write(); pPass->Write();
  hRejFrac.Write();
  hRejVsEta->Write(); hRejVsPhi->Write();
  hPassVsEta->Write();hPassVsPhi->Write();
  delete hRejVsEta; delete hRejVsPhi;
  delete hPassVsEta;delete hPassVsPhi;
  fqa->Close();
  std::cout<<"[Chi2QA]  QA histograms saved.\n";
}


// ---------------------------------------------------------------------------
//  Truth-vertex-Z spectrum  *with verbose QA*
//  • Drawn as blue full circles
//  • Y-axis is auto–scaled to the histogram maximum
// ---------------------------------------------------------------------------
void PlotVertexZTruthOnly(TH1F* h_truth_vz,
                          const char* outDir)
/* ------------------------------------------------------------------------ */
{
  /* ───────────────────────── sanity check ──────────────────────────── */
  if (!h_truth_vz)
  {
    std::cerr << "\n[vtxZ]  ERROR  null pointer:  h_truth_vz=" << h_truth_vz << '\n';
    return;
  }

  /* ───────────────────────── detach from any file ──────────────────── */
  h_truth_vz->SetDirectory(nullptr);

  /* ───────────────────────── global stats ──────────────────────────── */
  const int    firstBin   = h_truth_vz->FindFirstBinAbove(0.0);
  const int    lastBin    = h_truth_vz->FindLastBinAbove (0.0);
  const double integral   = h_truth_vz->Integral();

  std::cout << "[vtxZ] Truth  Entries=" << h_truth_vz->GetEntries()
            << "  Integral="  << integral
            << "  Mean="      << h_truth_vz->GetMean()
            << "  RMS="       << h_truth_vz->GetRMS()
            << "  FirstBin="  << firstBin
            << "  LastBin="   << lastBin
            << '\n';

  /* ───────────────────────── style ─────────────────────────────────── */
  h_truth_vz->SetMarkerStyle(20);
  h_truth_vz->SetMarkerSize (0.8);
  h_truth_vz->SetMarkerColor(kBlue+1);
  h_truth_vz->SetLineColor  (kBlue+1);
  h_truth_vz->SetLineWidth  (2);

  /* ───────────────────────── y-axis range ──────────────────────────── */
  const double yMax = 1.30 * h_truth_vz->GetMaximum();
  std::cout << "[vtxZ] Y-axis will span 0 … " << yMax << '\n';

  h_truth_vz->GetYaxis()->SetRangeUser(0.0, yMax);
  h_truth_vz->GetYaxis()->SetTitle("Counts");
  h_truth_vz->GetXaxis()->SetTitle("z_{vtx}  (cm)");

  /* ───────────────────────── draw & save ───────────────────────────── */
  TCanvas cV("VertexZ_Truth","Truth vertex-Z", 880, 620);
  gStyle->SetOptStat(0);

  h_truth_vz->Draw("E1");      // blue filled circles

  TString outName = Form("%s/VertexZ_Truth.png", outDir);
  cV.SaveAs(outName);

  std::cout << "[vtxZ] wrote " << outName
            << "  (Truth max=" << h_truth_vz->GetMaximum() << ")\n\n";
}


namespace
{
// ─────────────────────────────────────────────────────────────────────────────
//  Formatting constants (shared by φ & η)
// ─────────────────────────────────────────────────────────────────────────────
constexpr std::array<Color_t,5>  kCol { kGreen+2, kBlue, kBlack, kRed, kMagenta+1 };
constexpr std::array<Style_t,5>  kMk  { 20,        20,        20,     20,    20   };
constexpr std::array<const char*,5> kLab{
        "no corr, scratch",
        "b(E) corr, scratch",
        "no corr, cluster",
        "CorrectPosition, cluster",
        "b(E) corr, cluster"
};

// ─────────────────────────────────────────────────────────────────────────────
//  Generic legend helper – one definition for all canvases
// ─────────────────────────────────────────────────────────────────────────────
struct LegendPos   { double x1{0.16}, y1{0.67}, x2{0.45}, y2{0.90}; };

inline TLegend*
addVariantLegend(const std::vector<int>& shown,
                 const LegendPos& lp             = {},
                 const char*      header         = nullptr,
                 double           txtSize        = 0.028 )
{
    auto* lg = new TLegend(lp.x1, lp.y1, lp.x2, lp.y2);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(txtSize);
    if(header) lg->SetHeader(header,"C");

    for(int v : shown)
    {
        auto* mk = new TMarker(0,0,kMk[v]);
        mk->SetMarkerColor(kCol[v]);
        lg->AddEntry(mk,kLab[v],"p");
    }
    lg->Draw();
    return lg;                // caller keeps ownership
}


// ─────────────────────────────────────────────────────────────────────────────
//  Small helpers
// ─────────────────────────────────────────────────────────────────────────────
inline void ensureDir(const std::string& p)
{
    gSystem->mkdir(p.c_str(), /*recurse=*/true);
}

struct FitRes { double mu{}, dmu{}, sg{}, dsg{}; };

/// Robust one‑dim Gaussian fit (quartile pre‑seed, IRQ window)
inline FitRes robustGaussianFit(TH1* h)
{
    if(!h || h->Integral()==0) return {};

    const double q[3] = {0.25,0.5,0.75};
    double quart[3]; h->GetQuantiles(3,quart,const_cast<double*>(q));
    const double med = quart[1];
    const double iqr = quart[2]-quart[0];

    double lo = med-2.5*iqr , hi = med+2.5*iqr;
    if(lo>=hi){ lo=h->GetMean()-2*h->GetRMS(); hi=h->GetMean()+2*h->GetRMS(); }

    TF1 f("g","gaus",lo,hi);
    f.SetParameters(h->GetMaximum(),med,0.5*h->GetRMS());
    const auto r = h->Fit(&f,"QNR0S");
    if(!r.Get() || !r->IsValid())
    {
        const double rms = h->GetRMS();
        const double err = rms/std::sqrt(2.*std::max(1.,h->GetEntries()-1.));
        return { h->GetMean(), err, rms, err };
    }
    return { f.GetParameter(1), f.GetParError(1),
             std::fabs(f.GetParameter(2)), f.GetParError(2) };
}

/// Clone + normalise histogram (probability density)
inline TH1F* cloneAndNormPdf(const TH1F* src, int v)
{
    auto* h = static_cast<TH1F*>(src->Clone());
    h->SetDirectory(nullptr);
    if(h->Integral()>0) h->Scale(1./h->Integral());

    h->SetMarkerColor(kCol[v]); h->SetLineColor(kCol[v]);
    h->SetMarkerStyle(kMk[v]);  h->SetMarkerSize(0.9);
    return h;
}



// ─────────────────────────────────────────────────────────────────────────────
//  One‑slice overlay
// ─────────────────────────────────────────────────────────────────────────────
struct SliceResult
{
    std::array<FitRes,5> fit{};
    double yMax = 0.0;
};

SliceResult drawSlice(const std::array<std::vector<TH1F*>,5>& h,
                      int                           iE,
                      const std::vector<int>&       shown,
                      const std::pair<double,double>& eEdge,
                      TCanvas&                      cSlice,
                      const char*                   deltaSym)
{
    SliceResult out{};
    cSlice.cd();

    for(int v : shown)
    {
        if(!h[v][iE]) continue;
        auto* hh  = cloneAndNormPdf(h[v][iE], v);
        out.fit[v]= robustGaussianFit(hh);
        double localMax = 0.0;
        for (int b = 1; b <= hh->GetNbinsX(); ++b) {
            const double y = hh->GetBinContent(b) + hh->GetBinError(b);
            if (y > localMax) localMax = y;
        }
        out.yMax = std::max(out.yMax, localMax);

        hh->SetTitle(Form("%s  [%g, %g)  GeV", deltaSym,
                          eEdge.first, eEdge.second));
        hh->GetXaxis()->SetTitle(Form("%s = reco #minus truth  [rad]", deltaSym));
        hh->GetYaxis()->SetTitle("Probability density");
        hh->GetYaxis()->SetRangeUser(0,1.25*out.yMax);
        hh->Draw( (v==shown.front()) ? "E" : "E SAME" );
    }
    const double ymax = 1.10 * out.yMax;       // 10 % head‑room
    for (auto *pr : *cSlice.GetListOfPrimitives()) {
        if (auto *h = dynamic_cast<TH1*>(pr))
            h->GetYaxis()->SetRangeUser(0., ymax);
    }
    cSlice.Modified();     // force pad to recompute axes
    cSlice.Update();

    addVariantLegend(shown, LegendPos{});
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Energy‑trend “summary sheet”   (μ & σ stacked pads)
// ─────────────────────────────────────────────────────────────────────────────
void drawMuSigmaSummary(const std::string&     fileStem,
                        const std::vector<double>&                eCtr,
                        const std::array<std::vector<double>,5>&  MU,
                        const std::array<std::vector<double>,5>&  dMU,
                        const std::array<std::vector<double>,5>&  SG,
                        const std::array<std::vector<double>,5>&  dSG,
                        const std::vector<int>&                   shown)
{
    const double halfBinLo = (eCtr.size()>1) ? 0.5*(eCtr[1]-eCtr[0]) : 1.0;
    const double halfBinHi = (eCtr.size()>1) ? 0.5*(eCtr.back()-eCtr[eCtr.size()-2]) : 1.0;
    const double xAxisMin  = 0.0;
    const double xAxisMax  = eCtr.back() + halfBinHi;

    auto frameRange=[&](bool wantMu)
    {
        double lo=1e30, hi=-1e30;
        for(int v:shown)
            for(size_t i=0;i<eCtr.size();++i)
            {
                const auto& val = wantMu ? MU[v] : SG[v];
                const auto& err = wantMu ? dMU[v] : dSG[v];
                lo = std::min(lo, val[i]-err[i]);
                hi = std::max(hi, val[i]+err[i]);
            }
        const double pad = 0.2*(hi-lo);
        return std::pair{lo-pad, hi+pad};
    };

    auto makeGraph=[&](const std::vector<double>& Y,
                       const std::vector<double>& dY,
                       int v)
    {
        std::vector<double> ex(eCtr.size(),0.);
        return new TGraphErrors((int)eCtr.size(),
                                eCtr.data(), Y.data(),
                                ex.data(),  dY.data());
    };

    TCanvas c("cMuSig","μ & σ vs E",900,800);
    c.Divide(1,2,0,0);

    // --- μ pad ----------------------------------------------------------
    auto padU = static_cast<TPad*>(c.cd(1));
    padU->SetLeftMargin(0.15); padU->SetBottomMargin(0.04);

    const auto [muLo,muHi] = frameRange(true);
    TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
    frU.SetMinimum(muLo); frU.SetMaximum(muHi); frU.Draw("AXIS");
    TLine l0(xAxisMin,0,xAxisMax,0); l0.SetLineStyle(2); l0.Draw();

    for(int v:shown){
        auto* g = makeGraph(MU[v],dMU[v],v);
        g->SetMarkerStyle(kMk[v]); g->SetMarkerColor(kCol[v]);
        g->SetLineColor  (kCol[v]); g->SetMarkerSize(1.1);
        g->Draw("P SAME");
    }

    padU->cd();   // make sure we are inside the μ-pad
    addVariantLegend(shown,
                     LegendPos{0.75, 0.72, 0.90, 0.92},   // upper-right corner
                     nullptr,                            // no legend header
                     0.035);

    // --- σ pad ----------------------------------------------------------
    auto padL = static_cast<TPad*>(c.cd(2));
    padL->SetLeftMargin(0.15); padL->SetTopMargin(0.06); padL->SetBottomMargin(0.35);

    const auto [sgLo,sgHi] = frameRange(false);
    TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
    frL.SetMinimum(std::max(0.,sgLo)); frL.SetMaximum(sgHi); frL.Draw("AXIS");

    for(int v:shown){
        auto* g = makeGraph(SG[v],dSG[v],v);
        g->SetMarkerStyle(kMk[v]); g->SetMarkerColor(kCol[v]);
        g->SetLineColor(kCol[v]);  g->SetMarkerSize(1.1);
        g->Draw("P SAME");
    }

    const std::string outPng =
        std::string(fileStem) + "/" +
        gSystem->BaseName(fileStem.c_str()) + "_MeanSigmaVsE.png";
    c.SaveAs(outPng.c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
//  RMSE, σ/|μ| and width‑ratio helpers
// ─────────────────────────────────────────────────────────────────────────────
std::array<std::vector<double>,5> buildFromFormula(
        const std::array<std::vector<double>,5>& MU,
        const std::array<std::vector<double>,5>& SG,
        const std::function<double(double,double)>& f)
{
    std::array<std::vector<double>,5> res{};
    for(int v=0; v<5; ++v)
    {
        res[v].resize(MU[v].size());
        for(size_t i=0;i<MU[v].size();++i)
            res[v][i] = f(MU[v][i],SG[v][i]);
    }
    return res;
}

void drawDiag(const std::string&                 outPng,
              const std::string&                 yTitle,
              const std::vector<double>&         eCtr,
              const std::array<std::vector<double>,5>& Y,
              const std::vector<int>&            shown)
{
    const double xMin = eCtr.front(), xMax = eCtr.back();

    double yLo=1e30,yHi=-1e30;
    for(int v:shown)
        for(double y:Y[v]){ yLo=std::min(yLo,y); yHi=std::max(yHi,y); }
    const double pad=0.07*(yHi-yLo); yLo-=pad; yHi+=pad;
    if(yLo<0) yLo=0;

    TCanvas c(outPng.c_str(),outPng.c_str(),900,640);
    c.SetLeftMargin(0.15); c.SetRightMargin(0.06);

    const std::string hTitle = ";E [GeV];" + yTitle;
    TH1F fr("fr",hTitle.c_str(),1,xMin,xMax);
    fr.SetMinimum(yLo); fr.SetMaximum(yHi); fr.Draw("AXIS");

    std::vector<double> ex(eCtr.size(),0.);
    for(int v:shown)
    {
        auto* g = new TGraphErrors((int)eCtr.size(),
                                   eCtr.data(),Y[v].data(),ex.data(),ex.data());
        g->SetMarkerStyle(kMk[v]); g->SetMarkerColor(kCol[v]);
        g->SetLineColor(kCol[v]);  g->SetMarkerSize(1.1);
        g->Draw("P SAME");
    }
    addVariantLegend(shown, LegendPos{}, nullptr, 0.025);
    c.SaveAs(outPng.c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
//  μ vs lnE diagnostic (shared by φ & η)
// ─────────────────────────────────────────────────────────────────────────────
void drawMuVsLnE(const std::string&                 baseDir,
                 const std::vector<double>&         eCtr,
                 const std::array<std::vector<double>,5>& MU,
                 const std::array<std::vector<double>,5>& dMU,
                 const std::vector<int>&            shown,
                 const char*                         fileStem)
{
    if(eCtr.size()<2) return;

    std::vector<double> lnE(eCtr.size());
    std::transform(eCtr.begin(),eCtr.end(),lnE.begin(),
                   [](double e){return std::log(e);} );

    const double xLo = *std::min_element(lnE.begin(),lnE.end()) - 0.05;
    const double xHi = *std::max_element(lnE.begin(),lnE.end()) + 0.05;

    double yLo=1e30,yHi=-1e30;
    for(int v:shown){
        yLo=std::min(yLo,*std::min_element(MU[v].begin(),MU[v].end()));
        yHi=std::max(yHi,*std::max_element(MU[v].begin(),MU[v].end()));
    }
    const double pad=0.15*(yHi-yLo); yLo-=pad; yHi+=0.35*(yHi-yLo);

    TCanvas c("cLn","μ vs lnE",900,640);
    c.SetLeftMargin(0.15); c.SetRightMargin(0.06);

    TH1F fr("",";ln E  [GeV];#mu  [rad]",1,xLo,xHi);
    fr.SetMinimum(yLo); fr.SetMaximum(yHi); fr.Draw("AXIS");

    std::ofstream fout(baseDir + std::string("/") + fileStem + "_MuVsLogE_fit.txt");
    fout<<"# variant   p0(rad)   p1(rad)\n"<<std::scientific<<std::setprecision(6);

    TLegend lg(0.16,0.75,0.90,0.90);
    lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.02);

    std::vector<double> ex(eCtr.size(), 0.);
    for (int v : shown)
    {
        auto* g = new TGraphErrors((int)eCtr.size(), lnE.data(), MU[v].data(),
                                   ex.data(), dMU[v].data());
        g->SetMarkerStyle(kMk[v]);
        g->SetMarkerColor(kCol[v]);
        g->SetLineColor  (kCol[v]);
        g->SetMarkerSize(1.1);
        g->Draw("P SAME");

        TF1 f("f", "pol1", xLo, xHi);
        f.SetLineColor (kCol[v]);
        f.SetLineStyle (2);
        f.SetLineWidth (2);
        g->Fit(&f, "Q");
        f.Draw("SAME");

        fout << std::left << std::setw(28) << kLab[v]
             << f.GetParameter(0) << "   " << f.GetParameter(1) << '\n';

        auto* mk = new TMarker(0, 0, kMk[v]);
        mk->SetMarkerColor(kCol[v]);
        lg.AddEntry(mk,
                    Form("%s  (a=%.2e, b=%.2e)", kLab[v],
                         f.GetParameter(0), f.GetParameter(1)),
                    "p");
    }
    addVariantLegend(shown, LegendPos{}, nullptr, 0.02);

    /* ------------------------------------------
     *  finally write the PNG next to the *.txt
     * ------------------------------------------ */
    const std::string outPng =
        baseDir + std::string("/") + fileStem + std::string("_MuVsLogE.png");
    c.SaveAs(outPng.c_str());
    c.Close();
}

// ─────────────────────────────────────────────────────────────────────────────
//  Generic residual‑plot driver (φ or η)
// ─────────────────────────────────────────────────────────────────────────────
void makeResidualPlots(const std::vector<std::pair<double,double>>& eEdges,
                       EBinningMode                                binMode,
                       const std::array<std::vector<TH1F*>,5>&    h,
                       const char*                                 outDirRaw,
                       const std::vector<int>&                     variantsToPlot)
{
    if(eEdges.empty()){ std::cerr<<"[MakeResidualPlots] eEdges empty – abort\n"; return; }

    const std::string outDir(outDirRaw);
    ensureDir(outDir);

    /* ------------------------------------------------------------------ *
     *  Decide once which residual we handle by looking at the directory
     *  path that the wrapper provides:
     *                        “…/deltaPhi”  →  Δφ
     *                        “…/deltaEta”  →  Δη
     * ------------------------------------------------------------------ */
    const bool  isPhi     = (outDir.find("deltaPhi") != std::string::npos);
    const char* deltaSym = isPhi ? "#Delta#phi" : "#Delta#eta";
    const char* fileStem = isPhi ? "DeltaPhi"   : "DeltaEta";

    // ─── containers for summary graphs ─────────────────────────────────
    const int N_E = static_cast<int>(eEdges.size());
    std::array<std::vector<double>,5> MU,dMU, SG,dSG;
    for(int v=0; v<5; ++v){ MU[v].reserve(N_E); dMU[v].reserve(N_E);
                            SG[v].reserve(N_E); dSG[v].reserve(N_E); }

    std::vector<double> eCtr; eCtr.reserve(N_E);

    // ─── per‑slice overlays ────────────────────────────────────────────
    TCanvas cSlice("cSlice","slice",900,700);
    const std::string sliceDir = outDir + "/" + deltaSym + "_slices";
    ensureDir(sliceDir);

    /* optional z‑vertex label, e.g.  "|z_{vtx}| ∈ [0, 5) cm" */
    std::string vtxLabel;
    {
        const std::string probe = "/vertexDependent/vz";
        auto pos = outDir.find(probe);
        if(pos != std::string::npos){
            std::string tag = outDir.substr(pos + probe.size());  // "0_5"
            float vzLo = 0.f, vzHi = 0.f;
            if(std::sscanf(tag.c_str(),"%f_%f",&vzLo,&vzHi) == 2){
                vtxLabel = Form("|z_{vtx}| #in [%.0f, %.0f) cm", vzLo, vzHi);
            }
        }
    }
    for(int iE=0;iE<N_E;++iE)
    {
        // -------------------------------------------------------------
        // 1) overlaid slice (unchanged)
        // -------------------------------------------------------------
        eCtr.push_back(0.5*(eEdges[iE].first+eEdges[iE].second));

        const auto res = drawSlice(h,iE,variantsToPlot,eEdges[iE],
                                   cSlice,deltaSym);

        if(!vtxLabel.empty()){
            cSlice.cd();
            TLatex tz;
            tz.SetNDC(); tz.SetTextFont(42); tz.SetTextSize(0.028);
            tz.SetTextAlign(11);               // left‑top corner
            tz.DrawLatex(0.19,0.9, vtxLabel.c_str());
        }

        cSlice.SaveAs(Form("%s/%s_%02d.png",sliceDir.c_str(),fileStem,iE));
        cSlice.Clear();

        // collect fit results for the summary graphs
        for(int v:variantsToPlot){
            MU[v].push_back(res.fit[v].mu );  dMU[v].push_back(res.fit[v].dmu);
            SG[v].push_back(res.fit[v].sg );  dSG[v].push_back(res.fit[v].dsg);
        }

        // -------------------------------------------------------------
        // 2) NEW: individual components in the same colours
        //      one sub‑folder per energy‑bin
        // -------------------------------------------------------------
        const std::string gridDir =
              Form("%s/%s_%02d_components", sliceDir.c_str(), fileStem, iE);
        ensureDir(gridDir);

        TCanvas cGrid("cGrid","components_grid",1200,800);
        cGrid.Divide(3,2,0,0);          // 3 columns × 2 rows

        std::array<FitRes,5> fitVar{};  // remember fit results
        int padIdx = 1;

        for(int v : variantsToPlot)
        {
            if(!h[v][iE]) { ++padIdx; continue; }

            cGrid.cd(padIdx++);

            auto*   hh = cloneAndNormPdf(h[v][iE], v);
            fitVar[v] = robustGaussianFit(hh);

            // dynamic y‑range that includes the error bars
            double localMax = 0.0;
            for(int b=1;b<=hh->GetNbinsX();++b){
                const double y = hh->GetBinContent(b) + hh->GetBinError(b);
                if(y>localMax) localMax = y;
            }

            hh->SetTitle(kLab[v]);
            hh->GetXaxis()->SetTitle(
                   Form("%s = reco #minus truth  [rad]", deltaSym));
            hh->GetYaxis()->SetTitle("Probability density");
            hh->GetYaxis()->SetRangeUser(0,1.10*localMax);
            hh->Draw("E");

            // annotate energy‑bin
            TLatex txt;
            txt.SetNDC(); txt.SetTextFont(42);
            txt.SetTextSize(0.035); txt.SetTextAlign(13);
            txt.SetTextColor(kCol[v]);
            txt.DrawLatex(0.15,0.86,
                Form("[%.1f, %.1f) GeV",
                     eEdges[iE].first, eEdges[iE].second));
            if(!vtxLabel.empty())
                txt.DrawLatex(0.15,0.78, vtxLabel.c_str());
        }

        // ---------- summary pad (bottom‑right) ------------------------
        cGrid.cd(6);
        gPad->SetLeftMargin(0.02); gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.02);  gPad->SetBottomMargin(0.02);
        gPad->Clear();

        TLatex info;
        info.SetNDC(); info.SetTextFont(42);
        info.SetTextSize(0.032); info.SetTextAlign(13);

        double yPos = 0.92, dy = 0.14;
        for(int v : variantsToPlot)
        {
            info.SetTextColor(kCol[v]);
            info.DrawLatex(0.045, yPos,
                Form("%s :  #mu = %.3g #pm %.3g   #sigma_{RMS} = %.3g #pm %.3g",
                     kLab[v],
                     fitVar[v].mu,  fitVar[v].dmu,
                     fitVar[v].sg,  fitVar[v].dsg));
            yPos -= dy;
        }

        cGrid.SaveAs(
            Form("%s/%s_%02d_grid.png",
                 gridDir.c_str(), fileStem, iE));
        cGrid.Close();
    }

    // ─── μ & σ summary ────────────────────────────────────────────────
    drawMuSigmaSummary(outDir, eCtr, MU,dMU, SG,dSG, variantsToPlot);

    // ─── RMSE & σ/|μ| diagnostics ─────────────────────────────────────
    const auto RMSE = buildFromFormula(MU,SG,
                       [](double mu,double sg){return std::hypot(mu,sg);} );
    const auto Frac = buildFromFormula(MU,SG,
                       [](double mu,double sg){return sg/std::max(1e-9,std::fabs(mu));} );

    drawDiag(outDir + "/" + fileStem + "_RMSE.png",
             "#sqrt{#mu^{2}+#sigma_{RMS}^{2}}  [rad]", eCtr, RMSE, variantsToPlot);

    drawDiag(outDir + "/" + fileStem + "_FracRes.png",
             "#sigma_{RMS} / |#mu|", eCtr, Frac, variantsToPlot);

    // ─── width ratio (vs v=2) ─────────────────────────────────────────
    const int ref=2;
    std::array<std::vector<double>,5> WR{};
    for(int v:variantsToPlot)
    {
        WR[v].resize(SG[v].size());
        for(size_t i=0;i<SG[v].size();++i)
            WR[v][i] = SG[ref][i]>0 ? SG[v][i]/SG[ref][i] : 0.;
    }
    drawDiag(outDir + "/" + fileStem + "_WidthRatioVsE.png",
             "#sigma_{RMS}_{variant}/#sigma_{RMS}_{baseline}", eCtr, WR, variantsToPlot);

    // ─── μ vs lnE ──────────────────────────────────────────────────────
    drawMuVsLnE(outDir, eCtr, MU,dMU, variantsToPlot, fileStem);

    std::cout << "[MakeResidualPlots] wrote plots to  " << outDir << '\n';
}

} // end anonymous namespace

// ============================================================================
//  PUBLIC WRAPPERS  (φ, η  and combined φ + η)  ––– VERBOSE / DEBUG EDITION
// ============================================================================

// ————————————————————————————————————————————————————————————————
//  Helpers & state needed for the combined plots
// ————————————————————————————————————————————————————————————————
namespace
{
    //-------------------------------------------------------------------
    //  Utility: coloured banner to make log‑messages stand out
    //-------------------------------------------------------------------
    inline void dbgBanner(const std::string& head,
                          const std::string& msg,
                          char fill = '=')
    {
        const std::size_t W = 78;
        std::string line(W, fill);
        std::cout << "\n" << line << "\n"
                  << "[" << head << "] " << msg << "\n"
                  << line << "\n";
    }

    //-------------------------------------------------------------------
    //  Global scratch storage – filled by MakeDeltaPhiPlots(…)
    //-------------------------------------------------------------------
    static std::vector<std::pair<double,double>>     gEdgePhi;
    static std::array<std::vector<TH1F*>,5>          gHPhiGlobal{};
    static EBinningMode                              gBinModePhi  = EBinningMode::kRange;
    static std::vector<int>                          gVariantsPhi = {0,1,2,3,4};

    //-------------------------------------------------------------------
    //  Sanity helper: check that vectors are well‑formed
    //-------------------------------------------------------------------
    bool checkViewIntegrity(const std::array<std::vector<TH1F*>,5>& view,
                            const char* tag)
    {
        bool ok = true;
        for (int v = 0; v < 5; ++v)
        {
            for (std::size_t i = 0; i < view[v].size(); ++i)
            {
                if (!view[v][i])
                {
                    std::cerr << "[WARN] <" << tag << ">  v=" << v
                              << "  iE=" << i << "  → nullptr histogram\n";
                    ok = false;
                }
            }
        }
        return ok;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Produce Δη + Δφ overlays
    //    • one PNG per energy slice   (unchanged)
    //    • one 2×4 “table” per variant (NEW)
    // ─────────────────────────────────────────────────────────────────────────────
    void makeCombinedDeltaEtaPhiPlots(
            const std::vector<std::pair<double,double>>& eEdges,
            EBinningMode                                 /*binMode*/,
            const std::array<std::vector<TH1F*>,5>&      hEta,
            const std::array<std::vector<TH1F*>,5>&      hPhi,
            const char*                                  outDir,
            const std::vector<int>&                      variantsToPlot)
    {
        dbgBanner("makeCombinedDeltaEtaPhiPlots",
                  Form("Entering – eEdges=%zu  outDir=%s", eEdges.size(), outDir));

        if (eEdges.empty()) {
            std::cerr << "[ERROR] eEdges vector is empty – nothing to plot.\n";
            return;
        }

        checkViewIntegrity(hEta,"η‑view");
        checkViewIntegrity(hPhi,"φ‑view");
        ensureDir(outDir);

        const char* fileStem = "DeltaEtaPhi";

        // ─────────────────────────────────────────────────────────────────
        //  PART 1 – per‑energy‑slice grids  (fixed legend & labels)
        // ─────────────────────────────────────────────────────────────────
        for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
        {
            TCanvas cGrid(Form("cGrid_%zu", iE), "combined_eta_phi", 1200, 800);
            cGrid.Divide(3, 2, 0, 0);                      // 3 × 2 layout

            int padIdx = 1;
            for (int v : variantsToPlot)
            {
                if (!hEta[v][iE] || !hPhi[v][iE]) { ++padIdx; continue; }
                cGrid.cd(padIdx++);

                // ───────────────────────────────────────────────────
                // 1) draw the two PDFs
                // ───────────────────────────────────────────────────
                auto* hEtaPDF = cloneAndNormPdf(hEta[v][iE], v);
                auto* hPhiPDF = cloneAndNormPdf(hPhi[v][iE], v);
                hPhiPDF->SetMarkerStyle(kMk[v] + 4);   // open symbol for Δφ
                hPhiPDF->SetLineStyle(2);              // dashed line for Δφ

                const double yMax = 1.10 *
                    std::max(hEtaPDF->GetMaximum(), hPhiPDF->GetMaximum());

                hEtaPDF->SetTitle("");                 // labels are added manually
                hEtaPDF->GetXaxis()->SetTitle("#Delta  [rad]");
                hEtaPDF->GetYaxis()->SetTitle("Probability density");
                hEtaPDF->GetYaxis()->SetRangeUser(0., yMax);

                hEtaPDF->Draw("E");                    // first:  Δη
                hPhiPDF->Draw("E SAME");               // then:   Δφ

                // ───────────────────────────────────────────────────
                // 2) statistical numbers (robust Gaussian fit)
                // ───────────────────────────────────────────────────
                const FitRes fEta = robustGaussianFit(hEta[v][iE]);
                const FitRes fPhi = robustGaussianFit(hPhi[v][iE]);

                // anchor for all TLatex / legend in pad coordinates
                const double xLeft = gPad->GetLeftMargin() + 0.04;

                // ───────────────────────────────────────────────────
                // 3) legend (bottom‑right) – tells which marker is which residual
                //     • solid  =  Δη
                //     • open   =  Δφ
                // ───────────────────────────────────────────────────
                const double padR   = 1.0 - gPad->GetRightMargin() - 0.02;  // 2 % in from pad edge
                const double legW   = 0.26;                                 // legend width
                const double xLegL  = padR - legW;                          // left x of legend

                /*  keep the legend safely above the pad’s bottom‑margin:
                    0.28 works for both the top row (small margin) and the bottom
                    row (large margin required for X–axis labels).                   */
                const double yLegB  = 0.28;                 // bottom y of legend
                const double yLegT  = yLegB + 0.16;         // ~16 % tall

                auto *leg = new TLegend(xLegL, yLegB, padR, yLegT);
                leg->SetBorderSize(0);          // no frame
                leg->SetFillStyle(0);           // transparent
                leg->SetTextSize(0.03);

                // dummy markers that mimic the real styles
                auto *mEta = new TMarker(0, 0, kMk[v]);     // solid   – Δη
                mEta->SetMarkerColor(kCol[v]); mEta->SetMarkerSize(1.1);

                auto *mPhi = new TMarker(0, 0, kMk[v] + 4); // open    – Δφ
                mPhi->SetMarkerColor(kCol[v]); mPhi->SetMarkerSize(1.1);

                leg->AddEntry(mEta, "#Delta#eta", "p");
                leg->AddEntry(mPhi, "#Delta#phi", "p");
                leg->Draw("same");              // keep inside pad’s primitive list

                gPad->Modified();               // guarantee proper layout

                // ───────────────────────────────────────────────────
                // 4) variant label – coloured, just below the legend
                // ───────────────────────────────────────────────────
                TLatex vLab;  vLab.SetNDC();  vLab.SetTextColor(kCol[v]);
                vLab.SetTextSize(0.034);  vLab.SetTextAlign(11);
                vLab.DrawLatex(xLeft, 0.71, kLab[v]);

                // ───────────────────────────────────────────────────
                // 5) Gaussian μ ± σ (same x‑anchor as legend)
                // ───────────────────────────────────────────────────
                TLatex stats;  stats.SetNDC();  stats.SetTextSize(0.035);
                stats.SetTextAlign(11);
                stats.DrawLatex(xLeft, 0.64,
                    Form("#Delta#eta : %.2g #pm %.2g", fEta.mu, fEta.sg));
                stats.DrawLatex(xLeft, 0.58,
                    Form("#Delta#phi : %.2g #pm %.2g", fPhi.mu, fPhi.sg));

                // ───────────────────────────────────────────────────
                // 6) energy‑range label (upper‑right corner)
                // ───────────────────────────────────────────────────
                TLatex eLab;  eLab.SetNDC();  eLab.SetTextSize(0.04);
                eLab.SetTextAlign(33);
                eLab.DrawLatex(0.88, 0.92,
                    Form("[%.0f, %.0f) GeV",
                         eEdges[iE].first, eEdges[iE].second));

                // ------------------------------------------------------------------
                // 6)  save this pad as an **individual PNG**
                //       <outDir>/<variant>/<fileStem>_<variant>_E<slice>.png
                //      (define the variant‑specific folder locally)
                // ------------------------------------------------------------------
                std::string vNameLocal = kLab[v];                    // human‑readable label
                for (char& c : vNameLocal)                           // sanitise for file‑system
                    if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';

                std::string vDirLocal = std::string(outDir) + "/" + vNameLocal;
                ensureDir(vDirLocal);                                // create if missing

                gPad->Modified();                                    // finalise layout
                gPad->Update();
                gPad->SaveAs(Form("%s/%s_%s_E%02zu.png",
                                  vDirLocal.c_str(),                 // …/variant/…
                                  fileStem,                          // e.g. DeltaEtaPhi
                                  vNameLocal.c_str(),                // sanitised variant label
                                  iE));                              // energy‑slice index
            }

            cGrid.cd(6);               // sixth pad kept blank
            gPad->Clear();

            cGrid.SaveAs(Form("%s/%s_%02zu_grid.png", outDir, fileStem, iE));
            cGrid.Close();
        }


        // ─────────────────────────────────────────────────────────────────
        //  PART 2 – per‑variant 2×4 tables  (NEW)
        // ─────────────────────────────────────────────────────────────────
        const int nCol = 4, nRow = 2, maxPads = nCol * nRow;

        for (int v : variantsToPlot)
        {
            // create folder “…/combinedDeltaEtaPhi/<variant>/”
            std::string vName  = kLab[v];                  // human‑readable label
            for (char& c : vName)                          // sanitise for FS path
                if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';

            std::string vDir = std::string(outDir) + "/" + vName;
            ensureDir(vDir);

            TCanvas cTab(Form("cTab_%d", v), kLab[v], 1600, 900);
            cTab.SetTopMargin(0.12);           // gives the title room to breathe
            cTab.Divide(nCol, nRow, 0, 0);

            int pad = 1;
            for (std::size_t iE = 0; iE < eEdges.size() && pad <= maxPads; ++iE)
            {
                if (!hEta[v][iE] || !hPhi[v][iE]) continue;

                cTab.cd(pad++);                      // go to next pad

                // ────────────────────────────────────────────────
                // 1) prepare & draw the PDFs
                // ────────────────────────────────────────────────
                auto* hEtaPDF = cloneAndNormPdf(hEta[v][iE], v);
                auto* hPhiPDF = cloneAndNormPdf(hPhi[v][iE], v);
                hPhiPDF->SetMarkerStyle(kMk[v] + 4);       // open symbol for Δφ
                hPhiPDF->SetLineStyle(2);                  // dashed line for Δφ

                const double yMax = 1.10 *
                    std::max(hEtaPDF->GetMaximum(), hPhiPDF->GetMaximum());

                hEtaPDF->SetTitle("");                     // we write labels manually
                // ─── X axis ──────────────────────────────────────────────────────────
                hEtaPDF->GetXaxis()->SetTitle("#Delta  [rad]");
                TAxis* axX = hEtaPDF->GetXaxis();
                axX->SetLabelSize(0.03);      // smaller tick‑label font
                axX->SetTitleSize(0.035);     // keep title readable
                axX->SetTitleOffset(1.2);     // lift title a bit
                axX->SetNdivisions(505);      // 5 major ticks, no minor → avoids clutter

                // ─── Y axis ──────────────────────────────────────────────────────────
                hEtaPDF->GetYaxis()->SetTitle("Probability density");
                TAxis* axY = hEtaPDF->GetYaxis();
                axY->SetLabelSize(0.03);      // match Y‑axis label size to X
                axY->SetTitleSize(0.035);
                axY->SetTitleOffset(1.5);

                hEtaPDF->GetYaxis()->SetRangeUser(0., yMax);

                hEtaPDF->Draw("E");                        // first:  Δη
                hPhiPDF->Draw("E SAME");                   // then:   Δφ

                // ────────────────────────────────────────────────
                // 2) statistics (robust Gaussian fit)
                // ────────────────────────────────────────────────
                const FitRes fEta = robustGaussianFit(hEta[v][iE]);
                const FitRes fPhi = robustGaussianFit(hPhi[v][iE]);

                // ────────────────────────────────────────────────
                // 3) legend   (upper‑left) – x‑position adapts to pad margin
                //     • solid symbol  →  Δη
                //     • open  symbol  →  Δφ
                // ────────────────────────────────────────────────
                const double xLeft   = gPad->GetLeftMargin() + 0.04;   // stay clear of Y‑axis
                const double legW    = 0.28;                           // fixed legend width
                const double yTop    = 0.92;
                const double yBottom = 0.78;

                auto* leg = new TLegend(xLeft, yBottom, xLeft + legW, yTop);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.05);

                /* create dummy markers that carry the exact style we used */
                auto* mEta = new TMarker(0, 0, kMk[v]);       // solid symbol = Δη
                mEta->SetMarkerColor(kCol[v]);
                mEta->SetMarkerSize(1.2);

                auto* mPhi = new TMarker(0, 0, kMk[v] + 4);   // open symbol  = Δφ
                mPhi->SetMarkerColor(kCol[v]);
                mPhi->SetMarkerSize(1.2);

                leg->AddEntry(mEta, "#Delta#eta", "p");
                leg->AddEntry(mPhi, "#Delta#phi", "p");
                leg->Draw();      // the pad now owns the legend – it will stay until SaveAs


                // ────────────────────────────────────────────────
                // 4) μ ± σ numbers just below the legend – same x anchor
                // ────────────────────────────────────────────────
                TLatex stats;  stats.SetNDC();  stats.SetTextSize(0.038);  stats.SetTextAlign(11);
                stats.DrawLatex(xLeft, 0.45,
                    Form("#Delta#eta : %.2g #pm %.2g", fEta.mu, fEta.sg));
                stats.DrawLatex(xLeft, 0.38,
                    Form("#Delta#phi : %.2g #pm %.2g", fPhi.mu, fPhi.sg));

                // ────────────────────────────────────────────────
                // 5) energy‑range label   (upper‑right)
                // ────────────────────────────────────────────────
                TLatex eLab;  eLab.SetNDC();  eLab.SetTextSize(0.055);  eLab.SetTextAlign(33);
                eLab.DrawLatex(0.88, 0.92,
                    Form("[%.0f, %.0f) GeV",
                         eEdges[iE].first, eEdges[iE].second));

            }

            /* ------------------------------------------------------------------
             *  (A)  finish the 2×4 “component” table  – identical to old code
             * ------------------------------------------------------------------ */
            cTab.cd(0);
            TLatex head; head.SetNDC(); head.SetTextSize(0.042);
            head.SetTextAlign(22);

            const char* tabTitle =
                (!strcmp(kLab[v],"CorrectPosition, cluster"))
                    ? "Coresoftware Code w/ CorrectPosition"
                : (!strcmp(kLab[v],"b(E) corr, scratch"))
                    ? "From Scratch Code w/ Position Correction"
                    : kLab[v];

            head.DrawLatex(0.5,0.97, tabTitle);
            cTab.SaveAs(Form("%s/%s_%s.png",
                             vDir.c_str(), fileStem, vName.c_str()));
            cTab.Close();

            /* ------------------------------------------------------------------
             *  (B)  NEW –  μ(E)  &  σ(E)  summary for this *single* variant
             *        – top pad:  μ vs E (Δη solid, Δφ open)
             *        – bottom :  σ vs E (same symbols, dashed for φ)
             * ------------------------------------------------------------------ */

            /* gather the per–slice fit results that we just produced */
            std::vector<double>  eCtr, muEta, dmuEta, sgEta, dsgEta,
                                 muPhi, dmuPhi, sgPhi, dsgPhi;
            for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
            {
                if (!hEta[v][iE] || !hPhi[v][iE]) continue;

                const double eMid = 0.5 * (eEdges[iE].first + eEdges[iE].second);

                const FitRes fEta = robustGaussianFit(hEta[v][iE]);
                const FitRes fPhi = robustGaussianFit(hPhi[v][iE]);

                eCtr .push_back(eMid);

                muEta.push_back(fEta.mu);  dmuEta.push_back(fEta.dmu);
                sgEta.push_back(fEta.sg);  dsgEta.push_back(fEta.dsg);

                muPhi.push_back(fPhi.mu);  dmuPhi.push_back(fPhi.dmu);
                sgPhi.push_back(fPhi.sg);  dsgPhi.push_back(fPhi.dsg);
            }

            if (eCtr.empty())   // nothing to draw → skip
                continue;

            /* ---------- cosmetics shared by both pads ------------------- */
            const double xMin = 0.0;
            const double xMax = eCtr.back() +
                                0.5*(eEdges.back().second - eEdges.back().first);
            std::vector<double> ex(eCtr.size(), 0.0);

            /* ---------- create the canvas -------------------------------- */
            TCanvas cVS(Form("cVS_%d", v), "mu_sigma_vs_E", 900, 800);
            cVS.Divide(1, 2, 0, 0);

            /* ===== top pad :  μ(E) ====================================== */
            cVS.cd(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.04);
            gPad->SetTopMargin(0.10);          // extra head-room for the canvas title

            const double muLo = std::min( *std::min_element(muEta.begin(), muEta.end()),
                                          *std::min_element(muPhi.begin(), muPhi.end()) );
            const double muHi = std::max( *std::max_element(muEta.begin(), muEta.end()),
                                          *std::max_element(muPhi.begin(), muPhi.end()) );
            const double padMu = 0.2 * (muHi - muLo);

            TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMin,xMax);
            frU.SetMinimum(muLo - padMu);  frU.SetMaximum(muHi + padMu);
            frU.Draw("AXIS");
            
            
            /*  dashed y = 0 reference (μ-pad only) */
            TLine ref0(xMin, 0.0, xMax, 0.0);
            ref0.SetLineStyle(2);   // dashed
            ref0.SetLineWidth(1);
            ref0.Draw();


            TGraphErrors gMuEta(eCtr.size(), eCtr.data(), muEta.data(),
                                ex.data(), dmuEta.data());
            gMuEta.SetMarkerStyle(kMk[v]);       gMuEta.SetMarkerColor(kCol[v]);
            gMuEta.SetLineColor  (kCol[v]);      gMuEta.Draw("P SAME");

            TGraphErrors gMuPhi(eCtr.size(), eCtr.data(), muPhi.data(),
                                ex.data(), dmuPhi.data());
            gMuPhi.SetMarkerStyle(kMk[v]+4);     gMuPhi.SetMarkerColor(kCol[v]);
            gMuPhi.SetLineColor  (kCol[v]);      gMuPhi.SetLineStyle(2);
            gMuPhi.Draw("P SAME");

            /* legend – show only on μ-pad, top-right corner */
            TLegend legU(0.88,0.74,0.97,0.9);  // x1,y1,x2,y2  (NDC)
            legU.SetBorderSize(0);  legU.SetFillStyle(0); legU.SetTextSize(0.048);
            legU.AddEntry(&gMuEta,"#Delta#eta","p");
            legU.AddEntry(&gMuPhi,"#Delta#phi","p");
            legU.Draw();

            /* ===== bottom pad :  σ(E) =================================== */
            cVS.cd(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetTopMargin(0.06);          // keep room for x-axis title

            const double sgHi = std::max( *std::max_element(sgEta.begin(), sgEta.end()),
                                          *std::max_element(sgPhi.begin(), sgPhi.end()) );

            TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xMin,xMax);
            frL.SetMinimum(0.0);           frL.SetMaximum(1.15 * sgHi);
            frL.Draw("AXIS");

            TGraphErrors gSgEta(eCtr.size(), eCtr.data(), sgEta.data(),
                                ex.data(), dsgEta.data());
            gSgEta.SetMarkerStyle(kMk[v]);       gSgEta.SetMarkerColor(kCol[v]);
            gSgEta.SetLineColor  (kCol[v]);      gSgEta.Draw("P SAME");

            TGraphErrors gSgPhi(eCtr.size(), eCtr.data(), sgPhi.data(),
                                ex.data(), dsgPhi.data());
            gSgPhi.SetMarkerStyle(kMk[v]+4);     gSgPhi.SetMarkerColor(kCol[v]);
            gSgPhi.SetLineColor  (kCol[v]);      gSgPhi.SetLineStyle(2);
            gSgPhi.Draw("P SAME");               //  ← no legend here


            cVS.cd(0);
            TLatex hSum;  hSum.SetNDC();  hSum.SetTextAlign(22);  hSum.SetTextSize(0.036);

            std::string sumTitle =
                (!strcmp(kLab[v],"CorrectPosition, cluster"))
                    ? "Coresoftware Code w/ CorrectPosition  #mu/#sigma_{RMS} Summary"
                : (!strcmp(kLab[v],"b(E) corr, scratch"))
                    ? "From Scratch Code w/ Position Correction  #mu/#sigma_{RMS} Summary"
                    : std::string(kLab[v]) + "  #mu/#sigma_{RMS} Summary";
            hSum.DrawLatex(0.5,0.98, sumTitle.c_str());   // slightly higher & smaller

            cVS.SaveAs(Form("%s/%s_%s_MuSigmaVsE.png",
                            vDir.c_str(), fileStem, vName.c_str()));
            cVS.Close();
        }

        /* ------------------------------------------------------------------
         *  PART 3 –  NEW comparison folders
         *
         *      Two‑variant overlays
         *      ─────────────────────
         *          • compareClusterizerWithWithout
         *                2  (“no corr, cluster”)
         *                3  (“CorrectPosition, cluster”)
         *          • compareFromScratchWithWithout
         *                0  (“no corr, scratch”)
         *                1  (“b(E) corr, scratch”)
         *          • compareCPfromScratchVsClusterizer
         *                3  (“CorrectPosition, cluster”)
         *                1  (“b(E) corr, scratch”)
         *          • compareRawfromScratchVsClusterizer
         *                0  (“no corr, scratch”)
         *                2  (“no corr, cluster”)
         *
         *      Four‑variant overlay
         *      ────────────────────
         *          • fourWayCompareCorrectUncorr
         *                0, 1, 2, 3  (all four baseline / corrected flavours)
         *
         *      For every folder we create, **per residual type (Δφ & Δη)**:
         *          – a 2 × 4 table (or 4‑way table) for all energy slices
         *          – a matching μ(E) & σ(E) summary sheet
         * ------------------------------------------------------------------ */
        struct CompSpec
        {
            const char*            subdir;
            std::vector<int>       v;         // any length ≥2
        };
        const std::vector<CompSpec> comps = {
            { "compareClusterizerWithWithout",    {2,3} },
            { "compareFromScratchWithWithout",    {0,1} },
            { "compareCPfromScratchVsClusterizer",{3,1} },
            { "compareRawfromScratchVsClusterizer",{0,2} },
            { "fourWayCompareCorrectUncorr",      {0,1,2,3} }
        };

        const int nColC   = 4;
        const int nRowC   = 2;
        const int maxPadsC= nColC*nRowC;           // fits all 8 energy slices

        constexpr Color_t kCmpCol[2]={ kBlue+1 , kRed+1 };

        auto markerFor=[&](std::size_t /*vIdx*/,std::size_t /*which*/){
            return 20;                      // filled circle for every curve
        };
    
        for (const auto& C : comps)
        {
            const std::string baseDir = std::string(outDir) + "/" + C.subdir;
            ensureDir(baseDir);

            for (int resType = 0; resType < 2; ++resType)   // 0 = Δφ , 1 = Δη
            {
                const bool      isPhi  = (resType == 0);
                const char*     resTag = isPhi ? "DeltaPhi" : "DeltaEta";
                const char*     sym    = isPhi ? "#Delta#phi" : "#Delta#eta";
                const auto&     hView  = isPhi ? hPhi : hEta;

                TCanvas cTbl(Form("cTbl_%s_%s", resTag, C.subdir), "cmp_table", 1600, 900);
                cTbl.Divide(nColC, nRowC, 0, 0);

                /* give every pad a roomier top margin */
                for(int ip=1; ip<=nColC*nRowC; ++ip){
                    cTbl.cd(ip);
                    gPad->SetTopMargin(0.15);      // was the default 0.05
                }

                int padId = 1;
                std::vector<double> eCtr;
                std::vector<std::vector<double>> MU, dMU, SG, dSG;
                const std::size_t Nvar = C.v.size();
                MU .resize(Nvar); dMU.resize(Nvar);
                SG .resize(Nvar); dSG.resize(Nvar);

                for (std::size_t iE = 0; iE < eEdges.size() && padId <= maxPadsC; ++iE)
                {
                    bool haveAll = true;
                    for (int vidx : C.v)
                        if (!hView[vidx][iE]) { haveAll=false; break; }
                    if (!haveAll) continue;

                    cTbl.cd(padId++);

                    double localYmax = 0.0;
                    for (std::size_t j = 0; j < Nvar; ++j)
                    {
                        const int vIdx = C.v[j];
                        auto* h = cloneAndNormPdf(hView[vIdx][iE], vIdx);
                        if (j>0){                        // use open symbol & dashed line
                            h->SetMarkerStyle(20);
                            h->SetMarkerColor(kCmpCol[j<2 ? j : 0]);
                            h->SetLineColor  (kCmpCol[j<2 ? j : 0]);
                            h->SetLineStyle(2);
                        }
                        localYmax = std::max(localYmax, h->GetMaximum());
                        h->GetXaxis()->SetTitle(Form("%s  [rad]", sym));
                        h->GetYaxis()->SetTitle("Probability density");
                        h->GetYaxis()->SetRangeUser(0.,1.30*localYmax);
                        h->Draw( (j==0) ? "E" : "E SAME" );

                        const FitRes fr = robustGaussianFit(hView[vIdx][iE]);
                        MU [j].push_back(fr.mu ); dMU[j].push_back(fr.dmu);
                        SG [j].push_back(fr.sg ); dSG[j].push_back(fr.dsg);
                    }

                    /* energy‑range label, one per pad */
                    TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.04); eLab.SetTextAlign(33);
                    eLab.DrawLatex(0.88,0.92,
                        Form("[%.0f, %.0f) GeV", eEdges[iE].first, eEdges[iE].second));

                    /* legend */
                    TLegend lg(0.58,0.72,0.92,0.92);
                    lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.03);
                    for (std::size_t j=0;j<Nvar;++j){
                        const int vIdx=C.v[j];
                        auto* m=new TMarker(0,0,markerFor(vIdx,j));
                        m->SetMarkerColor(kCmpCol[j<2 ? j : 0]); m->SetMarkerSize(1.1);
                        lg.AddEntry(m,kLab[vIdx],"p");
                    }
                    lg.Draw();

                    if (eCtr.size() < iE+1)             // push eCtr only once per slice
                        eCtr.push_back(0.5*(eEdges[iE].first+eEdges[iE].second));
                }

                /* title & write PNG */
                cTbl.cd(0);
                TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextSize(0.045);

                std::string tMain;
                if      (!strcmp(C.subdir,"compareClusterizerWithWithout"))
                    tMain = Form("%s  Coresoftware Code With/Without CorrectPosition", sym);
                else if (!strcmp(C.subdir,"compareFromScratchWithWithout"))
                    tMain = Form("%s  From Scratch Code With/Without Position Correction", sym);
                else if (!strcmp(C.subdir,"compareCPfromScratchVsClusterizer"))
                    tMain = Form("%s  From Scratch Code vs Clusterizer With Position Correction", sym);
                else if (!strcmp(C.subdir,"compareRawfromScratchVsClusterizer"))
                    tMain = Form("%s  From Scratch Code vs Clusterizer Without Position Correction", sym);
                else
                    tMain = Form("%s, %zu way comparison", sym, Nvar);

                head.DrawLatex(0.5,0.96, tMain.c_str());
                cTbl.SaveAs(Form("%s/%s_Table.png", baseDir.c_str(), resTag));
                cTbl.Close();

                /* --- (B) μ(E) & σ(E) summary sheet --------------------------- */
                if (eCtr.empty()) continue;
                const double xMax = eCtr.back() + 0.5*(eEdges.back().second - eEdges.back().first);
                std::vector<double> ex(eCtr.size(),0.0);

                TCanvas cSum(Form("cSum_%s_%s", resTag, C.subdir), "cmp_summary", 900, 800);
                cSum.Divide(1,2,0,0);

                /* μ pad ------------------------------------------------------- */
                cSum.cd(1);
                gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.04);
                double muLow=1e30, muHigh=-1e30;
                for (std::size_t j=0;j<Nvar;++j)
                    for (std::size_t k=0;k<MU[j].size();++k){
                        muLow = std::min(muLow, MU[j][k]-dMU[j][k]);
                        muHigh= std::max(muHigh,MU[j][k]+dMU[j][k]);
                    }
                TH1F frU("frU",";E [GeV];#mu [rad]",1,0.0,xMax);
                frU.SetMinimum(muLow-0.1*(muHigh-muLow));
                frU.SetMaximum(muHigh+0.1*(muHigh-muLow));
                frU.Draw("AXIS");
                TLine ref0(0.0,0.0,xMax,0.0); ref0.SetLineStyle(2); ref0.Draw();

                TLegend legU(0.18,0.75,0.42,0.90);
                legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.03);

                for (std::size_t j=0;j<Nvar;++j){
                    int vIdx=C.v[j];
                    TGraphErrors* g=new TGraphErrors(eCtr.size(),
                                                     eCtr.data(), MU[j].data(),
                                                     ex.data(),  dMU[j].data());
                    g->SetMarkerStyle(20);
                    g->SetMarkerColor(kCmpCol[j<2 ? j : 0]);
                    g->SetLineColor  (kCmpCol[j<2 ? j : 0]);
                    if (j>0) g->SetLineStyle(2);
                    g->Draw("P SAME");
                    legU.AddEntry(g,kLab[vIdx],"p");
                }
                legU.Draw();

                /* σ pad ------------------------------------------------------- */
                cSum.cd(2);
                gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);
                double sgHigh=-1e30;
                for (std::size_t j=0;j<Nvar;++j)
                    for (double vsg:SG[j]) sgHigh=std::max(sgHigh,vsg);
                TH1F frL("frL",";E [GeV];#sigma_{RMS} [rad]",1,0.0,xMax);
                frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHigh); frL.Draw("AXIS");

                for (std::size_t j=0;j<Nvar;++j){
                    int vIdx=C.v[j];
                    TGraphErrors* g=new TGraphErrors(eCtr.size(),
                                                     eCtr.data(), SG[j].data(),
                                                     ex.data(),  dSG[j].data());
                    g->SetMarkerStyle(20);                            // filled circle
                    g->SetMarkerColor(kCmpCol[j < 2 ? j : 0]);       // blue (j=0) or red (j=1)
                    g->SetLineColor  (kCmpCol[j < 2 ? j : 0]);
                    if (j>0) g->SetLineStyle(2);                     // dashed for the 2nd curve
                    g->Draw("P SAME");
                }

                /* add summary title line */
                cSum.cd(0);
                TLatex hS;  hS.SetNDC(); hS.SetTextAlign(22); hS.SetTextSize(0.045);

                std::string tSum = tMain + "  Gaussian #mu/#sigma_{RMS} Summary";
                hS.DrawLatex(0.5,0.97, tSum.c_str());

                cSum.SaveAs(Form("%s/%s_MuSigmaVsE.png", baseDir.c_str(), resTag));
                cSum.Close();
            }   // ── end loop over residual type (Δφ / Δη)
            
            /* ------------------------------------------------------------
             *  EXTRA summary for the 4-way overlay:  draw *both* residuals
             *  (Δη filled, Δφ open) on one μ-pad and one σ-pad.
             *  Only produced for the “fourWayCompareCorrectUncorr” folder.
             * ------------------------------------------------------------ */
            if (!strcmp(C.subdir,"fourWayCompareCorrectUncorr"))
            {
                /* ----- gather fit results for all variants & both residuals --- */
                std::vector<double> eCtr;
                std::array<std::vector<double>,4> MUe, MUf, sGe, sGf, dMUe, dMUf, dsGe, dsGf;

                for (std::size_t iE=0;iE<eEdges.size();++iE)
                {
                    bool ok=true;
                    for (int vIdx : {0,1,2,3})
                        if (!hEta[vIdx][iE] || !hPhi[vIdx][iE]) { ok=false; break; }
                    if (!ok) continue;

                    const double eMid=0.5*(eEdges[iE].first+eEdges[iE].second);
                    eCtr.push_back(eMid);

                    for (int j=0;j<4;++j)
                    {
                        const FitRes fE = robustGaussianFit(hEta[j][iE]);
                        const FitRes fF = robustGaussianFit(hPhi[j][iE]);

                        MUe[j].push_back(fE.mu);  dMUe[j].push_back(fE.dmu);
                        sGe[j].push_back(fE.sg);  dsGe[j].push_back(fE.dsg);

                        MUf[j].push_back(fF.mu);  dMUf[j].push_back(fF.dmu);
                        sGf[j].push_back(fF.sg);  dsGf[j].push_back(fF.dsg);
                    }
                }
                if (!eCtr.empty())
                {
                    const double xMax=eCtr.back()+0.5*(eEdges.back().second-eEdges.back().first);
                    std::vector<double> ex(eCtr.size(),0.0);

                    TCanvas c4("cSum_BothRes_fourWay","combined_mu_sigma",900,800);
                    c4.Divide(1,2,0,0);

                    /* === μ-pad ============================================ */
                    c4.cd(1);
                    gPad->SetLeftMargin(0.15);   // keep generous room for the y-axis
                    gPad->SetRightMargin(0.02);  // slim buffer on the right
                    gPad->SetBottomMargin(0.06); // a little more room for the x-axis ticks
                    gPad->SetTopMargin(0.12);    // ↑ extra head-room for title & legend

                    double muLo= 1e30, muHi=-1e30;
                    for(int j=0;j<4;++j)
                        for(std::size_t k=0;k<MUe[j].size();++k){
                            muLo=std::min(muLo,std::min(MUe[j][k]-dMUe[j][k],
                                                        MUf[j][k]-dMUf[j][k]));
                            muHi=std::max(muHi,std::max(MUe[j][k]+dMUe[j][k],
                                                        MUf[j][k]+dMUf[j][k]));
                        }
                    TH1F frU("frU",";E [GeV];#mu [rad]",1,0.0,xMax);
                    /* widen the vertical range: 25 % padding and force the axis to include y=0 */
                    const double padFrac = 0.25;                      // 25 % extra space
                    double yMin = muLo - padFrac * (muHi - muLo);
                    double yMax = muHi + padFrac * (muHi - muLo);
                    if (yMin > 0.0)  yMin = 0.0;                      // keep 0 within view
                    if (yMax < 0.0)  yMax = 0.0;
                    frU.SetMinimum(yMin);
                    frU.SetMaximum(yMax);

                    frU.Draw("AXIS");
                    auto* ref0 = new TLine(0.0, 0.0, xMax, 0.0);
                    ref0->SetLineStyle(2);          // dashed
                    ref0->SetLineColor(kBlack);     // black
                    ref0->Draw();

                    /* ---------- colour & marker conventions we enforce ---------------- */
                    static const Color_t varCol[4]   = { kBlue+1 , kBlue+1 , kRed+1  , kRed+1  };
                    static const Style_t openCirc    = 24 ,  closedCirc  = 20;
                    static const Style_t openSquare  = 25 ,  closedSquare= 21;

                    auto chooseStyle = [&](int j, bool isPhi){
                        const bool wantSquare = (j==1 || j==3);           // “with correction” → square
                        return wantSquare ? (isPhi ? closedSquare : openSquare)
                                          : (isPhi ? closedCirc   : openCirc );
                    };

                    /* ----- FIRST: draw all data points (μ-pad already current) -------- */
                    static const int ord[4] = {0,2,1,3};                  // preferred display order
                    for(int idx = 0; idx < 4; ++idx)
                    {
                        const int j = ord[idx];

                        /* Δη (open marker) */
                        auto* gEta = new TGraphErrors((int)eCtr.size(),
                                                      eCtr.data(), MUe[j].data(),
                                                      ex.data(),  dMUe[j].data());
                        gEta->SetMarkerStyle( chooseStyle(j,false) );
                        gEta->SetMarkerColor( varCol[j] );
                        gEta->SetLineColor  ( varCol[j] );
                        gEta->SetMarkerSize(0.8);
                        gEta->Draw("P SAME");

                        /* Δφ (closed marker) */
                        auto* gPhi = new TGraphErrors((int)eCtr.size(),
                                                      eCtr.data(), MUf[j].data(),
                                                      ex.data(),  dMUf[j].data());
                        gPhi->SetMarkerStyle( chooseStyle(j,true) );
                        gPhi->SetMarkerColor( varCol[j] );
                        gPhi->SetLineColor  ( varCol[j] );
                        gPhi->SetMarkerSize(0.8);
                        gPhi->Draw("P SAME");
                    }

                    /* ---------- legend:  open-icon / closed-icon  text -------------- */
                    TLegend legU(0.7,0.65,0.98,0.875);           // upper-right corner
                    legU.SetBorderSize(0);
                    legU.SetFillStyle(0);
                    legU.SetTextSize(0.034);

                    legU.SetNColumns(4);                          // ∘ | / | ● | text
                    legU.SetColumnSeparation(0.002);              // bring icons & slash closer
                    legU.SetMargin(0.02);

                    for (int idx = 0; idx < 4; ++idx)
                    {
                        const int j = ord[idx];                   // display order 0-2-1-3

                        auto* mkEta = new TMarker(0,0, chooseStyle(j,false));   // open symbol (η)
                        auto* mkPhi = new TMarker(0,0, chooseStyle(j,true));    // filled symbol (φ)
                        mkEta->SetMarkerColor(varCol[j]);
                        mkPhi->SetMarkerColor(varCol[j]);

                        const char* dummy = " ";        // non‑empty → ROOT prints no class‑name

                        /* col‑0 : open icon (η) */
                        legU.AddEntry(mkEta, dummy, "p");

                        /* col‑1 : slash separator */
                        legU.AddEntry((TObject*)nullptr, " / ", "");

                        /* col‑2 : closed icon (φ) */
                        legU.AddEntry(mkPhi, dummy, "p");

                        /* col-3 : descriptive text */
                        const char* lbl =
                            (j==0)? "#Delta#eta/#Delta#phi  From Scratch, no correction" :
                            (j==2)? "#Delta#eta/#Delta#phi  Coresoftware, no correction" :
                            (j==1)? "#Delta#eta/#Delta#phi  From Scratch, with correction" :
                                    "#Delta#eta/#Delta#phi  Coresoftware, with correction";
                        legU.AddEntry((TObject*)nullptr, lbl, "");
                    }
                    legU.Draw();


                    /* === σ-pad ============================================ */
                    c4.cd(2);
                    gPad->SetLeftMargin(0.15);
                    gPad->SetTopMargin(0.06);          // keep room for axis title

                    double sgHi = -1e30;
                    for(int j=0;j<4;++j){
                        for(double v : sGe[j]) sgHi = std::max(sgHi,v);
                        for(double v : sGf[j]) sgHi = std::max(sgHi,v);
                    }

                    TH1F frL("frL",";E [GeV];#sigma_{RMS} [rad]",1,0.0,xMax);
                    frL.SetMinimum(0.0);
                    frL.SetMaximum(1.15*sgHi);
                    frL.Draw("AXIS");

                    for(int j=0;j<4;++j){
                        /* Δη : open marker --------------------------------- */
                        TGraphErrors* gSe = new TGraphErrors(eCtr.size(),
                                                             eCtr.data(), sGe[j].data(),
                                                             ex.data(),  dsGe[j].data());
                        gSe->SetMarkerStyle( chooseStyle(j,false) );   // same rule as μ-pad
                        gSe->SetMarkerColor( varCol[j] );
                        gSe->SetLineColor  ( varCol[j] );
                        gSe->SetMarkerSize(0.8);                       // smaller markers
                        gSe->Draw("P SAME");

                        /* Δφ : closed marker ------------------------------- */
                        TGraphErrors* gSf = new TGraphErrors(eCtr.size(),
                                                             eCtr.data(), sGf[j].data(),
                                                             ex.data(),  dsGf[j].data());
                        gSf->SetMarkerStyle( chooseStyle(j,true) );
                        gSf->SetMarkerColor( varCol[j] );
                        gSf->SetLineColor  ( varCol[j] );
                        gSf->SetMarkerSize(0.8);
                        gSf->Draw("P SAME");
                    }

                    /* canvas title – full 4‑way overlay --------------------------------------- */
                    c4.cd(0);
                    TLatex tC; tC.SetNDC(); tC.SetTextSize(0.032); tC.SetTextAlign(22);
                    tC.DrawLatex(0.5,0.97,
                        "#Delta#eta / #Delta#phi, From Scratch/Coresoftware Comparison Gaussian #mu/#sigma_{RMS} Summary");

                    c4.SaveAs(Form("%s/DeltaEtaPhi_Combined_MuSigmaVsE.png", baseDir.c_str()));
                    c4.Close();
            
                    
                    /* -------------------------------------------------------------------------
                     *  ADDITIONAL 2‑way overlays – uncorrected (0,2)   and   corrected (1,3)
                     *  Each canvas shows 4 points per slice (η open, φ filled).
                     * ------------------------------------------------------------------------- */
                    struct DuoSpec { const char* tag; std::array<int,2> v; };
                    const std::array<DuoSpec,2> duos = {
                        DuoSpec{ "Uncorrected", {0,2} },     // “no corr, scratch”  +  “no corr, cluster”
                        DuoSpec{ "Corrected",   {1,3} }      // “b(E) corr, scratch”+“CorrectPosition, cluster”
                    };
                    for (const DuoSpec& D : duos)
                    {
                        /* --- create a fresh canvas ----------------------------------------- */
                        TCanvas c2(Form("c2_%s",D.tag),"mu_sigma_duo",900,800);
                        c2.Divide(1,2,0,0);

                        /* ---------- helper lambdas ---------------------------------------- */
                        auto getLoHi = [&](bool wantMu){
                            double lo=1e30, hi=-1e30;
                            for(int idx=0;idx<2;++idx){
                                int j=D.v[idx];
                                const auto& V  = wantMu ? MUe[j] : sGe[j];
                                const auto& dV = wantMu ? dMUe[j]: dsGe[j];
                                const auto& Vf = wantMu ? MUf[j] : sGf[j];
                                const auto& dVf= wantMu ? dMUf[j]: dsGf[j];
                                for(size_t k=0;k<V.size();++k){
                                    lo=std::min(lo,std::min(V[k]-dV[k],Vf[k]-dVf[k]));
                                    hi=std::max(hi,std::max(V[k]+dV[k],Vf[k]+dVf[k]));
                                }
                            }
                            return std::pair{lo,hi};
                        };
                        auto drawGraphs=[&](bool wantMu){
                            for(int idx=0;idx<2;++idx){
                                int j=D.v[idx];
                                /* η → open marker ------------------------------------------ */
                                TGraphErrors* gE=new TGraphErrors(eCtr.size(),
                                                                  eCtr.data(),
                                                                  (wantMu?MUe[j]:sGe[j]).data(),
                                                                  ex.data(),
                                                                  (wantMu?dMUe[j]:dsGe[j]).data());
                                gE->SetMarkerStyle( chooseStyle(j,false) );
                                gE->SetMarkerColor( varCol[j] ); gE->SetLineColor( varCol[j] );
                                gE->SetMarkerSize(0.9); gE->Draw("P SAME");

                                /* φ → filled marker ---------------------------------------- */
                                TGraphErrors* gF=new TGraphErrors(eCtr.size(),
                                                                  eCtr.data(),
                                                                  (wantMu?MUf[j]:sGf[j]).data(),
                                                                  ex.data(),
                                                                  (wantMu?dMUf[j]:dsGf[j]).data());
                                gF->SetMarkerStyle( chooseStyle(j,true) );
                                gF->SetMarkerColor( varCol[j] ); gF->SetLineColor( varCol[j] );
                                gF->SetMarkerSize(0.9); gF->Draw("P SAME");
                            }
                        };

                        /* === μ‑pad ======================================================== */
                        c2.cd(1);
                        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.05); gPad->SetTopMargin(0.12);
                        auto [muLo,muHi]=getLoHi(true);
                        TH1F frU("frU",";E [GeV];#mu [rad]",1,0.0,xMax);
                        frU.SetMinimum(muLo-0.1*(muHi-muLo)); frU.SetMaximum(muHi+0.1*(muHi-muLo));
                        frU.Draw("AXIS");
                        auto* ref0 = new TLine(0.0, 0.0, xMax, 0.0);
                        ref0->SetLineStyle(2);
                        ref0->SetLineColor(kBlack);
                        ref0->Draw();
                        drawGraphs(true);

                        /* legend – two variants, slash separator --------------------------- */
                        TLegend lg(0.75,0.78,0.99,0.89); lg.SetBorderSize(0); lg.SetFillStyle(0);
                        lg.SetTextSize(0.035);
                        lg.SetNColumns(4);                          // ∘ | / | ● | text
                        lg.SetColumnSeparation(0.002);              // bring icons & slash closer
                        lg.SetMargin(0.02);
                        
                        lg.SetTextSize(0.028); lg.SetNColumns(4); lg.SetColumnSeparation(0.002);
                        for(int idx=0;idx<2;++idx){
                            int j=D.v[idx];
                            auto* mkE = new TMarker(0,0, chooseStyle(j,false));   // open  (η)
                            auto* mkF = new TMarker(0,0, chooseStyle(j,true));    // filled (φ)
                            mkE->SetMarkerColor(varCol[j]);
                            mkF->SetMarkerColor(varCol[j]);

                            const char* dummy = " ";

                            lg.AddEntry(mkE, dummy, "p");               // η  icon
                            lg.AddEntry((TObject*)nullptr," / ","");    // slash
                            lg.AddEntry(mkF, dummy, "p");               // φ  icon

                            const char* txt =
                                (j==0)? "#Delta#eta/#Delta#phi  From Scratch" :
                                (j==2)? "#Delta#eta/#Delta#phi  Coresoftware"  :
                                (j==1)? "#Delta#eta/#Delta#phi  From Scratch"  :
                                        "#Delta#eta/#Delta#phi  Coresoftware";
                            lg.AddEntry((TObject*)nullptr,txt,"");
                        }
                        lg.Draw();

                        /* === σ‑pad ======================================================== */
                        c2.cd(2);
                        gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);
                        auto [sgLo,sgHi]=getLoHi(false);
                        TH1F frL("frL",";E [GeV];#sigma_{RMS} [rad]",1,0.0,xMax);
                        frL.SetMinimum(std::max(0.0,sgLo-0.05*sgHi)); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");
                        drawGraphs(false);

                        /* title & save ----------------------------------------------------- */
                        c2.cd(0);
                        TLatex tit; tit.SetNDC(); tit.SetTextSize(0.032); tit.SetTextAlign(22);
                        tit.DrawLatex(0.5,0.97,
                            Form("#Delta#eta / #Delta#phi  %s Comparison Gaussian #mu/#sigma_{RMS}", D.tag));
                        c2.SaveAs(Form("%s/DeltaEtaPhi_%s_MuSigmaVsE.png", baseDir.c_str(), D.tag));
                        c2.Close();
                    }
                }
            }   // ── end extra 4-way summary
        }       // ── end loop over folders / comparisons

        std::cout << "[makeCombinedDeltaEtaPhiPlots] FINISHED – output dir: "
                  << outDir << '\n';
    }

} // namespace (anonymous, verbose)

// ————————————————————————————————————————————————————————————————
//  φ residuals (global only) – now verbose
// ————————————————————————————————————————————————————————————————
inline void MakeDeltaPhiPlots(const std::vector<std::pair<double,double>>& eEdges,
                              EBinningMode                                binMode,
                              const std::array<std::vector<TH1F*>,5>&    hPhi,
                              const char*                                 outDir,
                              const std::vector<int>& variantsToPlot = {0,1,2,3,4})
{
    dbgBanner("MakeDeltaPhiPlots","Entering");

    std::cout << "  * energy slices  : " << eEdges.size() << '\n'
              << "  * variants       : ";
    for (int v : variantsToPlot) std::cout << v << ' ';
    std::cout << '\n';

    // Cache global φ information for the later η‑call
    gEdgePhi     = eEdges;
    gHPhiGlobal  = hPhi;
    gBinModePhi  = binMode;
    gVariantsPhi = variantsToPlot;

    const std::string dir = std::string(outDir)+"/deltaPhi";
    ensureDir(dir);

    makeResidualPlots(eEdges, binMode, hPhi, dir.c_str(), variantsToPlot);

    dbgBanner("MakeDeltaPhiPlots","Done");
}

// ─────────────────────────────────────────────────────────────────────────────
//  From‑Scratch vs Clusterizer comparison plots
//      (ex‑NEW ➌ block, lifted verbatim)
//      • keeps identical behaviour – only the wrapper name is new
// ─────────────────────────────────────────────────────────────────────────────
void makeScratchVsClusterizerCompare(
        const std::vector<std::pair<double,double>>&               eEdges,
        const std::array<std::vector<std::vector<TH1F*>>,5>&       hEtaVz,
        const char*                                                outDir)
{
    /* original NEW ➌ code starts here – 100 % unchanged */
    gStyle->SetOptTitle(0);

    struct CompSpec { const char* tag; int vA; int vB; };
    const std::array<CompSpec,4> comps = {
        CompSpec{ "Scratch_RAW_vs_Cluster_RAW" , 0, 2 },
        CompSpec{ "Scratch_CORR_vs_Cluster_CORR", 1, 3 },
        CompSpec{ "Scratch_RAW_vs_Scratch_CORR" , 0, 1 },
        CompSpec{ "Cluster_RAW_vs_Cluster_CORR", 2, 3 }
    };

    const std::string cmpRoot =
        std::string(outDir) + "/fromScratchVsClusterizerCompare";
    ensureDir(cmpRoot);

    constexpr int nCol = 4, nRow = 2, maxPads = nCol * nRow;

    const int N_Vz = N_VzBins;

    for (int iVz = 0; iVz < N_Vz; ++iVz)
    {
        const std::string vzDir = Form("%s/vz_%g_%g",
                                       cmpRoot.c_str(),
                                       vzEdge[iVz], vzEdge[iVz+1]);
        ensureDir(vzDir);

        for (const auto& C : comps)
        {
            /* ①  2×4 residual‑shape table – code identical to original */
            {
                TCanvas c(Form("cEtaCmp_%s_vz%d",C.tag,iVz), "", 1600,900);
                c.SetTopMargin(0.12);
                c.Divide(nCol, nRow, 0, 0);

                int pad = 1;
                for (std::size_t iE = 0;
                     iE < eEdges.size() && pad <= maxPads; ++iE)
                {
                    if (   iE >= hEtaVz[C.vA].size()
                        || iVz >= (int)hEtaVz[C.vA][iE].size()
                        || !hEtaVz[C.vA][iE][iVz] )
                        if (   iE >= hEtaVz[C.vB].size()
                            || iVz >= (int)hEtaVz[C.vB][iE].size()
                            || !hEtaVz[C.vB][iE][iVz] )
                            continue;

                    c.cd(pad++);
                    double yMax = 0.0;
                    auto* lg = new TLegend(0.15,0.8,0.55,0.97);
                    lg->SetBorderSize(0); lg->SetFillStyle(0);
                    lg->SetTextSize(0.045);

                    for (int which = 0; which < 2; ++which)
                    {
                        const int v = (which==0 ? C.vA : C.vB);
                        if (   iE >= hEtaVz[v].size()
                            || iVz >= (int)hEtaVz[v][iE].size()
                            || !hEtaVz[v][iE][iVz]
                            || hEtaVz[v][iE][iVz]->Integral()==0 )
                            continue;

                        TH1F* src = hEtaVz[v][iE][iVz];
                        auto* h = cloneAndNormPdf(src, v);

                        Style_t mk = 20;  Color_t col = kBlack;
                        switch (v) {
                            case 0:  mk = 25; col = kRed;     break;
                            case 1:  mk = 21; col = kRed;     break;
                            case 2:  mk = 24; col = kBlue+1;  break;
                            case 3:  mk = 20; col = kBlue+1;  break;
                        }
                        h->SetMarkerColor(col); h->SetLineColor(col);
                        h->SetMarkerStyle(mk);  h->SetMarkerSize(1.1);

                        yMax = std::max(yMax, h->GetMaximum()*1.15);

                        if (which==0)
                        {
                            h->SetTitle("");
                            h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                            h->GetYaxis()->SetTitle("Probability density");
                            h->Draw("E");
                        } else h->Draw("E SAME");

                        lg->AddEntry(h, kLab[v], "lep");
                    }

                    for (auto* obj : *gPad->GetListOfPrimitives())
                        if (auto* hh = dynamic_cast<TH1*>(obj))
                            hh->GetYaxis()->SetRangeUser(0., yMax);

                    TLatex vzlab; vzlab.SetNDC();
                    vzlab.SetTextSize(0.037); vzlab.SetTextAlign(33);
                    vzlab.DrawLatex(0.88,0.92,
                        Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]));
                    lg->Draw();
                }
                c.SaveAs(Form("%s/EtaCmp_%s_Table.png", vzDir.c_str(), C.tag));
                c.Close();
            }

            /* ②  μ(E) / σ(E) summary – identical logic to original */
            {
                std::vector<double> eCtr;  eCtr.reserve(eEdges.size());
                std::array<std::vector<double>,2> MU,dMU,SG,dSG;

                for (std::size_t iE=0;iE<eEdges.size();++iE)
                {
                    const double eMid=0.5*(eEdges[iE].first+eEdges[iE].second);
                    eCtr.push_back(eMid);
                    for(int k=0;k<2;++k){
                        MU[k].push_back(0.0); dMU[k].push_back(0.0);
                        SG[k].push_back(0.0); dSG[k].push_back(0.0);
                    }

                    for(int k=0;k<2;++k)
                    {
                        const int v = (k==0?C.vA:C.vB);
                        if (iE>=hEtaVz[v].size() ||
                            iVz>=(int)hEtaVz[v][iE].size() ||
                            !hEtaVz[v][iE][iVz] ||
                            hEtaVz[v][iE][iVz]->Integral()==0) continue;

                        const FitRes fr = robustGaussianFit(hEtaVz[v][iE][iVz]);
                        MU[k].back() = fr.mu;  dMU[k].back() = fr.dmu;
                        SG[k].back() = fr.sg;  dSG[k].back() = fr.dsg;
                    }
                }

                const double xMax = eCtr.back() +
                     0.5*(eEdges.back().second-eEdges.back().first);
                std::vector<double> ex(eCtr.size(),0.0);

                TCanvas c(Form("cEtaCmpMS_%s_vz%d",C.tag,iVz), "", 900,800);
                c.Divide(1,2,0,0);

                /* μ‑pad */
                c.cd(1);
                gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.06);

                double muLo=1e30, muHi=-1e30;
                for(int k=0;k<2;++k)
                    for (std::size_t j=0;j<MU[k].size();++j) {
                        muLo=std::min(muLo,MU[k][j]-dMU[k][j]);
                        muHi=std::max(muHi,MU[k][j]+dMU[k][j]);
                    }
                TH1F frU("frU",";E [GeV];#mu [rad]",1,0.0,xMax);
                frU.SetMinimum(std::min(0.0, muLo-0.1*(muHi-muLo)));
                frU.SetMaximum(muHi+0.1*(muHi-muLo));
                frU.Draw("AXIS");
                TLine ref0(0,0,xMax,0); ref0.SetLineStyle(2); ref0.Draw("same");

                TLegend legU(0.70,0.75,0.92,0.92);
                legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.04);

                for(int k=0;k<2;++k){
                    const int v = (k==0?C.vA:C.vB);
                    Style_t mk=20; Color_t col=kBlack;
                    switch(v){
                        case 0: mk=25; col=kRed;     break;
                        case 1: mk=21; col=kRed;     break;
                        case 2: mk=24; col=kBlue+1;  break;
                        case 3: mk=20; col=kBlue+1;  break;
                    }

                    TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                       eCtr.data(), MU[k].data(),
                                                       ex.data(),  dMU[k].data());
                    g->SetMarkerColor(col); g->SetLineColor(col);
                    g->SetMarkerStyle(mk);  g->SetMarkerSize(1.1);
                    g->Draw("P SAME");

                    auto* m = new TMarker(0,0,mk);
                    m->SetMarkerColor(col); m->SetMarkerSize(1.1);
                    legU.AddEntry(m,kLab[v],"p");
                }
                legU.Draw();

                /* σ‑pad */
                c.cd(2);
                gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);

                double sgHi=-1e30;
                for(int k=0;k<2;++k)
                    for(double s:SG[k]) sgHi=std::max(sgHi,s);
                TH1F frL("frL",";E [GeV];#sigma_{RMS} [rad]",1,0.0,xMax);
                frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");

                for(int k=0;k<2;++k){
                    const int v=(k==0?C.vA:C.vB);
                    Style_t mk=20; Color_t col=kBlack;
                    switch(v){
                        case 0: mk=25; col=kRed;     break;
                        case 1: mk=21; col=kRed;     break;
                        case 2: mk=24; col=kBlue+1;  break;
                        case 3: mk=20; col=kBlue+1;  break;
                    }
                    TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                       eCtr.data(), SG[k].data(),
                                                       ex.data(),  dSG[k].data());
                    g->SetMarkerColor(col); g->SetLineColor(col);
                    g->SetMarkerStyle(mk);  g->SetMarkerSize(1.1);
                    g->Draw("P SAME");
                }

                c.SaveAs(Form("%s/EtaCmp_%s_MuSigma.png", vzDir.c_str(), C.tag));
                c.Close();
            }
        }   /* end CompSpec */
    }       /* end iVz loop */
}           /* end helper */


void makeThreeSliceOverlays(
        int                                                     v,          // variant index
        const std::string&                                      vName,      // sanitised variant label
        const std::string&                                      vDir,       // base dir for this variant
        const std::vector<std::pair<double,double>>&            eEdges,
        const std::array<std::vector<std::vector<TH1F*>>,5>&    hEtaVz,
        const char*                                              deltaAxis,  // "#Delta#eta" or "#Delta#phi"
        char                                                     signTag = '\0')   // 'P' (+z), 'N' (-z), or '\0' (|z|)
{
    const int N_Vz = N_VzBins;                  // use global count
    const std::array<int,3> pick3 = {0, N_Vz/2, N_Vz-1};

    // pick a consistent basename for output files
    const std::string baseName =
        (std::string(deltaAxis).find("#phi") != std::string::npos) ? "DeltaPhi" : "DeltaEta";

    // palette-based color for any iVz
    auto vzColor = [&](int iVz)->Color_t {
        const int ncol = gStyle->GetNumberOfColors();
        const int idx  = std::max(0, std::min(ncol-1,
                        (int)std::floor( (double)iVz / std::max(1, N_Vz-1) * (ncol-1) )));
        return TColor::GetColorPalette(idx);
    };

    std::string vDir3 = vDir + "/threeSlices";
    ensureDir(vDir3);

    // ───────────────────── 1 × 3 canvas helper (Elo / Ehi) ──────────────────
    auto drawThreeTable =
    [&](const char* tag, std::size_t startIdx)
    {
        TCanvas cT3(Form("cVz3_%s_%d",tag,v),"vz_three",3200,1200);
        cT3.SetTopMargin(0.07); cT3.SetBottomMargin(0.06);
        cT3.Divide(3,1,0.002,0.0);

        for (int pad=1; pad<=3; ++pad)
        {
            const std::size_t iE = startIdx + pad-1;
            if (iE >= eEdges.size()) continue;

            bool have=false;
            for (int iVz: pick3)
                if (iVz < (int)hEtaVz[v][iE].size() &&
                    hEtaVz[v][iE][iVz] &&
                    hEtaVz[v][iE][iVz]->Integral()>0){ have=true; break; }
            if (!have) continue;

            cT3.cd(pad);
            gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.08);  gPad->SetBottomMargin(0.12);

            double yMax = 0.0;
            auto* leg = new TLegend(0.13,0.73,0.50,0.93);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
            leg->AddEntry((TObject*)nullptr,
                          Form("E: [%.0f, %.0f) GeV",
                               eEdges[iE].first, eEdges[iE].second), "");

            bool axesDrawn=false;
            for (int slot=0; slot<3; ++slot)
            {
                const int iVz = pick3[slot];
                if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                TH1F* src = hEtaVz[v][iE][iVz];
                if (!src || src->Integral()==0) continue;

                auto* h = cloneAndNormPdf(src, v);
                const Color_t col = vzColor(iVz);
                h->SetMarkerColor(col); h->SetLineColor(col);
                h->SetMarkerStyle(20);  h->SetMarkerSize(1.2);
                yMax = std::max(yMax, h->GetMaximum()*1.20);

                if (!axesDrawn) {
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle(Form("%s  [rad]", deltaAxis));
                    h->GetYaxis()->SetTitle("Probability density");
                    h->GetXaxis()->SetTitleSize(0.045);
                    h->GetXaxis()->SetLabelSize(0.035);
                    h->GetXaxis()->SetTitleOffset(1.10);
                    h->GetYaxis()->SetTitleSize(0.045);
                    h->GetYaxis()->SetLabelSize(0.035);
                    h->GetYaxis()->SetTitleOffset(1.40);
                    h->Draw("E");
                    axesDrawn = true;
                } else {
                    h->Draw("E SAME");
                }

                // legend label with optional ± sign
                if (signTag == 'P')
                    leg->AddEntry(h,
                        Form("z: +[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
                else if (signTag == 'N')
                    leg->AddEntry(h,
                        Form("z: -[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
                else
                    leg->AddEntry(h,
                        Form("|z_{vtx}|: [%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
            }

            for (auto* obj : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(obj))
                    h->GetYaxis()->SetRangeUser(0., yMax);

            leg->Draw();
            gPad->Modified(); gPad->Update();
        }

        cT3.cd(0);
        TLatex hd; hd.SetNDC(); hd.SetTextAlign(22); hd.SetTextSize(0.05);
        hd.DrawLatex(0.5,0.96,
            (std::string(kLab[v])+"  , low / mid / high slices  ("+tag+')').c_str());

        cT3.Modified(); cT3.Update();
        cT3.SaveAs(Form("%s/%s_Vtx3Slices_%s_%s.png",
                        vDir3.c_str(), baseName.c_str(), vName.c_str(), tag));
        cT3.Close();
    };

    if (eEdges.size() >= 1) drawThreeTable("Elo", 0);
    if (eEdges.size() >= 4) drawThreeTable("Ehi", eEdges.size() >= 3 ? eEdges.size()-3 : 0);

    // μ(E) & σ(E) summary for the same three slices
    std::vector<double> eCtr3; eCtr3.reserve(eEdges.size());
    std::array<std::vector<double>,3> MU3,dMU3,SG3,dSG3;

    for (std::size_t iE=0;iE<eEdges.size();++iE){
        const double eMid=0.5*(eEdges[iE].first+eEdges[iE].second);
        eCtr3.push_back(eMid);
        for(int slot=0;slot<3;++slot){
            MU3[slot].push_back(0.0); dMU3[slot].push_back(0.0);
            SG3[slot].push_back(0.0); dSG3[slot].push_back(0.0);
            int iVz=pick3[slot];
            if(iE>=hEtaVz[v].size()||iVz>=(int)hEtaVz[v][iE].size()) continue;
            TH1F* src=hEtaVz[v][iE][iVz];
            if(!src||src->Integral()==0) continue;
            FitRes fr=robustGaussianFit(src);
            MU3[slot].back()=fr.mu;  dMU3[slot].back()=fr.dmu;
            SG3[slot].back()=fr.sg;  dSG3[slot].back()=fr.dsg;
        }
    }

    const double xMax3 = eCtr3.back() +
        0.5*(eEdges.back().second-eEdges.back().first);
    std::vector<double> ex3(eCtr3.size(),0.0);

    TCanvas cMS3(Form("cVtx3MuSig_%d",v),"vtx3_mu_sigma",900,800);
    cMS3.Divide(1,2,0,0);

    // μ‑pad
    cMS3.cd(1);
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.04);
    double muLo=1e30,muHi=-1e30;
    for(int slot=0;slot<3;++slot)
        for(std::size_t k=0;k<MU3[slot].size();++k){
            muLo=std::min(muLo,MU3[slot][k]-dMU3[slot][k]);
            muHi=std::max(muHi,MU3[slot][k]+dMU3[slot][k]);
        }
    TH1F frU3("frU3",Form(";E  [GeV];%s  [rad]",deltaAxis),1,0.0,xMax3);
    frU3.SetMinimum(muLo-0.1*(muHi-muLo));
    frU3.SetMaximum(muHi+0.1*(muHi-muLo));
    frU3.Draw("AXIS");
    TLine ref0(0,0,xMax3,0); ref0.SetLineStyle(2); ref0.Draw("same");

    TLegend lu3(0.70,0.60,0.92,0.90);
    lu3.SetBorderSize(0); lu3.SetFillStyle(0); lu3.SetTextSize(0.035);

    for(int slot=0;slot<3;++slot){
        int iVz=pick3[slot];
        TGraphErrors* g=new TGraphErrors(eCtr3.size(), eCtr3.data(),
                                         MU3[slot].data(),
                                         ex3.data(),  dMU3[slot].data());
        const Color_t col = vzColor(iVz);
        g->SetMarkerColor(col); g->SetLineColor(col);
        g->SetMarkerStyle(20);  g->SetMarkerSize(1.1);
        g->Draw("P SAME");
        if (signTag == 'P')
            lu3.AddEntry(g,Form("+[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
        else if (signTag == 'N')
            lu3.AddEntry(g,Form("-[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
        else
            lu3.AddEntry(g,Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
    }
    lu3.Draw();

    // σ‑pad
    cMS3.cd(2);
    gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);
    double sgHi=-1e30;
    for(int slot=0;slot<3;++slot)
        for(double vsg:SG3[slot]) sgHi=std::max(sgHi,vsg);
    TH1F frL3("frL3",";E  [GeV];#sigma_{RMS}  [rad]",1,0.0,xMax3);
    frL3.SetMinimum(0.0); frL3.SetMaximum(1.15*sgHi); frL3.Draw("AXIS");

    for(int slot=0;slot<3;++slot){
        int iVz=pick3[slot];
        TGraphErrors* g=new TGraphErrors(eCtr3.size(), eCtr3.data(),
                                         SG3[slot].data(),
                                         ex3.data(),  dSG3[slot].data());
        const Color_t col = vzColor(iVz);
        g->SetMarkerColor(col); g->SetLineColor(col);
        g->SetMarkerStyle(20);  g->SetMarkerSize(1.1);
        g->Draw("P SAME");
    }

    cMS3.cd(0);
    TLatex ttl; ttl.SetNDC(); ttl.SetTextAlign(22); ttl.SetTextSize(0.04);
    ttl.DrawLatex(0.5,0.97,
        (std::string(kLab[v])+"  #mu/#sigma_{RMS} vs E  (low/mid/high slices)").c_str());
    cMS3.SaveAs(Form("%s/%s_Vtx3_MuSigmaVsE_%s.png",
                     vDir3.c_str(), baseName.c_str(), vName.c_str()));
    cMS3.Close();
}


void makeGlobalMatrixSummary(
        const std::vector<std::pair<double,double>>&            eEdges,
        const std::array<std::vector<std::vector<TH1F*>>,5>&    hEtaVz,
        const char*                                             outDir,
        const std::vector<int>&                                 variantsToPlot,
        const char*                                             deltaAxis)   // "#Delta#eta" or "#Delta#phi"
{
    const int N_Vz = N_VzBins;

    const std::string mDir = std::string(outDir) + "/vertexCondensedSummary";
    ensureDir(mDir);

    const std::string baseName =
        (std::string(deltaAxis).find("#phi") != std::string::npos) ? "DeltaPhi" : "DeltaEta";

    auto vzColor = [&](int iVz)->Color_t {
        const int ncol = gStyle->GetNumberOfColors();
        const int idx  = std::max(0, std::min(ncol-1,
                        (int)std::floor( (double)iVz / std::max(1, N_Vz-1) * (ncol-1) )));
        return TColor::GetColorPalette(idx);
    };

    const std::size_t nEPlot = std::min<std::size_t>(3, eEdges.size());
    const int nColM = static_cast<int>(nEPlot);
    const int nRowM = static_cast<int>(variantsToPlot.size());

    TCanvas cMat("cVzMatrix", "vtx_matrix", 4500, 3000);
    cMat.SetTopMargin(0.03);  cMat.SetBottomMargin(0.03);
    cMat.SetLeftMargin(0.02); cMat.SetRightMargin (0.02);
    cMat.Divide(nColM, nRowM, 0.008, 0.008);

    int padMat = 1;
    for (int r = 0; r < nRowM; ++r)
    {
        const int v = variantsToPlot[r];

        for (std::size_t c = 0; c < nEPlot; ++c, ++padMat)
        {
            cMat.cd(padMat);

            if (c >= hEtaVz[v].size() ||
                hEtaVz[v][c].empty() ||
                (!hEtaVz[v][c][0] && !hEtaVz[v][c].back()))
                continue;

            double ymax = 0.0;
            std::array<FitRes,2> fr{};

            auto* leg = new TLegend(0.70, 0.22, 0.88, 0.42);
            leg->SetBorderSize(0);  leg->SetFillStyle(0);  leg->SetTextSize(0.055);

            for (int k = 0; k < 2; ++k)
            {
                const int iVz = (k == 0 ? 0 : N_Vz - 1);   // low / high
                if (iVz >= (int)hEtaVz[v][c].size()) continue;
                TH1F* src = hEtaVz[v][c][iVz];
                if (!src || src->Integral() == 0)   continue;

                auto* h = cloneAndNormPdf(src, v);
                const Color_t col = vzColor(iVz);
                h->SetMarkerColor(col);
                h->SetLineColor(col);
                h->SetMarkerStyle(k == 0 ? 20 : 24);
                h->SetMarkerSize(1.2);
                ymax = std::max(ymax, h->GetMaximum() * 1.15);

                if (k == 0) {
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle(Form("%s  [rad]", deltaAxis));
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                } else
                    h->Draw("E SAME");

                fr[k] = robustGaussianFit(src);
                leg->AddEntry(h,
                              Form("[%.0f, %.0f) cm",
                                   vzEdge[iVz], vzEdge[iVz + 1]), "p");
            }

            gPad->Update();
            for (auto* o : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(o))
                    h->GetYaxis()->SetRangeUser(0., ymax);

            leg->Draw();

            TLatex vLab; vLab.SetNDC();
            vLab.SetTextSize(0.055); vLab.SetTextAlign(11);
            vLab.DrawLatex(gPad->GetLeftMargin() + 0.03, 0.86, kLab[v]);

            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.055); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.91, 0.88,
                Form("[%.0f, %.0f) GeV",
                     eEdges[c].first, eEdges[c].second));

            TLatex tbl; tbl.SetNDC(); tbl.SetTextFont(132); tbl.SetTextAlign(13);
            tbl.SetTextSize(0.052);
            double ytbl = 0.33;
            for (int k = 0; k < 2; ++k)
            {
                int iVz = (k == 0 ? 0 : N_Vz - 1);
                tbl.DrawLatex(gPad->GetLeftMargin() + 0.04, ytbl,
                    Form("#mu=%.2g  #sigma_{RMS}=%.2g  z:%g-%g",
                         fr[k].mu, fr[k].sg,
                         vzEdge[iVz], vzEdge[iVz + 1]));
                ytbl -= 0.06;
            }
        }
    }

    cMat.SaveAs(Form("%s/%s_Vtx_Matrix.png", mDir.c_str(), baseName.c_str()));
    cMat.Close();
}


void makeEtaVertexTables(const std::vector<std::pair<double,double>>&               eEdges,
                         const std::array<std::vector<std::vector<TH1F*>>,5>&       hEtaVz,
                         const char*                                                outDir,             // “…/deltaEta|deltaPhi/vertexDependent”
                         const std::vector<int>&                                    variantsToPlot = {0,1,2,3,4},
                         const char*                                                deltaAxis = "#Delta#eta",   // or "#Delta#phi"
                         char                                                       signTag = '\0')            // 'P' | 'N' | '\0'
{
    constexpr int nCol = 4, nRow = 2, maxPads = nCol * nRow;
    const int N_Vz = N_VzBins;

    const std::string baseName =
        (std::string(deltaAxis).find("#phi") != std::string::npos) ? "DeltaPhi" : "DeltaEta";

    auto vzColor = [&](int iVz)->Color_t {
        const int ncol = gStyle->GetNumberOfColors();
        const int idx  = std::max(0, std::min(ncol-1,
                        (int)std::floor( (double)iVz / std::max(1, N_Vz-1) * (ncol-1) )));
        return TColor::GetColorPalette(idx);
    };

    // iterate over reconstruction variants
    for (int v : variantsToPlot)
    {
        // sanitise variant name for file‑system usage
        std::string vName = kLab[v];
        for (char& c : vName)
            if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';

        const std::string vDir = std::string(outDir) + "/" + vName;
        ensureDir(vDir);

        TCanvas cTab(Form("cVzTab_%d", v), kLab[v], 1600, 900);
        cTab.SetTopMargin(0.12);
        cTab.Divide(nCol, nRow, 0, 0);

        int pad = 1;
        for (std::size_t iE = 0; iE < eEdges.size() && pad <= maxPads; ++iE)
        {
            // skip if no vertex histograms exist for this energy bin
            bool haveSlice = false;
            for (int iVz = 0; iVz < N_Vz && !haveSlice; ++iVz)
                if (iVz < (int)hEtaVz[v][iE].size() && hEtaVz[v][iE][iVz])
                    haveSlice = true;
            if (!haveSlice) continue;

            int curPad = pad++;
            cTab.cd(curPad);

            // containers sized to N_Vz
            std::vector<FitRes> fitInfo(N_Vz);
            std::vector<bool>   sliceUsed(N_Vz,false);

            double yMax = 0.0;

            const double xLeg1 = gPad->GetLeftMargin() + 0.05;
            const double xLeg2 = xLeg1 + 0.26;
            TLegend* leg = new TLegend(xLeg1, 0.73, xLeg2, 0.93);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);

            // overlay all vertex slices for this energy bin
            for (int iVz = 0; iVz < N_Vz; ++iVz)
            {
                if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                TH1F* src = hEtaVz[v][iE][iVz];
                if (!src || src->Integral() == 0) continue;

                auto* h = cloneAndNormPdf(src, v);
                const Color_t col = vzColor(iVz);
                h->SetMarkerColor(col);
                h->SetLineColor  (col);
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.9);

                yMax = std::max(yMax, h->GetMaximum() * 1.15);

                if (leg->GetNRows() == 0)
                {
                    h->SetTitle(Form("%s  [%.0f, %.0f) GeV",
                                     deltaAxis, eEdges[iE].first, eEdges[iE].second));
                    h->GetXaxis()->SetTitle(Form("%s  [rad]", deltaAxis));
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                }
                else
                    h->Draw("E SAME");

                // robust Gaussian fit
                const FitRes fr = robustGaussianFit(src);
                fitInfo[iVz]   = fr;
                sliceUsed[iVz] = true;

                // legend label with optional ± sign
                if (signTag == 'P')
                    leg->AddEntry(h,
                        Form("z: +[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                        "lep");
                else if (signTag == 'N')
                    leg->AddEntry(h,
                        Form("z: -[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                        "lep");
                else
                    leg->AddEntry(h,
                        Form("|z_{vtx}|: [%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                        "lep");
            }

            gPad->Update();
            for (auto* obj : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(obj))
                    h->GetYaxis()->SetRangeUser(0., yMax);

            leg->Draw();

            // energy‑range label (top‑right)
            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.037); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.88, 0.92,
                Form("[%.0f, %.0f) GeV", eEdges[iE].first, eEdges[iE].second));

            // μ | σ | z  table
            const int rowIdx = (curPad-1) / nCol;
            double yRow = (rowIdx == 0 ? 0.40 : 0.50);

            const double xTbl1 = gPad->GetLeftMargin() + 0.04;
            const double xTbl2 = xTbl1 + 0.12;
            const double xTbl3 = xTbl2 + 0.14;

            TLatex tbl; tbl.SetNDC(); tbl.SetTextAlign(13);
            tbl.SetTextSize(0.032);
            tbl.SetTextFont(132);
            tbl.DrawLatex(xTbl1, yRow, "#mu");
            tbl.DrawLatex(xTbl2, yRow, "#sigma_{RMS}");
            tbl.DrawLatex(xTbl3, yRow, "z\\,[cm]");
            tbl.SetTextSize(0.026);
            yRow -= 0.050;

            for (int iVz = 0; iVz < N_Vz; ++iVz)
            {
                if (!sliceUsed[iVz]) continue;
                tbl.SetTextColor(vzColor(iVz));
                tbl.DrawLatex(xTbl1, yRow, Form("%.2g", fitInfo[iVz].mu));
                tbl.DrawLatex(xTbl2, yRow, Form("%.2g", fitInfo[iVz].sg));
                if (signTag == 'P')
                    tbl.DrawLatex(xTbl3, yRow, Form("+%.0f..+%.0f", vzEdge[iVz], vzEdge[iVz+1]));
                else if (signTag == 'N')
                    tbl.DrawLatex(xTbl3, yRow, Form("-%.0f..-%.0f", vzEdge[iVz], vzEdge[iVz+1]));
                else
                    tbl.DrawLatex(xTbl3, yRow, Form("%.0f-%.0f", vzEdge[iVz], vzEdge[iVz+1]));
                yRow -= 0.045;
            }
            tbl.SetTextColor(kBlack);

            gPad->Modified();
            gPad->Update();
            gPad->SaveAs(Form("%s/%s_VtxSlices_%s_E%02zu.png",
                               vDir.c_str(), baseName.c_str(), vName.c_str(), iE));
        }

        // canvas-level title
        cTab.cd(0);
        TLatex head; head.SetNDC(); head.SetTextAlign(22);
        head.SetTextSize(0.045);
        head.DrawLatex(0.5, 0.975, kLab[v]);

        // full 2×4 table
        cTab.SaveAs(Form("%s/%s_VtxSlices_%s.png",
                         vDir.c_str(), baseName.c_str(), vName.c_str()));
        cTab.Close();

        // μ(E) & σ(E) summary (six curves → now N_Vz curves)
        {
            std::vector<double> eCtr;  eCtr.reserve(eEdges.size());
            std::vector<std::vector<double>> MUvz(N_Vz), dMUvz(N_Vz), SGvz(N_Vz), dSGvz(N_Vz);

            for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
            {
                const double eMid = 0.5 * (eEdges[iE].first + eEdges[iE].second);
                eCtr.push_back(eMid);

                if (iE >= hEtaVz[v].size()) continue;
                for (int iVz = 0; iVz < N_Vz; ++iVz)
                {
                    MUvz[iVz].push_back(0.0);  dMUvz[iVz].push_back(0.0);
                    SGvz[iVz].push_back(0.0);  dSGvz[iVz].push_back(0.0);

                    if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                    TH1F* src = hEtaVz[v][iE][iVz];
                    if (!src || src->Integral() == 0)     continue;

                    const FitRes fr = robustGaussianFit(src);
                    MUvz[iVz].back()  = fr.mu;  dMUvz[iVz].back() = fr.dmu;
                    SGvz[iVz].back()  = fr.sg;  dSGvz[iVz].back() = fr.dsg;
                }
            }

            TCanvas cMS(Form("cVtxMuSig_%d", v), "vtx_mu_sigma", 900, 800);
            cMS.Divide(1, 2, 0, 0);
            std::vector<double> ex(eCtr.size(), 0.0);

            // μ pad
            cMS.cd(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.06);
            gPad->SetBottomMargin(0.04);

            double muMin =  1e30, muMax = -1e30;
            for (int iVz = 0; iVz < N_Vz; ++iVz)
                for (std::size_t j = 0; j < MUvz[iVz].size(); ++j) {
                    muMin = std::min(muMin, MUvz[iVz][j] - dMUvz[iVz][j]);
                    muMax = std::max(muMax, MUvz[iVz][j] + dMUvz[iVz][j]);
                }
            const double xMax = eCtr.back() +
                                0.5 * (eEdges.back().second - eEdges.back().first);

            TH1F frU("frU",Form(";E  [GeV];%s  [rad]",deltaAxis),1,0.0,xMax);
            frU.SetMinimum(muMin - 0.05*(muMax-muMin));
            frU.SetMaximum(muMax + 0.05*(muMax-muMin));
            frU.Draw("AXIS");
            TLine ref0(0.0, 0.0, xMax, 0.0); ref0.SetLineStyle(2); ref0.Draw("same");

            for (int iVz = 0; iVz < N_Vz; ++iVz) {
                auto* g = new TGraphErrors((int)eCtr.size(),
                                           eCtr.data(), MUvz[iVz].data(),
                                           ex.data(),  dMUvz[iVz].data());
                const Color_t col = vzColor(iVz);
                g->SetMarkerColor(col);
                g->SetLineColor  (col);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.1);
                g->Draw("P SAME");
            }

            // σ pad
            cMS.cd(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.06);
            gPad->SetTopMargin(0.06);
            gPad->SetBottomMargin(0.35);

            double sgMax = -1e30;
            for (int iVz = 0; iVz < N_Vz; ++iVz)
                for (double vsg : SGvz[iVz]) sgMax = std::max(sgMax, vsg);

            TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,0.0,xMax);
            frL.SetMinimum(0.0);
            frL.SetMaximum(sgMax * 1.05);
            frL.Draw("AXIS");

            for (int iVz = 0; iVz < N_Vz; ++iVz) {
                auto* g = new TGraphErrors((int)eCtr.size(),
                                           eCtr.data(), SGvz[iVz].data(),
                                           ex.data(),  dSGvz[iVz].data());
                const Color_t col = vzColor(iVz);
                g->SetMarkerColor(col);
                g->SetLineColor  (col);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.1);
                g->Draw("P SAME");
            }

            cMS.SaveAs(Form("%s/%s_Vtx_MuSigmaVsE_%s.png",
                            vDir.c_str(), baseName.c_str(), vName.c_str()));
            cMS.Close();
        }

        // three-slice overlays (low/mid/high)
        makeThreeSliceOverlays(v, vName, vDir, eEdges, hEtaVz, deltaAxis, signTag);
    }

    // “condensed” summary – one canvas per energy slice (low/high only)
    const std::string summaryDir = std::string(outDir) + "/vertexCondensedSummary";
    ensureDir(summaryDir);

    constexpr int nColC   = 2, nRowC = 3, maxPadsC = nColC * nRowC;
    const std::array<int,2> pickVz = {0, N_Vz-1};

    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        TCanvas cSum(Form("cCond_%zu", iE), "condensed_vtx", 2400, 2400);
        cSum.SetTopMargin   (0.06);
        cSum.SetBottomMargin(0.04);
        cSum.SetLeftMargin  (0.02);
        cSum.SetRightMargin (0.02);
        cSum.Divide(nColC, nRowC, 0.005, 0.005);

        for (int ip = 1; ip <= maxPadsC; ++ip) {
            cSum.cd(ip);
            gPad->SetLeftMargin  (0.12);
            gPad->SetRightMargin (0.04);
            gPad->SetTopMargin   (0.08);
            gPad->SetBottomMargin(0.12);
        }

        int padC = 1;
        for (int v : variantsToPlot)
        {
            if (padC > maxPadsC) break;
            if (iE >= hEtaVz[v].size()) continue;

            cSum.cd(padC++);

            double ymax = 0.0;
            std::array<FitRes,2> fr{};
            auto* l = new TLegend(0.65, 0.25, 0.88, 0.4);
            l->SetBorderSize(0);  l->SetFillStyle(0);  l->SetTextSize(0.042);

            for (int k = 0; k < 2; ++k)
            {
                int iVz = pickVz[k];
                if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                TH1F* src = hEtaVz[v][iE][iVz];
                if (!src || src->Integral()==0)    continue;

                auto* h = cloneAndNormPdf(src, v);
                const Color_t col = vzColor(iVz);
                h->SetMarkerColor(col);
                h->SetLineColor  (col);
                h->SetMarkerStyle( k==0 ? 20 : 24 );
                h->SetMarkerSize(1.0);

                ymax = std::max(ymax, h->GetMaximum()*1.15);

                if (k==0) {
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle(Form("%s  [rad]", deltaAxis));
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                } else
                    h->Draw("E SAME");

                fr[k] = robustGaussianFit(src);

                if (signTag == 'P')
                    l->AddEntry(h, Form("+[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]), "p");
                else if (signTag == 'N')
                    l->AddEntry(h, Form("-[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]), "p");
                else
                    l->AddEntry(h, Form("[%.0f, %.0f) cm",  vzEdge[iVz], vzEdge[iVz+1]), "p");
            }

            gPad->Update();
            for (auto* o : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(o))
                    h->GetYaxis()->SetRangeUser(0., ymax);

            l->Draw();

            TLatex vtxt; vtxt.SetNDC();
            vtxt.SetTextSize(0.045); vtxt.SetTextAlign(11);
            vtxt.DrawLatex(gPad->GetLeftMargin()+0.04, 0.85, kLab[v]);

            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.045); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.88, 0.87,
                Form("[%.0f, %.0f) GeV", eEdges[iE].first, eEdges[iE].second));

            TLatex t; t.SetNDC(); t.SetTextFont(132); t.SetTextAlign(13);
            t.SetTextSize(0.035);

            const double x0 = gPad->GetLeftMargin() + 0.04;
            double y0 = 0.30;
            for (int k = 0; k < 2; ++k)
            {
                int iVz = pickVz[k];
                t.SetTextColor(vzColor(iVz));
                t.DrawLatex(x0, y0,
                    Form("#mu=%.2g  #sigma_{RMS}=%.2g  z:%g-%g",
                         fr[k].mu, fr[k].sg,
                         vzEdge[iVz], vzEdge[iVz+1]));
                y0 -= 0.05;
            }
            t.SetTextColor(kBlack);
        }

        cSum.SaveAs(Form("%s/%s_Condensed_E%02zu.png",
                         summaryDir.c_str(), baseName.c_str(), iE));
        cSum.Close();
    }

    // global matrix summary (first 3 energy bins)
    makeGlobalMatrixSummary(eEdges, hEtaVz, outDir, variantsToPlot, deltaAxis);
}



// ————————————————————————————————————————————————————————————————
//  η residuals – global + |z_vtx| resolved, plus φ + η overlay – verbose
// ————————————————————————————————————————————————————————————————
inline void MakeDeltaEtaPlots(
        const std::vector<std::pair<double,double>>&               eEdges,
        EBinningMode                                               binMode,
        const std::array<std::vector<TH1F*>,5>&                    hEtaGlobal,
        const std::array<std::vector<std::vector<TH1F*>>,5>&       hEtaVz,
        const char*                                                outDir,
        const std::vector<int>& variantsToPlot = {0,1,2,3,4})
{
    dbgBanner("MakeDeltaEtaPlots","Entering");

    //-------------------------------------------------------------------
    // 1) global Δη
    //-------------------------------------------------------------------
    const std::string etaDir = std::string(outDir)+"/deltaEta";
    ensureDir(etaDir);

    std::cout << "[Info] Calling makeResidualPlots for global Δη\n";
    makeResidualPlots(eEdges, binMode, hEtaGlobal,
                      etaDir.c_str(), variantsToPlot);

    //-------------------------------------------------------------------
    // 2) φ + η combined
    //-------------------------------------------------------------------
    if (!gHPhiGlobal.empty())
    {
        std::cout << "[Info] gHPhiGlobal is populated – building combined plots.\n";

        // quick integrity check
        if (gEdgePhi.size() != eEdges.size())
        {
            std::cerr << "[WARN] Number of φ energy slices ("
                      << gEdgePhi.size()
                      << ") ≠ η energy slices (" << eEdges.size()
                      << ").  Proceeding with min size.\n";
        }

        const std::string comboDir =
            std::string(outDir)+"/combinedDeltaEtaPhi";

        makeCombinedDeltaEtaPhiPlots(eEdges, binMode,
                                     hEtaGlobal, gHPhiGlobal,
                                     comboDir.c_str(), variantsToPlot);
    }
    else
    {
        std::cerr << "[WARN] gHPhiGlobal is empty – combined φ+η plots skipped.\n";
    }

    //-------------------------------------------------------------------
    // 3) vertex‑dependent Δη slices
    //-------------------------------------------------------------------
    static constexpr std::array<float,7> vzEdge = {0,5,10,15,20,25,30};
    constexpr int N_Vz = vzEdge.size() - 1;

    std::cout << "[Info] Producing vertex‑dependent Δη plots (N_Vz="
              << N_Vz << ")\n";

    for (int iVz = 0; iVz < N_Vz; ++iVz)
    {
        std::array<std::vector<TH1F*>,5> view{};   // [variant][iE]

        for (int v = 0; v < 5; ++v)
        {
            view[v].resize(eEdges.size(), nullptr);
            for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
            {
                if (iVz < (int)hEtaVz[v][iE].size())
                    view[v][iE] = hEtaVz[v][iE][iVz];
            }
        }

        const std::string tag =
            Form("vz%.0f_%.0f", vzEdge[iVz], vzEdge[iVz+1]);
        const std::string dir = etaDir + "/vertexDependent/" + tag;
        ensureDir(dir);

        std::cout << "  • Vz‑bin " << tag
                  << "  → calling makeResidualPlots\n";

        makeResidualPlots(eEdges, binMode, view,
                          dir.c_str(), variantsToPlot);
    }

    // ─── summary tables: Δη over all |z_vtx| slices, per variant ───
    makeEtaVertexTables(eEdges, hEtaVz,
                        (etaDir + "/vertexDependent").c_str(),
                        variantsToPlot);
    dbgBanner("MakeDeltaEtaPlots","Done");
}





/* ===================================================================== *
 *  PlotBcompare … overlay & interpolate b‑values, add Δb/b panel + stats
 * ===================================================================== */
void PlotBcompare(const std::map<double,double>& bRMS,
                  const std::map<double,double>& bPhi,
                  const char*                    outDir = "./Bcompare",
                  /* optional 1 σ errors; pass {} to skip */
                  const std::map<double,double>& eRMS  = {},
                  const std::map<double,double>& ePhi  = {})
{
  /* 0) General set‑up ------------------------------------------------- */
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gSystem->mkdir(outDir,true);

  if (bRMS.empty() || bPhi.empty())
  { std::cerr << "[PlotBcompare] ERROR: empty input map – abort\n"; return; }

  /* 1) Harvest the common energy points ------------------------------ */
  std::vector<double> eCtr, vRMS, vPhi, sRMS, sPhi;
  for (const auto& kv : bRMS)
  {
    const auto itP = bPhi.find(kv.first);
    if (itP == bPhi.end()) continue;
    if (!std::isfinite(kv.second) || !std::isfinite(itP->second)) continue;

    eCtr .push_back(kv.first);
    vRMS .push_back(kv.second);
    vPhi .push_back(itP->second);
    sRMS .push_back( (eRMS.count(kv.first)? eRMS.at(kv.first) : 0.0) );
    sPhi .push_back( (ePhi.count(kv.first)? ePhi.at(kv.first) : 0.0) );
  }
  if (eCtr.size()<2)
  { std::cerr << "[PlotBcompare] <2 common points – nothing to do\n"; return; }

  /* 2) Point graphs with 7 % horizontal jitter ----------------------- */
  const double dE = (eCtr.back()-eCtr.front()) /
                    std::max<std::size_t>(1,eCtr.size()-1);
  const double dx = 0.07*dE;

  auto gR = std::make_unique<TGraphErrors>(eCtr.size());
  auto gP = std::make_unique<TGraphErrors>(eCtr.size());
  for (std::size_t i=0;i<eCtr.size();++i)
  {
    gR->SetPoint(i,eCtr[i]-dx,vRMS[i]); gR->SetPointError(i,0.,sRMS[i]);
    gP->SetPoint(i,eCtr[i]+dx,vPhi[i]); gP->SetPointError(i,0.,sPhi[i]);
  }
  gR->SetMarkerStyle(24); gR->SetMarkerSize(1.2);
  gR->SetMarkerColor(kRed+1);  gR->SetLineColor(kRed+1);
  gP->SetMarkerStyle(20); gP->SetMarkerSize(1.2);
  gP->SetMarkerColor(kBlue+1); gP->SetLineColor(kBlue+1);

  auto sR = std::make_unique<TSpline3>("sRMS", gR.get());
  auto sP = std::make_unique<TSpline3>("sPhi", gP.get());
  sR->SetLineColor(kRed+1);  sR->SetLineWidth(2);
  sP->SetLineColor(kBlue+1); sP->SetLineWidth(2);

  /* 3) Fine‑grid Δb/b curve ----------------------------------------- */
  const int nFine = 300;
  auto gD  = std::make_unique<TGraph>(nFine);
  for (int i=0;i<nFine;++i)
  {
    const double x  = eCtr.front() + (eCtr.back()-eCtr.front())*i/(nFine-1);
    const double dr = (sP->Eval(x)-sR->Eval(x))/sR->Eval(x);
    gD->SetPoint(i,x,dr);
  }
  gD->SetLineColor(kBlack); gD->SetLineWidth(2);

  /* 3a) Global metrics ---------------------------------------------- */
  double sum=0., sum2=0., maxAbs=0.; int iMax=-1;
  for (int i=0;i<nFine;++i)
  { double x,d; gD->GetPoint(i,x,d); sum+=d; sum2+=d*d;
    if (std::fabs(d)>maxAbs){ maxAbs=std::fabs(d); iMax=i; } }
  const double mean = sum/nFine;
  const double rms  = std::sqrt(sum2/nFine - mean*mean);

  /* χ² with 2 % default uncertainty if none supplied */
  double chi2=0.; std::size_t ndf=0;
  for (std::size_t i=0;i<eCtr.size();++i)
  {
    const double sig = sRMS[i]>0? sRMS[i] : 0.02*vRMS[i];
    const double d   = vPhi[i]-vRMS[i];
    chi2 += d*d/(sig*sig); ++ndf;
  }

  /* 3b) Pretty terminal table --------------------------------------- */
  std::cout << "\n┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n"
            << "┃  Comparison of Ash‑b extraction methods                       ┃\n"
            << "┣━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━┫\n"
            << "┃  E   ┃  b_opt   ┃  b_phi   ┃ Δb/b [%] ┃ σ_opt ┃ σ_phi ┃\n"
            << "┣━━━━━━╋━━━━━━━━━━╋━━━━━━━━━━╋━━━━━━━━━━╋━━━━━━━━╋━━━━━━━━┫\n";
  for (std::size_t i=0;i<eCtr.size();++i)
    printf("┃%6.1f┃%10.5f┃%10.5f┃%+10.3f┃%8.4f┃%8.4f┃\n",
           eCtr[i], vRMS[i], vPhi[i], 100.*(vPhi[i]-vRMS[i])/vRMS[i],
           sRMS[i], sPhi[i]);
  std::cout << "┗━━━━━━┻━━━━━━━━━━┻━━━━━━━━━━┻━━━━━━━━━━┻━━━━━━━━┻━━━━━━━━┛\n";

  const double xMaxDev = (iMax>=0? (eCtr.front()+
                        (eCtr.back()-eCtr.front())*iMax/(nFine-1)) : 0.);
  std::cout << "\nSummary:\n"
            << "   ⟨Δb/b⟩          = " << std::setw(8) << std::fixed
            << std::setprecision(4) << 100*mean << "  %\n"
            << "   RMS(Δb/b)      = " << std::setw(8)
            << 100*rms << "  %\n"
            << "   χ² / ndf       = " << chi2 << " / " << ndf
            << "  = " << chi2/ndf << "\n"
            << "   Max |Δb/b|     = " << 100*maxAbs << "  %"
            << "  at E ≈ " << xMaxDev << " GeV\n"
            << "   Interpretation : mean within "
            << (std::fabs(100*mean)<1 ? "1 %" : "few %")
            << "; spread ~" << 100*rms << " %.\n\n";

  /* 4) Build canvas --------------------------------------------------- */
  TCanvas c("cBcompare","b comparison",900,950); c.SetTicks();

  /* Robust y‑range: median ±4×MAD  +15 % headroom                     */
  std::vector<double> tmp=vRMS; tmp.insert(tmp.end(),vPhi.begin(),vPhi.end());
  std::nth_element(tmp.begin(),tmp.begin()+tmp.size()/2,tmp.end());
  const double med = tmp[tmp.size()/2];
  std::vector<double> madV; madV.reserve(tmp.size());
  for (double v:tmp) madV.push_back(std::fabs(v-med));
  std::nth_element(madV.begin(),madV.begin()+madV.size()/2,madV.end());
  const double yLo = 0.14;
  const double yHi = 0.19;   // +15 %

  /* Upper pad (no grid) ---------------------------------------------- */
  TPad pT("pT","",0,0.32,1,1);
  pT.SetLeftMargin(0.14); pT.SetBottomMargin(0.02);
  pT.Draw(); pT.cd();

  TH2F fr("fr",";E_{slice centre}  [GeV];b  (tower units)",
           100,eCtr.front()-dE,eCtr.back()+dE, 100,yLo,yHi);
  fr.GetXaxis()->SetLabelSize(0); fr.Draw();

  sR->Draw("CSAME"); sP->Draw("CSAME");
  gR->Draw("P SAME"); gP->Draw("P SAME");

  /* ───── stats box – bottom-left of the top pad ───────────────────── */
  TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.028);

  tl.DrawLatex(0.16,0.26,
                 Form("#LT#Delta b/b#GT = %+.3f %% ", 100.*mean));

  tl.DrawLatex(0.16,0.22,
                 Form("RMS(#Delta b/b)  =  %.3f %%", 100.*rms));

  tl.DrawLatex(0.16,0.18,
                 Form("#chi^{2}/ndf      =  %.1f / %zu  = %.2f",
                      chi2, ndf, chi2/ndf));

  tl.DrawLatex(0.16,0.14,
                 Form("max |#Delta b/b|  =  %.3f %%  at  %.1f GeV",
                      100.*maxAbs, xMaxDev));

  TLegend lg(0.40,0.78,0.60,0.92);
  lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.032);
  lg.AddEntry(gR.get(),"b_{opt} (RMS)","lp");
  lg.AddEntry(gP.get(),"b_{#varphi} (fit)","lp");
  lg.Draw();

  /* Lower pad – relative diff ---------------------------------------- */
  c.cd();
  TPad pB("pB","",0,0,1,0.32);
  pB.SetLeftMargin(0.14); pB.SetTopMargin(0.02);
  pB.SetBottomMargin(0.30); pB.Draw(); pB.cd();

  TH2F fr2("fr2",";E_{slice centre}  [GeV];(b_{#varphi}-b_{opt}) / b_{opt}",
            100,eCtr.front()-dE,eCtr.back()+dE, 100,-0.25,0.25);
  fr2.Draw();

  TBox band(eCtr.front()-dE,-0.05,eCtr.back()+dE,0.05);
  band.SetFillStyle(3004); band.SetFillColor(kGray+1); band.Draw("same");

  gD->Draw("L SAME");

  /* 5) Save to file --------------------------------------------------- */
  TString png = TString::Format("%s/b_compare_RMS_vs_fit.png",outDir);
  c.SaveAs(png);
  std::cout << "[PlotBcompare]  canvas saved to  " << png << '\n';

  TFile fout(Form("%s/b_compare.root",outDir),"RECREATE");
  gR->Write("gRMS"); gP->Write("gPhi"); gD->Write("gRelDiff");
  sR->Write("sRMS"); sP->Write("sPhi");  fout.Close();
}

/* ========================================================================== */
/*  drawLego3D – single view  + clean 2×2 multi‑view                          */
/*     • ORIGINAL png     :  lego_<tag>.png                                   */
/*     • NEW four‑view png:  lego_<tag>_views.png                              */
/* ========================================================================== */
void drawLego3D(TH3F*        h,
                const char*  tag,          // “unc”, “cor”, …
                const char*  hdr,          // header text
                const char*  outDir,       // directory for PNGs
                double       fTtl,         // title‑font size
                double       fLbl)         // axis‑label font size
{
    if (!h) return;

    /* ------------------------------------------------------------------ *
     * 0. Cosmetics common to ALL drawings
     * ------------------------------------------------------------------ */
    h->SetMinimum(h->GetMinimum());
    h->SetMaximum(h->GetMaximum());

    h->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    h->GetYaxis()->SetTitle("block #phi_{local, 2#times 2}");
    h->GetZaxis()->SetTitle("Cluster E  [GeV]");

    h->GetXaxis()->SetLabelSize(fLbl - 0.002);
    h->GetYaxis()->SetLabelSize(fLbl - 0.002);
    h->GetZaxis()->SetLabelSize(fLbl - 0.002);

    /* ****************************************************************** *
     * 1.  ORIGINAL single‑view canvas (unchanged)
     * ****************************************************************** */
    TCanvas cMain(Form("c_%s", tag), "", 1400, 1080);
    cMain.SetRightMargin(0.15);
    cMain.SetBottomMargin(0.14);

    h->Draw("LEGO2Z0");
    cMain.Update();                                                // palette

    /* ------------------------------------------------------------------ *
     * helper: format palette labels as ×10^e  (robust against ROOT 6.30+)
     * ------------------------------------------------------------------ */
    auto prettifyPalette = [&](TPaletteAxis* pal)
    {
        if (!pal) return;

        /* tight bar ----------------------------------------------------- */
        pal->SetX1NDC(0.905);  pal->SetX2NDC(0.925);
        pal->SetBorderSize(0);  pal->SetFillStyle(0);

        /* scientific‑notation header ------------------------------------ */
        double zMax = h->GetMaximum();        // ROOT‑version‑proof
        int    e    = static_cast<int>(std::floor(std::log10(zMax)));
        double scale = std::pow(10., e);

        pal->GetAxis()->SetTitle( Form("counts ×10^{%d}", e) );
        pal->GetAxis()->CenterTitle();
        pal->GetAxis()->SetTitleOffset(0.60);

        /* relabel individual ticks -------------------------------------- *
         * Works with both TAxis (≤6.28) and TGaxis (≥6.30).               */
        if (auto* axTA = dynamic_cast<TAxis*>(pal->GetAxis()))
        {
            int nb = axTA->GetNbins();
            for (int i = 1; i <= nb; ++i)
            {
                double val = axTA->GetBinLowEdge(i) / scale;
                axTA->ChangeLabel(i, -1, -1, -1, -1, -1, Form("%.0f", val));
            }
        }
        else if (auto* axTG = dynamic_cast<TGaxis*>(pal->GetAxis()))
        {
            /* TGaxis no longer exposes bin‑edge helpers; rebuild labels   *
             * from axis world limits.                                     */
            double wmin = axTG->GetWmin();
            double wmax = axTG->GetWmax();
            int    nDiv = axTG->GetNdiv()%100;       // major divisions
            double step = (wmax - wmin) / nDiv;

            for (int i = 0; i <= nDiv; ++i)
            {
                double valLab = (wmin + i*step) / scale;
                axTG->ChangeLabel(i+1, -1, -1, -1, -1, -1,
                                  Form("%.0f", valLab));
            }
        }
    };

    if (auto* pal = dynamic_cast<TPaletteAxis*>
          (h->GetListOfFunctions()->FindObject("palette")))
    {
        pal->SetLabelSize(fLbl + 0.005);
        prettifyPalette(pal);     // <-- new, version‑proof helper
    }

    TLatex tit; tit.SetNDC(); tit.SetTextFont(42); tit.SetTextAlign(13);
    tit.SetTextSize(fTtl);
    tit.DrawLatex(0.12, 0.94, Form("#bf{%s}", hdr));
    tit.SetTextSize(fTtl - 0.010);
    tit.DrawLatex(0.12, 0.90, Form("entries: %.0f", h->GetEntries()));

    cMain.Print(Form("%s/lego_%s.png", outDir, tag), "png 600");

    /* ********************************************************************** *
     * 2.  E I G H T – V I E W   S H E E T   (4 × 2 pads)
     *
     *    ┌────────────────┬────────────────┬────────────────┬────────────────┐
     *    │  φ =   0°      │  φ =  90°      │  φ = 180°      │  φ = 270°      │ ← row 1
     *    │  (tilted)      │  (tilted)      │  (tilted)      │  (tilted)      │
     *    ├────────────────┼────────────────┼────────────────┼────────────────┤
     *    │  Front face    │  Right face    │  Back face     │  Left face     │ ← row 2
     *    │  θ = 12°       │  θ = 12°       │  θ = 12°       │  θ = 12°       │
     *    └────────────────┴────────────────┴────────────────┴────────────────┘
     * ********************************************************************** */

    /* ---- reference camera from the single‑view canvas (for row 1) ---- */
    cMain.cd();
    const double thetaTilt = gPad->GetTheta();   // usually ≈ 30°
    const double phiTilt   = gPad->GetPhi();     // usually ≈ 45°

    /* ---- camera angles for the grid ---------------------------------- */
    const double phiAbs[4]  = {   0.0,   90.0,  180.0,  270.0 };
    const char*  faceTxt[4] = { "Front", "Right", "Back", "Left" };
    const double thetaFace  = 12.0;              // gentle “bow‑down”

    /* ---- canvas ------------------------------------------------------- */
    TCanvas c8(Form("c_%s_views", tag),
               Form("%s – tilted & side‑on views", hdr),
               2100, 1150);
    c8.Divide(4, 2, 0.004, 0.010);               // almost invisible gutters

    /* permanent clone – avoids touching the original histogram ---------- */
    std::unique_ptr<TH3F> hMaster( static_cast<TH3F*>(h->Clone()) );
    hMaster->SetDirectory(nullptr);

    for (int row = 0; row < 2; ++row)
    for (int col = 0; col < 4; ++col)
    {
        const int ip = 1 + row*4 + col;          // pad index (1‑based)
        c8.cd(ip);

        gPad->SetRightMargin (0.13);
        gPad->SetBottomMargin(0.12);
        gPad->SetBorderMode  (0);

        /* ---------------- draw ---------------------------------------- */
        auto* hDraw = static_cast<TH3F*>(hMaster->Clone());
        hDraw->SetDirectory(nullptr);
        hDraw->Draw("LEGO2Z0");

        /* ---------------- camera -------------------------------------- */
        const double th = (row == 0) ? thetaTilt : thetaFace;
        const double ph = (row == 0)
                          ? phiTilt + phiAbs[col]   // keep ROOT’s default skew
                          : phiAbs[col];            // strict front / side / back
        gPad->SetTheta(th);
        gPad->SetPhi  (std::fmod(ph + 360.0, 360.0));
        gPad->Modified(); gPad->Update();

        /* keep colour bar only in pad 1 -------------------------------- */
        if (ip != 1)
            for (auto* obj : *gPad->GetListOfPrimitives())
                if (obj->InheritsFrom(TPaletteAxis::Class()))
                { obj->Delete(); break; }
    }

    /* ---------------- save -------------------------------------------- */
    c8.Print(Form("%s/lego_%s_views.png", outDir, tag), "png 600");
}



/* ======================================================================= */
/*  auditResidual – residual |ΔE|/⟨E⟩ map + η / φ profiles                 */
/*                                                                          *
 *  Z‑value plotted: |ΔE| / ⟨E⟩ • 100 %                                    *
 *  – ΔE  =  mean cluster‑energy in block − global mean                     *
 *  – ⟨E⟩ =  global mean cluster‑energy                                     *                                          *
 *  ‣ Stores the three finished objects (2‑D, φ‑profile, η‑profile)         *
 *    in a static cache keyed by tag.                                       *
 *  ‣ As soon as BOTH “UNCORRECTED” and “CORRECTED” are present, it          *
 *    automatically writes three two‑pad comparison canvases.               */
/* ======================================================================= */
void auditResidual(TH3F*  h,
                   const char* tag,          // "UNCORRECTED" | "CORRECTED"
                   const char* outDir,
                   double fTtl,
                   double fLbl)
{
  if (!h) return;

/* ---------------------------------------------------------------------- */
/* 0.  small helper to put a tag in the top‑right of the *current pad*    */
/* ---------------------------------------------------------------------- */
  auto drawTag = [&](const char* txt)
  {
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(fTtl);
    t.SetTextAlign(33);         // top‑right
    t.DrawLatex(0.85,0.92,txt);
  };

/* ---------------------------------------------------------------------- */
/* 1. build residual map                                                  */
/* ---------------------------------------------------------------------- */
  std::unique_ptr<TH2D> meanMap( h->Project3DProfile("xy") );
  const double mu = meanMap->GetMean();

  std::unique_ptr<TH2D> res(
      static_cast<TH2D*>( meanMap->Clone(Form("res_%s",tag)) ) );

  const int nx = res->GetXaxis()->GetNbins();
  const int ny = res->GetYaxis()->GetNbins();
  for (int ix=1; ix<=nx; ++ix)
    for (int iy=1; iy<=ny; ++iy)
      res->SetBinContent(ix,iy,
        100.*std::fabs(meanMap->GetBinContent(ix,iy)-mu)/mu);

  res->SetMinimum(0.0);
  res->SetMaximum(res->GetMaximum());      // just to freeze the range

    /* ---------------------------------------------------------------------- */
    /* 2. 2-D colour map (single histogram)                                   */
    /* ---------------------------------------------------------------------- */

    /* ---------- 2.a  colour-scale bookkeeping ----------------------------- *
     * We lock the palette range to the first (“UNCORRECTED”) call and reuse  *
     * it for the later “CORRECTED” call, so both plots share one scale.      */
    static double uncorrZmax = -1.;                       // remembered once
    const bool   isUnc = (std::strcmp(tag,"UNCORRECTED")==0);
    const bool   isCor = (std::strcmp(tag,"CORRECTED"  )==0);

    double zMax;
    if (isUnc) {                                         // first call
      zMax        = res->GetMaximum();
      uncorrZmax  = zMax;                                // lock-in
    } else {                                             // later call
      zMax = (uncorrZmax > 0 ? uncorrZmax : res->GetMaximum());
    }
    res->SetMinimum(0.0);
    res->SetMaximum(zMax);                               // unified scale

    /* ---------- 2.b  canvas & draw ---------------------------------------- */
    TPaletteAxis *pal = nullptr;

    TCanvas c2(Form("c_res_%s",tag),"",1000,850);
    c2.SetLeftMargin(0.18);  c2.SetBottomMargin(0.16);
    c2.SetRightMargin(0.12); c2.SetTopMargin(0.08);      // extra headroom

    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);

    res->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    res->GetYaxis()->SetTitle("block #phi_{local, 2#times 2}");
    res->GetZaxis()->SetTitle("");                      // blank vertical title
    res->Draw("COLZ");
    gPad->Update();                                     // palette now exists

    pal = dynamic_cast<TPaletteAxis*>(
             res->GetListOfFunctions()->FindObject("palette"));

    if (pal) {
      pal->SetX1NDC(0.92);  pal->SetX2NDC(0.945);
      pal->GetAxis()->SetTitle("");
      pal->SetLabelSize((fLbl-0.002)*0.75);
      pal->SetBorderSize(0); pal->SetFillStyle(0);

      /* ---------- horizontal Z-axis label (bold, larger) ------------------ */
      const double xMid = 0.5*(pal->GetX1NDC() + pal->GetX2NDC());
      TLatex tz; tz.SetNDC();
      tz.SetTextFont(62);        // 62 = Helvetica-Bold
      tz.SetTextSize(0.028);     // larger
      tz.SetTextAlign(21);       // centred
      tz.DrawLatex(xMid, pal->GetY2NDC() + 0.025,
                   "|#DeltaE| / #LT E #GT  [%]");
    }

    drawTag(tag);
    c2.Print(Form("%s/residual2_%s.png",outDir,tag),"png 600");

/* ---------------------------------------------------------------------- */
/* 3. build φ and η profiles                                              */
/* ---------------------------------------------------------------------- */
  std::unique_ptr<TProfile> pPhi( res->ProfileX(Form("pPhi_%s",tag)) );
  std::unique_ptr<TProfile> pEta( res->ProfileY(Form("pEta_%s",tag)) );

  auto drawProfile = [&](TProfile* p,const char* ax)
  {
    TCanvas cp(Form("c_%s_%s",ax,tag),"",900,600);
    cp.SetLeftMargin(0.18); cp.SetBottomMargin(0.15);
    p->SetLineWidth(2); p->SetMarkerStyle(kFullCircle); p->SetMarkerSize(1.1);
    p->GetYaxis()->SetTitle("|#DeltaE| / #LT E #GT  [%]");
    if (ax[0]=='p')
      p->GetXaxis()->SetTitle("block #phi_{local, 2#times 2}");
    else
      p->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    p->Draw("P");
    /* mean line */
    TLine L(p->GetXaxis()->GetXmin(),p->GetMean(2),
            p->GetXaxis()->GetXmax(),p->GetMean(2));
    L.SetLineColor(kGray+2); L.Draw();
    drawTag(tag);
    cp.Print(Form("%s/residual_%s_%s.png",outDir,ax,tag),"png 600");
  };

  drawProfile(pPhi.get(),"phi");
  drawProfile(pEta.get(),"eta");

    /* ------------------------------------------------------------------ */
    /* 4.  STATIC cache + on‑the‑fly comparisons                           */
    /* ------------------------------------------------------------------ */
    struct CacheEntry
    {
      std::unique_ptr<TH2D>     res2D;
      std::unique_ptr<TProfile> pPhi, pEta;
    };
    static std::map<std::string,CacheEntry> cache;

    /* --- move the freshly‑built histograms into the cache -------------- */
    CacheEntry ce;
    ce.res2D.swap(res);
    ce.pPhi .swap(pPhi);
    ce.pEta .swap(pEta);
    cache[tag] = std::move(ce);

    /* need both tags before we can draw any comparison ------------------ */
    if (cache.count("UNCORRECTED")==0 || cache.count("CORRECTED")==0) return;

    /* ================================================================== */
    /*  helper – side‑by‑side 2‑D residual maps + stand‑alone palette     */
    /*  (identical look to the single‑histogram “COLZ” palette)           */
    /* ================================================================== */
    auto makeSide2D = [&](TH2 *unc, TH2 *cor, const char *outPng)
    {
       std::cout << "[makeSide2D] >>> start, writing " << outPng << '\n';

       /* ---------- 0. sanity -------------------------------------------- */
       if (!unc || !cor)                   { std::cerr << "[ERROR] null hist\n"; return; }
       if (unc->GetNbinsX()==0||cor->GetNbinsX()==0){ std::cerr<<"[ERROR] empty\n";return;}

       const double zMax = std::max(unc->GetMaximum(), cor->GetMaximum());
       if (!std::isfinite(zMax))           { std::cerr << "[ERROR] bad zMax\n"; return; }
       for (TH2 *h : {unc,cor}) { h->SetMinimum(0.0); h->SetMaximum(zMax); }

       /* ---------- 1. canvas & pads ------------------------------------- */
       const int W = 1900 , H = 850 ;
       TCanvas c("cSide2D","side‑by‑side residual maps",W,H);

       TPad pUnc("pUnc","UNC",0.00,0.00,0.46,1.00);
       TPad pCor("pCor","COR",0.46,0.00,0.92,1.00);
       TPad pPal("pPal","PAL",0.92,0.00,1.00,1.00);
       pUnc.Draw(); pCor.Draw(); pPal.Draw();

       gStyle->SetPalette(kBird);
       gStyle->SetNumberContours(255);

       /* ---------- 2. UNCORRECTED map ----------------------------------- */
       pUnc.cd();  pUnc.SetLeftMargin(.18); pUnc.SetRightMargin(.02); pUnc.SetBottomMargin(.16);
       unc->Draw("COL");  drawTag("UNCORRECTED");

       /* ---------- 3. CORRECTED map ------------------------------------- */
       pCor.cd();  pCor.SetLeftMargin(.18); pCor.SetRightMargin(.02); pCor.SetBottomMargin(.16);
       cor->Draw("COL");  drawTag("CORRECTED");

       /* ----------------------------------------------------------------- */
       /* 4. palette pad – manufacture palette with one invisible pixel     */
       /* ----------------------------------------------------------------- */
       pPal.cd();
       pPal.SetLeftMargin  (0.00);
       pPal.SetRightMargin (0.00);
       pPal.SetTopMargin   (0.12);
       pPal.SetBottomMargin(0.18);
       pPal.SetFrameLineWidth(0);          // no frame at all
       pPal.SetTickx(0); pPal.SetTicky(0); // no pad ticks

       /* --- helper histogram: 1×1 bin, same z‑range, fully transparent --- */
       TH2D hDummy("hDummy","",1,0,1,1,0,1);
       hDummy.SetMinimum(0.0); hDummy.SetMaximum(zMax);
       hDummy.SetStats(0);
       hDummy.SetLineWidth(0); hDummy.SetFillStyle(0);
       hDummy.GetXaxis()->SetTickLength(0);
       hDummy.GetYaxis()->SetTickLength(0);
       hDummy.GetXaxis()->SetLabelSize(0);
       hDummy.GetYaxis()->SetLabelSize(0);

       hDummy.Draw("COLZ0");              // palette generated, no frame/axes
       gPad->Update();

       /* retrieve palette ------------------------------------------------- */
       auto *pal = dynamic_cast<TPaletteAxis*>(
                     hDummy.GetListOfFunctions()->FindObject("palette"));
       if (!pal) { std::cerr << "[ERROR] palette axis not found\n"; return; }

       /* --- 4.a resize palette to exactly 0.92–0.945 canvas X ---------- */
       const double xAbs1 = 0.92 , xAbs2 = 0.945;      // absolute canvas NDC
       const double padX1 = pPal.GetXlowNDC();
       const double padW  = pPal.GetWNDC();
       const double x1    = (xAbs1 - padX1)/padW;       // → pPal’s local NDC
       const double x2    = (xAbs2 - padX1)/padW;
       pal->SetX1NDC(x1);  pal->SetX2NDC(x2);

       /* cosmetic clone from single‑plot palette ------------------------- */
       pal->GetAxis()->SetTitle("");
       pal->SetLabelSize((fLbl-0.002)*0.75);
       pal->SetBorderSize(0);
       pal->SetFillStyle(0);

       /* --- 4.b external title (same size/offset) ----------------------- */
       const double xMid = 0.5*(pal->GetX1NDC()+pal->GetX2NDC());
       TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.022); t.SetTextAlign(21);
       t.DrawLatex(xMid, pal->GetY2NDC() + 0.018, "|#DeltaE| / #LT E #GT  [%]");

       gPad->Modified(); gPad->Update();

       /* ---------- 5. export -------------------------------------------- */
       c.SaveAs(outPng);
       if (!gSystem->AccessPathName(outPng))
             std::cout << "[makeSide2D] <<< finished OK – file written\n";
       else std::cerr << "[ERROR] SaveAs produced no file!\n";
    };


    /* ================================================================== */
    /*  helper #2 – overlay of 1‑D profiles                               */
    /* ================================================================== */
    auto overlay1D = [&](TProfile* unc, TProfile* cor,
                         const char* xTit, const char* outPng)
    {
      TCanvas c("cOver1D","",900,700);
      c.SetLeftMargin(0.18);  c.SetBottomMargin(0.16);

      /* aesthetics ----------------------------------------------------- */
      unc->SetLineWidth(2);   unc->SetMarkerStyle(kFullCircle);
      unc->SetMarkerColor(kBlue+1);     unc->SetLineColor(kBlue+1);

      cor->SetLineWidth(2);   cor->SetMarkerStyle(kFullSquare);
      cor->SetMarkerColor(kRed+1);      cor->SetLineColor(kRed+1);

      unc->GetXaxis()->SetTitle(xTit);
      unc->GetYaxis()->SetTitle("|#DeltaE| / #LT E #GT  [%]");

      /* dynamic Y‑range to fit both graphs ----------------------------- */
      const double yMin = std::min(unc->GetMinimum(),cor->GetMinimum())*0.95;
      const double yMax = std::max(unc->GetMaximum(),cor->GetMaximum())*1.05;
      unc->GetYaxis()->SetRangeUser(std::max(0.,yMin),yMax);

      unc->Draw("P");     cor->Draw("P SAME");

      TLegend L(0.8,0.85,0.9,0.92);
      L.SetBorderSize(0);  L.SetFillStyle(0); L.SetTextSize(0.022);
      L.AddEntry(unc,"UNCORRECTED","lp");
      L.AddEntry(cor,"CORRECTED"  ,"lp");
      L.Draw();

      c.Print(outPng,"png 600");
    };

    /* ================================================================== */
    /* 5.  produce comparison files                                       */
    /* ================================================================== */
    makeSide2D(cache["UNCORRECTED"].res2D.get(),
               cache["CORRECTED" ].res2D.get(),
               Form("%s/residual2_SIDE.png",outDir));

    overlay1D(cache["UNCORRECTED"].pPhi.get(),
              cache["CORRECTED" ].pPhi.get(),
              "block #phi_{local, 2#times 2}",
              Form("%s/residual_phi_OVER.png",outDir));

    overlay1D(cache["UNCORRECTED"].pEta.get(),
              cache["CORRECTED" ].pEta.get(),
              "block #eta_{local, 2#times 2}",
              Form("%s/residual_eta_OVER.png",outDir));

}
/* ======================================================================= */


/* ════════════════════════════════════════════════════════════════════ *
 *  High‑quality LEGO3D spin‑gif generator  –  lightweight output      *
 * ════════════════════════════════════════════════════════════════════ */
void makeLegoGifHD(TH3          *h1,                /* mandatory          */
                    const char   *tag1,
                    const char   *hdr1,
                    TH3          *h2      = nullptr, /* optional ‑ morph   */
                    const char   *tag2    = "",
                    const char   *hdr2    = "",
                    const char   *outDir  = ".",
                    int   nFrames         = 180,
                    double theta0         = 28.0,
                    double theta1         = 38.0,
                    double phi0           =   0.0,
                    double phiArc         = 360.0,
                    int   SS              = 2,
                    int   W0              = 1280,
                    int   H0              = 720,
                    int   FPS             = 25)
 {
     /* ------------------------------------------------------------------ */
     /* 0. defend against misuse                                           */
     /* ------------------------------------------------------------------ */
     if (!h1) { Error("makeLegoGifHD","h1 null"); return; }
     if (nFrames < 2) { Error("makeLegoGifHD","need ≥2 frames"); return; }
     if (h2 && ( h2->GetNbinsX()!=h1->GetNbinsX() ||
                 h2->GetNbinsY()!=h1->GetNbinsY() ||
                 h2->GetNbinsZ()!=h1->GetNbinsZ() ))
     {
         Error("makeLegoGifHD","h1 / h2 binning mismatch – no morph GIF");
         h2 = nullptr;
     }

 #if ROOT_VERSION_CODE >= ROOT_VERSION(6,30,0)
     ROOT::EnableImplicitMT(0);
 #endif
     gROOT->SetBatch(kTRUE);
     gStyle->SetOptStat(0);
     gStyle->SetNumberContours(140);
     gStyle->SetCanvasPreferGL(true);
     gErrorIgnoreLevel = kWarning;

     /* small helper – everything that renders one GIF ------------------- */
     auto renderOne = [&](TH3 *h, const char *tag, const char *hdr,
                          const char *suffix, bool downScale)
     {
         const int W = W0*SS, H = H0*SS;
         std::unique_ptr<TCanvas> c(
             new TCanvas(Form("c_%s_%s",tag,suffix), "", W, H));
         c->SetRightMargin(.15);  c->SetBottomMargin(.13);
         c->SetFillColorAlpha(kWhite,0);

         /* static scene -------------------------------------------------- */
         h->GetXaxis()->SetTitleOffset(1.25);
         h->GetYaxis()->SetTitleOffset(1.55);
         h->GetZaxis()->SetTitleOffset(0.90);
         h->Draw("LEGO2Z0");
         gPad->Modified(); gPad->Update();

         if (auto *pal =
             dynamic_cast<TPaletteAxis*>(h->GetListOfFunctions()
                                               ->FindObject("palette")))
         {
             pal->SetX1NDC(0.86); pal->SetX2NDC(0.89);
             pal->SetBorderSize(0); pal->SetFillStyle(0);
             pal->SetLabelSize(0.022*SS/2.0);
         }
         const Long64_t nEnt = h->GetEntries();
         TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextAlign(13);
         t.SetTextSize(0.05*SS/2.0);
         t.DrawLatex(.04,.965,Form("#bf{%s  (N = %lld)}",hdr,(Long64_t)nEnt));
         gPad->Modified(); gPad->Update();

         /* dirs & filenames ---------------------------------------------- */
         TString base   = Form("%s/lego_%s%s", outDir, tag, suffix);
         TString pngDir = base + "_png"; TString gifOut = base + ".gif";
         gSystem->mkdir(pngDir,kTRUE);

         /* frames -------------------------------------------------------- */
         const char spin[4]={'|','/','-','\\'}; const TDatime t0;
         for(int i=0;i<nFrames;++i)
         {
             double u   = double(i)/(nFrames-1);
             gPad->SetPhi (phi0 + phiArc*u);
             gPad->SetTheta(theta0 + (theta1-theta0)*
                            0.5*(1-std::cos(u*TMath::Pi())));
             gPad->Modified(); gPad->Update();
             c->Print(Form("%s/f%04d.png",pngDir.Data(),i),"png");
             if(i%std::max(1,nFrames/100)==0||i==nFrames-1)
             {
                 int e=TDatime().Convert()-t0.Convert();
                 int eta=e*(nFrames-i-1)/std::max(1,i+1);
                 printf("\r  %c %3d %%  %4d s ETA",
                        spin[i&3], int(100.*(i+1)/nFrames), eta); fflush(stdout);
             }
         } puts("\r  ✔ frame rendering done.          ");

         /* encoding ------------------------------------------------------ */
         TString ff=gSystem->Which(nullptr,"ffmpeg");
         if(ff.IsNull()
            && !gSystem->AccessPathName("/opt/homebrew/bin/ffmpeg",kExecutePermission))
             ff="/opt/homebrew/bin/ffmpeg";
         if(ff.IsNull()
            && !gSystem->AccessPathName("/usr/local/bin/ffmpeg",kExecutePermission))
             ff="/usr/local/bin/ffmpeg";
         if(ff.IsNull()){
             Error("makeLegoGifHD","ffmpeg not found"); return;
         }

         TString pal = pngDir+"/pal.png";
         TString cmd = Form("%s -loglevel error -y -i %s/f%%04d.png "
                            "-vf \"scale=%d:%d:flags=lanczos,palettegen=max_colors=128\" %s",
                            ff.Data(), pngDir.Data(), W0, H0, pal.Data());
         if(gSystem->Exec(cmd)) { Error("ffmpeg","palettegen failed"); return; }

         cmd = Form("%s -loglevel error -y -framerate %d -i %s/f%%04d.png -i %s "
                    "-lavfi \"scale=%d:%d:flags=lanczos[x];[x][1:v]paletteuse"
                    "=dither=sierra2_4a\" -gifflags +transdiff -r %d %s",
                    ff.Data(), FPS, pngDir.Data(), pal.Data(),
                    W0, H0, FPS, gifOut.Data());
         if(gSystem->Exec(cmd)){ Error("ffmpeg","paletteuse failed"); return; }

         gSystem->Exec(Form("rm -rf %s",pngDir.Data()));
         gSystem->Unlink(pal);
         Printf("  ↪  %s", gifOut.Data());
     };

     /* ------------------------------------------------------------------ *
      * 1. GIF #1  (h1)                                                    *
      * ------------------------------------------------------------------ */
     renderOne(h1, tag1, hdr1, "_spin", true);

     /* ------------------------------------------------------------------ *
      * 2. GIF #2  (h2) – only if supplied                                 *
      * ------------------------------------------------------------------ */
     if (h2) renderOne(h2, tag2, hdr2, "_spin", true);

    /* ------------------------------------------------------------------ *
     * 3.  MORPH GIF  :  raw → corrected                                  *
     * ------------------------------------------------------------------ */
    if (!h2) return;            // user did not request a morph

    /* --- 3.a  common Z‑range ------------------------------------------- */
    const double zMin = std::min(h1->GetMinimum(), h2->GetMinimum());
    const double zMax = std::max(h1->GetMaximum(), h2->GetMaximum());

    /* --- 3.b  template histogram --------------------------------------- */
    std::unique_ptr<TH3> hMix( static_cast<TH3*>(h1->Clone("h_mix")) );
    hMix->SetDirectory(nullptr);
    hMix->SetMinimum(zMin);
    hMix->SetMaximum(zMax);

    /* --- 3.c  canvas & folder ------------------------------------------ */
    const int W = W0*SS , H = H0*SS;
    std::unique_ptr<TCanvas> c(
        new TCanvas(Form("c_%s_%s_morph",tag1,tag2),"",W,H));
    c->SetRightMargin(.15); c->SetBottomMargin(.13);
    c->SetFillColorAlpha(kWhite,0);

    TString base   = Form("%s/lego_%s-%s_morph", outDir, tag1, tag2);
    TString pngDir = base + "_png";
    gSystem->mkdir(pngDir,kTRUE);

    /* --- 3.d  helper: encode a PNG folder → GIF ------------------------ */
    auto encodeGif = [&](const TString &dir, const TString &gifOut)
    {
        TString ff = gSystem->Which(nullptr,"ffmpeg");
        if (ff.IsNull() &&
            !gSystem->AccessPathName("/opt/homebrew/bin/ffmpeg",kExecutePermission))
            ff="/opt/homebrew/bin/ffmpeg";
        if (ff.IsNull() &&
            !gSystem->AccessPathName("/usr/local/bin/ffmpeg",kExecutePermission))
            ff="/usr/local/bin/ffmpeg";
        if (ff.IsNull()) { Error("makeLegoGifHD","ffmpeg not found"); return; }

        TString pal = dir + "/pal.png";
        TString cmd = Form("%s -loglevel error -y -i %s/f%%04d.png "
                           "-vf \"scale=%d:%d:flags=lanczos,palettegen=max_colors=128\" %s",
                           ff.Data(), dir.Data(), W0, H0, pal.Data());
        if (gSystem->Exec(cmd)) { Error("ffmpeg","palettegen failed"); return; }

        cmd = Form("%s -loglevel error -y -framerate %d -i %s/f%%04d.png -i %s "
                   "-lavfi \"scale=%d:%d:flags=lanczos[x];[x][1:v]paletteuse="
                   "dither=sierra2_4a\" -gifflags +transdiff -r %d %s",
                   ff.Data(), FPS, dir.Data(), pal.Data(),
                   W0, H0, FPS, gifOut.Data());
        if (gSystem->Exec(cmd)) { Error("ffmpeg","paletteuse failed"); return; }

        gSystem->Unlink(pal);
        gSystem->Exec(Form("rm -rf %s", dir.Data()));
    };

    /* --- 3.e  render frames -------------------------------------------- */
    const char spin[4] = {'|','/','-','\\'};
    const TDatime t0;

    for (int i=0; i<nFrames; ++i)
    {
        const double u  = double(i)/(nFrames-1);   // 0 … 1
        const double a  = u;                       // linear blend weight

        /* blend h1 → h2 --------------------------------------------------- */
        const int nX=h1->GetNbinsX(), nY=h1->GetNbinsY(), nZ=h1->GetNbinsZ();
        for (int ix=1; ix<=nX; ++ix)
          for (int iy=1; iy<=nY; ++iy)
            for (int iz=1; iz<=nZ; ++iz)
            {
                const double v = (1.0-a)*h1->GetBinContent(ix,iy,iz)
                               +        a *h2->GetBinContent(ix,iy,iz);
                hMix->SetBinContent(ix,iy,iz, v);
            }

        /* first frame: draw & keep palette -------------------------------- */
        if (i==0)
            hMix->Draw("LEGO2Z0");
        else
            hMix->Draw("LEGO2Z0 SAME");          // re‑uses existing palette

        if (i==0)                                // palette cosmetics once
        {
            if (auto *pal = dynamic_cast<TPaletteAxis*>(
                    hMix->GetListOfFunctions()->FindObject("palette")))
            {
                pal->SetX1NDC(0.86); pal->SetX2NDC(0.89);
                pal->SetBorderSize(0); pal->SetFillStyle(0);
                pal->SetLabelSize(0.022*SS/2.0);
            }
        }

        /* camera ---------------------------------------------------------- */
        const double phi = phi0 + phiArc*u;
        const double th  = theta0 + (theta1-theta0)*
                           0.5*(1-std::cos(u*TMath::Pi()));
        gPad->SetPhi(phi);  gPad->SetTheta(th);

        /* title ----------------------------------------------------------- */
        TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextAlign(13);
        tl.SetTextSize(0.05*SS/2.0);
        tl.DrawLatex(0.04,0.965, Form("#bf{MORPH  %3.0f %%}", 100.*a));

        gPad->Modified(); gPad->Update();
        c->Print(Form("%s/f%04d.png", pngDir.Data(), i), "png");

        if (i%std::max(1,nFrames/100)==0 || i==nFrames-1)
        {
            int e   = TDatime().Convert() - t0.Convert();
            int eta = e*(nFrames-i-1)/std::max(1,i+1);
            printf("\r  %c %3d %%  %4d s ETA",
                   spin[i&3], int(100.*(i+1)/nFrames), eta); fflush(stdout);
        }
    }
    puts("\r  ✔ morph frames rendered.            ");

    /* --- 3.f  encode & clean‑up ---------------------------------------- */
    encodeGif(pngDir, base + ".gif");
    Printf("  ↪  %s.gif", base.Data());
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
    const std::string dirMeanSigmaGaussCompare= std::string(outDir) + "/mean_sigmaGauss_HistCompare";
    const std::string dirIndivPerVar          = std::string(outDir) + "/invididualPerVariantSummary";
    const std::string dirEtaPhi2D             = std::string(outDir) + "/2DetaDphiDep";
    const std::string dirFirstBinOverlaysOnly = std::string(outDir) + "/firstBinOverlaysOnly";
    ensureDir(dirOverlayNoCorr.c_str());
    ensureDir(dirSeparateNoCorr.c_str());
    ensureDir(dirRawVsCorr.c_str());
    ensureDir(dirScratchRawVsCorr.c_str());
    ensureDir(dirFSvsClusNoCorr.c_str());
    ensureDir(dirFSvsClusWithCorr.c_str());
    ensureDir(dirMeanSigmaGaussCompare.c_str());
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


    // ----------------------- CSV helpers ----------------------------

    // Replace “.png” (or any extension) with “.csv”
    auto pngToCsvPath = [](const std::string& png)->std::string {
        const auto pos = png.find_last_of('.');
        return (pos == std::string::npos) ? (png + ".csv")
                                          : (png.substr(0,pos) + ".csv");
    };

    // Quote CSV fields that may contain commas
    auto csvQuote = [](const std::string& s)->std::string{
        std::string t = s;
        for (std::size_t i=0; i<t.size(); ++i)
            if (t[i] == '"') { t.insert(i, 1, '"'); ++i; }
        return "\"" + t + "\"";
    };

    // Find energy-bin edges for a given center (nearest)
    auto edgesFromCenter = [&](double ctr)->std::pair<double,double>{
        if (eEdges.empty()) return {0.0,0.0};
        double best = 1e99;
        std::pair<double,double> out = eEdges.front();
        for (const auto& pr : eEdges){
            const double c = 0.5*(pr.first + pr.second);
            const double d = std::fabs(c - ctr);
            if (d < best){ best = d; out = pr; }
        }
        return out;
    };

    // Detect residual kind from a vector of hists (“phi” / “eta”)
    auto detectResidual = [](const std::vector<TH1F*>& V)->std::string{
        for (auto* h : V){
            if (!h) continue;
            const TString n = h->GetName();
            if (n.Contains("phi")) return "phi";
            if (n.Contains("eta")) return "eta";
        }
        return "residual";
    };

    // Variant mapping from histogram name (used by FS vs Clus overlays)
    auto variantFromHistName = [](const std::string& n)->std::string{
        if (n.find("cpBcorr") != std::string::npos) return "b(E) corr, cluster";
        if (n.find("cpCorr")  != std::string::npos) return "CorrectPosition, cluster";
        if (n.find("cpRaw")   != std::string::npos) return "no corr, cluster";
        if (n.find("_corr_")  != std::string::npos) return "b(E) corr, scratch";
        if (n.find("_raw_")   != std::string::npos) return "no corr, scratch";
        return "variant";
    };

    auto firstHistName = [](const std::vector<TH1F*>& V)->std::string{
        for (auto* h : V) if (h) return std::string(h->GetName());
        return {};
    };

    // Core writer: one row per point per series
    auto saveMuSigmaCSV = [&](const std::string& pngPath,
                              const std::vector<std::string>& variants,
                              const std::vector<std::string>& residuals,
                              const std::vector<std::vector<double>>& eCtrS,
                              const std::vector<std::vector<double>>& muS,
                              const std::vector<std::vector<double>>& dmuS,
                              const std::vector<std::vector<double>>& sgS,
                              const std::vector<std::vector<double>>& dsgS)
    {
        const std::string csv = pngToCsvPath(pngPath);
        std::ofstream out(csv.c_str());
        if (!out){
            std::cerr << "[Playground][ERROR] Cannot open CSV for writing: " << csv << "\n";
            return;
        }
        out << "variant,residual,E_lo,E_hi,E_ctr,mu,mu_err,sigma,sigma_err\n";
        out.setf(std::ios::fixed);
        out << std::setprecision(8);

        for (std::size_t s=0; s<variants.size(); ++s){
            const auto& E   = eCtrS[s];
            const auto& MU  = muS[s];
            const auto& DMU = dmuS[s];
            const auto& SG  = sgS[s];
            const auto& DSG = dsgS[s];

            const std::string vq = csvQuote(variants[s]);
            const std::string rq = csvQuote(residuals[s]);

            for (std::size_t i=0; i<E.size() && i<MU.size() && i<SG.size(); ++i){
                // Skip obviously empty placeholders (both errors 0 and values 0)
                if ((DMU[i]==0.0 && DSG[i]==0.0) && (MU[i]==0.0 && SG[i]==0.0))
                    continue;

                const auto pr = edgesFromCenter(E[i]);
                out << vq << "," << rq << ","
                    << pr.first  << "," << pr.second << ","
                    << E[i]      << ","
                    << MU[i]     << "," << DMU[i] << ","
                    << SG[i]     << "," << DSG[i] << "\n";
            }
        }
        out.close();
        std::cout << "[Playground] Wrote " << csv << "\n";
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

    // -------------------- NEW: CSV (Gaussian) reader + overlay helpers --------------------
    struct GaussSeries {
        std::vector<double> E, mu, dmu, sg, dsg;
    };
    using GaussDB = std::map<std::string, std::map<std::string, GaussSeries>>; // variant -> residual -> series

    auto unquote = [](std::string s)->std::string {
        if (s.size()>=2 && s.front()=='"' && s.back()=='"') {
            s = s.substr(1, s.size()-2);
            // unescape ""
            std::string t; t.reserve(s.size());
            for (size_t i=0;i<s.size();++i){
                if (s[i]=='"' && i+1<s.size() && s[i+1]=='"') { t.push_back('"'); ++i; }
                else t.push_back(s[i]);
            }
            return t;
        }
        return s;
    };

    auto splitCSVLine = [&](const std::string& line)->std::vector<std::string>{
        std::vector<std::string> out;
        std::string cur; bool inQ=false;
        for(char ch : line){
            if (ch=='"'){ inQ=!inQ; cur.push_back(ch); }
            else if (ch==',' && !inQ){ out.push_back(cur); cur.clear(); }
            else cur.push_back(ch);
        }
        out.push_back(cur);
        for (auto& f : out) f = unquote(f);
        return out;
    };

    auto parseGaussCSV = [&](const std::string& path, GaussDB& db){
        std::ifstream fin(path.c_str());
        if (!fin.good()){
            std::cerr << "[Playground][WARN] Cannot open Gaussian CSV: " << path << "\n";
            return;
        }
        std::string line;
        // header
        if (!std::getline(fin, line)) return;

        // to avoid duplicates across multiple CSVs
        std::set<std::tuple<std::string,std::string,double,double>> seen;
        while (std::getline(fin, line)){
            if (line.empty()) continue;
            auto tok = splitCSVLine(line);
            if (tok.size() < 9) continue;

            const std::string variant  = tok[0];
            const std::string residual = tok[1];
            const double Elo  = atof(tok[2].c_str());
            const double Ehi  = atof(tok[3].c_str());
            const double Ectr = atof(tok[4].c_str());
            const double mu   = atof(tok[5].c_str());
            const double dmu  = atof(tok[6].c_str());
            const double sg   = atof(tok[7].c_str());
            const double dsg  = atof(tok[8].c_str());

            auto key = std::make_tuple(variant, residual, Elo, Ehi);
            if (seen.count(key)) continue;  // dedupe
            seen.insert(key);

            GaussSeries& S = db[variant][residual];
            S.E.push_back(Ectr);
            S.mu.push_back(mu);  S.dmu.push_back(dmu);
            S.sg.push_back(sg);  S.dsg.push_back(dsg);
        }
    };

    auto safeName = [](std::string s)->std::string{
        for (char& c : s) if (c==' ' || c==',' || c=='/' || c=='(' || c==')') c = '_';
        return s;
    };

    auto makeBuiltInSeries = [&](const std::vector<TH1F*>& H,
                                 const std::vector<std::pair<double,double>>& edges,
                                 double lo, double hi)
                                 -> GaussSeries {
        GaussSeries S;
        for (std::size_t i=0;i<edges.size();++i){
            if (!H[i] || H[i]->Integral()==0) continue;
            double m,me,r,re;
            std::tie(m,me,r,re) = statsFromHistRange(H[i], lo, hi);
            S.E .push_back(0.5*(edges[i].first + edges[i].second));
            S.mu.push_back(m);  S.dmu.push_back(me);
            S.sg.push_back(r);  S.dsg.push_back(re);
        }
        return S;
    };

    auto drawMeanRMSvsGaussOne = [&](const std::string& variantLabel,
                                     const std::string& residual,           // "phi" or "eta"
                                     const std::vector<TH1F*>& H,           // built-in hist set for this residual
                                     const GaussDB& db,
                                     const std::string& outDirPNG,
                                     const std::vector<std::pair<double,double>>& edges,
                                     double lo, double hi){
        if (H.empty()) return;

        // built-in (Mean/RMS)
        GaussSeries A = makeBuiltInSeries(H, edges, lo, hi);

        // gauss overlay (from CSVs) — translate display labels to CSV keys
        auto csvVariantKey = [&](const std::string& v)->std::string {
            if (v == "PDC raw")        return "no corr, scratch";
            if (v == "PDC corrected")  return "b(E) corr, scratch";
            // clusterizer labels already match CSV
            return v; // "no corr, cluster", "CorrectPosition, cluster"
        };

        GaussSeries B;
        const std::string key = csvVariantKey(variantLabel);
        auto iv = db.find(key);
        if (iv != db.end()){
            auto ir = iv->second.find(residual);   // "phi" or "eta"
            if (ir != iv->second.end()) B = ir->second;
        }
        // Collapse any duplicate Gaussian rows at the same E (keep the last)
        {
            const long long SCALE = 1000000; // 1e-6 GeV tolerance
            std::unordered_map<long long,size_t> lastPos;
            for (size_t j=0; j<B.E.size(); ++j) lastPos[ llround(B.E[j]*SCALE) ] = j;
            GaussSeries Buniq;
            Buniq.E.reserve(lastPos.size());
            for (const auto& kv : lastPos){
                size_t j = kv.second;
                Buniq.E .push_back(B.E[j]);
                Buniq.mu.push_back(B.mu[j]);   Buniq.dmu.push_back(B.dmu[j]);
                Buniq.sg.push_back(B.sg[j]);   Buniq.dsg.push_back(B.dsg[j]);
            }
            // keep B in increasing E order for tidy plots
            std::vector<size_t> ord(Buniq.E.size());
            std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](size_t i,size_t k){ return Buniq.E[i] < Buniq.E[k]; });
            GaussSeries Bsorted;
            for (size_t idx : ord){
                Bsorted.E .push_back(Buniq.E[idx]);
                Bsorted.mu.push_back(Buniq.mu[idx]);  Bsorted.dmu.push_back(Buniq.dmu[idx]);
                Bsorted.sg.push_back(Buniq.sg[idx]);  Bsorted.dsg.push_back(Buniq.dsg[idx]);
            }
            B = Bsorted;
        }

        // Keep only bins that exist in BOTH A (Mean/RMS from this run) and B (CSV)
        {
            const long long SCALE = 1000000;
            std::unordered_map<long long,size_t> posB;
            for (size_t j=0; j<B.E.size(); ++j) posB[ llround(B.E[j]*SCALE) ] = j;

            GaussSeries Akeep, Bkeep;
            for (size_t i=0; i<A.E.size(); ++i){
                long long key = llround(A.E[i]*SCALE);
                auto it = posB.find(key);
                if (it == posB.end()) continue;

                size_t j = it->second;
                // matched A (Mean/RMS)
                Akeep.E .push_back(A.E[i]);
                Akeep.mu.push_back(A.mu[i]);   Akeep.dmu.push_back(A.dmu[i]);
                Akeep.sg.push_back(A.sg[i]);   Akeep.dsg.push_back(A.dsg[i]);
                // matched B (Gaussian from CSV)
                Bkeep.E .push_back(B.E[j]);
                Bkeep.mu.push_back(B.mu[j]);   Bkeep.dmu.push_back(B.dmu[j]);
                Bkeep.sg.push_back(B.sg[j]);   Bkeep.dsg.push_back(B.dsg[j]);
            }
            A = Akeep;
            B = Bkeep;
        }

        // If nothing matches, skip the overlay safely
        if (A.E.empty() || B.E.empty()){
            std::cerr << "[Playground][WARN] No matched bins for " << variantLabel
                      << " / " << residual << " – skipping overlay.\n";
            return;
        }

        const Color_t col = (residual=="phi") ? (kRed+1) : (kBlue+1);

        // x-range from matched bins
        double xAxisMin = 0.0;
        double lastCtr  = std::max(A.E.back(), B.E.back());
        double binW     = edges.empty()? 0.0 : 0.5*(edges.back().second - edges.back().first);
        double xAxisMax = lastCtr + binW;


        // y-ranges (μ)
        auto minMu = [](const GaussSeries& S)->double{
            double v=+1e30; for(size_t i=0;i<S.mu.size();++i) v = std::min(v, S.mu[i]-S.dmu[i]); return v;
        };
        auto maxMu = [](const GaussSeries& S)->double{
            double v=-1e30; for(size_t i=0;i<S.mu.size();++i) v = std::max(v, S.mu[i]+S.dmu[i]); return v;
        };
        double muLo = +1e30, muHi = -1e30;
        if (!A.mu.empty()){ muLo = std::min(muLo, minMu(A)); muHi = std::max(muHi, maxMu(A)); }
        if (!B.mu.empty()){ muLo = std::min(muLo, minMu(B)); muHi = std::max(muHi, maxMu(B)); }
        if (muLo>muHi){ muLo=-1, muHi=1; }
        double padMu = 0.25*(muHi-muLo);

        // y-range (σ)
        auto maxSg = [](const GaussSeries& S)->double{
            double v=0; for(double x:S.sg) v=std::max(v,x); return v;
        };
        double sgHi = std::max(maxSg(A), maxSg(B));
        if (sgHi<=0) sgHi=1.0;

        std::vector<double> exA(A.E.size(),0.0), exB(B.E.size(),0.0);

        const std::string cname = "cGauss_" + safeName(variantLabel) + "_" + residual;

        // If a canvas with this name exists from a prior call, close it so nothing persists.
        if (auto* old = dynamic_cast<TCanvas*>(gROOT->FindObject(cname.c_str()))) {
            old->Close();
            gSystem->ProcessEvents();
        }

        TCanvas c(cname.c_str(), "Mean/RMS vs Gauss", 900, 800);
        c.Clear();

        auto pads = makeEqualHeightTwinPads(c, cname, 0.15,0.05, 0.14,0.02, 0.04,0.16);
        TPad* pTop = pads.first;
        TPad* pBot = pads.second;


        // μ(E)
        pTop->cd();
        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
        frU.GetXaxis()->SetLabelSize(0.0);
        frU.GetXaxis()->SetTitleSize(0.0);
        frU.GetXaxis()->SetTickLength(0.0);
        frU.Draw("AXIS");
        TLine l0(xAxisMin,0.0,xAxisMax,0.0); l0.SetLineStyle(2); l0.Draw();

        TGraphErrors gMuA((int)A.E.size(), A.E.data(), A.mu.data(), exA.data(), A.dmu.data());
        gMuA.SetMarkerStyle(20); // circles (requirement)
        gMuA.SetMarkerColor(col); gMuA.SetLineColor(col); gMuA.SetMarkerSize(1.05);
        gMuA.Draw("P SAME");
        TGraphErrors gMuB((int)B.E.size(), B.E.data(), B.mu.data(), exB.data(), B.dmu.data());
        gMuB.SetMarkerStyle(24); // open circles for Gauss overlay
        gMuB.SetMarkerColor(col); gMuB.SetLineColor(col); gMuB.SetMarkerSize(1.05);
        gMuB.Draw("P SAME");

        TLegend legU(0.62,0.7,0.93,0.85);
        legU.SetBorderSize(0); legU.SetFillStyle(0); legU.SetTextSize(0.035);
        legU.AddEntry(&gMuA, "Mean/RMS", "p");
        legU.AddEntry(&gMuB, "Gaussian fit", "p");
        legU.Draw();

        // σ(E)
        pBot->cd();
        TH1F frL("frL",";E  [GeV];#sigma_{RMS}  [rad]",1,xAxisMin,xAxisMax);
        frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");
        TGraphErrors gSgA((int)A.E.size(), A.E.data(), A.sg.data(), exA.data(), A.dsg.data());
        gSgA.SetMarkerStyle(20);
        gSgA.SetMarkerColor(col); gSgA.SetLineColor(col); gSgA.SetMarkerSize(1.05);
        gSgA.Draw("P SAME");
        TGraphErrors gSgB((int)B.E.size(), B.E.data(), B.sg.data(), exB.data(), B.dsg.data());
        gSgB.SetMarkerStyle(24); // open circles
        gSgB.SetMarkerColor(col); gSgB.SetLineColor(col); gSgB.SetMarkerSize(1.05);
        gSgB.Draw("P SAME");

        // Header
        c.cd();
        TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.038);
        h.DrawLatex(0.50,0.985,
            Form("%s  -  #Delta%s  Mean/RMS vs E  with Gaussian overlay",
                 variantLabel.c_str(), (residual=="phi" ? "#phi" : "#eta")));

        const std::string outPng = outDirPNG + "/Delta" + (residual=="phi"?"Phi_":"Eta_")
                                 + safeName(variantLabel) + "_MeanRMS_vs_Gauss.png";
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    
    // ---------------------------------------------------------------------
    // Unified residual-panel drawer: single curve OR overlay
    // ---------------------------------------------------------------------

    // Usage:
    //   std::vector<ResidualSpec> specs = {
    //     {h1, "#Delta#phi", kRed+1, 20, 1},
    //     {h2, "#Delta#eta", kBlue+1, 20, 1}
    //   };
    //   drawResidualPanel(gPad, "Clusterizer Output, No Position Correction",
    //                     eLo, eHi, specs, xMin, xMax, 1.0);
    //
    // 'textScale' can be <1 for pads inside tables.
    // ---------------------------------------------------------------------
    struct ResidualSpec {
        TH1F*       h;         // histogram to draw (original; function clones internally)
        const char* label;     // legend label & what appears in stats/nEvents text
        Color_t     col;       // color for markers/fit curve/text
        Style_t     mk;        // marker style (e.g. 20)
        int         ls;        // line style for Gaussian curve (1=solid, 2=dashed, …)
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

//        // Upper-right middle: stats for each curve (multi-line, always black text)
//        TLatex st; st.SetNDC(); st.SetTextFont(42);
//        st.SetTextAlign(33); st.SetTextSize(0.025 * textScale);
//        st.SetTextColor(kBlack);
//        double yStats = 0.78;
//        const double dy = 0.028 * textScale;         // vertical step between lines
//        const double padBetweenCurves = 0.018 * textScale; // extra gap after each curve block
//
//        for (std::size_t i=0; i<specs.size(); ++i) {
//            if (!specs[i].h || specs[i].h->Integral()==0) continue;
//            double m, me, r, re;
//            std::tie(m, me, r, re) = statsFromHistRange(specs[i].h, xMin, xMax);
//
//            // Line 1: label + Mean ± err
//            st.DrawLatex(0.92, yStats, Form("%s:  #mu = %.6g #pm %.2g", specs[i].label, m, me));
//            yStats -= dy;
//
//            // Line 2:        RMS ± err
//            st.DrawLatex(0.92, yStats, Form("       #sigma_{RMS} = %.6g #pm %.2g", r, re));
//            yStats -= (dy + padBetweenCurves);
//        }
//        // Middle-left: nEvents per curve (always black text)
//        TLatex nev; nev.SetNDC(); nev.SetTextFont(42);
//        nev.SetTextAlign(13); nev.SetTextSize(0.025 * textScale);
//        nev.SetTextColor(kBlack);
//        double yNev = 0.58;
//        for (const auto& s : specs) {
//            const double N = s.h ? s.h->GetEntries() : 0.0;
//            nev.DrawLatex(0.16, yNev, Form("nEvents %s = %.0f", s.label, N));
//            yNev -= (0.05 * textScale);
//        }

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

        const std::array<const char*,5> phiPat = {
            "h_phi_diff_raw_%s",
            "h_phi_diff_corr_%s",
            "h_phi_diff_cpRaw_%s",
            "h_phi_diff_cpCorr_%s",
            "h_phi_diff_cpBcorr_%s"
        };
        const std::array<const char*,5> phiLab = {
            "PDC raw",
            "PDC corrected",
            "no corr, cluster",
            "CorrectPosition, cluster",
            "b(E) corr, cluster"
        };

        std::array<std::vector<TH1F*>,5> PHI{};
        for (int v=0; v<5; ++v) PHI[v] = loadSet(phiPat[v]);

        std::vector<double> eCtr; eCtr.reserve(eEdges.size());
        for (std::size_t iE=0; iE<eEdges.size(); ++iE)
            eCtr.push_back(0.5*(eEdges[iE].first + eEdges[iE].second));

        std::array<std::vector<double>,5> MU, dMU, SG, dSG;
        std::vector<int> shown; shown.reserve(5);

        for (int v=0; v<5; ++v) {
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
                g->SetMarkerStyle(kMk[v]);
                g->SetMarkerColor(kCol[v]);
                g->SetLineColor  (kCol[v]);
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
                g->SetMarkerStyle(kMk[v]);
                g->SetMarkerColor(kCol[v]);
                g->SetLineColor  (kCol[v]);
                g->SetMarkerSize(1.0);
                g->Draw("P SAME");
                gSg.emplace_back(std::move(g));
            }

            c.cd();
            TLatex h; h.SetNDC(); h.SetTextAlign(22); h.SetTextSize(0.040);
            h.DrawLatex(0.50,0.985,"#Delta#phi  –  Mean / RMS vs E  (variants)");
            const std::string outP = dirPhiLog + "/DeltaPhi_MeanSigmaVsE_5Variants.png";
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

        // μ vs ln(E) fits (all available variants)
        if (eCtr.size() >= 2 && !shown.empty()) {
            std::vector<double> lnE(eCtr.size());
            std::transform(eCtr.begin(), eCtr.end(), lnE.begin(),
                           [](double e){ return std::log(e); });

            const double xLo = *std::min_element(lnE.begin(),lnE.end()) - 0.05;
            const double xHi = *std::max_element(lnE.begin(),lnE.end()) + 0.05;

            double yLo=1e30,yHi=-1e30;
            for (int v : shown) {
                yLo = std::min(yLo, *std::min_element(MU[v].begin(), MU[v].end()));
                yHi = std::max(yHi, *std::max_element(MU[v].begin(), MU[v].end()));
            }
            const double pad=0.15*(yHi-yLo); yLo-=pad; yHi+=0.35*(yHi-yLo);

            TCanvas c("cLnAll","μ vs lnE (Δφ variants)",900,640);
            c.SetLeftMargin(0.15); c.SetRightMargin(0.06);

            TH1F fr("",";ln E  [GeV];#mu  [rad]",1,xLo,xHi);
            fr.SetMinimum(yLo); fr.SetMaximum(yHi); fr.Draw("AXIS");

            std::ofstream fout(dirPhiLog + "/DeltaPhi_MuVsLogE_fit.txt");
            fout << "# variant                         p0(rad)          p1(rad)\n";

            std::vector<double> ex(eCtr.size(), 0.0);
            TLegend lg(0.16,0.75,0.90,0.90);
            lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.028);

            for (int v : shown) {
                auto* g = new TGraphErrors((int)eCtr.size(), lnE.data(), MU[v].data(),
                                           ex.data(), dMU[v].data());
                g->SetMarkerStyle(kMk[v]);
                g->SetMarkerColor(kCol[v]);
                g->SetLineColor  (kCol[v]);
                g->SetMarkerSize(1.1);
                g->Draw("P SAME");

                TF1 f("f","pol1",xLo,xHi);
                f.SetLineColor (kCol[v]);
                f.SetLineStyle (2);
                f.SetLineWidth (2);
                g->Fit(&f,"Q");
                f.Draw("SAME");

                fout << Form("%-28s % .6e   % .6e\n",
                             phiLab[v], f.GetParameter(0), f.GetParameter(1));

                auto* mk = new TMarker(0, 0, kMk[v]);
                mk->SetMarkerColor(kCol[v]);
                lg.AddEntry(mk,
                            Form("%s  (a=%.2e, b=%.2e)", phiLab[v],
                                 f.GetParameter(0), f.GetParameter(1)),
                            "p");
            }
            lg.Draw();

            c.SaveAs((dirPhiLog + "/DeltaPhi_MuVsLogE_All.png").c_str());
            c.Close();
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
    // ---------------------------------------------------------------------
    // Mean/RMS vs Gaussian-fit overlays per variant (8 PNGs total)
    //       Output folder: mean_sigmaGauss_HistCompare
    //       Gaussian sources:
    //        - Clusterizer:   .../Clusterizer_RawVsCorr_EtaPhi/DeltaEtaPhi_cpRaw_vs_cpCorr_MeanSigmaVsE_Overlay_4Series.csv
    //        - From-scratch:  .../fromScratch_RawVsCorr_EtaPhi/DeltaEtaPhi_raw_vs_corr_MeanSigmaVsE_Overlay_4Series.csv
    // ---------------------------------------------------------------------
    {
        const std::string gaussBase = "/Users/patsfan753/Desktop/scratchPDC_fitGuassianOutputAfterPhiTiltCorr";
        const std::string fClus  = gaussBase + "/Clusterizer_RawVsCorr_EtaPhi/DeltaEtaPhi_cpRaw_vs_cpCorr_MeanSigmaVsE_Overlay_4Series.csv";
        const std::string fScratch = gaussBase + "/fromScratch_RawVsCorr_EtaPhi/DeltaEtaPhi_raw_vs_corr_MeanSigmaVsE_Overlay_4Series.csv";

        GaussDB gdb;
        parseGaussCSV(fClus,   gdb);
        parseGaussCSV(fScratch,gdb);

        // Variant -> (phi-hists, eta-hists)
        struct Vspec { const char* label; const std::vector<TH1F*>* Hphi; const std::vector<TH1F*>* Heta; };
        std::vector<Vspec> V = {
            { "PDC raw",               &HphiScratchRaw, &HetaScratchRaw },
            { "PDC corrected",         &HphiScratchCorr,&HetaScratchCorr},
            { "no corr, cluster",      &Hphi,           &Heta           },
            { "CorrectPosition, cluster",&HphiCorr,     &HetaCorr       }
        };

        for (const auto& v : V){
            if (v.Hphi && !v.Hphi->empty()) drawMeanRMSvsGaussOne(v.label, "phi", *v.Hphi, gdb, dirMeanSigmaGaussCompare, eEdges, xMin, xMax);
            if (v.Heta && !v.Heta->empty()) drawMeanRMSvsGaussOne(v.label, "eta", *v.Heta, gdb, dirMeanSigmaGaussCompare, eEdges, xMin, xMax);
        }

        auto makeSingleResidualTable = [&](const std::vector<TH1F*>& H,
                                          const char* sym, Color_t col,
                                          const std::string& header,
                                          const std::string& outPng)
        {
            const int nCol=4, nRow=2, maxPads=nCol*nRow;
            TCanvas c(Form("cIndiv_%s", outPng.c_str()), "per-variant residual table", 1600, 900);
            c.SetTopMargin(0.10);
            c.Divide(nCol, nRow, 0, 0);

            int pad=1;
            for (std::size_t iE=0; iE<eEdges.size() && pad<=maxPads; ++iE)
            {
                if (!H[iE] || H[iE]->Integral()==0) { ++pad; continue; }
                c.cd(pad++);
                setPadMargins();

                std::vector<ResidualSpec> specs = {
                    { H[iE], sym, col, (Style_t)20, 1 }
                };

                // smaller text on table pads for neatness
                drawResidualPanel(gPad, header,
                                  eEdges[iE].first, eEdges[iE].second,
                                  specs, xMin, xMax, 0.85);
            }
            c.cd(0);
            drawHeaderNDC(header.c_str(), 0.975, 0.045);
            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // PDC raw
        makeSingleResidualTable(HphiScratchRaw, "#Delta#phi", kRed+1,
            "From scratch no correction  -  #Delta#phi  (All energy bins)",
            dirIndivPerVar + "/DeltaPhi_PDC_raw_AllBins_Table.png");
        makeSingleResidualTable(HetaScratchRaw, "#Delta#eta", kBlue+1,
            "From scratch no correction  -  #Delta#eta  (All energy bins)",
            dirIndivPerVar + "/DeltaEta_PDC_raw_AllBins_Table.png");

        // PDC corrected
        makeSingleResidualTable(HphiScratchCorr, "#Delta#phi", kRed+1,
            "From scratch with correction  -  #Delta#phi  (All energy bins)",
            dirIndivPerVar + "/DeltaPhi_PDC_corrected_AllBins_Table.png");
        makeSingleResidualTable(HetaScratchCorr, "#Delta#eta", kBlue+1,
            "From scratch with correction  -  #Delta#eta  (All energy bins)",
            dirIndivPerVar + "/DeltaEta_PDC_corrected_AllBins_Table.png");

        // Clusterizer no corr
        makeSingleResidualTable(Hphi, "#Delta#phi", kRed+1,
            "Clusterizer without Correction  -  #Delta#phi  (All energy bins)",
            dirIndivPerVar + "/DeltaPhi_no_corr__cluster_AllBins_Table.png");
        makeSingleResidualTable(Heta, "#Delta#eta", kBlue+1,
            "Clusterizer without Correction  -  #Delta#eta  (All energy bins)",
            dirIndivPerVar + "/DeltaEta_no_corr__cluster_AllBins_Table.png");

        // Clusterizer corrected
        makeSingleResidualTable(HphiCorr, "#Delta#phi", kRed+1,
            "Clusterizer w/ Correction  -  #Delta#phi  (All energy bins)",
            dirIndivPerVar + "/DeltaPhi_CorrectPosition__cluster_AllBins_Table.png");
        makeSingleResidualTable(HetaCorr, "#Delta#eta", kBlue+1,
            "Clusterizer w/ Correction  -  #Delta#eta  (All energy bins)",
            dirIndivPerVar + "/DeltaEta_CorrectPosition__cluster_AllBins_Table.png");

        // -----------------------------------------------------------------
        // NEW: Single-bin (6–8 GeV) overlays for clusterizer raw vs corrected
        //       • φ: red;  open circle = raw, closed circle = corrected
        //       • η: blue; open circle = raw, closed circle = corrected
        // -----------------------------------------------------------------
        auto findSliceIndex = [&](double Elo, double Ehi)->int {
            for (std::size_t i=0; i<eEdges.size(); ++i){
                if (std::fabs(eEdges[i].first - Elo) < 1e-6 &&
                    std::fabs(eEdges[i].second - Ehi) < 1e-6)
                    return static_cast<int>(i);
            }
            return -1;
        };

        auto makeSingleBinClusterizerOverlay = [&](int iE,
                                                   const std::vector<TH1F*>& Hraw,
                                                   const std::vector<TH1F*>& Hcorr,
                                                   const char* sym,
                                                   Color_t col,
                                                   const std::string& outPng)
        {
            if (iE < 0) return;
            if (iE >= static_cast<int>(Hraw.size()) || iE >= static_cast<int>(Hcorr.size())) return;
            if (!Hraw[iE] || !Hcorr[iE]) return;
            if (Hraw[iE]->Integral()<=0 || Hcorr[iE]->Integral()<=0) return;

            TCanvas c(Form("cOneBin_%s", outPng.c_str()), "single-bin overlay", 1000, 750);

            std::vector<ResidualSpec> specs = {
                // RAW (no correction) = closed circle (20)
                { Hraw[iE],  Form("%s without Correction", sym), col, (Style_t)20, 1 },
                // CORR (with correction) = open circle (24)
                { Hcorr[iE], Form("%s with Correction",     sym), col, (Style_t)24, 1 }
            };

            const double eLo = eEdges[iE].first;
            const double eHi = eEdges[iE].second;

            drawResidualPanel(&c,
                              "Clusterizer (single slice)",
                              eLo, eHi,
                              specs, xMin, xMax, 1.0);

            // --- draw vertical guides at x = ±0.025 rising ~20% of the y-range, with labels ---
            c.cd(); // ensure we are on the current canvas/pad that has the frame
            const double yMin = gPad->GetUymin();
            const double yMax = gPad->GetUymax();
            const double yTop = yMin + 0.20*(yMax - yMin);     // 20% of visible range
            const double yLbl = yMin + 0.22*(yMax - yMin);     // slightly above the line

            // guide at +0.025
            TLine vpos( 0.025, yMin, 0.025, yTop );
            vpos.SetLineColor(kGray+2);
            vpos.SetLineStyle(2);
            vpos.SetLineWidth(2);
            vpos.Draw();

            // guide at -0.025
            TLine vneg(-0.025, yMin, -0.025, yTop );
            vneg.SetLineColor(kGray+2);
            vneg.SetLineStyle(2);
            vneg.SetLineWidth(2);
            vneg.Draw();

            // label(s) above the lines
            TLatex lab;
            lab.SetTextAlign(21);         // center horizontally, bottom vertically
            lab.SetTextSize(0.030);       // user-coordinate text (scaled by pad)
            lab.SetTextColor(kGray+2);

            // You can label both, but labeling the +0.025 is usually enough:
            lab.DrawLatex( 0.025, yLbl, "one tower width ~ 0.025" );
            // If you also want the negative side labeled, uncomment the next line:
            // lab.DrawLatex(-0.025, yLbl, "one tower width ~ 0.025");

            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // -----------------------------------------------------------------
        // NEW: Single-bin (6–8 GeV) overlays for PDC (from-scratch) raw vs corrected
        //       • same rendering as clusterizer overlays, but PDC labels/header
        // -----------------------------------------------------------------
        auto makeSingleBinScratchOverlay = [&](int iE,
                                               const std::vector<TH1F*>& Hraw,
                                               const std::vector<TH1F*>& Hcorr,
                                               const char* sym,
                                               Color_t col,
                                               const std::string& outPng)
        {
            if (iE < 0) return;
            if (iE >= static_cast<int>(Hraw.size()) || iE >= static_cast<int>(Hcorr.size())) return;
            if (!Hraw[iE] || !Hcorr[iE]) return;
            if (Hraw[iE]->Integral()<=0 || Hcorr[iE]->Integral()<=0) return;

            TCanvas c(Form("cOneBin_%s", outPng.c_str()), "single-bin overlay (PDC)", 1000, 750);

            std::vector<ResidualSpec> specs = {
                // RAW (no corr) = closed circle (20)
                { Hraw[iE],  Form("%s from scratch (no corr)", sym), col, (Style_t)20, 1 },
                // CORR (with corr) = open circle (24)
                { Hcorr[iE], Form("%s from scratch (corr)",    sym), col, (Style_t)24, 1 }
            };

            const double eLo = eEdges[iE].first;
            const double eHi = eEdges[iE].second;

            drawResidualPanel(&c,
                              "From-scratch (single slice)",
                              eLo, eHi,
                              specs, xMin, xMax, 1.0);

            // vertical guides at ±0.025 with label (same as clusterizer overlay)
            c.cd();
            const double yMin = gPad->GetUymin();
            const double yMax = gPad->GetUymax();
            const double yTop = yMin + 0.20*(yMax - yMin);
            const double yLbl = yMin + 0.22*(yMax - yMin);

            TLine vpos( 0.025, yMin, 0.025, yTop ); vpos.SetLineColor(kGray+2); vpos.SetLineStyle(2); vpos.SetLineWidth(2); vpos.Draw();
            TLine vneg(-0.025, yMin, -0.025, yTop ); vneg.SetLineColor(kGray+2); vneg.SetLineStyle(2); vneg.SetLineWidth(2); vneg.Draw();

            TLatex lab; lab.SetTextAlign(21); lab.SetTextSize(0.030); lab.SetTextColor(kGray+2);
            lab.DrawLatex( 0.025, yLbl, "one tower width ~ 0.025" );

            c.SaveAs(outPng.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << outPng << "\n";
        };

        // 6–8 GeV slice index
        const int iE_6_8 = findSliceIndex(6.0, 8.0);

        // φ overlay (clusterizer raw vs corrected)
        makeSingleBinClusterizerOverlay(iE_6_8, Hphi, HphiCorr,
            "#Delta#phi", kRed+1,
            dirIndivPerVar + "/DeltaPhi_clusterizer_6to8GeV_RawVsCorr_Overlay.png");

        // η overlay (clusterizer raw vs corrected) — CLUS should be RED
        makeSingleBinClusterizerOverlay(iE_6_8, Heta, HetaCorr,
            "#Delta#eta", kRed+1,
            dirIndivPerVar + "/DeltaEta_clusterizer_6to8GeV_RawVsCorr_Overlay.png");

        // φ overlay (PDC raw vs corrected) — PDC should be BLUE
        makeSingleBinScratchOverlay(iE_6_8, HphiScratchRaw, HphiScratchCorr,
            "#Delta#phi", kBlue+1,
            dirIndivPerVar + "/DeltaPhi_fromScratch_6to8GeV_RawVsCorr_Overlay.png");

        // η overlay (PDC raw vs corrected)
        makeSingleBinScratchOverlay(iE_6_8, HetaScratchRaw, HetaScratchCorr,
            "#Delta#eta", kBlue+1,
            dirIndivPerVar + "/DeltaEta_fromScratch_6to8GeV_RawVsCorr_Overlay.png");
        }


    // -----------------------------------------------------------------
    // NEW: Save Δφ vs η 2D maps per energy bin & variant (only if present)
    // -----------------------------------------------------------------
    {
        // variant key in ROOT file  ->  pretty name for PNG
        struct V2D { const char* key; const char* label; };
        const std::vector<V2D> v2d = {
            {"h2_dphi_vs_eta_RAW_%s",     "Clusterizer_RAW"},
            {"h2_dphi_vs_eta_CP_%s",      "Clusterizer_CP"},
            {"h2_dphi_vs_eta_BCORR_%s",   "Clusterizer_BCORR"},
            {"h2_dphi_vs_eta_PDCraw_%s",  "PDC_RAW"},
            {"h2_dphi_vs_eta_PDCCorr_%s", "PDC_CORR"}
        };

        auto drawAndSave2D = [&](TH2F* h, const std::string& pngPath,
                                 double eLo, double eHi, const char* label)
        {
            if (!h || h->GetEntries() <= 0) return;

            const std::string cname = "c2D_" + safeName(pngPath);
            if (auto* old = dynamic_cast<TCanvas*>(gROOT->FindObject(cname.c_str()))) {
                old->Close();
                gSystem->ProcessEvents();
            }

            TCanvas c(cname.c_str(), "Δφ vs η", 900, 800);
            c.SetLeftMargin(0.12); c.SetRightMargin(0.14);
            c.SetTopMargin(0.08);  c.SetBottomMargin(0.12);

            gPad->SetRightMargin(0.14);
            h->SetContour(50);
            h->SetStats(0);
            h->Draw("COLZ");

            // Header & energy text
            TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextSize(0.04);
            head.DrawLatex(0.50, 0.96, Form("#Delta#phi vs #eta  –  %s", label));

            TLatex ebox; ebox.SetNDC(); ebox.SetTextAlign(33); ebox.SetTextSize(0.035);
            ebox.DrawLatex(0.96, 0.92, Form("[%.0f, %.0f) GeV", eLo, eHi));

            c.SaveAs(pngPath.c_str());
            c.Close();
            std::cout << "[Playground] Wrote " << pngPath << "\n";
        };

        for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
        {
            const double eLo = eEdges[iE].first;
            const double eHi = eEdges[iE].second;
            const TString tag = makeSliceTag(iE, eEdges[iE]); // matches booking tag

            for (const auto& spec : v2d)
            {
                TH2F* h2 = dynamic_cast<TH2F*>( fIn->Get(Form(spec.key, tag.Data())) );
                if (!h2) continue;

                const std::string outPng = dirEtaPhi2D + "/DeltaPhiVsEta_"
                                         + safeName(spec.label) + "_"
                                         + tag.Data() + ".png";
                drawAndSave2D(h2, outPng, eLo, eHi, spec.label);
            }
        }
    }

    // =====================================================================
    // NEW: 4WAY summary overlays (first non-empty lowest E bin) + 4-curve μ/σ(E)
    //       Outputs into: <outDir>/4WAYneededSUMMARIES
    // =====================================================================
    const std::string dir4Way = std::string(outDir) + "/4WAYneededSUMMARIES";
    ensureDir(dir4Way.c_str());

    // Find the first energy slice where BOTH hists in a pair are non-empty
    auto firstNonEmptyPair = [&](const std::vector<TH1F*>& A,
                                 const std::vector<TH1F*>& B)->int {
        for (std::size_t i=0; i<eEdges.size(); ++i) {
            if (i<A.size() && i<B.size() && A[i] && B[i] &&
                A[i]->Integral()>0 && B[i]->Integral()>0) return static_cast<int>(i);
        }
        return -1;
    };

    // Draw one two-series overlay for a given residual (phi/eta)
    auto overlayTwoSeriesFirstBin = [&](const char* baseLabel,   // "Variable-b (Benchmark)" or "Constant-b (Production)"
                                        Color_t color,           // kBlue+1 or kRed+1
                                        bool isPhi,              // true: phi, false: eta
                                        const std::vector<TH1F*>& Vraw,
                                        const std::vector<TH1F*>& Vcor,
                                        const std::string& outPngStem)
    {
        const int iE = firstNonEmptyPair(Vraw, Vcor);
        if (iE < 0) {
            std::cerr << "[Playground][WARN] 4WAY overlay skipped for " << baseLabel
                      << (isPhi? " (phi)" : " (eta)") << " – no non-empty first bin.\n";
            return;
        }

        TCanvas c(Form("c4Way_%s_%s", baseLabel, isPhi? "phi" : "eta"),
                  "4WAY first-bin overlay", 1000, 750);

        const char* symShort = isPhi ? "#phi" : "#eta";
        std::vector<ResidualSpec> specs = {
            // Uncorrected = filled circle (20)
            { Vraw[iE], Form("%s Uncorr, %s", baseLabel, symShort),  color, (Style_t)20, 1 },
            // Corrected   = open circle (24)
            { Vcor[iE], Form("%s Corr, %s",   baseLabel, symShort),  color, (Style_t)24, 1 }
        };

        // Dynamic header: "<Δφ or Δη> Uncorrected/Corrected Overlay for <Constant b-Correction | Variable b-Correction>"
        const char* residTitle = isPhi ? "#Delta#phi" : "#Delta#eta";
        const bool isConstB    = (std::string(baseLabel).find("Constant-b") != std::string::npos);
        const char* familyTitle= isConstB ? "Constant b-Correction" : "Variable b-Correction";
        std::string header     = Form("%s Uncorrected/Corrected Overlay for %s", residTitle, familyTitle);

        drawResidualPanel(&c,
                          header,
                          eEdges[iE].first, eEdges[iE].second,
                          specs, xMin, xMax, 1.0);

        const std::string outPng = dir4Way + "/" + outPngStem + ".png";
        c.SaveAs(outPng.c_str());
        c.Close();
        std::cout << "[Playground] Wrote " << outPng << "\n";
    };

    // ---- The 4 requested first-bin overlays ----
    // PDC (2×2 Block-Local), φ and η
    overlayTwoSeriesFirstBin("Variable-b (Benchmark)", kBlue+1,  true,
                             HphiScratchRaw,  HphiScratchCorr,
                             "PDC_phi_RawVsCorr_FirstBin");
    overlayTwoSeriesFirstBin("Variable-b (Benchmark)", kBlue+1,  false,
                             HetaScratchRaw,  HetaScratchCorr,
                             "PDC_eta_RawVsCorr_FirstBin");

    // CLUS (Constant-b (Production)), φ and η
    overlayTwoSeriesFirstBin("Constant-b (Production)",    kRed+1,   true,
                             Hphi,             HphiCorr,
                             "CLUS_phi_RawVsCorr_FirstBin");
    overlayTwoSeriesFirstBin("Constant-b (Production)",    kRed+1,   false,
                             Heta,             HetaCorr,
                             "CLUS_eta_RawVsCorr_FirstBin");

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

    auto drawMuSigmaFour = [&](bool isPhi,
                               const MuSigSeries& CLUSraw,
                               const MuSigSeries& CLUScor,
                               const MuSigSeries& PDCraw,
                               const MuSigSeries& PDCcor,
                               const std::string& outPng)
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
        if (muLo > muHi) { muLo = -1; muHi = 1; }
        const double padMu = 0.25*(muHi - muLo);

        TH1F frU("frU",";E  [GeV];#mu  [rad]",1,xMinE,xMaxE);
        frU.SetMinimum(muLo - padMu);
        frU.SetMaximum(muHi + padMu);
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

        // ---------- σ(E)
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

        // Legend in the μ(E) panel (upper-right), 2 columns
        pTop->cd();  // ensure the legend is drawn in the top pad
        TLegend lg(0.37, 0.7, 0.92, 0.86);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.041);
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

    // Build series and draw for φ
    MuSigSeries CLUSraw_phi = makeMuSigSeries(Hphi);
    MuSigSeries CLUScor_phi = makeMuSigSeries(HphiCorr);
    MuSigSeries PDCraw_phi  = makeMuSigSeries(HphiScratchRaw);
    MuSigSeries PDCcor_phi  = makeMuSigSeries(HphiScratchCorr);
    drawMuSigmaFour(true,  CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi,
                    dir4Way + "/FourVariants_MeanSigmaVsE_Phi.png");

    // Build series and draw for η
    MuSigSeries CLUSraw_eta = makeMuSigSeries(Heta);
    MuSigSeries CLUScor_eta = makeMuSigSeries(HetaCorr);
    MuSigSeries PDCraw_eta  = makeMuSigSeries(HetaScratchRaw);
    MuSigSeries PDCcor_eta  = makeMuSigSeries(HetaScratchCorr);
    drawMuSigmaFour(false, CLUSraw_eta, CLUScor_eta, PDCraw_eta, PDCcor_eta,
                    dir4Way + "/FourVariants_MeanSigmaVsE_Eta.png");

    // Parallel PNGs with μ(E) clamped to ±0.025 rad (top panel only)
    drawMuSigmaFour_muClamp(true,
                            CLUSraw_phi, CLUScor_phi, PDCraw_phi, PDCcor_phi,
                            dir4Way + "/FourVariants_MeanSigmaVsE_Phi_muClamp_pm0p025.png");

    drawMuSigmaFour_muClamp(false,
                            CLUSraw_eta, CLUScor_eta, PDCraw_eta, PDCcor_eta,
                            dir4Way + "/FourVariants_MeanSigmaVsE_Eta_muClamp_pm0p025.png");
    // ---------------------------------------------------------------------
    std::cout << "[Playground] Completed Δφ & Δη outputs into: " << outDir << "\n";
}




void PDCanalysis()
{
  // 1) Style / stat
  gROOT->LoadMacro("sPhenixStyle.C");
  SetsPhenixStyle();
    
  gStyle->SetOptStat(0);
    //PositionDep_sim_ALL_withUpdatedPhiTiltCorr_pinkValuesOnly
    //
    
  const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root";

  // 2) Open input
  std::cout << "[INFO] Opening file: " << inFile << "\n";
  TFile* fIn = TFile::Open(inFile, "READ");
  if(!fIn || fIn->IsZombie())
  {
    std::cerr << "[ERROR] Could not open file: " << inFile << std::endl;
    return;
  }
  std::cout << "[INFO] Successfully opened file: " << inFile << "\n";
    EBinningMode binMode = EBinningMode::kRange;

    if (fIn->Get("h3_blockCoord_E_range") == nullptr &&   // no “range” histo
        fIn->Get("h3_blockCoord_E_disc")  != nullptr)     // but “disc” exists
    {
        binMode = EBinningMode::kDiscrete;
        std::cout << "[PDCanalysis]  Detected DISCRETE energy binning\n";
    }
    else
    {
        std::cout << "[PDCanalysis]  Detected RANGE energy binning\n";
    }

    /* build the tag exactly like in the producer code */
    const char* modeTag = (binMode == EBinningMode::kRange ? "range" : "disc");
    
    // 6) Output directory + bValFile
    const char* outDir = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput";
    gSystem->mkdir(outDir, true);

    const char* outDirAll = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/allHistOutput";
    gSystem->mkdir(outDirAll, true);
    
    std::cout << std::left
              << std::setw(30) << "HISTOGRAM NAME"
              << std::setw(20) << "CLASS TYPE"
              << std::setw(15) << "ENTRIES"
              << std::endl;
    std::cout << std::string(65, '=') << std::endl;

  // We'll store the list of keys for reuse
  TList* keyList = fIn->GetListOfKeys();
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
    
//  std::cout << "\n[INFO] Now saving *every* histogram as a PNG => " << outDirAll << std::endl;
//
//  TIter nextAll(fIn->GetListOfKeys());
//  TKey* keyAll;
//  while( (keyAll = (TKey*) nextAll()) )
//  {
//      TClass* cl = gROOT->GetClass(keyAll->GetClassName());
//      if(!cl) continue;
//
//      TObject* obj = keyAll->ReadObj();
//      if(!obj) continue;
//
//      // We'll skip non-histogram objects
//      if(!cl->InheritsFrom("TH1") && !cl->InheritsFrom("TProfile")) continue;
//
//      TH1* htmp = dynamic_cast<TH1*>(obj);
//      if(!htmp) continue;
//
//      // Build output PNG name
//      TString histName = htmp->GetName();
//      TString outPNG = TString(outDirAll) + "/" + histName + ".png";
//
//      TCanvas ctemp("ctemp","",800,600);
//      ctemp.SetLeftMargin(0.12);
//      ctemp.SetRightMargin(0.15);
//      ctemp.SetBottomMargin(0.12);
//
//      // Decide how to draw
//      if( htmp->InheritsFrom("TH2") ) {
//        htmp->Draw("COLZ");
//      }
//      else if(htmp->InheritsFrom("TH3")) {
//        // 3D is tricky to visualize. We'll do "BOX".
//        htmp->Draw("BOX");
//      }
//      else if(htmp->InheritsFrom("TProfile")) {
//        htmp->Draw("E1");
//      }
//      else {
//        // presumably TH1
//        htmp->Draw("E");
//      }
//
//      ctemp.SaveAs(outPNG);
//      delete htmp;
//  }

    /* ------------------------------------------------------------------ */
    /* 2) fetch the 3-D histograms using the detected tag                 */
    /* ------------------------------------------------------------------ */
    const TString hNameUnc = Form("h3_blockCoord_E_%s",      modeTag);
    const TString hNameCor = Form("h3_blockCoord_Ecorr_%s",  modeTag);

    TH3F* hUnc3D = dynamic_cast<TH3F*>( fIn->Get(hNameUnc) );
    if (!hUnc3D)
    {
      std::cerr << "[ERROR] Histogram '" << hNameUnc << "' not found – abort.\n";
      fIn->Close();
      return;
    }
    hUnc3D->Sumw2();

    TH3F* hCor3D = dynamic_cast<TH3F*>( fIn->Get(hNameCor) );   // may be null
    if (hCor3D) hCor3D->Sumw2();
    

  // 5) Our 8 energy slices => [2..4), [4..6), etc.
  std::vector<std::pair<double,double>> eEdges = {
    {2.0,4.0},{4.0,6.0},{6.0,8.0},{8.0,10.0},
    {10.0,12.0},{12.0,15.0},{15.0,20.0},{20.0,30.0}
  };

    
  std::string bValFile = std::string(outDir) + "/bValues.txt";
  std::ofstream bOut(bValFile);
  if(!bOut.is_open())
  {
    std::cerr << "[WARN] Cannot open " << bValFile << " => won't save b-values.\n";
  }
  else
  {
    bOut << "# E range    best-b  (PHI or ETA)\n";
  }


  MakeDeltaPhiEtaPlayground(
        "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL.root",
        "/Users/patsfan753/Desktop/scratchPDC",
        -0.035, 0.035
  );
    
  Plot2DBlockEtaPhi(hUnc3D, hCor3D, isFirstPass, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/2DPlots");

  FitLocalPhiEta(hUnc3D,           // uncorrected 3-D histogram
                   hCor3D,           // corrected 3-D histogram (may be nullptr)
                   isFirstPass,      // same toggle you already use
                   eEdges,           // vector with the eight E-ranges
                    "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits",           // where PNGs will be written
                   bOut);            // SAME open ofstream instance

  OverlayUncorrPhiEta(hUnc3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
  PlotPhiShiftAndWidth(hUnc3D, hCor3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
    /* ----- build the input views ------------------------------------------------- */
    /* ------------------------------------------------------------ *
     *  Re‑create the five histogram vectors for every E‑slice
     *  by loading them directly from the ROOT file.  The Tag that
     *  the producer used is identical to the one in the Fit code:
     *        RANGE  →  "low_high"   (e.g.  2_4)
     *        DISC   →  "E<low>"     (e.g.  E2)
     * ------------------------------------------------------------ */

    auto makeSliceTag = [binMode](std::size_t iE,
                                  const std::pair<double,double>& edge)
    {
        return (binMode == EBinningMode::kRange)
               ? TString::Format("%.0f_%.0f", edge.first, edge.second)
               : TString::Format("E%.0f",     edge.first);
    };

    /* φ‑residual histograms ------------------------------------------------- */
    std::array<std::vector<TH1F*>,5> phiView;
    for (int v = 0; v < 5; ++v) phiView[v].resize(eEdges.size(), nullptr);

    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        const TString tag = makeSliceTag(iE, eEdges[iE]);

        phiView[0][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_raw_%s"    , tag.Data())) );
        phiView[1][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_corr_%s"   , tag.Data())) );
        phiView[2][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpRaw_%s"  , tag.Data())) );
        phiView[3][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpCorr_%s" , tag.Data())) );
        phiView[4][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpBcorr_%s", tag.Data())) );
    }

    /* η‑residual histograms  (global |z_vtx|-mixed) ------------------------- */
    std::array<std::vector<TH1F*>,5> etaGlobal;
    for (int v = 0; v < 5; ++v) etaGlobal[v].resize(eEdges.size(), nullptr);

    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        const TString tag = makeSliceTag(iE, eEdges[iE]);

        etaGlobal[0][iE] = static_cast<TH1F*>( fIn->Get(Form("h_eta_diff_raw_%s"    , tag.Data())) );
        etaGlobal[1][iE] = static_cast<TH1F*>( fIn->Get(Form("h_eta_diff_corr_%s"   , tag.Data())) );
        etaGlobal[2][iE] = static_cast<TH1F*>( fIn->Get(Form("h_eta_diff_cpRaw_%s"  , tag.Data())) );
        etaGlobal[3][iE] = static_cast<TH1F*>( fIn->Get(Form("h_eta_diff_cpCorr_%s" , tag.Data())) );
        etaGlobal[4][iE] = static_cast<TH1F*>( fIn->Get(Form("h_eta_diff_cpBcorr_%s", tag.Data())) );
    }

    /* φ‑residual histograms  (global |z_vtx|-mixed) ------------------------- */
    std::array<std::vector<TH1F*>,5> phiGlobalSlices;
    for (int v = 0; v < 5; ++v) phiGlobalSlices[v].resize(eEdges.size(), nullptr);

    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        const TString tag = makeSliceTag(iE, eEdges[iE]);

        phiGlobalSlices[0][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_raw_%s"    , tag.Data())) );
        phiGlobalSlices[1][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_corr_%s"   , tag.Data())) );
        phiGlobalSlices[2][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpRaw_%s"  , tag.Data())) );
        phiGlobalSlices[3][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpCorr_%s" , tag.Data())) );
        phiGlobalSlices[4][iE] = static_cast<TH1F*>( fIn->Get(Form("h_phi_diff_cpBcorr_%s", tag.Data())) );
    }


    /* legacy |z| and signed‑z containers: η */
    std::array<std::vector<std::vector<TH1F*>>,5> etaVzAbs;   // legacy |z|
    std::array<std::vector<std::vector<TH1F*>>,5> etaVzPos;   // z>0
    std::array<std::vector<std::vector<TH1F*>>,5> etaVzNeg;   // z<0

    /* legacy |z| and signed‑z containers: φ */
    std::array<std::vector<std::vector<TH1F*>>,5> phiVzAbs;   // legacy |z|
    std::array<std::vector<std::vector<TH1F*>>,5> phiVzPos;   // z>0
    std::array<std::vector<std::vector<TH1F*>>,5> phiVzNeg;   // z<0

    for (int v = 0; v < 5; ++v)
    {
        etaVzAbs[v].resize(eEdges.size());
        etaVzPos[v].resize(eEdges.size());
        etaVzNeg[v].resize(eEdges.size());

        phiVzAbs[v].resize(eEdges.size());
        phiVzPos[v].resize(eEdges.size());
        phiVzNeg[v].resize(eEdges.size());

        for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
        {
            etaVzAbs[v][iE].resize(N_VzBins,  nullptr);  // absolute |z|
            etaVzPos[v][iE].resize(N_VzBins,  nullptr);  // z > 0
            etaVzNeg[v][iE].resize(N_VzBins,  nullptr);  // z < 0

            phiVzAbs[v][iE].resize(N_VzBins,  nullptr);  // absolute |z|
            phiVzPos[v][iE].resize(N_VzBins,  nullptr);  // z > 0
            phiVzNeg[v][iE].resize(N_VzBins,  nullptr);  // z < 0
        }
    }

    /* –– fill (η & φ) –– */
    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        const TString eTag = makeSliceTag(iE, eEdges[iE]);   // e.g. "2_4"

        for (int iVz = 0; iVz < N_VzBins; ++iVz)
        {
            const TString vzTag = Form("vz%d_%d",
                                       int(vzEdge[iVz]), int(vzEdge[iVz+1]));
            /* η */
            etaVzAbs[0][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_raw_%s_%s"    , eTag.Data(), vzTag.Data())) );
            etaVzAbs[1][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_corr_%s_%s"   , eTag.Data(), vzTag.Data())) );
            etaVzAbs[2][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpRaw_%s_%s"  , eTag.Data(), vzTag.Data())) );
            etaVzAbs[3][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpCorr_%s_%s" , eTag.Data(), vzTag.Data())) );
            etaVzAbs[4][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpBcorr_%s_%s", eTag.Data(), vzTag.Data())) );
            /* φ */
            phiVzAbs[0][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_phi_diff_raw_%s_%s"    , eTag.Data(), vzTag.Data())) );
            phiVzAbs[1][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_phi_diff_corr_%s_%s"   , eTag.Data(), vzTag.Data())) );
            phiVzAbs[2][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_phi_diff_cpRaw_%s_%s"  , eTag.Data(), vzTag.Data())) );
            phiVzAbs[3][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_phi_diff_cpCorr_%s_%s" , eTag.Data(), vzTag.Data())) );
            phiVzAbs[4][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_phi_diff_cpBcorr_%s_%s", eTag.Data(), vzTag.Data())) );
        }

        for (int iVz = 0; iVz < N_VzBins; ++iVz)
        {
            const TString vzTagP = Form("%s_vzP%d_%d", eTag.Data(),
                                        int(vzEdge[iVz]), int(vzEdge[iVz+1]));
            /* η */
            etaVzPos[0][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_raw_%s"   , vzTagP.Data()) ));
            etaVzPos[1][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_corr_%s"  , vzTagP.Data()) ));
            etaVzPos[2][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpRaw_%s", vzTagP.Data()) ));
            etaVzPos[3][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpCorr_%s",vzTagP.Data()) ));
            etaVzPos[4][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpBcorr_%s",vzTagP.Data()) ));
            /* φ */
            phiVzPos[0][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_raw_%s"   , vzTagP.Data()) ));
            phiVzPos[1][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_corr_%s"  , vzTagP.Data()) ));
            phiVzPos[2][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpRaw_%s", vzTagP.Data()) ));
            phiVzPos[3][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpCorr_%s",vzTagP.Data()) ));
            phiVzPos[4][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpBcorr_%s",vzTagP.Data()) ));
        }

        for (int iVz = 0; iVz < N_VzBins; ++iVz)
        {
            const TString vzTagN = Form("%s_vzN%d_%d", eTag.Data(),
                                        int(vzEdge[iVz]), int(vzEdge[iVz+1]));
            /* η */
            etaVzNeg[0][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_raw_%s"   , vzTagN.Data()) ));
            etaVzNeg[1][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_corr_%s"  , vzTagN.Data()) ));
            etaVzNeg[2][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpRaw_%s", vzTagN.Data()) ));
            etaVzNeg[3][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpCorr_%s", vzTagN.Data()) ));
            etaVzNeg[4][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_eta_diff_cpBcorr_%s", vzTagN.Data()) ));
            /* φ */
            phiVzNeg[0][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_raw_%s"   , vzTagN.Data()) ));
            phiVzNeg[1][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_corr_%s"  , vzTagN.Data()) ));
            phiVzNeg[2][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpRaw_%s", vzTagN.Data()) ));
            phiVzNeg[3][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpCorr_%s", vzTagN.Data()) ));
            phiVzNeg[4][iE][iVz] = static_cast<TH1F*>( fIn->Get( Form("h_phi_diff_cpBcorr_%s", vzTagN.Data()) ));
        }
    }

  /* ----- produce all plots ----------------------------------------------------- */
  const char* OUTROOT = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits";

  /* global φ (unchanged) */
  MakeDeltaPhiPlots(eEdges, binMode, phiView, OUTROOT);

  /* legacy absolute-|z| (η) */
  MakeDeltaEtaPlots(eEdges, binMode, etaGlobal, etaVzAbs, OUTROOT);

  /* η (positive/negative z) */
  const std::string etaPosDir = std::string(OUTROOT) + "/deltaEta/vertexDependent/positiveVertexDependent";
  const std::string etaNegDir = std::string(OUTROOT) + "/deltaEta/vertexDependent/negativeVertexDependent";
  ensureDir(etaPosDir);
  ensureDir(etaNegDir);
  makeEtaVertexTables(eEdges, etaVzPos, etaPosDir.c_str(), {0,1,2,3,4}, "#Delta#eta", 'P');
  makeEtaVertexTables(eEdges, etaVzNeg, etaNegDir.c_str(), {0,1,2,3,4}, "#Delta#eta", 'N');

  /* φ (positive/negative z) */
  const std::string phiPosDir = std::string(OUTROOT) + "/deltaPhi/vertexDependent/positiveVertexDependent";
  const std::string phiNegDir = std::string(OUTROOT) + "/deltaPhi/vertexDependent/negativeVertexDependent";
  ensureDir(phiPosDir);
  ensureDir(phiNegDir);
  makeEtaVertexTables(eEdges, phiVzPos, phiPosDir.c_str(), {0,1,2,3,4}, "#Delta#phi", 'P');
  makeEtaVertexTables(eEdges, phiVzNeg, phiNegDir.c_str(), {0,1,2,3,4}, "#Delta#phi", 'N');

    
  const char* out2DDir = "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/2DPlots";
  gSystem->mkdir(out2DDir, true);
  drawLego3D(hUnc3D, "unc", "UNCORRECTED",
               out2DDir, 0.045, 0.035);
  drawLego3D(hCor3D, "cor", "CORRECTED",
               out2DDir, 0.045, 0.035);

  // residual QA
  auditResidual(hUnc3D, "UNCORRECTED", out2DDir,
                  0.045, 0.035);
  auditResidual(hCor3D, "CORRECTED",  out2DDir,
                  0.045, 0.035);

//
//  makeLegoGifHD(hUnc3D, "unc", "UNCORRECTED",
//                  hCor3D, "cor", "CORRECTED",
//                  out2DDir);



  // after opening the file …
  TH1F* h_truth_vz = static_cast<TH1F*>( fIn->Get("h_truth_vz") );

  if(!h_truth_vz){
      std::cerr << "[ERROR] vertex-Z histograms not found!\n";
  } else {
      PlotVertexZTruthOnly(h_truth_vz, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/QA");
  }

    PlotBvaluesVsEnergy(hUnc3D, eEdges,
                        "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
//    // Retrieve the two maps
//    auto bRMSMap  = plotAshLogRMS_sideBySide(inFile);                // RMS-optimised Ash-b
//    auto bPhiMap  = PlotBvaluesVsEnergy(hUnc3D, eEdges,
//                        "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
//
//    // Overlay & export
//    PlotBcompare(bRMSMap, bPhiMap,
//                 "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput");
////
//    
//  std::string chi2Dir = std::string(outDir) + "/QA/Chi2_QA";
//  PlotChi2QA(inFile, chi2Dir.c_str());  // after opening the file …
//  TH1F* h_truth_vz = static_cast<TH1F*>( fIn->Get("h_truth_vz") );
//
//  if(!h_truth_vz){
//      std::cerr << "[ERROR] vertex-Z histograms not found!\n";
//  } else {
//      PlotVertexZTruthOnly(h_truth_vz, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/QA");
//  }
//
//
//    
//  std::string chi2Dir = std::string(outDir) + "/QA/Chi2_QA";
//  PlotChi2QA(inFile, chi2Dir.c_str());
  // (final) close text file if open
  if(bOut.is_open())
  {
    bOut.close();
    std::cout << "[INFO] Wrote b-values => " << bValFile << "\n"
              << "[INFO] (Check file for lines labeled PHI vs ETA.)\n";
  }

  // close input
  fIn->Close();
  delete fIn;

  std::cout << "[INFO] LocalPhiOnly2by2() completed successfully.\n\n";
}
