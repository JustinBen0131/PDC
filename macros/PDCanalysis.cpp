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
    lat.DrawLatex(0.58, 0.74, Form("#sigma=%.3f", bestSigma));

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
    gAshRef->GetYaxis()->SetTitle("#sigma_{x}  (local tower units)");

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
        Form("%.0f #leq E < %.0f GeV  (b=%.2f,  #sigma=%.3f)",
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
    gLogRef->GetYaxis()->SetTitle("#sigma_{x}  (local tower units)");

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
        Form("%.0f #leq E < %.0f GeV  (w_{0}=%.2f,  #sigma=%.3f)",
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
                                "#sigma_{x}  (local tower units)");
  makeSinglePanel(logRes, "Log", "w_{0}",
                                "#sigma_{x}  (local tower units)");
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
  gAsh.GetYaxis()->SetTitle("#sigma_{x}  (mm)");
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
               Form("#sigma_{x}=%.2f + %.2f E^{-1/2}  (mm)",
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
          pad->SetLeftMargin (0.12);
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


/******************************************************************************
 * OverlayDeltaPhiSlices  (v4 – rock‑solid fitter + extra diagnostics)
 *
 *  • keeps every file/PNG the old code wrote – paths are unchanged
 *  • adds      Δφ_RMSErrorVsE.png      and      Δφ_FracResVsE.png
 *
 * 2025‑06‑12  –  dramatically improved convergence & outlier resilience
 ******************************************************************************/
void OverlayDeltaPhiSlices(const std::vector<std::pair<double,double>>& eEdges,
                           EBinningMode  binMode,
                           bool          isFirstPass,
                           const char*   outDir)
{
  /* ------------------------------------------------------------------ *
   * 0)   House‑keeping
   * ------------------------------------------------------------------ */
  gSystem->mkdir(outDir, /*recursive=*/true);
  gStyle->SetOptStat(0);

  const int N_E = static_cast<int>(eEdges.size());
  if (!N_E) { std::cerr << "[Δφ] eEdges empty – abort\n"; return; }

  /* energy‑slice labels exactly like the producer -------------------- */
  auto makeLabel = [&](int i)->std::string
  {
      if (binMode == EBinningMode::kRange)
          return Form("%.0f_%.0f", eEdges[i].first, eEdges[i].second);
      return Form("E%.0f", eEdges[i].first);              // discrete‑bin tag
  };

  /* ------------------------------------------------------------------ *
   * 1)   Bullet‑proof 1‑D Gaussian fitter
   * ------------------------------------------------------------------ *
   *  – adaptive fit window (inter‑quartile ±2.5×IQR)
   *  – multi‑start (4 amplitude/width seeds)   – “S” flag → TFitResult
   *  – iterative 3σ clip (≤3 passes, abort if N<50)
   *  – fall‑back to robust RMS if every fit fails
   * ------------------------------------------------------------------ */
  struct GRes { double N, mu, muErr, sg, sgErr; };

  auto robustGauss = [](TH1* h)->GRes
  {
      if (!h || h->Integral() == 0) return {0,0,0,0,0};

      /* -------- initial robust estimates -------------------------------- */
      const double q[3] = {0.25, 0.50, 0.75};
      double quart[3]; h->GetQuantiles(3, quart, const_cast<double*>(q));
      const double med = quart[1];
      const double iqr = quart[2] - quart[0];

      double winLo = med - 2.5 * iqr;
      double winHi = med + 2.5 * iqr;
      if (winLo >= winHi) {          // fall‑back to ±2×RMS
          winLo = h->GetMean() - 2.0*h->GetRMS();
          winHi = h->GetMean() + 2.0*h->GetRMS();
      }

      TF1 fG("fG","gaus", winLo, winHi);

      /* -------- multi‑start seeds --------------------------------------- */
      const double A = h->GetMaximum();
      const double S = std::max(1e-4, 0.5 * h->GetRMS());
      std::vector<std::tuple<double,double,double>> seeds = {
          { A       , med         , S       },
          { 1.5*A   , med         , 1.5*S   },
          { 0.5*A   , med + 0.3*S , 0.7*S   },
          { 1.2*A   , med - 0.3*S , 1.2*S   }
      };

      int    nIter   = 0;
      double bestχ2  = 1e99;
      GRes   best    {h->GetEntries(), 0,0, 0,0};

      const char* fitOpt = "QNR0S";        // Q‑quiet, N‑no draw, R‑range, 0‑do not store in hist, S‑store
      while (++nIter <= 3)                 // ≤3 σ‑clipping cycles
      {
          for (auto s : seeds) {
              fG.SetParameters(std::get<0>(s), std::get<1>(s), std::get<2>(s));
              TFitResultPtr r = h->Fit(&fG, fitOpt);
              if (!r.Get() || !r->IsValid()) continue;
              if (r->Chi2() < bestχ2 && std::fabs(fG.GetParameter(2))>1e-6) {
                  bestχ2   = r->Chi2();
                  best.mu  = fG.GetParameter(1);
                  best.muErr = fG.GetParError(1);
                  best.sg  = std::fabs(fG.GetParameter(2));
                  best.sgErr= fG.GetParError(2);
              }
          }
          if (bestχ2 < 1e98) break;       // success

          /* tighten window and repeat ------------------------------------ */
          winLo +=  0.5 * fG.GetParameter(2);
          winHi -=  0.5 * fG.GetParameter(2);
          if (winLo >= winHi || h->Integral() < 50) break;
          fG.SetRange(winLo, winHi);
      }

      /* fall‑back: robust RMS ------------------------------------------- */
      if (bestχ2 >= 1e98) {
          const double rms = h->GetRMS();
          const double rmsErr = rms / std::sqrt(2. * (h->GetEntries()-1));
          best.mu = h->GetMean();
          best.muErr = rmsErr;
          best.sg = rms;
          best.sgErr = rmsErr;
      }
      return best;
  };

    /* ------------------------------------------------------------------ *
     * 2)   Master canvas for per‑slice overlays
     *       – wider pads (2000 px total)
     *       – Y‑axis title pushed further left
     *       – smaller X‑axis tick labels
     *       – bold, colour‑coded “RAW / CORR” text, shifted right
     * ------------------------------------------------------------------ */
    const int nCols = 4;
    const int nRows = (binMode == EBinningMode::kDiscrete ? 4 : 2);

    /* • make the canvas wider (2000 px instead of 1600 px) */
    TCanvas c4x2(Form("DeltaPhiOverlay_%dx%d", nCols, nRows),
                 "#Delta#phi raw vs corrected (normalised)",
                 2000,                       /*  ⇦  width  */
                 600 * nRows);               /*      height */
    c4x2.Divide(nCols, nRows);

    /* containers for the later summary plots ---------------------------------- */
    std::vector<double> eCtr,
                        muRaw,  muRawErr,  sgRaw,  sgRawErr,
                        muCor,  muCorErr,  sgCor,  sgCorErr;

    /* keep objects alive until the very end ----------------------------------- */
    std::vector<TObject*> keep;

    /* ------------------------------------------------------------------ *
     * 3)   loop over energy slices
     * ------------------------------------------------------------------ */
    for (int iE = 0; iE < N_E; ++iE)
    {
        const double  eLo  = eEdges[iE].first;
        const double  eHi  = eEdges[iE].second;
        const TString tag  = (binMode==EBinningMode::kRange)
                           ? Form("%.0f_%.0f", eLo, eHi)
                           : Form("E%.0f",     eLo);

        /* ------------- fetch histograms ---------------------------------- */
        TH1F *hRaw  = static_cast<TH1F*>( gROOT->FindObject(
                            Form("h_phi_diff_raw_%s",  tag.Data()) ));
        TH1F *hCorr = (!isFirstPass)
                      ? static_cast<TH1F*>( gROOT->FindObject(
                            Form("h_phi_diff_corr_%s", tag.Data()) ))
                      : nullptr;

        if (!hRaw) {
            std::cerr << "[Δφ] missing RAW hist for tag " << tag << '\n';
            continue;
        }

        /* ------------- clone & re‑normalise ------------------------------ */
        hRaw  = static_cast<TH1F*>( hRaw ->Clone(Form("hRaw_%d", iE)) );
        hRaw ->SetDirectory(nullptr);   keep.push_back(hRaw);

        if (hCorr) {
            hCorr = static_cast<TH1F*>( hCorr->Clone(Form("hCor_%d", iE)) );
            hCorr->SetDirectory(nullptr); keep.push_back(hCorr);
        }

        if (hRaw ->Integral()>0) hRaw ->Scale(1./hRaw ->Integral());
        if (hCorr && hCorr->Integral()>0) hCorr->Scale(1./hCorr->Integral());

        /* ------------- robust Gaussian fits ------------------------------ */
        GRes R = robustGauss(hRaw);
        GRes C = (hCorr ? robustGauss(hCorr) : GRes{0,0,0,0,0});

        /* ------------- book‑keeping for the summary plots ---------------- */
        eCtr    .push_back( 0.5*(eLo+eHi) );
        muRaw   .push_back( R.mu );   muRawErr.push_back( R.muErr );
        sgRaw   .push_back( R.sg );   sgRawErr.push_back( R.sgErr );
        muCor   .push_back( C.mu );   muCorErr.push_back( C.muErr );
        sgCor   .push_back( C.sg );   sgCorErr.push_back( C.sgErr );

        /* ------------- draw on the pad ---------------------------------- */
        TPad* pad = static_cast<TPad*>( c4x2.cd(iE+1) );
        pad->SetLeftMargin  (0.25);        // ⇦ more space for Y‑title
        pad->SetBottomMargin(0.19);
        pad->SetTopMargin   (0.10);
        pad->SetRightMargin (0.06);

        /* raw = black ● ,  corr = red ● */
        hRaw ->SetMarkerStyle(20); hRaw ->SetMarkerColor(kBlack); hRaw ->SetLineColor(kBlack);
        if (hCorr) {
            hCorr->SetMarkerStyle(20); hCorr->SetMarkerColor(kRed); hCorr->SetLineColor(kRed);
        }

        /* axis cosmetics --------------------------------------------------- */
        const double yMax = std::max( hRaw->GetMaximum(),
                                      hCorr ? hCorr->GetMaximum() : 0.0 );
        hRaw->SetTitle   ( Form("#Delta#phi  [%.1f, %.1f)  GeV", eLo, eHi) );
        hRaw->GetXaxis()->SetTitle("#Delta#phi  (reco - truth)  [rad]");
        hRaw->GetYaxis()->SetTitle("Normalised counts");

        hRaw->GetYaxis()->SetTitleOffset( 2.15 );   // ⇦ shift Y‑title further left
        hRaw->GetXaxis()->SetLabelSize ( 0.028 );   // ⇦ smaller tick labels
        hRaw->GetYaxis()->SetRangeUser( 0, 1.30*yMax );

        /* draw ------------------------------------------------------------- */
        hRaw->Draw("E");
        if (hCorr) hCorr->Draw("E SAME");

        /* legend ----------------------------------------------------------- */
        TLegend* leg = new TLegend(0.28,0.72,0.54,0.88);   // (x1,y1,x2,y2) NDC
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42);  leg->SetTextSize(0.038);
        leg->AddEntry(hRaw ,  "RAW"       , "lp");
        if (hCorr) leg->AddEntry(hCorr,  "CORRECTED", "lp");
        leg->Draw();  keep.push_back(leg);

        /* bold, colour‑coded stats block ---------------------------------- */
        {
            /* RAW line – black, lower */
            TLatex txR;  txR.SetNDC();
            txR.SetTextFont(42);              // base font
            txR.SetTextSize(0.042);
            txR.SetTextColor(kBlack);
            txR.SetTextAlign(33);             // right‑aligned
            txR.DrawLatex(0.93, 0.84,
                Form("#bf{RAW:  #mu = %.4f,  #sigma = %.4f}", R.mu, R.sg));

            /* CORR line – red, above RAW */
            if (hCorr)
            {
                TLatex txC;  txC.SetNDC();
                txC.SetTextFont(42);
                txC.SetTextSize(0.042);
                txC.SetTextColor(kRed);
                txC.SetTextAlign(33);
                txC.DrawLatex(0.93, 0.90,
                    Form("#bf{CORR: #mu = %.4f, #sigma = %.4f}", C.mu, C.sg));
            }
        }
    } /* -------- end slice loop ---------------------------------------------- */


  /* ------------------------------------------------------------------ *
   * 4)   summary graphs (µ vs E, σ vs E) – unchanged appearance
   * ------------------------------------------------------------------ */
  const int nPts = static_cast<int>(eCtr.size());
  if (nPts == 0) return;

  auto makeG = [&](const std::vector<double>& y,
                   const std::vector<double>& yErr,
                   double dx, Color_t col, int mStyle)->TGraphErrors*
  {
      std::vector<double> x(nPts), ex(nPts,0.);
      for (int i=0;i<nPts;++i) x[i]=eCtr[i]+dx;
      auto g=new TGraphErrors(nPts, x.data(), y.data(), ex.data(), yErr.data());
      g->SetMarkerColor(col); g->SetLineColor(col);
      g->SetMarkerStyle(mStyle); g->SetMarkerSize(1.1);
      return g;
  };

  const double dx = 0.00;
  auto gMuRaw = makeG(muRaw, muRawErr, -dx, kBlack, 20);
  auto gMuCor = makeG(muCor, muCorErr, +dx, kRed  , 20);
  auto gSiRaw = makeG(sgRaw, sgRawErr, -dx, kBlack, 20);
  auto gSiCor = makeG(sgCor, sgCorErr, +dx, kRed  , 20);

  keep.insert(keep.end(),{gMuRaw,gMuCor,gSiRaw,gSiCor});

    /* -- summary canvas (µ top, σ bottom) --------------------------------- */
    {
        /* ----- axis range: pad one GeV on the left, go to the true upper edge ---- */
        const double xMin = eEdges.front().first  - 1.0;   // lower edge of first slice – 1 GeV
        const double xMax = eEdges.back ().second;         // upper edge of last slice  (30 GeV)

        TCanvas cSum("cMuSi","#Delta#phi mean / sigma vs E",860,760);
        TPad pT("pT","",0,0.34,1,1);  pT.Draw();
        TPad pB("pB","",0,0   ,1,0.34); pB.Draw();

        /* ---------- helper lambdas -------------------------------------- */
        auto fullMin = [](const std::vector<double>& v1,const std::vector<double>& e1,
                          const std::vector<double>& v2,const std::vector<double>& e2)->double
        {
            double m =  1e30;
            for (std::size_t i=0;i<v1.size();++i) m = std::min(m,v1[i]-(i<e1.size()?e1[i]:0.));
            for (std::size_t i=0;i<v2.size();++i) m = std::min(m,v2[i]-(i<e2.size()?e2[i]:0.));
            return m;
        };
        auto fullMax = [](const std::vector<double>& v1,const std::vector<double>& e1,
                          const std::vector<double>& v2,const std::vector<double>& e2)->double
        {
            double M = -1e30;
            for (std::size_t i=0;i<v1.size();++i) M = std::max(M,v1[i]+(i<e1.size()?e1[i]:0.));
            for (std::size_t i=0;i<v2.size();++i) M = std::max(M,v2[i]+(i<e2.size()?e2[i]:0.));
            return M;
        };

        /* ---------------- μ(E) panel ------------------------------------ */
        const double muMin = fullMin(muRaw,muRawErr,muCor,muCorErr);
        const double muMax = fullMax(muRaw,muRawErr,muCor,muCorErr);

        pT.cd();
        pT.SetLeftMargin(0.15);
        pT.SetBottomMargin(0.03);                       // hairline gap

        TH1F fMu("fMu","; ;#mu_{Gauss}  [rad]",1,xMin,xMax);
        fMu.SetMinimum(muMin - 0.10*std::fabs(muMin));
        fMu.SetMaximum(muMax + 0.10*std::fabs(muMax));

        /* hide the x-axis on the upper pad */
        TAxis* axT = fMu.GetXaxis();
        axT->SetLabelSize(0);
        axT->SetTickLength(0);
        axT->SetTitleSize(0);

        fMu.Draw("AXIS");
        gMuRaw->Draw("P SAME");
        gMuCor->Draw("P SAME");
        
        TLine* l0 = new TLine(xMin, 0.0, xMax, 0.0);
        l0->SetLineStyle(2);          // dashed
        l0->SetLineWidth(2);
        l0->SetLineColor(kGray+2);
        l0->Draw();

        /* keep the object until the function exits */
        keep.push_back(l0);           // <- change ‘guard’ to ‘keep’

        TLegend legMu(0.18,0.79,0.43,0.91);
        legMu.SetBorderSize(0); legMu.SetFillStyle(0); legMu.SetTextSize(0.035);
        legMu.AddEntry(gMuRaw,"RAW","lp");
        legMu.AddEntry(gMuCor,"CORRECTED","lp");
        legMu.Draw();

        /* ---------------- σ(E) panel ------------------------------------ */
        const double siMin = fullMin(sgRaw,sgRawErr,sgCor,sgCorErr);
        const double siMax = fullMax(sgRaw,sgRawErr,sgCor,sgCorErr);

        pB.cd();
        pB.SetLeftMargin(0.15);
        pB.SetTopMargin(0.06);
        pB.SetBottomMargin(0.38);                      // room for labels

        TH1F fSi("fSi",";E_{slice centre}  [GeV];#sigma_{Gauss}  [rad]",1,xMin,xMax);
        fSi.SetMinimum(std::max(0.0, siMin - 0.10*std::fabs(siMin)));
        fSi.SetMaximum(siMax + 0.10*std::fabs(siMax));

        TAxis* axB = fSi.GetXaxis();
        axB->SetNdivisions(505);          // 5 major, 5 minor
        axB->SetTickLength(0.05);
        axB->SetTitleOffset(1.1);

        fSi.Draw("AXIS");
        gSiRaw->Draw("P SAME");
        gSiCor->Draw("P SAME");

        /* ---------------------------------------------------------------- */
        cSum.SaveAs(TString(outDir)+"/DeltaPhi_MeanSigmaVsE.png");
    }

    /* ------------------------------------------------------------------ *
     * 5)   EXTRA diagnostics  (quadrature RMSE and fractional resolution)
     * ------------------------------------------------------------------ */
    {
        /* helper: returns a TGraphErrors with a horizontal offset “dxOff”   */
        auto gShift = [&](const std::vector<double>& y,
                          const std::vector<double>& yErr,
                          double dxOff, Color_t col) -> TGraphErrors*
        {
            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x(N), ex(N, 0.);
            for (int i = 0; i < N; ++i) x[i] = eCtr[i] + dxOff;

            auto *g = new TGraphErrors(N,
                                       x.data(),  y.data(),
                                       ex.data(), yErr.data());
            g->SetMarkerStyle(20);       // solid circle  ●
            g->SetMarkerSize (0.9);      // slightly smaller
            g->SetMarkerColor(col);
            g->SetLineColor  (col);
            keep.push_back(g);
            return g;
        };

        /* generic drawer used twice (RMSE & frac-res) -------------------- */
        auto makeDiag = [&](const char* cname,       // ROOT canvas / frame name
                            const char* yTit,        // y-axis label
                            const std::vector<double>& yR,
                            const std::vector<double>& yRE,
                            const std::vector<double>& yC,
                            const std::vector<double>& yCE,
                            const char* png)
        {
            const double dxOff = 0.12;               // RAW left, CORR right

            /* ---------- dynamic y-range : include ±error bars ------------------ */
            double yMin =  1e30,  yMax = -1e30;
            for (size_t i = 0; i < yR.size(); ++i){
                yMin = std::min({yMin, yR[i] - yRE[i], yC[i] - yCE[i]});
                yMax = std::max({yMax, yR[i] + yRE[i], yC[i] + yCE[i]});
            }
            /* make the bottom margin tighter – only 2 % below min, 5 % above max */
            const double span = yMax - yMin;
            yMin -= 0.02 * span;
            yMax += 0.05 * span;

            /* ---------- canvas & axis frame ----------------------------------- */
            TCanvas c(cname, cname, 860, 620);
            c.SetLeftMargin(0.15);  c.SetRightMargin(0.06);

            TH1F frame("f",
                       Form(";E_{slice}  [GeV];%s", yTit),
                       1, eCtr.front() - 1.0, eEdges.back().second);   // x-axis out to 30 GeV
            frame.SetMinimum(yMin);
            frame.SetMaximum(yMax);
            frame.Draw("AXIS");

            /* ---------- graphs ------------------------------------------------- */
            auto gRaw = gShift(yR, yRE, -dxOff, kBlack);
            auto gCor = gShift(yC, yCE, +dxOff, kRed  );
            gRaw->Draw("P SAME");  gCor->Draw("P SAME");

            /* ---------- legend placement & styling -------------------------------- */
            const bool isFracPlot = (std::strstr(png, "FracRes") != nullptr);

            /*  ── coordinates in NDC ───────────────────────────────────────────────
             *     RMSE  → bottom-left   : (0.15,0.12) – (0.45,0.32)
             *     Frac  → top-right     : (0.70,0.72) – (0.90,0.92)
             */
            double lx1 = isFracPlot ? 0.70 : 0.7;
            double ly1 = isFracPlot ? 0.72 : 0.4;
            double lx2 = isFracPlot ? 0.90 : 0.8;
            double ly2 = isFracPlot ? 0.92 : 0.52;

            TLegend leg(lx1, ly1, lx2, ly2);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextSize(0.026);
            leg.AddEntry(gRaw, "RAW",        "lp");
            leg.AddEntry(gCor, "CORRECTED",  "lp");
            leg.Draw();


            c.SaveAs(TString(outDir) + "/" + png);
        };

        /* (i) RMSE = √(μ²+σ²) ------------------------------------------- */
        std::vector<double> rmseR(nPts), rmseRE(nPts),
                            rmseC(nPts), rmseCE(nPts);
        for (int i = 0; i < nPts; ++i) {
            rmseR[i]  = std::hypot(muRaw[i], sgRaw[i]);
            rmseRE[i] = rmseR[i] * std::sqrt(
                           std::pow(muRawErr[i]/muRaw[i],2) +
                           std::pow(sgRawErr[i]/sgRaw[i],2));

            rmseC[i]  = std::hypot(muCor[i], sgCor[i]);
            rmseCE[i] = (rmseC[i] > 0.0)
                        ? rmseC[i] * std::sqrt(
                              std::pow(muCorErr[i]/muCor[i],2) +
                              std::pow(sgCorErr[i]/sgCor[i],2))
                        : 0.0;
        }
        makeDiag("cRMSE",
                 "#sqrt{#mu^{2} + #sigma^{2}}  [rad]",
                 rmseR, rmseRE, rmseC, rmseCE,
                 "DeltaPhi_RMSErrorVsE.png");

        /* (ii)  σ / |μ| --------------------------------------------------- */
        std::vector<double> fracR(nPts), fracRE(nPts),
                            fracC(nPts), fracCE(nPts);
        for (int i = 0; i < nPts; ++i) {
            if (std::fabs(muRaw[i]) > 1e-6) {
                fracR[i]  = sgRaw[i] / std::fabs(muRaw[i]);
                fracRE[i] = fracR[i] * std::sqrt(
                               std::pow(sgRawErr[i]/sgRaw[i],2) +
                               std::pow(muRawErr[i]/muRaw[i],2));
            }
            if (std::fabs(muCor[i]) > 1e-6) {
                fracC[i]  = sgCor[i] / std::fabs(muCor[i]);
                fracCE[i] = fracC[i] * std::sqrt(
                               std::pow(sgCorErr[i]/sgCor[i],2) +
                               std::pow(muCorErr[i]/muCor[i],2));
            }
        }
        makeDiag("cFrac",
                 "#sigma / |#mu|",
                 fracR, fracRE, fracC, fracCE,
                 "DeltaPhi_FracResVsE.png");
    }
    /* ------------------------------------------------------------------ */
  TString outAll = TString(outDir)+"/DeltaPhiCompare_AllOutput.png";
  c4x2.SaveAs(outAll);
  std::cout << "[Δφ] wrote summary → " << outAll << '\n';

  for (auto* o: keep) delete o;
}

void OverlayDeltaPhiClusterizerCP(const std::vector<std::pair<double,double>>& eEdges,
                                  EBinningMode  binMode,
                                  bool          isFirstPass,
                                  const char*   outDir)
{
  /* ---------- 0. boiler‑plate ---------- */
  gSystem->mkdir(outDir, /*recurse=*/true);
  gStyle->SetOptStat(0);

  const int N_E = static_cast<int>(eEdges.size());
  if (!N_E) { std::cerr << "[Δφ‑CP] eEdges empty – abort\n"; return; }

  auto tagOf = [&](int i)->std::string
  {
    return (binMode==EBinningMode::kRange)
         ? Form("%.0f_%.0f",eEdges[i].first,eEdges[i].second)
         : Form("E%.0f",     eEdges[i].first);
  };

  /* ---------- helper: robust Gaussian fit ---------- */
  struct GRes { double mu,muErr,sg,sgErr; };
  auto robust = [](TH1* h)->GRes
  {
    if(!h||h->Integral()==0) return {0,0,0,0};
    double q[3]={0.25,0.50,0.75}, quart[3]; h->GetQuantiles(3,quart,q);
    double med=quart[1],iqr=quart[2]-quart[0];
    double lo=med-2.5*iqr,hi=med+2.5*iqr;
    if(lo>=hi){ lo=h->GetMean()-2*h->GetRMS(); hi=h->GetMean()+2*h->GetRMS(); }
    TF1 f("g","gaus",lo,hi);  f.SetParameters(h->GetMaximum(),med,0.5*h->GetRMS());
    TFitResultPtr r=h->Fit(&f,"QNR0S");
    if(!r.Get()||!r->IsValid()) {            // fallback – RMS
       double rms=h->GetRMS(),err=rms/std::sqrt(2.*(h->GetEntries()-1));
       return {h->GetMean(),err,rms,err};
    }
    return {f.GetParameter(1),f.GetParError(1),
            std::fabs(f.GetParameter(2)),f.GetParError(2)};
  };

  /* ---------- 1. master canvas ---------- */
  const int nCols=4 , nRows=(binMode==EBinningMode::kDiscrete?4:2);
  TCanvas cMain("cΔφCP","Δφ comparison – clusteriser",1600,600*nRows);
  cMain.Divide(nCols,nRows);

  /* containers for summary graphs */
  std::vector<double> eCtr,
    muRaw,muCorr,muB , muRE,muCE,muBE,
    sgRaw,sgCorr,sgB , sgRE,sgCE,sgBE;

  std::vector<TObject*> guard;   // keep clones alive

  /* ---------- 2. per‑slice loop ---------- */
  for(int i=0;i<N_E;++i)
  {
    const double eLo=eEdges[i].first, eHi=eEdges[i].second;
    const std::string tag=tagOf(i);

    TH1F* hR = (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpRaw_%s"  ,tag.c_str()));
    TH1F* hC = (!isFirstPass)
               ? (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpCorr_%s",tag.c_str()))
               : nullptr;
    TH1F* hB = (!isFirstPass)
               ? (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpBcorr_%s",tag.c_str()))
               : nullptr;

    if(!hR) { std::cerr<<"[Δφ‑CP] missing RAW hist "<<tag<<"\n"; continue; }

    /* clone & normalise */
    auto clone = [&](TH1F* src,const char* nm,Color_t col,int mStyle)->TH1F*
    {
      if(!src) return nullptr;
      TH1F* c=(TH1F*)src->Clone(nm); c->SetDirectory(nullptr);
      if(c->Integral()>0) c->Scale(1./c->Integral());
      c->SetMarkerColor(col); c->SetLineColor(col);
      c->SetMarkerStyle(mStyle); c->SetMarkerSize(0.9);
      guard.push_back(c); return c;
    };

    TH1F* r = clone(hR,Form("r%d",i),kBlack,20);
    TH1F* c = clone(hC,Form("c%d",i),kRed  ,21);
    TH1F* b = clone(hB,Form("b%d",i),kBlue+1,22);

    /* robust fits */
    GRes R=robust(r), C=robust(c), B=robust(b);

    eCtr.push_back(0.5*(eLo+eHi));
    muRaw.push_back(R.mu); muRE.push_back(R.muErr);
    sgRaw.push_back(R.sg); sgRE.push_back(R.sgErr);

    muCorr.push_back(C.mu); muCE.push_back(C.muErr);
    sgCorr.push_back(C.sg); sgCE.push_back(C.sgErr);

    muB.push_back(B.mu);   muBE.push_back(B.muErr);
    sgB.push_back(B.sg);   sgBE.push_back(B.sgErr);

    /* draw overlay */
    TPad* pad=(TPad*)cMain.cd(i+1);
    pad->SetLeftMargin(0.22); pad->SetBottomMargin(0.18);

    double yMax=r->GetMaximum();
    if(c) yMax=std::max(yMax,c->GetMaximum());
    if(b) yMax=std::max(yMax,b->GetMaximum());

    r->SetTitle(Form("#Delta#phi  [%.1f, %.1f) GeV",eLo,eHi));
    r->GetXaxis()->SetTitle("#Delta#phi = #phi_{reco} - #phi_{truth}  [rad]");
    r->GetYaxis()->SetTitle("Probability density");
    r->GetYaxis()->SetRangeUser(0,1.25*yMax);

    r->Draw("E");
    if(c) c->Draw("E SAME");
    if(b) b->Draw("E SAME");

    TLegend* leg=new TLegend(0.24,0.70,0.57,0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    leg->AddEntry(r,"base, no CorrectPosition","lp");
    if(c) leg->AddEntry(c,"+ CorrectPosition"    ,"lp");
    if(b) leg->AddEntry(b,"+ updated b‑corr"     ,"lp");
    leg->Draw(); guard.push_back(leg);
  } /* slice loop */

  /* ---------- helper: make TGraphErrors ---------- */
  auto makeG=[&](const std::vector<double>& y,const std::vector<double>& ye,
                 double dx,Color_t col,int m)->TGraphErrors*
  {
    const int N=(int)eCtr.size();
    std::vector<double> x(N), ex(N,0.);
    for(int i=0;i<N;++i) x[i]=eCtr[i]+dx;
    auto g=new TGraphErrors(N,x.data(),y.data(),ex.data(),ye.data());
    g->SetMarkerColor(col); g->SetLineColor(col);
    g->SetMarkerStyle(m); g->SetMarkerSize(1.);
    guard.push_back(g); return g;
  };

  const double dx=0.0;
  auto gMuR = makeG(muRaw ,muRE,-dx,kBlack  ,20);
  auto gMuC = makeG(muCorr,muCE, 0,kRed    ,20);
  auto gMuB = makeG(muB   ,muBE,+dx,kBlue+1,20);

  auto gSiR = makeG(sgRaw ,sgRE,-dx,kBlack  ,20);
  auto gSiC = makeG(sgCorr,sgCE, 0,kRed    ,20);
  auto gSiB = makeG(sgB   ,sgBE,+dx,kBlue+1,20);

  /* ---------- 3. summary canvas (μ & σ) ---------- */
  {
    const double xMin=eCtr.front()-1;
    const double xMax = eEdges.back ().second;
    TCanvas cSum("cMuSigma","#Delta#phi mean / sigma vs E",860,780);
    TPad pT("pT","",0,0.37,1,1); pT.Draw();
    TPad pB("pB","",0,0  ,1,0.37); pB.Draw();

    /* μ(E) ------------------------------------------------ */
    pT.cd();
    pT.SetLeftMargin(0.15);
    pT.SetBottomMargin(0.03);          // hair-line gap to the lower pad

    // --- 1. compute inclusive [min,max] over all three curves, incl. ±Δμ ---
    double muMin =  1e30,  muMax = -1e30;
    auto upd = [&](double m, double e){          // helper to update extrema
          muMin = std::min(muMin, m - e);
          muMax = std::max(muMax, m + e);
    };
    for (size_t i=0; i<eCtr.size(); ++i){
          upd(muRaw [i], muRE[i]);
          upd(muCorr[i], muCE[i]);
          upd(muB   [i], muBE[i]);
    }

    // --- 2. add a symmetric 5 % margin so the end caps never touch the frame
    const double pad = 0.05 * (muMax - muMin);
    muMin -= pad;  muMax += pad;

    // --- 3. create the dummy frame with the new limits ---------------------
    TH1F fMu("fMu","; ;Gaussian mean  #mu  [rad]", 1, xMin, xMax);
    fMu.SetMinimum(muMin);
    fMu.SetMaximum(muMax);
    TAxis* axT = fMu.GetXaxis();       // hide x-labels on the top pad
    axT->SetLabelSize(0);  axT->SetTickLength(0);  axT->SetTitleSize(0);
    fMu.Draw("AXIS");

    // --- 4. overlay the three TGraphErrors ---------------------------------
    gMuR->Draw("P SAME");
    gMuC->Draw("P SAME");
    gMuB->Draw("P SAME");
      
      /* dashed reference line at μ = 0 ------------------------------------- */
      TLine* l0 = new TLine(xMin, 0.0, xMax, 0.0);   //  <-- create it!
      l0->SetLineStyle(2);          // kDashed
      l0->SetLineWidth(2);
      l0->SetLineColor(kGray+2);    // subtle grey
      l0->Draw();
      guard.push_back(l0);          // keep it alive until the end

    TLegend lg(0.38,0.1,0.82,0.35);
    lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.034);
    lg.AddEntry(gMuR,"Base Clusterizer, no BEmcRecCEMC::correctPosition","lp");
    lg.AddEntry(gMuC,"With BEmcRecCEMC::correctPosition","lp");
    lg.AddEntry(gMuB,"With Energy-Dep b-correction","lp"); lg.Draw();

    /* σ(E) ------------------------------------------------ */
    pB.cd(); pB.SetLeftMargin(0.15); pB.SetTopMargin(0.07);
    pB.SetBottomMargin(0.35);
    TH1F fSi("fSi",";Energy slice centre  [GeV];Gaussian width  #sigma  [rad]",
             1,xMin,xMax);
    fSi.SetMinimum(0);
    fSi.SetMaximum(1.3* *std::max_element(sgRaw.begin(),sgRaw.end()));
    fSi.Draw("AXIS");
    gSiR->Draw("P SAME"); gSiC->Draw("P SAME"); gSiB->Draw("P SAME");

    cSum.SaveAs(TString(outDir)+"/DeltaPhiCP_MeanSigmaVsE.png");
  }

    /* ---------- 4. diagnostics : RMSE  &  σ / |μ|  ----------------------- */
    auto diag = [&](const char* cname,            // ROOT canvas / pad title
                    const char* yTit,             // y-axis label
                    const std::vector<double>& yR, const std::vector<double>& yRE,
                    const std::vector<double>& yC, const std::vector<double>& yCE,
                    const std::vector<double>& yB, const std::vector<double>& yBE,
                    const char* filePNG)          // PNG filename
    {
        /* --- helper to build a graph with an x-offset --------------------- */
        auto gOf = [&](const std::vector<double>& y,
                       const std::vector<double>& ye,
                       double dx , Color_t col)->TGraphErrors*
        {
            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x(N), ex(N, 0.);
            for (int i = 0; i < N; ++i) x[i] = eCtr[i] + dx;
            auto *g = new TGraphErrors(N, x.data(), y.data(), ex.data(), ye.data());
            g->SetMarkerStyle(20); g->SetMarkerSize (0.8);
            g->SetMarkerColor(col); g->SetLineColor(col);
            guard.push_back(g); return g;
        };

        const double dxOff = 0.18;                     // black-red-blue spacing
        auto gR = gOf(yR, yRE, -dxOff, kBlack  );
        auto gC = gOf(yC, yCE,   0.0 , kRed    );
        auto gB = gOf(yB, yBE, +dxOff, kBlue+1);

        /* ---------- dynamic y-range incl. ±error -------------------------- */
        double yMin =  1e30, yMax = -1e30;
        for (size_t i = 0; i < yR.size(); ++i){
            yMin = std::min({yMin, yR[i]-yRE[i], yC[i]-yCE[i], yB[i]-yBE[i]});
            yMax = std::max({yMax, yR[i]+yRE[i], yC[i]+yCE[i], yB[i]+yBE[i]});
        }
        if (yMin > 0) yMin = 0;                       // keep zero on screen
        const double pad = 0.07 * (yMax - yMin);      // 7 % breathing room
        yMin -= pad;  yMax += pad;

        /* ---------- canvas & axis frame ----------------------------------- */
        TCanvas c(cname, cname, 860, 620);
        c.SetLeftMargin(0.15);  c.SetRightMargin(0.06);

        TH1F frame("f", TString::Format(";Energy slice centre  [GeV];%s", yTit).Data(),
                   1, eCtr.front() - 1.0,  eEdges.back().second); // full edge (to 30 GeV)
        frame.SetMinimum(yMin); frame.SetMaximum(yMax);
        frame.Draw("AXIS");

        /* ---------- draw the three curves --------------------------------- */
        gR->Draw("P SAME");  gC->Draw("P SAME");  gB->Draw("P SAME");

        /* ---------- legend placement & styling ---------------------------- */
        const bool isFrac = (std::strstr(filePNG, "FracRes") != nullptr);
        double lx1 = 0.15, lx2 = 0.60;
        double ly1 = isFrac ? 0.75 : 0.2;            // top-left vs bottom-left
        double ly2 = isFrac ? 0.92 : 0.4;

        TLegend lg(lx1, ly1, lx2, ly2);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.030);
        lg.AddEntry(gR, "Base Clusterizer, no BEmcRecCEMC::correctPosition", "lp");
        lg.AddEntry(gC, "With BEmcRecCEMC::correctPosition",                "lp");
        lg.AddEntry(gB, "With Energy-Dep b-correction",                     "lp");
        lg.Draw();

        c.SaveAs(TString(outDir) + "/" + filePNG);
    };

  /* RMSE */
  std::vector<double> rmR,rmRE, rmC,rmCE, rmB,rmBE;
  for(size_t i=0;i<eCtr.size();++i){
    auto rm=[&](double mu,double sg,double dmu,double dsg,
                double& err){
        double r=std::hypot(mu,sg);
        err=r*std::sqrt(std::pow(dmu/mu,2)+std::pow(dsg/sg,2));
        return r;
    };
    double e;
    rmR.push_back(rm(muRaw[i],sgRaw[i],muRE[i],sgRE[i],e)); rmRE.push_back(e);
    rmC.push_back(rm(muCorr[i],sgCorr[i],muCE[i],sgCE[i],e)); rmCE.push_back(e);
    rmB.push_back(rm(muB   [i],sgB   [i],muBE[i],sgBE[i],e)); rmBE.push_back(e);
  }
  diag("cRMSE","#sqrt{#mu^{2}+#sigma^{2}}  [rad]",
       rmR,rmRE, rmC,rmCE, rmB,rmBE,
       "DeltaPhiCP_RMSErrorVsE.png");

  /* σ/|μ| */
  std::vector<double> frR,frRE, frC,frCE, frB,frBE;
  for(size_t i=0;i<eCtr.size();++i){
    auto frac=[&](double mu,double sg,double dmu,double dsg,double& err){
        double f=sg/std::fabs(mu);
        err=f*std::sqrt(std::pow(dsg/sg,2)+std::pow(dmu/mu,2)); return f; };
    double e;
    frR.push_back(frac(muRaw[i],sgRaw[i],muRE[i],sgRE[i],e)); frRE.push_back(e);
    frC.push_back(frac(muCorr[i],sgCorr[i],muCE[i],sgCE[i],e)); frCE.push_back(e);
    frB.push_back(frac(muB   [i],sgB   [i],muBE[i],sgBE[i],e)); frBE.push_back(e);
  }
  diag("cFrac","#sigma / |#mu|",
       frR,frRE, frC,frCE, frB,frBE,
       "DeltaPhiCP_FracResVsE.png");

  /* ---------- 5. save overlay sheet ---------- */
  cMain.SaveAs(TString(outDir)+"/DeltaPhiCPCompare_AllOutput.png");
  std::cout<<"[Δφ‑CP] wrote overlays & summaries to "<<outDir<<'\n';

  for(auto* o:guard) delete o;
}


/******************************************************************************
 * OverlayDeltaPhiFiveWays  (2025‑06‑15)
 *
 *  – overlays 5 Δφ flavours per energy‑slice
 *  – writes   • per‑slice canvas
 *             • μ(E), σ(E)
 *             • RMSE(E) = √(μ²+σ²)
 *             • σ/|μ|(E)
 *  – axis titles and legend keys are explicit & unambiguous
 ******************************************************************************/
void OverlayDeltaPhiFiveWays(const std::vector<std::pair<double,double>>& eEdges,
                             EBinningMode  binMode,
                             bool          isFirstPass,
                             const char*   outDir)
{
  /* ---------- 0. boiler‑plate ---------- */
  gSystem->mkdir(outDir, /*recurse=*/true);
  gStyle->SetOptStat(0);

  static const Color_t  colArr[5]   = {kBlack, kRed, kBlue+1, kOrange+7, kMagenta+1};
  static const Style_t  mksArr[5]   = {20,     21,   22,      33,        29};
  static const char*    legTxtArr[5]= {
        "no correction, from scratch coord transforms",
        "energy dep b correction, from scratch coord transforms",
        "no correction, clusterizer coord transforms",
        "BEmcRecCEMC::CorrectPosition, clusterizer coord transforms",
        "energy dep b correction, clusterizer coord transforms"
  };
    
  const int N_E = static_cast<int>(eEdges.size());
  if (!N_E) { std::cerr << "[Δφ‑5] eEdges empty – abort\n"; return; }

  auto tagOf=[&](int i)->std::string
  {
    return (binMode==EBinningMode::kRange)
         ? Form("%.0f_%.0f",eEdges[i].first,eEdges[i].second)
         : Form("E%.0f",eEdges[i].first);
  };

  /* ---------- helper: robust Gaussian fit ---------- */
  struct GRes { double mu,muErr,sg,sgErr; };
  auto robust=[&](TH1* h)->GRes
  {
    if(!h||h->Integral()==0) return {0,0,0,0};
    const double q[3]={0.25,0.50,0.75};
    double quart[3]; h->GetQuantiles(3,quart,const_cast<double*>(q));
    double med=quart[1],iqr=quart[2]-quart[0];
    double lo=med-2.5*iqr , hi=med+2.5*iqr;
    if(lo>=hi){ lo=h->GetMean()-2*h->GetRMS(); hi=h->GetMean()+2*h->GetRMS(); }
    TF1 g("g","gaus",lo,hi);
    g.SetParameters(h->GetMaximum(),med,0.5*h->GetRMS());
    TFitResultPtr r=h->Fit(&g,"QNR0S");
    if(!r.Get()||!r->IsValid()){
      double rms=h->GetRMS(),err=rms/std::sqrt(2.*(h->GetEntries()-1));
      return {h->GetMean(),err,rms,err};
    }
    return {g.GetParameter(1),g.GetParError(1),
            std::fabs(g.GetParameter(2)),g.GetParError(2)};
  };

  /* ---------- 1. master canvas ---------- */
  const int nCols=4 , nRows=(binMode==EBinningMode::kDiscrete?4:2);
  TCanvas cMain("cΔφ5","Δφ five‑way comparison",1800,600*nRows);
  cMain.Divide(nCols,nRows);

  /* storage for summary graphs -------------------------------------- */
  std::vector<double> eCtr;
  enum {kRaw,kRawB,kCP,kCPcp,kCPb};
  std::array<std::vector<double>,5> mu, muE, sg, sgE;

  std::vector<TObject*> guard;   // prevent premature deletion

  /* ---------- 2. per‑slice loop ---------- */
  for(int i=0;i<N_E;++i)
  {
    const std::string tag=tagOf(i);

    TH1F *h[kCPb+1]={};
    h[kRaw ] = (TH1F*)gROOT->FindObject(Form("h_phi_diff_raw_%s"   ,tag.c_str()));
    h[kRawB] = (TH1F*)gROOT->FindObject(Form("h_phi_diff_corr_%s"  ,tag.c_str()));
    if(!isFirstPass){
      h[kCP ] = (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpRaw_%s" ,tag.c_str()));
      h[kCPcp]= (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpCorr_%s",tag.c_str()));
      h[kCPb] = (TH1F*)gROOT->FindObject(Form("h_phi_diff_cpBcorr_%s",tag.c_str()));
    }

    if(!h[kRaw]) { std::cerr<<"[Δφ‑5] missing RAW hist "<<tag<<"\n"; continue; }


    /* clone, normalise, fit, store ---------------------------------- */
    GRes gRes[5]={{}};
    for(int v=0;v<=kCPb;++v){
      if(!h[v]) continue;
      TH1F* c=(TH1F*)h[v]->Clone(Form("h_%d_%d",v,i));
      c->SetDirectory(nullptr);
      if(c->Integral()>0) c->Scale(1./c->Integral());
      c->SetMarkerColor(colArr[v]); c->SetLineColor(colArr[v]);
      c->SetMarkerStyle(mksArr[v]); c->SetMarkerSize(0.8);
      guard.push_back(c); h[v]=c;               // replace by clone
      gRes[v]=robust(c);
    }

    /* save results for summary plots ------------------------------- */
    eCtr.push_back(0.5*(eEdges[i].first+eEdges[i].second));
    auto push=[&](int v,const GRes& r){
      mu [v].push_back(r.mu ); muE[v].push_back(r.muErr);
      sg [v].push_back(r.sg ); sgE[v].push_back(r.sgErr);
    };
    for(int v=0;v<=kCPb;++v) if(h[v]) push(v,gRes[v]);

    /* ---------- draw overlay ---------- */
    TPad* pad=(TPad*)cMain.cd(i+1);
    pad->SetLeftMargin(0.23); pad->SetBottomMargin(0.18);

    /* find y‑max */
    double yMax=0; for(int v=0;v<=kCPb;++v) if(h[v]) yMax=std::max(yMax,h[v]->GetMaximum());

    h[kRaw]->SetTitle(Form("#Delta#phi  [%g, %g) GeV",
                           eEdges[i].first,eEdges[i].second));
    h[kRaw]->GetXaxis()->SetTitle("#Delta#phi = #phi_{reco} - #phi_{truth}  [rad]");
    h[kRaw]->GetYaxis()->SetTitle("Probability density");
    h[kRaw]->GetYaxis()->SetRangeUser(0,1.25*yMax);

    for(int v=0;v<=kCPb;++v)
      if(h[v]) h[v]->Draw(v==kRaw? "E":"E SAME");

    TLegend* lg=new TLegend(0.24,0.66,0.78,0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028);
    for(int v=0;v<=kCPb;++v) if(h[v]) lg->AddEntry(h[v],legTxtArr[v],"lp");
    lg->Draw(); guard.push_back(lg);
  } /* slice loop */

  /* ---------- helper: make graph ---------- */
  auto makeG=[&](int v,double dx,Color_t col,Style_t m)->TGraphErrors*
  {
    const int N=(int)eCtr.size();
    std::vector<double> x(N),ex(N,0.);
    for(int i=0;i<N;++i) x[i]=eCtr[i]+dx;
    TGraphErrors* g=new TGraphErrors(N,x.data(),mu[v].data(),ex.data(),muE[v].data());
    g->SetMarkerColor(col); g->SetLineColor(col);
    g->SetMarkerStyle(m); g->SetMarkerSize(1.1);
    guard.push_back(g); return g;
  };

    /* ---------- 3. summary canvas μ & σ ---------- */
    {
        /* offset between the 5 curves (-2 … +2) */
        const double dx = 0.0;

        /* ------------------------------------------------------------------ *
         *  Helper → return {min,max} over *all* curves, including ±σErr,
         *  then pad the range by 5 % so markers never touch the frame
         * ------------------------------------------------------------------ */
        auto rangeOf = [&](bool useMu)->std::pair<double,double>
        {
            double lo =  1e30, hi = -1e30;
            for (int v = 0; v <= kCPb; ++v) {
                if (mu[v].empty()) continue;
                const auto &val =  useMu ?  mu[v] :  sg[v];
                const auto &err =  useMu ?  muE[v] : sgE[v];
                for (size_t i = 0; i < val.size(); ++i) {
                    lo = std::min(lo, val[i] - err[i]);
                    hi = std::max(hi, val[i] + err[i]);
                }
            }
            const double pad = 0.05 * (hi - lo);
            return {lo - pad, hi + pad};
        };

        /* ------------------------------------------------------------------ *
         *  Canvas + two stacked pads
         * ------------------------------------------------------------------ */
        TCanvas cSum("cMuSigma5", "#Delta#phi five-way summary", 900, 780);
        TPad pT("pT","",0,0.37,1,1);  pT.Draw();
        TPad pB("pB","",0,0   ,1,0.37); pB.Draw();

        const double xMin =  eCtr.front() - 1.0;          // pad left by 1 GeV
        const double xMax =  eEdges.back().second;        // true upper E-edge

        /* ====================  μ(E) top pad  ==================== */
        pT.cd();
        pT.SetLeftMargin(0.15);
        pT.SetBottomMargin(0.03);               // hair-line gap

        auto [muMin, muMax] = rangeOf(true);    // auto-scaled y-range
        TH1F fMu("fMu","; ;Gaussian mean  #mu  [rad]",1,xMin,xMax);
        fMu.SetMinimum(muMin);
        fMu.SetMaximum(muMax);

        // hide x–axis labels/ticks on the top pad
        TAxis* axT = fMu.GetXaxis();
        axT->SetLabelSize(0);  axT->SetTickLength(0);  axT->SetTitleSize(0);
        fMu.Draw("AXIS");

        /* dashed reference line at μ = 0 */
        TLine *l0 = new TLine(xMin, 0.0, xMax, 0.0);
        l0->SetLineStyle(2); l0->SetLineWidth(2); l0->SetLineColor(kGray+2);
        l0->Draw();                 guard.push_back(l0);

        /* plot all five curves – circles, coloured, x-offset by (v-2)*dx */
        for (int v = 0; v <= kCPb; ++v) {
            if (mu[v].empty()) continue;
            TGraphErrors *g = makeG(v, (v-2)*dx, colArr[v], 20);   // 20 = solid circle
            g->Draw("P SAME");
        }

        /* legend (use invisible TMarker clones to show colour) */
        TLegend lg(0.4,0.08,0.82,0.28);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.032);
        for (int v = 0; v <= kCPb; ++v) {
            if (mu[v].empty()) continue;
            auto m = new TMarker(0,0,20);           // circle
            m->SetMarkerColor(colArr[v]); m->SetMarkerSize(1.1);
            guard.push_back(m);
            lg.AddEntry(m, legTxtArr[v], "p");
        }
        lg.Draw();

        /* ====================  σ(E) bottom pad  ==================== */
        pB.cd();
        pB.SetLeftMargin(0.15);
        pB.SetTopMargin(0.07);
        pB.SetBottomMargin(0.35);

        auto [siMin, siMax] = rangeOf(false);
        TH1F fSi("fSi",
                 ";Energy slice centre  [GeV];Gaussian width  #sigma  [rad]",
                 1, xMin, xMax);
        fSi.SetMinimum(std::max(0.0, siMin));
        fSi.SetMaximum(siMax);
        fSi.Draw("AXIS");

        for (int v = 0; v <= kCPb; ++v) {
            if (sg[v].empty()) continue;

            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x(N), ex(N,0.);
            for (int i = 0; i < N; ++i) x[i] = eCtr[i] + (v-2)*dx;

            TGraphErrors *g = new TGraphErrors(N, x.data(), sg[v].data(),
                                               ex.data(), sgE[v].data());
            g->SetMarkerStyle(20);                    // circle
            g->SetMarkerColor(colArr[v]);
            g->SetLineColor  (colArr[v]);
            g->SetMarkerSize(1.1);
            guard.push_back(g);
            g->Draw("P SAME");
        }

        /* ------------------------------------------------------------------ */
        cSum.SaveAs(TString(outDir) + "/DeltaPhi5_MeanSigmaVsE.png");
    }


    /* ---------- 4. diagnostics :  RMSE  &  σ / |μ|  ---------- */
    auto diag = [&](const char* file,                 // PNG stem
                    const char* yTit,                 // y–axis label
                    const std::array<std::vector<double>,5>& y,
                    const std::array<std::vector<double>,5>& yE)
    {
        const double dxShift = 0.12;                  // ±x-offset between curves
        TCanvas c(file, file, 900, 640);
        c.SetLeftMargin(0.15);  c.SetRightMargin(0.06);

        /* ---------- robust y-range : include ±error & add 7 % padding ----- */
        double yMin =  1e30,  yMax = -1e30;
        for (int v = 0; v <= kCPb; ++v)
            for (size_t i = 0; i < y[v].size(); ++i) {
                yMin = std::min(yMin, y[v][i] - yE[v][i]);
                yMax = std::max(yMax, y[v][i] + yE[v][i]);
            }
        if (yMin > 0.0) yMin = 0.0;                   // RMSE / frac-res are ≥ 0
        const double pad = 0.07 * (yMax - yMin);      // 7 % breathing room
        yMin -= pad;   yMax += pad;

        /* ---------- frame -------------------------------------------------- */
        TH1F frame("fr",
                   Form(";Energy slice centre  [GeV];%s", yTit),
                   1, eCtr.front() - 1.0,  eEdges.back().second);   // true upper edge (30 GeV)
        frame.SetMinimum(yMin);
        frame.SetMaximum(yMax);
        frame.Draw("AXIS");

        /* ---------- legend (bottom-left) ----------------------------------- */
        /* ---------- legend --------------------------------------------------- */
        /*  – RMSE stays bottom-left
         *  – FracRes goes to *top-left* with a slightly smaller font            */
        const bool isFracRes = (std::strcmp(file, "DeltaPhi5_FracRes") == 0);

        const double x1 = 0.15;                      // left edge identical
        const double x2 = 0.55;                      // right edge identical

        double y1, y2;                               // decide vertical placement
        if (isFracRes) {    /* top-left */
            y1 = 0.75;                               // NDC
            y2 = 0.92;
        } else {            /* keep RMSE bottom-left */
            y1 = 0.2;
            y2 = 0.4;
        }

        TLegend lg(x1, y1, x2, y2);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.024);       // a touch smaller than before
        lg.SetTextFont(42);          // plain font for readability

        /* ---------- draw all five curves – circles only -------------------- */
        for (int v = 0; v <= kCPb; ++v)
        {
            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x (N), ex(N, 0.);
            for (int i = 0; i < N; ++i) x[i] = eCtr[i] + (v - 2) * dxShift;

            TGraphErrors* g = new TGraphErrors(
                                 N, x.data(), y[v].data(),
                                 ex.data(),   yE[v].data());

            g->SetMarkerStyle(20);          // solid circle
            g->SetMarkerSize (1.1);
            g->SetMarkerColor(colArr[v]);
            g->SetLineColor  (colArr[v]);

            guard.push_back(g);
            g->Draw("P SAME");

            lg.AddEntry(g, legTxtArr[v], "p");
        }
        lg.Draw();

        c.SaveAs( (TString(outDir) + "/" + file + ".png").Data() );
    };

    /* ---------- prepare RMSE and  σ/|μ|  ---------------------------------- */
    std::array<std::vector<double>,5> rm , rmE ,
                                      fr , frE ;

    for (int v = 0; v <= kCPb; ++v) {
        rm [v].resize(eCtr.size());  rmE[v].resize(eCtr.size());
        fr [v].resize(eCtr.size());  frE[v].resize(eCtr.size());
    }

    for (size_t i = 0; i < eCtr.size(); ++i) {
        for (int v = 0; v <= kCPb; ++v) {

            /* grab the central values and their errors -------------------- */
            double muVal   =  mu [v][i];
            double sgVal   =  sg [v][i];
            double dmuVal  = muE[v][i];
            double dsgVal  = sgE[v][i];

            /* ---------- RMSE = √(μ² + σ²) -------------------------------- */
            double rmse    = std::hypot(muVal, sgVal);
            double rmseErr = rmse * std::sqrt( std::pow(dmuVal/muVal,2) +
                                               std::pow(dsgVal/sgVal,2) );
            rm [v][i] = rmse;
            rmE[v][i] = rmseErr;

            /* ---------- fractional resolution  σ / |μ| ------------------- */
            double frac    = sgVal / std::fabs(muVal);
            double fracErr = frac  * std::sqrt( std::pow(dsgVal/sgVal,2) +
                                                std::pow(dmuVal/muVal,2) );
            fr [v][i] = frac;
            frE[v][i] = fracErr;
        }
    }


    /* ---------- draw and save the two summary plots ----------------------- */
    diag("DeltaPhi5_RMSE",    "#sqrt{#mu^{2} + #sigma^{2}}  [rad]", rm, rmE);
    diag("DeltaPhi5_FracRes", "#sigma / |#mu|",                     fr, frE);


  /* ---------- 5. save slice sheet ---------- */
  cMain.SaveAs(TString(outDir)+"/DeltaPhi5_Compare_AllSlices.png");
  std::cout<<"[Δφ‑5] overlays & summaries written to "<<outDir<<'\n';

  for(auto* o:guard) delete o;
}



/******************************************************************************
 * OverlayDeltaEtaSlices  (v1 – cloned from OverlayDeltaPhiSlices)
 *
 *  – Consumes the Δη histograms produced by
 *      • fillDEtaRawAndCorrected()
 *      – names:  h_eta_diff_raw_<tag>,  h_eta_diff_corr_<tag>
 *
 *  – Produces exactly the same suite of plots that OverlayDeltaPhiSlices
 *    makes for Δφ, just with “η” everywhere:
 *        ΔηOverlay_….png           (per–slice overlay sheet)
 *        DeltaEta_MeanSigmaVsE.png (μ/σ summary, zero‑line + two pads)
 *        DeltaEta_RMSErrorVsE.png
 *        DeltaEta_FracResVsE.png
 *        DeltaEtaCompare_AllOutput.png   (φ‑style ‘all‑output’ sheet)
 *
 *  Everything else – fit strategy, zero‑line, automatic axis scaling,
 *  marker styles, legend placement – is byte‑for‑byte identical.
 *
 *  2025‑06‑15  – first public version
 ******************************************************************************/
void OverlayDeltaEtaSlices(const std::vector<std::pair<double,double>>& eEdges,
                           EBinningMode  binMode,
                           bool          isFirstPass,
                           const char*   outDir)
{
    /* ------------------------------------------------------------------ *
     * 0)  Boiler‑plate
     * ------------------------------------------------------------------ */
    gSystem->mkdir(outDir, /*recursive=*/true);
    gStyle->SetOptStat(0);

    const int N_E = static_cast<int>(eEdges.size());
    if (!N_E) { std::cerr << "[Δη] eEdges empty – abort\n"; return; }

    /* reproducible energy‑slice label ----------------------------------- */
    auto makeLabel = [&](int i)->std::string
    {
        return (binMode==EBinningMode::kRange)
             ? Form("%.0f_%.0f", eEdges[i].first, eEdges[i].second)
             : Form("E%.0f",      eEdges[i].first);
    };

    /* ------------------------------------------------------------------ *
     * 1)  Bullet‑proof 1‑D Gaussian fitter (identical to Δφ routine)
     * ------------------------------------------------------------------ */
    struct GRes { double N, mu, muErr, sg, sgErr; };

    auto robustGauss = [](TH1* h)->GRes
    {
        if (!h || h->Integral()==0) return {0,0,0,0,0};

        const double q[3]={0.25,0.50,0.75};
        double quart[3]; h->GetQuantiles(3,quart,const_cast<double*>(q));
        const double med=quart[1], iqr=quart[2]-quart[0];

        double lo=med-2.5*iqr, hi=med+2.5*iqr;
        if (lo>=hi){ lo=h->GetMean()-2.*h->GetRMS(); hi=h->GetMean()+2.*h->GetRMS(); }

        TF1 f("g","gaus",lo,hi);

        const double A=h->GetMaximum(), S=std::max(1e-4,0.5*h->GetRMS());
        std::vector<std::tuple<double,double,double>> seed={
            {A      ,med      ,S    },{1.5*A,med,1.5*S},
            {0.5*A  ,med+0.3*S,0.7*S},{1.2*A,med-0.3*S,1.2*S} };

        int nIter=0; double best=1e99;
        GRes res{h->GetEntries(),0,0,0,0};
        while(++nIter<=3)
        {
            for(auto s:seed){
                f.SetParameters(std::get<0>(s),std::get<1>(s),std::get<2>(s));
                TFitResultPtr r=h->Fit(&f,"QNR0S");
                if(!r.Get()||!r->IsValid()) continue;
                if(r->Chi2()<best && std::fabs(f.GetParameter(2))>1e-6){
                    best=r->Chi2();
                    res.mu=f.GetParameter(1); res.muErr=f.GetParError(1);
                    res.sg=std::fabs(f.GetParameter(2)); res.sgErr=f.GetParError(2);
                }
            }
            if(best<1e98) break;
            lo+=0.5*f.GetParameter(2); hi-=0.5*f.GetParameter(2);
            if(lo>=hi||h->Integral()<50) break;
            f.SetRange(lo,hi);
        }
        if(best>=1e98){                      // fallback robust RMS
            const double rms=h->GetRMS();
            const double err=rms/std::sqrt(2.*(h->GetEntries()-1));
            res.mu=h->GetMean(); res.muErr=err; res.sg=rms; res.sgErr=err;
        }
        return res;
    };

    /* ------------------------------------------------------------------ *
     * 2)  Per‑slice overlay canvas
     * ------------------------------------------------------------------ */
    const int nCols=4, nRows=(binMode==EBinningMode::kDiscrete ? 4 : 2);
    TCanvas c4x2(Form("DeltaEtaOverlay_%dx%d",nCols,nRows),
                 "#Delta#eta raw vs corrected (normalised)",
                 1600,600*nRows);
    c4x2.Divide(nCols,nRows);

    /* book‑keeping containers ------------------------------------------ */
    std::vector<double> eCtr,
                        muRaw,muRawErr, sgRaw,sgRawErr,
                        muCor,muCorErr, sgCor,sgCorErr;

    std::vector<TObject*> keep;      // prevent premature deletion

    /* ------------------------------------------------------------------ *
     * 3)  Loop over slices – identical drawing/fit logic
     * ------------------------------------------------------------------ */
    for(int iE=0;iE<N_E;++iE)
    {
        const double eLo=eEdges[iE].first,  eHi=eEdges[iE].second;
        const std::string tag=makeLabel(iE);

        TH1F* hRaw = dynamic_cast<TH1F*>( gROOT->FindObject(
                         Form("h_eta_diff_raw_%s",tag.c_str()) ));
        TH1F* hCor = (!isFirstPass)
                     ? dynamic_cast<TH1F*>( gROOT->FindObject(
                         Form("h_eta_diff_corr_%s",tag.c_str()) ))
                     : nullptr;
        if(!hRaw){ std::cerr<<"[Δη] missing RAW hist "<<tag<<'\n'; continue; }

        hRaw = (TH1F*)hRaw->Clone(Form("hRaw_%d",iE));  hRaw->SetDirectory(nullptr);
        keep.push_back(hRaw);
        if(hCor){ hCor=(TH1F*)hCor->Clone(Form("hCor_%d",iE)); hCor->SetDirectory(nullptr); keep.push_back(hCor); }

        if(hRaw->Integral()>0) hRaw->Scale(1./hRaw->Integral());
        if(hCor && hCor->Integral()>0) hCor->Scale(1./hCor->Integral());

        GRes R=robustGauss(hRaw);
        GRes C=hCor?robustGauss(hCor):GRes{0,0,0,0,0};

        eCtr.push_back(0.5*(eLo+eHi));
        muRaw.push_back(R.mu); muRawErr.push_back(R.muErr);
        sgRaw.push_back(R.sg); sgRawErr.push_back(R.sgErr);
        muCor.push_back(C.mu); muCorErr.push_back(C.muErr);
        sgCor.push_back(C.sg); sgCorErr.push_back(C.sgErr);

        /* ---------- drawing ------------------------------------------- */
        TPad* pad=(TPad*)c4x2.cd(iE+1);
        pad->SetLeftMargin(0.22); pad->SetBottomMargin(0.18);
        pad->SetRightMargin(0.06); pad->SetTopMargin(0.10);

        hRaw->SetMarkerStyle(20); hRaw->SetMarkerColor(kBlack); hRaw->SetLineColor(kBlack);
        if(hCor){ hCor->SetMarkerStyle(20); hCor->SetMarkerColor(kRed); hCor->SetLineColor(kRed); }

        double yMax=hRaw->GetMaximum(); if(hCor) yMax=std::max(yMax,hCor->GetMaximum());
        hRaw->GetYaxis()->SetRangeUser(0,1.25*yMax);
        hRaw->SetTitle(Form("#Delta#eta  [%.1f, %.1f) GeV",eLo,eHi));
        hRaw->GetXaxis()->SetTitle("#Delta#eta  (reco - truth)");
        hRaw->GetYaxis()->SetTitle("Normalised counts");

        hRaw->Draw("E"); if(hCor) hCor->Draw("E SAME");

        TLegend* lg=new TLegend(0.24,0.70,0.48,0.86);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
        lg->AddEntry(hRaw,"RAW","lp"); if(hCor) lg->AddEntry(hCor,"CORRECTED","lp");
        lg->Draw(); keep.push_back(lg);

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030); tx.SetTextAlign(33);
        tx.DrawLatex(0.92,0.88,Form("RAW  #mu=%.4f  #sigma=%.4f",R.mu,R.sg));
        if(hCor) tx.DrawLatex(0.92,0.82,Form("CORR #mu=%.4f  #sigma=%.4f",C.mu,C.sg));
    } /* slice loop */

    /* ------------------------------------------------------------------ *
     * 4)  μ(E) & σ(E) summary (identical styling incl. zero‑line)
     * ------------------------------------------------------------------ */
    const int nPts=(int)eCtr.size(); if(!nPts) return;
    auto makeG=[&](const std::vector<double>& y,const std::vector<double>& ye,
                   double dx,Color_t col)->TGraphErrors*
    {
        std::vector<double> x(nPts), ex(nPts,0.);
        for(int i=0;i<nPts;++i) x[i]=eCtr[i]+dx;
        auto g=new TGraphErrors(nPts,x.data(),y.data(),ex.data(),ye.data());
        g->SetMarkerStyle(20); g->SetMarkerSize(1.1);
        g->SetMarkerColor(col); g->SetLineColor(col);
        keep.push_back(g); return g;
    };
    const double dx=0.08;
    auto gMuR=makeG(muRaw,muRawErr,-dx,kBlack);
    auto gMuC=makeG(muCor,muCorErr,+dx,kRed);
    auto gSiR=makeG(sgRaw,sgRawErr,-dx,kBlack);
    auto gSiC=makeG(sgCor,sgCorErr,+dx,kRed);

    /* summary canvas --------------------------------------------------- */
    {
        const double xMin=eEdges.front().first-1., xMax=eEdges.back().second;
        TCanvas cSum("cMuSiEta","#Delta#eta mean / sigma vs E",860,760);
        TPad pT("pT","",0,0.34,1,1); pT.Draw();
        TPad pB("pB","",0,0   ,1,0.34); pB.Draw();

        /* μ(E) top pad -------------------------------------------------- */
        auto fullMin=[&](const std::vector<double>& a,const std::vector<double>& ae,
                         const std::vector<double>& b,const std::vector<double>& be)
        {
            double m=1e30; for(size_t i=0;i<a.size();++i) m=std::min(m,a[i]-ae[i]);
            for(size_t i=0;i<b.size();++i) m=std::min(m,b[i]-be[i]); return m;
        };
        auto fullMax=[&](const std::vector<double>& a,const std::vector<double>& ae,
                         const std::vector<double>& b,const std::vector<double>& be)
        {
            double M=-1e30; for(size_t i=0;i<a.size();++i) M=std::max(M,a[i]+ae[i]);
            for(size_t i=0;i<b.size();++i) M=std::max(M,b[i]+be[i]); return M;
        };
        const double muMin=fullMin(muRaw,muRawErr,muCor,muCorErr);
        const double muMax=fullMax(muRaw,muRawErr,muCor,muCorErr);

        pT.cd(); pT.SetLeftMargin(0.15); pT.SetBottomMargin(0.03);
        TH1F fMu("fMu","; ;#mu_{Gauss}  [rad]",1,xMin,xMax);
        fMu.SetMinimum(muMin-0.1*std::fabs(muMin)); fMu.SetMaximum(muMax+0.1*std::fabs(muMax));
        fMu.GetXaxis()->SetLabelSize(0); fMu.GetXaxis()->SetTickLength(0); fMu.GetXaxis()->SetTitleSize(0);
        fMu.Draw("AXIS");
        gMuR->Draw("P SAME"); gMuC->Draw("P SAME");

        TLine* l0=new TLine(xMin,0.,xMax,0.); l0->SetLineStyle(2); l0->SetLineWidth(2); l0->SetLineColor(kGray+2);
        l0->Draw(); keep.push_back(l0);

        TLegend legMu(0.18,0.79,0.43,0.91);
        legMu.SetBorderSize(0); legMu.SetFillStyle(0); legMu.SetTextSize(0.035);
        legMu.AddEntry(gMuR,"RAW","lp"); legMu.AddEntry(gMuC,"CORRECTED","lp"); legMu.Draw();

        /* σ(E) bottom pad ---------------------------------------------- */
        const double siMin=fullMin(sgRaw,sgRawErr,sgCor,sgCorErr);
        const double siMax=fullMax(sgRaw,sgRawErr,sgCor,sgCorErr);

        pB.cd(); pB.SetLeftMargin(0.15); pB.SetTopMargin(0.06); pB.SetBottomMargin(0.38);
        TH1F fSi("fSi",";E_{slice centre}  [GeV];#sigma_{Gauss}  [rad]",1,xMin,xMax);
        fSi.SetMinimum(std::max(0.,siMin-0.1*std::fabs(siMin)));
        fSi.SetMaximum(siMax+0.1*std::fabs(siMax));
        fSi.GetXaxis()->SetNdivisions(505); fSi.GetXaxis()->SetTickLength(0.05); fSi.GetXaxis()->SetTitleOffset(1.1);
        fSi.Draw("AXIS");
        gSiR->Draw("P SAME"); gSiC->Draw("P SAME");

        cSum.SaveAs(TString(outDir)+"/DeltaEta_MeanSigmaVsE.png");
    }

    /* ------------------------------------------------------------------ *
     * 5)  EXTRA diagnostics  (RMSE & σ/|μ|)
     * ------------------------------------------------------------------ */
    {
        const int n=nPts;
        std::vector<double> rmR(n),rmRE(n), rmC(n),rmCE(n),
                            frR(n),frRE(n), frC(n),frCE(n);
        for(int i=0;i<n;++i){
            rmR[i]=std::hypot(muRaw[i],sgRaw[i]);
            rmRE[i]=rmR[i]*std::sqrt(std::pow(muRawErr[i]/muRaw[i],2)+
                                     std::pow(sgRawErr[i]/sgRaw[i],2));
            rmC[i]=std::hypot(muCor[i],sgCor[i]);
            rmCE[i]=(rmC[i]>0)?rmC[i]*std::sqrt(std::pow(muCorErr[i]/muCor[i],2)+
                                                std::pow(sgCorErr[i]/sgCor[i],2)):0;

            frR[i]=sgRaw[i]/std::fabs(muRaw[i]);
            frRE[i]=frR[i]*std::sqrt(std::pow(sgRawErr[i]/sgRaw[i],2)+
                                     std::pow(muRawErr[i]/muRaw[i],2));
            frC[i]=sgCor[i]/std::fabs(muCor[i]);
            frCE[i]=frC[i]*std::sqrt(std::pow(sgCorErr[i]/sgCor[i],2)+
                                     std::pow(muCorErr[i]/muCor[i],2));
        }

        auto makeG=[&](const std::vector<double>& y,const std::vector<double>& ye,
                       double dx,Color_t col)->TGraphErrors*
        {
            std::vector<double> x(n),ex(n,0.); for(int i=0;i<n;++i) x[i]=eCtr[i]+dx;
            auto g=new TGraphErrors(n,x.data(),y.data(),ex.data(),ye.data());
            g->SetMarkerStyle(20); g->SetMarkerSize(0.9);
            g->SetMarkerColor(col); g->SetLineColor(col); keep.push_back(g); return g;
        };

        auto drawDiag=[&](const char* name,const char* yTit,
                          const std::vector<double>& yR,const std::vector<double>& yRE,
                          const std::vector<double>& yC,const std::vector<double>& yCE,
                          const char* png,bool topRight)
        {
            const double dxOff=0.12;
            auto gR=makeG(yR,yRE,-dxOff,kBlack);
            auto gC=makeG(yC,yCE,+dxOff,kRed);

            double yMin=1e30,yMax=-1e30;
            for(size_t i=0;i<yR.size();++i){
                yMin=std::min({yMin,yR[i]-yRE[i],yC[i]-yCE[i]});
                yMax=std::max({yMax,yR[i]+yRE[i],yC[i]+yCE[i]});
            }
            const double span=yMax-yMin; yMin-=0.02*span; yMax+=0.05*span;

            TCanvas c(name,name,860,620); c.SetLeftMargin(0.15); c.SetRightMargin(0.06);
            TH1F fr("fr",Form(";E_{slice} [GeV];%s",yTit),1,eCtr.front()-1.0,eEdges.back().second);
            fr.SetMinimum(yMin); fr.SetMaximum(yMax); fr.Draw("AXIS");
            gR->Draw("P SAME"); gC->Draw("P SAME");

            double lx1= topRight?0.70:0.15 , ly1= topRight?0.72:0.12;
            double lx2= topRight?0.90:0.45 , ly2= topRight?0.92:0.32;
            TLegend lg(lx1,ly1,lx2,ly2);
            lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.028);
            lg.AddEntry(gR,"RAW","lp"); lg.AddEntry(gC,"CORRECTED","lp"); lg.Draw();

            c.SaveAs(TString(outDir)+"/"+png);
        };

        drawDiag("cRMSEeta","#sqrt{#mu^{2}+#sigma^{2}}  [rad]",
                 rmR,rmRE, rmC,rmCE, "DeltaEta_RMSErrorVsE.png", /*topRight=*/false);

        drawDiag("cFracEta","#sigma / |#mu|",
                 frR,frRE, frC,frCE, "DeltaEta_FracResVsE.png",   /*topRight=*/true);
    }

    /* ------------------------------------------------------------------ *
     * 6)  all‑output sheet
     * ------------------------------------------------------------------ */
    TString outAll=TString(outDir)+"/DeltaEtaCompare_AllOutput.png";
    c4x2.SaveAs(outAll);
    std::cout<<"[Δη] wrote summary → "<<outAll<<'\n';

    for(auto* o:keep) delete o;
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
 *  – ⟨E⟩ =  global mean cluster‑energy                                     *
 *                                                                          *
 *  NEW (24‑Jun‑2025)                                                       *
 *  ─────────────────                                                       *
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
    t.DrawLatex(0.97,0.96,txt);
  };

/* ---------------------------------------------------------------------- */
/* 1. build residual map                                                  */
/* ---------------------------------------------------------------------- */
  std::unique_ptr<TH2D> meanMap( h->Project3DProfile("yx") );
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
/* 2. 2‑D colour map (single histogram)                                   */
/* ---------------------------------------------------------------------- */
  TCanvas c2(Form("c_res_%s",tag),"",1000,850);
  c2.SetLeftMargin(0.18); c2.SetBottomMargin(0.16); c2.SetRightMargin(0.12);

  gStyle->SetPalette(kBird);     gStyle->SetNumberContours(255);

  res->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
  res->GetYaxis()->SetTitle("block #phi_{local, 2#times 2}");
  res->GetZaxis()->SetTitle("|#DeltaE| / #LT E #GT  [%]");
  res->Draw("COLZ");  c2.Update();

  /* palette cosmetics */
  if (auto* pal = dynamic_cast<TPaletteAxis*>(
        res->GetListOfFunctions()->FindObject("palette")))
  {
    pal->SetX1NDC(0.92); pal->SetX2NDC(0.945);
    pal->GetAxis()->SetTitle("");  pal->SetLabelSize((fLbl-0.002)*0.75);
    pal->SetBorderSize(0); pal->SetFillStyle(0);

    const double xMid=0.5*(pal->GetX1NDC()+pal->GetX2NDC());
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.022);
    t.SetTextAlign(21);
    t.DrawLatex(xMid,pal->GetY2NDC()+0.018,"|#DeltaE| / #LT E #GT  [%]");
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
    /*  helper #1 – two‑pad 2‑D comparison (single palette)               */
    /* ================================================================== */
    auto makeSide2D = [&](TH2* unc, TH2* cor, const char* outPng)
    {
      /* ---------- canvas skeleton ------------------------------------ */
      TCanvas c("cSide2D","",1800,850);
      c.Divide(2,1,0,0);                // equal‑sized pads

      /* identical Z‑range --------------------------------------------- */
      const double zMax = std::max(unc->GetMaximum(),cor->GetMaximum());
      unc->SetMaximum(zMax);  cor->SetMaximum(zMax);

      /* ---------- left pad – UNCORRECTED ----------------------------- */
      c.cd(1);
      gPad->SetLeftMargin(0.18);
      gPad->SetRightMargin(0.06);       // slim, but equal to right pad
      gPad->SetBottomMargin(0.16);

      unc->Draw("COL");                 // no palette in this pad
      drawTag("UNCORRECTED");

      /* ---------- right pad – CORRECTED (+ palette) ------------------ */
      c.cd(2);
      gPad->SetLeftMargin(0.18);
      gPad->SetRightMargin(0.06);       // same absolute size as pad 1
      gPad->SetBottomMargin(0.16);

      cor->Draw("COLZ");
      drawTag("CORRECTED");

      /* shift palette so its ticks are fully visible ------------------ */
      gPad->Update();
      if (auto* pal = dynamic_cast<TPaletteAxis*>(
            cor->GetListOfFunctions()->FindObject("palette")))
      {
        pal->SetX1NDC(0.93);            // snug to the pad’s right edge
        pal->SetX2NDC(0.97);
      }

      c.Print(outPng,"png 600");
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

      TLegend L(0.65,0.80,0.90,0.92);
      L.SetBorderSize(0);  L.SetFillStyle(0);
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


void makeLegoGifHD(TH3 *h,
                   const char *tag,
                   const char *hdr,
                   const char *outDir,
                   int   nFrames = 180,       // 2°/frame → one full orbit in 7.2 s @ 25 fps
                   double theta0  = 28.0,     // inclination at φ = φ0
                   double theta1  = 38.0,     // inclination at φ = φ0 + φArc
                   double phi0    =   0.0,    // start azimuth (deg)
                   double phiArc  = 360.0,    // azimuth sweep (deg)
                   int   SS       = 2,        // supersampling factor (×2 = good AA, small IO)
                   int   W0       = 1280,     // output width  (px)
                   int   H0       = 720,      // output height (px)
                   int   FPS      = 25)       // frame‑rate (Google Slides default)
{
    /*----------------------------------------------------------------------*/
    /* 0. sanity checks                                                     */
    /*----------------------------------------------------------------------*/
    if (!h)       { Error("makeLegoGifHD", "null histogram pointer"); return; }
    if (nFrames<2){ Error("makeLegoGifHD", "nFrames must be ≥ 2");     return; }

    /*----------------------------------------------------------------------*/
    /* 1. ensure batch mode and thread safety                               */
    /*----------------------------------------------------------------------*/
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,30,0)
    // Avoid the implicit multi‑threading ROOT enables since 6.30
    ROOT::EnableImplicitMT(0);
#endif
    gROOT->SetBatch(kTRUE);

    /*----------------------------------------------------------------------*/
    /* 2. one‑off ROOT style                                                */
    /*----------------------------------------------------------------------*/
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(140);
    gStyle->SetCanvasPreferGL(true);      // use OpenGL → MSAA
    gErrorIgnoreLevel = kWarning;         // suppress harmless colour warnings

    /*----------------------------------------------------------------------*/
    /* 3. supersampled canvas                                               */
    /*----------------------------------------------------------------------*/
    const int W = W0 * SS,  H = H0 * SS;
    TString canvName = Form("c_%s_gif", tag);
    std::unique_ptr<TCanvas> c(new TCanvas(canvName, "", W, H));
    c->SetRightMargin (0.15);
    c->SetBottomMargin(0.13);
    c->SetFillColorAlpha(kWhite, 0.0);

    /*----------------------------------------------------------------------*/
    /* 4. static scene (draw only once)                                     */
    /*----------------------------------------------------------------------*/
    h->GetXaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetTitleOffset(1.55);
    h->GetZaxis()->SetTitleOffset(0.90);
    h->Draw("LEGO2Z0");
    gPad->Modified(); gPad->Update();

    if (auto *pal =
        dynamic_cast<TPaletteAxis*>(h->GetListOfFunctions()->FindObject("palette")))
    {
        pal->SetX1NDC(0.86); pal->SetX2NDC(0.89);
        pal->SetBorderSize(0);
        pal->SetFillStyle(0);
        pal->SetLabelSize(0.022 * SS / 2.0);
    }

    TLatex head; head.SetNDC(); head.SetTextFont(42);
    head.SetTextAlign(13);
    head.SetTextSize(0.05 * SS / 2.0);
    head.DrawLatex(0.04, 0.965, Form("#bf{%s}", hdr));
    gPad->Modified(); gPad->Update();

    /*----------------------------------------------------------------------*/
    /* 5. filenames & directories                                           */
    /*----------------------------------------------------------------------*/
    TString base   = Form("%s/lego_%s_spin", outDir, tag);
    TString pngDir = base + "_png";
    TString gifOut = base + ".gif";
    gSystem->mkdir(pngDir, kTRUE);

    /*----------------------------------------------------------------------*/
    /* 6. render PNG frames (8‑bit palette)                                 */
    /*----------------------------------------------------------------------*/
    Printf("\n[makeLegoGifHD]  %s", hdr);
    Printf("  output dir  : %s", outDir);
    Printf("  resolution  : %d × %d  (supersampled ×%d)", W0, H0, SS);
    Printf("  frames/fps  : %d @ %d fps\n------------------------------", nFrames, FPS);

    const char spin[4] = {'|','/','-','\\'};
    const TDatime tStart;

    for (int i=0;i<nFrames;++i)
    {
        double u   = double(i)/(nFrames-1);
        double phi = phi0 + phiArc*u;
        double th  = theta0 + (theta1-theta0)*
                     0.5*(1-std::cos(u*TMath::Pi()));   // cosine‑ease

        gPad->SetPhi(phi); gPad->SetTheta(th);
        gPad->Modified();  gPad->Update();

        TString fn = Form("%s/f%04d.png", pngDir.Data(), i);
        c->Print(fn, "png");  // 8‑bit, ~1 MB at 1280×720

        if (i%(nFrames/100<1?1:nFrames/100)==0||i==nFrames-1)
        {
            int elapsed = TDatime().Convert()-tStart.Convert();
            int eta     = elapsed*(nFrames-i-1)/(i+1);
            printf("\r  %c %3d %%  %4d s ETA", spin[i%4], int(100.*(i+1)/nFrames), eta);
            fflush(stdout);
        }
    }
    printf("\r  ✔ frame rendering done.                    \n");

    /*----------------------------------------------------------------------*/
    /* 7. encode GIF via ffmpeg > ImageMagick; abort if neither present     */
    /*----------------------------------------------------------------------*/
    TString ffmpeg  = gSystem->Which(nullptr,"ffmpeg");
    TString convert = gSystem->Which(nullptr,"convert");

    /* Fallback #1 – common Home-brew prefix (Apple/Intel) */
    if (ffmpeg .IsNull() && !gSystem->AccessPathName("/opt/homebrew/bin/ffmpeg", kExecutePermission))
        ffmpeg  = "/opt/homebrew/bin/ffmpeg";

    if (convert.IsNull() && !gSystem->AccessPathName("/opt/homebrew/bin/convert", kExecutePermission))
        convert = "/opt/homebrew/bin/convert";

    /* Fallback #2 – classic /usr/local for older Intel Macs */
    if (ffmpeg .IsNull() && !gSystem->AccessPathName("/usr/local/bin/ffmpeg", kExecutePermission))
        ffmpeg  = "/usr/local/bin/ffmpeg";

    if (convert.IsNull() && !gSystem->AccessPathName("/usr/local/bin/convert", kExecutePermission))
        convert = "/usr/local/bin/convert";

    if (!ffmpeg.IsNull())   /*------------ ffmpeg ––––––––––––––––––*/
    {
        TString pal = pngDir+"/palette.png";
        TString cmd  = Form("%s -loglevel error -y -i %s/f%%04d.png "
                            "-vf palettegen=max_colors=256 %s",
                            ffmpeg.Data(), pngDir.Data(), pal.Data());
        if (gSystem->Exec(cmd)!=0) { Error("makeLegoGifHD","ffmpeg palettegen failed"); return; }

        cmd = Form("%s -loglevel error -y -framerate %d -i %s/f%%04d.png -i %s "
                   "-lavfi paletteuse -gifflags +transdiff -r %d %s",
                   ffmpeg.Data(), FPS, pngDir.Data(), pal.Data(), FPS, gifOut.Data());
        if (gSystem->Exec(cmd)!=0) { Error("makeLegoGifHD","ffmpeg paletteuse failed"); return; }

        gSystem->Unlink(pal);
        Printf("  ✔ encoded with ffmpeg (%s)", ffmpeg.Data());
    }
    else if (!convert.IsNull())   /*------- ImageMagick convert ––––*/
    {
        int delay = lround(100.0/FPS);  // centiseconds
        TString cmd = Form("%s -delay %d -loop 0 %s/f*.png "
                           "-layers Optimize %s",
                           convert.Data(), delay, pngDir.Data(), gifOut.Data());
        if (gSystem->Exec(cmd)!=0) { Error("makeLegoGifHD","convert failed"); return; }
        Printf("  ✔ encoded with ImageMagick (%s)", convert.Data());
    }
    else
    {
        Error("makeLegoGifHD",
              "neither ffmpeg nor ImageMagick found – install one to produce GIFs");
        return;
    }

    /*----------------------------------------------------------------------*/
    /* 8. cleanup                                                           */
    /*----------------------------------------------------------------------*/
    gSystem->Exec(Form("rm -rf %s", pngDir.Data()));
    Printf("  ↪  %s\n  ▸ ready for Google Slides / Keynote (loops ∞)\n", gifOut.Data());
}






void PDCanalysis()
{
  // 1) Style / stat
  gROOT->LoadMacro("sPhenixStyle.C");
  SetsPhenixStyle();
    
  gStyle->SetOptStat(0);
    //_withVirgilesChange
  const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL_withoutPhiTiltCorr.root";

  // 2) Open input
  std::cout << "[INFO] Opening file: " << inFile << "\n";
  TFile* fIn = TFile::Open(inFile, "READ");
  if(!fIn || fIn->IsZombie())
  {
    std::cerr << "[ERROR] Could not open file: " << inFile << std::endl;
    return;
  }
  std::cout << "[INFO] Successfully opened file: " << inFile << "\n";


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


  Plot2DBlockEtaPhi(hUnc3D, hCor3D, isFirstPass, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/2DPlots");

  FitLocalPhiEta(hUnc3D,           // uncorrected 3-D histogram
                   hCor3D,           // corrected 3-D histogram (may be nullptr)
                   isFirstPass,      // same toggle you already use
                   eEdges,           // vector with the eight E-ranges
                    "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits",           // where PNGs will be written
                   bOut);            // SAME open ofstream instance

  OverlayUncorrPhiEta(hUnc3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
  PlotPhiShiftAndWidth(hUnc3D, hCor3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
  OverlayDeltaPhiSlices(eEdges, binMode, isFirstPass, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits/deltaPhiFromScratch");
    
  /* ------------------------------------------------------------------ */
  /*  NEW: clusteriser overlays  (3‑curve)                              */
  /* ------------------------------------------------------------------ */
  OverlayDeltaPhiClusterizerCP(eEdges, binMode, isFirstPass,
        "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits/deltaPhiClusterizerOnly");

  /* ------------------------------------------------------------------ */
  /*  NEW: five‑way master overlay  (bespoke ±b, CP, CP+CP‑b)           */
  /* ------------------------------------------------------------------ */
  OverlayDeltaPhiFiveWays(eEdges, binMode, isFirstPass,
        "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits/deltaPhiClusterizerVsFromScratch");
    
  OverlayDeltaEtaSlices( eEdges,          // vector< pair<double,double> >
                           binMode,         // EBinningMode  (kRange | kDiscrete)
                           isFirstPass,     // bool – true on pass-1
       "/Users/patsfan753/Desktop/PositionDependentCorrection/"
       "SimOutput/1DplotsAndFits/deltaEtaFromScratch" );
    
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


// makeLegoGifHD(hCor3D, "cor", "CORRECTED",   out2DDir);
// makeLegoGifHD(hUnc3D, "unc", "UNCORRECTED", out2DDir);


  // after opening the file …
  TH1F* h_truth_vz = static_cast<TH1F*>( fIn->Get("h_truth_vz") );

  if(!h_truth_vz){
      std::cerr << "[ERROR] vertex-Z histograms not found!\n";
  } else {
      PlotVertexZTruthOnly(h_truth_vz, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/QA");
  }

    
    // Retrieve the two maps
    auto bRMSMap  = plotAshLogRMS_sideBySide(inFile);                // RMS-optimised Ash-b
    auto bPhiMap  = PlotBvaluesVsEnergy(hUnc3D, eEdges,
                        "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");

    // Overlay & export
    PlotBcompare(bRMSMap, bPhiMap,
                 "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput");

    
  std::string chi2Dir = std::string(outDir) + "/QA/Chi2_QA";
  PlotChi2QA(inFile, chi2Dir.c_str());
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
