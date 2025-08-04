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
    TH1F frU("frU",";E_{ctr}  [GeV];#mu  [rad]",1,xAxisMin,xAxisMax);
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
    TH1F frL("frL",";E_{ctr}  [GeV];#sigma  [rad]",1,xAxisMin,xAxisMax);
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

    const std::string hTitle = ";E_{ctr} [GeV];" + yTitle;
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
                Form("%s :  #mu = %.3g #pm %.3g   #sigma = %.3g #pm %.3g",
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
             "#sqrt{#mu^{2}+#sigma^{2}}  [rad]", eCtr, RMSE, variantsToPlot);

    drawDiag(outDir + "/" + fileStem + "_FracRes.png",
             "#sigma / |#mu|", eCtr, Frac, variantsToPlot);

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
             "#sigma_{variant}/#sigma_{baseline}", eCtr, WR, variantsToPlot);

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

            TH1F frU("frU",";E_{ctr}  [GeV];#mu  [rad]",1,xMin,xMax);
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

            TH1F frL("frL",";E_{ctr}  [GeV];#sigma  [rad]",1,xMin,xMax);
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
                    ? "Coresoftware Code w/ CorrectPosition  #mu/#sigma Summary"
                : (!strcmp(kLab[v],"b(E) corr, scratch"))
                    ? "From Scratch Code w/ Position Correction  #mu/#sigma Summary"
                    : std::string(kLab[v]) + "  #mu/#sigma Summary";
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
                TH1F frU("frU",";E_{ctr} [GeV];#mu [rad]",1,0.0,xMax);
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
                TH1F frL("frL",";E_{ctr} [GeV];#sigma [rad]",1,0.0,xMax);
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

                std::string tSum = tMain + "  Gaussian #mu/#sigma Summary";
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
                    TH1F frU("frU",";E_{ctr} [GeV];#mu [rad]",1,0.0,xMax);
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

                    TH1F frL("frL",";E_{ctr} [GeV];#sigma [rad]",1,0.0,xMax);
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
                        "#Delta#eta / #Delta#phi, From Scratch/Coresoftware Comparison Gaussian #mu/#sigma Summary");

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
                        TH1F frU("frU",";E_{ctr} [GeV];#mu [rad]",1,0.0,xMax);
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
                        TH1F frL("frL",";E_{ctr} [GeV];#sigma [rad]",1,0.0,xMax);
                        frL.SetMinimum(std::max(0.0,sgLo-0.05*sgHi)); frL.SetMaximum(1.15*sgHi); frL.Draw("AXIS");
                        drawGraphs(false);

                        /* title & save ----------------------------------------------------- */
                        c2.cd(0);
                        TLatex tit; tit.SetNDC(); tit.SetTextSize(0.032); tit.SetTextAlign(22);
                        tit.DrawLatex(0.5,0.97,
                            Form("#Delta#eta / #Delta#phi  %s Comparison Gaussian #mu/#sigma", D.tag));
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
//   per‑variant 2×4 tables that overlay the six |z_vtx| Δη slices
//      • one sub‑folder  “…/deltaEta/vertexDependent/<variant>/”
//      • file name       “DeltaEta_VtxSlices_<variant>.png”
// ─────────────────────────────────────────────────────────────────────────────
void makeEtaVertexTables(const std::vector<std::pair<double,double>>&               eEdges,
                         const std::array<std::vector<std::vector<TH1F*>>,5>&       hEtaVz,
                         const char*                                                outDir,     // “…/deltaEta/vertexDependent”
                         const std::vector<int>&                                    variantsToPlot = {0,1,2,3,4})
{
    /* canvas geometry */
    constexpr int nCol = 4, nRow = 2, maxPads = nCol * nRow;

    /* vertex‑bin edges */
    static constexpr std::array<float,7> vzEdge = {0,5,10,15,20,25,30};
    constexpr int N_Vz = vzEdge.size() - 1;

    /* colours for each vertex slice (all markers are style 20) */
    static constexpr std::array<Color_t,6> kVzCol =
        { kRed+1, kOrange+1, kYellow+2, kGreen+2, kAzure+2, kViolet+1 };
    static constexpr std::array<Style_t,6> kVzMk  = { 20,20,20,20,20,20 };

    /* iterate over reconstruction variants */
    for (int v : variantsToPlot)
    {
        /* sanitise variant name for file‑system usage */
        std::string vName = kLab[v];
        for (char& c : vName)
            if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';

        const std::string vDir = std::string(outDir) + "/" + vName;
        ensureDir(vDir);

        TCanvas cTab(Form("cEtaVzTab_%d", v), kLab[v], 1600, 900);
        cTab.SetTopMargin(0.12);                 // breathing room for title
        cTab.Divide(nCol, nRow, 0, 0);

        int pad = 1;
        for (std::size_t iE = 0; iE < eEdges.size() && pad <= maxPads; ++iE)
        {
            /* skip if no vertex histograms exist for this energy bin */
            bool haveSlice = false;
            for (int iVz = 0; iVz < N_Vz && !haveSlice; ++iVz)
                if (iVz < (int)hEtaVz[v][iE].size() && hEtaVz[v][iE][iVz])
                    haveSlice = true;
            if (!haveSlice) continue;

            int curPad = pad++;          // remember pad index before increment
            cTab.cd(curPad);

            /* containers for later summary table */
            std::array<FitRes,6> fitInfo{};   // one per vz slice
            std::array<bool,   6> sliceUsed{}; sliceUsed.fill(false);

            double yMax = 0.0;

            /* legend position: push right by pad’s left margin */
            const double xLeg1 = gPad->GetLeftMargin() + 0.05;
            const double xLeg2 = xLeg1 + 0.26;
            TLegend* leg = new TLegend(xLeg1, 0.73, xLeg2, 0.93);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);

            /* overlay all vertex slices for this energy bin */
            for (int iVz = 0; iVz < N_Vz; ++iVz)
            {
                if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                TH1F* src = hEtaVz[v][iE][iVz];
                if (!src || src->Integral() == 0) continue;

                auto* h = cloneAndNormPdf(src, v);       // colour replaced below
                h->SetMarkerColor(kVzCol[iVz]);
                h->SetLineColor  (kVzCol[iVz]);
                h->SetMarkerStyle(kVzMk[iVz]);
                h->SetMarkerSize(0.9);

                yMax = std::max(yMax, h->GetMaximum() * 1.15);

                if (leg->GetNRows() == 0)
                {
                    h->SetTitle(Form("Δη  [%.0f, %.0f) GeV",
                                     eEdges[iE].first, eEdges[iE].second));
                    h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                }
                else
                    h->Draw("E SAME");

                /* robust Gaussian fit */
                const FitRes fr = robustGaussianFit(src);
                fitInfo[iVz]   = fr;
                sliceUsed[iVz] = true;

                leg->AddEntry(h,
                    Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                    "lep");
            }

            /* fix y‑range on all primitives */
            gPad->Update();
            for (auto* obj : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(obj))
                    h->GetYaxis()->SetRangeUser(0., yMax);

            leg->Draw();

            /* energy‑range label (top‑right) */
            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.037); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.88, 0.92,
                Form("[%.0f, %.0f) GeV", eEdges[iE].first, eEdges[iE].second));

            /* decide table Y‑start: raise in bottom row to avoid clipping */
            const int rowIdx = (curPad-1) / nCol;     // 0 = top, 1 = bottom
            double yRow = (rowIdx == 0 ? 0.40 : 0.50);

            /* μ | σ | z  table – now bottom‑left, smaller font, proper ROOT underline */
            const double xTbl1 = gPad->GetLeftMargin() + 0.04;   // μ column
            const double xTbl2 = xTbl1 + 0.12;                   // σ column
            const double xTbl3 = xTbl2 + 0.14;                   // z column

            TLatex tbl; tbl.SetNDC(); tbl.SetTextAlign(13);

            /* column headers – escape the blank in “z [cm]” with \,   */
            tbl.SetTextSize(0.032);
            /* underline works only if you **turn on TeX-style parsing first**   */
            tbl.SetTextFont(132);               // Helvetica-like, TeX parser ON

            tbl.DrawLatex(xTbl1, yRow, "#mu");
            tbl.DrawLatex(xTbl2, yRow, "#sigma");
            tbl.DrawLatex(xTbl3, yRow, "z\\,[cm]");   // escaped thin-space

            tbl.SetTextSize(0.026);           // body of the table a bit smaller
            yRow -= 0.050;                    // first data row

            for (int iVz = 0; iVz < N_Vz; ++iVz)
            {
                if (!sliceUsed[iVz]) continue;

                tbl.SetTextColor(kVzCol[iVz]);
                tbl.DrawLatex(xTbl1, yRow, Form("%.2g", fitInfo[iVz].mu));
                tbl.DrawLatex(xTbl2, yRow, Form("%.2g", fitInfo[iVz].sg));
                tbl.DrawLatex(xTbl3, yRow, Form("%.0f-%.0f", vzEdge[iVz], vzEdge[iVz+1])); // ASCII hyphen
                yRow -= 0.045;                // tight but safe vertical spacing
            }
            tbl.SetTextColor(kBlack);         // restore default colour

            // ───────────────────────────────────────────────────────────
            //  NEW — write this overlay (one energy‑bin) as its own PNG
            //      • file lives in the same <variant> directory
            //      • file name: “DeltaEta_VtxSlices_<variant>_E<slice>.png”
            // ───────────────────────────────────────────────────────────
            gPad->Modified();                       // finalise layout
            gPad->Update();
            gPad->SaveAs(Form("%s/DeltaEta_VtxSlices_%s_E%02zu.png",
                               vDir.c_str(),        // “…/deltaEta/vertexDependent/<variant>/”
                               vName.c_str(),       // sanitised variant label
                               iE));                // energy‑slice index
        }

        /* canvas-level title */
        cTab.cd(0);
        TLatex head; head.SetNDC(); head.SetTextAlign(22);

        /* make the font a bit smaller for the long Clusterizer label */
        head.SetTextSize(0.045);

        /* special wording for the “CorrectPosition__cluster” variant (v == 2) */
        const char* labTxt =
            (v == 2 ? "Clusterizer Code with Correction, Bins of z-vertex"
                    : kLab[v]);

        head.DrawLatex(0.5, 0.975, labTxt);

        /* write PNG for the full 6‑slice table */
        cTab.SaveAs(Form("%s/DeltaEta_VtxSlices_%s.png",
                         vDir.c_str(), vName.c_str()));
        cTab.Close();

        // ───────────────────────────────────────────────────────────────
        //  NEW ➊  ––  per‑vertex‑slice μ(E) & σ(E) summary
        //             • top pad  : μ vs E
        //             • bottom   : σ vs E
        //             • six curves = six |z_vtx| bins (kVzCol colours)
        //             • output    : “…/DeltaEta_Vtx_MuSigmaVsE_<variant>.png”
        // ───────────────────────────────────────────────────────────────
        {
            /* accumulate centre‑of‑energy and fit results */
            std::vector<double> eCtr;  eCtr.reserve(eEdges.size());
            std::array<std::vector<double>,6> MUvz, dMUvz, SGvz, dSGvz;

            for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
            {
                const double eMid = 0.5 * (eEdges[iE].first + eEdges[iE].second);
                eCtr.push_back(eMid);

                /* ensure all slices have a slot for this energy */
                for (int iVz = 0; iVz < N_Vz; ++iVz) {
                    MUvz[iVz].push_back(0.0);  dMUvz[iVz].push_back(0.0);
                    SGvz[iVz].push_back(0.0);  dSGvz[iVz].push_back(0.0);
                }

                if (iE >= hEtaVz[v].size()) continue;
                for (int iVz = 0; iVz < N_Vz; ++iVz)
                {
                    if (iVz >= (int)hEtaVz[v][iE].size()) continue;
                    TH1F* src = hEtaVz[v][iE][iVz];
                    if (!src || src->Integral() == 0)     continue;

                    const FitRes fr = robustGaussianFit(src);
                    MUvz[iVz].back()  = fr.mu;  dMUvz[iVz].back() = fr.dmu;
                    SGvz[iVz].back()  = fr.sg;  dSGvz[iVz].back() = fr.dsg;
                }
            }

            /* prepare canvas */
            TCanvas cMS(Form("cVtxMuSig_%d", v), "vtx_mu_sigma", 900, 800);
            cMS.Divide(1, 2, 0, 0);
            std::vector<double> ex(eCtr.size(), 0.0);

            // ---------- μ pad --------------------------------------------------
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

            TH1F frU("frU",";E_{ctr}  [GeV];#mu  [rad]",1,0.0,xMax);
            frU.SetMinimum(muMin - 0.05*(muMax-muMin));
            frU.SetMaximum(muMax + 0.05*(muMax-muMin));
            frU.Draw("AXIS");

            /* dashed reference line at μ = 0 */
            TLine ref0(0.0, 0.0, xMax, 0.0);
            ref0.SetLineStyle(2); ref0.Draw("same");

            /* draw curves */
            for (int iVz = 0; iVz < N_Vz; ++iVz) {
                auto* g = new TGraphErrors((int)eCtr.size(),
                                           eCtr.data(), MUvz[iVz].data(),
                                           ex.data(),  dMUvz[iVz].data());
                g->SetMarkerColor(kVzCol[iVz]);
                g->SetLineColor  (kVzCol[iVz]);
                g->SetMarkerStyle(kVzMk[iVz]);
                g->SetMarkerSize(1.1);
                g->Draw("P SAME");
            }

            /* legend – top‑right of upper pad only */
            TLegend legU(0.70, 0.6, 0.92, 0.9);
            legU.SetBorderSize(0);
            legU.SetFillStyle(0);
            legU.SetTextSize(0.035);
            for (int iVz = 0; iVz < N_Vz; ++iVz) {
                auto* m = new TMarker(0,0,kVzMk[iVz]); m->SetMarkerColor(kVzCol[iVz]);
                legU.AddEntry(m,Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),"p");
            }
            legU.Draw();

            // ---------- σ pad --------------------------------------------------
            cMS.cd(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.06);
            gPad->SetTopMargin(0.06);
            gPad->SetBottomMargin(0.35);

            double sgMax = -1e30;
            for (int iVz = 0; iVz < N_Vz; ++iVz)
                for (double vsg : SGvz[iVz]) sgMax = std::max(sgMax, vsg);

            TH1F frL("frL",";E_{ctr}  [GeV];#sigma  [rad]",1,0.0,xMax);
            frL.SetMinimum(0.0);                     // start exactly at σ = 0
            frL.SetMaximum(sgMax * 1.05);            // 5 % head‑room
            frL.Draw("AXIS");

            for (int iVz = 0; iVz < N_Vz; ++iVz) {
                auto* g = new TGraphErrors((int)eCtr.size(),
                                           eCtr.data(), SGvz[iVz].data(),
                                           ex.data(),  dSGvz[iVz].data());
                g->SetMarkerColor(kVzCol[iVz]);
                g->SetLineColor  (kVzCol[iVz]);
                g->SetMarkerStyle(kVzMk[iVz]);
                g->SetMarkerSize(1.1);
                g->Draw("P SAME");
            }

            /* save canvas */
            cMS.SaveAs(Form("%s/DeltaEta_Vtx_MuSigmaVsE_%s.png",
                            vDir.c_str(), vName.c_str()));
            cMS.Close();
        }
        /* ════════════════════════════════════════════════════════════════════════
         *  NEW ➋ – low / mid / high |z_vtx| overlays
         *         • folder   “…/<variant>/threeSlices/”
         *         • 2×4 table for every energy bin (3 slices only)
         *         • μ(E) & σ(E) summary for the same three slices
         * ══════════════════════════════════════════════════════════════════════ */
        {
            const std::array<int,3> pick3 = {0, N_Vz/2, N_Vz-1};      // low–mid–high
            std::string vDir3 = vDir + "/threeSlices";
            ensureDir(vDir3);

            /* ────────────────────────────────────────────────────────────────────────────
             *  per-energy-bin OVERLAYS – now two **very-large** 1 × 3 canvases
             *      • “Elo” = first three energy-bins   (bins 0-1-2)
             *      • “Ehi” = last three energy-bins    (bins N-3, N-2, N-1)
             *      • canvas geometry  ➔  3200 × 1200 px  (plenty of vertical space)
             *      • every pad shows:
             *            – three |z_vtx| slices (low/mid/high) in colour
             *            – legend with the  |z_vtx|  ranges **and** the exact E-bin
             * ────────────────────────────────────────────────────────────────────────── */
            auto drawThreeTable =
            [&](const char* tag, std::size_t startIdx)
            {
                TCanvas cT3(Form("cEtaVz3_%s_%d",tag,v),"eta_vtx_three",3200,1200);
                cT3.SetTopMargin   (0.07);   // more head-room for big title
                cT3.SetBottomMargin(0.06);
                cT3.Divide(3,1,0.002,0.0);   // 1 row × 3 cols, 0.2 % gutters

                for (int pad=1; pad<=3; ++pad)
                {
                    const std::size_t iE = startIdx + pad-1;
                    if (iE >= eEdges.size()) continue;

                    /* verify that at least one of the 3 slices exists */
                    bool have=false;
                    for (int iVz: pick3)
                        if (iVz < (int)hEtaVz[v][iE].size() &&
                            hEtaVz[v][iE][iVz] &&
                            hEtaVz[v][iE][iVz]->Integral()>0){ have=true; break; }
                    if (!have) continue;

                    cT3.cd(pad);
                    gPad->SetLeftMargin(0.11);
                    gPad->SetRightMargin(0.03);
                    gPad->SetTopMargin(0.08);
                    gPad->SetBottomMargin(0.12);

                    double yMax = 0.0;

                    auto* leg = new TLegend(0.13, 0.73, 0.50, 0.93);      // ← upper‑left corner
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.042);                              // same font size
                    leg->AddEntry(static_cast<TObject*>(nullptr),         // energy label only
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
                        h->SetMarkerColor(kVzCol[iVz]);
                        h->SetLineColor  (kVzCol[iVz]);
                        h->SetMarkerStyle(kVzMk[iVz]);
                        h->SetMarkerSize(1.2);            // slightly bigger
                        yMax = std::max(yMax, h->GetMaximum()*1.20);

                        if (!axesDrawn) {
                            h->SetTitle("");

                            /* main-axis titles */
                            h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                            h->GetYaxis()->SetTitle("Probability density");

                            /* shrink fonts & pull them closer to the axes */
                            h->GetXaxis()->SetTitleSize(0.045);   // smaller axis-title font
                            h->GetXaxis()->SetLabelSize(0.035);   // smaller numeric labels
                            h->GetXaxis()->SetTitleOffset(1.10);  // move title closer to axis ticks

                            h->GetYaxis()->SetTitleSize(0.045);
                            h->GetYaxis()->SetLabelSize(0.035);
                            h->GetYaxis()->SetTitleOffset(1.40);  // shift Y-title toward axis

                            h->Draw("E");                         // first slice draws axes
                            axesDrawn = true;
                        } else {
                            h->Draw("E SAME");
                        }


                        leg->AddEntry(h,
                                      Form("|z_{vtx}|: [%.0f, %.0f) cm",
                                           vzEdge[iVz], vzEdge[iVz+1]),
                                      "p");
                    }

                    /* uniform Y-range for all histograms in this pad */
                    for (auto* obj : *gPad->GetListOfPrimitives())
                        if (auto* h = dynamic_cast<TH1*>(obj))
                            h->GetYaxis()->SetRangeUser(0., yMax);

                    leg->Draw();
                    gPad->Modified();   // make sure legend becomes part of pad
                    gPad->Update();
                }

                /* big canvas-level heading */
                cT3.cd(0);
                TLatex hd; hd.SetNDC(); hd.SetTextAlign(22); hd.SetTextSize(0.05);
                hd.DrawLatex(0.5,0.96,
                    (std::string(kLab[v]) +
                     "  , low / mid / high |z_{vtx}|  (" + tag + ')').c_str());

                cT3.Modified();   // finalise layout, then write PNG
                cT3.Update();
                cT3.SaveAs(Form("%s/DeltaEta_Vtx3Slices_%s_%s.png",
                                vDir3.c_str(), vName.c_str(), tag));
                cT3.Close();
            };

            /* build “Elo”  (bins 0-1-2)  and  “Ehi”  (bins N-3, N-2, N-1) */
            if (eEdges.size() >= 1) drawThreeTable("Elo", 0);
            if (eEdges.size() >= 4) drawThreeTable("Ehi",
                                                   eEdges.size() >= 3 ? eEdges.size()-3 : 0);


            /* ---------- μ(E) & σ(E) summary for the same three slices -------- */
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

            /* μ pad */
            cMS3.cd(1);
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.04);
            double muLo=1e30,muHi=-1e30;
            for(int slot=0;slot<3;++slot)
                for(std::size_t k=0;k<MU3[slot].size();++k){
                    muLo=std::min(muLo,MU3[slot][k]-dMU3[slot][k]);
                    muHi=std::max(muHi,MU3[slot][k]+dMU3[slot][k]);
                }
            TH1F frU3("frU3",";E_{ctr}  [GeV];#mu  [rad]",1,0.0,xMax3);
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
                g->SetMarkerColor(kVzCol[iVz]); g->SetLineColor(kVzCol[iVz]);
                g->SetMarkerStyle(kVzMk[iVz]);  g->SetMarkerSize(1.1);
                g->Draw("P SAME");
                lu3.AddEntry(g,Form("[%.0f, %.0f) cm",
                                    vzEdge[iVz], vzEdge[iVz+1]),"p");
            }
            lu3.Draw();

            /* σ pad */
            cMS3.cd(2);
            gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);
            double sgHi=-1e30;
            for(int slot=0;slot<3;++slot)
                for(double vsg:SG3[slot]) sgHi=std::max(sgHi,vsg);
            TH1F frL3("frL3",";E_{ctr}  [GeV];#sigma  [rad]",1,0.0,xMax3);
            frL3.SetMinimum(0.0); frL3.SetMaximum(1.15*sgHi); frL3.Draw("AXIS");

            for(int slot=0;slot<3;++slot){
                int iVz=pick3[slot];
                TGraphErrors* g=new TGraphErrors(eCtr3.size(), eCtr3.data(),
                                                 SG3[slot].data(),
                                                 ex3.data(),  dSG3[slot].data());
                g->SetMarkerColor(kVzCol[iVz]); g->SetLineColor(kVzCol[iVz]);
                g->SetMarkerStyle(kVzMk[iVz]);  g->SetMarkerSize(1.1);
                g->Draw("P SAME");
            }

            /* title & save */
            cMS3.cd(0);
            TLatex ttl; ttl.SetNDC(); ttl.SetTextAlign(22); ttl.SetTextSize(0.04);
            ttl.DrawLatex(0.5,0.97,
                (std::string(kLab[v])+"  #mu/#sigma vs E  (low/mid/high |z_{vtx}|)").c_str());
            cMS3.SaveAs(Form("%s/DeltaEta_Vtx3_MuSigmaVsE_%s.png",
                             vDir3.c_str(), vName.c_str()));
            cMS3.Close();
        }
        /* ─── keep the original brace that ends the per‑variant loop ─── */
    }   // ————— end of per-variant loop ———————————————————————————


    /* ═══════════════════════════════════════════════════════════════════════════
     *  NEW ➌  –  From‑Scratch vs Clusterizer comparisons   (complete rewrite)
     *           • four overlay scenarios
     *                ① scratch‑RAW      vs cluster‑RAW          (0 ↔ 2)
     *                ② scratch‑CORR     vs cluster‑CORR         (1 ↔ 3)
     *                ③ scratch  RAW ↔ CORR                     (0 ↔ 1)
     *                ④ cluster  RAW ↔ CORR                     (2 ↔ 3)
     *           • output directory hierarchy
     *                   …/fromScratchVsClusterizerCompare/
     *                       vz_0_5/          ← per‑vertex bin
     *                           EtaCmp_<tag>_Table.png
     *                           EtaCmp_<tag>_MuSigma.png
     *                       vz_5_10/
     *                       …
     *           • marker convention
     *                   variant A  → filled circle  (kMarker=20)
     *                   variant B  → open   circle  (kMarker=24)
     * ═════════════════════════════════════════════════════════════════════════*/
    {
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
        ensureDir(cmpRoot);          // top‑level folder

        /* geometry reused from the main routine --------------------------------*/
        constexpr int nCol = 4, nRow = 2, maxPads = nCol * nRow;

        /* loop over vertex‑slice index (six bins → six sub‑folders) -------------*/
        for (int iVz = 0; iVz < N_Vz; ++iVz)
        {
            const std::string vzDir = Form("%s/vz_%g_%g",
                                           cmpRoot.c_str(),
                                           vzEdge[iVz], vzEdge[iVz+1]);
            ensureDir(vzDir);

            /* ————— one comparison canvas per CompSpec ————————————————*/
            for (const auto& C : comps)
            {
                /* ①  2×4 residual‑shape table ----------------------------------*/
                {
                    TCanvas c(Form("cEtaCmp_%s_vz%d",C.tag,iVz),          // name
                              "",                                         // ← **empty** title
                              1600,900);
                    c.SetTopMargin(0.12);
                    c.Divide(nCol, nRow, 0, 0);

                    int pad = 1;
                    for (std::size_t iE = 0;
                         iE < eEdges.size() && pad <= maxPads; ++iE)
                    {
                        /* skip if *both* variants miss this (E, z) histogram */
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

                        /* two variants, A first (filled), B second (open) */
                        for (int which = 0; which < 2; ++which)
                        {
                            const int v = (which==0 ? C.vA : C.vB);

                            /* safety: histogram may be absent */
                            if (   iE >= hEtaVz[v].size()
                                || iVz >= (int)hEtaVz[v][iE].size()
                                || !hEtaVz[v][iE][iVz]
                                || hEtaVz[v][iE][iVz]->Integral()==0 )
                                continue;

                            TH1F* src = hEtaVz[v][iE][iVz];
                            auto* h = cloneAndNormPdf(src, v);
                            /* --- unified colour / marker convention ----------------------------------
                             * scratch : v0 RAW  -> open  red  square   (mk 25)
                             *           v1 CORR -> solid red  square   (mk 21)
                             * cluster : v2 RAW  -> open  blue circle   (mk 24)
                             *           v3 CORR -> solid blue circle   (mk 20)
                             * ------------------------------------------------------------------------ */
                            Style_t mk   = 20;          // default (filled circle)
                            Color_t col  = kBlack;      // fallback

                            switch (v) {
                                case 0:  mk = 25; col = kRed;       break;   // scratch-RAW
                                case 1:  mk = 21; col = kRed;       break;   // scratch-CORR
                                case 2:  mk = 24; col = kBlue+1;    break;   // cluster-RAW
                                case 3:  mk = 20; col = kBlue+1;    break;   // cluster-CORR
                            }

                            h->SetMarkerColor(col);
                            h->SetLineColor  (col);
                            h->SetMarkerStyle(mk);
                            h->SetMarkerSize(1.1);          // slightly larger marker


                            yMax = std::max(yMax, h->GetMaximum()*1.15);

                            if (which==0)
                            {
                                h->SetTitle("");
                                h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                                h->GetYaxis()->SetTitle("Probability density");
                                h->Draw("E");
                            }
                            else  h->Draw("E SAME");

                            lg->AddEntry(h, kLab[v], "lep");
                        }

                        /* common Y‑range */
                        for (auto* obj : *gPad->GetListOfPrimitives())
                            if (auto* hh = dynamic_cast<TH1*>(obj))
                                hh->GetYaxis()->SetRangeUser(0., yMax);

                        /* vertex‑range label (top‑right) */
                        TLatex vzlab; vzlab.SetNDC();
                        vzlab.SetTextSize(0.037); vzlab.SetTextAlign(33);
                        vzlab.DrawLatex(0.88,0.92,
                            Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]));

                        lg->Draw();
                    }

                    c.SaveAs(Form("%s/EtaCmp_%s_Table.png",
                                  vzDir.c_str(), C.tag));
                    c.Close();
                }

                /* ②  μ(E) / σ(E) summary for this vertex‑slice ----------------*/
                {
                    /* accumulate μ,σ for the chosen vertex‑slice only */
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

                    TCanvas c(Form("cEtaCmpMS_%s_vz%d",C.tag,iVz),
                              "",                                         // ← **empty** title
                              900,800);
                    c.Divide(1,2,0,0);

                    /* —— μ pad —— */
                    c.cd(1);
                    gPad->SetLeftMargin(0.15);
                    gPad->SetBottomMargin(0.06);

                    /* determine Y‑range and force it to straddle 0 so the baseline is visible */
                    double muLo =  1e30, muHi = -1e30;
                    for (int k = 0; k < 2; ++k)
                        for (std::size_t j = 0; j < MU[k].size(); ++j) {
                            muLo = std::min(muLo, MU[k][j] - dMU[k][j]);
                            muHi = std::max(muHi, MU[k][j] + dMU[k][j]);
                        }
                    const double yMin = std::min(0.0, muLo - 0.10*(muHi - muLo));   // always ≤ 0
                    const double yMax = muHi + 0.10*(muHi - muLo);                  // head‑room

                    TH1F frU("frU","",1,0.0,xMax);              // absolutely empty title string
                    frU.GetXaxis()->SetTitle("E_{ctr} [GeV]");
                    frU.GetYaxis()->SetTitle("#mu [rad]");
                    frU.SetMinimum(yMin);
                    frU.SetMaximum(yMax);
                    frU.Draw("AXIS");

                    /* solid, dashed black baseline drawn *after* the axes so it is visible */
                    TLine ref0(0.0, 0.0, xMax, 0.0);
                    ref0.SetLineColor(kBlack);
                    ref0.SetLineStyle(2);                       // ROOT style 2 = dashed
                    ref0.Draw("same");

                    TLegend legU(0.70,0.750,0.92,0.92);
                    legU.SetBorderSize(0); legU.SetFillStyle(0);
                    legU.SetTextSize(0.04);

                    for(int k=0; k<2; ++k){
                        const int v = (k==0 ? C.vA : C.vB);
                        /* --- apply the same unified style table as above ------------------------ */
                        Style_t mk   = 20;
                        Color_t col  = kBlack;

                        switch (v) {
                            case 0:  mk = 25; col = kRed;       break;   // scratch-RAW
                            case 1:  mk = 21; col = kRed;       break;   // scratch-CORR
                            case 2:  mk = 24; col = kBlue+1;    break;   // cluster-RAW
                            case 3:  mk = 20; col = kBlue+1;    break;   // cluster-CORR
                        }

                        TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                           eCtr.data(), MU[k].data(),
                                                           ex.data(),  dMU[k].data());
                        g->SetMarkerColor(col);  g->SetLineColor(col);
                        g->SetMarkerStyle(mk);   g->SetMarkerSize(1.1);
                        g->Draw("P SAME");

                        /* legend entry – marker only, no error bar */
                        auto* m = new TMarker(0,0,mk);
                        m->SetMarkerColor(col);  m->SetMarkerSize(1.1);
                        legU.AddEntry(m, kLab[v], "p");

                    }

                    legU.Draw();

                    /* —— σ pad —— */
                    c.cd(2);
                    gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.06);
                    double sgHi=-1e30;
                    for(int k=0;k<2;++k)
                        for(double s:SG[k]) sgHi=std::max(sgHi,s);
                    TH1F frL("frL",";E_{ctr} [GeV];#sigma [rad]",1,0.0,xMax);
                    frL.SetMinimum(0.0); frL.SetMaximum(1.15*sgHi);
                    frL.Draw("AXIS");
                    for(int k=0; k<2; ++k){
                        const int  v  = (k==0 ? C.vA : C.vB);

                        /* unified colour / marker convention */
                        Style_t mk   = 20;
                        Color_t col  = kBlack;
                        switch (v) {
                            case 0:  mk = 25; col = kRed;     break;   // scratch-RAW  (open red square)
                            case 1:  mk = 21; col = kRed;     break;   // scratch-CORR (filled red square)
                            case 2:  mk = 24; col = kBlue+1;  break;   // cluster-RAW  (open blue circle)
                            case 3:  mk = 20; col = kBlue+1;  break;   // cluster-CORR (filled blue circle)
                        }

                        TGraphErrors* g = new TGraphErrors((int)eCtr.size(),
                                                           eCtr.data(), SG[k].data(),
                                                           ex.data(),  dSG[k].data());
                        g->SetMarkerColor(col);  g->SetLineColor(col);
                        g->SetMarkerStyle(mk);   g->SetMarkerSize(1.1);
                        g->Draw("P SAME");
                    }

                    /* ---------- finally write the PNG before closing ---------- */
                    c.SaveAs(Form("%s/EtaCmp_%s_MuSigma.png",
                                  vzDir.c_str(), C.tag));
                    c.Close();
                }
            }   /* —— end of CompSpec loop —— */
        }       /* —— end of vertex‑bin loop —— */
    }           /* —— end of NEW ➌ block —— */





    // ─────────────────────────────────────────────────────────────────
    //  “condensed” summary – one canvas **per energy slice**
    //      • layout: 4 rows × 2 columns  (max 8 variants shown)
    //      • overlays **only** the first (iVz = 0) and last
    //        (iVz = N_Vz – 1) |z_vtx| bins
    //      • output dir:  “…/vertexDependent/vertexCondensedSummary/”
    // ─────────────────────────────────────────────────────────────────
    const std::string summaryDir =
        std::string(outDir) + "/vertexCondensedSummary";
    ensureDir(summaryDir);

    constexpr int nColC   = 2, nRowC = 3, maxPadsC = nColC * nRowC;
    const std::array<int,2> pickVz = {0, N_Vz-1};          // low & high slices

    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        /* roomier square canvas – 2400 × 2400 px */
        TCanvas cSum(Form("cEtaCond_%zu", iE), "condensed_eta_vtx",
                     2400, 2400);

        /* slim global margins – pads fill almost the entire canvas */
        cSum.SetTopMargin   (0.06);
        cSum.SetBottomMargin(0.04);
        cSum.SetLeftMargin  (0.02);
        cSum.SetRightMargin (0.02);

        /* divide into 3 rows × 2 columns, tiny 0.5 % gutters */
        cSum.Divide(nColC, nRowC, 0.005, 0.005);

        /* give every pad sane internal margins */
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
            if (padC > maxPadsC) break;            // only 8 pads available
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

                auto* h = cloneAndNormPdf(src, v);          // recolour
                h->SetMarkerColor(kCol[v]);                 // colour = variant
                h->SetLineColor  (kCol[v]);
                h->SetMarkerStyle( k==0 ? 20   /* filled */   // low‑|z_vtx|
                                    : 24 ); /* open */        // high‑|z_vtx|
                h->SetMarkerSize(1.0);

                ymax = std::max(ymax, h->GetMaximum()*1.15);

                if (k==0) {
                    h->SetTitle("");                        // custom labels below
                    h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                } else
                    h->Draw("E SAME");

                fr[k] = robustGaussianFit(src);
                l->AddEntry(h,
                            Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                            "p");
            }

            gPad->Update();
            for (auto* o : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(o))
                    h->GetYaxis()->SetRangeUser(0., ymax);

            l->Draw();

            // variant label (colour‑coded)
            TLatex vtxt; vtxt.SetNDC(); vtxt.SetTextColor(kCol[v]);
            vtxt.SetTextSize(0.045); vtxt.SetTextAlign(11);
            vtxt.DrawLatex(gPad->GetLeftMargin()+0.04, 0.85, kLab[v]);

            // energy‑range label (upper‑right)
            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.045); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.88, 0.87,
                Form("[%.0f, %.0f) GeV", eEdges[iE].first, eEdges[iE].second));

            // μ / σ table (bottom‑left)
            TLatex t; t.SetNDC(); t.SetTextFont(132); t.SetTextAlign(13);
            t.SetTextSize(0.035);

            const double x0 = gPad->GetLeftMargin() + 0.04;
            double y0 = 0.30;
            for (int k = 0; k < 2; ++k)
            {
                int iVz = pickVz[k];
                t.SetTextColor(kVzCol[iVz]);
                t.DrawLatex(x0, y0,
                    Form("#mu=%.2g  #sigma=%.2g  z:%g-%g",
                         fr[k].mu, fr[k].sg,
                         vzEdge[iVz], vzEdge[iVz+1]));
                y0 -= 0.05;
            }
            t.SetTextColor(kBlack);
        }

        cSum.SaveAs(Form("%s/DeltaEta_Condensed_E%02zu.png",
                         summaryDir.c_str(), iE));
        cSum.Close();
    }
    
    // ─────────────────────────────────────────────────────────────────
    //  global matrix summary
    //      rows    → flavour variants      (up to 5, one per row)
    //      columns → the five lowest E‑bins
    //      every cell: overlay low‑|z_vtx| (filled)  +
    //                              high‑|z_vtx| (open)
    // ─────────────────────────────────────────────────────────────────
    /* show only the 3 lowest‑E slices → wider pads */
    const std::size_t nEPlot = std::min<std::size_t>(3, eEdges.size());
    const int nColM = static_cast<int>(nEPlot);                  // 3 columns
    const int nRowM = static_cast<int>(variantsToPlot.size());   // 5 rows

    const std::string mDir =
        std::string(outDir) + "/vertexCondensedSummary";
    ensureDir(mDir);

    /* wide landscape canvas – 4500 px × 3000 px */
    TCanvas cMat("cEtaVzMatrix", "eta_vtx_matrix", 4500, 3000);
    cMat.SetTopMargin   (0.03);
    cMat.SetBottomMargin(0.03);
    cMat.SetLeftMargin  (0.02);
    cMat.SetRightMargin (0.02);

    /* 3 cols × 5 rows, 0.8 % gutters */
    cMat.Divide(nColM, nRowM, 0.008, 0.008);

    /* iterate:   outer → variant row,   inner → energy‑bin column */
    int padMat = 1;
    for (int r = 0; r < nRowM; ++r)
    {
        const int v = variantsToPlot[r];

        for (std::size_t c = 0; c < nEPlot; ++c, ++padMat)
        {
            cMat.cd(padMat);

            /* skip pad completely if no histograms exist */
            if (c >= hEtaVz[v].size() ||
                hEtaVz[v][c].empty()     ||
                (!hEtaVz[v][c][0] && !hEtaVz[v][c].back()) )
            {
                continue;
            }

            double ymax = 0.0;
            std::array<FitRes,2> fr{};

            /* pad‑local legend (bottom‑right of each cell) */
            auto* leg = new TLegend(0.7, 0.22, 0.88, 0.42);
            leg->SetBorderSize(0);  leg->SetFillStyle(0);  leg->SetTextSize(0.055);

            /* loop over the two chosen |z_vtx| slices: low (filled) & high (open) */
            for (int k = 0; k < 2; ++k)
            {
                const int iVz = (k == 0 ? 0 : N_Vz-1);
                if (iVz >= (int)hEtaVz[v][c].size()) continue;
                TH1F* src = hEtaVz[v][c][iVz];
                if (!src || src->Integral() == 0)  continue;

                auto* h = cloneAndNormPdf(src, v);
                h->SetMarkerColor(kCol[v]);   h->SetLineColor(kCol[v]);
                h->SetMarkerStyle( k==0 ? 20 : 24 );   // filled / open
                h->SetMarkerSize(1.2);

                ymax = std::max(ymax, h->GetMaximum()*1.15);

                if (k == 0) {
                    /* first slice draws axes */
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle("#Delta#eta  [rad]");
                    h->GetYaxis()->SetTitle("Probability density");
                    h->Draw("E");
                } else
                    h->Draw("E SAME");

                fr[k] = robustGaussianFit(src);
                leg->AddEntry(h,
                              Form("[%.0f, %.0f) cm", vzEdge[iVz], vzEdge[iVz+1]),
                              "p");
            }

            /* harmonise Y‑range */
            gPad->Update();
            for (auto* o : *gPad->GetListOfPrimitives())
                if (auto* h = dynamic_cast<TH1*>(o))
                    h->GetYaxis()->SetRangeUser(0., ymax);

            leg->Draw();

            /* variant label – colour‑coded, left‑top inside pad */
            TLatex vLab; vLab.SetNDC(); vLab.SetTextColor(kCol[v]);
            vLab.SetTextSize(0.055); vLab.SetTextAlign(11);
            vLab.DrawLatex(gPad->GetLeftMargin()+0.03, 0.86, kLab[v]);

            /* energy‑bin label – right‑top */
            TLatex eLab; eLab.SetNDC(); eLab.SetTextSize(0.055); eLab.SetTextAlign(33);
            eLab.DrawLatex(0.91, 0.88,
                Form("[%.0f, %.0f) GeV", eEdges[c].first, eEdges[c].second));

            /* μ / σ table – bottom‑left */
            TLatex tbl; tbl.SetNDC(); tbl.SetTextFont(132); tbl.SetTextAlign(13);
            tbl.SetTextSize(0.052);

            double ytbl = 0.33;
            for (int k = 0; k < 2; ++k)
            {
                int iVz = (k == 0 ? 0 : N_Vz-1);
                tbl.SetTextColor(k == 0 ? kCol[v] : kCol[v]);
                tbl.DrawLatex(gPad->GetLeftMargin()+0.04, ytbl,
                    Form("#mu=%.2g  #sigma=%.2g  z:%g-%g",
                         fr[k].mu, fr[k].sg,
                         vzEdge[iVz], vzEdge[iVz+1]));
                ytbl -= 0.06;
            }
            tbl.SetTextColor(kBlack);
        }
    }

    /* save matrix PNG */
    cMat.SaveAs(Form("%s/DeltaEta_Vtx_Matrix.png", mDir.c_str()));
    cMat.Close();
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

    /* vertex‑resolved Δη vectors ------------------------------------------- */
    static constexpr std::array<float,7> vzEdge = {0,5,10,15,20,25,30};
    constexpr int N_Vz = vzEdge.size() - 1;          // = 6 bins

    std::array<std::vector<std::vector<TH1F*>>,5> etaVz;

    /* –– allocate –– */
    for (int v = 0; v < 5; ++v)
    {
        etaVz[v].resize(eEdges.size());              // energy dimension
        for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
            etaVz[v][iE].resize(N_Vz, nullptr);      // vz‑dimension
    }

    /* –– fill –– */
    for (std::size_t iE = 0; iE < eEdges.size(); ++iE)
    {
        const TString eTag = makeSliceTag(iE, eEdges[iE]);   // e.g. “2_4”

        for (int iVz = 0; iVz < N_Vz; ++iVz)
        {
            const TString vzTag = Form("vz%d_%d",
                                       int(vzEdge[iVz]), int(vzEdge[iVz+1])); // e.g. “vz5_10”

            /* variant‑selector ------------------------------------------------ */
            etaVz[0][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_raw_%s_%s"    , eTag.Data(), vzTag.Data())) );
            etaVz[1][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_corr_%s_%s"   , eTag.Data(), vzTag.Data())) );
            etaVz[2][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpRaw_%s_%s"  , eTag.Data(), vzTag.Data())) );
            etaVz[3][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpCorr_%s_%s" , eTag.Data(), vzTag.Data())) );
            etaVz[4][iE][iVz] = static_cast<TH1F*>( fIn->Get(
                                    Form("h_eta_diff_cpBcorr_%s_%s", eTag.Data(), vzTag.Data())) );
        }
    }

    /* ----- produce all plots ----------------------------------------------------- */
    /* ‘binMode’ (or whatever local/argument name you use) is already in scope */
    MakeDeltaPhiPlots(eEdges, binMode, phiView,
                      "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");

    MakeDeltaEtaPlots(eEdges, binMode,
                      etaGlobal, etaVz,
                      "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");

    
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
//
//

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
