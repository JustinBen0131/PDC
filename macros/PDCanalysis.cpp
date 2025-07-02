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


constexpr double kMaxAbsMu   =  0.05;   // [rad]  forbid |μ| above this
constexpr double kMinSigma   =  1e-4;   // [rad]  floor for σ
constexpr int    kMaxRetries =     6;   // extra passes if initial fit fails

/* --------------------------------------------------------------- *
 *  doFit  – one Gaussian fit with parameter limits + smart seeds  *
 * --------------------------------------------------------------- */
static bool
doFit(TH1* h, TF1& f, double winLo, double winHi,
      double& outMu, double& outMuErr,
      double& outSi, double& outSiErr, double& outChi2)
{
    /* enforce physical ranges on μ and σ ------------------------- */
    f.SetParLimits(1, -kMaxAbsMu,  kMaxAbsMu);
    f.SetParLimits(2,  kMinSigma,  1.0);

    /* four diverse seeds around the IQR centre ------------------- */
    const double A = h->GetMaximum();
    const double S = std::max(kMinSigma, 0.5*h->GetRMS());
    const double M = 0.5*(winLo + winHi);

    const double s[4][3] = { {A,     M      , S   },
                             {1.5*A, M      , 1.5*S},
                             {0.7*A, M+0.3*S, 0.8*S},
                             {1.2*A, M-0.3*S, 1.2*S} };

    double bestChi2 = 1e99;

    for (auto& p : s)
    {
        f.SetParameters(p[0], p[1], p[2]);
        TFitResultPtr r = h->Fit(&f, "QNR0S");
        if (!r.Get() || !r->IsValid()) continue;

        if (r->Chi2() < bestChi2) {
            bestChi2   = r->Chi2();
            outMu      = f.GetParameter(1);
            outMuErr   = f.GetParError (1);
            outSi      = f.GetParameter(2);
            outSiErr   = f.GetParError (2);
            outChi2 = r->Chi2() /
                      std::max<double>(1.0, static_cast<double>(r->Ndf()));
        }
    }
    return bestChi2 < 1e98;      // success flag
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

    /* --------------------------------------------------------------- *
     *  robustGauss – iterative, window‑shrinking Gaussian fit         *
     * --------------------------------------------------------------- */
    auto robustGauss = [](TH1* h) -> GRes
    {
        if (!h || h->Integral() == 0) return {0,0,0,0,0};

        /* IQR‑based start window ------------------------------------ */
        const double q[3] = {0.25, 0.50, 0.75};
        double quart[3];  h->GetQuantiles(3, quart, const_cast<double*>(q));
        double winLo = quart[0] - 2.5*(quart[2]-quart[0]);
        double winHi = quart[2] + 2.5*(quart[2]-quart[0]);

        TF1 fG("fG", "gaus", winLo, winHi);

        GRes g     {h->GetEntries(),0,0,0,0};
        double chi2 = 1e99;

        for (int pass = 0; pass <= kMaxRetries; ++pass)
        {
            if (doFit(h,fG,winLo,winHi,
                      g.mu,g.muErr,g.sg,g.sgErr,chi2))
                break;                          // good fit

            /* shrink window by 30 % on each side and retry ------------ */
            const double span = winHi-winLo;
            winLo += 0.30*span;
            winHi -= 0.30*span;
            if (winHi - winLo < 6.*kMinSigma) break;   // window too tight
            fG.SetRange(winLo, winHi);
        }

        /* ultimate fall‑back: robust RMS ----------------------------- */
        if (!std::isfinite(g.sg) || g.sg < kMinSigma) {
            g.mu     = h->GetMean();
            g.sg     = h->GetRMS();
            g.muErr  = g.sg / std::sqrt(2.*std::max(1., g.N-1.));
            g.sgErr  = g.muErr;
        }
        return g;
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
        
        /* ---------------------------------------------------------- *
         *  Diagnostic & safety net for the corrected fit             *
         * ---------------------------------------------------------- */
        {
            /* -------- ANSI helpers ---------------------------------- */
            const char *BOLD = "\033[1m", *RED  = "\033[1;31m",
                       *YEL  = "\033[1;33m", *GRN  = "\033[1;32m",
                       *RST  = "\033[0m";

            auto printRes = [&](const char* lbl, const GRes& g,
                                Color_t termCol = kBlack)
            {
                const char* clr = (termCol==kRed  ? RED :
                                   termCol==kGreen? GRN : "");
                std::cout << Form("    %s%-4s%s  N=%6.0f  μ=%+9.6f ±%8.6f"
                                  "   σ=%8.6f ±%8.6f\n",
                                  clr,lbl,RST,
                                  g.N, g.mu, g.muErr, g.sg, g.sgErr);
            };

            /* always print RAW result */
            printRes("RAW", R, kBlack);

            if (hCorr) {            /* only if corrected histogram exists */
                /* -------- sanity indicators -------------------------- */
                /* -------- sanity check : keep the fit unless it is truly invalid ---- */
                const bool badSigma   = (C.sg < kMinSigma)               || !std::isfinite(C.sg);
                const bool badRelErr  = (C.sgErr/C.sg > 1.0);            // >100 % relative σ‐err
                const bool badFit     = badSigma || badRelErr;

                printRes("CORR", C, badFit ? kRed : kMagenta+1);

                if (badFit)
                {
                    std::cerr << "\033[1;31m[Δφ] slice " << tag
                              << " – fit rejected (σ=" << C.sg
                              << ", σ_err/σ=" << C.sgErr/C.sg << ")\033[0m\n";

                    /* robust fall‑back: mean / RMS ------------------------------------ */
                    C.mu     = hCorr->GetMean();
                    C.sg     = hCorr->GetRMS();
                    C.muErr  = C.sg / std::sqrt(2.*std::max(1.0,C.N-1.0));
                    C.sgErr  = C.muErr;

                    printRes("CORR*", C, kGreen);
                }
            }
        }


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

        hRaw ->SetMarkerStyle(20); hRaw ->SetMarkerColor(kGreen + 2); hRaw ->SetLineColor(kGreen + 1);
        if (hCorr) {
            hCorr->SetMarkerStyle(20); hCorr->SetMarkerColor(kMagenta + 1); hCorr->SetLineColor(kMagenta + 2);
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

        double x  = 0.91;
        double dy = 0.04;          // vertical offset between lines

        // ---------- RAW (green) -------------------------------------------
        {
           TLatex tx; tx.SetNDC();
           tx.SetTextFont(42);
           tx.SetTextSize(0.036);
           tx.SetTextColor(kGreen+2);
           tx.SetTextAlign(33);

           double y0 = 0.88;       // shift block down to avoid overlap
           tx.DrawLatex(x, y0,        Form("#bf{#mu_{RAW} = %.4f}",   R.mu));
           tx.DrawLatex(x, y0 - dy,   Form("#bf{#sigma_{RAW} = %.4f}",R.sg));
        }

        
        // ---------- CORR (magenta) ----------------------------------------
        if (hCorr) {
           TLatex tx; tx.SetNDC();
           tx.SetTextFont(42);
           tx.SetTextSize(0.036);
           tx.SetTextColor(kMagenta+1);
           tx.SetTextAlign(33);

           double y0 = 0.8;       // first line (µ)
           tx.DrawLatex(x, y0,        Form("#bf{#mu_{CORR} = %.4f}",  C.mu));
           tx.DrawLatex(x, y0 - dy,   Form("#bf{#sigma_{CORR} = %.4f}",C.sg));
        }
        
        {
            /* e.g.  outDir/DeltaPhiSlice_8_12.png  or  .../DeltaPhiSlice_E9.png */
            TString slicePNG = TString(outDir) +
                               Form("/DeltaPhiSlice_%s.png", tag.Data());

            /* make sure everything is painted before printing           */
            pad->Modified();
            pad->Update();

            /* ROOT 5/6‑safe: write exactly this pad to file             */
            pad->Print(slicePNG, "png");

            std::cout << "           wrote slice → " << slicePNG << '\n';
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
  auto gMuRaw = makeG(muRaw, muRawErr, -dx, kGreen + 2, 20);
  auto gMuCor = makeG(muCor, muCorErr, +dx, kMagenta + 1  , 20);
  auto gSiRaw = makeG(sgRaw, sgRawErr, -dx, kGreen + 2, 20);
  auto gSiCor = makeG(sgCor, sgCorErr, +dx, kMagenta + 1, 20);

  keep.insert(keep.end(),{gMuRaw,gMuCor,gSiRaw,gSiCor});

    /* -- summary canvas (µ top, σ bottom) --------------------------------- */
    {
        /* ----- X‑range : pad 1 GeV on the left, go to upper edge of last slice --- */
        const double xMin = eEdges.front().first  - 1.0;
        const double xMax = eEdges.back ().second;

        TCanvas cSum("cMuSi","#Delta#phi mean / sigma vs E",860,760);
        TPad pT("pT","",0,0.34,1,1);  pT.Draw();
        TPad pB("pB","",0,0   ,1,0.34); pB.Draw();

        /* handy helpers to get global min / max including error bars ------------- */
        auto fullMin = [](const std::vector<double>& v,const std::vector<double>& e,
                          const std::vector<double>& w,const std::vector<double>& f)
        {
            double m =  1e30;
            for (size_t i=0;i<v.size();++i) m = std::min(m, v[i]-(i<e.size()?e[i]:0.));
            for (size_t i=0;i<w.size();++i) m = std::min(m, w[i]-(i<f.size()?f[i]:0.));
            return m;
        };
        auto fullMax = [](const std::vector<double>& v,const std::vector<double>& e,
                          const std::vector<double>& w,const std::vector<double>& f)
        {
            double M = -1e30;
            for (size_t i=0;i<v.size();++i) M = std::max(M, v[i]+(i<e.size()?e[i]:0.));
            for (size_t i=0;i<w.size();++i) M = std::max(M, w[i]+(i<f.size()?f[i]:0.));
            return M;
        };

        /* ---------------- μ(E) panel ------------------------------------------ */
        const double muMin = fullMin(muRaw,muRawErr,muCor,muCorErr);
        const double muMax = fullMax(muRaw,muRawErr,muCor,muCorErr);

        pT.cd();
        pT.SetLeftMargin(0.15);
        pT.SetBottomMargin(0.03);          // hairline gap

        TH1F fMu("fMu","; ;#mu_{Gauss}  [rad]",1,xMin,xMax);
        fMu.SetMinimum(muMin - 0.10*std::fabs(muMin));
        fMu.SetMaximum(muMax + 0.10*std::fabs(muMax));

        /* hide x–axis on upper pad -------------------------------------------- */
        fMu.GetXaxis()->SetLabelSize(0);
        fMu.GetXaxis()->SetTickLength(0);
        fMu.GetXaxis()->SetTitleSize(0);

        fMu.Draw("AXIS");
        gMuRaw->Draw("P SAME");
        gMuCor->Draw("P SAME");

        TLine *l0 = new TLine(xMin,0.0,xMax,0.0);
        l0->SetLineStyle(2); l0->SetLineWidth(2); l0->SetLineColor(kGray+2);
        l0->Draw();   keep.push_back(l0);

        TLegend legMu(0.18,0.79,0.43,0.91);
        legMu.SetBorderSize(0); legMu.SetFillStyle(0); legMu.SetTextSize(0.035);
        legMu.AddEntry(gMuRaw,"RAW","lp");
        legMu.AddEntry(gMuCor,"CORRECTED","lp");
        legMu.Draw();

        /* ---------------- σ(E) panel ----------------------------------------- */
        const double siMin = fullMin(sgRaw,sgRawErr,sgCor,sgCorErr);
        const double siMax = fullMax(sgRaw,sgRawErr,sgCor,sgCorErr);

        pB.cd();
        pB.SetLeftMargin(0.15);
        pB.SetTopMargin(0.06);
        pB.SetBottomMargin(0.38);

        TH1F fSi("fSi",";E_{slice centre}  [GeV];#sigma_{Gauss}  [rad]",1,xMin,xMax);
        fSi.SetMinimum(std::max(0.0, siMin - 0.10*std::fabs(siMin)));
        fSi.SetMaximum(siMax + 0.10*std::fabs(siMax));

        fSi.GetXaxis()->SetNdivisions(505);
        fSi.GetXaxis()->SetTickLength(0.05);
        fSi.GetXaxis()->SetTitleOffset(1.1);

        fSi.Draw("AXIS");
        gSiRaw->Draw("P SAME");
        gSiCor->Draw("P SAME");

        cSum.SaveAs(TString(outDir)+"/DeltaPhi_MeanSigmaVsE.png");
    }

    /* ------------------------------------------------------------------ *
     * 5)   EXTRA diagnostics  (quadrature RMSE and fractional resolution)
     * ------------------------------------------------------------------ */
    {
        /* helper: graph builder with horizontal offset ---------------------- */
        auto gShift = [&](const std::vector<double>& y,
                          const std::vector<double>& yErr,
                          double dxOff, Color_t col) -> TGraphErrors*
        {
            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x(N), ex(N,0.);
            for (int i=0;i<N;++i) x[i] = eCtr[i] + dxOff;

            auto *g = new TGraphErrors(N,
                                       x.data(),  y.data(),
                                       ex.data(), yErr.data());
            g->SetMarkerStyle(20);
            g->SetMarkerSize (0.9);
            g->SetMarkerColor(col);
            g->SetLineColor  (col);
            keep.push_back(g);
            return g;
        };

        /* generic drawer (used twice) --------------------------------------- */
        auto makeDiag =
            [&](const char* cname, const char* yTit,
                const std::vector<double>& yR,const std::vector<double>& yRE,
                const std::vector<double>& yC,const std::vector<double>& yCE,
                const char* png)
        {
            const double dxOff = 0.12;

            double yMin =  1e30, yMax = -1e30;
            for (size_t i=0;i<yR.size();++i){
                yMin = std::min({yMin,yR[i]-yRE[i],yC[i]-yCE[i]});
                yMax = std::max({yMax,yR[i]+yRE[i],yC[i]+yCE[i]});
            }
            const double span = yMax-yMin;
            yMin -= 0.02*span;   yMax += 0.05*span;

            TCanvas c(cname,cname,860,620);
            c.SetLeftMargin(0.15); c.SetRightMargin(0.06);

            TH1F frame("f",
                       Form(";E_{slice}  [GeV];%s",yTit),
                       1, eCtr.front()-1.0, eEdges.back().second);
            frame.SetMinimum(yMin); frame.SetMaximum(yMax);
            frame.Draw("AXIS");

            auto gRaw = gShift(yR,yRE,-dxOff,kGreen+1);
            auto gCor = gShift(yC,yCE,+dxOff,kMagenta+2);
            gRaw->Draw("P SAME"); gCor->Draw("P SAME");

            const bool isFrac = std::strstr(png,"FracRes")!=nullptr;
            double lx1=isFrac?0.70:0.15, ly1=isFrac?0.72:0.12,
                   lx2=isFrac?0.90:0.45, ly2=isFrac?0.92:0.32;

            TLegend leg(lx1,ly1,lx2,ly2);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.026);
            leg.AddEntry(gRaw,"RAW","lp");
            leg.AddEntry(gCor,"CORRECTED","lp");
            leg.Draw();

            c.SaveAs(TString(outDir)+"/"+png);
        };

        /* (i) RMSE = sqrt(mu²+sigma²) --------------------------------------- */
        std::vector<double> rmseR(nPts), rmseRE(nPts),
                            rmseC(nPts), rmseCE(nPts);

        for (int i=0;i<nPts;++i) {
            rmseR[i]  = std::hypot(muRaw[i],sgRaw[i]);
            rmseRE[i] = rmseR[i]*std::sqrt( std::pow(muRawErr[i]/muRaw[i],2) +
                                            std::pow(sgRawErr[i]/sgRaw[i],2) );

            rmseC[i]  = std::hypot(muCor[i],sgCor[i]);
            rmseCE[i] = (rmseC[i]>0)
                        ? rmseC[i]*std::sqrt( std::pow(muCorErr[i]/muCor[i],2) +
                                               std::pow(sgCorErr[i]/sgCor[i],2) )
                        : 0.0;
        }
        makeDiag("cRMSE",
                 "#sqrt{#mu^{2} + #sigma^{2}}  [rad]",
                 rmseR,rmseRE, rmseC,rmseCE,
                 "DeltaPhi_RMSErrorVsE.png");

        /* ================================================================== *
         *  Fractional resolution  σ / |μ|                                    *
         *  – only plot slices with |μ| > 3·σ_μ (otherwise μ ≈ 0)             *
         * ================================================================== */
        std::vector<double> fracR(nPts,0), fracRE(nPts,0),
                            fracC(nPts,0), fracCE(nPts,0);

        for (int i=0;i<nPts;++i)
        {
            /* RAW ----------------------------------------------------------- */
            if (std::fabs(muRaw[i]) > 3*muRawErr[i]) {
                fracR[i]  = sgRaw[i] / std::fabs(muRaw[i]);
                fracRE[i] = fracR[i] *
                            std::sqrt( std::pow(sgRawErr[i]/sgRaw[i],2) +
                                       std::pow(muRawErr[i]/muRaw[i],2) );
            }
            /* CORRECTED ----------------------------------------------------- */
            if (std::fabs(muCor[i]) > 3*muCorErr[i]) {
                fracC[i]  = sgCor[i] / std::fabs(muCor[i]);
                fracCE[i] = fracC[i] *
                            std::sqrt( std::pow(sgCorErr[i]/sgCor[i],2) +
                                       std::pow(muCorErr[i]/muCor[i],2) );
            }
        }

        makeDiag("cFrac",
                 "#sigma / |#mu|   (only if |#mu|>3#sigma_{#mu})",
                 fracR,fracRE, fracC,fracCE,
                 "DeltaPhi_FracResVsE.png");
    }
    /* ------------------------------------------------------------------ */
    /* ================================================================= *
     *  EXTRA summary – width ratio  R(E) = σ_CORR / σ_RAW               *
     *  ( R < 1  ⇒  improved resolution )                                *
     * ================================================================= */
    {
        std::vector<double> ratio(nPts,0.0), ratioErr(nPts,0.0);

        /* build R(E) and its uncertainty -------------------------------- */
        for (int i=0;i<nPts;++i) {
            if (sgRaw[i] > 0.0) {
                ratio[i]    = sgCor[i] / sgRaw[i];
                ratioErr[i] = ratio[i] *
                              std::sqrt( std::pow(sgCorErr[i]/sgCor[i],2) +
                                         std::pow(sgRawErr[i]/sgRaw[i],2) );
            }
        }

        /* determine y‑range with padding -------------------------------- */
        double yMin =  1e30, yMax = -1e30;
        for (int i=0;i<nPts;++i) {
            yMin = std::min(yMin, ratio[i]-ratioErr[i]);
            yMax = std::max(yMax, ratio[i]+ratioErr[i]);
        }
        yMin = std::max(0.0, yMin - 0.05);
        yMax = std::min(1.2,  yMax + 0.05);

        /* canvas & axis -------------------------------------------------- */
        TCanvas cRatio("cRatio","width ratio", 860, 620);
        cRatio.SetLeftMargin(0.15);
        cRatio.SetRightMargin(0.06);

        TH1F fR("fR",
                ";E_{slice}  [GeV];#sigma_{corr} / #sigma_{raw}",
                1, eCtr.front()-1.0, eEdges.back().second);
        fR.SetMinimum(yMin);
        fR.SetMaximum(yMax);
        fR.Draw("AXIS");

        /* graph with error bars ----------------------------------------- */
        auto gR = new TGraphErrors(nPts);
        for (int i=0;i<nPts;++i) {
            gR->SetPoint      (i, eCtr[i], ratio[i]);
            gR->SetPointError (i,     0.0, ratioErr[i]);
        }
        gR->SetMarkerStyle(20);
        gR->SetMarkerColor(kBlue+2);
        gR->SetLineColor  (kBlue+2);
        gR->Draw("PE SAME");              // P = points, E = error bars
        keep.push_back(gR);

        /* unity reference line ------------------------------------------ */
        TLine l1(fR.GetXaxis()->GetXmin(), 1.0,
                 fR.GetXaxis()->GetXmax(), 1.0);
        l1.SetLineStyle(2);
        l1.SetLineWidth(2);
        l1.SetLineColor(kGray+2);
        l1.Draw();

        /* legend (top‑right) -------------------------------------------- */
        TLegend leg(0.70, 0.80, 0.94, 0.92);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.032);
        leg.AddEntry(gR, "width ratio", "lp");
        leg.Draw();

        /* save ----------------------------------------------------------- */
        cRatio.SaveAs(TString(outDir)+"/DeltaPhi_WidthRatioVsE.png");
    }
    /* ================================================================= */


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
 * OverlayDeltaPhiFiveWays
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

  static const Color_t  colArr[5]   = {kGreen + 2, kMagenta + 1, kBlack, kRed, kBlue};
  static const Style_t  mksArr[5]   = {20,     20,   20,      20,        20};
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
    // ---------- cache parameters for later summary ----------
    static std::array<double,5> parA {{0.}},   // intercept  a
                                 parB {{0.}};  // slope      b
    static std::array<bool ,5>  hasFit{{false}};
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

    /* ------------------------------------------------------------------ *
     *  Δφ summary: Gaussian μ & σ vs energy, five reconstruction variants
     * ------------------------------------------------------------------ *
     *  This version adds a single “changePositions” switch:
     *    – changePositions = true  → legend moves to the top‑right of the μ‑pad
     *                                and extra head‑room is added so the legend
     *                                never hides points;
     *    – changePositions = false → original legend position and the y‑range is
     *                                left exactly as the automatic calculation
     *                                produced it.
     * ------------------------------------------------------------------ */

    /* ============================================================
     * user options
     * ============================================================
     */
    const bool   changePositions   = true;   // ← toggle legend & y‑range
    const double extraHeadRoomFrac = 0.30;   // add 25 % head‑room if TRUE
    const double manualMuMax       = 0.0;    // >0 ⇒ explicit upper limit

    /* ---------- 3. summary canvas μ & σ ---------- */
    {
        /* horizontal offset between the 5 TGraph points (-2 … +2) */
        const double dx = 0.0;

        /* ------------------------------------------------------------------ *
         *  Helper → return {min,max} over *all* curves, including ±σErr,
         *  then pad the range by 5 % so markers never touch the frame.
         * ------------------------------------------------------------------ */
        auto rangeOf = [&](bool wantMu) -> std::pair<double,double>
        {
            double lo =  1e30, hi = -1e30;
            for (int v = 0; v <= kCPb; ++v) {
                if (mu[v].empty()) continue;
                const auto& val =  wantMu ?  mu [v] :  sg [v];
                const auto& err =  wantMu ?  muE[v] :  sgE[v];
                for (size_t i = 0; i < val.size(); ++i) {
                    lo = std::min(lo, val[i] - err[i]);
                    hi = std::max(hi, val[i] + err[i]);
                }
            }
            const double pad = 0.05 * (hi - lo);
            return {lo - pad, hi + pad};
        };

        /* ------------------------------------------------------------------ *
         *  Canvas with two vertically stacked pads
         * ------------------------------------------------------------------ */
        TCanvas cSum("cMuSigma5", "#Delta#phi five‑way summary", 900, 780);
        TPad pT("pT","",0,0.37,1,1);  pT.Draw();
        TPad pB("pB","",0,0   ,1,0.37); pB.Draw();

        const double xMin =  eCtr.front() - 1.0;        // pad left by 1 GeV
        const double xMax =  eEdges.back().second;      // right edge is true E‑edge

        /* ====================  μ(E) top pad  ==================== */
        pT.cd();
        pT.SetLeftMargin(0.15);
        pT.SetBottomMargin(0.03);                       // hair‑line gap

        auto [muMin, muMax] = rangeOf(/*wantMu=*/true);

        /* optional extra head‑room when legend sits at the top */
        if (changePositions) {
            muMax += extraHeadRoomFrac * (muMax - muMin);
            if (manualMuMax > 0.0) muMax = manualMuMax;   // explicit override
        }

        TH1F fMu("fMu","; ;Gaussian mean  #mu  [rad]",1,xMin,xMax);
        fMu.SetMinimum(muMin);
        fMu.SetMaximum(muMax);

        /* hide x‑axis on the upper pad */
        TAxis* axTop = fMu.GetXaxis();
        axTop->SetLabelSize(0);
        axTop->SetTickLength(0);
        axTop->SetTitleSize(0);

        fMu.Draw("AXIS");

        /* dashed reference line μ = 0 */
        auto* l0 = new TLine(xMin, 0.0, xMax, 0.0);
        l0->SetLineStyle(2); l0->SetLineWidth(2); l0->SetLineColor(kGray+2);
        l0->Draw();             guard.push_back(l0);

        /* plot the five μ‑curves (if present) */
        for (int v = 0; v <= kCPb; ++v) {
            if (mu[v].empty()) continue;
            TGraphErrors* g = makeG(v, (v-2)*dx, colArr[v], 20);
            g->Draw("P SAME");
        }

        /* legend – coordinates depend on changePositions */
        double x1, y1, x2, y2;
        if (changePositions) {           // top‑right
            x1 = 0.38; y1 = 0.75; x2 = 0.85; y2 = 0.93;
        } else {                         // original position (bottom‑right)
            x1 = 0.40; y1 = 0.08; x2 = 0.82; y2 = 0.28;
        }

        TLegend lg(x1, y1, x2, y2);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.SetTextSize(0.032);
        for (int v = 0; v <= kCPb; ++v) {
            if (mu[v].empty()) continue;
            auto* m = new TMarker(0,0,20);           // invisible anchor
            m->SetMarkerColor(colArr[v]);
            m->SetMarkerSize(1.1);
            guard.push_back(m);
            lg.AddEntry(m, legTxtArr[v], "p");
        }
        lg.Draw();

        /* ====================  σ(E) bottom pad  ==================== */
        pB.cd();
        pB.SetLeftMargin(0.15);
        pB.SetTopMargin(0.07);
        pB.SetBottomMargin(0.35);

        auto [sgMin, sgMax] = rangeOf(/*wantMu=*/false);

        TH1F fSi("fSi",
                 ";Energy slice centre  [GeV];Gaussian width  #sigma  [rad]",
                 1, xMin, xMax);
        fSi.SetMinimum(std::max(0.0, sgMin));
        fSi.SetMaximum(sgMax);
        fSi.Draw("AXIS");

        for (int v = 0; v <= kCPb; ++v) {
            if (sg[v].empty()) continue;

            const int N = static_cast<int>(eCtr.size());
            std::vector<double> x(N), ex(N,0.);
            for (int i = 0; i < N; ++i) x[i] = eCtr[i] + (v-2)*dx;

            TGraphErrors* g = new TGraphErrors(
                      N, x.data(), sg[v].data(),
                      ex.data(),  sgE[v].data());

            g->SetMarkerStyle(20);
            g->SetMarkerSize (1.1);
            g->SetMarkerColor(colArr[v]);
            g->SetLineColor  (colArr[v]);
            guard.push_back(g);
            g->Draw("P SAME");
        }

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

    /* ================================================================= *
     *  width‑ratio summary – σ_variant / σ_baseline                      *
     *  baseline variant = kCP  ( no‑corr, clusterizer transforms )       *
     * ================================================================= */
    {
        const int base = kCP;                       // index of the reference
        const int N    = static_cast<int>(eCtr.size());

        /* build R(E) and its uncertainty for every variant -------------- */
        std::array<std::vector<double>,5> R , dR ;
        for (int v = 0; v <= kCPb; ++v) {
            R [v].resize(N, 0.0);
            dR[v].resize(N, 0.0);

            if (sg[v].empty()) continue;            // variant absent

            for (int i = 0; i < N; ++i) {
                const double si   = sg [v][i];
                const double dsi  = sgE[v][i];
                const double sRef = sg [base][i];
                const double dsRef= sgE[base][i];

                if (sRef <= 0.0) continue;          // avoid divide‑by‑zero

                R [v][i]  =  si / sRef;
                dR[v][i]  =  R[v][i] *
                             std::sqrt( std::pow(dsi /si ,2) +
                                        std::pow(dsRef/sRef,2) );
            }
        }

        /* auto y‑range --------------------------------------------------- */
        double yMin =  1e30, yMax = -1e30;
        for (int v=0; v<=kCPb; ++v)
            for (int i=0;i<N;++i){
                yMin = std::min(yMin, R[v][i]-dR[v][i]);
                yMax = std::max(yMax, R[v][i]+dR[v][i]);
            }
        yMin = std::max(0.0, yMin-0.05);
        yMax = std::min(1.5, yMax+0.05);            // clamp to 1.5

        /* canvas & frame ------------------------------------------------- */
        TCanvas cRat("cRat","width ratios", 900,640);
        cRat.SetLeftMargin(0.15);  cRat.SetRightMargin(0.06);

        TH1F frame("fRat",
                   ";Energy slice centre  [GeV];#sigma_{variant} / #sigma_{baseline}",
                   1, eCtr.front()-1.0,  eEdges.back().second);
        frame.SetMinimum(yMin);
        frame.SetMaximum(yMax);
        frame.Draw("AXIS");

        /* dashed unity line --------------------------------------------- */
        TLine l1(frame.GetXaxis()->GetXmin(), 1.0,
                 frame.GetXaxis()->GetXmax(), 1.0);
        l1.SetLineStyle(2); l1.SetLineWidth(2); l1.SetLineColor(kGray+2);
        l1.Draw();

        /* horizontal separation between variants ------------------------ */
        const double dxOff = 0.14;

        /* legend – top‑right -------------------------------------------- */
        TLegend lg(0.42,0.45,0.69,0.65);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.025);

        for (int v = 0; v <= kCPb; ++v)
        {
            if (v == base || sg[v].empty()) continue;       // skip baseline curve

            std::vector<double> x(N), ex(N,0.);
            for (int i=0;i<N;++i) x[i] = eCtr[i] + (v-base)*dxOff;

            TGraphErrors* g = new TGraphErrors(N,
                                               x.data(), R [v].data(),
                                               ex.data(), dR[v].data());
            g->SetMarkerStyle(20);
            g->SetMarkerColor(colArr[v]);
            g->SetLineColor  (colArr[v]);
            g->Draw("PE SAME");                             // with error bars
            guard.push_back(g);

            lg.AddEntry(g, legTxtArr[v], "lp");
        }
        lg.Draw();

        cRat.SaveAs(TString(outDir)+"/DeltaPhi5_WidthRatioVsE.png");
    }
    /* ================================================================= */


  /* ---------- 5. save slice sheet ---------- */
  cMain.SaveAs(TString(outDir)+"/DeltaPhi5_Compare_AllSlices.png");
  std::cout<<"[Δφ‑5] overlays & summaries written to "<<outDir<<'\n';
    
    /* ---------- 6.  log-E diagnostic : μ  versus  ln E  (five variants) ------ *
     *  • overlays μ(E) for all 5 reconstruction flavours in ln E space          *
     *  • fits   μ(E) = p0 + p1·lnE   (natural log, E in GeV)                    *
     *  • writes the coefficients to  …/DeltaPhi5_MuVsLogE_fit.txt               *
     * ------------------------------------------------------------------------ */
    {
        const int N = static_cast<int>(eCtr.size());
        if (N < 2){
            std::cerr << "[Δφ-5] log-E plot skipped (need ≥2 slices)\n";
            /* nothing else to do */
        } else {

            /* ------------- build ln E array once ----------------------------- */
            std::vector<double> lnE(N);
            for (int i = 0; i < N; ++i) lnE[i] = std::log(eCtr[i]);

            /* ------------- helper → graph builder --------------------------- */
            auto makeGraph = [&](int v)->TGraphErrors*{
                std::vector<double> ex(N,0.);
                auto *g = new TGraphErrors(
                            N, lnE.data(),   mu [v].data(),
                            ex.data(),       muE[v].data());
                g->SetMarkerStyle(mksArr[v]);
                g->SetMarkerColor(colArr[v]);
                g->SetLineColor  (colArr[v]);
                guard.push_back(g);
                return g;
            };

            std::array<TGraphErrors*,5> g{};
            for (int v = 0; v <= kCPb; ++v)
                if (!mu[v].empty()) g[v] = makeGraph(v);

            /* ------------- canvas & axes ------------------------------------ */
            TCanvas cLn("cMuVsLogE5",
                        "#Delta#phi Gaussian #mu  versus lnE  –  five reconstructions",
                        900,640);
            cLn.SetLeftMargin(0.15);  cLn.SetRightMargin(0.06);
            cLn.SetTopMargin (0.08);  cLn.SetBottomMargin(0.12);

            const double xLo = *std::min_element(lnE.begin(), lnE.end()) - 0.05;
            const double xHi = *std::max_element(lnE.begin(), lnE.end()) + 0.05;

            double yLo =  1e30 , yHi = -1e30;
            for (int v = 0; v <= kCPb; ++v)
                if (!mu[v].empty()){
                    yLo = std::min(yLo, *std::min_element(mu[v].begin(), mu[v].end()));
                    yHi = std::max(yHi, *std::max_element(mu[v].begin(), mu[v].end()));
                }
            const double pad = 0.15 * (yHi - yLo);     // +15 % extra head-room
            yLo -= pad;  yHi += 1.25*pad;              // more room on top

            TH1F frame("fLn5",";ln E  [GeV];Gaussian mean  #mu  [rad]",1,xLo,xHi);
            frame.SetMinimum(yLo);
            frame.SetMaximum(yHi);
            frame.Draw("AXIS");

            /* ------------- draw points & dashed linear fits ------------------ */
            std::array<TF1*,5> fit{};
            
            for (int v = 0; v <= kCPb; ++v)
                if (g[v]){
                    g[v]->Draw("P SAME");

                    fit[v]  = new TF1(Form("f%d",v),"pol1",xLo,xHi);
                    fit[v]->SetLineColor (colArr[v]);
                    fit[v]->SetLineStyle (2);          // thin dashed
                    fit[v]->SetLineWidth (2);
                    g[v]->Fit(fit[v],"Q");             // quiet
                    fit[v]->Draw("SAME");
                    parA[v]   = fit[v]->GetParameter(0);   // save intercept
                    parB[v]   = fit[v]->GetParameter(1);   // save slope
                    hasFit[v] = true;                      // mark this variant as fitted
                    guard.push_back(fit[v]);
                }
            
            /* ---------- derive one (a,b) pair for the tilt correction ---------- */
            double sumA = 0.0 , sumB = 0.0;     // running totals of p0 , p1
            int    nFit = 0;                    // number of successful fits

            for (int v = 0; v <= kCPb; ++v)
            {
              if (!fit[v]) continue;
              sumA += fit[v]->GetParameter(0);   // p0  from  μ(E) ≈ p0 + p1 lnE
              sumB += fit[v]->GetParameter(1);   // p1
              ++nFit;
            }

            if (nFit > 0)
            {
              /* average over all reconstruction variants */
              const double a_mean = sumA / nFit;          // still in  μ(E)  convention
              const double b_mean = sumB / nFit;

              /* convert to the φ = a – b lnE convention used in CorrectShowerDepth:        *
               * we need  φ(E) = –μ(E)  ⇒  a_tilt = –a_mean ,  b_tilt =  +b_mean            */
              const double a_tilt = -a_mean;              // positive ≈ +1.7 mrad
              const double b_tilt =  b_mean;              // ≈ +0.83 mrad

              std::cout << "\n[Δφ‑5]  recommended azimuth‑tilt constants  (detailed geometry)\n"
                        << "        a  = " << std::scientific << a_tilt << "  rad\n"
                        << "        b  = " << std::scientific << b_tilt << "  rad\n\n";
            }
            else
            {
              std::cerr << "[Δφ‑5]  WARNING: no valid fits – cannot derive tilt constants\n";
            }


            /* ------------- compact legend (bottom-right) --------------------- */
            TLegend lg(0.16, 0.75,   // x-min, y-min   (NDC)
                       0.90, 0.90);  // x-max, y-max
            lg.SetBorderSize(0);
            lg.SetFillStyle(0);
            lg.SetTextFont(42);
            lg.SetTextSize(0.02);
            lg.SetMargin(0.10);       // <<< shrink gap between marker and text
            

            /* add an entry only if both the graph and its fit exist – use a marker clone
             * instead of the TGraphErrors itself to keep the legend compact (no wings)   */
            for (int v = 0; v <= kCPb; ++v)
                if (g[v] && fit[v]) {
                    auto *mk = new TMarker(0, 0, mksArr[v]);   // marker only
                    mk->SetMarkerColor(colArr[v]);
                    mk->SetMarkerSize(1.2);                    // match the plot size
                    guard.push_back(mk);                       // survive until clean-up

                    lg.AddEntry(mk,
                                Form("%s  (a = %.2e,  b = %.2e)",
                                     legTxtArr[v],
                                     fit[v]->GetParameter(0),
                                     fit[v]->GetParameter(1)),
                                "p");                          // “p” → marker only
                }

            lg.Draw();
            
            /* ------------ explanatory label for the fit --------------------- */
            TLatex info;
            info.SetNDC();
            info.SetTextFont(42);
            info.SetTextSize(0.025);
            info.SetTextAlign(13);          // left-top
            /* place just below the legend – tweak y-coordinate if needed */
            info.DrawLatex(0.5, 0.18,
                           "#mu(E) = a - b lnE   (  a: intercept  |  b: slope )");

            /* ------------- save PNG & ASCII ---------------------------------- */
            TString pngName = TString(outDir) + "/DeltaPhi5_MuVsLogE.png";
            cLn.SaveAs(pngName);

            std::ofstream fout((std::string(outDir)+"/DeltaPhi5_MuVsLogE_fit.txt").c_str());
            fout << "# Linear fit  mu(E) = p0 + p1*lnE   (E in GeV)\n"
                 << "# columns: variant   p0(rad)   p1(rad)\n"
                 << std::scientific << std::setprecision(6);

            for (int v = 0; v <= kCPb; ++v)
                if (fit[v])
                    fout << std::left << std::setw(36) << legTxtArr[v]
                         << fit[v]->GetParameter(0) << "   "
                         << fit[v]->GetParameter(1) << "\n";
            fout.close();

            std::cout << "[Δφ-5] wrote cleaned-up lnE diagnostic plot → "
                      << pngName << '\n';
        }
    }
    
    /* ===================================================================== *
     *  EXTRA: summary of fit parameters  a  (intercept)  and  b  (slope)    *
     *         – one point per reconstruction flavour + Virgile’s reference  *
     *         – produces  .../DeltaPhi5_a_b_Summary.png                     *
     * ===================================================================== */
    try
    {
        std::cout << "[Δφ‑5:ab] building coefficient summary …\n";

        /* ---------- 1. gather coefficients --------------------------------- */
        struct Coeff {
            TString name;   double a;   double b;   Color_t col;  Style_t mk;
        };
        std::vector<Coeff> coeffs;

        /* publication‑friendly labels, same order as everywhere else -------- */
        static const char* shortLbl[5] = {
            "scratch, none", "scratch, b(E)",
            "cluster, none", "cluster, CP", "cluster, b(E)"
        };

        for (int v = 0; v <= kCPb; ++v)
            if (hasFit[v])               // keep only variants that were fitted
                coeffs.push_back({ shortLbl[v],
                                   parA[v], parB[v],
                                   colArr[v], mksArr[v] });

        /* add Virgile’s reference tilt (grey star) -------------------------- */
        coeffs.push_back({ "Virgile's Output",
                           +3.3e-3, +9.9e-4,
                           static_cast<Color_t>(kGray+2), 22 });

        if (coeffs.empty())
            throw std::runtime_error("[Δφ‑5:ab] No coefficients collected – cannot draw summary.");

        const int NV = static_cast<int>(coeffs.size());
        std::vector<double> x(NV);                     // x‑positions 1 … NV
        for (int i = 0; i < NV; ++i) x[i] = i + 1.0;

        /* ---------- helper: draw one pad (a or b) -------------------------- */
        auto drawPad = [&](TPad* p, const char* yLab, bool wantA)
        {
            if (!p) throw std::runtime_error("[Δφ‑5:ab] Null pad pointer.");

            /* ---- generic pad cosmetics ----------------------------------- */
            p->SetLeftMargin  (0.12);
            p->SetRightMargin (0.04);
            p->SetBottomMargin(0.16);
            p->SetTopMargin   (0.05);
            p->SetGrid(1,1);
            p->Draw();
            p->cd();

            /* ---- y‑range : global min/max over all points ---------------- */
            double ymin =  1e30, ymax = -1e30;
            for (const auto& c : coeffs) {
                const double y = wantA ? c.a : c.b;
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
            double pad = 0.15 * (ymax - ymin);
            if (pad == 0) pad = std::fabs(ymax) * 0.15 + 1.e-6;
            ymin -= pad;  ymax += pad;

            /* ---- axis frame (kept by the pad, no manual delete) ---------- */
            TH1F* fr = new TH1F(Form("fr_%s", yLab),
                                ";Variant index; ", NV+2, 0.5, NV + 1.5);
            fr->SetDirectory(nullptr);
            fr->SetMinimum(ymin);  fr->SetMaximum(ymax);
            fr->GetYaxis()->SetTitle(yLab);
            fr->GetXaxis()->SetLabelSize(0);                 // hide numbers
            fr->Draw("AXIS");

            /* dashed 0‑line for intercept panel ---------------------------- */
            if (wantA) {
                TLine* l0 = new TLine(fr->GetXaxis()->GetXmin(), 0.0,
                                      fr->GetXaxis()->GetXmax(), 0.0);
                l0->SetLineStyle(2); l0->SetLineWidth(2); l0->SetLineColor(kGray+2);
                l0->Draw();
            }

            /* ---- draw a one‑point TGraph per coefficient ----------------- */
            std::vector<TGraph*> keepAlive;      // local, *no* manual delete
            for (int i = 0; i < NV; ++i)
            {
                const double y = wantA ? coeffs[i].a : coeffs[i].b;
                double ex = 0., ey = 0.;

                auto* g = new TGraphErrors(1, &x[i], &y, &ex, &ey);
                g->SetMarkerStyle(coeffs[i].mk);
                g->SetMarkerColor(coeffs[i].col);
                g->SetMarkerSize (1.2);
                g->Draw("P SAME");

                keepAlive.push_back(g);          // survive until pad/canvas dies
            }

            /* ---- x‑axis labels under the frame --------------------------- */
            TLatex lx; lx.SetTextFont(42); lx.SetTextSize(0.035); lx.SetTextAlign(22);
            for (int i = 0; i < NV; ++i)
                lx.DrawLatex(x[i], ymin - 0.07*(ymax - ymin), coeffs[i].name);
        };

        /* ---------- build the two‑pad canvas ----------------------------- */
        TCanvas cAB("cAB", "Fit‑parameter summary", 900, 700);
        cAB.Divide(1, 2, 0, 0);

        drawPad( static_cast<TPad*>(cAB.cd(1)), "Intercept  a  [rad]", true  );
        drawPad( static_cast<TPad*>(cAB.cd(2)), "Slope  b  [rad]"    , false );

        /* ---------- save -------------------------------------------------- */
        TString outAB = TString(outDir) + "/DeltaPhi5_a_b_Summary.png";
        cAB.SaveAs(outAB);
        std::cout << "[Δφ‑5:ab] wrote coefficient summary  →  "
                  << outAB << '\n';
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n[Δφ‑5:ab] FATAL: " << e.what() << "\n";
        throw;                                // let caller decide what to do
    }



  for(auto* o:guard) delete o;
}

/**********************************************************************
 *  OverlayDeltaEtaFiveWays  –  robust & talkative version
 *  ------------------------------------------------------------------
 *  Author : ChatGPT rewrite, 2025‑06‑28
 *  Usage  : drop this function in **exactly** as‑is.  Before running
 *           you may tune the global   gOverlayEtaDbg   variable from
 *           the ROOT prompt, e.g.
 *
 *              root [0] gOverlayEtaDbg = 3;  // max verbosity
 *
 *           Verbosity levels
 *              0  silent except fatal errors
 *              1  entry / exit + file‑written messages          (default)
 *              2  per‑slice progress and main ranges
 *              3  full dump of μ/σ tables
 *********************************************************************/

namespace {
  /* ------------------------------------------------------------------ *
   *  DEBUG / ASSERT helpers                                            *
   * ------------------------------------------------------------------ */
  int gOverlayEtaDbg = 1;                     // run‑time adjustable

  #define DBG(lvl,msg)                                                             \
      do{ if(gOverlayEtaDbg >= lvl)                                                \
              std::cout << "[Δη‑DBG] " << msg << std::endl; }while(0)

  #define CHK(cond,msg)                                                            \
      do{ if(!(cond)){                                                             \
              std::ostringstream oss; oss << msg;                                  \
              throw std::runtime_error(oss.str()); } }while(0)
}

/* ====================================================================== *
 *  OverlayDeltaEtaFiveWays (robust variant)
 * ====================================================================== */
void OverlayDeltaEtaFiveWays(const std::vector<std::pair<double,double>>& eEdges,
                             EBinningMode  binMode,
                             bool          isFirstPass,
                             const char*   outDir = "plots/DeltaEta5")
{
  /* wrap the *entire* body so that any exception is caught once      */
  try{
    DBG(1,"ENTER  – isFirstPass="<<isFirstPass
           <<"  slices="<<eEdges.size()
           <<"  outDir="<<outDir);

    /* -------- 0. Book‑keeping / style ------------------------------ */
    gSystem->mkdir(outDir,true);
    gStyle->SetOptStat(0);

    static const Color_t  colArr[5] = {kGreen+2,kMagenta+1,kBlack,kRed,kBlue};
    static const Style_t  mksArr[5] = {20,20,20,20,20};
    static const char*    legTxt[5] = {
      "scratch, none","scratch, b(E)",
      "cluster, none","cluster, CP","cluster, b(E)"
    };

    const int nSlices = static_cast<int>(eEdges.size());
    CHK(nSlices>0,"eEdges is empty – nothing to do");

    auto tagOf = [&](int i)->std::string{
      return (binMode==EBinningMode::kRange)
           ? Form("%.0f_%.0f",eEdges[i].first,eEdges[i].second)
           : Form("E%.0f",eEdges[i].first);
    };

    /* ---------------------------------------------------------------- *
     *  tiny helper – robust Gaussian fit with IQR pre‑clipping         *
     * ---------------------------------------------------------------- */
    struct GRes { double mu,muE,sg,sgE; };
    auto robust = [&](TH1* h)->GRes
    {
      if(!h || h->Integral()==0) return {0,0,0,0};

      const double q[3]={0.25,0.50,0.75};  double quart[3];
      h->GetQuantiles(3,quart,(double*)q);
      double med=quart[1], iqr=quart[2]-quart[0];
      double lo = med-2.5*iqr , hi = med+2.5*iqr;
      if(lo>=hi){ lo=h->GetMean()-2*h->GetRMS(); hi=h->GetMean()+2*h->GetRMS(); }

      TF1 g("g","gaus",lo,hi);
      g.SetParameters(h->GetMaximum(),med,0.5*h->GetRMS());

      /* turn *off* stats box updates – they do heavy access in TList */
      g.SetParNames("A","mu","sigma");
      TFitResultPtr fr = h->Fit(&g,"QNR0S");

      if(!fr || !fr->IsValid()){
        const double rms=h->GetRMS();
        const double err=rms/std::sqrt(std::max(1.,h->GetEntries()-1.));
        return {h->GetMean(),err,rms,err};
      }
      return { g.GetParameter(1),g.GetParError(1),
               std::fabs(g.GetParameter(2)),g.GetParError(2)};
    };

    /* -------- 1. containers --------------------------------------- */
    enum {kRAW,kRAWb,kCP,kCPcp,kCPb};            // histogram flavours

    std::array<std::vector<double>,5> mu , muE , sg , sgE ;
    std::vector<double> eCtr;

    std::vector<std::unique_ptr<TH1>>  clones;   // we own *only* our clones
    std::vector<std::unique_ptr<TObject>>   keep;/* graphs, legends, lines */

    /* -------- 2. master canvas with overlays ---------------------- */
    const int nCols=4, nRows=(binMode==EBinningMode::kDiscrete?4:2);
    auto cMain = std::make_unique<TCanvas>("cDEta5","Δη five‑way",1800,600*nRows);
    cMain->Divide(nCols,nRows);

    for(int is=0; is<nSlices; ++is)
    {
      const std::string tag = tagOf(is);
      DBG(2,"… slice "<<is<<" / "<<nSlices<<"  tag="<<tag);

      /* locate histograms already sitting in memory ---------------- */
      TH1F* src[5] = {};
      src[kRAW ] = (TH1F*)gROOT->FindObject(Form("h_eta_diff_raw_%s",tag.c_str()));
      src[kRAWb] = (TH1F*)gROOT->FindObject(Form("h_eta_diff_corr_%s",tag.c_str()));

      if(!isFirstPass){
        src[kCP  ] = (TH1F*)gROOT->FindObject(Form("h_eta_diff_cpRaw_%s" ,tag.c_str()));
        src[kCPcp] = (TH1F*)gROOT->FindObject(Form("h_eta_diff_cpCorr_%s",tag.c_str()));
        src[kCPb ] = (TH1F*)gROOT->FindObject(Form("h_eta_diff_cpBcorr_%s",tag.c_str()));
      }

      CHK(src[kRAW],"Missing mandatory histogram h_eta_diff_raw_"<<tag);

      /* ---- clone + style + fit ---------------------------------- */
      GRes res[5]={{}};
      TH1* hPlot[5]={};

      for(int v=0; v<=kCPb; ++v)
        if(src[v]){
          auto c = std::unique_ptr<TH1F>((TH1F*)src[v]->Clone(
                                  Form("h_eta_%d_%d",v,is)));
          c->SetDirectory(nullptr);
          if(c->Integral()>0) c->Scale(1./c->Integral());

          c->SetMarkerStyle(mksArr[v]); c->SetMarkerSize(0.8);
          c->SetMarkerColor(colArr[v]); c->SetLineColor(colArr[v]);

          res[v]  = robust(c.get());
          hPlot[v]= c.get();
          clones.emplace_back(std::move(c));
        }

      eCtr.push_back(0.5*(eEdges[is].first+eEdges[is].second));
      for(int v=0; v<=kCPb; ++v){
        mu [v].push_back(res[v].mu );  muE[v].push_back(res[v].muE);
        sg [v].push_back(res[v].sg );  sgE[v].push_back(res[v].sgE);
      }

      /* ---- draw overlay ---------------------------------------- */
      TPad* pad = (TPad*)cMain->cd(is+1);
      CHK(pad,"null pad pointer");
      pad->SetLeftMargin(0.23); pad->SetBottomMargin(0.18);

      double yMax=0; for(int v=0; v<=kCPb; ++v) if(hPlot[v])
        yMax = std::max(yMax, hPlot[v]->GetMaximum());

      hPlot[kRAW]->SetTitle(Form("#Delta#eta  [%.0f,%.0f) GeV",
                                 eEdges[is].first,eEdges[is].second));
      hPlot[kRAW]->GetXaxis()->SetTitle("#eta_{reco} - #eta_{truth}");
      hPlot[kRAW]->GetYaxis()->SetTitle("probability density");
      hPlot[kRAW]->GetYaxis()->SetRangeUser(0,1.25*yMax);

      for(int v=0; v<=kCPb; ++v)
        if(hPlot[v]) hPlot[v]->Draw(v==kRAW? "E":"E SAME");

      auto lg = std::make_unique<TLegend>(0.24,0.66,0.78,0.90);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028);
      for(int v=0; v<=kCPb; ++v)
        if(hPlot[v]) lg->AddEntry(hPlot[v],legTxt[v],"lp");
      lg->Draw();  keep.emplace_back(std::move(lg));
    } /* slice loop */

    /* write the big sheet early so we know it succeeded ---------- */
    {
      TString fn = TString::Format("%s/DeltaEta5_Compare_AllSlices.png",outDir);
      cMain->SaveAs(fn);
      DBG(1,"Wrote "<<fn);
    }

    /* -------- 3. build summary tables --------------------------- */
    const int N = static_cast<int>(eCtr.size());
    CHK(N>0,"No valid slices accumulated – summary aborted");

    auto makeGraph =
      [&](const std::vector<double>& Y,const std::vector<double>& dY,int v)
          -> std::unique_ptr<TGraphErrors>
      {
        std::vector<double> X(N),eX(N,0.);
        for(int i=0;i<N;++i) X[i]=eCtr[i];
        auto g = std::make_unique<TGraphErrors>(N,X.data(),Y.data(),
                                                eX.data(),dY.data());
        g->SetMarkerStyle(mksArr[v]); g->SetMarkerSize(1.1);
        g->SetMarkerColor(colArr[v]); g->SetLineColor(colArr[v]);
        return g;
      };

    /* -------- 3a. μ / σ canvas ---------------------------------- */
    {
      TCanvas c("cMuSig","Δη μ/σ",900,780);

      /* --- top (μ) pad --- */
      TPad pT("pT","",0,0.37,1,1); pT.SetLeftMargin(0.15);
      pT.SetBottomMargin(0.03); pT.Draw(); pT.cd();

      double yLo=1e30,yHi=-1e30;
      for(int v=0;v<=kCPb;++v) for(int i=0;i<N;++i){
        yLo=std::min(yLo,mu[v][i]-muE[v][i]);
        yHi=std::max(yHi,mu[v][i]+muE[v][i]);
      }
      yLo-=0.05*(yHi-yLo); yHi+=0.25*(yHi-yLo);   // head‑room for legend

      TH1F frame("f","; ;Gaussian μ  (η‑units)",1,
                 eCtr.front()-1.0,eCtr.back()+1.0);
      frame.SetMinimum(yLo); frame.SetMaximum(yHi); frame.Draw("AXIS");

      TLine l0(frame.GetXaxis()->GetXmin(),0,
               frame.GetXaxis()->GetXmax(),0);
      l0.SetLineStyle(2); l0.SetLineColor(kGray+2); l0.Draw();

      auto lg = std::make_unique<TLegend>(0.38,0.78,0.85,0.93);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);

      for(int v=0;v<=kCPb;++v){
        auto g = makeGraph(mu[v],muE[v],v);
        if(g->GetN()==0) continue;
        g->Draw("P SAME");
        lg->AddEntry(g.get(),legTxt[v],"p");
        keep.emplace_back(std::move(g));
      }
      lg->Draw(); keep.emplace_back(std::move(lg));

      /* --- bottom (σ) pad --- */
      c.cd();
      TPad pB("pB","",0,0,1,0.37); pB.SetLeftMargin(0.15);
      pB.SetTopMargin(0.07); pB.SetBottomMargin(0.35);
      pB.Draw(); pB.cd();

      yLo=1e30; yHi=-1e30;
      for(int v=0;v<=kCPb;++v) for(int i=0;i<N;++i){
        yLo=std::min(yLo,sg[v][i]-sgE[v][i]);
        yHi=std::max(yHi,sg[v][i]+sgE[v][i]);
      }
      yLo = std::max(0.0,yLo-0.05*(yHi-yLo)); yHi+=0.05*(yHi-yLo);

      TH1F frame2("f2",";E slice centre [GeV];Gaussian σ  (η‑units)",1,
                  eCtr.front()-1.0,eCtr.back()+1.0);
      frame2.SetMinimum(yLo); frame2.SetMaximum(yHi); frame2.Draw("AXIS");

      for(int v=0;v<=kCPb;++v){
        auto g = makeGraph(sg[v],sgE[v],v);
        if(g->GetN()==0) continue;
        g->Draw("P SAME");
        keep.emplace_back(std::move(g));
      }

      TString fn = TString::Format("%s/DeltaEta5_MeanSigmaVsE.png",outDir);
      c.SaveAs(fn);
      DBG(1,"Wrote "<<fn);
    }

    /* -------- 3b. derived metrics (RMSE, frac‑res) -------------- */
    std::array<std::vector<double>,5> rm , rmE , fr , frE ;
    rm.fill({}); rmE.fill({}); fr.fill({}); frE.fill({});
    for(int v=0;v<=kCPb;++v){
      rm [v].resize(N); rmE[v].resize(N);
      fr [v].resize(N); frE[v].resize(N);

      for(int i=0;i<N;++i){
        const double m=mu[v][i], s=sg[v][i];
        const double dm=std::max(1e-12,muE[v][i]);
        const double ds=std::max(1e-12,sgE[v][i]);

        const double R  = std::hypot(m,s);
        const double dR = R? R*std::hypot(dm/m,ds/s):0.;

        const double F  = (std::fabs(m)>1e-12)? s/std::fabs(m):0.;
        const double dF = F? F*std::hypot(ds/s,dm/m):0.;

        rm [v][i]=R;  rmE[v][i]=dR;
        fr [v][i]=F;  frE[v][i]=dF;
      }
    }

    auto makeDiag =
      [&](const char* stem,const char* yTit,
          const std::array<std::vector<double>,5>& Y,
          const std::array<std::vector<double>,5>& dY)
    {
      TCanvas c(stem,stem,900,640); c.SetLeftMargin(0.15); c.SetRightMargin(0.06);

      double yLo=1e30,yHi=-1e30;
      for(int v=0;v<=kCPb;++v) for(int i=0;i<N;++i){
        yLo=std::min(yLo,Y[v][i]-dY[v][i]);
        yHi=std::max(yHi,Y[v][i]+dY[v][i]);
      }
      if(yLo>0) yLo=0; const double pad=0.07*(yHi-yLo); yLo-=pad; yHi+=pad;

      TH1F frame("f",Form(";E slice centre [GeV];%s",yTit),1,
                 eCtr.front()-1.0,eCtr.back()+1.0);
      frame.SetMinimum(yLo); frame.SetMaximum(yHi); frame.Draw("AXIS");

      auto lg = std::make_unique<TLegend>(0.15,0.75,0.55,0.93);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.024);

      for(int v=0;v<=kCPb;++v){
        auto g = makeGraph(Y[v],dY[v],v);
        if(g->GetN()==0) continue;
        g->Draw("P SAME");
        lg->AddEntry(g.get(),legTxt[v],"p");
        keep.emplace_back(std::move(g));
      }
      lg->Draw(); keep.emplace_back(std::move(lg));

      TString fn = TString::Format("%s/%s.png",outDir,stem);
      c.SaveAs(fn);
      DBG(1,"Wrote "<<fn);
    };

    makeDiag("DeltaEta5_RMSE"   ,"#sqrt{#mu^{2}+ #sigma^{2}}" ,rm ,rmE);
    makeDiag("DeltaEta5_FracRes","#sigma / |#mu|"             ,fr ,frE);

    DBG(1,"EXIT – OK");
  } /* try body */

  /* -------- global catch : print & re‑throw so ROOT prompt survives ---- */
  catch(const std::exception& e){
    std::cerr << "[Δη‑FATAL] "<<e.what()<<std::endl;
    throw;        // so you get the full cling back‑trace if desired
  }
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
    /* 2. 2‑D colour map (single histogram)                                   */
    /* ---------------------------------------------------------------------- */

    /* --- NEW: remember the maximum of the UNCORRECTED plot ---------------- */
    static double uncorrZmax = -1.;      // sentinel – will stay ≥ 0 once filled
    const bool   isUnc       = (std::strcmp(tag,"UNCORRECTED")==0);
    const bool   isCor       = (std::strcmp(tag,"CORRECTED"  )==0);

    /* keep the histogram self–contained, but postpone fixing the upper limit
       until we know the UNCORRECTED reference value                           */
    if (isUnc)
      uncorrZmax = res->GetMaximum();          // remember for later calls

    /* always use the UNCORRECTED range once it is known -------------------- */
    const double zMax =
      (uncorrZmax>0 ? uncorrZmax             // already defined → use it
                    : res->GetMaximum());    // first call & happens to be COR

    res->SetMinimum(0.0);
    res->SetMaximum(zMax);                    // <-- unified palette range

    /* ---------------------------------------------------------------------- */
    TCanvas c2(Form("c_res_%s",tag),"",1000,850);
    c2.SetLeftMargin(0.18); c2.SetBottomMargin(0.16); c2.SetRightMargin(0.12);

    gStyle->SetPalette(kBird);     gStyle->SetNumberContours(255);

    res->GetXaxis()->SetTitle("block #eta_{local, 2#times 2}");
    res->GetYaxis()->SetTitle("block #phi_{local, 2#times 2}");
    res->GetZaxis()->SetTitle("|#DeltaE| / #LT E #GT  [%]");
    res->Draw("COLZ");

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
    //_withVirgilesChange
  const char* inFile = "/Users/patsfan753/Desktop/PositionDependentCorrection/PositionDep_sim_ALL_test.root";

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
    
  OverlayDeltaEtaFiveWays(eEdges, binMode, isFirstPass,
          "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits/deltaEtaClusterizerVsFromScratch");
    
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


//  makeLegoGifHD(hUnc3D, "unc", "UNCORRECTED",
//                  hCor3D, "cor", "CORRECTED",
//                  out2DDir);                       // ← morph included



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
