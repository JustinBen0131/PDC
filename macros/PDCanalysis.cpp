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

/**
 * \brief Toggle: If isFirstPass=false, we overlay corrected histograms.
 *                If true, we do uncorrected only.
 */
static bool isFirstPass = false;

/**
 * \brief 8 energy bins: E_edges[] = {2,4,6,8,10,12,15,20,30}
 *        => N_E=8
 */
constexpr double E_edges[] = {2,4,6,8,10,12,15,20,30};
constexpr int    N_E = (sizeof(E_edges)/sizeof(E_edges[0])) - 1;

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
    const std::vector<double>& bScan_cm,    // in actual cm steps (like 0.00..1.60)
    double cellSize,                        // 5.55 cm
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir
)
{
  ScanResults results;
  results.tgVec.resize(N_E);

  std::cout << "\n[DEBUG] doAshScan => ENTER. suffix='"<<suffix<<"', #Ebins="<<N_E
            <<", #bScan="<<bScan_cm.size()<<", cellSize="<<cellSize<<"\n";

  for(int iE=0; iE<N_E; iE++)
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
      TString hName = Form("h_dx_ash_b%.4f_E%d", bVal, iE);
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
    int col = 1 + (iE % 10);
    g->SetMarkerStyle(20 + (iE % 10));
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
    const std::vector<double>& w0Scan,
    const std::function<double(TH1*,const TString&)>& sigmaFunc,
    const TString& suffix,
    const TString& histOutDir
)
{
  ScanResults results;
  results.tgVec.resize(N_E);
  results.minY=DBL_MAX;
  results.maxY=-DBL_MAX;

  std::cout<<"\n[DEBUG] doLogScan => ENTER. suffix='"<<suffix<<"', #Ebins="<<N_E
           <<" w0Scan.size="<<w0Scan.size()<<"\n";

  for(int iE=0; iE<N_E; iE++)
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
      TString hN  = Form("h_dx_log_w0%.2f_E%d", wVal, iE);

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
    int col = 2 + (iE % 10);
    g->SetMarkerStyle(21 + (iE % 10));
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


/**
 * \brief drawAshLogSideBySide
 *   Draws two sets of TGraphs:
 *    - left = asinh scan
 *    - right= log scan
 *   Each TGraph is for a different E-slice
 *   Also prints debug info about minY, maxY, best param
 */
void drawAshLogSideBySide(
    const ScanResults& ashRes,
    const ScanResults& logRes,
    const char* methodName,
    const double* E_edges,
    int N_E,
    const TString& baseDir
)
{
  std::cout << "\n[DEBUG] drawAshLogSideBySide => method='"<<methodName<<"', N_E="
            <<N_E<<", baseDir='"<<baseDir<<"'\n";

  // 1) Canvas with 2 pads
  TString cName = Form("cSideBySide_%s", methodName);
  TCanvas cSide(cName, cName, 1600,600);
  cSide.Divide(2,1);

  // 2) LEFT pad => Ash
  cSide.cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  double asMin = ashRes.minY;
  double asMax = ashRes.maxY;
  std::cout<<"[DEBUG] => Ash => global minY="<<asMin<<", maxY="<<asMax<<"\n";

  if(!ashRes.tgVec.empty() && ashRes.tgVec[0])
  {
    TGraph* firstGr = ashRes.tgVec[0];
    std::cout<<"   => first Ash TGraph has N="<<firstGr->GetN()<<" points.\n";

    firstGr->Draw("ALP");
    firstGr->GetXaxis()->SetTitle("b (cm)");
    firstGr->GetYaxis()->SetTitle("#sigma_{x} (cm)");

    double yLo = asMin - 0.1*fabs(asMin);
    double yHi = asMax + 0.1*fabs(asMax);
    if(yLo<0) yLo=0.;
    firstGr->GetYaxis()->SetRangeUser(yLo, yHi);

    // overlay the other E-bins
    for(size_t i=1; i<ashRes.tgVec.size(); i++)
    {
      TGraph* gA = ashRes.tgVec[i];
      if(gA && gA->GetN()>0) gA->Draw("LP SAME");
      else std::cout<<"   [WARN] => Ash i="<<i<<" => TGraph missing or empty!\n";
    }
  }
  else
  {
    std::cout<<"[WARN] => No Ash TGraph to draw => skipping.\n";
  }

  // build a legend for left pad
  TLegend legA(0.13,0.65,0.48,0.88);
  legA.SetBorderSize(0);
  legA.SetFillStyle(0);
  for(int iE=0; iE<N_E; iE++)
  {
    double bBest=(iE<(int)ashRes.bestParam.size()? ashRes.bestParam[iE]: 0.);
    double sBest=(iE<(int)ashRes.bestSigma.size()? ashRes.bestSigma[iE]: 0.);
    TGraph* gr= (iE<(int)ashRes.tgVec.size()? ashRes.tgVec[iE]: nullptr);

    if(!gr || gr->GetN()==0)
    {
      legA.AddEntry((TObject*)nullptr,
        Form("%.1f< E<%.1f => NO data", E_edges[iE], E_edges[iE+1]), "");
      continue;
    }
    legA.AddEntry(gr,
       Form("%.1f< E<%.1f (b=%.3g, #sigma=%.3f)", E_edges[iE], E_edges[iE+1],
            bBest, sBest),
       "lp");
  }
  legA.Draw();

  TLatex latA; latA.SetNDC(true);
  latA.SetTextSize(0.04);
  latA.DrawLatex(0.15,0.93, Form("Ash scan [%s]", methodName));

  // 3) RIGHT pad => Log
  cSide.cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();

  double lgMin = logRes.minY;
  double lgMax = logRes.maxY;
  std::cout<<"[DEBUG] => Log => global minY="<<lgMin<<", maxY="<<lgMax<<"\n";

  if(!logRes.tgVec.empty() && logRes.tgVec[0])
  {
    TGraph* firstGr= logRes.tgVec[0];
    std::cout<<"   => first Log TGraph => N="<<firstGr->GetN()<<" points.\n";

    firstGr->Draw("ALP");
    firstGr->GetXaxis()->SetTitle("w_{0}");
    firstGr->GetYaxis()->SetTitle("#sigma_{x} (cm)");

    double yLo = lgMin - 0.1*fabs(lgMin);
    double yHi = lgMax + 0.1*fabs(lgMax);
    if(yLo<0) yLo=0.;
    firstGr->GetYaxis()->SetRangeUser(yLo, yHi);

    // overlay for other E
    for(size_t i=1; i<logRes.tgVec.size(); i++)
    {
      TGraph* gL=logRes.tgVec[i];
      if(gL && gL->GetN()>0) gL->Draw("LP SAME");
      else std::cout<<"   [WARN] => Log i="<<i<<" => TGraph missing or empty!\n";
    }
  }
  else
  {
    std::cout<<"[WARN] => No Log TGraph to draw => skipping.\n";
  }

  // legend right
  TLegend legB(0.6,0.65,0.85,0.88);
  legB.SetBorderSize(0);
  legB.SetFillStyle(0);
  for(int iE=0; iE<N_E; iE++)
  {
    double wBest=(iE<(int)logRes.bestParam.size()? logRes.bestParam[iE]: 0.);
    double sBest=(iE<(int)logRes.bestSigma.size()? logRes.bestSigma[iE]: 0.);
    TGraph* gr= (iE<(int)logRes.tgVec.size()? logRes.tgVec[iE]: nullptr);

    if(!gr || gr->GetN()==0)
    {
      legB.AddEntry((TObject*)nullptr,
        Form("%.1f< E<%.1f => NO data", E_edges[iE], E_edges[iE+1]), "");
      continue;
    }
    legB.AddEntry(gr,
      Form("%.1f< E<%.1f (w0=%.2f, #sigma=%.3f)", E_edges[iE], E_edges[iE+1],
           wBest, sBest),
      "lp");
  }
  legB.Draw();

  TLatex latB; latB.SetNDC(true);
  latB.SetTextSize(0.04);
  latB.DrawLatex(0.15,0.93, Form("Log scan [%s]", methodName));

  // 4) Save
  TString outName = Form("%s/SideBySide_%s.png", baseDir.Data(), methodName);
  cSide.SaveAs(outName);
  std::cout<<"\n[INFO] => wrote side-by-side Ash/Log => "<<outName<<"\n";
}


/**
 * \brief Main function that orchestrates the scanning and plotting
 *        of Ash-b and Log-w0 histograms.
 *
 * Usage:   root -l plotAshLogRMS_sideBySide.cpp\(\"PositionDep_sim_ALL.root\"\)
 */
void plotAshLogRMS_sideBySide(const char* infile="PositionDep_sim_ALL.root")
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
    return;
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
  auto ashFit = doAshScan(fIn.get(),
                          N_E,
                          E_edges,
                          bScan_cm,
                          5.55, // cellSize
                          coreGaussianSigma, // function pointer
                          "fit",
                          histOutDir
                          );
  auto logFit = doLogScan(fIn.get(),
                          N_E,
                          E_edges,
                          w0Scan,
                          coreGaussianSigma,
                          "fit",
                          histOutDir
                          );
  drawAshLogSideBySide(ashFit, logFit, "FIT", E_edges, N_E, baseDir);

  // 5B) "RMS"
  std::cout<<"\n=== [STEP] Doing *RMS* approach with rawRMS ===\n";
  // Wrap rawRMS in a lambda if you like, or use directly
  auto ashRMS = doAshScan(fIn.get(),
                          N_E,
                          E_edges,
                          bScan_cm,
                          5.55,
                          rawRMS,
                          "rms",
                          histOutDir
                          );
  auto logRMS = doLogScan(fIn.get(),
                          N_E,
                          E_edges,
                          w0Scan,
                          rawRMS,
                          "rms",
                          histOutDir
                          );
  drawAshLogSideBySide(ashRMS, logRMS, "RMS", E_edges, N_E, baseDir);

  std::cout<<"\n[INFO] => plotAshLogRMS_sideBySide completed all tasks.\n\n";
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
        hUnc->GetXaxis()->SetTitle("block #phi_{CG}");
        hUnc->GetYaxis()->SetTitle("counts");
        hUnc->SetMarkerStyle(20);
        hUnc->Draw("E");
        vPhiUnc.emplace_back(std::move(hUnc));

        // ---------- 2)  asinh fit ------------------------------------------
        TF1 trial("asinhφ",asinhModel,-0.5,0.5,2);
        trial.SetParNames("Norm","bVal");
        trial.SetParLimits(0,1e-9,1e9);
        trial.SetParLimits(1,1e-5,2.0);

        double bestChi2 = 1e50, bestB=0.; int bestNdf=0;
        std::unique_ptr<TF1> bestFit;

        for(auto g : {std::pair{0.1,0.10}, {0.2,0.14},
                      std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(g.first,g.second);
            TFitResultPtr r = vPhiUnc.back()->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2() < bestChi2){
                bestChi2=r->Chi2(); bestB=trial.GetParameter(1); bestNdf=r->Ndf();
                bestFit = std::make_unique<TF1>(trial);
                bestFit->SetName(Form("bestF_phi_E%d",i));
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
                                Form("asinh fit (b=%.3g)",bestB),"l");
        lg.Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042);
        tl.DrawLatex(0.14,0.84,Form("b = %.3g",bestB));
        tl.DrawLatex(0.14,0.76,Form("#chi^{2}_{LL}/NDF = %.2f",
                        bestNdf>0? bestChi2/bestNdf : 0.));
        tl.SetTextAlign(32);
        tl.DrawLatex(0.88,0.85,Form("Uncorr: %.0f",nUnc));
        if(hCorRaw) tl.DrawLatex(0.88,0.82,Form("Corr: %.0f",nCor));

        if(bOut.is_open() && bestB>0)
            bOut << Form("PHI [%.1f,%.1f)  %.6f\n",eLo,eHi,bestB);
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
        hUnc->GetXaxis()->SetTitle("block #eta_{CG}");
        hUnc->GetYaxis()->SetTitle("counts");
        hUnc->SetMarkerStyle(20);
        hUnc->Draw("E");
        vEtaUnc.emplace_back(std::move(hUnc));

        TF1 trial("asinhη",asinhModel,-0.5,0.5,2);
        trial.SetParNames("Norm","bVal");
        trial.SetParLimits(0,1e-9,1e9);
        trial.SetParLimits(1,1e-5,2.0);

        double bestChi2=1e50,bestB=0.; int bestNdf=0;
        std::unique_ptr<TF1> bestFit;

        for(auto g : {std::pair{0.1,0.10}, {0.2,0.14},
                      std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(g.first,g.second);
            TFitResultPtr r = vEtaUnc.back()->Fit(&trial,"RQL0S","",-0.5,0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2()<bestChi2){
                bestChi2=r->Chi2(); bestB=trial.GetParameter(1); bestNdf=r->Ndf();
                bestFit=std::make_unique<TF1>(trial);
                bestFit->SetName(Form("bestF_eta_E%d",i));
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
                                Form("asinh fit (b=%.3g)",bestB),"l");
        lg.Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042);
        tl.DrawLatex(0.14,0.84,Form("b = %.3g",bestB));
        tl.DrawLatex(0.14,0.76,Form("#chi^{2}_{LL}/NDF = %.2f",
                        bestNdf>0? bestChi2/bestNdf : 0.));
        tl.SetTextAlign(32);
        tl.DrawLatex(0.88,0.85,Form("Uncorr: %.0f",nUnc));
        if(hCorRaw) tl.DrawLatex(0.88,0.82,Form("Corr: %.0f",nCor));

        if(bOut.is_open() && bestB>0)
            bOut << Form("ETA [%.1f,%.1f)  %.6f\n",eLo,eHi,bestB);
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
  lat.DrawLatex(0.16,0.93,"#bf{(a)}  Mean shift");

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
  leg.AddText("#color[kRed+1]{#square}  #Delta#LT#varphi#GT");
  leg.AddText("#color[kAzure+2]{#square}  RMS ratio");
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
          h2_unc->GetZaxis()->CenterTitle();          // optional
          h2_unc->GetZaxis()->SetTitleOffset(2.1);    // move label away
          h2_unc->Draw("COLZ");
          // do NOT delete here – wait until after SaveAs()
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
            h2_cor->GetXaxis()->SetTitle("blockEta_{CG}");
            h2_cor->GetYaxis()->SetTitle("blockPhi_{CG}");
            h2_cor->GetZaxis()->SetTitle("Energy [GeV]");
            h2_cor->GetZaxis()->SetTitleOffset(2.1);
            h2_cor->GetZaxis()->CenterTitle();
            h2_cor->Draw("COLZ");
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

        double bestChi2 = 1e50, bestB = 0.0;
        TF1*   bestFit  = nullptr;

        for(auto seed : { std::pair{0.1,0.10}, {0.2,0.14},
                          std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(seed.first, seed.second);
            TFitResultPtr r = h->Fit(&trial, "RQL0S", "", -0.5, 0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2() < bestChi2){
                bestChi2 = r->Chi2();
                bestB    = trial.GetParameter(1);
                if(bestFit) delete bestFit;
                bestFit   = new TF1(trial);
            }
        }
        if(bestFit){
            bestFit->SetLineColor(col);
            bestFit->SetLineWidth(2);
        }
        return {bestB, bestFit};
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
        //  draw – set a common y-range
        // -----------------------------------------------------------------
        const double ymax = std::max(hPhi->GetMaximum(), hEta->GetMaximum());
        hEta->GetYaxis()->SetRangeUser(0.0, 1.30 * std::max(1.0, ymax));

        hEta->Draw("E");            // blue first
        hPhi->Draw("E SAME");       // then red
        if(fEta) fEta->Draw("SAME");
        if(fPhi) fPhi->Draw("SAME");

        TPad* pad = (TPad*) gPad;                  // current pad
        TLegend* leg = new TLegend(0.15, 0.75,     // (x1,y1,x2,y2) in NDC
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

    for(auto obj : keepAlive) delete obj;
}

// ---------------------------------------------------------------------------
//  Scatter-plot of best-fit 𝗯 values vs. energy slice
//  • red squares = φ    (Y-projection)   b_{φ}
//  • blue circles = η   (X-projection)   b_{η}
//  • y-axis limits chosen automatically to include *all* points
// ---------------------------------------------------------------------------
void PlotBvaluesVsEnergy(TH3F*  hUnc3D,
                         const std::vector<std::pair<double,double>>& eEdges,
                         const char* outDir)
/* ------------------------------------------------------------------------ */
{
    if(!hUnc3D){
        std::cerr << "[bPlot] null hUnc3D – abort\n";
        return;
    }

    gStyle->SetOptStat(0);

    const int N_E = static_cast<int>(eEdges.size());

    std::vector<double> vX, vBphi, vBeta;    // store results
    vX.reserve(N_E);  vBphi.reserve(N_E);  vBeta.reserve(N_E);

    // ───────────────────────────────────────────────────────────────────
    // helper: identical to the fitter used elsewhere
    // ───────────────────────────────────────────────────────────────────
    auto fitAsinh = [&](TH1D* h)->double
    {
        TF1 trial("trial", asinhModel, -0.5, 0.5, 2);
        trial.SetParNames("Norm", "b");
        trial.SetParLimits(0, 1e-9, 1e9);
        trial.SetParLimits(1, 1e-5, 2.0);

        double bestChi2 = 1e50, bestB = 0.0;

        for(auto seed : { std::pair{0.1,0.10}, {0.2,0.14},
                          std::pair{0.3,0.20}, {0.5,0.08}})
        {
            trial.SetParameters(seed.first, seed.second);
            TFitResultPtr r = h->Fit(&trial, "RQL0S", "", -0.5, 0.5);
            if(!r.Get() || !r->IsValid()) continue;
            if(r->Chi2() < bestChi2){
                bestChi2 = r->Chi2();
                bestB    = trial.GetParameter(1);
            }
        }
        return bestB;
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

        const double bPhi = fitAsinh(hPhi);
        const double bEta = fitAsinh(hEta);

        vX   .push_back( 0.5*(eLo+eHi) );   // use bin centre
        vBphi.push_back( bPhi );
        vBeta.push_back( bEta );

        keepAlive.insert(keepAlive.end(), {hPhi, hEta});
    }

    // ───────────────────────────────────────────────────────────────────
    //  build graphs
    // ───────────────────────────────────────────────────────────────────
    TGraph* gPhi = new TGraph( (int)vX.size(), &vX[0], &vBphi[0] );
    TGraph* gEta = new TGraph( (int)vX.size(), &vX[0], &vBeta[0] );

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

    TLegend leg(0.15,0.78,0.40,0.90);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
    leg.AddEntry(gPhi,"b_{#varphi}","lp");
    leg.AddEntry(gEta,"b_{#eta}","lp");
    leg.Draw();

    keepAlive.insert(keepAlive.end(), {gPhi,gEta,hFrame});

    TString outName = Form("%s/bValues_vs_E.png", outDir);
    cB.SaveAs(outName);
    std::cout << "[bPlot] wrote " << outName << '\n';

    for(auto obj : keepAlive) delete obj;
}


/***************************************************************************
 * PlotChi2QA
 * ----------
 *  • inFile : full path to a ROOT file that contains
 *      h2_chi2_tot_etaPhi      (TH2F)   – all clusters, before χ² cut
 *      h2_chi2_rej_etaPhi      (TH2F)   – clusters rejected by χ² cut
 *      p_chi2_pass_etaPhi      (TProfile2D) – ⟨pass flag⟩ per tower
 *  • outDir : directory where the PNGs will be written (created if absent)
 *
 *  The function produces:
 *    (1) Chi2_Total.png        – total-occupancy heat-map
 *    (2) Chi2_Rejected.png     – rejected-occupancy heat-map
 *    (3) Chi2_PassFraction.png – ⟨pass⟩   (profile)   0 … 1
 *    (4) Chi2_RejectFraction.png – 1−⟨pass⟩  for convenience
 *
 *  Author: <you>, 2025-06-03
 ***************************************************************************/
void PlotChi2QA(const char* inFile,
                const char* outDir = "./Chi2_QA")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, /*recursive=*/true);

  /* ------------------------------------------------------------------ */
  /* 1) open file & fetch histograms                                    */
  /* ------------------------------------------------------------------ */
  std::unique_ptr<TFile> fin( TFile::Open(inFile,"READ") );
  if(!fin || fin->IsZombie()){
      std::cerr << "[Chi2QA]  ERROR – cannot open " << inFile << '\n';
      return;
  }
  auto* hTot = dynamic_cast<TH2F*>     ( fin->Get("h2_chi2_tot_etaPhi") );
  auto* hRej = dynamic_cast<TH2F*>     ( fin->Get("h2_chi2_rej_etaPhi") );
  auto* pPass= dynamic_cast<TProfile2D*>( fin->Get("p_chi2_pass_etaPhi") );

  if(!hTot || !hRej || !pPass){
      std::cerr << "[Chi2QA]  ERROR – one or more χ² QA histograms "
                   "missing in file:\n"
                << "          " << inFile << '\n';
      return;
  }

  /* detach from the file so we can close it right away */
  hTot ->SetDirectory(nullptr);
  hRej ->SetDirectory(nullptr);
  pPass->SetDirectory(nullptr);
  fin->Close();

  /* ------------------------------------------------------------------ */
  /* 2) build reject-fraction TH2F  =  (1 – pass)                       */
  /* ------------------------------------------------------------------ */
    // construct an empty TH2D that has exactly the same axes as pPass
  TH2D hRejFrac("h_chi2_rejFrac_etaPhi",
                  "Reject fraction after #chi^{2} cut;"
                  "Tower #eta index;Tower #varphi index;Reject fraction",
                  pPass->GetNbinsX(), pPass->GetXaxis()->GetXmin(), pPass->GetXaxis()->GetXmax(),
                  pPass->GetNbinsY(), pPass->GetYaxis()->GetXmin(), pPass->GetYaxis()->GetXmax());


  hRejFrac.SetName ("h_chi2_rejFrac_etaPhi");
  hRejFrac.SetTitle("Reject fraction after #chi^{2} cut;"
                    "Tower #eta index;Tower #varphi index;Reject fraction");

  const int nX = hRejFrac.GetNbinsX();
  const int nY = hRejFrac.GetNbinsY();
  for(int ix=1; ix<=nX; ++ix)
     for(int iy=1; iy<=nY; ++iy){
         double pass = pPass->GetBinContent(ix,iy);
         if(pass<0) pass = 0;          // profile may be empty
         hRejFrac.SetBinContent(ix,iy, 1.0 - pass);
     }

  /* colour palette nicer than default */
  gStyle->SetPalette(kViridis);

  /* helper lambda – draw & save one heat-map ------------------------- */
  auto drawSave = [&](TH2* h,
                      const char* cname,
                      const char* png)
  {
      TCanvas c(cname,cname,1000,760);
      c.SetRightMargin(0.14);
      h->Draw("COLZ");
      c.SaveAs(Form("%s/%s",outDir,png));
      std::cout << "[Chi2QA] wrote  " << outDir << '/' << png << '\n';
  };

  drawSave(hTot , "cChi2Tot" , "Chi2_Total.png");
  drawSave(hRej , "cChi2Rej" , "Chi2_Rejected.png");
  drawSave(pPass, "cChi2Pass", "Chi2_PassFraction.png");
  drawSave(&hRejFrac,"cChi2Frac","Chi2_RejectFraction.png");

  /* ------------------------------------------------------------------ */
  /* 3) global numbers                                                  */
  /* ------------------------------------------------------------------ */
  const double totEntries = hTot->GetEntries();
  const double rejEntries = hRej->GetEntries();
  std::cout << "\n[Chi2QA]  GLOBAL SUMMARY\n"
            << "          clusters before χ²  : " << totEntries << '\n'
            << "          clusters rejected   : " << rejEntries << '\n'
            << "          overall reject frac : "
            << (totEntries>0 ? rejEntries/totEntries : 0) << "\n\n";
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
 * OverlayDeltaPhiSlices  (v3 – legend bottom-left, compact stats)
 ******************************************************************************/
void OverlayDeltaPhiSlices(const std::vector<std::pair<double,double>>& eEdges,
                           bool   isFirstPass,
                           const char* outDir)
{
    gSystem->mkdir(outDir, /*recursive=*/true);
    gStyle->SetOptStat(0);

    const int N_E = static_cast<int>(eEdges.size());
    if(N_E==0){ std::cerr<<"[DeltaPhiOverlay] eEdges empty – abort\n"; return; }

    /* helper for Gaussian fit → N, μ, σ ---------------------------------- */
    struct Info { double N, mu, sigma; };
    auto getInfo = [](TH1* h)->Info
    {
        if(!h || h->Integral()==0) return {0,0,0};
        double m = h->GetMean(), r = h->GetRMS();
        TF1 g("g","gaus", m-2.5*r, m+2.5*r);
        g.SetParameters(h->GetMaximum(), m, r/2.);
        h->Fit(&g,"RQ0");
        return { h->GetEntries(), g.GetParameter(1), fabs(g.GetParameter(2)) };
    };

    /* summary canvas ------------------------------------------------------ */
    TCanvas c4x2("DeltaPhiOverlay_4by2",
                 "#Delta#phi raw vs corrected (normalised)",1600,1200);
    c4x2.Divide(4,2);

    std::vector<TObject*> keep;   // prevent premature deletion

    for(int i=0;i<N_E;++i)
    {
        const double eLo=eEdges[i].first, eHi=eEdges[i].second;
        TString rawName  = Form("h_phi_diff_raw_%.0f_%.0f", eLo, eHi);
        TString corrName = Form("h_phi_diff_corr_%.1f_%.1f", eLo, eHi);

        TH1F* hRawSrc  = dynamic_cast<TH1F*>( gROOT->FindObject(rawName)  );
        TH1F* hCorrSrc = (!isFirstPass) ? dynamic_cast<TH1F*>( gROOT->FindObject(corrName) )
                                        : nullptr;
        if(!hRawSrc){ std::cerr<<"[DeltaPhiOverlay] missing "<<rawName<<"\n"; continue; }

        /* clone (safe to modify) */
        TH1F* hRaw  = static_cast<TH1F*>(hRawSrc ->Clone(Form("hRaw_%d",i)));
        TH1F* hCorr = hCorrSrc ? static_cast<TH1F*>(hCorrSrc->Clone(Form("hCorr_%d",i))) : nullptr;
        hRaw ->SetDirectory(nullptr); if(hCorr) hCorr->SetDirectory(nullptr);
        keep.push_back(hRaw); if(hCorr) keep.push_back(hCorr);

        /* area-normalise --------------------------------------------------- */
        if(hRaw ->Integral()>0) hRaw ->Scale(1./hRaw ->Integral());
        if(hCorr && hCorr->Integral()>0) hCorr->Scale(1./hCorr->Integral());

        /* stats (after normalisation – entries unchanged) */
        Info R = getInfo(hRaw);
        Info C = hCorr ? getInfo(hCorr) : Info{0,0,0};

        /* basic style ------------------------------------------------------ */
        hRaw ->SetLineColor(kBlack); hRaw ->SetMarkerColor(kBlack);
        hRaw ->SetMarkerStyle(20);   hRaw ->SetMarkerSize(0.8);

        if(hCorr){
            hCorr->SetLineColor(kRed); hCorr->SetMarkerColor(kRed);
            hCorr->SetMarkerStyle(21); hCorr->SetMarkerSize(0.8);
        }

        /* ---- (A) summary pad -------------------------------------------- */
        c4x2.cd(i+1);
        TString ttl=Form("#Delta#phi  [%.1f, %.1f)  GeV",eLo,eHi);
        hRaw->SetTitle(ttl);
        hRaw->GetXaxis()->SetTitle("#Delta#phi (reco - truth)  [rad]");
        hRaw->GetYaxis()->SetTitle("Normalised counts");
        hRaw->GetYaxis()->SetTitleOffset(1.3);

        double ymax=hRaw->GetMaximum(); if(hCorr) ymax=std::max(ymax,hCorr->GetMaximum());
        hRaw->GetYaxis()->SetRangeUser(0,1.25*ymax);
        hRaw->Draw("E"); if(hCorr) hCorr->Draw("E SAME");

        /* legend – bottom-left */
        TLegend* legS = new TLegend(0.15,0.2,0.45,0.3);
        legS->SetBorderSize(0); legS->SetFillStyle(0); legS->SetTextSize(0.03);
        legS->AddEntry(hRaw ,"RAW"      ,"lp");
        if(hCorr) legS->AddEntry(hCorr,"CORRECTED","lp");
        legS->Draw(); keep.push_back(legS);

        /* text (compact) */
        TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.024);
        tl.SetTextAlign(13);
        tl.DrawLatex(0.15,0.85,Form("RAW : N=%d,  #mu=%.4f,  #sigma=%.4f",
                                    (int)R.N, R.mu, R.sigma));
        if(hCorr)
            tl.DrawLatex(0.15,0.81,Form("CORR: N=%d,  #mu=%.4f,  #sigma=%.4f",
                                        (int)C.N, C.mu, C.sigma));

        /* ---- (B) individual canvas -------------------------------------- */
        TString cNm = Form("cDeltaPhi_%d",i);
        TString cTi = Form("#Delta#phi overlay  %.1f–%.1f GeV",eLo,eHi);
        TCanvas cInd(cNm,cTi,800,600);
        cInd.SetLeftMargin(0.13); cInd.SetRightMargin(0.06); cInd.SetBottomMargin(0.12);

        hRaw->Draw("E"); hRaw->GetYaxis()->SetRangeUser(0,1.25*ymax);
        if(hCorr) hCorr->Draw("E SAME");

        TLegend legI(0.15,0.08,0.45,0.23);
        legI.SetBorderSize(0); legI.SetFillStyle(0); legI.SetTextSize(0.03);
        legI.AddEntry(hRaw ,"RAW"      ,"lp");
        if(hCorr) legI.AddEntry(hCorr,"CORRECTED","lp");
        legI.Draw();

        TLatex tlI; tlI.SetNDC(); tlI.SetTextFont(42); tlI.SetTextSize(0.030);
        tlI.SetTextAlign(13);
        tlI.DrawLatex(0.15,0.85,Form("RAW : N=%d,  #mu=%.4f,  #sigma=%.4f",
                                     (int)R.N, R.mu, R.sigma));
        if(hCorr)
            tlI.DrawLatex(0.15,0.81,Form("CORR: N=%d,  #mu=%.4f,  #sigma=%.4f",
                                         (int)C.N, C.mu, C.sigma));

        TString outPNG=TString(outDir)+Form("/DeltaPhiCompare_%.0fto%.0f.png",eLo,eHi);
        cInd.SaveAs(outPNG);
        std::cout<<"   ↳ wrote "<<outPNG<<'\n';
    }

    TString sumPNG=TString(outDir)+"/DeltaPhiCompare_AllOutput.png";
    c4x2.SaveAs(sumPNG);
    std::cout<<"[DeltaPhiOverlay] wrote summary → "<<sumPNG<<'\n';

    for(auto* obj:keep) delete obj;
}


/**
 * \brief LocalPhiOnly2by2
 * Reads two TH3Fs:
 *   "h2_cluster_block_cord_E"           (UNCORRECTED)
 *   "h2_cluster_block_cord_E_corrected" (CORRECTED, if isFirstPass==false)
 *
 * Each has axes:
 *   X => blockEta  [-0.5 .. +1.5], 14 bins
 *   Y => blockPhi  [-0.5 .. +1.5], 14 bins
 *   Z => cluster E, 8 bins with edges {2,4,6,8,10,12,15,20,30}
 *
 * Produces:
 * (A) 2D block‐eta vs block‐phi plots (4×2)
 * (B) 1D local‐φ distributions => asinh fits
 * (C) 1D local‐η distributions => asinh fits
 *
 * Then writes best-fit b-values to a text file (“bValues.txt”).
 *
 * Additional note:
 *   To properly overlay corrected local‐φ or local‐η, remember to
 *   set the Z axis to your E slice and reset the *unused* local axis
 *   to its full range.
 */
void PDCanalysis()
{
  // 1) Style / stat
  gStyle->SetOptStat(0);
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
    
  std::cout << "\n[INFO] Now saving *every* histogram as a PNG => " << outDirAll << std::endl;

  TIter nextAll(fIn->GetListOfKeys());
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
      TString outPNG = TString(outDirAll) + "/" + histName + ".png";

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

  // 3) Grab uncorrected TH3F
  TH3F* hUnc3D = dynamic_cast<TH3F*>( fIn->Get("h2_cluster_block_cord_E") );
  if(!hUnc3D)
  {
    std::cerr << "[ERROR] 'h2_cluster_block_cord_E' not found. Aborting.\n";
    fIn->Close();
    return;
  }
  hUnc3D->Sumw2();
  std::cout << "[DEBUG] Uncorrected TH3F => hUnc3D\n";

  // 4) If isFirstPass==false => also get corrected TH3F
  TH3F* hCor3D = nullptr;
  if(!isFirstPass)
  {
    hCor3D = dynamic_cast<TH3F*>( fIn->Get("h2_cluster_block_cord_E_corrected") );
    if(!hCor3D)
    {
      std::cerr << "[WARN] isFirstPass=false but no 'h2_cluster_block_cord_E_corrected' found.\n";
    }
    else
    {
      hCor3D->Sumw2();
      std::cout << "[DEBUG] Corrected TH3F => hCor3D\n";
    }
  }

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
    
  PlotBvaluesVsEnergy(hUnc3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
  PlotPhiShiftAndWidth(hUnc3D, hCor3D, eEdges, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    
  OverlayDeltaPhiSlices(eEdges, isFirstPass, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/1DplotsAndFits");
    

  {
      gStyle->SetNumberContours(50);      // smooth colours

      auto beautify = [](TH3F* h)
      {
        // X axis
        h->GetXaxis()->SetTitle("block #eta_{CG}");
        h->GetXaxis()->CenterTitle();
        h->GetXaxis()->SetTitleOffset(1.45);
        h->GetXaxis()->SetLabelOffset(0.007);

        // Y axis
        h->GetYaxis()->SetTitle("block #phi_{CG}");
        h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetTitleOffset(1.90);
        h->GetYaxis()->SetLabelOffset(0.007);

        // Z axis
        h->GetZaxis()->SetTitle("Cluster E  [GeV]");
        h->GetZaxis()->CenterTitle();
        h->GetZaxis()->SetTitleOffset(1.35);
        h->GetZaxis()->SetLabelOffset(0.007);
      };

        auto drawWithEntries = [&beautify]               // <-- add &
                (TH3F* h,
                 const char* canvName,
                 const char* pngPath)
        {
          if (!h) return;

          TCanvas c(canvName, canvName, 1200, 900);
          beautify(h);                                   // now allowed
          h->Draw("LEGO2Z");

          TLatex txt;
          txt.SetNDC(kTRUE);
          txt.SetTextFont(42);
          txt.SetTextSize(0.035);
          txt.SetTextAlign(13);
          txt.DrawLatex(0.12, 0.92,
                        Form("#bf{Entries: %.0f}", h->GetEntries()));

          c.SaveAs(pngPath);
          std::cout << "[INFO] wrote 3-D plot → " << pngPath << '\n';
        };


      // --------------------------- un-corrected ------------------------------
      drawWithEntries(hUnc3D,
                      "c3D_unc",
                      Form("%s/2DPlots/h2_cluster_block_cord_E.png", outDir));

      // ----------------------------- corrected -------------------------------
      drawWithEntries(hCor3D,
                      "c3D_cor",
                      Form("%s/2DPlots/h2_cluster_block_cord_E_corrected.png", outDir));
  } // end LEGO block

  // after opening the file …
  TH1F* h_truth_vz = static_cast<TH1F*>( fIn->Get("h_truth_vz") );

  if(!h_truth_vz){
      std::cerr << "[ERROR] vertex-Z histograms not found!\n";
  } else {
      PlotVertexZTruthOnly(h_truth_vz, "/Users/patsfan753/Desktop/PositionDependentCorrection/SimOutput/QA");
  }

    
  plotAshLogRMS_sideBySide(inFile);
    
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
