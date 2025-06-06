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
#include <TFitResult.h>
#include <TLatex.h>
#include <vector>
#include <utility>
#include <cmath>
#include <memory>
#include <cfloat>  // for DBL_MAX

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
//  1) Projects local-φ or local-η from h3 in slices of pT (here E-slices).
//  2) Overlays all E-bins on one canvas => LocalPhiFits.png or LocalEtaFits.png
//  3) If coord=phi, also does single‐bin overlays of raw vs corrected Δφ
//  4) Finally, we produce a canvas with each bin's local‐φ distribution & fit.
//
// **Changed**: now it's a 4×2 canvas for 8 bins (instead of 2×2).
////////////////////////////////////////////////////////////////////////////////
void doLocalPhiEtaFits(TH3F* h3,
                       const double pT_bins_low[/*N_E*/],
                       const double pT_bins_high[/*N_E*/],
                       std::ofstream& bResults,
                       const TString& outDirPhi)
{
  // We have N_E=8 bins
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

  // Style arrays for the bins
  int colors[N_E]   = { kBlack, kRed+1, kBlue+1, kGreen+2,
                        kMagenta+1, kOrange+2, kAzure+2, kSpring+9 };
  int markers[N_E]  = {20, 20, 20, 20, 21, 21, 21, 21};

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
                 << "# E_low  E_high  b_phi\n";
      }
      else
      {
        bResults << "\n# -- Now processing local eta fits --\n"
                 << "# E_low  E_high  b_eta\n";
      }
    }

    // 1) Create a canvas for the “overlay” of all N_E=8 bins
    TCanvas cCoord(Form("cCoord_%s", coordName),
                   Form("Local %s distributions", coordName),
                   900, 700);
    cCoord.SetLeftMargin(0.12);
    cCoord.SetRightMargin(0.05);
    cCoord.SetBottomMargin(0.12);
    cCoord.cd();

    TLegend leg(0.18, 0.68, 0.55, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.028);

    // We'll store the 1D projections and best-fit b-values
    std::vector<TH1D*> hLocalVec(N_E, nullptr);
    std::vector<double> bestBvals(N_E, 0.0);

    // 2) Loop over the 8 E-slices => single overlay
    for(int i=0; i<N_E; i++)
    {
      // Determine Z range for that E slice
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
                  << " for E bin i=" << i << std::endl;
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
                  << coordName << ", E bin i=" << i
                  << " => bVal=0.\n";
      }
      bestBvals[i] = bVal;

      // Legend entry
      TString legEntry = Form("E=[%.1f,%.1f] GeV : b=%.3g",
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

    // Save that “overlay”
    TString mainOut = outDirPhi + "/" + outPNGname;
    cCoord.SaveAs(mainOut);
    std::cout << "[INFO] Wrote overlay " << coordName
              << " => " << mainOut << std::endl;

    // -----------------------------------------------------------------------
    // If it's local φ, also produce:
    //   (a) single‐bin overlay of raw vs. corrected Δφ
    //   (b) a **4×2** layout for 8 bin local‐φ distributions + asinh fits
    // -----------------------------------------------------------------------
    if (isPhi)
    {
        //----------------------------------------------------------------------
        // (a) Single‐bin overlay: raw vs corrected, if those histograms exist
        //----------------------------------------------------------------------
        for(int i=0; i<N_E; i++)
        {
            double ptLo = pT_bins_low[i];
            double ptHi = pT_bins_high[i];
            
            // 1) Retrieve the RAW Δφ histogram
            TString rawName = Form("h_phi_diff_raw_%.0f_%.0f", ptLo, ptHi);
            TH1F* hRaw = dynamic_cast<TH1F*>(gROOT->FindObject(rawName));
            if(!hRaw)
            {
                std::cerr << "[WARNING] No raw Δphi hist '" << rawName
                          << "' => skip overlay for E bin i=" << i << std::endl;
                continue;
            }
            
            // 2) Retrieve the CORRECTED Δφ histogram
            TString corrName = Form("h_phi_diff_corr_%.1f_%.1f", ptLo, ptHi);
            TH1F* hCorr = dynamic_cast<TH1F*>(gROOT->FindObject(corrName));
            if(!hCorr)
            {
                std::cerr << "[WARNING] No corrected Δphi hist '" << corrName
                          << "' => skip overlay for E bin i=" << i << std::endl;
                continue;
            }
            
            // 3) Make a canvas & style each histogram
            TString cName  = Form("cDeltaPhiCompare_%.1fto%.1f", ptLo, ptHi);
            TString cTitle = Form("#Delta#phi overlay [%.1f < E < %.1f]",
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
                               Form("Raw #Delta#phi (%.1f< E<%.1f)", ptLo, ptHi),
                               "lp");
            legSingle.AddEntry(hCorr,
                               Form("Corrected #Delta#phi (%.1f< E<%.1f)", ptLo, ptHi),
                               "lp");
            legSingle.Draw();
            
            // 6) Save
            TString singleOut = Form("%s/DeltaPhiCompare_%.1fto%.1f.png",
                                     outDirPhi.Data(), ptLo, ptHi);
            cSingle.SaveAs(singleOut);
            
            std::cout << "[INFO] Single-bin Δphi overlay E=["
            << ptLo << "," << ptHi
            << "] => " << singleOut << std::endl;
        } // end loop over i
        
        //----------------------------------------------------------------------
        // (b) A **4×2** canvas showing each bin's local‐φ distribution & best-fit,
        //     possibly overlaying the corrected local‐φ projection
        //----------------------------------------------------------------------
        TCanvas c2by2("cLocalPhi4by2", "Local φ fits per E bin", 1600, 1000);
        c2by2.Divide(4, 2);  // 4 columns, 2 rows = 8 pads

        // Attempt to retrieve the corrected TH3
        TH3F* h3corr = dynamic_cast<TH3F*>(gROOT->FindObject("h2_cluster_block_cord_E_corrected"));
        if(!h3corr)
        {
            std::cout << "[WARNING] 'h2_cluster_block_cord_E_corrected' not found; no overlay.\n";
        }
        
        for (int i=0; i<N_E; i++)
        {
            c2by2.cd(i+1);
            // from the big overlay
            TH1D* histLocal = hLocalVec[i];
            if (!histLocal) continue;  // skip if empty
            
            // Put E range + b param in the title
            double bVal = 0.0;
            TF1* bf = dynamic_cast<TF1*>(
                       histLocal->GetListOfFunctions()->FindObject(
                         Form("bestFunc_phi_%d", i)
                       )
                     );
            if(bf) bVal = bf->GetParameter(1);

            TString binTitle = Form("Local #phi: E=[%.1f,%.1f], b=%.3g",
                                    pT_bins_low[i], pT_bins_high[i], bVal);
            histLocal->SetTitle(binTitle);
            histLocal->GetXaxis()->SetTitle("local #phi_{CG}");
            histLocal->GetYaxis()->SetTitle("counts (normalized)");
            histLocal->Draw("E");
            
            // Overlay the best function
            if(bf) bf->Draw("SAME");
            
            //-------------------------------------------------------------
            // Overlap the corrected local-φ from h3corr, if present
            //-------------------------------------------------------------
            if (h3corr)
            {
                // Define the pT bin range
                int zLoCorr = h3corr->GetZaxis()->FindFixBin(pT_bins_low[i] + 1e-9);
                int zHiCorr = h3corr->GetZaxis()->FindFixBin(pT_bins_high[i] - 1e-9);
                if(zLoCorr < 1) zLoCorr=1;
                if(zHiCorr> h3corr->GetNbinsZ()) zHiCorr= h3corr->GetNbinsZ();
                
                // Project in the φ dimension => "y"
                h3corr->GetXaxis()->SetRange(1, h3corr->GetNbinsX());
                h3corr->GetZaxis()->SetRange(zLoCorr, zHiCorr);
                
                TH1D* hCorrLocal = (TH1D*) h3corr->Project3D("y");
                if (hCorrLocal)
                {
                    double integC = hCorrLocal->Integral();
                    if (integC > 1e-9) hCorrLocal->Scale(1./integC);
                    
                    hCorrLocal->SetLineColor(kMagenta+1);
                    hCorrLocal->SetMarkerColor(kMagenta+1);
                    hCorrLocal->SetMarkerStyle(25);
                    hCorrLocal->SetMarkerSize(1.0);
                    hCorrLocal->SetTitle("");
                    
                    hCorrLocal->Draw("SAME E");
                    
                    // Possibly add a small local legend
                    TLegend legCorr(0.50, 0.65, 0.88, 0.82);
                    legCorr.SetBorderSize(0);
                    legCorr.SetFillStyle(0);
                    legCorr.SetTextSize(0.035);
                    legCorr.AddEntry(histLocal, "Uncorrected", "lp");
                    legCorr.AddEntry(hCorrLocal, "Corrected", "lp");
                    legCorr.Draw("SAME");
                }
                else
                {
                    std::cout << "[WARNING] Could not project corrected local-φ for bin " << i << std::endl;
                }
            } // end if(h3corr)
        }
        
        // Save the 4×2 layout
        TString out2by2 = outDirPhi + "/LocalPhiFits_4by2.png";
        c2by2.SaveAs(out2by2);
        std::cout << "[INFO] Wrote 4×2 local φ bin-by-bin fits => "
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
// PDCanalysis()
//
//  1) Lists all histograms (like the original code).
//  2) Makes a 2D block-eta vs block-phi plot from h2_cluster_block_cord_E
//  3) Makes local-phi distributions for E slices [2,4], [4,6], [6,8], [8,10],
//     [10,12], [12,15], [15,20], [20,30], fits them, and records b-values.
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
  // 2) Grab the 3D histogram for block coords vs E
  ///////////////////////////////////////////////////////////////////////
  TH3F* h3 = (TH3F*) f->Get("h2_cluster_block_cord_E");
    
  if(!h3)
  {
    std::cerr << "[ERROR] 'h2_cluster_block_cord_E' not found => can't do block-eta vs. block-phi.\n";
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
  // 3) Create 2D plot from h3 by integrating over entire E
  ///////////////////////////////////////////////////////////////////////
    {
      // First, restore the Z-axis range to include all E
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
      std::cout << "[INFO] Wrote 2D block-eta vs block-phi (all E) => "
                << (outDir2D + "/BlockEtaPhi_2D.png") << std::endl;
    }

    ///////////////////////////////////////////////////////////////////
    // 3b) Also produce 2D plots for each E slice
    ///////////////////////////////////////////////////////////////////
    {
      for(int i=0; i<N_E; i++)
      {
        // Determine Z-range in the 3D histogram (E axis)
        int zLo = h3->GetZaxis()->FindFixBin( E_edges[i] + 1e-9 );
        int zHi = h3->GetZaxis()->FindFixBin( E_edges[i+1] - 1e-9 );
        if(zLo < 1) zLo = 1;
        if(zHi > h3->GetNbinsZ()) zHi = h3->GetNbinsZ();

        // Apply that Z-range, then project onto X-Y => 2D
        h3->GetZaxis()->SetRange(zLo, zHi);
        TH2F* h2_slice = (TH2F*) h3->Project3D("xy");

        // Name & title
        h2_slice->SetName( Form("h2_blockCoord2D_E%.0fto%.0f", E_edges[i], E_edges[i+1]) );
        h2_slice->SetTitle( Form("Block #eta vs Block #phi (%.1f < E < %.1f);block #eta;block #phi",
                                   E_edges[i], E_edges[i+1]) );

        // Draw
        TCanvas c2dSlice( Form("c2d_E%.0fto%.0f", E_edges[i], E_edges[i+1]),
                          "Block Eta-Phi (E slice)", 900, 700);
        c2dSlice.SetLeftMargin(0.12);
        c2dSlice.SetRightMargin(0.15);
        c2dSlice.SetBottomMargin(0.12);

        h2_slice->SetStats(false);
        h2_slice->Draw("COLZ");

        // Save a PNG for each slice
        TString outName = Form("/BlockEtaPhi_2D_E%.0fto%.0f.png",
                               E_edges[i], E_edges[i+1]);
        c2dSlice.SaveAs(outDir2D + outName);

        std::cout << "[INFO] Wrote 2D block-eta vs block-phi => "
                  << (outDir2D + outName)
                  << " for E=[" << E_edges[i] << ", " << E_edges[i+1] << "]" << std::endl;
      }
    }

  ///////////////////////////////////////////////////////////////////////
  // 4) Local φ slices in E bins => do 1D projections, fits
  ///////////////////////////////////////////////////////////////////////

  // We'll record b-values to a text file
  std::ofstream bResults(Form("%s/bValues.txt", baseDir.Data()));
  if(!bResults.is_open()) {
    std::cerr << "[WARNING] Could not open bValues.txt for writing.\n";
  } else {
    bResults << "# E_low  E_high   bValue\n";
  }

  double pT_bins_low[N_E], pT_bins_high[N_E];
  for(int j=0; j<N_E; j++){
      pT_bins_low[j]  = E_edges[j];
      pT_bins_high[j] = E_edges[j+1];
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

  // Finally call the side‐by‐side RMS vs b/w0
  plotAshLogRMS_sideBySide(filename);

  ///////////////////////////////////////////////////////////////////////
  // 6) Done
  ///////////////////////////////////////////////////////////////////////
  f->Close();
  delete f;

  std::cout << "[INFO] PDCanalysis() completed successfully.\n";
}
