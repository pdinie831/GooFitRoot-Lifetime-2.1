//
//
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
//#include <TCint.h>
#include <TGenericClassInfo.h> 
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TLegend.h>
#include <TMinuit.h>

#include "TRandom.h" 
#include  <TStopwatch.h>
#include "TH1F.h"
#include "TH2F.h"			// unused?
#include "TStyle.h"
#include "TCanvas.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/PDFs/combine/MappedPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/mypdf/RGaussianPdf.h>
#include <goofit/PDFs/mypdf/ErfcMassPdf.h>
#include <goofit/PDFs/mypdf/SimpleCheby2Pdf.h> 
#include <goofit/PDFs/mypdf/ExpGausProdBPdf.h>
#include <goofit/PDFs/mypdf/ErfcPolyPdf.h>
#include <goofit/PDFs/mypdf/ErfEffiBpPdf.h>
//#include "ExpGausPEEPdf.h" 
//#include "ExpGausMPdf.h" 
//#include "ExpGausWithIntPdf.h"
//#include "ExpGausPEEfixSigmaPdf.h" 
//#include "ExpGausProdBPdf.h"
//#include "ExpGausProdEffiBPdf.h"
//#include "ExpGausPEESigmaBPdf.h" 
//#include "PolyEffiPdf.h" 
//#include "ErfcPolyPdf.h"
//#include "ErfcMassPdf.h"
//#include "SigmoidBpPdf.h"
//#include <goofit/PDFs/mypdf/SigmoidBpPdf.h>
//#include "ErfEffiBpPdf.h"
//#include "SigmoidGausPdf.h"
//#include "GooFit/BivarGaussianConstrPdf.h"
//#include "GooFit/TrivarGaussianConstrPdf.h"
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

using namespace std; 
using namespace GooFit;
using namespace ROOT;
void fitTauSBModel();
GooFit::Application *app_ptr;
bool minuit1;


int main (int argc, char** argv) {

  TApplication tapp("TApp",&argc, argv);
  GooFit::Application app("testGoofit3DBp-2016 fit example", argc, argv);
  app_ptr = &app;
//  app.require_subcommand();

//  app.add_flag("--minuit1", minuit1, "Use Minuit 1 instead of Minuit 2");
 
  TStopwatch TimeWatch;
  TimeWatch.Start();

  fitTauSBModel(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  GOOFIT_PARSE(app);
  return 0 ;
}


void fitTauSBModel(){


  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);

   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gROOT->ForceStyle();
   gStyle->SetOptStat(000000);
   gStyle->SetOptFit(000000);


   TCanvas* c1 = new TCanvas("c1","Mass PLOTS",200,10,900,780);
   TCanvas* c2 = new TCanvas("c2","cTau PLOTS",200,10,900,780);
   TCanvas* c3 = new TCanvas("c3","STau PLOTS",200,10,900,780);
//   TCanvas*c2 = new TCanvas("c2","PLOTS",200,10,900,780);
   //TPad* pad1 = (TPad*)c1->GetPad(0);
//   TPad* pad2 = (TPad*)c2->GetPad(0);
   //pad1->SetLeftMargin(0.15); 
//   pad2->SetLeftMargin(0.15); 




//  Char_t    InputFileName[300] = "rooFit/Ntuple2012BPJpsiK_RD-SaraCuts-2016-MC.root";
//    Char_t    InputFileName[300] = "rooFit/test-Run2016G-Charmonium-v1.root";
//Char_t    InputFileName[300] = "rooFit/test-Run2016H-Charmonium-v2.root";
//  Char_t    InputFileName[300] = "testproof2DCut-Bp-SaraCuts-NewCutPDGMass.root";
//  Char_t    InputFileName[300] = "testproof2DCut-Run2016-SaraCuts-NewCutPDGMass-TriggerMatchMuons.root";
//      Char_t    InputFileName[300] = "testproof2DCut-Run2016-SaraCuts-NewCutPDGMass-Els3-testMu4.root";
      Char_t    InputFileName[300] = "testproof2DCut-Run2016-SaraCuts-NewCutPDGMass-Els3-testMu4-TrigMatch.root";
//      Char_t    InputFileName[300] = "testproof2DCut-Run2016-SaraCuts-NewCutPDGMass.root";
//    Char_t    InputFileName[300] = "testproof2DCut-RunHCDF2016-SaraCuts-NewCutPDGMass.root
  Char_t    InputTauBpTreeName[10]   = "TauBpTree";
  TFile*InputFile = TFile::Open(InputFileName,"READ","ROOT file");
  
   Char_t    OutFileName[300] =  "testGoofit3DBp-2016.root";
   gSystem->Exec(Form("mv %s %s.tmp",OutFileName,OutFileName));
   TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
   
   int PlotLineWidth = 1;
//   float PlotLineWidth = 1.2;
   float MarkerSize    = 0.35;
 
  double xBpMass;
//  double xBpTau;
  double xBpcTau;
  double xSBpcTau;
//  double c_const       = 0.0299792458;
//     double XMinSign = 5.12;
//     double XMaxSign = 5.44;
  double XMinSign = 5.1;
  double XMaxSign = 5.45;
  double BpMass   = 5.2784;
//double BpSigma  = 0.020;
//  double BpSigma  = 0.022;
  double BpSigma  = 0.028;

//  double NSigmaSB = 3.6;
  double NSigmaSBL = -2.;
  double NSigmaSBR = 0;
  double BiasSB   = 6;
 
//      double XMinSign = 5.12;
//      double XMaxSign = 5.44;
//double XMinSign = 5.15;
//double XMaxSign = 5.40;
//   double XMinSBL = 4.879000;
//   double XMaxSBL = 5.159000;
//   double XMinSBR = 5.399000;
//   double XMaxSBR = 5.679000;
  double XMinSBL = XMinSign - (BiasSB+NSigmaSBL)*BpSigma;
  double XMaxSBL = BpMass   -  BiasSB           *BpSigma;
//  double XMaxSBL = BpMass - BiasSB          *BpSigma;
  double XMinSBR = BpMass   +  BiasSB           *BpSigma;
  double XMaxSBR = XMaxSign + (BiasSB+NSigmaSBR)*BpSigma;
//  double XMaxSBR = BpMass +(BiasSB+NSigmaSB)*BpSigma;
  double XMin = 0.004;
  double XMax = 0.30;
  double SXMin = 0.0003;
  double SXMax = 0.005;
//  double SXMax = 0.0079;
//  double SXMax = 0.0079;
  double XStepSign = 0.001;
  double XStepcTau = 0.001;
//  double XStepScTau = 0.00002;
  double XStepMinuit = 0.00001;
  double XHScale = 4;
  
  double c_const       = 0.0299792458;

  char strbuffer[1000];


  printf("(xBpMass>%8f && xBpMass<%8f || xBpMass>%8f && xBpMass<%8f)\n",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR);

//  GooFit::Observable* xMass  = new GooFit::Observable("xMass",XMinSign, XMaxSign); 
  GooFit::Observable xMass("xMass",XMinSign, XMaxSign); 
  xMass.setNumBins( (XMaxSign -XMinSign)/XStepSign );
  TH1F HxMass( "HxMass" , "B^{+} Mass"    ,	     xMass.getNumBins(), xMass.getLowerLimit(), xMass.getUpperLimit());
  TH1F pdfHist("pdfHist", "B^{+} Mass Fit",  XHScale*xMass.getNumBins(), xMass.getLowerLimit(), xMass.getUpperLimit());
  TH1F sigHist("sigHist", "B^{+} Mass Fit",  XHScale*xMass.getNumBins(), xMass.getLowerLimit(), xMass.getUpperLimit());
  TH1F bkgHist("bgkHist", "B^{+} Mass Fit",  XHScale*xMass.getNumBins(), xMass.getLowerLimit(), xMass.getUpperLimit());

//  GooFit::Observable* xcTau  = new GooFit::Observable("xcTau",XMin, XMax); 
  GooFit::Observable xcTau("xcTau",XMin, XMax); 
  int NBINS = (XMax -XMin)/XStepcTau;
  xcTau.setNumBins(NBINS);
  GooFit::Observable xScTau("xScTau",SXMin, SXMax); 
  xScTau.setNumBins(NBINS);
  std::cout<<"xMass.getNumBins() = "<<xMass.getNumBins()<<std::endl;
  std::cout<<"xcTau.getNumBins() = "<<xcTau.getNumBins()<<std::endl;
  std::cout<<"xScTau.getNumBins() = "<<xScTau.getNumBins()<<std::endl;
  std::cout<<"xScTau.getNumBins() = "<<xScTau.getNumBins()<<std::endl;
 
  TH1F HxcTau(    "HxcTau"   , "B^{+} cTau",          xcTau.getNumBins(), xcTau.getLowerLimit(),  xcTau.getUpperLimit());
  TH1F HxcTauSB(  "HxcTauSB" , "B^{+} cTau SB",       xcTau.getNumBins(), xcTau.getLowerLimit(),  xcTau.getUpperLimit());
  TH1F HxScTau(   "HxScTau"  , "B^{+} cTau Sigma",    xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
  TH1F HxScTauSB( "HxScTauSB", "B^{+} cTau Sigma SB", xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//   TH1F pdf_cTau_Hist( "pdf_cTau_Hist" , "B+ cTau model   pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());
//   TH1F sig_cTau_Hist( "sig_cTau_Hist" , "B+ cTau signal  pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());
//   TH1F bkg_cTau_Hist( "bkg_cTau_Hist" , "B+ cTau bckg    pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());
//   TH1F pdf_STau_Hist( "pdf_STau_Hist" , "B+ ScTau model  pdf",    xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//   TH1F sig_STau_Hist( "sig_STau_Hist" , "B+ ScTau signal pdf",    xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//   TH1F bkg_STau_Hist( "bkg_STau_Hist" , "B+ ScTau bckg   pdf",    xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
// 
//   TH1F pdf_cTau_Hist2D( "pdf_cTau_Hist2D" , "B+ cTau model   pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());
//   TH1F sig_cTau_Hist2D( "sig_cTau_Hist2D" , "B+ cTau signal  pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());
//   TH1F bkg_cTau_Hist2D( "bkg_cTau_Hist2D" , "B+ cTau bckg    pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit());              

//  TH2F pdf_cTauSTau_Hist2D( "pdf_cTauSTau_Hist2D" , "B+ cTau model   pdf",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//  TH2F sig_cTauSTau_Hist2D( "sig_cTauSTau_Hist2D" , "B+ cTau model   sig",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//  TH2F bkg_cTauSTau_Hist2D( "bkg_cTauSTau_Hist2D" , "B+ cTau model   bkg",    xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
  TH2F pdf_cTauSTau_Hist2D( "pdf_cTauSTau_Hist2D" , "Bc+ cTau model   pdf",  XHScale*xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), XHScale*xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
  TH2F sig_cTauSTau_Hist2D( "sig_cTauSTau_Hist2D" , "Bc+ cTau model   sig",  XHScale*xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), XHScale*xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
  TH2F bkg_cTauSTau_Hist2D( "bkg_cTauSTau_Hist2D" , "Bc+ cTau model   bkg",  XHScale*xcTau.getNumBins(), xcTau.getLowerLimit(), xcTau.getUpperLimit(), XHScale*xScTau.getNumBins(), xScTau.getLowerLimit(), xScTau.getUpperLimit());
//
// Mass Spectrum
//
  GooFit::Variable mean  ("mean"  ,5.2785,XStepMinuit, 5., 5.5);
  GooFit::Variable sigma1("sigma1",0.0139,XStepMinuit, 0., 1.);
  GooFit::Variable sigma2("sigma2",0.0228,XStepMinuit, 0., 1.);
  GooFit::Variable sigma3("sigma3",0.0601,XStepMinuit, 0., 1.);

//   GaussianPdf* gauss1 = new GaussianPdf("gauss1", xMass, mean, sigma1);
//   GaussianPdf* gauss2 = new GaussianPdf("gauss2", xMass, mean, sigma2);
//   GaussianPdf* gauss3 = new GaussianPdf("gauss3", xMass, mean, sigma3);
  RGaussianPdf* gauss1 = new RGaussianPdf("gauss1", xMass, mean, sigma1);
  RGaussianPdf* gauss2 = new RGaussianPdf("gauss2", xMass, mean, sigma2);
  RGaussianPdf* gauss3 = new RGaussianPdf("gauss3", xMass, mean, sigma3);

/*   GooFit::Variable* meanBckgBp   = new GooFit::Variable("meanBckgBp" ,5.360,0.00001, 5., 5.5);
  GooFit::Variable* sigmaBckgBp  = new GooFit::Variable("sigmaBckgBp",0.030,0.00001, 0. , 1. );
  GooFit::Variable* meanBckgB0   = new GooFit::Variable("meanBckgB0" ,5.090,0.00001, 5. , 5.2);
  GooFit::Variable* sigmaBckgB0  = new GooFit::Variable("sigmaBckgB0",0.025,0.00001, 0. , 1.);
 */  
//   GooFit::Variable* meanBckgBp   = new GooFit::Variable("meanBckgBp" ,5.37 ,0,  5.2, 5.7);
//   GooFit::Variable* sigmaBckgBp  = new GooFit::Variable("sigmaBckgBp",0.033,0,  0.013 , 0.05 );
//   GooFit::Variable* meanBckgB0   = new GooFit::Variable("meanBckgB0" ,5.090,0,  5.0 , 5.15);
//   GooFit::Variable* sigmaBckgB0  = new GooFit::Variable("sigmaBckgB0",0.025,0,  0.01 , 0.05);
//   GooFit::Variable meanBckgBp ("meanBckgBp" ,5.35129e+00,0., 5.25, 5.45);
//   GooFit::Variable sigmaBckgBp("sigmaBckgBp",2.10790e-02,0., 0.013 , 0.04 );
  GooFit::Variable meanBckgBp ("meanBckgBp" ,5.35129e+00);
  GooFit::Variable sigmaBckgBp("sigmaBckgBp",2.10790e-02);
//  GooFit::Variable* meanBckgB0   = new GooFit::Variable("meanBckgB0" ,5.010,0., 5.0 , 5.15);
//   GooFit::Variable meanBckgB0 ("meanBckgB0" ,5.10,0., 5.0. , 5.15);
//   GooFit::Variable sigmaBckgB0("sigmaBckgB0",0.029,0., 0.02 , 0.04);
  GooFit::Variable meanBckgB0 ("meanBckgB0" ,5.10 );
  GooFit::Variable sigmaBckgB0("sigmaBckgB0",0.029);
  GooFit::Variable	   wb1("wb1",1.75628e-01, 0., 1.);
  GooFit::Variable	   wb2("wb2",1.17246e-01, 0., 1.);
  GooFit::Variable	   wb3("wb3",0.1, 0., 1.);

  RGaussianPdf* gaussBckgBp = new RGaussianPdf("gaussBckgBp", xMass, meanBckgBp, sigmaBckgBp);
  RGaussianPdf* gaussBckgB0 = new RGaussianPdf("gaussBckgB0", xMass, meanBckgB0, sigmaBckgB0);
  

  GooFit::Variable* wg1 = new GooFit::Variable("wg1",0.44, 0., 1.);
  GooFit::Variable* wg2 = new GooFit::Variable("wg2",0.5 , 0., 1.);
  GooFit::Variable* signalYield = new GooFit::Variable("signalYield",8.88646e+05,      100000., 2000000.);
  GooFit::Variable* bckgYield   = new GooFit::Variable("bckgYield"  ,1.25481e+05 ,      10000., 400000.);

 
//   GooFit::Variable* constaCoef = new GooFit::Variable("constaCoef", 70., XStepMinuit, 20., 1000); 
//   GooFit::Variable* linearCoef = new GooFit::Variable("linearCoef", 0.1, XStepMinuit, -3.5, 10.); 
//   GooFit::Variable* secondCoef = new GooFit::Variable("secondCoef", 0.1, XStepMinuit, 0, 10);
//   GooFit::Variable* thirdCoef  = new GooFit::Variable("thirdCoef" , 0.1, XStepMinuit, 0, 10);
//  GooFit::Variable* constaCoef = new GooFit::Variable("constaCoef", 1. ,XStepMinuit,0.,1000. ); 
//  GooFit::Variable* linearCoef = new GooFit::Variable("linearCoef", 0.001,XStepMinuit,0,10 ); 

//  GooFit::Variable* p0   = new GooFit::Variable("p0", -2.31800e-01,-10.,10. ); 
//  GooFit::Variable  p1("p1",-2.78404e-01,-10.,10. ); 
  GooFit::Variable p0("p0",-3.08433e-01,-10.,10. ); 
  GooFit::Variable p1("p1"          ,0.,-10.,10. ); 
  GooFit::Variable VMinSign("VMinSign",XMinSign ); 
  GooFit::Variable VMaxSign("VMaxSign",XMaxSign ); 
  SimpleCheby2Pdf* SimpleCheby2  = new SimpleCheby2Pdf("SimpleCheby2", xMass, p0, p1,VMinSign,VMaxSign);

// double  fullRange = XMaxSign - XMinSign;
// double  minScaled = -1. + 2. * (XMinSign - xminfull) / fullRange;
// 
// double  maxScaled = +1. - 2. * (xmaxfull - XMaxSign)) / fullRange; 

//  GooFit::Variable* aslope     = new GooFit::Variable("slope", -1.);
  //GooFit::Variable* aslope     = new GooFit::Variable("slope", 0.39, -10, 10);
  //GooFit::Variable* apower     = new GooFit::Variable("apower", 6, 0, 10);
//  GooFit::Variable* apower     = new GooFit::Variable("apower", 1.18, XStepMinuit, 0.9, 15.);
//  GooFit::Variable* apower     = new GooFit::Variable("apower", 1.18, XStepMinuit, 0.9, 6.);
//  GooFit::Variable* apower     = new GooFit::Variable("apower", 1.18, 0.001, 0.9, 5.);
  //GooFit::Variable* treshold   = new GooFit::Variable("treshold" ,5.168,XStepMinuit, 5.02, 6.);
//  GooFit::Variable* treshold   = new GooFit::Variable("treshold" ,5.33,0, 5.04, 6.);

 
 
  std::vector<GooFit::Variable> weightsSignalMass;
  weightsSignalMass.push_back(*wg1);
//  weightsSignalMass.push_back(wg2); 

  std::vector<PdfBase*> compsSignalMass;
  compsSignalMass.push_back(gauss1);
  compsSignalMass.push_back(gauss2);
//  compsSignalMass.push_back(gauss3);
  string str = "signalMass";
  sprintf(strbuffer, "signalMass");
  GooFit::AddPdf signalMass(str, weightsSignalMass, compsSignalMass); 
//  GooFit::AddPdf * signalMass = new GooFit::AddPdf(str, weightsSignalMass, compsSignalMass); 
//  signalMass.addSpecialMask(PdfBase::ForceCommonNorm) ;

//  vector<GooFit::Variable*> weightsPoly;
//  weightsPoly.push_back(constaCoef);
//  weightsPoly.push_back(linearCoef);
//  weightsPoly.push_back(secondCoef);
//  weightsPoly.push_back(thirdCoef);

  
//  PolynomialPdf* polyTmp = new PolynomialPdf("polyTmp", xMass, weightsPoly); 
//   std::vector<PdfBase*> compsPoly2;
//   compsSignalMass.push_back(polyTmp);
//   compsSignalMass.push_back(polyTmp);
//   
//   ProdPdf* poly      = new ProdPdf("poly"  ,compsPoly2 );
 
//  PolynomialPdf* poly = new PolynomialPdf("poly", xMass, weightsPoly); 

//   vector<GooFit::Variable*> weightsErfcMassBckg;
//   weightsErfcMassBckg.push_back(constaCoef);
//   weightsErfcMassBckg.push_back(linearCoef);
//   weightsErfcMassBckg.push_back(secondCoef);
//   weightsErfcMassBckg.push_back(thirdCoef);
  GooFit::Variable ps0("ps0",4.44421e+01,  8.5 , 120.); 
  GooFit::Variable ps1("ps1",5.13630e+00,   5.1,5.2);
  GooFit::Variable ps2("ps2",1.); 
  GooFit::Variable ps3("ps3",0.);
  ErfcMassPdf* ErfcMassBckg = new ErfcMassPdf("ErfcMassBckg",xMass,ps0,ps1,ps2,ps3);
  
  std::vector<GooFit::Variable> weightsBckgMass;
  weightsBckgMass.push_back(wb1);
  weightsBckgMass.push_back(wb2);
//  weightsBckgMass.push_back(wb3);

//  ArgusPdf* argus = new  ArgusPdf("argus", xMass, treshold, aslope, true, apower);  

  std::vector<PdfBase*> compsBckgMass;
  compsBckgMass.push_back(gaussBckgBp);
//  compsBckgMass.push_back(gaussBckgB0);
//    compsBckgMass.push_back(argus);
//  compsBckgMass.push_back(poly);
  compsBckgMass.push_back(SimpleCheby2);
  compsBckgMass.push_back(ErfcMassBckg);
 
  GooFit::AddPdf bckgMass("bckgMass", weightsBckgMass, compsBckgMass);
//  bckgMass.addSpecialMask(PdfBase::ForceCommonNorm) ;

  
//==============================================================================
//==============================================================================
//==============================================================================
// Lifetime
//==============================================================================
//==============================================================================
//==============================================================================

     GooFit::Variable cTau ("cTau"  ,1./(1.638 *c_const), 10., 30.);
//     GooFit::Variable tauSB1("tauSB1",1./(1.620 *c_const), 0., 1000.);
//     GooFit::Variable tauSB2("tauSB2",1./(0.400 *c_const), 0., 1000.);

     GooFit::Variable tauSB1("tauSB1",1/6.31009e-03, 0., 2000.);
     GooFit::Variable tauSB2("tauSB2",1/4.72861e-02, 0., 1000.);
     GooFit::Variable tauSB3("tauSB3",1/4.71086e-02, 0., 1000.);
     

//GooFit::Variable cTau  ("cTau"  ,1./( 1.638 *c_const),0.00, 1000.);
//GooFit::Variable tauSB1("tauSB1",1./( 1.440 *c_const),0.00, 1000.);
//GooFit::Variable tauSB2("tauSB2",1./( 1.600 *c_const),0.00, 1000.);

  GooFit::Variable meanResSign ("meanResSign"  ,0.);

  GooFit::Variable meanResBckg1("meanResBckg1" ,0.);
  GooFit::Variable meanResBckg2("meanResBckg2" ,0.);
//  GooFit::Variable meanResBckg2("meanResBckg2" ,0.,0,0.01);
  GooFit::Variable meanResBckg3("meanResBckg3" ,0.);

  GooFit::Variable sigmaRes    ("sigmaRes"	,0.0003,0,001);

  GooFit::Variable sigmaResBckg("sigmaResBckg"  ,0.0003,0,001);



//  GooFit::Variable* meanLandauErrSign      = new GooFit::Variable( "meanLandauErrSign"        ,0.0015 ,SXMin, SXMax);
//  GooFit::Variable* sigmaLandauErrorSign   = new GooFit::Variable( "sigmaLandauErrorSign"     ,0.0002 ,SXMin, SXMax);

//  GooFit::Variable* meanLandauErrBckg      = new GooFit::Variable( "meanLandauErrBckg"        ,0.0015 ,SXMin, SXMax);
//  GooFit::Variable* sigmaLandauErrorBckg   = new GooFit::Variable( "sigmaLandauErrorBckg"     ,0.0002 ,SXMin, SXMax);
 
//  GooFit::Variable* meanGaussianErrSign    = new GooFit::Variable( "meanGaussianErrSign"      ,0.0013  ,SXMin, SXMax);
//  GooFit::Variable* sigmaGaussianErrorSign = new GooFit::Variable( "sigmaGaussianErrorSign"   ,0.0003  ,0.00001, SXMax);

//  GooFit::Variable* meanGaussianErrBckg1   = new GooFit::Variable( "meanGaussianErrBckg1"      ,0.0013  ,0.,SXMin, SXMax);
//  GooFit::Variable* sigmaGaussianErrorBckg1= new GooFit::Variable( "sigmaGaussianErrorBckg1"   ,0.0003  ,0.,0.00001, SXMax);

//  GooFit::Variable* meanBifurGErrSign      = new GooFit::Variable( "meanBifurGErrSign"        ,0.0015 ,0., SXMax);
//  GooFit::Variable* sigmaLBifurGErrSign    = new GooFit::Variable( "sigmaLBifurGErrSign"      ,0.0003 ,0.00001, SXMax);
//  GooFit::Variable* sigmaRBifurGErrSign    = new GooFit::Variable( "sigmaRBifurGErrSign"      ,0.0009 ,0.00001, SXMax);

//  GooFit::Variable* meanBifurGErrBckg      = new GooFit::Variable( "meanBifurGErrBckg"	      ,0.0015 ,0., SXMax);
//  GooFit::Variable* sigmaLBifurGErrBckg    = new GooFit::Variable( "sigmaLBifurGErrBckg"      ,0.0003 ,0.00001, SXMax);
//  GooFit::Variable* sigmaRBifurGErrBckg    = new GooFit::Variable( "sigmaRBifurGErrBckg"      ,0.0009 ,0.00001, SXMax);
  
//   GooFit::Variable* tauErrSign   = new GooFit::Variable("tauErrSign",2100,0.,0., 10000.);
//   GooFit::Variable* tauErrBckg1  = new GooFit::Variable("tauErrBckg1",2100,0.,0., 10000.);

//   GooFit::Variable* meanGaussianErrSign1    = new GooFit::Variable( "meanGaussianErrSign1"     ,1.20736e-03  ,0.0001, SXMax);
//   GooFit::Variable* meanGaussianErrSign2    = new GooFit::Variable( "meanGaussianErrSign2"     ,1.66720e-03  ,0.0001, SXMax);
//   GooFit::Variable* meanGaussianErrSign3    = new GooFit::Variable( "meanGaussianErrSign3"     ,1.72100e-03  ,0.0001, SXMax);
  GooFit::Variable meanGaussianErrSign1( "meanGaussianErrSign1"     ,1.19914e-03  ,0.0001, SXMax);
  GooFit::Variable meanGaussianErrSign2( "meanGaussianErrSign2"     ,1.66008e-03  ,0.0001, SXMax);
  GooFit::Variable meanGaussianErrSign3( "meanGaussianErrSign3"     ,1.72100e-03  ,0.0001, SXMax);

//   GooFit::Variable* sigmaGaussianErrorSign1 = new GooFit::Variable( "sigmaGaussianErrorSign1"  ,2.57463e-04  ,0.0001, 0.001);
//   GooFit::Variable* sigmaGaussianErrorSign2 = new GooFit::Variable( "sigmaGaussianErrorSign2"  ,2.82726e-04  ,0.0001, 0.001);
//   GooFit::Variable* sigmaGaussianErrorSign3 = new GooFit::Variable( "sigmaGaussianErrorSign3"  ,4.55213e-04  ,0.0001, 0.001);
  GooFit::Variable sigmaGaussianErrorSign1( "sigmaGaussianErrorSign1"  ,2.83235e-04  ,0.00001, 0.001);
  GooFit::Variable sigmaGaussianErrorSign2( "sigmaGaussianErrorSign2"  ,2.83263e-04  ,0.00001, 0.001);
  GooFit::Variable sigmaGaussianErrorSign3( "sigmaGaussianErrorSign3"  ,4.55213e-04  ,0.00001, 0.001);

  GooFit::Variable meanGaussianErrBckg1   ( "meanGaussianErrBckg1"	,1.51119e-03 ,0.00005, 0.005);
  GooFit::Variable meanGaussianErrBckg2   ( "meanGaussianErrBckg2"	,1.51338e-03 ,0.00005, 0.005);
  GooFit::Variable meanGaussianErrBckg3   ( "meanGaussianErrBckg3"	,1.86837e-03 ,0.00005, 0.005);

  GooFit::Variable sigmaGaussianErrorBckg1( "sigmaGaussianErrorBckg1"	,3.65925e-04,0.0001, 0.001);
  GooFit::Variable sigmaGaussianErrorBckg2( "sigmaGaussianErrorBckg2"	,3.67023e-04,0.0001, 0.001);
  GooFit::Variable sigmaGaussianErrorBckg3( "sigmaGaussianErrorBckg3"	,3.17983e-04,0.0001, 0.001);


//   GooFit::Variable* tauErrSign1  = new GooFit::Variable("tauErrSign1" ,1.49315e+03,1000., 20000.);
//   GooFit::Variable* tauErrSign2  = new GooFit::Variable("tauErrSign2" ,3.39595e+03,1000., 11000.);
  GooFit::Variable tauErrSign1("tauErrSign1" ,3.38691e+03,1000., 22000.);
  GooFit::Variable tauErrSign2("tauErrSign2" ,1.49315e+04,1000., 22000.);
  GooFit::Variable tauErrSign3("tauErrSign3" ,1.02715e+04,1000., 11000.);

  GooFit::Variable tauErrBckg1("tauErrBckg1",1/3.78093e-04, 2000., 18000.);
  GooFit::Variable tauErrBckg2("tauErrBckg2",1/3.71851e-04, 1600., 20000.);
  GooFit::Variable tauErrBckg3("tauErrBckg3",1/2.76171e-04, 2000., 10000);

  GooFit::Variable wt1("wt1",3.17124e-02, 0., 1.);
  GooFit::Variable wt2("wt2",0.05, 0., 1.);
//Break
//GenTrue break 5 sigma 0.01-0.35
//  double ef0Param = 8.86226e-02       ;
//  double ef1Param = 2.45898e+02       ;
//  double ef2Param =-5.56651e-02       ;
//  double ef3Param =-2.09280e-01       ;
//  double ef4Param = 3.23734e+00       ;
//  double ef5Param =-2.48160e+01       ;
//  double ef6Param = 7.36517e+01       ;
//  double ef7Param =-7.46837e+01       ;
//GenTrue break 5 sigma 0.006-0.35 TriggerMatchMuons
// 	 double ef0Param =   	1.66702e-01     ;
// 	 double ef1Param =   	3.02277e+02     ;
// 	 double ef2Param =     -1.34226e-01     ;
// 	 double ef3Param =     -1.89685e-01     ;
// 	 double ef4Param =   	2.93533e+00     ;
// 	 double ef5Param =     -2.27782e+01     ;
// 	 double ef6Param =   	6.74533e+01     ;
// 	 double ef7Param =     -6.77444e+01     ;
// 	     double ef0Param =  1.66612e-01  	;
// 	     double ef1Param =  3.02190e+02  	;
// 	     double ef2Param = -1.34140e-01  	;
// 	     double ef3Param = -1.89753e-01  	;
// 	     double ef4Param =  2.93743e+00  	;
// 	     double ef5Param = -2.27969e+01  	;
// 	     double ef6Param =  6.75209e+01  	;
// 	     double ef7Param = -6.78281e+01  	;
//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.015-0.35 Els>3
//
// 	   double ef0Param =   2.62140e-02     ;
// 	   double ef1Param =   8.93091e+01     ;
// 	   double ef2Param =   1.34002e-02     ;
// 	   double ef3Param =  -4.58456e-01     ;
// 	   double ef4Param =   6.60342e+00     ;
// 	   double ef5Param =  -4.54769e+01     ;
// 	   double ef6Param =   1.32014e+02     ;
// 	   double ef7Param =  -1.36165e+02     ;
// 	   double sef0Param =  2.18070e-03      ;
// 	   double sef1Param =  8.29254e+00      ;
//            double rhoParam  =  0.714	        ;
//  	   double sef2Param =  2.35149e-03      ;
//   GooFit::Variable* sef11 = new GooFit::Variable("sef11", 4.755e-06 ); 
//   GooFit::Variable* sef22 = new GooFit::Variable("sef22", 6.877e+01 ); 
//   GooFit::Variable* sef33 = new GooFit::Variable("sef33", 5.529e-06 ); 
//   GooFit::Variable* sef12 = new GooFit::Variable("sef12", 1.292e-02 ); 
//   GooFit::Variable* sef13 = new GooFit::Variable("sef13",-4.998e-06 ); 
//   GooFit::Variable* sef23 = new GooFit::Variable("sef23",-1.649e-02 ); 
//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35 Els>3 abs(dY)<1.2
// 
// 	   double ef0Param =   1.29837e-01  	;
// 	   double ef1Param =   3.11998e+02  	;
// 	   double ef2Param =  -1.05559e-01  	;
// 	   double ef3Param =  -1.28769e-01  	;
// 	   double ef4Param =   2.00568e+00  	;
// 	   double ef5Param =  -1.55316e+01  	;
// 	   double ef6Param =   4.48051e+01  	;
// 	   double ef7Param =  -4.31350e+01  	;
//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35 Els>3 abs(dY)<1.6

// 	   double ef0Param =   1.57339e-01  	;
// 	   double ef1Param =   3.06289e+02  	;
// 	   double ef2Param =  -1.27415e-01  	;
// 	   double ef3Param =  -1.77634e-01  	;
// 	   double ef4Param =   2.77372e+00  	;
// 	   double ef5Param =  -2.15556e+01  	;
// 	   double ef6Param =   6.40929e+01  	;
// 	   double ef7Param =  -6.47876e+01  	;
//    	   double sef0Param =	 3.30840e-03	;
//    	   double sef1Param =	 3.85188e+00	;
//      	   double rhoParam  =	 0.893  	;
// 	   double sef2Param =   3.35169e-03     ;
// 
//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35 Els>3
// 	   double ef0Param =  1.69773e-01	  ;
// 	   double ef1Param =  3.04705e+02	  ;
// 	   double ef2Param = -1.37350e-01	  ;
// 	   double ef3Param = -1.85806e-01	  ;
// 	   double ef4Param =  2.88660e+00	  ;
// 	   double ef5Param = -2.25131e+01	  ;
// 	   double ef6Param =  6.68356e+01	  ;
// 	   double ef7Param = -6.72583e+01	  ;
//..Erfc
	   double ef0Param =  1.31463e+02  ;
	   double ef1Param =  3.89558e-03  ;
	   double ef2Param =  3.15785e-02  ;
	   double ef3Param = -1.48208e-01  ;
	   double ef4Param =  2.29056e+00  ;
	   double ef5Param = -1.84320e+01  ;
	   double ef6Param =  5.44423e+01  ;
	   double ef7Param = -5.35583e+01  ;
   	   double sef0Param =  1.44466e+00   ;
   	   double sef1Param =  3.90037e-03   ;
	   double sef2Param =  3.15401e-02   ;
     	   double rhoParam  =  0.587	     ;
 	   GooFit::Variable sef11("sef11", 2.087e+00 );
 	   GooFit::Variable sef22("sef22", 3.261e-10 );
 	   GooFit::Variable sef33("sef33", 2.511e-08 );
 	   GooFit::Variable sef12("sef12", 1.531e-05 );
 	   GooFit::Variable sef13("sef13",-1.406e-04 );
 	   GooFit::Variable sef23("sef23",-3.606e-10 );
// 	   double ef0Param =  1.70441e-01	  ;
// 	   double ef1Param =  3.05321e+02	  ;
// 	   double ef2Param = -1.38062e-01	  ;
// 	   double ef3Param = -1.86543e-01	  ;
// 	   double ef4Param =  2.90080e+00	  ;
// 	   double ef5Param = -2.26120e+01	  ;
// 	   double ef6Param =  6.71608e+01	  ;
// 	   double ef7Param = -6.76570e+01	  ;
//    	   double sef0Param =	 3.44552e-03	;
//    	   double sef1Param =	 3.05321e+02	;
// 	   double sef2Param =    3.49278e-03    ;
//      	   double rhoParam  =	 0.895  	;
//   GooFit::Variable* sef11 = new GooFit::Variable("sef11", 1.187e-05 ); 
//   GooFit::Variable* sef22 = new GooFit::Variable("sef22", 1.379e+01 ); 
//   GooFit::Variable* sef33 = new GooFit::Variable("sef33", 1.220e-05 ); 
//   GooFit::Variable* sef12 = new GooFit::Variable("sef12", 1.145e-02 ); 
//   GooFit::Variable* sef13 = new GooFit::Variable("sef13",-1.202e-05 ); 
//   GooFit::Variable* sef23 = new GooFit::Variable("sef23",-1.181e-02 ); 
//	   
// 	   double ef0Param =  1.70820e-01	  ;
// 	   double ef1Param =  3.05391e+02	  ;
// 	   double ef2Param = -1.38397e-01	  ;
// 	   double ef3Param = -1.85847e-01	  ;
// 	   double ef4Param =  2.88704e+00	  ;
// 	   double ef5Param = -2.25124e+01	  ;
// 	   double ef6Param =  6.68205e+01	  ;
// 	   double ef7Param = -6.72272e+01	  ;
//    	   double sef0Param =	 3.47581e-03	;
//    	   double sef1Param =	 3.72946e+00	;
//      	   double rhoParam  =	 0.893  	;
// 	   double sef2Param =    3.52303e-03    ;
// 	   double sef0Param =	 4.39679e-05	;
// 	   double sef1Param =	 1.07946e+00	;
// 	   double rhoParam  =	 -0.433  	;
//   GooFit::Variable* sef11 = new GooFit::Variable("sef11", 1.208e-05 ); 
//   GooFit::Variable* sef22 = new GooFit::Variable("sef22", 1.391e+01 ); 
//   GooFit::Variable* sef33 = new GooFit::Variable("sef33", 1.241e-05 ); 
//   GooFit::Variable* sef12 = new GooFit::Variable("sef12", 1.162e-02 ); 
//   GooFit::Variable* sef13 = new GooFit::Variable("sef13",-1.223e-05 ); 
//   GooFit::Variable* sef23 = new GooFit::Variable("sef23",-1.198e-02 ); 

//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35 Els>3 ADAPTIVE
//
 
// 
// I Option i effi fit adaptive 
// 
 
//  
// 	   double ef0Param =    1.91776e-01 	;
// 	   double ef1Param =    3.23179e+02 	;
// 	   double ef2Param =   -1.59890e-01 	;
// 	   double ef3Param =   -1.56484e-01 	;
// 	   double ef4Param =    2.35446e+00 	;
// 	   double ef5Param =   -1.83996e+01 	;
// 	   double ef6Param =    5.29587e+01 	;
// 	   double ef7Param =   -5.05694e+01 	;
//  	   double sef0Param =	4.25355e-03	;
//   	   double sef1Param =	3.55237e+00	;
// 	   double rhoParam  =	0.845		;

// 	   double ef0Param =   2.07891e-02     ;
// 	   double ef1Param =   2.39721e+02     ;
// 	   double ef2Param =   1.08529e-02     ;
// 	   double ef3Param =  -1.55834e-01     ;
// 	   double ef4Param =   2.35394e+00     ;
// 	   double ef5Param =  -1.85136e+01     ;
// 	   double ef6Param =   5.23539e+01     ;
// 	   double ef7Param =  -4.69369e+01     ;
//  	   double sef0Param =  3.03492e-03     ;
//   	   double sef1Param =  2.98088e+00     ;
// 	   double rhoParam  =  0.883	       ;

//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35 Els>3 Moro Sigmoid
// 
// 
// 	   double ef0Param =     1.70820e-01  	  ;
// 	   double ef1Param =     3.05391e+02  	  ;
// 	   double ef2Param =    -8.10192e-01  	  ;
// 	   double ef3Param =    -1.08797e+00  	  ;
// 	   double ef4Param =     1.69011e+01  	  ;
// 	   double ef5Param =    -1.31790e+02  	  ;
// 	   double ef6Param =     3.91175e+02  	  ;
// 	   double ef7Param =    -3.93556e+02  	  ;
//   	   double sef0Param =	 4.93457e-03	;
//   	   double sef1Param =	 4.60564e+00	;
// 	   double rhoParam  =	 0.933  	;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//GenTrue break 5 sigma 0.006-0.35
//		       double ef0Param =     1.66983e-01    ;
//		       double ef1Param =     3.02531e+02    ;
//		       double ef2Param =    -1.34503e-01    ;
//		       double ef3Param =    -1.88564e-01    ;
//		       double ef4Param =     2.93013e+00    ;
//		       double ef5Param =    -2.28080e+01    ;
//		       double ef6Param =     6.77249e+01    ;
//		       double ef7Param =    -6.82386e+01    ;
////		       double sef0Param =    3.42514e-03    ;
////		       double sef1Param =    3.76065e+00    ;
////		       double rhoParam  =    0.896	    ;
//		       double sef0Param =    4.40192e-05    ;
//		       double sef1Param =    3.02532e+02    ;
// 		     double rhoParam  =    -0.434	  ;
//       	     double ef0Param = 1.66895e-01	 ;
//       	     double ef1Param = 3.02444e+02	 ;
//       	     double ef2Param =-1.34418e-01	 ;
//       	     double ef3Param =-1.88642e-01	 ;
//       	     double ef4Param = 2.93237e+00	 ;
//       	     double ef5Param =-2.28278e+01	 ;
//       	     double ef6Param = 6.77957e+01	 ;
//       	     double ef7Param =-6.83260e+01	 ;
// //////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//  GooFit::Variable* ef0 = new GooFit::Variable("ef0",  5.06489e-02); 
//  GooFit::Variable* uno = new GooFit::Variable("uno",  1); 
//   GooFit::Variable* ef0 = new GooFit::Variable("ef0",  ef0Param,ef0Param- 3.41501e-03,ef0Param+ 3.41501e-03); 
//   GooFit::Variable* ef1 = new GooFit::Variable("ef1",  ef1Param,ef1Param-3.67208e+00,ef1Param+3.67208e+00); 
 double NSigEffi0  =    0.5;
 double NSigEffi1  =    1.0;
 double NSigEffi2  =    0.5;
  GooFit::Variable efx("efx",  ef0Param,ef0Param - NSigEffi0*sef0Param,ef0Param + NSigEffi0*sef0Param); 
  GooFit::Variable efy("efy",  ef1Param,ef1Param - NSigEffi1*sef1Param,ef1Param + NSigEffi1*sef1Param); 
  GooFit::Variable efz("efz",  ef2Param,ef2Param - NSigEffi2*sef2Param,ef2Param + NSigEffi2*sef2Param); 
//  GooFit::Variable* efx = new GooFit::Variable("efx",  ef0Param); 
//  GooFit::Variable* efy = new GooFit::Variable("efy",  ef1Param); 
  GooFit::Variable ef0("ef0",  ef0Param); 
  GooFit::Variable ef1("ef1",  ef1Param); 
  GooFit::Variable ef2("ef2",  ef2Param); 
  GooFit::Variable ef3("ef3",  ef3Param); 
  GooFit::Variable ef4("ef4",  ef4Param);  
  GooFit::Variable ef5("ef5",  ef5Param); 
  GooFit::Variable ef6("ef6",  ef6Param); 
  GooFit::Variable ef7("ef7",  ef7Param); 

  GooFit::Variable sef0("sef0",  sef0Param); 
  GooFit::Variable sef1("sef1",  sef1Param); 
  GooFit::Variable rho ("rho" ,  rhoParam); 


//  coeffEffi.push_back(uno);
//   vector<GooFit::Variable*> coeffEffi;
//   coeffEffi.push_back(ef0);
//   coeffEffi.push_back(ef1);
//   coeffEffi.push_back(ef2);
//  coeffEffi.push_back(ef3);

 
//   PolyEffiPdf* Effi = new PolyEffiPdf("Effi", xcTau, coeffEffi);
// PolynomialPdf* Effi = new PolynomialPdf("Effi", xcTau, coeffEffi);

//  ErfcPolyPdf  *Effi = new ErfcPolyPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4);
//      SigmoidBpMoroPdf  *Effi = new SigmoidBpMoroPdf("Effi", xcTau, efx,efy,ef2,ef3,ef4,ef5,ef6,ef7);
//      SigmoidBpPdf  *Effi    = new SigmoidBpPdf("Effi", xcTau, efx,efy,efz,ef3,ef4,ef5,ef6,ef7);
//      SigmoidBpPdf  *Effi = new SigmoidBpPdf("Effi", xcTau, efx,efy,ef2,ef3,ef4,ef5,ef6,ef7);
//      SigmoidBpPdf  *Effi = new SigmoidBpPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7);
//      SigmoidGausPdf  *Effi = new SigmoidGausPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6);

//      SigmoidBpPdf  *EffiFIX = new SigmoidBpPdf("EffiFIX", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7);
      ErfEffiBpPdf  *EffiFIX = new ErfEffiBpPdf("EffiFIX", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7);


     GooFit::Variable XMinV ("XMinV"  ,XMin );
     GooFit::Variable XMaxV ("XMaxV"  ,XMax );
     GooFit::Variable SXMinV("SXMinV" ,SXMin);
     GooFit::Variable SXMaxV("SXMaxV" ,SXMax);
//  Effi->setParameterConstantness(true); 
//  pdfFitBckg1 ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//  Effi        ->addSpecialMask(PdfBase::ForceSeparateNorm); 
   
//    ExpGausProdBPdf* DecayBp	= new ExpGausProdBPdf("DecayBp"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign1,
//    XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp1	= new ExpGausProdBPdf("DecayBp1"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign1, meanGaussianErrSign1,tauErrSign1,
      XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* DecayBp1	= new ExpGausProdEffiBPdf("DecayBp1"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign1,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//     ExpGausProdBPdf* DecayBp2	= new ExpGausProdBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign2,
//     XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp2	= new ExpGausProdBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign1,
      XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp3 = new ExpGausProdBPdf("DecayBp3"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign3, meanGaussianErrSign3,tauErrSign1,
      XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* DecayBp2	= new ExpGausProdEffiBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign2,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* pdfFitBckg1 = new ExpGausProdEffiBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//        ExpGausProdBPdf* pdfFitBckg1 = new ExpGausProdBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
//        XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* pdfFitBckg1 = new ExpGausProdBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg2, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
      XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg2, tauSB2, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
      XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* pdfFitBckg3 = new ExpGausProdBPdf("pdfFitBckg3", xcTau, xScTau, meanResBckg2, tauSB3, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
      XMinV,XMaxV,SXMinV,SXMaxV);
//        ExpGausProdBPdf* pdfFitBckg3 = new ExpGausProdBPdf("pdfFitBckg3", xcTau, xScTau, meanResBckg3, tauSB3, sigmaGaussianErrorBckg3, meanGaussianErrBckg3,tauErrBckg3,
//        XMinV,XMaxV,SXMinV,SXMaxV);
 // ExpGausProdEffiBPdf* pdfFitBckg2 = new ExpGausProdEffiBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
 // XMinV,XMaxV,SXMinV,SXMaxV);
//       ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
//       XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg2,
//  XMinV,XMaxV,SXMinV,SXMaxV);

// ExpPdf* DecayBp     = new ExpPdf("DecayBp"    ,  xcTau, cTau  );
// ExpPdf* pdfFitBckg1 = new ExpPdf("pdfFitBckg1", xcTau, tauSB1);
// ExpPdf* pdfFitBckg2 = new ExpPdf("pdfFitBckg2", xcTau, tauSB2);


//  ExpGausPEEPdf* DecayBp     = new ExpGausPEEPdf("DecayBp"    , xcTau, xScTau, meanResSign   , cTau  );
//  ExpGausPEEPdf* pdfFitBckg1 = new ExpGausPEEPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1);
//  ExpGausPEEPdf* pdfFitBckg2 = new ExpGausPEEPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2);

  
//  ExpGausPdf* DecayBp	    = new ExpGausPdf("DecayBp"	, xcTau, meanResSign   , sigmaRes, cTau  );
//  ExpGausPdf* pdfFitBckg1 = new ExpGausPdf("pdfFitBckg1", xcTau, meanResBckg, sigmaResBckg, tauSB1);
//  ExpGausPdf* pdfFitBckg2 = new ExpGausPdf("pdfFitBckg2", xcTau, meanResBckg, sigmaResBckg, tauSB2);

//    std::vector<PdfBase*> compspdfFitBp1;
//    compspdfFitBp1.push_back(DecayBp1);
//    compspdfFitBp1.push_back(Effi);
//    ProdPdf* pdfFitBp1	 = new ProdPdf("pdfFitBp1"  , compspdfFitBp1);
// 
//    std::vector<PdfBase*> compspdfFitBp2;
//    compspdfFitBp2.push_back(DecayBp2);
//    compspdfFitBp2.push_back(Effi);
//    ProdPdf* pdfFitBp2	 = new ProdPdf("pdfFitBp2"  , compspdfFitBp2);
// 
//   std::vector<PdfBase*> compspdfDecayBpAdd;
//   compspdfDecayBpAdd.push_back(pdfFitBp1);
//   compspdfDecayBpAdd.push_back(pdfFitBp2);
//   GooFit::AddPdf *pdfFitBp = new GooFit::AddPdf("pdfFitBp", weightsSignalMass, compspdfDecayBpAdd);
// 
   std::vector<GooFit::Variable> weightsSignalTau;
   weightsSignalTau.push_back(wt1);
//   weightsSignalTau.push_back(wt2); 
   std::vector<PdfBase*> compspdfDecayBpAdd;
   compspdfDecayBpAdd.push_back(DecayBp1);
   compspdfDecayBpAdd.push_back(DecayBp2);
//   compspdfDecayBpAdd.push_back(DecayBp3);
 GooFit::AddPdf DecayBp("DecayBp", weightsSignalMass, compspdfDecayBpAdd);
//     GooFit::AddPdf DecayBp("DecayBp", weightsSignalTau, compspdfDecayBpAdd);
//  GooFit::AddPdf* pdfFitBp =  new GooFit::AddPdf("pdfFitBp", weightsSignalMass, compspdfDecayBpAdd);
  
 
//  DecayBp.addSpecialMask(PdfBase::ForceSeparateNorm); 
 
 
 
    std::vector<PdfBase*> compspdfFitBp;
    compspdfFitBp.push_back(&DecayBp);
//    compspdfFitBp.push_back(DecayBp1);
    compspdfFitBp.push_back(EffiFIX);
    ProdPdf* pdfFitBp	 = new ProdPdf("pdfFitBp"  , compspdfFitBp);

  GooFit::Variable b1("b1",1.21140e-01,0., 1.);
  GooFit::Variable b2("b2",5.38694e-02,0., 1.);
  std::vector<GooFit::Variable> weightspdfFitBckg;
  weightspdfFitBckg.push_back(b1);
//  weightspdfFitBckg.push_back(b2);

  
  std::vector<PdfBase*> compspdfFitBckgAdd;
  compspdfFitBckgAdd.push_back(pdfFitBckg1);
  compspdfFitBckgAdd.push_back(pdfFitBckg2);
//  compspdfFitBckgAdd.push_back(pdfFitBckg3);
//  GooFit::AddPdf pdfFitBckgAdd("pdfFitBckgAdd", weightspdfFitBckg, compspdfFitBckgAdd);
//  GooFit::AddPdf pdfFitBckgTmp("pdfFitBckgTmp",weightsBckgMass , compspdfFitBckgAdd);
  GooFit::AddPdf pdfFitBckgTmp("pdfFitBckgTmp", weightspdfFitBckg, compspdfFitBckgAdd);
//      GooFit::AddPdf pdfFitBckg("pdfFitBckg", weightspdfFitBckg, compspdfFitBckgAdd);
  
   std::vector<PdfBase*> compspdfFitBpBckg;
   compspdfFitBpBckg.push_back(&pdfFitBckgTmp);
//   compspdfFitBpBckg.push_back(pdfFitBckg1);
   compspdfFitBpBckg.push_back(EffiFIX);
   ProdPdf pdfFitBckg("pdfFitBckg"  , compspdfFitBpBckg);
   
////////////////////////////////////////////////////   
////////////////////////////////////////////////////   
// bivariate gaussian constraint
////////////////////////////////////////////////////   
////////////////////////////////////////////////////   

// BivarGaussianConstrPdf* EffiConstr = new BivarGaussianConstrPdf("EffiConstr",xcTau, efx,efy, ef0,sef0, ef1,sef1, rho);   
// TrivarGaussianConstrPdf* EffiConstr = new TrivarGaussianConstrPdf("EffiConstr",xcTau, efx,efy,efz,ef0,ef1,ef2,sef11,sef22,sef33,sef12,sef13,sef23);   

////////////////////////////////////////////////////   
////////////////////////////////////////////////////   
// bivariate gaussian constraint
////////////////////////////////////////////////////   
////////////////////////////////////////////////////   
   
   
 
//  std::vector<PdfBase*> compspdfFitBckg;
//  compspdfFitBckg.push_back(pdfFitBckg1);
//  compspdfFitBckg.push_back(&pdfFitBckgAdd);
//  compspdfFitBckg.push_back(Effi);
 
//  ProdPdf* pdfFitBckg   = new ProdPdf("pdfFitBckg", compspdfFitBckg);

// Res Models   
//  ExpGausPEESigmaBPdf* ExpGauSign = new ExpGausPEESigmaBPdf("ExpGauSig" , xScTau,  sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign1,
//    SXMinV,SXMaxV);
//  ExpGausPEESigmaBPdf* ExpGauBckg = new ExpGausPEESigmaBPdf("ExpGauBckg", xScTau,  sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
//    SXMinV,SXMaxV);
//  ExpGausPEESigmaBPdf* ExpGauBckg2 = new ExpGausPEESigmaBPdf("ExpGauBckg2", xScTau,  sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
//    SXMinV,SXMaxV);


//  GooPdf* LandauErrorSign = new LandauPdf("LandauErrorSign", xScTau, meanLandauErrSign, sigmaLandauErrorSign);
//  GooPdf* LandauErrorBckg = new LandauPdf("LandauErrorBckg", xScTau, meanLandauErrBckg, sigmaLandauErrorBckg);

//  GooPdf* GaussianErrorSign = new GaussianPdf("GaussianErrorSign", xScTau, meanGaussianErrSign, sigmaGaussianErrorSign);
//  GooPdf* GaussianErrorBckg = new GaussianPdf("GaussianErrorBckg", xScTau, meanGaussianErrBckg1, sigmaGaussianErrorBckg);

//  GooPdf* BifurGErrorSign = new BifurGaussPdf("BifurGErrorSign", xScTau, meanBifurGErrSign, sigmaLBifurGErrSign,sigmaRBifurGErrSign);
//  GooPdf* BifurGErrorBckg = new BifurGaussPdf("BifurGErrorBckg", xScTau, meanBifurGErrBckg, sigmaLBifurGErrBckg,sigmaRBifurGErrBckg);

//  ExpGausPdf* ExpGauSign = new ExpGausPdf("ExpGauSig" , xScTau, meanGaussianErrSign, sigmaGaussianErrorSign, tauErrSign);
//  ExpGausPdf* ExpGauBckg = new ExpGausPdf("ExpGauBckg", xScTau, meanGaussianErrBckg1, sigmaGaussianErrorBckg1, tauErrBckg);

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//
// 2DFit
//  
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
  

  std::vector<PdfBase*> compsSignalLife;
  compsSignalLife.push_back(&signalMass);
  compsSignalLife.push_back(pdfFitBp);
//  compsSignalLife.push_back( GaussianErrorSign);
  //compsSignalLife.push_back(ExpGauSign);
//  compsSignalLife.push_back( LandauErrorSign);
  ///compsSignalLife.push_back(BifurGErrorSign);
//  compsSignalLife.push_back(DecayBp);

  std::vector<PdfBase*> compsBckgLife;
  compsBckgLife.push_back(&bckgMass);
  compsBckgLife.push_back(&pdfFitBckg);
//  compsBckgLife.push_back(pdfFitBckg1);
//  compsSignalLife.push_back( GaussianErrorBckg);
  //compsSignalLife.push_back(ExpGauBckg);
//  compsSignalLife.push_back( LandauErrorBckg);
  //compsSignalLife.push_back(BifurGErrorBckg);
//compsBckgLife.push_back(&pdfFitBckgAdd);

  ProdPdf* signalLife = new ProdPdf("signalLife", compsSignalLife);
  ProdPdf* bckgLife   = new ProdPdf("bckgLife  ", compsBckgLife);

  std::vector<GooFit::Variable> weightsYield;
  weightsYield.push_back(*signalYield);
  weightsYield.push_back(*bckgYield);
  

  std::vector<PdfBase*> compsModel;
  compsModel.push_back(signalLife);
//  compsModel.push_back(gaussBckgBp);
  compsModel.push_back(bckgLife);
  
  
//  compsModel.push_back(poly);
//  compsModel.push_back(argus);
  GooFit::AddPdf model("model", weightsYield, compsModel); 
//  GooFit::AddPdf model1("model1", weightsYield, compsModel); 

//  GooFit::Variable* uno = new GooFit::Variable("uno",0.5);
//  GooFit::Variable* L   = new GooFit::Variable("L",-100.,-100.,-100.);
//  std::vector<PdfBase*> compsModelConstr;
//    compsModelConstr.push_back(&model1);
//    compsModelConstr.push_back(EffiConstr);
//    compsModelConstr.push_back(signalLife);
//    compsModelConstr.push_back(bckgLife);
//  compsModelConstr.push_back(EffiConstr);

//  std::vector<GooFit::Variable*> weightsC;
//   weightsC.push_back(signalYield);
//   weightsC.push_back(bckgYield);
//  weightsC.push_back(uno);
//  weightsC.push_back(uno);

//  GooFit::AddPdf model("model", weightsC , compsModelConstr); 
//  ProdPdf model("model",  compsModelConstr); 


//  model.addSpecialMask(PdfBase::ForceCommonNorm) ;

//
// These are used for Plots....
//
  std::vector<PdfBase*> compsMass;
  compsMass.push_back(&signalMass);
  compsMass.push_back(&bckgMass);
  GooFit::AddPdf modelMass("modelMass", weightsYield, compsMass); 
  
  std::vector<PdfBase*> compscTau;
  compscTau.push_back(pdfFitBp);
  compscTau.push_back(&pdfFitBckg);
  GooFit::AddPdf model_cTau("model_cTau", weightsYield, compscTau); 
  
//  std::vector<PdfBase*> compsSTau;
//  compsSTau.push_back(LandauErrorSign);
//  compsSTau.push_back(LandauErrorBckg);
//  compsSTau.push_back(BifurGErrorSign);
//  compsSTau.push_back(BifurGErrorBckg);
//  compsSTau.push_back(GaussianErrorSign);
//  compsSTau.push_back(GaussianErrorBckg);
//  compsSTau.push_back(ExpGauSign);
//  compsSTau.push_back(ExpGauBckg);
//  GooFit::AddPdf model_STau("model_cTau", weightsYield, compsSTau); 
  
  

//
// Data
//
  std::vector<GooFit::Observable> dataVec;
  
  dataVec.push_back(xMass);
  dataVec.push_back(xcTau);
  dataVec.push_back(xScTau);
  UnbinnedDataSet* dataLife = new GooFit::UnbinnedDataSet(dataVec);
//
  if (!InputFile)
   {
     cout<<"File:"<<InputFileName<<" not found!!!"<<endl;
    exit(1);
   }
   InputFile->ls();
   
   TTree *TauBpTree    = (TTree*)InputFile->Get(InputTauBpTreeName);
   if(!TauBpTree ){
     cout<<"TTree cTau Data: "<< InputTauBpTreeName <<" not found!!!"<<endl;
     exit(1);
   }else{
     cout<<"TTree cTau Data: "<< InputTauBpTreeName <<" OK FOUND!!!"<<endl;
   }  
    
   TauBpTree->SetBranchAddress("xBpMass",&xBpMass);
//   TauBpTree->SetBranchAddress("xBpTau" ,&xBpTau);
   TauBpTree->SetBranchAddress("xBpcTau",&xBpcTau);
   TauBpTree->SetBranchAddress("xSBpcTau",&xSBpcTau);
   int nentries = (int)TauBpTree->GetEntries();
   
   for (Int_t i=0;i<nentries;i++) { 
    TauBpTree->GetEntry(i);
    if(xBpcTau>=XMin&&xBpcTau<=XMax&&xSBpcTau>=SXMin&&xSBpcTau<=SXMax){
     if(xBpMass>=XMinSign&&xBpMass<=XMaxSign){
      xMass.setValue(xBpMass)  ;
      xcTau.setValue(xBpcTau)   ;
      xScTau.setValue(xSBpcTau); 
      dataLife->addEvent();
      HxMass.Fill(xBpMass);
      HxcTau.Fill(xBpcTau);
      HxScTau.Fill(xSBpcTau);
     } 
     if(xBpMass>XMinSBL&&xBpMass<XMaxSBL){
        HxcTauSB.Fill(xBpcTau);
        HxScTauSB.Fill(xSBpcTau);
     } 
     if(xBpMass>XMinSBR&&xBpMass<XMaxSBR){ 
        HxcTauSB.Fill(xBpcTau);
        HxScTauSB.Fill(xSBpcTau);
     }  
    } 
   }
   char TXT [200];
   sprintf(TXT,"Mass        Entries = %7f",HxMass.GetEntries());
   cout<<"***********************************"<<endl;
   cout<<"***********************************\n"<<endl;
   cout<<"TauBpTree   Entries = "<<nentries<<endl;
   cout<<"Mass        Entries = "<<HxMass.GetEntries()<<endl;
   cout<<TXT<<endl;
   cout<<"SideBand    Entries = "<<HxcTauSB.GetEntries()<<endl;
   cout<<"\n***********************************"<<endl;
   cout<<"***********************************"<<endl;
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================
//  double arglis[10];
//  int ierflg= 0;
  model.setData(dataLife);
//  FitManager fitter(&modelC);
  GooFit::FitManagerMinuit1 fitter(&model);
  fitter.setMaxCalls(12000);
  cout<<"                  ===*** Start Fit ***=== "<<endl;
  cout<<"                  ===*** Start Fit ***=== "<<endl;
  cout<<"                  ===*** Start Fit ***=== "<<endl;
//  fitter.setupMinuit();
  Minuit1 * Minuit = fitter.getMinuitObject();
  Minuit->SetPrintLevel(1);
  Minuit->SetErrorDef(1);
//  Minuit->SetErrorDef(0.5);
  Double_t arglist[2]; int err = 0;
  arglist[0]= 12000; // maximum iterations
  arglist[1]= 1.0; 
  Minuit->mnexcm("MIGRAD",arglist,2,err);
  Minuit->mnexcm("HESSE",arglist, 0,err);

//  Minuit->SetErrorDef(0.5);

//   double arglis[10];
//   arglis[0]=8;
//   int ierflg= 0;
//   Minuit->mnexcm("FIX ",arglis,1,ierflg);
//   fitter.runCommand("MIGRAD");
//   for (int i=19;i<23;i++){
//    arglis[0]=i;
//    Minuit->mnexcm("FIX ",arglis,1,ierflg);
//   } 
//   fitter.runCommand("MIGRAD");
//   for (int i=1;i<37;i++){
//    arglis[0]=i;
//    Minuit->mnexcm("FIX ",arglis,1,ierflg);
//   }
//   for (int i=19;i<23;i++){
//    arglis[0]=i;
//    Minuit->mnexcm("REL ",arglis,1,ierflg); 
//   }
//   fitter.runCommand("MIGRAD");
//   for (int i=19;i<20;i++){
//    arglis[0]=i;
//    Minuit->mnexcm("REL ",arglis,1,ierflg); 
//   }
//      Minuit->FixParameter(0);
//      Minuit->FixParameter(1);
//      Minuit->FixParameter(2);
//      Minuit->FixParameter(3);
//      Minuit->FixParameter(4);
//      Minuit->FixParameter(5);
//      Minuit->FixParameter(6);

// 
//       Minuit->FixParameter(7);
// 
/* 
       Minuit->FixParameter(35);
       Minuit->FixParameter(36);
       Minuit->FixParameter(37);
       Minuit->FixParameter(38);
       Minuit->FixParameter(39);
       Minuit->FixParameter(40);
       Minuit->FixParameter(41);
       Minuit->FixParameter(42);
       Minuit->FixParameter(43);
       Minuit->FixParameter(44);
       Minuit->FixParameter(45);
       Minuit->FixParameter(46);
       Minuit->FixParameter(47);

 */
//      fitter.runCommand("MIGRAD");
   
//      Minuit->Release(7);
      
      
 //       Minuit->Release(36);
//        Minuit->Release(38);
//          Minuit->Release(39);
//          Minuit->Release(40);
//        Minuit->Release(41);
//        Minuit->Release(42);
//        Minuit->Release(43);
//        Minuit->Release(44);
//        Minuit->Release(45);
//        Minuit->Release(46);
// 
//
//       Minuit->Release(7);
//       fitter.runCommand("MIGRAD");

//      Minuit->FixParameter(18);
//      Minuit->FixParameter(19);
//      Minuit->FixParameter(20);
//      Minuit->FixParameter(21);
//Minuit->FixParameter(15);
//Minuit->FixParameter(16);
//Minuit->FixParameter(17);
//Minuit->FixParameter(18);
//     fitter.runCommand("MIGRAD"); 
//    for (int i=0;i<37;i++){
//     arglis[0]=i;
//     Minuit->mnexcm("FIX ",arglis,1,ierflg);
//    }
/*  for (int i=19;i<23;i++){
   arglis[0]=i;
   Minuit->mnexcm("REL ",arglis,1,ierflg);
  }
    arglis[0]=8;
      Minuit->mnexcm("REL ",arglis,1,ierflg);
*/
//    arglis[0]=8;
//    Minuit->mnexcm("REL ",arglis,1,ierflg);
//    arglis[0]=10;
//    Minuit->mnexcm("REL ",arglis,1,ierflg);
//    arglis[0]=16;
//    Minuit->mnexcm("REL ",arglis,1,ierflg);
//    arglis[0]=17;
//    Minuit->mnexcm("REL ",arglis,1,ierflg);
//     arglis[0]=34;
//     Minuit->mnexcm("REL ",arglis,1,ierflg);

//     Minuit->FixParameter(44);
//     Minuit->FixParameter(45);
//     Minuit->FixParameter(46);
//     Minuit->FixParameter(47);
    //fitter.runCommand("MIGRAD");
//    fitter.runMigrad(0.3);
 
//  fitter.runCommand("MINOS");
    //fitter.runCommand("HESSE");
//  fitter.fit();   
//  Minuit->SetPrintLevel(1);
//  Minuit->mnmigr();
//  Minuit->mnhess();
//  Minuit->mnmigr();
  std::vector<Variable> var; 
  double tmp_value, tmp_error;
  for(Variable &var : Minuit->getVaraibles()) {
      Minuit->GetParameter(var.getFitterIndex(), tmp_value, tmp_error);
      var.setValue(tmp_value);
      var.setError(tmp_error);
  }
//  fitter.getMinuitValues(); 
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================

//================================================================================
///PLOT

//XHScale=10;
// Mass
  UnbinnedDataSet gridMass(xMass);
  double totalDataMass = 0; 
  double NStepMass = XHScale*xMass.getNumBins();
  for (int i = 0; i < NStepMass; ++i) {
    double step = (xMass.getUpperLimit() - xMass.getLowerLimit())/NStepMass;
    xMass.setValue(xMass.getLowerLimit() + (i + 0.5) * step);
    gridMass.addEvent(); 
   totalDataMass++; 
  }

  modelMass.setData(&gridMass);
  std::vector<std::vector<double> > pdfValsMass = modelMass.getCompProbsAtDataPoints();
//  modelMass.getCompProbsAtDataPoints(pdfValsMass); 
  double totalPdfMass = 0; 
  for (int i = 0; i < gridMass.getNumEvents(); ++i) {
    gridMass.loadEvent(i); 
    pdfHist.Fill(xMass.getValue(), pdfValsMass[0][i]);
    sigHist.Fill(xMass.getValue(), pdfValsMass[1][i]);
    bkgHist.Fill(xMass.getValue(), pdfValsMass[2][i]);
    totalPdfMass += pdfValsMass[0][i]; 
  }
  
  
  pdfHist.Scale((signalYield->getValue()+bckgYield->getValue())/pdfHist.Integral()*XHScale);
  sigHist.Scale(signalYield->getValue()/sigHist.Integral()*XHScale);
  bkgHist.Scale(bckgYield->getValue()/bkgHist.Integral()*XHScale);
  std::cout<<"Signal Yield = "<< signalYield->getValue()<<std::endl;
  std::cout<<"Bckg   Yield = "<< bckgYield->getValue()<<std::endl;
  std::cout<<"(SB    Yield  = "<<HxcTauSB.GetEntries() <<")"<<std::endl;
  std::cout<<"Tot   Yield  = "<< signalYield->getValue()+bckgYield->getValue()<<std::endl;
//--------------------------------------------------  
// Tau


//XHScale=1;
  int NIntegral = 1;
  

//  vector<GooFit::Variable*> dataPlot;
//  dataPlot.push_back(xcTau);
//  dataPlot.push_back(xScTau);

  std::vector<GooFit::Observable> dataPlot2D;
  dataPlot2D.push_back(xcTau);
  dataPlot2D.push_back(xScTau);

//  vector<GooFit::Variable*> dataPlotS;
//  dataPlotS.push_back(xScTau);
  
  
//  UnbinnedDataSet grid_cTau(dataPlot);
  UnbinnedDataSet grid_cTau2D(dataPlot2D);
//  UnbinnedDataSet grid_STau(dataPlotS);
  
//  bool first = true;
//  UnbinnedDataSet grid_cTau(xcTau);
//  double totalData_cTau = 0; 
  double NStep  = XHScale*xcTau.getNumBins();
//  double NSStep   = XHScale*xcTau->getNumBins();
  double NSStep2D = NIntegral*XHScale*xScTau.getNumBins();
  double step  = (xcTau.getUpperLimit() - xcTau.getLowerLimit())/NStep;
//  double sstep = (xScTau.getUpperLimit() - xScTau.getLowerLimit())/NSStep;
  double sstep2D = (xScTau.getUpperLimit() - xScTau.getLowerLimit())/NSStep2D;
  for (int i = 0; i < NStep; ++i) {
    xcTau.setValue(xcTau.getLowerLimit()  + (i + 0.5) * step);
//    grid_cTau.addEvent(); 
//    totalData_cTau++; 
//    xScTau.getValue() = xScTau.getLowerLimit() + (i + 0.5) * sstep;
//    cout<<"X = "<<xcTau->getValue()<<" sx = "<<xScTau.getValue()<<endl;
//    grid_cTau.addEvent(); 
//    xcTau2D->getValue()  = xcTau2D ->getLowerLimit() + (i + 0.5) * step;
//   cout<<"======================================     \n"<<endl;
//   cout<<"======================================     \n"<<endl;
//   cout<<"======================================     \n"<<endl;
    for (int ii = 0; ii < NSStep2D; ++ii) {
     xScTau.setValue(xScTau.getLowerLimit() + (ii + 0.5) * sstep2D);
//     xScTau2D->getValue() = xScTau2D->getLowerLimit() + (ii + 0.5) * sstep2D;
//    cout<<"X = "<<xcTau->getValue()<<" sx = "<<xScTau.getValue()<<endl;
     grid_cTau2D.addEvent(); 
//     if (first) grid_STau.addEvent();
    }
//    first = false;
  }

//  model_cTau.setData(&grid_cTau);
//  vector<vector<double> > pdfVals_cTau;
//  model_cTau.getCompProbsAtDataPoints(pdfVals_cTau); 
//  double totalPdf_cTau = 0; 
//   for (int i = 0; i < grid_cTau.getNumEvents(); ++i) {
//     grid_cTau.loadEvent(i); 
//     pdf_cTau_Hist.Fill(xcTau->getValue() , pdfVals_cTau[0][i]);
//     sig_cTau_Hist.Fill(xcTau->getValue() , pdfVals_cTau[1][i]);
//     bkg_cTau_Hist.Fill(xcTau->getValue() , pdfVals_cTau[2][i]);
//     totalPdf_cTau += pdfVals_cTau[0][i]; 
//   }

//  double pdf_cTau_Integral2D = 0;
//  double sig_cTau_Integral2D = 0;
//  double bkg_cTau_Integral2D = 0;
//model.setData(&grid_cTau2D);
//vector<vector<double> > pdfVals_cTau2D;
//model.getCompProbsAtDataPoints(pdfVals_cTau2D);
     model_cTau.setData(&grid_cTau2D);
     vector<vector<double> > pdfVals_cTau2D = model_cTau.getCompProbsAtDataPoints();
//     model_cTau.getCompProbsAtDataPoints(pdfVals_cTau2D);
  for (int i = 0; i < grid_cTau2D.getNumEvents(); ++i) {
    grid_cTau2D.loadEvent(i); 
    pdf_cTauSTau_Hist2D.Fill(xcTau.getValue() ,xScTau.getValue() , pdfVals_cTau2D[0][i]);
    sig_cTauSTau_Hist2D.Fill(xcTau.getValue() ,xScTau.getValue() , pdfVals_cTau2D[1][i]);
    bkg_cTauSTau_Hist2D.Fill(xcTau.getValue() ,xScTau.getValue() , pdfVals_cTau2D[2][i]);
//     if (i%int(NSStep2D) == 1 && i>0){
//      pdf_cTau_Hist2D.Fill(xcTau->getValue() , pdf_cTau_Integral2D/step);
//      sig_cTau_Hist2D.Fill(xcTau->getValue() , sig_cTau_Integral2D/step);
//      bkg_cTau_Hist2D.Fill(xcTau->getValue() , bkg_cTau_Integral2D/step);
// //     cout<<"Int = "<<bkg_cTau_Integral2D<<endl;
//      pdf_cTau_Integral2D=0;
//      sig_cTau_Integral2D=0;
//      bkg_cTau_Integral2D=0;
// //     exit(0);
//     }else{
// //     cout<<"X = "<<xcTau->getValue()<<" sx = "<<xScTau.getValue()<<" NStep2D = "<<NStep2D<<endl;
//      pdf_cTau_Integral2D =+ pdfVals_cTau2D[0][i]*sstep2D; 
//      sig_cTau_Integral2D =+ pdfVals_cTau2D[1][i]*sstep2D; 
//      bkg_cTau_Integral2D =+ pdfVals_cTau2D[2][i]*sstep2D;
// //     cout<<"Int = "<<pdfVals_cTau2D[2][i]<<endl;
//     } 
  }
  
  TH1D * pdf_cTauSTau_X = pdf_cTauSTau_Hist2D.ProjectionX("pdf_cTauSTau_X");
  TH1D * pdf_cTauSTau_Y = pdf_cTauSTau_Hist2D.ProjectionY("pdf_cTauSTau_Y");

  TH1D * sig_cTauSTau_X = sig_cTauSTau_Hist2D.ProjectionX("sig_cTauSTau_X");
  TH1D * sig_cTauSTau_Y = sig_cTauSTau_Hist2D.ProjectionY("sig_cTauSTau_Y");

  TH1D * bkg_cTauSTau_X = bkg_cTauSTau_Hist2D.ProjectionX("bkg_cTauSTau_X");
  TH1D * bkg_cTauSTau_Y = bkg_cTauSTau_Hist2D.ProjectionY("bkg_cTauSTau_Y");

/*   vector<GooFit::Variable*> dataSPlot;
  dataSPlot.push_back(xScTau);
  dataSPlot.push_back(xcTau);
  UnbinnedDataSet grid_STau(dataSPlot);
//  UnbinnedDataSet grid_cTau(xcTau);
  double totalData_STau = 0; 
  NStep  = XHScale*xScTau.getNumBins();
  double NSStep = XHScale*xcTau->getNumBins();
  for (int i = 0; i < NSStep; ++i) {
    totalData_STau++; 
    double sstep = (xScTau.getUpperLimit() - xScTau.getLowerLimit())/NSStep;
    xScTau.getValue() = xScTau.getLowerLimit() + (i + 0.5) * sstep;
    grid_STau.addEvent(); 
  }
 */

//
// STau
//
//   model_STau.setData(&grid_STau);
//   vector<vector<double> > pdfVals_STau;
//   model_STau.getCompProbsAtDataPoints(pdfVals_STau); 
//   double totalPdf_STau = 0;  
//   for (int i = 0; i < grid_STau.getNumEvents(); ++i) {
//     grid_STau.loadEvent(i); 
//     pdf_STau_Hist.Fill(xScTau.getValue(), pdfVals_STau[0][i]);
//     sig_STau_Hist.Fill(xScTau.getValue(), pdfVals_STau[1][i]);
//     bkg_STau_Hist.Fill(xScTau.getValue(), pdfVals_STau[2][i]);
//     totalPdf_STau += pdfVals_STau[0][i]; 
//   }
  
//
// Models plot  
//   pdf_cTau_Hist.Scale((signalYield->getValue()+bckgYield->getValue())/pdf_cTau_Hist.Integral()*XHScale);
//   sig_cTau_Hist.Scale(signalYield->getValue()/sig_cTau_Hist.Integral()*XHScale);
//   bkg_cTau_Hist.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist.Integral()*XHScale);
// 
//   pdf_cTau_Hist2D.Scale((signalYield->getValue()+bckgYield->getValue())/pdf_cTau_Hist2D.Integral()*XHScale);
//   sig_cTau_Hist2D.Scale(signalYield->getValue()/sig_cTau_Hist2D.Integral()*XHScale);
//   bkg_cTau_Hist2D.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist2D.Integral()*XHScale);
    
  pdf_cTauSTau_X->Scale((signalYield->getValue()+bckgYield->getValue())/pdf_cTauSTau_X->Integral()*XHScale);
  sig_cTauSTau_X->Scale((signalYield->getValue())/sig_cTauSTau_X->Integral()*XHScale);
  bkg_cTauSTau_X->Scale((HxcTauSB.GetEntries())/bkg_cTauSTau_X->Integral()*XHScale);
//  bkg_cTauSTau_X->Scale((bckgYield->getValue())/bkg_cTauSTau_X->Integral()*XHScale);
  pdf_cTauSTau_Y->Scale((signalYield->getValue()+bckgYield->getValue())/pdf_cTauSTau_Y->Integral()*XHScale);
  sig_cTauSTau_Y->Scale((signalYield->getValue())/sig_cTauSTau_Y->Integral()*XHScale);
  bkg_cTauSTau_Y->Scale((HxcTauSB.GetEntries())/bkg_cTauSTau_Y->Integral()*XHScale);
//  bkg_cTauSTau_Y->Scale((bckgYield->getValue())/bkg_cTauSTau_Y->Integral()*XHScale);
  
  

//   sig_STau_Hist.Scale(signalYield->getValue()/sig_STau_Hist.Integral()*XHScale);
//   bkg_STau_Hist.Scale(HxScTauSB.GetEntries()/bkg_STau_Hist.Integral()*XHScale);
//   pdf_STau_Hist.Scale((signalYield->getValue()+bckgYield->getValue())/pdf_STau_Hist.Integral()*XHScale);
//  bkg_cTau_Hist.Scale(bckgYield->getValue()/bkg_cTau_Hist.Integral()*XHScale);
 
//   for (int i = 0; i < xMass.getNumBins(); ++i) {
//     double val = pdfHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= totalData;
//     pdfHist.SetBinContent(i+1, val); 
//     val = sigHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= sigFrac->getValue(); 
//     val *= totalData;
//     sigHist.SetBinContent(i+1, val); 
//     val = bkgHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= (1.0 - sigFrac->getValue());
//     val *= totalData;
//     bkgHist.SetBinContent(i+1, val); 
//   }

//  double sigmaw    = sigma1.getValue()*wg1->getValue()+ (1-wg1->getValue())*sigma2.getValue();
//  double sigmawErr = sqrt(sigma1.getError()*wg1->getValue()*sigma1.getError()*wg1->getValue()+ (1-wg1->getValue())*sigma2.getError()*(1-wg1->getValue())*sigma2.getError());
  c1->cd();
    TLegend* leg_sign = new TLegend(0.30,0.70,0.90,0.90);
    leg_sign->SetTextSize(0.025) ;
    leg_sign->SetTextAlign(31);
    leg_sign->SetBorderSize(0.);
    leg_sign->SetFillStyle(0);
    leg_sign->SetHeader("B^{+} mass spectrum  Fit Projection");
    if(signalYield->getError()!=0){
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f  #pm %5.0f",signalYield->getValue(),signalYield->getError()),"");
    }else{
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f Fixed",signalYield->getValue()),"");
    }
    if(bckgYield->getError()!=0){
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =    %5.0f  #pm  %5.0f",bckgYield->getValue(),bckgYield->getError()),"");
    }else{
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =    %5.0f  Fixed",bckgYield->getValue()),"");
    }
    
    if(mean.getError()!=0){
     leg_sign->AddEntry(&HxMass ,Form( "M_{B^{+}} =   %5.5f  #pm %5.5f",mean.getValue(),mean.getError()),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "M_{B^{+}} =   %5.5f Fixed",mean.getValue()),"");
     }
    if(sigma1.getError()!=0){
     leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{1}_{B^{+}} =   %5.5f  #pm %5.5f",sigma1.getValue(),sigma1.getError()),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{1}_{B^{+}} =   %5.5f Fixed",sigma1.getValue()),"");
    }
    if(sigma2.getError()!=0){
     leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{2}_{B^{+}} =   %5.5f  #pm %5.5f",sigma2.getValue(),sigma2.getError()),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{2}_{B^{+}} =   %5.5f Fixed",sigma2.getValue()),"");
    }
  HxMass.GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  HxMass.SetMarkerStyle(8);
  HxMass.SetMarkerSize(MarkerSize);
  HxMass.SetTitle("");
  HxMass.Draw("E1"); 
//  HxMass.Draw("p"); 
  pdfHist.SetLineWidth(PlotLineWidth); 
  pdfHist.SetLineColor(kBlue);
  pdfHist.Draw("same,HIST"); 
  sigHist.SetLineWidth(PlotLineWidth); 
  sigHist.SetLineColor(kMagenta);
  sigHist.SetLineStyle(kDashed); 
  sigHist.Draw("same,HIST"); 
  bkgHist.SetLineWidth(PlotLineWidth); 
  bkgHist.SetLineColor(kRed);
  bkgHist.SetLineStyle(kDashed); 
  bkgHist.Draw("same,HIST"); 
  leg_sign->Draw("same");
  HxMass.Write();
  pdfHist.Write();
  sigHist.Write();
  bkgHist.Write();
//  
  c2->cd();
  c2->SetLogy();
  TLegend* leg_pdfSB = new TLegend(0.50,0.65,0.90,0.90);
  leg_pdfSB->SetTextAlign(12);
  leg_pdfSB->SetHeader("B^{+} proper time Fit Projection");
  leg_pdfSB->SetTextSize(0.025) ;
  leg_pdfSB->SetBorderSize(0.);
  leg_pdfSB->SetFillStyle(0);
  leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[4]{#scale[1.5]{#tau}_{B^{+}}  =  %5.4f #pm %5.4f [ps]}",1/(c_const*cTau.getValue()),cTau.getError()/((c_const*cTau.getValue())*(cTau.getValue())))   ,"");
  if( b1.getError()!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b1   =  %5.4f #pm %5.4f     ",b1.getValue(),b1.getError())   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b1   =  %5.4f     Fixed     ",b1.getValue())   ,"");
  }   
//   if( b2->getError()!=0){
//       leg_pdfSB->AddEntry(&HxcTau ,Form( "b2   =  %5.3f #pm %5.3f     ",b2->getValue(),b2->getError())   ,"");
//   }else{      
//       leg_pdfSB->AddEntry(&HxcTau ,Form( "b2   =  %5.3f     Fixed     ",b2->getValue())   ,"");
//   }   
  if( tauSB1.getError()!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.4f #pm %5.4f     }",1/(tauSB1.getValue()),tauSB1.getError()/((tauSB1.getValue())*(tauSB1.getValue())))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.4f	 Fixed	  }",1/(tauSB1.getValue()))   ,"");
  }   
  if( tauSB2.getError()!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.4f #pm %5.4f     }",1/(tauSB2.getValue()),tauSB2.getError()/((tauSB2.getValue())*(tauSB2.getValue())))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.4f	 Fixed	  }",1/(tauSB2.getValue()))   ,"");
  }   
//   if( tauSB3->getError()!=0){
//       leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f #pm %5.3f     }",1/(tauSB3->getValue()),tauSB3->getError()/((tauSB3->getValue())*(tauSB3->getValue())))   ,"");
//   }else{      
//       leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f Fixed     }",1/(tauSB3->getValue()))   ,"");
//   }   
  leg_pdfSB->AddEntry(&HxcTau ,"#color[4]{#scale[0.8]{- Fit model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[6]{#scale[0.8]{- Signal model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[2]{#scale[0.8]{- Background model on SB}}"   ,"");
  HxcTau.SetMinimum(10);
  HxcTauSB.SetMinimum(10);
  HxcTau.GetXaxis()->SetTitle("ct (cm)");
  HxcTau.SetMarkerStyle(8);
  HxcTau.SetMarkerSize(MarkerSize);
  HxcTau.SetTitle("");
  HxcTau.Draw("E1");
  HxcTauSB.SetMarkerStyle(8);
  HxcTauSB.SetMarkerSize(0.5);
  HxcTauSB.SetMarkerColor(kRed);
  HxcTauSB.Draw("same,E1");
  leg_pdfSB->Draw("same");
/*   pdf_cTau_Hist2D.SetLineColor(kBlue);
  pdf_cTau_Hist2D.SetLineWidth(3); 
  pdf_cTau_Hist2D.Draw("same"); 
  sig_cTau_Hist2D.SetLineColor(kMagenta);
  sig_cTau_Hist2D.SetLineStyle(kDashed); 
  sig_cTau_Hist2D.SetLineWidth(2); 
  sig_cTau_Hist2D.Draw("same"); 
  bkg_cTau_Hist2D.SetLineColor(kRed);
  bkg_cTau_Hist2D.SetLineStyle(kDashed); 
  bkg_cTau_Hist2D.SetLineWidth(2); 
  bkg_cTau_Hist2D.Draw("same"); 
 */  
  pdf_cTauSTau_X->SetLineColor(kBlue);
  pdf_cTauSTau_X->SetLineWidth(PlotLineWidth);
  pdf_cTauSTau_X->Draw("same,HIST");
  sig_cTauSTau_X->SetLineColor(kMagenta);
  sig_cTauSTau_X->SetLineWidth(PlotLineWidth);
  sig_cTauSTau_X->SetLineStyle(kDashed);
  sig_cTauSTau_X->Draw("same,HIST");
  bkg_cTauSTau_X->SetLineColor(kRed);
  bkg_cTauSTau_X->SetLineWidth(PlotLineWidth);
  bkg_cTauSTau_X->SetLineStyle(kDashed);
  bkg_cTauSTau_X->Draw("same,HIST");
//   pdf_cTau_Hist.SetLineColor(kBlue);
//   pdf_cTau_Hist.SetLineWidth(3); 
//   pdf_cTau_Hist.Draw("same"); 
//   sig_cTau_Hist.SetLineColor(kMagenta);
//   sig_cTau_Hist.SetLineStyle(kDashed); 
//   sig_cTau_Hist.SetLineWidth(2); 
//   sig_cTau_Hist.Draw("same"); 
//   bkg_cTau_Hist.SetLineColor(kRed);
//   bkg_cTau_Hist.SetLineStyle(kDashed); 
//   bkg_cTau_Hist.SetLineWidth(2); 
//   bkg_cTau_Hist.Draw("same"); 
  
  
  c3->cd();
  TLegend* leg_pdfUncertainty = new TLegend(0.40,0.47,0.90,0.90);
  leg_pdfUncertainty->SetTextAlign(12);
  leg_pdfUncertainty->SetHeader("B^{+} Uncertainty Fit Projection");
  leg_pdfUncertainty->SetTextSize(0.025) ;
  leg_pdfUncertainty->SetBorderSize(0.);
  leg_pdfUncertainty->SetFillStyle(0);
//  leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "#color[4]{#scale[1.5]{#tau}_{B^{+}}  =  %5.3f #pm %5.3f     }",1/(c_const*cTau->getValue()),cTau->getError()/((c_const*cTau->getValue())*(cTau->getValue())))   ,"");
//   if( b1->getError()!=0){
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "b1   =  %5.3f #pm %5.3f     ",b1->getValue(),b1->getError())   ,"");
//   }else{      
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "b1   =  %5.3f     Fixed     ",b1->getValue())   ,"");
//   }   
//   if( tauUncertainty1->getError()!=0){
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Uncertainty1} =  %5.3f #pm %5.3f     }",1/(c_const*tauUncertainty1->getValue()),tauUncertainty1->getError()/((c_const*tauUncertainty1->getValue())*(tauUncertainty1->getValue())))   ,"");
//   }else{      
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Uncertainty1} =  %5.3f	Fixed	  }",1/(c_const*tauUncertainty1->getValue()))   ,"");
//   }   
//   if( tauUncertainty2->getError()!=0){
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Uncertainty2} =  %5.3f #pm %5.3f     }",1/(c_const*tauUncertainty2->getValue()),tauUncertainty2->getError()/((c_const*tauUncertainty2->getValue())*(tauUncertainty2->getValue())))   ,"");
//   }else{      
//       leg_pdfUncertainty->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Uncertainty2} =  %5.3f	Fixed	  }",1/(c_const*tauUncertainty2->getValue()))   ,"");
//   }   
  leg_pdfUncertainty->AddEntry(&HxScTau ,"#color[4]{#scale[0.8]{- Pdf model	     }}"   ,"");
  leg_pdfUncertainty->AddEntry(&HxScTau ,"#color[6]{#scale[0.8]{- Signal model   }}"   ,"");
  leg_pdfUncertainty->AddEntry(&HxScTau ,"#color[2]{#scale[0.8]{- Background model (on SB events)}}"   ,"");
  if( tauErrSign1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#tau}_{ErrSign }   =  %5.3e #pm %5.3e     }",1/tauErrSign1.getValue(),tauErrSign1.getError()/(tauErrSign1.getValue()*tauErrSign1.getValue()))   ,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#tau}_{ErrSign }   =  %5.3e   Fixed   }",1/tauErrSign1.getValue())   ,"");
  }   
  if( meanGaussianErrSign1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrSign1}  =  %5.3e #pm %5.3e     }",meanGaussianErrSign1.getValue(),meanGaussianErrSign1.getError())   ,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrSign1}  =  %5.3e   Fixed     }",meanGaussianErrSign1.getValue())   ,"");
  }   
  if( sigmaGaussianErrorSign1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrSign1}  =  %5.3e #pm %5.3e     }",sigmaGaussianErrorSign1.getValue(),sigmaGaussianErrorSign1.getError())	,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrSign1}  =  %5.3e	Fixed	}",sigmaGaussianErrorSign1.getValue())   ,"");
  }   
  if( meanGaussianErrSign2.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrSign2}  =  %5.3e #pm %5.3e     }",meanGaussianErrSign2.getValue(),meanGaussianErrSign2.getError())   ,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrSign2}  =  %5.3e   Fixed     }",meanGaussianErrSign2.getValue())   ,"");
  }   
  if( sigmaGaussianErrorSign2.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrSign2}  =  %5.3e #pm %5.3e     }",sigmaGaussianErrorSign2.getValue(),sigmaGaussianErrorSign2.getError())	,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrSign2}  =  %5.3e	Fixed	}",sigmaGaussianErrorSign2.getValue())   ,"");
  }   
  if( tauErrBckg1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#tau}_{ErrBckg }   =  %5.3e #pm %5.3e     }",1/tauErrBckg1.getValue(),tauErrBckg1.getError()/(tauErrBckg1.getValue()*tauErrBckg1.getValue()))   ,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#tau}_{ErrBckg }   =  %5.3e   Fixed   }",1/tauErrBckg1.getValue())   ,"");
  }   
  if( meanGaussianErrBckg1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrBckg1}  =  %5.3e #pm %5.3e     }",meanGaussianErrBckg1.getValue(),meanGaussianErrBckg1.getError())   ,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.0]{M}_{ErrBckg1}  =  %5.3e   Fixed     }",meanGaussianErrBckg1.getValue())   ,"");
  }   
  if( sigmaGaussianErrorBckg1.getError()!=0){
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrBckg1}  =  %5.3e #pm %5.3e     }",sigmaGaussianErrorBckg1.getValue(),sigmaGaussianErrorBckg1.getError())	,"");
  }else{      
      leg_pdfUncertainty->AddEntry(&HxScTau ,Form( "#color[2]{#scale[1.5]{#sigma}_{ErrBckg1}  =  %5.3e	Fixed	}",sigmaGaussianErrorBckg1.getValue())   ,"");
  }   
  HxScTau.Draw("E1"); 
  HxScTau.GetXaxis()->SetTitle("ct (cm)");
  HxScTau.SetTitle("");
  leg_pdfUncertainty->Draw("same");
  HxScTauSB.SetMarkerStyle(8);
  HxScTauSB.SetMarkerSize(MarkerSize);
  HxScTauSB.SetMarkerColor(kRed);
  HxScTauSB.Draw("same,E1");
  pdf_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  pdf_cTauSTau_Y->SetLineColor(kBlue);
  pdf_cTauSTau_Y->Draw("same,HIST");
  sig_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  sig_cTauSTau_Y->SetLineColor(kMagenta);
  sig_cTauSTau_Y->SetLineStyle(kDashed);
  sig_cTauSTau_Y->Draw("same,HIST");
  bkg_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  bkg_cTauSTau_Y->SetLineColor(kRed);
  bkg_cTauSTau_Y->SetLineStyle(kDashed);
  bkg_cTauSTau_Y->Draw("same,HIST");


//   pdf_STau_Hist.SetLineColor(kBlue);
//   pdf_STau_Hist.SetLineWidth(2); 
//   pdf_STau_Hist.Draw("same"); 
//   sig_STau_Hist.SetLineColor(kMagenta);
//   sig_STau_Hist.SetLineStyle(kDashed); 
//   sig_STau_Hist.SetLineWidth(2); 
//   sig_STau_Hist.Draw("same"); 
//   bkg_STau_Hist.SetLineColor(kRed);
//   bkg_STau_Hist.SetLineStyle(kDashed); 
//   bkg_STau_Hist.SetLineWidth(2); 
//   bkg_STau_Hist.Draw("same"); 
  
  HxcTau.Write();
  HxcTauSB.Write();
  HxScTau.Write();
  HxScTauSB.Write();
  pdf_cTauSTau_Hist2D.Write();
  sig_cTauSTau_Hist2D.Write();
  bkg_cTauSTau_Hist2D.Write();
//   pdf_cTau_Hist2D.Write();
//   sig_cTau_Hist2D.Write();
//   bkg_cTau_Hist2D.Write();
//   pdf_cTau_Hist.Write();
//   sig_cTau_Hist.Write();
//   bkg_cTau_Hist.Write();
//   pdf_STau_Hist.Write();
//   sig_STau_Hist.Write();
//   bkg_STau_Hist.Write();
  pdf_cTauSTau_X->Write();
  pdf_cTauSTau_Y->Write();
  sig_cTauSTau_X->Write();
  sig_cTauSTau_Y->Write();
  bkg_cTauSTau_X->Write();
  bkg_cTauSTau_Y->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  char PDFNameMass[50] = "Bp-Mass-2016.pdf";
  char PDFNamecTau[50] = "Bp-cTau-2016.pdf";
  char PDFNameReso[50] = "Bp-Reso-2016.pdf";
  char testo[130] ;
  sprintf(testo,"mv %s %s.tmp",PDFNameMass,PDFNameMass);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNamecTau,PDFNamecTau);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNameReso,PDFNameReso);
  gSystem->Exec(testo);
  
  std::cout<<"Tau  [ps]= "<<1/(c_const*cTau.getValue())<<"+/-"<<cTau.getError()/((c_const*cTau.getValue())*(cTau.getValue()))<<std::endl;
  std::cout<<"1/TauSB1 = "<<1/(tauSB1.getValue())<<"+/-"<<tauSB1.getError()/((tauSB1.getValue())*(tauSB1.getValue()))<<std::endl;
  std::cout<<"1/TauSB2 = "<<1/(tauSB2.getValue())<<"+/-"<<tauSB2.getError()/((tauSB2.getValue())*(tauSB2.getValue()))<<std::endl;
//  std::cout<<"1/TauSB3 = "<<1/(tauSB3->getValue())<<"+/-"<<tauSB3->getError()/((tauSB3->getValue())*(tauSB3->getValue()))<<std::endl;

  c1->Print(PDFNameMass);
  c2->Print(PDFNamecTau);
  c3->Print(PDFNameReso);
  OutFile->Close();
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;

  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;
}
