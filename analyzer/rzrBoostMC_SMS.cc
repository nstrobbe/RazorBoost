//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmd.h"
#include "utils.h"
#include "btagutils.h"
#include <math.h>

#include "TLorentzVector.h"

#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/pdg.h"
#else
#include "pdg.h"
#endif

using namespace std;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{

  // Get the trigger histogram:
  TFile* fhlt = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  if (!fhlt){
    fhlt = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  }
  if (!fhlt){
    cout << "Could not find trigger efficiency root file... Where did you put it??" << endl;
    return 1;
  }
  //TH2D* h_hlteff = (TH2D*)fhlt->Get("hBinValues");
  TH2D* h_hlteff = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph");
    
  // Get file list and histogram filename from command line
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples
  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read
  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

  // Select variables to be read
  selectVariables(stream);

  /*
	 Notes:
	
	 1. Use
	   ofile = outputFile(cmdline.outputfile, stream)
	
	 to skim events to output file in addition to writing out histograms.
	
	 2. Use
	   ofile.addEvent(event-weight)
	
	 to specify that the current event is to be added to the output file.
	 If omitted, the event-weight is defaulted to 1.
	
	 3. Use
	    ofile.count(cut-name, event-weight)
	
	 to keep track, in the count histogram, of the number of events
	 passing a given cut. If omitted, the event-weight is taken to be 1.
	 If you want the counts in the count histogram to appear in a given
	 order, specify the order, before entering the event loop, as in
	 the example below
	 
	    ofile.count("NoCuts", 0)
		ofile.count("GoodEvent", 0)
		ofile.count("Vertex", 0)
		ofile.count("MET", 0)
  */
  
  string fsample = cmdline.filelist;
  double xsect = cmdline.xsect;
  // ! totweight should contain the proper ISR weights if we want to do ISRreweighting !
  // ! totweight should contain the proper top pT weights if we want to do top Pt reweighting !
  double totweight = cmdline.totweight; 
  double lumi = cmdline.lumi;

  string sample = "";
  if ( argc > 6 )
    sample = string(argv[6]);
  string ISR = "";
  if ( argc > 7 )
    ISR = string(argv[7]);
  string TopPt = "";
  if ( argc > 8 )
    TopPt = string(argv[8]);
  string Pileup = "";
  if ( argc > 9 )
    Pileup = string(argv[9]);

  bool doISRreweighting = false;
  if (ISR == "ISR_True" 
      && (sample == "T2tt" || sample == "T1ttcc" || sample == "T1ttcc_old" || sample == "T1t1t"
	  || sample == "TTJets" || sample == "WJets" || sample == "ZJets" )
      ){
    doISRreweighting = true;
    cout << "Will do ISR reweighting" << endl;
  }
  bool doTopPtreweighting = false;
  if (sample == "TTJets" && TopPt == "TopPt_True"){
    doTopPtreweighting = true;
    doISRreweighting = false;
    cout << "Will do top pt reweighting" << endl;
  }
  bool doPileupReweighting = false;
  if (sample != "Data" && Pileup == "Pileup_True"){
    doPileupReweighting = true;
    cout << "Will do pileup reweighting" << endl;
  }

  // Get the pileup histogram:
  TString pileupname = "pileup_weights.root";
  if (sample == "T1ttcc_old" || sample == "T2tt" ){
    pileupname = "pileup_weights_sig52X.root";
  }
  if (sample == "T1ttcc" || sample == "T1t1t"){
    pileupname = "pileup_weights_sig53X.root";
  }
  TFile* fpileup = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname);
  if (!fpileup){
    fpileup = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/pileup/"+pileupname);
  }
  if (!fpileup){
    cout << "Could not find pileup weights root file... Where did you put it??" << endl;
    return 1;
  }
  TH1D* h_pileup = (TH1D*)fpileup->Get("pileup_weight");

  // ----------------------------------
  // -- Get the btag eff histograms: --
  // ----------------------------------

  TFile* fbeff;

  if (sample == "TTJets" or sample == "Top" or sample == "TTX") {
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_TTJets.root");
  } else if (sample == "T1ttcc" || sample == "T2tt") {
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_T1ttcc.root");
  } else {
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_QCD.root");
  }
  TH1D* h_pt_b_CSVMeff = (TH1D*)fbeff->Get("h_pt_b_CSVMeff");
  TH1D* h_pt_c_CSVMeff = (TH1D*)fbeff->Get("h_pt_c_CSVMeff");
  //TH1D* h_pt_l_CSVMeff = (TH1D*)fbeff->Get("h_pt_l_CSVMeff");
  TH1D* h_pt_lc_CSVMeff = (TH1D*)fbeff->Get("h_pt_lc_CSVMeff");

  TH1D* h_pt_b_CSVLeff = (TH1D*)fbeff->Get("h_pt_b_CSVLeff");
  TH1D* h_pt_c_CSVLeff = (TH1D*)fbeff->Get("h_pt_c_CSVLeff");
  //TH1D* h_pt_l_CSVLeff = (TH1D*)fbeff->Get("h_pt_l_CSVLeff");
  TH1D* h_pt_lc_CSVLeff = (TH1D*)fbeff->Get("h_pt_lc_CSVLeff");

  //TH2D* h_pt_eta_b_CSVMeff = (TH2D*)fbeff->Get("h_pt_eta_b_CSVMeff");
  //TH2D* h_pt_eta_c_CSVMeff = (TH2D*)fbeff->Get("h_pt_eta_c_CSVMeff");
  TH2D* h_pt_eta_l_CSVMeff = (TH2D*)fbeff->Get("h_pt_eta_l_CSVMeff");
  //TH2D* h_pt_eta_lc_CSVMeff = (TH2D*)fbeff->Get("h_pt_eta_lc_CSVMeff");

  //TH2D* h_pt_eta_b_CSVLeff = (TH2D*)fbeff->Get("h_pt_eta_b_CSVLeff");
  //TH2D* h_pt_eta_c_CSVLeff = (TH2D*)fbeff->Get("h_pt_eta_c_CSVLeff");
  TH2D* h_pt_eta_l_CSVLeff = (TH2D*)fbeff->Get("h_pt_eta_l_CSVLeff");
  //TH2D* h_pt_eta_lc_CSVLeff = (TH2D*)fbeff->Get("h_pt_eta_lc_CSVLeff");

  outputFile ofile(cmdline.outputfilename);

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------

  TH1::SetDefaultSumw2();

  // set the ranges and nbins according to which sample we are using
  // declare variables and give some random default values
  int nbins_stop = 10;
  int nbins_LSP = 10;
  int stop_min = 0; 
  int stop_max = 1000; 
  int LSP_min = 0; 
  int LSP_max = 500; 

  if (sample == "T1ttcc" || sample == "T1t1t"){
    // stop is actually gluino in this case
    nbins_stop = 15;
    nbins_LSP = 11;
    stop_min = 600; 
    stop_max = 1350; 
    LSP_min = 0; 
    LSP_max = 550; 
  } else if (sample == "T1ttcc_old"){
    nbins_stop = 113;
    nbins_LSP = 6;
    stop_min = 310; 
    stop_max = 880; 
    LSP_min = 300; 
    LSP_max = 900; 
  } else if (sample == "T2tt"){
    nbins_stop = 35;
    nbins_LSP = 37;
    stop_min = 150; 
    stop_max = 1025; 
    LSP_min = 0; 
    LSP_max = 925; 
  }

  TH1D* h_mstop = new TH1D("h_mstop","h_mstop",nbins_stop,stop_min,stop_max);
  TH1D* h_mLSP = new TH1D("h_mLSP","h_mLSP",nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP = new TH2D("h_mstop_mLSP","h_mstop_mLSP",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_nevents = new TH2D("h_mstop_mLSP_nevents","h_mstop_mLSP_nevents",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_NoCuts = new TH2D("h_mstop_mLSP_NoCuts","h_mstop_mLSP_NoCuts",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_Cleaning = new TH2D("h_mstop_mLSP_Cleaning","h_mstop_mLSP_Cleaning",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_Pileup = new TH2D("h_mstop_mLSP_Pileup","h_mstop_mLSP_Pileup",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_ISR = new TH2D("h_mstop_mLSP_ISR","h_mstop_mLSP_ISR",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_HCAL_noise = new TH2D("h_mstop_mLSP_HCAL_noise","h_mstop_mLSP_HCAL_noise",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_vertexg0 = new TH2D("h_mstop_mLSP_vertexg0","h_mstop_mLSP_vertexg0",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_njetge3 = new TH2D("h_mstop_mLSP_njetge3","h_mstop_mLSP_njetge3",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_HLT = new TH2D("h_mstop_mLSP_HLT","h_mstop_mLSP_HLT",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_jet1ptg200 = new TH2D("h_mstop_mLSP_jet1ptg200","h_mstop_mLSP_jet1ptg200",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_SIG = new TH2D("h_mstop_mLSP_SIG","h_mstop_mLSP_SIG",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_neleeq0 = new TH2D("h_mstop_mLSP_neleeq0","h_mstop_mLSP_neleeq0",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_nmueq0 = new TH2D("h_mstop_mLSP_nmueq0","h_mstop_mLSP_nmueq0",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_trackIso = new TH2D("h_mstop_mLSP_trackIso","h_mstop_mLSP_trackIso",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb0Ll = new TH2D("h_mstop_mLSP_g1Mb0Ll","h_mstop_mLSP_g1Mb0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W0Ll = new TH2D("h_mstop_mLSP_g1Mbg1W0Ll","h_mstop_mLSP_g1Mbg1W0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_1Mbg1W0Ll = new TH2D("h_mstop_mLSP_1Mbg1W0Ll","h_mstop_mLSP_1Mbg1W0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g2Mbg1W0Ll = new TH2D("h_mstop_mLSP_g2Mbg1W0Ll","h_mstop_mLSP_g2Mbg1W0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHat4 = new TH2D("h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHat4","h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHat4",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHatg4 = new TH2D("h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHatg4","h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHatg4",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb0Wg1uW0Ll = new TH2D("h_mstop_mLSP_g1Mb0Wg1uW0Ll","h_mstop_mLSP_g1Mb0Wg1uW0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
 
  TH2D* h_mstop_mLSP_0Lb0Ll = new TH2D("h_mstop_mLSP_0Lb0Ll","h_mstop_mLSP_0Lb0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1uW0Ll = new TH2D("h_mstop_mLSP_0Lbg1uW0Ll","h_mstop_mLSP_0Lbg1uW0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1uW0Ll_mdPhi0p3 = new TH2D("h_mstop_mLSP_0Lbg1uW0Ll_mdPhi0p3","h_mstop_mLSP_0Lbg1uW0Ll_mdPhi0p3",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat4 = new TH2D("h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat4","h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat4",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat5 = new TH2D("h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat5","h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat5",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1W0Ll = new TH2D("h_mstop_mLSP_0Lbg1W0Ll","h_mstop_mLSP_0Lbg1W0Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  
  TH2D* h_mstop_mLSP_1Ll = new TH2D("h_mstop_mLSP_1Ll","h_mstop_mLSP_1Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb1Ll = new TH2D("h_mstop_mLSP_g1Mb1Ll","h_mstop_mLSP_g1Mb1Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W1Ll = new TH2D("h_mstop_mLSP_g1Mbg1W1Ll","h_mstop_mLSP_g1Mbg1W1Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W1LlmT100 = new TH2D("h_mstop_mLSP_g1Mbg1W1LlmT100","h_mstop_mLSP_g1Mbg1W1LlmT100",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_1Mbg1W1LlmT100 = new TH2D("h_mstop_mLSP_1Mbg1W1LlmT100","h_mstop_mLSP_1Mbg1W1LlmT100",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g2Mbg1W1LlmT100 = new TH2D("h_mstop_mLSP_g2Mbg1W1LlmT100","h_mstop_mLSP_g2Mbg1W1LlmT100",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1W1LlmT = new TH2D("h_mstop_mLSP_g1Mbg1W1LlmT","h_mstop_mLSP_g1Mbg1W1LlmT",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_0Lb1Ll = new TH2D("h_mstop_mLSP_0Lb1Ll","h_mstop_mLSP_0Lb1Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y1Ll = new TH2D("h_mstop_mLSP_0Lbg1Y1Ll","h_mstop_mLSP_0Lbg1Y1Ll",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y1LlmT100 = new TH2D("h_mstop_mLSP_0Lbg1Y1LlmT100","h_mstop_mLSP_0Lbg1Y1LlmT100",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y1LlmT = new TH2D("h_mstop_mLSP_0Lbg1Y1LlmT","h_mstop_mLSP_0Lbg1Y1LlmT",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_2munoZmass = new TH2D("h_mstop_mLSP_2munoZmass","h_mstop_mLSP_2munoZmass",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2mu = new TH2D("h_mstop_mLSP_2mu","h_mstop_mLSP_2mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2mu0el = new TH2D("h_mstop_mLSP_2mu0el","h_mstop_mLSP_2mu0el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lb2mu0el = new TH2D("h_mstop_mLSP_0Lb2mu0el","h_mstop_mLSP_0Lb2mu0el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb2mu0el = new TH2D("h_mstop_mLSP_g1Mb2mu0el","h_mstop_mLSP_g1Mb2mu0el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y2mu0el = new TH2D("h_mstop_mLSP_0Lbg1Y2mu0el","h_mstop_mLSP_0Lbg1Y2mu0el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1Y2mu0el = new TH2D("h_mstop_mLSP_g1Mbg1Y2mu0el","h_mstop_mLSP_g1Mbg1Y2mu0el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_2elnoZmass = new TH2D("h_mstop_mLSP_2elnoZmass","h_mstop_mLSP_2elnoZmass",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2el = new TH2D("h_mstop_mLSP_2el","h_mstop_mLSP_2el",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2el0mu = new TH2D("h_mstop_mLSP_2el0mu","h_mstop_mLSP_2el0mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lb2el0mu = new TH2D("h_mstop_mLSP_0Lb2el0mu","h_mstop_mLSP_0Lb2el0mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb2el0mu = new TH2D("h_mstop_mLSP_g1Mb2el0mu","h_mstop_mLSP_g1Mb2el0mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y2el0mu = new TH2D("h_mstop_mLSP_0Lbg1Y2el0mu","h_mstop_mLSP_0Lbg1Y2el0mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1Y2el0mu = new TH2D("h_mstop_mLSP_g1Mbg1Y2el0mu","h_mstop_mLSP_g1Mbg1Y2el0mu",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);

  TH2D* h_mstop_mLSP_2lnoZmass = new TH2D("h_mstop_mLSP_2lnoZmass","h_mstop_mLSP_2lnoZmass",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2l = new TH2D("h_mstop_mLSP_2l","h_mstop_mLSP_2l",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_2l0ol = new TH2D("h_mstop_mLSP_2l0ol","h_mstop_mLSP_2l0ol",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lb2l0ol = new TH2D("h_mstop_mLSP_0Lb2l0ol","h_mstop_mLSP_0Lb2l0ol",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mb2l0ol = new TH2D("h_mstop_mLSP_g1Mb2l0ol","h_mstop_mLSP_g1Mb2l0ol",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_0Lbg1Y2l0ol = new TH2D("h_mstop_mLSP_0Lbg1Y2l0ol","h_mstop_mLSP_0Lbg1Y2l0ol",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);
  TH2D* h_mstop_mLSP_g1Mbg1Y2l0ol = new TH2D("h_mstop_mLSP_g1Mbg1Y2l0ol","h_mstop_mLSP_g1Mbg1Y2l0ol",nbins_stop,stop_min,stop_max,nbins_LSP,LSP_min,LSP_max);


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry)
    {
      // Read event into memory
      stream.read(entry);
            
      // Count events and get the total weight contibuted by the event
      double w = 1.;
      if (geneventinfoproduct_weight != 0) {
        w = geneventinfoproduct_weight*w;
      }
      
      // Write every ith event:
      if (entry % 10000 == 0) cout << entry << endl;
      
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file rzrBTanalyzer.h to find out what structs
      // are available. Each struct contains the field "selected", which
      // can be set as needed. Call saveSelectedObjects() before a call to
      // addEvent if you wish to save only the selected objects.
      
      fillObjects();

      // now find out which mass point this event belongs to
      double mg = lheeventproducthelper_mg;
      double mt = lheeventproducthelper_mt1;
     
      // fill out histograms using mt1
      double mt1 = mt;
      if (sample == "T1ttcc" || sample == "T1t1t")
	mt1 = mg;

      double mz1 = lheeventproducthelper_mz1;
      //cout << "mgluino = " << mg << ", mstop = " << mt << ", mLSP = " << mz1 << endl;

      if (mz1 == 0) continue; // mLSP=0 should be rejected

      h_mstop->Fill(mt1);
      h_mLSP->Fill(mz1);
      h_mstop_mLSP->Fill(mt1,mz1,1);
      h_mstop_mLSP_nevents->Fill(mt1,mz1,w);

      //cout << "will fill ofile with weight " << w << endl;
      h_mstop_mLSP_NoCuts->Fill(mt1,mz1,w);

      // Get rid of the noise in data before you start filling ANY histogram
      // by applying the filters:
      if (eventhelper_isRealData==1) {
        if (!(triggerresultshelper1_EcalDeadCellTriggerPrimitiveFilterPath==1)) continue;
        if (!(triggerresultshelper1_hcalLaserEventFilterPath==1)) continue;
        if (!(triggerresultshelper1_trackingFailureFilterPath==1)) continue;
        if (!(triggerresultshelper1_CSCTightHaloFilterPath==1)) continue;
        if (!(triggerresultshelper1_HBHENoiseFilterPath==1)) continue;
        if (!(triggerresultshelper1_primaryVertexFilterPath==1)) continue;
        if (!(triggerresultshelper1_noscrapingFilterPath==1)) continue;
        if (!(triggerresultshelper1_metNoiseCleaningPath==1)) continue;
        if (!(triggerresultshelper1_eeBadScFilterPath==1)) continue;
        if (!(triggerresultshelper1_trkPOGFiltersPath==1)) continue;
        if (!(sint_hcallasereventfilter2012_value==1)) continue;
      }
      // Get rid of the events with wrong kinematics in the MadGraph samples
      if (eventhelper_isRealData!=1) {
        if (triggerresultshelper2_totalKinematicsFilterPath==0) continue;
      }

      h_mstop_mLSP_Cleaning->Fill(mt1,mz1,w);

      // do pileup reweighting
      double num_vertices = pileupsummaryinfo[0].getTrueNumInteractions;
      // get the bin number in the pileup histogram
      int pileup_bin = (int)ceil(num_vertices);
      double w_pileup = 1.;
      if(doPileupReweighting)
	w_pileup = h_pileup->GetBinContent(pileup_bin);      

      w = w*w_pileup;

      h_mstop_mLSP_Pileup->Fill(mt1,mz1,w);

      // ----------------------
      // -- object selection --
      // ----------------------

      // General reference:
      // https://twiki.cern.ch/twiki//bin/view/CMS/Internal/ApprovedObjects

      // vertices - selected:
      std::vector<vertex_s> svertex;
      for (unsigned int i=0; i<vertex.size(); i++) {
	if (vertex[i].isFake) continue;
	if (!(vertex[i].ndof > 4) ) continue;
	if (!(fabs(vertex[i].z) < 24) ) continue;
	if (!(fabs(vertex[i].position_Rho) < 2) ) continue;
	svertex.push_back(vertex[i]);
      }


      // jets - selected:
      // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
      std::vector<cmgpfjet_s> sjet;
      std::vector<TLorentzVector> LVsjet;
      std::vector<cmgpfjet_s> sbjet;
      std::vector<cmgpfjet_s> slbjet;
      double HT = 0;
      // btag probabilities
      double PCSVLsim = 1.0;
      double PCSVLdata = 1.0;
      double PCSVMsim = 1.0;
      double PCSVMdata = 1.0;
      double sigmaSFFl = 0.0;
      double sigmaSFFs = 0.0;
      for (unsigned int i=0; i<cmgpfjet.size(); i++) {
	if (!(cmgpfjet[i].pt > 30) ) continue;
	if (!(fabs(cmgpfjet[i].eta) < 2.4) ) continue;
	//if (!(cmgpfjet[i].neutralHadronEnergyFraction < 0.99) ) continue;
	if (!(cmgpfjet[i].component_5_fraction + cmgpfjet[i].component_6_fraction < 0.99) ) continue;
	//if (!(cmgpfjet[i].neutralEmEnergyFraction < 0.99) ) continue;
	if (!(cmgpfjet[i].component_4_fraction < 0.99) ) continue;
	if (!(cmgpfjet[i].nConstituents > 1) ) continue;
	//if (fabs(cmgpfjet[i].eta) < 2.4) {
	if (!(cmgpfjet[i].component_1_fraction > 0) ) continue;
	if (!(cmgpfjet[i].component_1_number > 0) ) continue;
	if (!(cmgpfjet[i].component_2_fraction < 0.99) ) continue;
	//}
	sjet.push_back(cmgpfjet[i]);
	HT += cmgpfjet[i].pt;

	// btag SF mess
        double SFCSVMFl, dSFCSVMFl, SFCSVMFs, dSFCSVMFs;
        double SFCSVLFl, dSFCSVLFl, SFCSVLFs, dSFCSVLFs;
	double pt = cmgpfjet[i].pt;
	double eta = cmgpfjet[i].eta;
	double partonFlavour = cmgpfjet[i].partonFlavour;
        btagCSVMEEFull(partonFlavour, pt, eta, SFCSVMFl, dSFCSVMFl);
        btagCSVMEEFast(partonFlavour, pt, eta, SFCSVMFs, dSFCSVMFs);

        btagCSVLEEFull(partonFlavour, pt, eta, SFCSVLFl, dSFCSVLFl);
        btagCSVLEEFast(partonFlavour, pt, eta, SFCSVLFs, dSFCSVLFs);

	double eCSVM = 0;
	double eCSVL = 0;
	double SFCSVM, SFCSVL;
	// FastSim:
	if (sample == "T1ttcc" || sample == "T2tt") {
	  if (fabs(partonFlavour) == 5) {
	    eCSVM = geteff1D(h_pt_b_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_b_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) == 4) {
	    eCSVM = geteff1D(h_pt_c_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_c_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) != 4 && fabs(partonFlavour) != 5) {
	    eCSVM = geteff2D(h_pt_eta_l_CSVMeff, pt, fabs(eta));
	    eCSVL = geteff2D(h_pt_eta_l_CSVLeff, pt, fabs(eta));
	  }
	  if (eCSVM == 0 || eCSVL == 0)
	    cout << "Problem is still there" << endl;
	  SFCSVL = (SFCSVLFl + sigmaSFFl*dSFCSVLFl)*(SFCSVLFs + sigmaSFFs*dSFCSVLFs);
	  SFCSVM = (SFCSVMFl + sigmaSFFl*dSFCSVMFl)*(SFCSVMFs + sigmaSFFs*dSFCSVMFs);
	} else { // FullSim
	  if (fabs(partonFlavour) == 5) {
	    eCSVM = geteff1D(h_pt_b_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_b_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) != 5) {
	    eCSVM = geteff1D(h_pt_lc_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_lc_CSVLeff, pt);
	  }
	  SFCSVL = (SFCSVLFl + sigmaSFFl*dSFCSVLFl);
	  SFCSVM = (SFCSVMFl + sigmaSFFs*dSFCSVMFl);
	}

	// CSVM
	if (cmgpfjet[i].combinedSecondaryVertexBJetTags > 0.679) {
	  sbjet.push_back(cmgpfjet[i]);	  
	  PCSVMsim *= eCSVM;
	  PCSVMdata *= (eCSVM * SFCSVM);
	} else {
	  PCSVMsim *= (1 - eCSVM);
	  PCSVMdata *= (1 - eCSVM * SFCSVM);
	}
	// CSVL
	if (cmgpfjet[i].combinedSecondaryVertexBJetTags > 0.244) {
	  slbjet.push_back(cmgpfjet[i]);
	  PCSVLsim *= eCSVL;
	  PCSVLdata *= (eCSVL * SFCSVL);
	} else {
	  PCSVLsim *= (1 - eCSVL);
	  PCSVLdata *= (1 - eCSVL * SFCSVL);
	}

	TLorentzVector jl;
	jl.SetPtEtaPhiE(cmgpfjet[i].pt, cmgpfjet[i].eta,
			cmgpfjet[i].phi, cmgpfjet[i].energy);
	LVsjet.push_back(jl);
      }

      // event weights to be applied according to the selection later
      double wCSVM = PCSVMdata / PCSVMsim;
      double wCSVL = PCSVLdata / PCSVLsim;

      // CA8
      // W selection:
      std::vector<jethelper4_s> sjet2;
      std::vector<jethelper4_s> sW;
      std::vector<jethelper4_s> aW;
      std::vector<jethelper4_s> sY;
      for (unsigned int i=0; i<jethelper4.size(); i++) {
        if (!(jethelper4[i].pt > 30) ) continue;
        if (!(fabs(jethelper4[i].eta) < 3) ) continue;

	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
	sY.push_back(jethelper4[i]);
        //cout << jethelper4[i].pt << endl;
        sjet2.push_back(jethelper4[i]);
        double massdrop = 1;
        double daughmass = -9;
        if (jethelper4[i].daughter_0_mass > jethelper4[i].daughter_1_mass) {
          daughmass = jethelper4[i].daughter_0_mass;
        } else {
          daughmass = jethelper4[i].daughter_1_mass;
        };
        massdrop = daughmass / jethelper4[i].mass;
        double massdrop2 = 1;
        double daughmass2 = -9;
        if (jethelper4[i].daughter_0_pt > jethelper4[i].daughter_1_pt) {
          daughmass2 = jethelper4[i].daughter_0_mass;
        } else {
          daughmass2 = jethelper4[i].daughter_1_mass;
        };
        massdrop2 = daughmass2 / jethelper4[i].mass;
        double dRd1d2 = fdeltaR(jethelper4[i].daughter_0_eta,
                                jethelper4[i].daughter_0_phi,
                                jethelper4[i].daughter_1_eta,
                                jethelper4[i].daughter_1_phi
                                );
        double yasym = (TMath::Min(pow(jethelper4[i].daughter_0_pt, 2),
                                  pow(jethelper4[i].daughter_1_pt, 2))*
                       pow(dRd1d2,2))/
          pow(jethelper4[i].mass,2);

        // Match with the unpruned:
        double prjmatch = 0;
        int jpr = -1;
        double dRmn = 100;
        for (unsigned int j=0; j<jethelper5.size(); j++) {
          double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, jethelper4[i].eta, jethelper4[i].phi);
          if (dR < 0.7 && dR < dRmn) {
            dRmn = dR;
            prjmatch = 1;
            jpr = j;
            break;
          }
        }
	if (!(prjmatch==1)) continue;
	double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
	double tau3 = jethelper5[jpr].tau3;
        //if (!(massdrop < 0.31)) continue;
	if (tau21 >= 0.46 || tau3 >= 0.135) {
          aW.push_back(jethelper4[i]);
        }
	if (!(tau21 < 0.46) ) continue;
	if (!(tau3 < 0.135) ) continue;
        sW.push_back(jethelper4[i]);
      }

      // W selection:
      std::vector<jethelper4_s> sWAK5;
      for (unsigned int i=0; i<jethelper4.size(); i++) {
        if (!(jethelper4[i].pt > 30) ) continue;
        if (!(fabs(jethelper4[i].eta) < 3) ) continue;

	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
        double massdrop = 1;
        double daughmass = -9;
        if (jethelper4[i].daughter_0_mass > jethelper4[i].daughter_1_mass) {
          daughmass = jethelper4[i].daughter_0_mass;
        } else {
          daughmass = jethelper4[i].daughter_1_mass;
        };
        massdrop = daughmass / jethelper4[i].mass;
        if (!(massdrop < 0.31)) continue;
        sWAK5.push_back(jethelper4[i]);
      }



      // Muons - veto:
      // From https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon
      std::vector<cmgmuon_s> vmuon;
      for (unsigned int i=0; i<cmgmuon.size(); i++) {
	if (!(cmgmuon[i].pt > 5) ) continue;
	if (!(fabs(cmgmuon[i].eta) < 2.4) ) continue;
	vmuon.push_back(cmgmuon[i]);
      }
      // Muons - tight
      std::vector<cmgmuon1_s> smuon;
      std::vector<TVector3> V3mu;
      std::vector<TLorentzVector> LVmu;
      for (unsigned int i=0; i<cmgmuon1.size(); i++) {
	if (!(cmgmuon1[i].pt > 10) ) continue; // lowered from 25
	if (!(fabs(cmgmuon1[i].eta) < 2.4) ) continue;
	smuon.push_back(cmgmuon1[i]);
        TVector3 lmu;
        lmu.SetPtEtaPhi(cmgmuon1[i].pt, 0, cmgmuon1[i].phi);
        V3mu.push_back(lmu);
	TLorentzVector lvmu;
	lvmu.SetPtEtaPhiE(cmgmuon1[i].pt, cmgmuon1[i].eta,
			  cmgmuon1[i].phi, cmgmuon1[i].energy);
        LVmu.push_back(lvmu);
      }

      // Electrons - veto:
      // Use the "veto" criteria
      // From 
      std::vector<cmgelectron_s> velectron;
      for (unsigned int i=0; i<cmgelectron.size(); i++) {
	if (!(cmgelectron[i].pt > 5) ) continue;
	if (!(fabs(cmgelectron[i].eta) < 2.5) ) continue;
	if (!(fabs(cmgelectron[i].eta) < 1.442 && fabs(cmgelectron[i].eta) < 1.556) ) continue;
	velectron.push_back(cmgelectron[i]);
      }
      // Electrons - tight
      std::vector<cmgelectron1_s> selectron;
      std::vector<TVector3> V3el;
      std::vector<TLorentzVector> LVel;
      for (unsigned int i=0; i<cmgelectron1.size(); i++) {
        if (!(cmgelectron1[i].pt > 10) ) continue;
        if (!(fabs(cmgelectron1[i].eta) < 2.5) ) continue;
        // veto 1.442 < |eta| < 1.556? --> is already done in the collection ??
        if (!(fabs(cmgelectron1[i].eta) < 1.442 && fabs(cmgelectron1[i].eta) < 1.556) ) continue;
        selectron.push_back(cmgelectron1[i]);
        TVector3 lel;
        lel.SetPtEtaPhi(cmgelectron1[i].pt, 0, cmgelectron1[i].phi);
        V3el.push_back(lel);
	TLorentzVector lvel;
	lvel.SetPtEtaPhiE(cmgelectron1[i].pt, cmgelectron1[i].eta,
			  cmgelectron1[i].phi, cmgelectron1[i].energy);
        LVel.push_back(lvel);

      }

      // Taus - veto
      // From supposed to behttps://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation#TauID
      // but since that twiki is so horrible, I just followed Will's multijet note.  Screw'em!
      std::vector<cmgtau_s> vtau;
      for (unsigned int i=0; i<cmgtau.size(); i++) {
	if (!(cmgtau[i].pt > 10) ) continue;
	if (!(fabs(cmgtau[i].eta) < 2.3) ) continue;
	//if (!(tau2[i].byMediumCombinedIsolationDeltaBetaCorr == 1) ) continue;
	vtau.push_back(cmgtau[i]);
      }
      // Taus - tight
      std::vector<cmgtau1_s> stau;
      for (unsigned int i=0; i<cmgtau1.size(); i++) {
        if (!(cmgtau1[i].pt > 10) ) continue;
        if (!(fabs(cmgtau1[i].eta) < 2.3) ) continue;
        //if (!(tau2[i].byMediumCombinedIsolationDeltaBetaCorr == 1) ) continue; 
        stau.push_back(cmgtau1[i]);
      }


      // look at number of Ws and bs before any selection (besides cleaning and trigger)
      double nWs = sW.size();
      double nYs = sY.size();
      double nbs = sbjet.size();

      // -------------------------------
      // -- Calculate Razor variables --
      // -------------------------------

      // Calculate MR and R2 ignoring muons
      TVector3 V3met;
      V3met.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      TLorentzVector met;
      met.SetPtEtaPhiE(cmgbasemet2[0].pt, 0, cmgbasemet2[0].phi, cmgbasemet2[0].energy);
      std::vector<TLorentzVector> LVhemis = CombineJets(LVsjet);

      double MR = -9999;
      double MTR = -9999;
      double R2 = -9999;
      if (LVhemis.size() == 2) {
	MR = CalcMR(LVhemis[0], LVhemis[1]);
	if (MR != MR) continue;
	MTR = CalcMTR(LVhemis[0], LVhemis[1], V3met);
	R2 = pow((MTR / MR),2);
      }
      
      // Calculate MR and R2 adding mus to MET
      TVector3 V3metmu;
      V3metmu.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      for (unsigned int i=0; i<V3mu.size(); i++) {
        V3metmu += V3mu[i];
      }

      double MTRmetmu = -9999;
      double R2metmu = -9999;
      if (LVhemis.size() == 2) {
        MTRmetmu = CalcMTR(LVhemis[0], LVhemis[1], V3metmu);
        R2metmu = pow((MTRmetmu / MR),2);
      }

      // Calculate MR and R2 adding electrons to MET
      TVector3 V3metel;
      V3metel.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      for (unsigned int i=0; i<V3el.size(); i++) {
        V3metel += V3el[i];
      }

      double MTRmetel = -9999;
      double R2metel = -9999;
      if (LVhemis.size() == 2) {
        MTRmetel = CalcMTR(LVhemis[0], LVhemis[1], V3metel);
        R2metel = pow((MTRmetel / MR),2);
      }


      // --------------------------------------------------------
      // -- Everything computed from generator level particles --
      // --------------------------------------------------------

 
      // *****************************************************
      // ***  ISR Reweighting recipe for Madgraph samples  ***
      // *****************************************************
      // per event weights to apply
      double w_ISR_nominal = 1.;
      double w_ISR_up = 1.; // always stays 1, i.e. no reweighting
      double w_ISR_down = 1.;
      if (doISRreweighting)
	{
	  // recipe can be found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMST2ccMadgraph8TeV
	  // find system recoiling against ISR: 
	  TLorentzVector recoilsystem(0,0,0,0); 
	  int ID_to_find = -1;
	  if (sample == "T2tt")
	    ID_to_find = 1000006;
	  if (sample == "T1ttcc" || sample == "T1ttcc_old" || sample == "T1t1t")
	    ID_to_find = 1000021;
	  if (sample == "TTJets")
	    ID_to_find = 6;
	  if (sample == "WJets")
	    ID_to_find = 24;
	  if (sample == "ZJets")
	    ID_to_find = 23;
	  for (unsigned int i=0; i<genparticlehelper.size(); i++) {
	    if (genparticlehelper[i].status != 3) continue;
	    if (fabs(genparticlehelper[i].pdgId) == ID_to_find){
	      TLorentzVector TLV_temp; 
	      TLV_temp.SetPtEtaPhiM(genparticlehelper[i].pt,genparticlehelper[i].eta
				    ,genparticlehelper[i].phi,genparticlehelper[i].mass);
	      recoilsystem += TLV_temp;
	    }
	  }
	  // get the pt of the recoil system, apply weights accordingly
	  double pt_ISR = recoilsystem.Pt();
	  if (pt_ISR <= 120){
	    w_ISR_nominal = 1.;
	    w_ISR_down = 1.;
	  } else if (pt_ISR <= 150){
	    w_ISR_nominal = 0.95;
	    w_ISR_down = 0.9;
	  } else if (pt_ISR <= 250){
	    w_ISR_nominal = 0.9;
	    w_ISR_down = 0.8;
	  } else {
	    w_ISR_nominal = 0.8;
	    w_ISR_down = 0.6;
	  }
	}

      if (doISRreweighting) {
	w = w*w_ISR_nominal;
      }

      h_mstop_mLSP_ISR->Fill(mt1,mz1,w);  


      // **************************************************************
      // ***  Top Pt Reweighting recipe for Madgraph TTbar samples  ***
      // **************************************************************
      // per event weights to apply
      double w_TopPt_nominal = 1.;
      double w_TopPt_up = 1.; 
      double w_TopPt_down = 1.; // always stays 1, i.e. no reweighting     
      if(doTopPtreweighting){
	for (unsigned int i=0; i<genparticlehelper.size(); i++) {
	  if (genparticlehelper[i].status != 3) continue;
	  if (fabs(genparticlehelper[i].pdgId) == 6){
	    double wtemp = GetTopPtScaleFactor(genparticlehelper[i].pt);
	    //cout << "wtemp: " << wtemp << endl;
	    w_TopPt_nominal = w_TopPt_nominal * wtemp;
	  }
	}
	w_TopPt_up = w_TopPt_nominal;
	w_TopPt_nominal = sqrt(w_TopPt_nominal);
      }



      // ---------------------
      // -- event selection --
      // ---------------------

      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;
      h_mstop_mLSP_HCAL_noise->Fill(mt1,mz1,w);

      // at least one good primary vertex
      if (!(svertex.size() > 0)) continue;
      h_mstop_mLSP_vertexg0->Fill(mt1,mz1,w);

      // at least three jets
      if (!(sjet.size() >= 3)) continue;
      h_mstop_mLSP_njetge3->Fill(mt1,mz1,w);
      
      // Calculate the HLT weight and include it in the total weight:
      double whlt = 1;
      if (eventhelper_isRealData==0) {
	for (int i=1; i<h_hlteff->GetNbinsX()+1; i++) {
	  double xmin = h_hlteff->GetXaxis()->GetBinLowEdge(i);
	  double xmax = h_hlteff->GetXaxis()->GetBinUpEdge(i);
	  if (!(HT >= xmin && HT < xmax)) continue;
	  for (int j=1; j<h_hlteff->GetNbinsY()+1; j++) {
	    double ymin = h_hlteff->GetYaxis()->GetBinLowEdge(j);
	    double ymax = h_hlteff->GetYaxis()->GetBinUpEdge(j);
	    if (sjet[0].pt >= ymin && sjet[0].pt < ymax) {
	      whlt = h_hlteff->GetBinContent(i, j);
	      //cout << xmin << " " << MR << " " << xmax << " " << ymin << " " << R2 << " " << ymax << " " << whlt << endl;
	      break;
	    }
	  }
	}
      }
      w = w*whlt;
      h_mstop_mLSP_HLT->Fill(mt1,mz1,w);

      // pt of first jet greater than 200 GeV
      if (!(sjet[0].pt > 200)) continue;
      h_mstop_mLSP_jet1ptg200->Fill(mt1,mz1,w);

      // count number of leptons
      int nlooseelectrons = velectron.size();
      int nloosemuons = vmuon.size(); 
      int nlooseleptons = nlooseelectrons + nloosemuons;
      int ntightmuons = smuon.size();
      int ntightelectrons = selectron.size();
      int ntightleptons = ntightelectrons + ntightmuons;
      
      // count number of b jets
      int nloosebs = slbjet.size();
      int nmediumbs = sbjet.size();
      
      // Only select events in MR-R2 SIG region 
      //if (!(MR > 800 && R2 > 0.08)) continue;
      if (MR > 800 && R2 > 0.08){
	h_mstop_mLSP_SIG->Fill(mt1,mz1,w);
	
	// Compute the minDeltaPhi variable, taking the first three jets into account
	// Compute the minDeltaPhiHat variable, taking the first three jets into account
	double minDeltaPhi = 99.;
	double minDeltaPhiHat = 99;
	for (int jet=0; jet<3; ++jet){
	  double mdphi = fdeltaPhi(sjet[jet].phi,V3met.Phi());
	  if (mdphi < minDeltaPhi)
	    minDeltaPhi = mdphi;
	  // compute resolution on T_i
	  double sigma_T_i_2 = 0;
	  for (int j=0; j<(signed)sjet.size(); ++j){
	    if(j==jet) continue;
	    double alpha = fdeltaPhi(sjet[jet].phi,sjet[j].phi); 
	    sigma_T_i_2 += pow( 0.1*sjet[j].pt*sin(alpha) ,2); // could also have to be sin(pi - alpha), but sin is the same for both anyways...
	  }
	  // compute resolution on dphi_i
	  double sigma_dphi_i = atan(sqrt(sigma_T_i_2)/V3met.Pt());
	  double mdphihat = mdphi/sigma_dphi_i;
	  if (mdphihat < minDeltaPhiHat)
	    minDeltaPhiHat = mdphihat;
	}

	// ----------------------------------------------------------------------------------------------------
	// 0 Lepton trajectory
	// ----------------------------------------------------------------------------------------------------
	if (nlooseelectrons == 0){
	  h_mstop_mLSP_neleeq0->Fill(mt1,mz1,w);
	  
	  if (nloosemuons == 0) {
	    h_mstop_mLSP_nmueq0->Fill(mt1,mz1,w);
	    
	    if (eventhelperextra_trackIso == 0){
	      h_mstop_mLSP_trackIso->Fill(mt1,mz1,w);
	      
	      if (nmediumbs > 0){
		w = w*wCSVM;
		h_mstop_mLSP_g1Mb0Ll->Fill(mt1,mz1,w);

		// g1Mb g1W 0Ll -- SIGNAL region
		if( sW.size() > 0){
		  h_mstop_mLSP_g1Mbg1W0Ll->Fill(mt1,mz1,w);
		  
		  if (nmediumbs == 1){
		    h_mstop_mLSP_1Mbg1W0Ll->Fill(mt1,mz1,w); 
		  } else {
		    h_mstop_mLSP_g2Mbg1W0Ll->Fill(mt1,mz1,w); 
		  }

		  if (minDeltaPhiHat > 4){
		    h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHatg4->Fill(mt1,mz1,w);
		  } else {
		    h_mstop_mLSP_g1Mbg1W0Ll_mdPhiHat4->Fill(mt1,mz1,w);
		  }
		  
		} // end of sW.size() > 0
		else if (aW.size() > 0){
		  h_mstop_mLSP_g1Mb0Wg1uW0Ll->Fill(mt1,mz1,w);
		}
	      } // end of nmediumbs > 0
	      
	      if (nloosebs == 0){
		w = w*wCSVL;
		h_mstop_mLSP_0Lb0Ll->Fill(mt1,mz1,w);
		
		// 0Lbg1uW0Ll -- QCD control region
		if( aW.size() > 0){
		  h_mstop_mLSP_0Lbg1uW0Ll->Fill(mt1,mz1,w);
		  
		  // cut on mindDeltaPhi
		  if (minDeltaPhi < 0.3){
		    h_mstop_mLSP_0Lbg1uW0Ll_mdPhi0p3->Fill(mt1,mz1,w);
		  } // end of minDeltaPhi < 0.3

		  if (minDeltaPhiHat < 4){
		    h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat4->Fill(mt1,mz1,w);
		  } // end of minDeltaPhiHat < 4

		  if (minDeltaPhiHat < 5){
		    h_mstop_mLSP_0Lbg1uW0Ll_mdPhiHat5->Fill(mt1,mz1,w);
		  } // end of minDeltaPhiHat < 5
		} // end of aW.size() > 0
		
		// 0Lbg1W0Ll
		if( sW.size() > 0){
		  h_mstop_mLSP_0Lbg1W0Ll->Fill(mt1,mz1,w);
		} // end of sW.size() > 0
		
	      } // end of nloosebs == 0
	    } // end veto iso track
	  } // end veto muon
	}  // end veto electron
	
	
	//---------------------------------------------------------------------------------
	// 1 Loose lepton trajectory
	//---------------------------------------------------------------------------------
	if(nlooseleptons == 1){
	  // Calculate mT 
	  TLorentzVector lepton;
	  if (nlooseelectrons == 1)
	    lepton.SetPtEtaPhiE(velectron[0].pt, velectron[0].eta, velectron[0].phi, velectron[0].energy);
	  else if (nloosemuons == 1)
	    lepton.SetPtEtaPhiE(vmuon[0].pt, vmuon[0].eta, vmuon[0].phi, vmuon[0].energy);
	  double mT = CalcMT(lepton,met);
	  
	  h_mstop_mLSP_1Ll->Fill(mt1,mz1,w);
	  
	  if (nmediumbs > 0){
	    w = w*wCSVM;
	    h_mstop_mLSP_g1Mb1Ll->Fill(mt1,mz1,w);
	    
	    if( sW.size() > 0 ){
	      h_mstop_mLSP_g1Mbg1W1Ll->Fill(mt1,mz1,w);

	      // TTJets Control region
	      if (mT < 100){
		h_mstop_mLSP_g1Mbg1W1LlmT100->Fill(mt1,mz1,w);

		if (nmediumbs == 1){
		  h_mstop_mLSP_1Mbg1W1LlmT100->Fill(mt1,mz1,w);
		} else {
		  h_mstop_mLSP_g2Mbg1W1LlmT100->Fill(mt1,mz1,w);
		}
		
		// mT window
		if (mT > 30){
		  h_mstop_mLSP_g1Mbg1W1LlmT->Fill(mt1,mz1,w);
		} // end mT > 30
	      } // end mT < 100
	    } // end sW.size()
	  } // end nmediumbs > 0


	  if (nloosebs == 0){
	    w = w*wCSVL;
	    h_mstop_mLSP_0Lb1Ll->Fill(mt1,mz1,w);
	    
	    if( sY.size() > 0 ){
	      h_mstop_mLSP_0Lbg1Y1Ll->Fill(mt1,mz1,w);
	      
	      // WJets Control Region
	      if (mT < 100){ 
		h_mstop_mLSP_0Lbg1Y1LlmT100->Fill(mt1,mz1,w);
		
		// mT window
		if (mT > 30){
		  h_mstop_mLSP_0Lbg1Y1LlmT->Fill(mt1,mz1,w);
		} // end mT > 30
	      } // end mT < 100
	    } // end sY.size()
	  } // end nloosebs > 0

	} // end nlooseleptons == 1
      } // end of MR>800 R2>0.08
      

      // -----------------------------------------------------------------------------
      // Dilepton trajectory
      // -----------------------------------------------------------------------------

      // Start the 2mu stuff here.
      // Need to use R2metmu to define the signal region
      if (MR > 800 && R2metmu > 0.08){
	TLorentzVector LVZcand;
	if (ntightmuons == 2 && nloosemuons == 2) {
	  // Make sure that the muons are opposite-signed:
	  if (!(smuon[0].charge == -smuon[1].charge)) continue;
	  // See if the 2mu makes a Z:
	  for (unsigned int m=0; m<LVmu.size(); m++) {
	    LVZcand += LVmu[m];
	  }
	  double Zmass = LVZcand.M();

	  h_mstop_mLSP_2munoZmass->Fill(mt1,mz1,w);
	  h_mstop_mLSP_2lnoZmass->Fill(mt1,mz1,w);
	  
	  if (!(Zmass >= 60 && Zmass <= 120)) continue;
	  h_mstop_mLSP_2mu->Fill(mt1,mz1,w);
	  h_mstop_mLSP_2l->Fill(mt1,mz1,w);
	  
	  if (nlooseelectrons == 0){
	    h_mstop_mLSP_2mu0el->Fill(mt1,mz1,w);
	    h_mstop_mLSP_2l0ol->Fill(mt1,mz1,w);
	    
	    if (nloosebs == 0){
	      w = w*wCSVL;
	      h_mstop_mLSP_0Lb2mu0el->Fill(mt1,mz1,w);
	      h_mstop_mLSP_0Lb2l0ol->Fill(mt1,mz1,w);
	      
	      if (nYs >= 1){ // Z no b, mu CR 
		h_mstop_mLSP_0Lbg1Y2mu0el->Fill(mt1,mz1,w);
		h_mstop_mLSP_0Lbg1Y2l0ol->Fill(mt1,mz1,w);
	      }// nYs >= 1
	    }// end nloosebs == 0

	    if (nmediumbs >= 1){
	      w = w*wCSVM;
	      h_mstop_mLSP_g1Mb2mu0el->Fill(mt1,mz1,w);
	      h_mstop_mLSP_g1Mb2l0ol->Fill(mt1,mz1,w);
	      
	      if (nYs >= 1){ // Z with b, mu CR
		h_mstop_mLSP_g1Mbg1Y2mu0el->Fill(mt1,mz1,w);
		h_mstop_mLSP_g1Mbg1Y2l0ol->Fill(mt1,mz1,w);
	      } // end nYs >= 1
	    } // end nmediumbs >= 1
	  } // end nlooseelectrons == 0
	} // end ntightmuons == 2
      } // end of MR>800 R2metmu>0.08

      // Start the 2el stuff here.
      // Need to use R2metel to define the signal region
      if (MR > 800 && R2metel > 0.08){
	TLorentzVector LVZcand;
	if (ntightelectrons == 2 && nlooseelectrons == 2) {
	  // Make sure that the electrons are opposite-signed:
	  if (!(selectron[0].charge == -selectron[1].charge)) continue;
	  // See if the 2el makes a Z:
	  for (unsigned int e=0; e<LVel.size(); e++) {
	    LVZcand += LVel[e];
	  }
	  double Zmass = LVZcand.M();

	  h_mstop_mLSP_2elnoZmass->Fill(mt1,mz1,w);
	  h_mstop_mLSP_2lnoZmass->Fill(mt1,mz1,w);

     	  if (!(Zmass >= 60 && Zmass <= 120)) continue;
	  h_mstop_mLSP_2el->Fill(mt1,mz1,w);
	  h_mstop_mLSP_2l->Fill(mt1,mz1,w);
	  
	  if (nloosemuons == 0){
	    h_mstop_mLSP_2el0mu->Fill(mt1,mz1,w);
	    h_mstop_mLSP_2l0ol->Fill(mt1,mz1,w);
	    
	    if (nloosebs == 0){
	      w = w*wCSVL;
	      h_mstop_mLSP_0Lb2el0mu->Fill(mt1,mz1,w);
	      h_mstop_mLSP_0Lb2l0ol->Fill(mt1,mz1,w);

	      if (nYs >= 1){
		h_mstop_mLSP_0Lbg1Y2el0mu->Fill(mt1,mz1,w);
		h_mstop_mLSP_0Lbg1Y2l0ol->Fill(mt1,mz1,w);
	      }// nYs >= 1
	    }// end nloosebs == 0
	    
	    if (nmediumbs >= 1){
	      w = w*wCSVM;
	      h_mstop_mLSP_g1Mb2el0mu->Fill(mt1,mz1,w);
	      h_mstop_mLSP_g1Mb2l0ol->Fill(mt1,mz1,w);
	    
	      if (nYs >= 1){
		h_mstop_mLSP_g1Mbg1Y2el0mu->Fill(mt1,mz1,w);
		h_mstop_mLSP_g1Mbg1Y2l0ol->Fill(mt1,mz1,w);
	      } // end nYs >= 1
	    } // end nmediumbs >= 1
	  } // end nlooseelectrons == 0
	} // end ntightelectrons == 2
      } // end of MR>800 R2metel>0.08


    } // end event loop
  
  fhlt->Close();
  fpileup->Close();
  fbeff->Close();
  stream.close();
  ofile.close();
  return 0;
}
