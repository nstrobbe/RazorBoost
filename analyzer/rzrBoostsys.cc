
//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmdsys.h"
#include "utils.h"
#include "systutils.h"
#include <math.h>
//#include "JetCorrectionUncertainty.h"
//#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"

#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/pdg.h"
#else
#include "pdg.h"
#endif

using namespace std;
//using namespace LHAPDF;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{

  // Get the trigger histogram:
  //TFile* fhlt = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  //if (!fhlt){
  TFile*  fhlt = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  //}
  if (!fhlt){
    cout << "Could not find trigger efficiency root file... Where did you put it??" << endl;
    return 1;
  }
  TH2D* h_hlteff = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph");
  TH2D* h_hlteff_up = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph_up");
  TH2D* h_hlteff_low = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph_low");

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

  cout << "Output file: " << cmdline.outputfilename.c_str() << endl;
  cout << "Systematics file: " << cmdline.systfilename.c_str() << endl;

  // Open the systematics file:
  ifstream systFile(cmdline.systfilename.c_str());
  if ( !systFile.good() ) error("unable to open systematics file: " + cmdline.systfilename);
  std::vector<double> vsyst;
  double syst;
  while (systFile >> syst) {
    cout << "Systematic: " << syst << endl;
    vsyst.push_back(syst);
  }

  // Assign the systematic sigmas:
  int nsyst = vsyst.size();
  double sigmaJEC = vsyst[0];
  double sigmaJECCA8 = vsyst[1];
  double sigmaHLT = vsyst[2];
  double sigmabtagFl = vsyst[3];
  double sigmabtagFs = vsyst[4];
  double sigmaW = vsyst[5];
  double sigmaeleVeto = vsyst[6];
  double sigmaPU = vsyst[7];
  double sigmaISR = vsyst[8];
  double sigmaTopPt = vsyst[9];

  string sample = "";
  if ( argc > 7 )
    sample = string(argv[7]);
  string ISR = "";
  if ( argc > 8 )
    ISR = string(argv[8]);
  string TopPt = "";
  if ( argc > 9 )
    TopPt = string(argv[9]);
  string Pileup = "";
  if ( argc > 10 )
    Pileup = string(argv[10]);
  string Runs = "";
  if ( argc > 11 )
    Runs = string(argv[11]);

  bool doISRreweighting = false;
  if (ISR == "ISR_True" 
      && (sample == "T2tt" || sample == "T1ttcc_old" || sample == "T1t1t"
	  || sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
	  || sample == "TTJets" || sample == "WJets" || sample == "ZJets" )
      ){
    doISRreweighting = true;
    cout << "Will do ISR reweighting" << endl;
  }
  bool doTopPtreweighting = false;
  if (sample == "TTJets" && TopPt == "TopPt_True"){
    doTopPtreweighting = true;
    doISRreweighting = false; // Never do both ISR and TopPt reweighting
    cout << "Will do top pt reweighting" << endl;
  }
  bool doPileupReweighting = false;
  if (sample != "Data" && Pileup == "Pileup_True"){
    doPileupReweighting = true;
    cout << "Will do pileup reweighting" << endl;
  }

  // ---------------------------------------
  // --- Get the correct pileup histogram --
  // ---------------------------------------

  TString pileupname = "pileup_weights.root";
  TString pileupname_up = "pileup_weights_up.root";
  TString pileupname_down = "pileup_weights_down.root";
  if (Runs == "AB"){
    pileupname = "pileup_weights_AB.root";
    cout << "Using pileup profile for runs A+B only" << endl;
  }
  if (sample == "T1ttcc_old" || sample == "T2tt"){
    pileupname = "pileup_weights_sig52X.root";
    pileupname_up = "pileup_weights_sig52X_up.root";
    pileupname_down = "pileup_weights_sig52X_down.root";
    if (Runs == "AB")
      pileupname = "pileup_weights_AB_sig52X.root";
  }
  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T1t1t"){
    pileupname = "pileup_weights_sig53X.root";
    pileupname_up = "pileup_weights_sig53X_up.root";
    pileupname_down = "pileup_weights_sig53X_down.root";
    if (Runs == "AB")
      pileupname = "pileup_weights_AB_sig53X.root";
  }
  //TFile* fpileup = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname);
  //TFile* fpileup_up = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname_up);
  //TFile* fpileup_down = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname_down);
  //if (!fpileup){
  TFile*  fpileup = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/pileup/"+pileupname);
  //}
  //if (!fpileup_up){
  TFile* fpileup_up = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/pileup/"+pileupname_up);
  //}
  //if (!fpileup_down){
  TFile* fpileup_down = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/pileup/"+pileupname_down);
  //}
  if (!fpileup || !fpileup_up || !fpileup_down){
    cout << "Could not find pileup weights root file... Where did you put it??" << endl;
    return 1;
  }
  TH1D* h_pileup = (TH1D*)fpileup->Get("pileup_weight");
  TH1D* h_pileup_up = (TH1D*)fpileup_up->Get("pileup_weight");
  TH1D* h_pileup_down = (TH1D*)fpileup_down->Get("pileup_weight");

  // ---------------------------------------
  // -- Get the btag eff histograms:
  // ---------------------------------------

  TString beff_name = "";
  if (sample == "TTJets" or sample == "Top" or sample == "TTX") {
    beff_name = "btageff_TTJets.root";
  } else if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
	     || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t") {
    beff_name = "btageff_T1ttcc.root";
  } else {
    beff_name = "btageff_QCD.root";
  }
  TFile* fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/"+beff_name);
  if (!fbeff){
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/"+beff_name);
  }
  if (!fbeff){
    cout << "Could not find btag root file... Where did you put it??" << endl;
    return 1;
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


  // JEC uncertainty class:
  //JetCorrectionUncertainty* jecUnc = new JetCorrectionUncertainty("MCcorrs/JEC/Summer13_V5_MC_Uncertainty_AK5PF.txt");
  //JetCorrectionUncertainty* jecUncCA8 = new JetCorrectionUncertainty("MCcorrs/JEC/Summer13_V5_MC_Uncertainty_AK7PFchs.txt");


  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights
  // -- The original MC weight will be divided by this quantity
  // --------------------------------------------------------------
  double weightnorm = 1.;
  if (xsect != -1 && totweight != -1 && lumi != -1) {
    weightnorm = (xsect*lumi)/totweight;
  }

  // for SMS's we take a different approach, and want histograms to be normalized to the efficiency
  // So normalization only uses "totweight" and not xsect or lumi
  TH2D* h_smscounts;
  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25"  || sample == "T1ttcc_DM80" 
      || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
    // open file with counts
    TFile* f_smscounts = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/smsinput/signal_counts.root");
    if (!f_smscounts)
      f_smscounts = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/smsinput/signal_counts.root");
    // get the proper histogram
    TString hname_smscounts = "";
    if (doISRreweighting){
      hname_smscounts = sample+"_ISR";
    } else {
      hname_smscounts = sample+"_noISR";
    }
    h_smscounts = (TH2D*)f_smscounts->Get(hname_smscounts);
    weightnorm = 1.;
  }

  cout << "lumi: " << lumi << endl;
  cout << "xsect: " << xsect << endl;
  cout << "totweight: " << totweight << endl;
  cout << "weightnorm: " << weightnorm << endl;

  outputFile ofile(cmdline.outputfilename);

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------
  TH1::SetDefaultSumw2();

  string binning = "variable";
  //string binning = "";

  // for fixed bin width
  int nbins_MR = 20;
  double MRmn = 0;
  double MRmx = 4000;
  Double_t* bins_MR = getFixedBinEdges(nbins_MR, MRmn, MRmx);

  int nbins_R2 = 25;
  double R2mn = 0;
  double R2mx = 1;
  Double_t* bins_R2 = getFixedBinEdges(nbins_R2, R2mn, R2mx);

  // for variable bin width
  if (binning == "variable"){
    nbins_MR = 7;
    nbins_R2 = 7;
    Double_t bins_MR_tmp[] = {0.,600.,800.,1000.,1200.,1600.,2000.,4000.};
    bins_MR = 0;
    bins_MR = getVariableBinEdges(nbins_MR+1,bins_MR_tmp);
    Double_t bins_R2_tmp[] = {0.,0.04,0.08,0.12,0.16,0.24,0.5,1.};
    bins_R2 = 0;
    bins_R2 = getVariableBinEdges(nbins_R2+1,bins_R2_tmp);
  }

  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);

  // g1Mb g1W 0Ll ; Signal box: >= 1 Mb; >= 1 W; 0 Ll
  TH1D * h_MR_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhiHatg4", "h_MR_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH2D("h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhiHat < 4 (RA2b value)
  TH1D * h_MR_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhiHat4", "h_MR_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhiHat4", "h_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhiHat4", "h_MR_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH2D("h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4", "h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // g1Mbg1W1LlmT ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll; mT<100
  TH1D * h_MR_g1Mbg1W1LlmT100 = new TH1D("h_MR_g1Mbg1W1LlmT100", "h_MR_g1Mbg1W1LlmT100", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100 = new TH1D("h_R2_g1Mbg1W1LlmT100", "h_R2_g1Mbg1W1LlmT100", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100", "h_MR_R2_g1Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W1LlmT100 = new TH2D("h_uw_MR_R2_g1Mbg1W1LlmT100", "h_uw_MR_R2_g1Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // 0Lbg1Y1LlmT ; Wj control region: == 0 Lb; >= 1 Y; 1 Ll; 30<mT<100
  TH1D * h_MR_0Lbg1Y1LlmT = new TH1D("h_MR_0Lbg1Y1LlmT", "h_MR_0Lbg1Y1LlmT", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT = new TH1D("h_R2_0Lbg1Y1LlmT", "h_R2_0Lbg1Y1LlmT", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT = new TH2D("h_MR_R2_0Lbg1Y1LlmT", "h_MR_R2_0Lbg1Y1LlmT", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1Y1LlmT = new TH2D("h_uw_MR_R2_0Lbg1Y1LlmT", "h_uw_MR_R2_0Lbg1Y1LlmT", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // Define the order of bins in the counts histogram:
  
  ofile.count("NoCuts", 0.0);
  ofile.count("Cleaning", 0.0);
  ofile.count("Pileup", 0.0);
  ofile.count("HCAL_noise", 0.0);
  ofile.count("vertexg0", 0.0);
  ofile.count("njetge3", 0.0);
  ofile.count("HLT", 0.0);
  ofile.count("jet1ptg200", 0.0);
  ofile.count("SIG", 0.0); // MR>800 R2>0.08

  ofile.count("neleeq0", 0.0);
  ofile.count("nmueq0", 0.0);
  ofile.count("trackIso", 0.0);
  ofile.count("g1Mb0Ll", 0.0);
  ofile.count("g1Mbg1W0Ll", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  ofile.count("0Lb0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhiHat4", 0.0);

  ofile.count("1Ll", 0.0);
  ofile.count("g1Mb1Ll", 0.0);
  ofile.count("g1Mbg1W1Ll", 0.0);
  ofile.count("g1Mbg1W1LlmT100", 0.0);
  ofile.count("0Lb1Ll", 0.0);
  ofile.count("0Lbg1Y1Ll", 0.0);
  ofile.count("0Lbg1Y1LlmT", 0.0);


  
  TH1D* TTallhad = new TH1D("counts_TTallhad","",1,0,1);
  TTallhad->SetBit(TH1::kCanRebin);
  TH1D* TTsemilep = new TH1D("counts_TTsemilep","",1,0,1);
  TTsemilep->SetBit(TH1::kCanRebin);
  TH1D* TTdilep = new TH1D("counts_TTdilep","",1,0,1);
  TTdilep->SetBit(TH1::kCanRebin);

  TTallhad->Fill("NoCuts", 0.0);
  TTallhad->Fill("Cleaning", 0.0);
  TTallhad->Fill("Pileup", 0.0);
  TTallhad->Fill("HCAL_noise", 0.0);
  TTallhad->Fill("vertexg0", 0.0);
  TTallhad->Fill("njetge3", 0.0);
  TTallhad->Fill("HLT", 0.0);
  TTallhad->Fill("jet1ptg200", 0.0);
  TTallhad->Fill("SIG", 0.0); // MR>800 R2>0.08
  TTallhad->Fill("neleeq0", 0.0);
  TTallhad->Fill("nmueq0", 0.0);
  TTallhad->Fill("trackIso", 0.0);
  TTallhad->Fill("g1Mb0Ll", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTallhad->Fill("0Lb0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTallhad->Fill("1Ll", 0.0);
  TTallhad->Fill("g1Mb1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100", 0.0);
  TTallhad->Fill("0Lb1Ll", 0.0);
  TTallhad->Fill("0Lbg1Y1Ll", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT", 0.0);

  TTsemilep->Fill("NoCuts", 0.0);
  TTsemilep->Fill("Cleaning", 0.0);
  TTsemilep->Fill("Pileup", 0.0);
  TTsemilep->Fill("HCAL_noise", 0.0);
  TTsemilep->Fill("vertexg0", 0.0);
  TTsemilep->Fill("njetge3", 0.0);
  TTsemilep->Fill("HLT", 0.0);
  TTsemilep->Fill("jet1ptg200", 0.0);
  TTsemilep->Fill("SIG", 0.0); // MR>800 R2>0.08
  TTsemilep->Fill("neleeq0", 0.0);
  TTsemilep->Fill("nmueq0", 0.0);
  TTsemilep->Fill("trackIso", 0.0);
  TTsemilep->Fill("g1Mb0Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTsemilep->Fill("0Lb0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTsemilep->Fill("1Ll", 0.0);
  TTsemilep->Fill("g1Mb1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100", 0.0);
  TTsemilep->Fill("0Lb1Ll", 0.0);
  TTsemilep->Fill("0Lbg1Y1Ll", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT", 0.0);


  TTdilep->Fill("NoCuts", 0.0);
  TTdilep->Fill("Cleaning", 0.0);
  TTdilep->Fill("Pileup", 0.0);
  TTdilep->Fill("HCAL_noise", 0.0);
  TTdilep->Fill("vertexg0", 0.0);
  TTdilep->Fill("njetge3", 0.0);
  TTdilep->Fill("HLT", 0.0);
  TTdilep->Fill("jet1ptg200", 0.0);
  TTdilep->Fill("SIG", 0.0); // MR>800 R2>0.08
  TTdilep->Fill("neleeq0", 0.0);
  TTdilep->Fill("nmueq0", 0.0);
  TTdilep->Fill("trackIso", 0.0);
  TTdilep->Fill("g1Mb0Ll", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTdilep->Fill("0Lb0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTdilep->Fill("1Ll", 0.0);
  TTdilep->Fill("g1Mb1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100", 0.0);
  TTdilep->Fill("0Lb1Ll", 0.0);
  TTdilep->Fill("0Lbg1Y1Ll", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT", 0.0);

  // -----------------------------------------------
  // -- extra histograms when running over signal --
  // -----------------------------------------------
  
  // Binning of the signal samples
  int nbins_mother = 10;
  int nbins_LSP = 10;
  int mother_min = 0; 
  int mother_max = 1000; 
  int LSP_min = 0; 
  int LSP_max = 500; 

  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T1t1t"){
    // mother is gluino
    nbins_mother = 15;
    nbins_LSP = 11;
    mother_min = 600; 
    mother_max = 1350; 
    LSP_min = 0; 
    LSP_max = 550; 
  } else if (sample == "T1ttcc_old"){
    // "mother" is stop as scan has fixed gluino mass
    nbins_mother = 113;
    nbins_LSP = 6;
    mother_min = 310; 
    mother_max = 880; 
    LSP_min = 300; 
    LSP_max = 900; 
  } else if (sample == "T2tt"){
    // mother is stop
    nbins_mother = 35;
    nbins_LSP = 37;
    mother_min = 150; 
    mother_max = 1025; 
    LSP_min = 0; 
    LSP_max = 925; 
  }

  // We need one histogram per region, per mass point
  TH2D* list_S[nbins_mother][nbins_LSP];
  TH2D* list_Q[nbins_mother][nbins_LSP];
  TH2D* list_T[nbins_mother][nbins_LSP];
  TH2D* list_W[nbins_mother][nbins_LSP];

  TH2D* list_S_uw[nbins_mother][nbins_LSP];
  TH2D* list_Q_uw[nbins_mother][nbins_LSP];
  TH2D* list_T_uw[nbins_mother][nbins_LSP];
  TH2D* list_W_uw[nbins_mother][nbins_LSP];

  int step_mother = (mother_max - mother_min)/nbins_mother;
  int step_LSP = (LSP_max - LSP_min)/nbins_LSP;
  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
      || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
    int counter_i = 0;
    int counter_j = 0;
    for(int i=mother_min; i<mother_max; i+=step_mother){
      counter_j = 0;
      for(int j=LSP_min; j<LSP_max; j+=step_LSP){
	TString nameS = "h_S_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameT = "h_T_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameQ = "h_Q_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameW = "h_W_" + sample + "_" + to_string(i) + "_" + to_string(j);
	list_S[counter_i][counter_j] = new TH2D(nameS,nameS,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_T[counter_i][counter_j] = new TH2D(nameT,nameT,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_Q[counter_i][counter_j] = new TH2D(nameQ,nameQ,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_W[counter_i][counter_j] = new TH2D(nameW,nameW,nbins_MR,bins_MR,nbins_R2,bins_R2);

	TString nameSuw = "h_uw_S_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameTuw = "h_uw_T_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameQuw = "h_uw_Q_" + sample + "_" + to_string(i) + "_" + to_string(j);
	TString nameWuw = "h_uw_W_" + sample + "_" + to_string(i) + "_" + to_string(j);
	list_S_uw[counter_i][counter_j] = new TH2D(nameSuw,nameSuw,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_T_uw[counter_i][counter_j] = new TH2D(nameTuw,nameTuw,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_Q_uw[counter_i][counter_j] = new TH2D(nameQuw,nameQuw,nbins_MR,bins_MR,nbins_R2,bins_R2);
	list_W_uw[counter_i][counter_j] = new TH2D(nameWuw,nameWuw,nbins_MR,bins_MR,nbins_R2,bins_R2);

	counter_j++;
      }
      counter_i++;
    }
  }


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  //  nevents = 100;
  for(int entry=0; entry < nevents; ++entry)
    {
      // Read event into memory
      stream.read(entry);
      cout << "Event: " << entry << endl;    
      
      // Count events and get the total weight contibuted by the event
      h_totweight->Fill(1, geneventinfoproduct_weight);
      double w = 1.;
      if (geneventinfoproduct_weight != 0) {
        w = geneventinfoproduct_weight*weightnorm;
      }
      
      // Write every ith event:
      if (entry % 10000 == 0) cout << entry << endl;
      
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file rzrBTanalyzer.h to find out what structs
      // are available. Each struct contains the field "selected", which
      // can be set as needed. Call saveSelectedObjects() before a call to
      // addEvent if you wish to save only the selected objects.
      
      fillObjects();

      // get mass point information for signal samples
      double mg = lheeventproducthelper_mg;
      double mt1 = lheeventproducthelper_mt1;
      double mz1 = lheeventproducthelper_mz1;
      double m_mother = mt1; // set mother to stop, works for T2tt and T1ttcc_old
      if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T1t1t")
	m_mother = mg;
      if (sample == "T2tt" && mz1 == 0) continue; // MLSP=0 should be rejected for this sample

      // get normalization for sms's
      if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T2tt" || sample == "T1t1t"){
	int bin_mother = (m_mother - mother_min)/step_mother;
	int bin_LSP = (mz1 - LSP_min)/step_LSP;
	w = w/h_smscounts->GetBinContent(bin_mother+1,bin_LSP+1);	
      }


      //cout << "will fill ofile with weight " << w << endl;
      ofile.count("NoCuts", w);

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

      ofile.count("Cleaning", w);

      // -------------------------
      // -- Trigger requirement --
      // -------------------------

      if (eventhelper_isRealData==1) { 
	// Run2012A:
	if ( eventhelper_run >= 190456 && eventhelper_run < 190762 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 1 ||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v5 == 1 ||
		 //triggerresultshelper_HLT_HT750_v1 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v3 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 0 &&
		 (triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v5 == 1) &&
		  //triggerresultshelper_HLT_HT750_v1 == 1 &&
		 triggerresultshelper_HLT_DiPFJetAve320_v3 == 0 
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 1 ||
		 triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v5 == 1 ||
		 //triggerresultshelper_HLT_HT750_v1 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v3 == 1
		 ) ) continue;
	  }
	}

	if ( eventhelper_run >= 190762 && eventhelper_run < 191512 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1 ||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v6 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v4 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 0 &&
		 (triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v6 == 1) &&
		 //triggerresultshelper_HLT_HT750_v2 == 1 &&
		 triggerresultshelper_HLT_DiPFJetAve320_v4 == 0
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1 ||
		 triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v6 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v4 == 1
		 ) ) continue;
	  }
	}

	if ( eventhelper_run >= 191512 && eventhelper_run < 193746 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 ||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v7 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v5 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 0 &&
		 (triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v7 == 1) &&
		 //triggerresultshelper_HLT_HT750_v2 == 1 &&
		 triggerresultshelper_HLT_DiPFJetAve320_v5 == 0
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 ||
		 triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v7 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 triggerresultshelper_HLT_DiPFJetAve320_v5 == 1
		 ) ) continue;
	  }
	}

	// 2012B
	if ( eventhelper_run >= 193746 && eventhelper_run < 196039) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5 == 1 ||
	       triggerresultshelper_HLT_PFHT650_v8 == 1 ||
	       //triggerresultshelper_HLT_HT750_v3 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v6 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 196039 && eventhelper_run < 197770 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6 == 1 ||
	       triggerresultshelper_HLT_PFHT650_v9 == 1 ||
	       //triggerresultshelper_HLT_HT750_v4 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v6 == 1
	       ) ) continue;
	}

	// Run2012C:
	if ( eventhelper_run >= 197770 && eventhelper_run < 199648 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet400_v6 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v1 == 1 ||
	       //triggerresultshelper_HLT_HT750_v5 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v7 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 199648 && eventhelper_run < 202807 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v8 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v3 == 1 ||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v9 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 202807 && eventhelper_run < 203734 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v4 == 1 ||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v10 == 1
	       ) ) continue;
	}

	// Run2012D:
	if ( eventhelper_run >= 203734 && eventhelper_run < 208940 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1 ||
	       triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v4 == 1 ||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       triggerresultshelper_HLT_DiPFJetAve400_v10 == 1
	       ) ) continue;
	}
      }


      // ------------------------
      // -- Pileup reweighting --
      // ------------------------

      double num_vertices = pileupsummaryinfo[0].getTrueNumInteractions;
      // get the bin number in the pileup histogram
      int pileup_bin = (int)ceil(num_vertices);
      // get the nominal, up and down weights
      double w_pileup = 1.;
      double w_pileup_up = 1.;
      double w_pileup_down = 1.;
      if(doPileupReweighting){
	w_pileup = h_pileup->GetBinContent(pileup_bin);      
	w_pileup_up = h_pileup_up->GetBinContent(pileup_bin);      
	w_pileup_down = h_pileup_down->GetBinContent(pileup_bin);      
	double dw_pileup_up = max(w_pileup_up, w_pileup_down) - w_pileup;
	double dw_pileup_down = w_pileup - min(w_pileup_up, w_pileup_down);
	//cout << "PU: " << w_pileup << " " << dw_pileup_up << " " << dw_pileup_down << endl;
	if (sigmaPU >= 0.) {
	  w_pileup += sigmaPU*dw_pileup_up; 
	} else {
	  w_pileup += sigmaPU*dw_pileup_down;
	}
      }
      // Compute the pileup weight according to the systematic variation considered
      // Use difference between nominal and up/down as 1 sigma variation 
      w = w*w_pileup;

      ofile.count("Pileup",w);


      cout << "xpu: " << endl;

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
      double METcorrfromJEC_px = 0.0;
      double METcorrfromJEC_py = 0.0;
      for (unsigned int i=0; i<cmgpfjet.size(); i++) {
	// begin JEC SF
	double jecUnc = 0.;
	AK5PFCHSJECunc(cmgpfjet[i].pt, cmgpfjet[i].eta, jecUnc);
	double jecSF = 1. + (sigmaJEC * fabs(jecUnc));
	//cout << cmgpfjet[i].pt << " " << cmgpfjet[i].eta << " " << jecUnc << endl;
	// Put the jet in a TLorentzVector and scale it with JEC SF
	TLorentzVector jlnojecSF;
	jlnojecSF.SetPtEtaPhiE(cmgpfjet[i].pt, cmgpfjet[i].eta,
			cmgpfjet[i].phi, cmgpfjet[i].energy);
	TLorentzVector jl;
	jl.SetPxPyPzE(jlnojecSF.Px()*jecSF, jlnojecSF.Py()*jecSF, jlnojecSF.Pz()*jecSF, jlnojecSF.E()*jecSF);
	//cout << "SF, jlpt, jleta, corrpt, correta: " << jecSF << " " << jlnojecSF.Pt() << " " << jlnojecSF.Eta() 
	//     << " " << jl.Pt() << " " << jl.Eta() << " " << jl.Pt()/jl.Pt() << endl;
	// get the MET corrections:
	METcorrfromJEC_px += (jl.Px() - jlnojecSF.Px());
	METcorrfromJEC_py += (jl.Py() - jlnojecSF.Py());
	cmgpfjet[i].pt = jl.Pt();
	// end of JEC SF
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
        btagCSVMSFFull(partonFlavour, pt, fabs(eta), SFCSVMFl, dSFCSVMFl);
        btagCSVMSFFast(partonFlavour, pt, fabs(eta), SFCSVMFs, dSFCSVMFs);

        btagCSVLSFFull(partonFlavour, pt, fabs(eta), SFCSVLFl, dSFCSVLFl);
        btagCSVLSFFast(partonFlavour, pt, fabs(eta), SFCSVLFs, dSFCSVLFs);

	double eCSVM = 0, eCSVL = 0;
	double SFCSVM, SFCSVL;
	// FastSim:
	if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
	    || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t") {
	  if (fabs(partonFlavour) == 5) {
	    eCSVM = geteff1D(h_pt_b_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_b_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) == 4) {
	    eCSVM = geteff1D(h_pt_c_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_c_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) != 4 and fabs(partonFlavour != 5)) {
	    eCSVM = geteff2D(h_pt_eta_l_CSVMeff, pt, fabs(eta));
	    eCSVL = geteff2D(h_pt_eta_l_CSVLeff, pt, fabs(eta));
	  }
	  SFCSVL = (SFCSVLFl + sigmabtagFl*dSFCSVLFl)*(SFCSVLFs + sigmabtagFs*dSFCSVLFs);
	  SFCSVM = (SFCSVMFl + sigmabtagFl*dSFCSVMFl)*(SFCSVMFs + sigmabtagFs*dSFCSVMFs);
	} else { // FullSim
	  if (fabs(partonFlavour) == 5) {
	    eCSVM = geteff1D(h_pt_b_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_b_CSVLeff, pt);
	  }
	  if (partonFlavour != 5) {
	    eCSVM = geteff1D(h_pt_lc_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_lc_CSVLeff, pt);
	  }
	  SFCSVL = (SFCSVLFl + sigmabtagFl*dSFCSVLFl);
	  SFCSVM = (SFCSVMFl + sigmabtagFs*dSFCSVMFl);
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
	//TLorentzVector jl;
	//jl.SetPtEtaPhiE(cmgpfjet[i].pt, cmgpfjet[i].eta,
	//		cmgpfjet[i].phi, cmgpfjet[i].energy);
	LVsjet.push_back(jl);
      }

      double wCSVM = PCSVMdata / PCSVMsim;
      double wCSVL = PCSVLdata / PCSVLsim;

      //cout << sbjet.size() << " " << wCSVM << " " << slbjet.size() << " " << wCSVL << endl;
      //continue;


      // CA8
      // W selection:
      std::vector<jethelper4_s> sjet2;
      std::vector<jethelper4_s> sW;
      std::vector<jethelper4_s> aW;
      std::vector<jethelper4_s> sY;
      for (unsigned int i=0; i<jethelper4.size(); i++) {
	// begin JEC SF
	double jecUnc = 0.;
	AK5PFCHSJECunc(jethelper4[i].pt, jethelper4[i].eta, jecUnc);
	double jecSFCA8 = 1. + (sigmaJECCA8 * fabs(jecUnc));
	// Put the jet in a TLorentzVector and scale it with JEC SF
	TLorentzVector jlCA8nojecSF;
	jlCA8nojecSF.SetPtEtaPhiE(jethelper4[i].pt, jethelper4[i].eta,
			jethelper4[i].phi, jethelper4[i].energy);
	TLorentzVector jlCA8;
	jlCA8.SetPxPyPzE(jlCA8nojecSF.Px()*jecSFCA8, jlCA8nojecSF.Py()*jecSFCA8, jlCA8nojecSF.Pz()*jecSFCA8, jlCA8nojecSF.E()*jecSFCA8);
	jethelper4[i].pt = jlCA8.Pt();
        if (!(jethelper4[i].pt > 30) ) continue;
        if (!(fabs(jethelper4[i].eta) < 2.4) ) continue;
	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
        //if (!(jethelper4[i].mass > 65 && jethelper4[i].mass < 105)) continue;
	sY.push_back(jethelper4[i]);
        sjet2.push_back(jethelper4[i]);
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
	//if (tau21 >= 0.46 || tau3 >= 0.135) {
	if (tau21 >= 0.50) {
          aW.push_back(jethelper4[i]);
        }
	if (!(tau21 < 0.50) ) continue;
	//if (!(tau21 < 0.46) ) continue;
	//if (!(tau3 < 0.135) ) continue;
        sW.push_back(jethelper4[i]);
      }

      // Wtag scale factors:
      double SFWtag = 0.86;
      double dSFWtag = 0.08;
      double w_Wtag = pow((SFWtag + sigmaW*dSFWtag), sW.size());


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
	// veto 1.442 < |eta| < 1.556?
	if (!(fabs(cmgelectron[i].eta) < 1.442 && fabs(cmgelectron[i].eta) < 1.556) ) continue;
	velectron.push_back(cmgelectron[i]);
      }
      // electron SFs:
      double SFeleVeto = 1.;
      double dSFeleVeto = 0.;
      if (velectron.size() == 1) {
	eleLooseSFEFull(velectron[0].pt, fabs(velectron[0].eta), SFeleVeto, dSFeleVeto);
      }
      double w_eleVeto = SFeleVeto + sigmaeleVeto*dSFeleVeto;
      //cout << "ele: " << SFe << " " << dSFe << endl;

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

      // We are not using taus, so we comment this out:
      /*
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
      */


      // ---------------------
      // -- Razor variables --
      // ---------------------

      // Calculate MR and R2 ignoring muons
      TVector3 V3metnojecSF;
      V3metnojecSF.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      TLorentzVector metnojecSF;
      metnojecSF.SetPtEtaPhiE(cmgbasemet2[0].pt, 0, cmgbasemet2[0].phi, cmgbasemet2[0].energy);
      // MET with JEC SF
      TVector3 V3met;
      V3met.SetXYZ(V3metnojecSF.Px()-METcorrfromJEC_px, V3metnojecSF.Py()-METcorrfromJEC_py, 0.0);
      TLorentzVector met;
      met.SetPxPyPzE(V3met.Px(), V3met.Py(), 0.0, V3met.Pt());
      std::vector<TLorentzVector> LVhemis = CombineJets(LVsjet);

      //cout << METcorrfromJEC_px << " " << METcorrfromJEC_py << " " 
      //   << V3metnojecSF.Pt() << " " << metnojecSF.E() << " " 
      //   << V3met.Pt() << " " << met.E() << endl;

      double MR = -9999;
      double MTR = -9999;
      double R2 = -9999;
      if (LVhemis.size() == 2) {
	MR = CalcMR(LVhemis[0], LVhemis[1]);
	if (MR != MR) continue;
	MTR = CalcMTR(LVhemis[0], LVhemis[1], V3met);
	R2 = pow((MTR / MR),2);
      }

      // We are not using these definitions which were once efined for HLT studies.
      /*      
      // Calculate MR and R2 adding mus to MET
      TVector3 V3metmu;
      V3metmu.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      //cout << V3metmu.Px() << " " << V3metmu.Py() << " " << V3metmu.Phi() << endl;
      for (unsigned int i=0; i<V3mu.size(); i++) {
        V3metmu += V3mu[i];
	//cout << i << " " << V3mu[i].Px() << " " << V3mu[i].Py() << " " << V3mu[i].Phi() << endl;
	//cout << i << "  " << V3metmu.Px() << " " << V3metmu.Py() << " " << V3metmu.Phi() << endl;
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
      */


      // --------------------------------------------------------
      // -- Everything computed from generator level particles --
      // --------------------------------------------------------

      std::vector<genparticlehelper_s> tops;
      std::vector<genparticlehelper_s> Ws;      
      for (unsigned int i=0; i<genparticlehelper.size(); i++) {
        if (genparticlehelper[i].status != 3) continue;
        if (fabs(genparticlehelper[i].pdgId) == 6) {
          std::vector<genparticlehelper_s> topdaughters;
          int id1 = genparticlehelper[i].firstDaughter;
          int id2 = genparticlehelper[i].lastDaughter;
          topdaughters.push_back(genparticlehelper[id1]);
          topdaughters.push_back(genparticlehelper[id2]);
          if ((fabs(topdaughters[0].pdgId) == 5 ||
               fabs(topdaughters[0].pdgId) == 24) &&
              (fabs(topdaughters[1].pdgId) == 5 ||
               fabs(topdaughters[1].pdgId) == 24)   
              ) {
            tops.push_back(genparticlehelper[i]);
            for (unsigned int j=0; j<topdaughters.size(); j++) {
              if (fabs(topdaughters[j].pdgId)==24) {
                Ws.push_back(topdaughters[j]);
	      }
            }
          }
	}
      }
 
      // also get the composition of the decays of ttbar
      bool isTTallhad = false;
      bool isTTsemilep = false;
      bool isTTdilep = false;

      if (Ws.size() == 2) {
	int nlep = 0;
	for (unsigned int i=0; i<Ws.size(); i++) {
	  int iWd1 = Ws[i].firstDaughter;
	  int iWd2 = Ws[i].lastDaughter;
	  if (fabs(genparticlehelper[iWd1].pdgId) < 5 ||
	      fabs(genparticlehelper[iWd2].pdgId) < 5 ) {
	  } else {
	    nlep++;
	  }
	}
	if (nlep==0)
	  isTTallhad = true;
	else if(nlep==1)
	  isTTsemilep = true;
	else if(nlep==2)
	  isTTdilep = true;
      }
      if(isTTallhad)
	TTallhad->Fill("Cleaning",w);
      else if(isTTsemilep)
	TTsemilep->Fill("Cleaning",w);
      else if(isTTdilep)
	TTdilep->Fill("Cleaning",w);


      // ------------------------------------
      // -- Pick up trigger weights for MC --
      // ------------------------------------

      double w_trigger = 1.;
      double hlterr = 0.;
      if (eventhelper_isRealData==0) {
	if (sjet.size() > 0){ // need at least one jet
	  for (int i=1; i<h_hlteff->GetNbinsX()+1; i++) {
	    double xmin = h_hlteff->GetXaxis()->GetBinLowEdge(i);
	    double xmax = h_hlteff->GetXaxis()->GetBinUpEdge(i);
	    if (!(HT >= xmin && HT < xmax)) continue;
	    for (int j=1; j<h_hlteff->GetNbinsY()+1; j++) {
	      double ymin = h_hlteff->GetYaxis()->GetBinLowEdge(j);
	      double ymax = h_hlteff->GetYaxis()->GetBinUpEdge(j);
	      if (sjet[0].pt >= ymin && sjet[0].pt < ymax) {
		w_trigger = h_hlteff->GetBinContent(i, j);
		if (sigmaHLT >= 0) { 
		  hlterr = sigmaHLT*h_hlteff_up->GetBinContent(i, j);
		} else {
		  hlterr = sigmaHLT*h_hlteff_low->GetBinContent(i, j);
		}
		w_trigger += hlterr;
		if (w_trigger > 1.) w_trigger = 1.;
		if (w_trigger < 0.) w_trigger = 0.;
		//cout << i << " " << j << " " << h_hlteff->GetBinContent(i, j) << " " << h_hlteff_up->GetBinContent(i, j) << " " << h_hlteff_low->GetBinContent(i, j) << endl;
		break;
	      }
	    }
	  }
	}
      }

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
	  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T1ttcc_old" || sample == "T1t1t")
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
	  double dw_ISR_up = w_ISR_up - w_ISR_nominal;
	  double dw_ISR_down = w_ISR_nominal - w_ISR_down;
	  if (sigmaISR >= 0.) {
	    w_ISR_nominal += sigmaISR*dw_ISR_up;
	  } else {
	    w_ISR_nominal += sigmaISR*dw_ISR_down;
	  }
	}

      // Need to think how to apply these weights with their systematic variations
      if (doISRreweighting)
	w = w*w_ISR_nominal;
     
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

	double dw_TopPt_up = w_TopPt_up - w_TopPt_nominal;
	double dw_TopPt_down = w_TopPt_nominal - w_TopPt_down;
	//cout << "toppt: " << w_TopPt_down << " " << w_TopPt_nominal << " " << w_TopPt_up << endl;
	if (sigmaTopPt >= 0.) {
	  w_TopPt_nominal += sigmaTopPt*dw_TopPt_up;
	} else {
	  w_TopPt_nominal += sigmaTopPt*dw_TopPt_down;
	}
      }

      if(doTopPtreweighting)
	w = w*w_TopPt_nominal;
      

      cout << "xbefore sel" << endl;

      // ---------------------
      // -- event selection --
      // ---------------------


      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;
      ofile.count("HCAL_noise", w);
      if(isTTallhad)
	TTallhad->Fill("HCAL_noise", w);
      else if(isTTsemilep)
	TTsemilep->Fill("HCAL_noise", w);
      else if(isTTdilep)
	TTdilep->Fill("HCAL noise", w);

      // at least one good primary vertex
      if (!(svertex.size() > 0)) continue;
      ofile.count("vertexg0", w);
      if(isTTallhad)
	TTallhad->Fill("vertexg0", w);
      else if(isTTsemilep)
	TTsemilep->Fill("vertexg0", w);
      else if(isTTdilep)
	TTdilep->Fill("vertexg0", w);

      // at least three jets
      if (!(sjet.size() >= 3)) continue;
      ofile.count("njetge3", w);
      if(isTTallhad)
	TTallhad->Fill("njetge3", w);
      else if(isTTsemilep)
	TTsemilep->Fill("njetge3", w);
      else if(isTTdilep)
	TTdilep->Fill("njetge3", w);
      
      // Apply the HLT weight and include it in the total weight:
      double w_nohlt = w;
      w = w*w_trigger;

      ofile.count("HLT", w);
      if(isTTallhad)
	TTallhad->Fill("HLT", w);
      else if(isTTsemilep)
	TTsemilep->Fill("HLT", w);
      else if(isTTdilep)
	TTdilep->Fill("HLT", w);


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

      
      // pt of first jet greater than 200 GeV
      if (!(sjet[0].pt > 200)) continue;
      ofile.count("jet1ptg200", w);
      if(isTTallhad)
	TTallhad->Fill("jet1ptg200", w);
      else if(isTTsemilep)
	TTsemilep->Fill("jet1ptg200", w);
      else if(isTTdilep)
	TTdilep->Fill("jet1ptg200", w);

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
	ofile.count("SIG", w);

	if(isTTallhad)
	  TTallhad->Fill("SIG", w);
	else if(isTTsemilep)
	  TTsemilep->Fill("SIG", w);
	else if(isTTdilep)
	  TTdilep->Fill("SIG", w);
	
	// ----------------------------------------------------------------------------------------------------
	// 0 Lepton trajectory
	// ----------------------------------------------------------------------------------------------------
	if (nlooseelectrons == 0){
	  ofile.count("neleeq0", w);
	  if(isTTallhad)
	    TTallhad->Fill("neleeq0", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("neleeq0", w);
	  else if(isTTdilep)
	    TTdilep->Fill("neleeq0", w);
	  
	  if (nloosemuons == 0) {
	    ofile.count("nmueq0", w);
	    if(isTTallhad)
	      TTallhad->Fill("nmueq0", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("nmueq0", w);
	    else if(isTTdilep)
	      TTdilep->Fill("nmueq0", w);
	    
	    if (eventhelperextra_trackIso == 0){
	      ofile.count("trackIso", w);
	      if(isTTallhad)
		TTallhad->Fill("trackIso", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("trackIso", w);
	      else if(isTTdilep)
		TTdilep->Fill("trackIso", w);
	      
	      if (nmediumbs > 0){
		if (eventhelper_isRealData!=1) {
		  w = w*wCSVM;
		}
		ofile.count("g1Mb0Ll", w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mb0Ll", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mb0Ll", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mb0Ll", w);

		// g1Mb g1W 0Ll -- SIGNAL region
		if( sW.size() > 0){
		  if (eventhelper_isRealData!=1) {
		    w = w*w_Wtag;
		  }
		  ofile.count("g1Mbg1W0Ll",w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W0Ll", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W0Ll", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W0Ll", w);

		  if ( minDeltaPhiHat > 4){
		    ofile.count("g1Mbg1W0Ll_mdPhiHatg4",w);
		    h_MR_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, R2, w);
		    h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, R2, 1.);
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);

		    if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
			|| sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
		      int bin_mother = (m_mother - mother_min)/step_mother;
		      int bin_LSP = (mz1 - LSP_min)/step_LSP;
		      list_S[bin_mother][bin_LSP]->Fill(MR,R2,w);
		      list_S_uw[bin_mother][bin_LSP]->Fill(MR,R2,1.);
		    }
		 		    
		  } // end of  minDeltaPhiHat > 4
		} // end of sW.size() > 0
	      } // end of nmediumbs > 0
	      
	      if (nloosebs == 0){
		if (eventhelper_isRealData!=1) {
		  w = w*wCSVL;
		}
		ofile.count("0Lb0Ll", w);
		if(isTTallhad)
		  TTallhad->Fill("0Lb0Ll", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lb0Ll", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lb0Ll", w);
		
		// 0Lbg1uW0Ll -- QCD control region
		if( aW.size() > 0){
		  ofile.count("0Lbg1uW0Ll",w);
		  if(isTTallhad)
		    TTallhad->Fill("0Lbg1uW0Ll", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("0Lbg1uW0Ll", w);
		  else if(isTTdilep)
		    TTdilep->Fill("0Lbg1uW0Ll", w);
		  
		  if (minDeltaPhiHat < 4){
		    ofile.count("0Lbg1uW0Ll_mdPhiHat4",w);
		    h_MR_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, R2, 1.);
		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat4", w);

		    if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
			|| sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
		      int bin_mother = (m_mother - mother_min)/step_mother;
		      int bin_LSP = (mz1 - LSP_min)/step_LSP;
		      list_Q[bin_mother][bin_LSP]->Fill(MR,R2,w);
		      list_Q_uw[bin_mother][bin_LSP]->Fill(MR,R2,1.);
		    }
		  } // end of minDeltaPhiHat < 4

		} // end of aW.size() > 0
				
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
	  if (nlooseelectrons == 1) {
	    lepton.SetPtEtaPhiE(velectron[0].pt, velectron[0].eta, velectron[0].phi, velectron[0].energy);
	    if (eventhelper_isRealData!=1) {
	      w = w*w_eleVeto;
	    }
	  } else if (nloosemuons == 1)
	    lepton.SetPtEtaPhiE(vmuon[0].pt, vmuon[0].eta, vmuon[0].phi, vmuon[0].energy);
	  double mT = CalcMT(lepton,met);
	  
	  ofile.count("1Ll",w);
	  if(isTTallhad)
	    TTallhad->Fill("1Ll", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("1Ll", w);
	  else if(isTTdilep)
	    TTdilep->Fill("1Ll", w);
	  
	  if (nmediumbs > 0){
	    if (eventhelper_isRealData!=1) {
	      w = w*wCSVM;
	    }
	    ofile.count("g1Mb1Ll",w);
	    if(isTTallhad)
	      TTallhad->Fill("g1Mb1Ll", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("g1Mb1Ll", w);
	    else if(isTTdilep)
	      TTdilep->Fill("g1Mb1Ll", w);
	    
	    if( sW.size() > 0 ){
	      if (eventhelper_isRealData!=1) {
		w = w*w_Wtag;
	      }
	      ofile.count("g1Mbg1W1Ll",w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mbg1W1Ll", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mbg1W1Ll", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mbg1W1Ll", w);

	      // TTJets Control region
	      if (mT < 100){
		ofile.count("g1Mbg1W1LlmT100",w);
		h_MR_g1Mbg1W1LlmT100->Fill(MR, w);
		h_R2_g1Mbg1W1LlmT100->Fill(R2, w);
		h_MR_R2_g1Mbg1W1LlmT100->Fill(MR, R2, w);
		h_uw_MR_R2_g1Mbg1W1LlmT100->Fill(MR, R2, 1.);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1W1LlmT100", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1W1LlmT100", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1W1LlmT100", w);
		
		if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
		    || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
		  int bin_mother = (m_mother - mother_min)/step_mother;
		  int bin_LSP = (mz1 - LSP_min)/step_LSP;
		  list_T[bin_mother][bin_LSP]->Fill(MR,R2,w);
		  list_T_uw[bin_mother][bin_LSP]->Fill(MR,R2,1.);
		}
	      } // end mT < 100
	    } // end sW.size()
	  } // end nmediumbs > 0


	  if (nloosebs == 0){
	    if (eventhelper_isRealData!=1) {
	      w = w*wCSVL;
	    }
	    ofile.count("0Lb1Ll",w);
	    if(isTTallhad)
	      TTallhad->Fill("0Lb1Ll", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("0Lb1Ll", w);
	    else if(isTTdilep)
	      TTdilep->Fill("0Lb1Ll", w);
	    
	    if( sY.size() > 0 ){
	      ofile.count("0Lbg1Y1Ll",w);
	      if(isTTallhad)
		TTallhad->Fill("0Lbg1Y1Ll", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lbg1Y1Ll", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lbg1Y1Ll", w);
	      
	      // WJets Control Region
	      if (mT < 100 && mT > 30){ 
		ofile.count("0Lbg1Y1LlmT",w);
		h_MR_0Lbg1Y1LlmT->Fill(MR, w);
		h_R2_0Lbg1Y1LlmT->Fill(R2, w);
		h_MR_R2_0Lbg1Y1LlmT->Fill(MR, R2, w);
		h_uw_MR_R2_0Lbg1Y1LlmT->Fill(MR, R2, 1.);
		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y1LlmT", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y1LlmT", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y1LlmT", w);

		if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80"
		    || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t"){
		  int bin_mother = (m_mother - mother_min)/step_mother;
		  int bin_LSP = (mz1 - LSP_min)/step_LSP;
		  list_W[bin_mother][bin_LSP]->Fill(MR,R2,w);
		  list_W_uw[bin_mother][bin_LSP]->Fill(MR,R2,1.);
		}
	      } // end mT < 100 && mT > 30
	    } // end sY.size()
	  } // end nloosebs > 0

	} // end nlooseleptons == 1
      } // end of MR>800 R2>0.08
      

    } // end event loop
  
  fhlt->Close();
  fpileup->Close();
  fbeff->Close();
  stream.close();
  ofile.close();
  return 0;
}
