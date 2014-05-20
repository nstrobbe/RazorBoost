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

  // ----------------------------------------------------------------
  // -- Get the trigger histogram: // old name: extr_eff0_sm2.root --
  // ----------------------------------------------------------------

  TFile* fhlt; 
  fhlt = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  //if (!fhlt){
    fhlt = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  //}
  if (!fhlt){
    cout << "Could not find trigger efficiency root file... Where did you put it??" << endl;
    return 1;
  }
  //TH2D* h_hlteff = (TH2D*)fhlt->Get("hBinValues");
  TH2D* h_hlteff = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph");
    
  // ------------------------------
  // -- Parse command line stuff --
  // ------------------------------

  // Get file list and histogram filename from command line
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples
  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) 
    error("unable to open ntuple file(s)");

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

  // parse rest of command line arguments
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
  string Runs = "";
  if ( argc > 10 )
    Runs = string(argv[10]);

  bool doISRreweighting = false;
  if (ISR == "ISR_True" 
      && (sample == "T2tt" || sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" 
	  || sample == "T1ttcc_DM80" || sample == "T1ttcc_old" || sample == "T1t1t"
	  || sample == "TTJets" || sample == "WJets" || sample == "ZJets" )
      ){
    doISRreweighting = true;
    cout << "Will do ISR reweighting" << endl;
  }
  bool doTopPtreweighting = false;
  if (sample == "TTJets" && TopPt == "TopPt_True"){
    doTopPtreweighting = true;
    doISRreweighting = false; // Never do both ISR and TopPt reweighting
    cout << "Will do top pt reweighting and no ISR reweighting" << endl;
  }
  bool doPileupReweighting = false;
  if (sample != "Data" && Pileup == "Pileup_True"){
    doPileupReweighting = true;
    cout << "Will do pileup reweighting" << endl;
  }

  // ---------------------------------------
  // --- Get the correct pileup histogram --
  // ---------------------------------------

  TString pileupname = "pileup_weights.root"; // default
  if (Runs == "AB"){
    pileupname = "pileup_weights_AB.root";
    cout << "Using pileup profile for runs A+B only" << endl;
  }
  if (sample == "T1ttcc_old" || sample == "T2tt"){
    pileupname = "pileup_weights_sig52X.root";
    if (Runs == "AB")
      pileupname = "pileup_weights_AB_sig52X.root";
  }
  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" || sample == "T1t1t"){
    pileupname = "pileup_weights_sig53X.root";
    if (Runs == "AB")
      pileupname = "pileup_weights_AB_sig53X.root";
  }
  TFile* fpileup;
  //fpileup = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname);
  //if (!fpileup){
    fpileup = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/pileup/"+pileupname);
  //}
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
    //fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_TTJets.root");
    fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_TTJets.root");
  } else if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
	     || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t") {
    //fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_T1ttcc.root");
    fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_T1ttcc.root");
  } else {
    //fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_QCD.root");
    fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_QCD.root");
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

  // ------------------------------
  // -- Get the systematics file --
  // ------------------------------
  /*
  ifstream systFile("systematics_nominal.txt");
  if ( !systFile.good() ) error("unable to open systematics file");
  std::vector<double> vsyst;
  double syst;
  while (systFile >> syst)
    cout << "Systematic: " << syst << endl;
    vsyst.push_back(syst);
  */

  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights --
  // -- The original MC weight will be divided by this quantity  --
  // --------------------------------------------------------------

  double weightnorm = 1.;
  if (xsect != -1 && totweight != -1 && lumi != -1) { // i.e. is not data
    weightnorm = (xsect*lumi)/totweight;
  }

  cout << "lumi: " << lumi << endl;
  cout << "xsect: " << xsect << endl;
  cout << "totweight: " << totweight << endl;
  cout << "weightnorm: " << weightnorm << endl;

  outputFile ofile(cmdline.outputfilename);

  // ---------------------------------------------------------------------------
  // -- Declare histograms                                                    --
  // ---------------------------------------------------------------------------

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

  // for variable bin width in MR, R2
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

  // variable binning for trigger histogram
  int nbn_HT = 19;
  int nbn_j1pt = 12;
  Double_t bn_HT_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,900.,1000.,2500.};
  Double_t* bn_HT = getVariableBinEdges(nbn_HT+1,bn_HT_tmp);
  Double_t bn_j1pt_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,700.,1000.};
  Double_t* bn_j1pt = getVariableBinEdges(nbn_j1pt+1,bn_j1pt_tmp);

  // Save the total weight for cross check purposes
  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);

  // histograms for total ISR weight: Sum(all events) w_ISR
  // needed to do the reweighting properly without changing the overall cross section. 
  TH1D* h_totalISRweight_nominal = new TH1D("h_totalISRweight_nominal", "h_totalISRweight_nominal", 1, 1, 2); // nominal ISR weight
  TH1D* h_totalISRweight_up = new TH1D("h_totalISRweight_up", "h_totalISRweight_up", 1, 1, 2); //  ISR up weight
  TH1D* h_totalISRweight_down = new TH1D("h_totalISRweight_down", "h_totalISRweight_down", 1, 1, 2); // ISR down weight

  // histograms for top pt reweighting
  TH1D* h_totalTopPTweight_nominal = new TH1D("h_totalTopPTweight_nominal", "h_totalTopPTweight_nominal", 1, 1, 2); // nominal ISR weight
  TH1D* h_totalTopPTweight_up = new TH1D("h_totalTopPTweight_up", "h_totalTopPTweight_up", 1, 1, 2); //  ISR up weight
  TH1D* h_totalTopPTweight_down = new TH1D("h_totalTopPTweight_down", "h_totalTopPTweight_down", 1, 1, 2); // ISR down weight

  // Histograms for pileup reweighting
  TH1D* h_TrueNumVertices = new TH1D("h_TrueNumVertices","h_TrueNumVertices",100,0,100);
  TH1D* h_TrueNumVertices_reweighted = new TH1D("h_TrueNumVertices_reweighted","h_TrueNumVertices_reweighted",100,0,100);
  TH1D* h_TrueNumVertices_reweighted_withtrigger = new TH1D("h_TrueNumVertices_reweighted_withtrigger","h_TrueNumVertices_reweighted_withtrigger",100,0,100);

  TH1D* h_PV_nosel = new TH1D("h_PV_nosel","h_PV_nosel",50,0,50);
  TH1D* h_PV_reweighted_nosel = new TH1D("h_PV_reweighted_nosel","h_PV_reweighted_nosel",50,0,50);
  TH1D* h_PV_1j = new TH1D("h_PV_1j","h_PV_1j",50,0,50);
  TH1D* h_PV_reweighted_1j = new TH1D("h_PV_reweighted_1j","h_PV_reweighted_1j",50,0,50);
  TH1D* h_PV_HLT = new TH1D("h_PV_HLT","h_PV_HLT",50,0,50);
  TH1D* h_PV_reweighted_HLT = new TH1D("h_PV_reweighted_HLT","h_PV_reweighted_HLT",50,0,50);
  TH1D* h_PV_SIG = new TH1D("h_PV_SIG","h_PV_SIG",50,0,50);
  TH1D* h_PV_reweighted_SIG = new TH1D("h_PV_reweighted_SIG","h_PV_reweighted_SIG",50,0,50);
  TH1D* h_PV_g1Mbg1W1LlmT100 = new TH1D("h_PV_g1Mbg1W1LlmT100","h_PV_g1Mbg1W1LlmT100",50,0,50);
  TH1D* h_PV_reweighted_g1Mbg1W1LlmT100 = new TH1D("h_PV_reweighted_g1Mbg1W1LlmT100","h_PV_reweighted_g1Mbg1W1LlmT100",50,0,50);

  TH2D* h_PV_HT_nosel = new TH2D("h_PV_HT_nosel","h_PV_HT_nosel",50,0,50, 50,0,2000);
  TH2D* h_PV_HT_reweighted_nosel = new TH2D("h_PV_HT_reweighted_nosel","h_PV_HT_reweighted_nosel",50,0,50, 50,0,2000);

  // Cross check plots for trigger weights
  TH2D * h_HT_jet1pt_nosel = new TH2D("h_HT_jet1pt_nosel", "h_HT_jet1pt_nosel", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  TH2D * h_HT_jet1pt_nosel_notrigw = new TH2D("h_HT_jet1pt_nosel_notrigw", "h_HT_jet1pt_nosel_notrigw", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  // W tagging plots
  TH2D* h_jmass_jpt = new TH2D("h_jmass_jpt", "h_jmass_jpt", 20, 0, 350, 20, 0, 1400);
  TH2D* h_d1pt_d2pt = new TH2D("h_d1pt_d2pt", "h_d1pt_d2pt", 20, 0, 1400, 20, 0, 1400);
  TH2D* h_d1m_d2m = new TH2D("h_d1m_d2m", "h_d1m_d2m", 20, 0, 200, 20, 0, 200);
  TH1D* h_jmass = new TH1D("h_jmass", "h_jmass", 50, 0, 250);
  TH2D* h_d1ptsel_d2ptsel = new TH2D("h_d1ptsel_d2ptsel", "h_d1ptsel_d2ptsel", 20, 0, 1400, 20, 0, 1400);
  TH2D* h_d1msel_d2msel = new TH2D("h_d1msel_d2msel", "h_d1msel_d2msel", 20, 0, 200, 20, 0, 200);
  TH1D* h_mdrp = new TH1D("h_mdrp", "h_mdrp", 50, 0, 1);
  TH2D* h_mdrp_jpt = new TH2D("h_mdrp_jpt", "h_mdrp_jpt", 20, 0, 1, 20, 0, 1400);
  TH1D* h_mdrp2 = new TH1D("h_mdrp2", "h_mdrp2", 50, 0, 1);
  TH2D* h_mdrp2_jpt = new TH2D("h_mdrp2_jpt", "h_mdrp2_jpt", 20, 0, 1, 20, 0, 1400);
  TH1D* h_yasym = new TH1D("h_yasym", "h_yasym", 50, 0, 1);
  TH2D* h_mdrp_yasym = new TH2D("h_mdrp_yasym", "h_mdrp_yasym", 20, 0, 1, 20, 0, 1);
  TH2D* h_mdrp2_yasym = new TH2D("h_mdrp2_yasym", "h_mdrp2_yasym", 20, 0, 1, 20, 0, 1);

  TH1D* h_pt_closestCA8jet = new TH1D("h_pt_closestCA8jet","h_pt_closestCA8jet",50,0,1000);
  TH1D* h_pt_closestCA8jet_Wmass = new TH1D("h_pt_closestCA8jet_Wmass","h_pt_closestCA8jet_Wmass",50,0,1000);
  TH1D* h_pt_closestCA8jet_Wmass_tagged = new TH1D("h_pt_closestCA8jet_Wmass_tagged","h_pt_closestCA8jet_Wmass_tagged",50,0,1000);

  // Plots of relevant quantities with tag for the cut flow
  TH1D* h_nW_Pileup = new TH1D("h_nW_Pileup", "h_nW_Pileup", 5, 0, 5);
  TH1D* h_nb_Pileup = new TH1D("h_nb_Pileup", "h_nb_Pileup", 5, 0, 5);
  TH2D* h_nW_nb_Pileup = new TH2D("h_nW_nb_Pileup", "h_nW_nb_Pileup", 5, 0, 5, 5, 0, 5);
  TH1D* h_nWAK5_Pileup = new TH1D("h_nWAK5_Pileup", "h_nWAK5_Pileup", 5, 0, 5);
  TH2D* h_met_R2_Pileup = new TH2D("h_met_R2_Pileup", "h_met_R2_Pileup", 20, 0, 1000, 25, 0, 1);
  TH1D* h_met_Pileup = new TH1D("h_met_Pileup", "h_met_Pileup", 20, 0, 1000);
  TH2D* h_metmu_R2metmu_Pileup = new TH2D("h_metmu_R2metmu_Pileup", "h_metmu_R2metmu_Pileup", 20, 0, 1000, 25, 0, 1);
  TH1D* h_metmu_Pileup = new TH1D("h_metmu_Pileup", "h_metmu_Pileup", 20, 0, 1000);
  TH2D* h_metel_R2metel_Pileup = new TH2D("h_metel_R2metel_Pileup", "h_metel_R2metel_Pileup", 20, 0, 1000, 25, 0, 1);
  TH1D* h_metel_Pileup = new TH1D("h_metel_Pileup", "h_metel_Pileup", 20, 0, 1000);
  TH1D* h_j1pt_Pileup = new TH1D("h_j1pt_Pileup", "h_j1pt_Pileup", 20, 0, 1000);
  TH1D* h_nj_Pileup = new TH1D("h_nj_Pileup", "h_nj_Pileup", 20, 0, 20);

  TH1D* h_nW_jet1ptg200 = new TH1D("h_nW_jet1ptg200", "h_nW_jet1ptg200", 5, 0, 5);
  TH1D* h_nb_jet1ptg200 = new TH1D("h_nb_jet1ptg200", "h_nb_jet1ptg200", 5, 0, 5);
  TH2D* h_nW_nb_jet1ptg200 = new TH2D("h_nW_nb_jet1ptg200", "h_nW_nb_jet1ptg200", 5, 0, 5, 5, 0, 5);
  TH1D* h_nWAK5_jet1ptg200 = new TH1D("h_nWAK5_jet1ptg200", "h_nWAK5_jet1ptg200", 5, 0, 5);

  // MR, R2 plots for the different steps in the selection
  // need at least two jets to be able to compute MR and R2

  TH1D* h_MR_Pileup = new TH1D("h_MR_Pileup", "h_MR_Pileup", nbins_MR, bins_MR);
  TH1D* h_R2_Pileup = new TH1D("h_R2_Pileup", "h_R2_Pileup", nbins_R2, bins_R2);
  TH2D* h_MR_R2_Pileup = new TH2D("h_MR_R2_Pileup", "h_MR_R2_Pileup", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_Pileup = new TH1D("h_PV_Pileup", "h_PV_Pileup", 50, 0, 50);

  TH1D* h_R2metmu_Pileup = new TH1D("h_R2metmu_Pileup", "h_R2metmu_Pileup", nbins_R2, bins_R2);
  TH2D* h_MR_R2metmu_Pileup = new TH2D("h_MR_R2metmu_Pileup", "h_MR_R2metmu_Pileup", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_R2metel_Pileup = new TH1D("h_R2metel_Pileup", "h_R2metel_Pileup", nbins_R2, bins_R2);
  TH2D* h_MR_R2metel_Pileup = new TH2D("h_MR_R2metel_Pileup", "h_MR_R2metel_Pileup", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D* h_MR_HCAL_noise = new TH1D("h_MR_HCAL_noise", "h_MR_HCAL_noise", nbins_MR, bins_MR);
  TH1D* h_R2_HCAL_noise = new TH1D("h_R2_HCAL_noise", "h_R2_HCAL_noise", nbins_R2, bins_R2);
  TH2D* h_MR_R2_HCAL_noise = new TH2D("h_MR_R2_HCAL_noise", "h_MR_R2_HCAL_noise", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_HCAL_noise = new TH1D("h_PV_HCAL_noise", "h_PV_HCAL_noise", 50, 0, 50);

  TH1D* h_MR_vertexg0 = new TH1D("h_MR_vertexg0", "h_MR_vertexg0", nbins_MR, bins_MR);
  TH1D* h_R2_vertexg0 = new TH1D("h_R2_vertexg0", "h_R2_vertexg0", nbins_R2, bins_R2);
  TH2D* h_MR_R2_vertexg0 = new TH2D("h_MR_R2_vertexg0", "h_MR_R2_vertexg0", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_vertexg0 = new TH1D("h_PV_vertexg0", "h_PV_vertexg0", 50, 0, 50);

  TH1D* h_MR_njetge3 = new TH1D("h_MR_njetge3", "h_MR_njetge3", nbins_MR, bins_MR);
  TH1D* h_R2_njetge3 = new TH1D("h_R2_njetge3", "h_R2_njetge3", nbins_R2, bins_R2);
  TH2D* h_MR_R2_njetge3 = new TH2D("h_MR_R2_njetge3", "h_MR_R2_njetge3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_njetge3 = new TH1D("h_PV_njetge3", "h_PV_njetge3", 50, 0, 50);

  TH1D* h_MR_HLT = new TH1D("h_MR_HLT", "h_MR_HLT", nbins_MR, bins_MR);
  TH1D* h_R2_HLT = new TH1D("h_R2_HLT", "h_R2_HLT", nbins_R2, bins_R2);
  TH2D* h_MR_R2_HLT = new TH2D("h_MR_R2_HLT", "h_MR_R2_HLT", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_minDeltaPhi_HLT = new TH1D("h_minDeltaPhi_HLT", "h_minDeltaPhi_HLT", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_HLT = new TH2D("h_MR_minDeltaPhi_HLT", "h_MR_minDeltaPhi_HLT", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_HLT = new TH2D("h_R2_minDeltaPhi_HLT", "h_R2_minDeltaPhi_HLT", nbins_R2, bins_R2, 50, 0, 5);
  TH1D * h_minDeltaPhiHat_HLT = new TH1D("h_minDeltaPhiHat_HLT", "h_minDeltaPhiHat_HLT", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_HLT = new TH2D("h_MR_minDeltaPhiHat_HLT", "h_MR_minDeltaPhiHat_HLT", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_HLT = new TH2D("h_R2_minDeltaPhiHat_HLT", "h_R2_minDeltaPhiHat_HLT", nbins_R2, bins_R2, 30, 0, 15);

  TH1D* h_MR_jet1ptg200 = new TH1D("h_MR_jet1ptg200", "h_MR_jet1ptg200", nbins_MR, bins_MR);
  TH1D* h_R2_jet1ptg200 = new TH1D("h_R2_jet1ptg200", "h_R2_jet1ptg200", nbins_R2, bins_R2);
  TH2D* h_MR_R2_jet1ptg200 = new TH2D("h_MR_R2_jet1ptg200", "h_MR_R2_jet1ptg200", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_jet1ptg200 = new TH1D("h_PV_jet1ptg200", "h_PV_jet1ptg200", 50, 0, 50);

  TH1D * h_minDeltaPhi_jet1ptg200 = new TH1D("h_minDeltaPhi_jet1ptg200", "h_minDeltaPhi_jet1ptg200", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_jet1ptg200 = new TH2D("h_MR_minDeltaPhi_jet1ptg200", "h_MR_minDeltaPhi_jet1ptg200", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_jet1ptg200 = new TH2D("h_R2_minDeltaPhi_jet1ptg200", "h_R2_minDeltaPhi_jet1ptg200", nbins_R2, bins_R2, 50, 0, 5); 
  TH1D * h_minDeltaPhiHat_jet1ptg200 = new TH1D("h_minDeltaPhiHat_jet1ptg200", "h_minDeltaPhiHat_jet1ptg200", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_jet1ptg200 = new TH2D("h_MR_minDeltaPhiHat_jet1ptg200", "h_MR_minDeltaPhiHat_jet1ptg200", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_jet1ptg200 = new TH2D("h_R2_minDeltaPhiHat_jet1ptg200", "h_R2_minDeltaPhiHat_jet1ptg200", nbins_R2, bins_R2, 30, 0, 15);

  TH1D * h_njets_jet1ptg200 = new TH1D("h_njets_jet1ptg200","h_njets_jet1ptg200",15,0,15);
  TH1D * h_nbjets_jet1ptg200 = new TH1D("h_nbjets_jet1ptg200","h_nbjets_jet1ptg200",6,0,6);
  TH1D * h_met_jet1ptg200 = new TH1D("h_met_jet1ptg200","h_met_jet1ptg200",20,0,1000);
  TH1D * h_jet1pt_jet1ptg200 = new TH1D("h_jet1pt_jet1ptg200","h_jet1pt_jet1ptg200",20,0,1000);
  TH1D * h_jet2pt_jet1ptg200 = new TH1D("h_jet2pt_jet1ptg200","h_jet2pt_jet1ptg200",20,0,1000);
  TH1D * h_jet3pt_jet1ptg200 = new TH1D("h_jet3pt_jet1ptg200","h_jet3pt_jet1ptg200",20,0,1000);

  TH1D * h_MR_SIG = new TH1D("h_MR_SIG", "h_MR_SIG", nbins_MR, bins_MR);
  TH1D * h_R2_SIG = new TH1D("h_R2_SIG", "h_R2_SIG", nbins_R2, bins_R2);
  TH2D * h_MR_R2_SIG = new TH2D("h_MR_R2_SIG", "h_MR_R2_SIG", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_minDeltaPhi_SIG = new TH1D("h_minDeltaPhi_SIG", "h_minDeltaPhi_SIG", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_SIG = new TH2D("h_MR_minDeltaPhi_SIG", "h_MR_minDeltaPhi_SIG", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_SIG = new TH2D("h_R2_minDeltaPhi_SIG", "h_R2_minDeltaPhi_SIG", nbins_R2, bins_R2, 50, 0, 5);
  TH1D * h_minDeltaPhiHat_SIG = new TH1D("h_minDeltaPhiHat_SIG", "h_minDeltaPhiHat_SIG", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_SIG = new TH2D("h_MR_minDeltaPhiHat_SIG", "h_MR_minDeltaPhiHat_SIG", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_SIG = new TH2D("h_R2_minDeltaPhiHat_SIG", "h_R2_minDeltaPhiHat_SIG", nbins_R2, bins_R2, 30, 0, 15);

  TH1D * h_HT_SIG = new TH1D("h_HT_SIG", "h_HT_SIG", 100, 0, 2000);
  TH2D * h_HT_jet1pt_SIG = new TH2D("h_HT_jet1pt_SIG", "h_HT_jet1pt_SIG", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  TH2D * h_HT_jet1pt_SIG_notrigw = new TH2D("h_HT_jet1pt_SIG_notrigw", "h_HT_jet1pt_SIG_notrigw", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  // 0 lepton trajectory
  TH1D * h_MR_neleeq0 = new TH1D("h_MR_neleeq0", "h_MR_neleeq0", nbins_MR, bins_MR);
  TH1D * h_R2_neleeq0 = new TH1D("h_R2_neleeq0", "h_R2_neleeq0", nbins_R2, bins_R2);
  TH2D * h_MR_R2_neleeq0 = new TH2D("h_MR_R2_neleeq0", "h_MR_R2_neleeq0", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_nmueq0 = new TH1D("h_MR_nmueq0", "h_MR_nmueq0", nbins_MR, bins_MR);
  TH1D * h_R2_nmueq0 = new TH1D("h_R2_nmueq0", "h_R2_nmueq0", nbins_R2, bins_R2);
  TH2D * h_MR_R2_nmueq0 = new TH2D("h_MR_R2_nmueq0", "h_MR_R2_nmueq0", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_trackIso = new TH1D("h_MR_trackIso", "h_MR_trackIso", nbins_MR, bins_MR);
  TH1D * h_R2_trackIso = new TH1D("h_R2_trackIso", "h_R2_trackIso", nbins_R2, bins_R2);
  TH2D * h_MR_R2_trackIso = new TH2D("h_MR_R2_trackIso", "h_MR_R2_trackIso", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_minDeltaPhi_0Ll = new TH1D("h_minDeltaPhi_0Ll", "h_minDeltaPhi_0Ll", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_0Ll = new TH2D("h_MR_minDeltaPhi_0Ll", "h_MR_minDeltaPhi_0Ll", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_0Ll = new TH2D("h_R2_minDeltaPhi_0Ll", "h_R2_minDeltaPhi_0Ll", nbins_R2, bins_R2, 50, 0, 5);
  TH1D * h_minDeltaPhiHat_0Ll = new TH1D("h_minDeltaPhiHat_0Ll", "h_minDeltaPhiHat_0Ll", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_0Ll = new TH2D("h_MR_minDeltaPhiHat_0Ll", "h_MR_minDeltaPhiHat_0Ll", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_0Ll = new TH2D("h_R2_minDeltaPhiHat_0Ll", "h_R2_minDeltaPhiHat_0Ll", nbins_R2, bins_R2, 30, 0, 15);
  TH1D* h_PV_0Ll = new TH1D("h_PV_0Ll", "h_PV_0Ll", 50, 0, 50);

  // g1Mb 0Ll
  TH1D * h_MR_g1Mb0Ll = new TH1D("h_MR_g1Mb0Ll", "h_MR_g1Mb0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb0Ll = new TH1D("h_R2_g1Mb0Ll", "h_R2_g1Mb0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb0Ll = new TH2D("h_MR_R2_g1Mb0Ll", "h_MR_R2_g1Mb0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_g1Mb0Ll = new TH1D("h_PV_g1Mb0Ll", "h_PV_g1Mb0Ll", 50, 0, 50);

  // g1Mb g1W 0Ll ; Signal box: >= 1 Mb; >= 1 W; 0 Ll
  TH1D * h_MR_g1Mbg1W0Ll = new TH1D("h_MR_g1Mbg1W0Ll", "h_MR_g1Mbg1W0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll = new TH1D("h_R2_g1Mbg1W0Ll", "h_R2_g1Mbg1W0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll = new TH2D("h_MR_R2_g1Mbg1W0Ll", "h_MR_R2_g1Mbg1W0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1W0Ll = new TH1D("h_njets_g1Mbg1W0Ll","h_njets_g1Mbg1W0Ll",15,0,15);
  TH1D * h_nbjets_g1Mbg1W0Ll = new TH1D("h_nbjets_g1Mbg1W0Ll","h_nbjets_g1Mbg1W0Ll",6,0,6);
  TH1D * h_met_g1Mbg1W0Ll = new TH1D("h_met_g1Mbg1W0Ll","h_met_g1Mbg1W0Ll",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W0Ll = new TH1D("h_jet1pt_g1Mbg1W0Ll","h_jet1pt_g1Mbg1W0Ll",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W0Ll = new TH1D("h_jet2pt_g1Mbg1W0Ll","h_jet2pt_g1Mbg1W0Ll",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W0Ll = new TH1D("h_jet3pt_g1Mbg1W0Ll","h_jet3pt_g1Mbg1W0Ll",20,0,1000);

  TH1D * h_dphimegajets_g1Mbg1W0Ll = new TH1D("h_dphimegajets_g1Mbg1W0Ll","h_dphimegajets_g1Mbg1W0Ll",40,0,4);
  TH1D * h_dphijet1jet2_g1Mbg1W0Ll = new TH1D("h_dphijet1jet2_g1Mbg1W0Ll","h_dphijet1jet2_g1Mbg1W0Ll",40,0,4);
  TH1D * h_HT_g1Mbg1W0Ll = new TH1D("h_HT_g1Mbg1W0Ll", "h_HT_g1Mbg1W0Ll", 100, 0, 2000);
  TH1D * h_PV_g1Mbg1W0Ll = new TH1D("h_PV_g1Mbg1W0Ll", "h_PV_g1Mbg1W0Ll", 50, 0, 50);
  TH1D * h_minDeltaPhiHat_g1Mbg1W0Ll = new TH1D("h_minDeltaPhiHat_g1Mbg1W0Ll", "h_minDeltaPhiHat_g1Mbg1W0Ll", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_g1Mbg1W0Ll = new TH2D("h_MR_minDeltaPhiHat_g1Mbg1W0Ll", "h_MR_minDeltaPhiHat_g1Mbg1W0Ll", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_g1Mbg1W0Ll = new TH2D("h_R2_minDeltaPhiHat_g1Mbg1W0Ll", "h_R2_minDeltaPhiHat_g1Mbg1W0Ll", nbins_R2, bins_R2, 30, 0, 15);

  TH1D * h_minDeltaPhi_g1Mbg1W0Ll = new TH1D("h_minDeltaPhi_g1Mbg1W0Ll", "h_minDeltaPhi_g1Mbg1W0Ll", 50, 0, 1);
  TH2D * h_MR_minDeltaPhi_g1Mbg1W0Ll = new TH2D("h_MR_minDeltaPhi_g1Mbg1W0Ll", "h_MR_minDeltaPhi_g1Mbg1W0Ll", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_g1Mbg1W0Ll = new TH2D("h_R2_minDeltaPhi_g1Mbg1W0Ll", "h_R2_minDeltaPhi_g1Mbg1W0Ll", nbins_R2, bins_R2, 50, 0, 5);

  // separate 1b and >=2b plots
  TH1D * h_MR_1Mbg1W0Ll = new TH1D("h_MR_1Mbg1W0Ll", "h_MR_1Mbg1W0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_1Mbg1W0Ll = new TH1D("h_R2_1Mbg1W0Ll", "h_R2_1Mbg1W0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_1Mbg1W0Ll = new TH2D("h_MR_R2_1Mbg1W0Ll", "h_MR_R2_1Mbg1W0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_g2Mbg1W0Ll = new TH1D("h_MR_g2Mbg1W0Ll", "h_MR_g2Mbg1W0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g2Mbg1W0Ll = new TH1D("h_R2_g2Mbg1W0Ll", "h_R2_g2Mbg1W0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g2Mbg1W0Ll = new TH2D("h_MR_R2_g2Mbg1W0Ll", "h_MR_R2_g2Mbg1W0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // signal box + mindeltaphihat > 4
  TH1D * h_MR_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhiHatg4", "h_MR_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4 = new TH2D("h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", "h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_njets_g1Mbg1W0Ll_mdPhiHatg4","h_njets_g1Mbg1W0Ll_mdPhiHatg4",15,0,15);
  TH1D * h_nbjets_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_nbjets_g1Mbg1W0Ll_mdPhiHatg4","h_nbjets_g1Mbg1W0Ll_mdPhiHatg4",6,0,6);
  TH1D * h_met_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_met_g1Mbg1W0Ll_mdPhiHatg4","h_met_g1Mbg1W0Ll_mdPhiHatg4",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_jet1pt_g1Mbg1W0Ll_mdPhiHatg4","h_jet1pt_g1Mbg1W0Ll_mdPhiHatg4",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_jet2pt_g1Mbg1W0Ll_mdPhiHatg4","h_jet2pt_g1Mbg1W0Ll_mdPhiHatg4",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W0Ll_mdPhiHatg4 = new TH1D("h_jet3pt_g1Mbg1W0Ll_mdPhiHatg4","h_jet3pt_g1Mbg1W0Ll_mdPhiHatg4",20,0,1000);

  // signal box + mindeltaphihat < 4
  TH1D * h_MR_g1Mbg1W0Ll_mdPhiHat4 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhiHat4", "h_MR_g1Mbg1W0Ll_mdPhiHat4", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhiHat4 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhiHat4", "h_R2_g1Mbg1W0Ll_mdPhiHat4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhiHat4 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhiHat4", "h_MR_R2_g1Mbg1W0Ll_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // signal box + mindeltaphi > 0.3
  TH1D * h_MR_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhig0p3", "h_MR_g1Mbg1W0Ll_mdPhig0p3", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhig0p3", "h_R2_g1Mbg1W0Ll_mdPhig0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhig0p3 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhig0p3", "h_MR_R2_g1Mbg1W0Ll_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p3 = new TH2D("h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p3", "h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_njets_g1Mbg1W0Ll_mdPhig0p3","h_njets_g1Mbg1W0Ll_mdPhig0p3",15,0,15);
  TH1D * h_nbjets_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_nbjets_g1Mbg1W0Ll_mdPhig0p3","h_nbjets_g1Mbg1W0Ll_mdPhig0p3",6,0,6);
  TH1D * h_met_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_met_g1Mbg1W0Ll_mdPhig0p3","h_met_g1Mbg1W0Ll_mdPhig0p3",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_jet1pt_g1Mbg1W0Ll_mdPhig0p3","h_jet1pt_g1Mbg1W0Ll_mdPhig0p3",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_jet2pt_g1Mbg1W0Ll_mdPhig0p3","h_jet2pt_g1Mbg1W0Ll_mdPhig0p3",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W0Ll_mdPhig0p3 = new TH1D("h_jet3pt_g1Mbg1W0Ll_mdPhig0p3","h_jet3pt_g1Mbg1W0Ll_mdPhig0p3",20,0,1000);

  // signal box + mindeltaphi < 0.3
  TH1D * h_MR_g1Mbg1W0Ll_mdPhi0p3 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhi0p3", "h_MR_g1Mbg1W0Ll_mdPhi0p3", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhi0p3 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhi0p3", "h_R2_g1Mbg1W0Ll_mdPhi0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhi0p3 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhi0p3", "h_MR_R2_g1Mbg1W0Ll_mdPhi0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // signal box + mindeltaphi > 0.5
  TH1D * h_MR_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhig0p5", "h_MR_g1Mbg1W0Ll_mdPhig0p5", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhig0p5", "h_R2_g1Mbg1W0Ll_mdPhig0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhig0p5 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhig0p5", "h_MR_R2_g1Mbg1W0Ll_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p5 = new TH2D("h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p5", "h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_njets_g1Mbg1W0Ll_mdPhig0p5","h_njets_g1Mbg1W0Ll_mdPhig0p5",15,0,15);
  TH1D * h_nbjets_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_nbjets_g1Mbg1W0Ll_mdPhig0p5","h_nbjets_g1Mbg1W0Ll_mdPhig0p5",6,0,6);
  TH1D * h_met_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_met_g1Mbg1W0Ll_mdPhig0p5","h_met_g1Mbg1W0Ll_mdPhig0p5",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_jet1pt_g1Mbg1W0Ll_mdPhig0p5","h_jet1pt_g1Mbg1W0Ll_mdPhig0p5",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_jet2pt_g1Mbg1W0Ll_mdPhig0p5","h_jet2pt_g1Mbg1W0Ll_mdPhig0p5",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W0Ll_mdPhig0p5 = new TH1D("h_jet3pt_g1Mbg1W0Ll_mdPhig0p5","h_jet3pt_g1Mbg1W0Ll_mdPhig0p5",20,0,1000);

  // signal box + mindeltaphi < 0.5
  TH1D * h_MR_g1Mbg1W0Ll_mdPhi0p5 = new TH1D("h_MR_g1Mbg1W0Ll_mdPhi0p5", "h_MR_g1Mbg1W0Ll_mdPhi0p5", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W0Ll_mdPhi0p5 = new TH1D("h_R2_g1Mbg1W0Ll_mdPhi0p5", "h_R2_g1Mbg1W0Ll_mdPhi0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W0Ll_mdPhi0p5 = new TH2D("h_MR_R2_g1Mbg1W0Ll_mdPhi0p5", "h_MR_R2_g1Mbg1W0Ll_mdPhi0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // g1Mb 0W g1uW 0Ll ; Signal-like box: >= 1 Mb; 0 W - 1 uW; 0 Ll
  TH1D * h_MR_g1Mb0Wg1uW0Ll     = new TH1D("h_MR_g1Mb0Wg1uW0Ll", "h_MR_g1Mb0Wg1uW0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb0Wg1uW0Ll     = new TH1D("h_R2_g1Mb0Wg1uW0Ll", "h_R2_g1Mb0Wg1uW0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb0Wg1uW0Ll  = new TH2D("h_MR_R2_g1Mb0Wg1uW0Ll", "h_MR_R2_g1Mb0Wg1uW0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mb0Wg1uW0Ll  = new TH1D("h_njets_g1Mb0Wg1uW0Ll","h_njets_g1Mb0Wg1uW0Ll",15,0,15);
  TH1D * h_nbjets_g1Mb0Wg1uW0Ll = new TH1D("h_nbjets_g1Mb0Wg1uW0Ll","h_nbjets_g1Mb0Wg1uW0Ll",6,0,6);
  TH1D * h_met_g1Mb0Wg1uW0Ll    = new TH1D("h_met_g1Mb0Wg1uW0Ll","h_met_g1Mb0Wg1uW0Ll",20,0,1000);
  TH1D * h_jet1pt_g1Mb0Wg1uW0Ll = new TH1D("h_jet1pt_g1Mb0Wg1uW0Ll","h_jet1pt_g1Mb0Wg1uW0Ll",20,0,1000);
  TH1D * h_jet2pt_g1Mb0Wg1uW0Ll = new TH1D("h_jet2pt_g1Mb0Wg1uW0Ll","h_jet2pt_g1Mb0Wg1uW0Ll",20,0,1000);
  TH1D * h_jet3pt_g1Mb0Wg1uW0Ll = new TH1D("h_jet3pt_g1Mb0Wg1uW0Ll","h_jet3pt_g1Mb0Wg1uW0Ll",20,0,1000);
  TH1D* h_PV_g1Mb0Wg1uW0Ll = new TH1D("h_PV_g1Mb0Wg1uW0Ll", "h_PV_g1Mb0Wg1uW0Ll", 50, 0, 50);

  // 0Lb 0Ll
  TH1D * h_MR_0Lb0Ll = new TH1D("h_MR_0Lb0Ll", "h_MR_0Lb0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_0Lb0Ll = new TH1D("h_R2_0Lb0Ll", "h_R2_0Lb0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lb0Ll = new TH2D("h_MR_R2_0Lb0Ll", "h_MR_R2_0Lb0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_0Lb0Ll = new TH1D("h_PV_0Lb0Ll", "h_PV_0Lb0Ll", 50, 0, 50);

  TH1D * h_minDeltaPhi_0Lb0Ll = new TH1D("h_minDeltaPhi_0Lb0Ll", "h_minDeltaPhi_0Lb0Ll", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_0Lb0Ll = new TH2D("h_MR_minDeltaPhi_0Lb0Ll", "h_MR_minDeltaPhi_0Lb0Ll", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_0Lb0Ll = new TH2D("h_R2_minDeltaPhi_0Lb0Ll", "h_R2_minDeltaPhi_0Lb0Ll", nbins_R2, bins_R2, 50, 0, 5);
  
  TH1D * h_minDeltaPhiHat_0Lb0Ll = new TH1D("h_minDeltaPhiHat_0Lb0Ll", "h_minDeltaPhiHat_0Lb0Ll", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_0Lb0Ll = new TH2D("h_MR_minDeltaPhiHat_0Lb0Ll", "h_MR_minDeltaPhiHat_0Lb0Ll", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_0Lb0Ll = new TH2D("h_R2_minDeltaPhiHat_0Lb0Ll", "h_R2_minDeltaPhiHat_0Lb0Ll", nbins_R2, bins_R2, 30, 0, 15);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll
  TH1D * h_MR_0Lbg1uW0Ll = new TH1D("h_MR_0Lbg1uW0Ll", "h_MR_0Lbg1uW0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll = new TH1D("h_R2_0Lbg1uW0Ll", "h_R2_0Lbg1uW0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll = new TH2D("h_MR_R2_0Lbg1uW0Ll", "h_MR_R2_0Lbg1uW0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_0Lbg1uW0Ll = new TH1D("h_PV_0Lbg1uW0Ll", "h_PV_0Lbg1uW0Ll", 50, 0, 50);

  TH1D * h_minDeltaPhi_0Lbg1uW0Ll = new TH1D("h_minDeltaPhi_0Lbg1uW0Ll", "h_minDeltaPhi_0Lbg1uW0Ll", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_0Lbg1uW0Ll = new TH2D("h_MR_minDeltaPhi_0Lbg1uW0Ll", "h_MR_minDeltaPhi_0Lbg1uW0Ll", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_0Lbg1uW0Ll = new TH2D("h_R2_minDeltaPhi_0Lbg1uW0Ll", "h_R2_minDeltaPhi_0Lbg1uW0Ll", nbins_R2, bins_R2, 50, 0, 5);
  
  TH1D * h_minDeltaPhiHat_0Lbg1uW0Ll = new TH1D("h_minDeltaPhiHat_0Lbg1uW0Ll", "h_minDeltaPhiHat_0Lbg1uW0Ll", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_0Lbg1uW0Ll = new TH2D("h_MR_minDeltaPhiHat_0Lbg1uW0Ll", "h_MR_minDeltaPhiHat_0Lbg1uW0Ll", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_0Lbg1uW0Ll = new TH2D("h_R2_minDeltaPhiHat_0Lbg1uW0Ll", "h_R2_minDeltaPhiHat_0Lbg1uW0Ll", nbins_R2, bins_R2, 30, 0, 15);
  
  TH1D * h_njets_0Lbg1uW0Ll = new TH1D("h_njets_0Lbg1uW0Ll","h_njets_0Lbg1uW0Ll",15,0,15);
  TH1D * h_met_0Lbg1uW0Ll = new TH1D("h_met_0Lbg1uW0Ll","h_met_0Lbg1uW0Ll",20,0,1000);
  TH1D * h_jet1pt_0Lbg1uW0Ll = new TH1D("h_jet1pt_0Lbg1uW0Ll","h_jet1pt_0Lbg1uW0Ll",20,0,1000);
  TH1D * h_jet2pt_0Lbg1uW0Ll = new TH1D("h_jet2pt_0Lbg1uW0Ll","h_jet2pt_0Lbg1uW0Ll",20,0,1000);
  TH1D * h_jet3pt_0Lbg1uW0Ll = new TH1D("h_jet3pt_0Lbg1uW0Ll","h_jet3pt_0Lbg1uW0Ll",20,0,1000);
  TH1D * h_HT_0Lbg1uW0Ll = new TH1D("h_HT_0Lbg1uW0Ll", "h_HT_0Lbg1uW0Ll", 100, 0, 2000);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhi < 0.3
  TH1D * h_MR_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhi0p3", "h_MR_0Lbg1uW0Ll_mdPhi0p3", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhi0p3", "h_R2_0Lbg1uW0Ll_mdPhi0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhi0p3 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhi0p3", "h_MR_R2_0Lbg1uW0Ll_mdPhi0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p3 = new TH2D("h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p3", "h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_njets_0Lbg1uW0Ll_mdPhi0p3","h_njets_0Lbg1uW0Ll_mdPhi0p3",15,0,15);
  TH1D * h_met_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_met_0Lbg1uW0Ll_mdPhi0p3","h_met_0Lbg1uW0Ll_mdPhi0p3",20,0,1000);
  TH1D * h_jet1pt_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_jet1pt_0Lbg1uW0Ll_mdPhi0p3","h_jet1pt_0Lbg1uW0Ll_mdPhi0p3",20,0,1000);
  TH1D * h_jet2pt_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_jet2pt_0Lbg1uW0Ll_mdPhi0p3","h_jet2pt_0Lbg1uW0Ll_mdPhi0p3",20,0,1000);
  TH1D * h_jet3pt_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_jet3pt_0Lbg1uW0Ll_mdPhi0p3","h_jet3pt_0Lbg1uW0Ll_mdPhi0p3",20,0,1000);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhi > 0.3
  TH1D * h_MR_0Lbg1uW0Ll_mdPhig0p3 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhig0p3", "h_MR_0Lbg1uW0Ll_mdPhig0p3", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhig0p3 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhig0p3", "h_R2_0Lbg1uW0Ll_mdPhig0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhig0p3 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhig0p3", "h_MR_R2_0Lbg1uW0Ll_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhi < 0.5
  TH1D * h_MR_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhi0p5", "h_MR_0Lbg1uW0Ll_mdPhi0p5", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhi0p5", "h_R2_0Lbg1uW0Ll_mdPhi0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhi0p5 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhi0p5", "h_MR_R2_0Lbg1uW0Ll_mdPhi0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p5 = new TH2D("h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p5", "h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_njets_0Lbg1uW0Ll_mdPhi0p5","h_njets_0Lbg1uW0Ll_mdPhi0p5",15,0,15);
  TH1D * h_met_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_met_0Lbg1uW0Ll_mdPhi0p5","h_met_0Lbg1uW0Ll_mdPhi0p5",20,0,1000);
  TH1D * h_jet1pt_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_jet1pt_0Lbg1uW0Ll_mdPhi0p5","h_jet1pt_0Lbg1uW0Ll_mdPhi0p5",20,0,1000);
  TH1D * h_jet2pt_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_jet2pt_0Lbg1uW0Ll_mdPhi0p5","h_jet2pt_0Lbg1uW0Ll_mdPhi0p5",20,0,1000);
  TH1D * h_jet3pt_0Lbg1uW0Ll_mdPhi0p5 = new TH1D("h_jet3pt_0Lbg1uW0Ll_mdPhi0p5","h_jet3pt_0Lbg1uW0Ll_mdPhi0p5",20,0,1000);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhi > 0.5
  TH1D * h_MR_0Lbg1uW0Ll_mdPhig0p5 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhig0p5", "h_MR_0Lbg1uW0Ll_mdPhig0p5", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhig0p5 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhig0p5", "h_R2_0Lbg1uW0Ll_mdPhig0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhig0p5 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhig0p5", "h_MR_R2_0Lbg1uW0Ll_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhiHat < 4 (RA2b value)
  TH1D * h_MR_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhiHat4", "h_MR_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhiHat4", "h_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhiHat4", "h_MR_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4 = new TH2D("h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4", "h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // also get the complement for closure tests
  TH1D * h_MR_0Lbg1uW0Ll_mdPhiHatg4 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhiHatg4", "h_MR_0Lbg1uW0Ll_mdPhiHatg4", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhiHatg4 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhiHatg4", "h_R2_0Lbg1uW0Ll_mdPhiHatg4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhiHatg4 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhiHatg4", "h_MR_R2_0Lbg1uW0Ll_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_njets_0Lbg1uW0Ll_mdPhiHat4","h_njets_0Lbg1uW0Ll_mdPhiHat4",15,0,15);
  TH1D * h_met_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_met_0Lbg1uW0Ll_mdPhiHat4","h_met_0Lbg1uW0Ll_mdPhiHat4",20,0,1000);
  TH1D * h_jet1pt_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_jet1pt_0Lbg1uW0Ll_mdPhiHat4","h_jet1pt_0Lbg1uW0Ll_mdPhiHat4",20,0,1000);
  TH1D * h_jet2pt_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_jet2pt_0Lbg1uW0Ll_mdPhiHat4","h_jet2pt_0Lbg1uW0Ll_mdPhiHat4",20,0,1000);
  TH1D * h_jet3pt_0Lbg1uW0Ll_mdPhiHat4 = new TH1D("h_jet3pt_0Lbg1uW0Ll_mdPhiHat4","h_jet3pt_0Lbg1uW0Ll_mdPhiHat4",20,0,1000);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhiHat < 5
  TH1D * h_MR_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhiHat5", "h_MR_0Lbg1uW0Ll_mdPhiHat5", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhiHat5", "h_R2_0Lbg1uW0Ll_mdPhiHat5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhiHat5 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhiHat5", "h_MR_R2_0Lbg1uW0Ll_mdPhiHat5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_njets_0Lbg1uW0Ll_mdPhiHat5","h_njets_0Lbg1uW0Ll_mdPhiHat5",15,0,15);
  TH1D * h_met_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_met_0Lbg1uW0Ll_mdPhiHat5","h_met_0Lbg1uW0Ll_mdPhiHat5",20,0,1000);
  TH1D * h_jet1pt_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_jet1pt_0Lbg1uW0Ll_mdPhiHat5","h_jet1pt_0Lbg1uW0Ll_mdPhiHat5",20,0,1000);
  TH1D * h_jet2pt_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_jet2pt_0Lbg1uW0Ll_mdPhiHat5","h_jet2pt_0Lbg1uW0Ll_mdPhiHat5",20,0,1000);
  TH1D * h_jet3pt_0Lbg1uW0Ll_mdPhiHat5 = new TH1D("h_jet3pt_0Lbg1uW0Ll_mdPhiHat5","h_jet3pt_0Lbg1uW0Ll_mdPhiHat5",20,0,1000);

  // QCD control region: 0 Lb; >= 1 W; 0 Ll
  TH1D * h_MR_0Lbg1W0Ll = new TH1D("h_MR_0Lbg1W0Ll", "h_MR_0Lbg1W0Ll", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1W0Ll = new TH1D("h_R2_0Lbg1W0Ll", "h_R2_0Lbg1W0Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1W0Ll = new TH2D("h_MR_R2_0Lbg1W0Ll", "h_MR_R2_0Lbg1W0Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);


  // 1 loose lepton trajectory
  TH1D * h_MR_1Ll = new TH1D("h_MR_1Ll", "h_MR_1Ll", nbins_MR, bins_MR);
  TH1D * h_R2_1Ll = new TH1D("h_R2_1Ll", "h_R2_1Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_1Ll = new TH2D("h_MR_R2_1Ll", "h_MR_R2_1Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_HT_1Ll = new TH1D("h_HT_1Ll", "h_HT_1Ll", 100, 0, 2000);
  TH1D* h_PV_1Ll = new TH1D("h_PV_1Ll", "h_PV_1Ll", 50, 0, 50);

  // g1Mb1Ll
  TH1D * h_MR_g1Mb1Ll = new TH1D("h_MR_g1Mb1Ll", "h_MR_g1Mb1Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb1Ll = new TH1D("h_R2_g1Mb1Ll", "h_R2_g1Mb1Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb1Ll = new TH2D("h_MR_R2_g1Mb1Ll", "h_MR_R2_g1Mb1Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_g1Mb1Ll = new TH1D("h_PV_g1Mb1Ll", "h_PV_g1Mb1Ll", 50, 0, 50);

  TH1D * h_mT_g1Mb1Ll = new TH1D("h_mT_g1Mb1Ll", "h_mT_g1Mb1Ll", 50, 0, 500);
  TH1D * h_njets_g1Mb1Ll = new TH1D("h_njets_g1Mb1Ll","h_njets_g1Mb1Ll",15,0,15);
  TH1D * h_nbjets_g1Mb1Ll = new TH1D("h_nbjets_g1Mb1Ll","h_nbjets_g1Mb1Ll",6,0,6);
  TH1D * h_met_g1Mb1Ll = new TH1D("h_met_g1Mb1Ll","h_met_g1Mb1Ll",20,0,1000);
  TH1D * h_jet1pt_g1Mb1Ll = new TH1D("h_jet1pt_g1Mb1Ll","h_jet1pt_g1Mb1Ll",20,0,1000);
  TH1D * h_jet2pt_g1Mb1Ll = new TH1D("h_jet2pt_g1Mb1Ll","h_jet2pt_g1Mb1Ll",20,0,1000);
  TH1D * h_jet3pt_g1Mb1Ll = new TH1D("h_jet3pt_g1Mb1Ll","h_jet3pt_g1Mb1Ll",20,0,1000);
  TH1D * h_leptonpt_g1Mb1Ll = new TH1D("h_leptonpt_g1Mb1Ll","h_leptonpt_g1Mb1Ll",25,0,500);

  // g1Mbg1W1Ll ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll
  TH1D * h_MR_g1Mbg1W1Ll = new TH1D("h_MR_g1Mbg1W1Ll", "h_MR_g1Mbg1W1Ll", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1Ll = new TH1D("h_R2_g1Mbg1W1Ll", "h_R2_g1Mbg1W1Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1Ll = new TH2D("h_MR_R2_g1Mbg1W1Ll", "h_MR_R2_g1Mbg1W1Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_mT_g1Mbg1W1Ll = new TH1D("h_mT_g1Mbg1W1Ll", "h_mT_g1Mbg1W1Ll", 50, 0, 500);
  TH2D * h_MR_mT_g1Mbg1W1Ll = new TH2D("h_MR_mT_g1Mbg1W1Ll", "h_MR_mT_g1Mbg1W1Ll", nbins_MR, bins_MR, 50, 0, 500);
  TH2D * h_R2_mT_g1Mbg1W1Ll = new TH2D("h_R2_mT_g1Mbg1W1Ll", "h_R2_mT_g1Mbg1W1Ll", nbins_R2, bins_R2, 50, 0, 500);

  TH1D * h_njets_g1Mbg1W1Ll = new TH1D("h_njets_g1Mbg1W1Ll","h_njets_g1Mbg1W1Ll",15,0,15);
  TH1D * h_nbjets_g1Mbg1W1Ll = new TH1D("h_nbjets_g1Mbg1W1Ll","h_nbjets_g1Mbg1W1Ll",6,0,6);
  TH1D * h_met_g1Mbg1W1Ll = new TH1D("h_met_g1Mbg1W1Ll","h_met_g1Mbg1W1Ll",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W1Ll = new TH1D("h_jet1pt_g1Mbg1W1Ll","h_jet1pt_g1Mbg1W1Ll",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W1Ll = new TH1D("h_jet2pt_g1Mbg1W1Ll","h_jet2pt_g1Mbg1W1Ll",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W1Ll = new TH1D("h_jet3pt_g1Mbg1W1Ll","h_jet3pt_g1Mbg1W1Ll",20,0,1000);
  TH1D * h_leptonpt_g1Mbg1W1Ll = new TH1D("h_leptonpt_g1Mbg1W1Ll","h_leptonpt_g1Mbg1W1Ll",25,0,500);

  // g1Mbg1W1LlmT ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll; mT<100
  TH1D * h_MR_g1Mbg1W1LlmT100 = new TH1D("h_MR_g1Mbg1W1LlmT100", "h_MR_g1Mbg1W1LlmT100", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100 = new TH1D("h_R2_g1Mbg1W1LlmT100", "h_R2_g1Mbg1W1LlmT100", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100", "h_MR_R2_g1Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W1LlmT100 = new TH2D("h_uw_MR_R2_g1Mbg1W1LlmT100", "h_uw_MR_R2_g1Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1W1LlmT100 = new TH1D("h_njets_g1Mbg1W1LlmT100","h_njets_g1Mbg1W1LlmT100",15,0,15);
  TH1D * h_nbjets_g1Mbg1W1LlmT100 = new TH1D("h_nbjets_g1Mbg1W1LlmT100","h_nbjets_g1Mbg1W1LlmT100",6,0,6);
  TH1D * h_met_g1Mbg1W1LlmT100 = new TH1D("h_met_g1Mbg1W1LlmT100","h_met_g1Mbg1W1LlmT100",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W1LlmT100 = new TH1D("h_jet1pt_g1Mbg1W1LlmT100","h_jet1pt_g1Mbg1W1LlmT100",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W1LlmT100 = new TH1D("h_jet2pt_g1Mbg1W1LlmT100","h_jet2pt_g1Mbg1W1LlmT100",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W1LlmT100 = new TH1D("h_jet3pt_g1Mbg1W1LlmT100","h_jet3pt_g1Mbg1W1LlmT100",20,0,1000);
  TH1D * h_leptonpt_g1Mbg1W1LlmT100 = new TH1D("h_leptonpt_g1Mbg1W1LlmT100","h_leptonpt_g1Mbg1W1LlmT100",25,0,500);
  TH1D * h_HT_g1Mbg1W1LlmT100 = new TH1D("h_HT_g1Mbg1W1LlmT100", "h_HT_g1Mbg1W1LlmT100", 100, 0, 2000);
  TH1D * h_minDeltaPhiHat_g1Mbg1W1LlmT100 = new TH1D("h_minDeltaPhiHat_g1Mbg1W1LlmT100", "h_minDeltaPhiHat_g1Mbg1W1LlmT100", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_g1Mbg1W1LlmT100 = new TH2D("h_MR_minDeltaPhiHat_g1Mbg1W1LlmT100", "h_MR_minDeltaPhiHat_g1Mbg1W1LlmT100", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_g1Mbg1W1LlmT100 = new TH2D("h_R2_minDeltaPhiHat_g1Mbg1W1LlmT100", "h_R2_minDeltaPhiHat_g1Mbg1W1LlmT100", nbins_R2, bins_R2, 30, 0, 15);
  TH1D * h_minDeltaPhi_g1Mbg1W1LlmT100 = new TH1D("h_minDeltaPhi_g1Mbg1W1LlmT100", "h_minDeltaPhi_g1Mbg1W1LlmT100", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_g1Mbg1W1LlmT100 = new TH2D("h_MR_minDeltaPhi_g1Mbg1W1LlmT100", "h_MR_minDeltaPhi_g1Mbg1W1LlmT100", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_g1Mbg1W1LlmT100 = new TH2D("h_R2_minDeltaPhi_g1Mbg1W1LlmT100", "h_R2_minDeltaPhi_g1Mbg1W1LlmT100", nbins_R2, bins_R2, 50, 0, 5);

  // == 1b
  TH1D * h_MR_1Mbg1W1LlmT100 = new TH1D("h_MR_1Mbg1W1LlmT100", "h_MR_1Mbg1W1LlmT100", nbins_MR, bins_MR);
  TH1D * h_R2_1Mbg1W1LlmT100 = new TH1D("h_R2_1Mbg1W1LlmT100", "h_R2_1Mbg1W1LlmT100", nbins_R2, bins_R2);
  TH2D * h_MR_R2_1Mbg1W1LlmT100 = new TH2D("h_MR_R2_1Mbg1W1LlmT100", "h_MR_R2_1Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // >= 2b
  TH1D * h_MR_g2Mbg1W1LlmT100 = new TH1D("h_MR_g2Mbg1W1LlmT100", "h_MR_g2Mbg1W1LlmT100", nbins_MR, bins_MR);
  TH1D * h_R2_g2Mbg1W1LlmT100 = new TH1D("h_R2_g2Mbg1W1LlmT100", "h_R2_g2Mbg1W1LlmT100", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g2Mbg1W1LlmT100 = new TH2D("h_MR_R2_g2Mbg1W1LlmT100", "h_MR_R2_g2Mbg1W1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // mdPhiHatg4
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhiHatg4", "h_MR_g1Mbg1W1LlmT100_mdPhiHatg4", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhiHatg4", "h_R2_g1Mbg1W1LlmT100_mdPhiHatg4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4", "h_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH2D("h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4", "h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_njets_g1Mbg1W1LlmT100_mdPhiHatg4","h_njets_g1Mbg1W1LlmT100_mdPhiHatg4",15,0,15);
  TH1D * h_nbjets_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_nbjets_g1Mbg1W1LlmT100_mdPhiHatg4","h_nbjets_g1Mbg1W1LlmT100_mdPhiHatg4",6,0,6);
  TH1D * h_met_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_met_g1Mbg1W1LlmT100_mdPhiHatg4","h_met_g1Mbg1W1LlmT100_mdPhiHatg4",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_jet1pt_g1Mbg1W1LlmT100_mdPhiHatg4","h_jet1pt_g1Mbg1W1LlmT100_mdPhiHatg4",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_jet2pt_g1Mbg1W1LlmT100_mdPhiHatg4","h_jet2pt_g1Mbg1W1LlmT100_mdPhiHatg4",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_jet3pt_g1Mbg1W1LlmT100_mdPhiHatg4","h_jet3pt_g1Mbg1W1LlmT100_mdPhiHatg4",20,0,1000);
  TH1D * h_leptonpt_g1Mbg1W1LlmT100_mdPhiHatg4 = new TH1D("h_leptonpt_g1Mbg1W1LlmT100_mdPhiHatg4","h_leptonpt_g1Mbg1W1LlmT100_mdPhiHatg4",25,0,500);

  // mdPhiHat4
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhiHat4 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhiHat4", "h_MR_g1Mbg1W1LlmT100_mdPhiHat4", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhiHat4 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhiHat4", "h_R2_g1Mbg1W1LlmT100_mdPhiHat4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhiHat4 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhiHat4", "h_MR_R2_g1Mbg1W1LlmT100_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // mdPhi0p3
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhi0p3 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhi0p3", "h_MR_g1Mbg1W1LlmT100_mdPhi0p3", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhi0p3 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhi0p3", "h_R2_g1Mbg1W1LlmT100_mdPhi0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p3 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p3", "h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  // mdPhig0p3
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhig0p3", "h_MR_g1Mbg1W1LlmT100_mdPhig0p3", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhig0p3", "h_R2_g1Mbg1W1LlmT100_mdPhig0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3", "h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3 = new TH2D("h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3", "h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_njets_g1Mbg1W1LlmT100_mdPhig0p3","h_njets_g1Mbg1W1LlmT100_mdPhig0p3",15,0,15);
  TH1D * h_nbjets_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_nbjets_g1Mbg1W1LlmT100_mdPhig0p3","h_nbjets_g1Mbg1W1LlmT100_mdPhig0p3",6,0,6);
  TH1D * h_met_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_met_g1Mbg1W1LlmT100_mdPhig0p3","h_met_g1Mbg1W1LlmT100_mdPhig0p3",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p3","h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p3",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p3","h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p3",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p3","h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p3",20,0,1000);
  TH1D * h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p3 = new TH1D("h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p3","h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p3",25,0,500);

  // mdPhi0p5
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhi0p5 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhi0p5", "h_MR_g1Mbg1W1LlmT100_mdPhi0p5", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhi0p5 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhi0p5", "h_R2_g1Mbg1W1LlmT100_mdPhi0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p5 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p5", "h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  // mdPhig0p5
  TH1D * h_MR_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_MR_g1Mbg1W1LlmT100_mdPhig0p5", "h_MR_g1Mbg1W1LlmT100_mdPhig0p5", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_R2_g1Mbg1W1LlmT100_mdPhig0p5", "h_R2_g1Mbg1W1LlmT100_mdPhig0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5 = new TH2D("h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5", "h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5 = new TH2D("h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5", "h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_njets_g1Mbg1W1LlmT100_mdPhig0p5","h_njets_g1Mbg1W1LlmT100_mdPhig0p5",15,0,15);
  TH1D * h_nbjets_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_nbjets_g1Mbg1W1LlmT100_mdPhig0p5","h_nbjets_g1Mbg1W1LlmT100_mdPhig0p5",6,0,6);
  TH1D * h_met_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_met_g1Mbg1W1LlmT100_mdPhig0p5","h_met_g1Mbg1W1LlmT100_mdPhig0p5",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p5","h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p5",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p5","h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p5",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p5","h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p5",20,0,1000);
  TH1D * h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p5 = new TH1D("h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p5","h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p5",25,0,500);

  // g1Mbg1W1LlmT ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll; 30<mT<100
  TH1D * h_MR_g1Mbg1W1LlmT = new TH1D("h_MR_g1Mbg1W1LlmT", "h_MR_g1Mbg1W1LlmT", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1W1LlmT = new TH1D("h_R2_g1Mbg1W1LlmT", "h_R2_g1Mbg1W1LlmT", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1W1LlmT = new TH2D("h_MR_R2_g1Mbg1W1LlmT", "h_MR_R2_g1Mbg1W1LlmT", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_mT_g1Mbg1W1LlmT = new TH1D("h_mT_g1Mbg1W1LlmT", "h_mT_g1Mbg1W1LlmT", 50, 0, 500);
  TH2D * h_MR_mT_g1Mbg1W1LlmT = new TH2D("h_MR_mT_g1Mbg1W1LlmT", "h_MR_mT_g1Mbg1W1LlmT", nbins_MR, bins_MR, 50, 0, 500);
  TH2D * h_R2_mT_g1Mbg1W1LlmT = new TH2D("h_R2_mT_g1Mbg1W1LlmT", "h_R2_mT_g1Mbg1W1LlmT", nbins_R2, bins_R2, 50, 0, 500);


  // 0Lb1Ll
  TH1D * h_MR_0Lb1Ll = new TH1D("h_MR_0Lb1Ll", "h_MR_0Lb1Ll", nbins_MR, bins_MR);
  TH1D * h_R2_0Lb1Ll = new TH1D("h_R2_0Lb1Ll", "h_R2_0Lb1Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lb1Ll = new TH2D("h_MR_R2_0Lb1Ll", "h_MR_R2_0Lb1Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_0Lb1Ll = new TH1D("h_PV_0Lb1Ll", "h_PV_0Lb1Ll", 50, 0, 50);

  TH1D * h_mT_0Lb1Ll = new TH1D("h_mT_0Lb1Ll", "h_mT_0Lb1Ll", 50, 0, 500);
  TH2D * h_MR_mT_0Lb1Ll = new TH2D("h_MR_mT_0Lb1Ll", "h_MR_mT_0Lb1Ll", nbins_MR, bins_MR, 50, 0, 500);
  TH2D * h_R2_mT_0Lb1Ll = new TH2D("h_R2_mT_0Lb1Ll", "h_R2_mT_0Lb1Ll", nbins_R2, bins_R2, 50, 0, 500);

  TH1D * h_met_0Lb1Ll = new TH1D("h_met_0Lb1Ll","h_met_0Lb1Ll",20,0,1000);
  TH1D * h_njets_0Lb1Ll = new TH1D("h_njets_0Lb1Ll","h_njets_0Lb1Ll",15,0,15);
  TH1D * h_jet1pt_0Lb1Ll = new TH1D("h_jet1pt_0Lb1Ll","h_jet1pt_0Lb1Ll",20,0,1000);
  TH1D * h_jet2pt_0Lb1Ll = new TH1D("h_jet2pt_0Lb1Ll","h_jet2pt_0Lb1Ll",20,0,1000);
  TH1D * h_jet3pt_0Lb1Ll = new TH1D("h_jet3pt_0Lb1Ll","h_jet3pt_0Lb1Ll",20,0,1000);
  TH1D * h_leptonpt_0Lb1Ll = new TH1D("h_leptonpt_0Lb1Ll","h_leptonpt_0Lb1Ll",25,0,500);

  // 0Lbg1Y1Ll ; Wj control region: == 0 Lb; >= 1 Y; 1 Ll
  TH1D * h_MR_0Lbg1Y1Ll = new TH1D("h_MR_0Lbg1Y1Ll", "h_MR_0Lbg1Y1Ll", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1Ll = new TH1D("h_R2_0Lbg1Y1Ll", "h_R2_0Lbg1Y1Ll", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1Ll = new TH2D("h_MR_R2_0Lbg1Y1Ll", "h_MR_R2_0Lbg1Y1Ll", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D* h_PV_0Lbg1Y1Ll = new TH1D("h_PV_0Lbg1Y1Ll", "h_PV_0Lbg1Y1Ll", 50, 0, 50);

  TH1D * h_mT_0Lbg1Y1Ll = new TH1D("h_mT_0Lbg1Y1Ll", "h_mT_0Lbg1Y1Ll", 50, 0, 500);
  TH2D * h_MR_mT_0Lbg1Y1Ll = new TH2D("h_MR_mT_0Lbg1Y1Ll", "h_MR_mT_0Lbg1Y1Ll", nbins_MR, bins_MR, 50, 0, 500);
  TH2D * h_R2_mT_0Lbg1Y1Ll = new TH2D("h_R2_mT_0Lbg1Y1Ll", "h_R2_mT_0Lbg1Y1Ll", nbins_R2, bins_R2, 50, 0, 500);

  TH1D * h_met_0Lbg1Y1Ll = new TH1D("h_met_0Lbg1Y1Ll","h_met_0Lbg1Y1Ll",20,0,1000);
  TH1D * h_njets_0Lbg1Y1Ll = new TH1D("h_njets_0Lbg1Y1Ll","h_njets_0Lbg1Y1Ll",15,0,15);
  TH1D * h_jet1pt_0Lbg1Y1Ll = new TH1D("h_jet1pt_0Lbg1Y1Ll","h_jet1pt_0Lbg1Y1Ll",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1Ll = new TH1D("h_jet2pt_0Lbg1Y1Ll","h_jet2pt_0Lbg1Y1Ll",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1Ll = new TH1D("h_jet3pt_0Lbg1Y1Ll","h_jet3pt_0Lbg1Y1Ll",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1Ll = new TH1D("h_leptonpt_0Lbg1Y1Ll","h_leptonpt_0Lbg1Y1Ll",25,0,500);

  // 0Lbg1Y1LlmT ; Wj control region: == 1 Lb; >= 1 Y; 1 Ll; mT<100
  TH1D * h_MR_0Lbg1Y1LlmT100 = new TH1D("h_MR_0Lbg1Y1LlmT100", "h_MR_0Lbg1Y1LlmT100", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT100 = new TH1D("h_R2_0Lbg1Y1LlmT100", "h_R2_0Lbg1Y1LlmT100", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT100 = new TH2D("h_MR_R2_0Lbg1Y1LlmT100", "h_MR_R2_0Lbg1Y1LlmT100", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1Y1LlmT100 = new TH1D("h_njets_0Lbg1Y1LlmT100","h_njets_0Lbg1Y1LlmT100",15,0,15);
  TH1D * h_met_0Lbg1Y1LlmT100 = new TH1D("h_met_0Lbg1Y1LlmT100","h_met_0Lbg1Y1LlmT100",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y1LlmT100 = new TH1D("h_jet1pt_0Lbg1Y1LlmT100","h_jet1pt_0Lbg1Y1LlmT100",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1LlmT100 = new TH1D("h_jet2pt_0Lbg1Y1LlmT100","h_jet2pt_0Lbg1Y1LlmT100",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1LlmT100 = new TH1D("h_jet3pt_0Lbg1Y1LlmT100","h_jet3pt_0Lbg1Y1LlmT100",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1LlmT100 = new TH1D("h_leptonpt_0Lbg1Y1LlmT100","h_leptonpt_0Lbg1Y1LlmT100",25,0,500);

  // 0Lbg1Y1LlmT ; Wj control region: == 0 Lb; >= 1 Y; 1 Ll; 30<mT<100
  TH1D * h_MR_0Lbg1Y1LlmT = new TH1D("h_MR_0Lbg1Y1LlmT", "h_MR_0Lbg1Y1LlmT", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT = new TH1D("h_R2_0Lbg1Y1LlmT", "h_R2_0Lbg1Y1LlmT", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT = new TH2D("h_MR_R2_0Lbg1Y1LlmT", "h_MR_R2_0Lbg1Y1LlmT", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1Y1LlmT = new TH2D("h_uw_MR_R2_0Lbg1Y1LlmT", "h_uw_MR_R2_0Lbg1Y1LlmT", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_mT_0Lbg1Y1LlmT = new TH1D("h_mT_0Lbg1Y1LlmT", "h_mT_0Lbg1Y1LlmT", 50, 0, 500);
  TH2D * h_MR_mT_0Lbg1Y1LlmT = new TH2D("h_MR_mT_0Lbg1Y1LlmT", "h_MR_mT_0Lbg1Y1LlmT", nbins_MR, bins_MR, 50, 0, 500);
  TH2D * h_R2_mT_0Lbg1Y1LlmT = new TH2D("h_R2_mT_0Lbg1Y1LlmT", "h_R2_mT_0Lbg1Y1LlmT", nbins_R2, bins_R2, 50, 0, 500);

  TH1D * h_njets_0Lbg1Y1LlmT = new TH1D("h_njets_0Lbg1Y1LlmT","h_njets_0Lbg1Y1LlmT",15,0,15);
  TH1D * h_met_0Lbg1Y1LlmT = new TH1D("h_met_0Lbg1Y1LlmT","h_met_0Lbg1Y1LlmT",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y1LlmT = new TH1D("h_jet1pt_0Lbg1Y1LlmT","h_jet1pt_0Lbg1Y1LlmT",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1LlmT = new TH1D("h_jet2pt_0Lbg1Y1LlmT","h_jet2pt_0Lbg1Y1LlmT",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1LlmT = new TH1D("h_jet3pt_0Lbg1Y1LlmT","h_jet3pt_0Lbg1Y1LlmT",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1LlmT = new TH1D("h_leptonpt_0Lbg1Y1LlmT","h_leptonpt_0Lbg1Y1LlmT",25,0,500);
  TH1D * h_HT_0Lbg1Y1LlmT = new TH1D("h_HT_0Lbg1Y1LlmT", "h_HT_0Lbg1Y1LlmT", 100, 0, 2000);
  TH1D * h_minDeltaPhiHat_0Lbg1Y1LlmT = new TH1D("h_minDeltaPhiHat_0Lbg1Y1LlmT", "h_minDeltaPhiHat_0Lbg1Y1LlmT", 30, 0, 15);
  TH2D * h_MR_minDeltaPhiHat_0Lbg1Y1LlmT = new TH2D("h_MR_minDeltaPhiHat_0Lbg1Y1LlmT", "h_MR_minDeltaPhiHat_0Lbg1Y1LlmT", nbins_MR, bins_MR, 30, 0, 15);
  TH2D * h_R2_minDeltaPhiHat_0Lbg1Y1LlmT = new TH2D("h_R2_minDeltaPhiHat_0Lbg1Y1LlmT", "h_R2_minDeltaPhiHat_0Lbg1Y1LlmT", nbins_R2, bins_R2, 30, 0, 15);
  TH1D * h_minDeltaPhi_0Lbg1Y1LlmT = new TH1D("h_minDeltaPhi_0Lbg1Y1LlmT", "h_minDeltaPhi_0Lbg1Y1LlmT", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_0Lbg1Y1LlmT = new TH2D("h_MR_minDeltaPhi_0Lbg1Y1LlmT", "h_MR_minDeltaPhi_0Lbg1Y1LlmT", nbins_MR, bins_MR, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_0Lbg1Y1LlmT = new TH2D("h_R2_minDeltaPhi_0Lbg1Y1LlmT", "h_R2_minDeltaPhi_0Lbg1Y1LlmT", nbins_R2, bins_R2, 50, 0, 5);

  // mdPhiHatg4 
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhiHatg4", "h_MR_0Lbg1Y1LlmT_mdPhiHatg4", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhiHatg4", "h_R2_0Lbg1Y1LlmT_mdPhiHatg4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4", "h_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4 = new TH2D("h_uw_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4", "h_uw_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_njets_0Lbg1Y1LlmT_mdPhiHatg4","h_njets_0Lbg1Y1LlmT_mdPhiHatg4",15,0,15);
  TH1D * h_met_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_met_0Lbg1Y1LlmT_mdPhiHatg4","h_met_0Lbg1Y1LlmT_mdPhiHatg4",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_jet1pt_0Lbg1Y1LlmT_mdPhiHatg4","h_jet1pt_0Lbg1Y1LlmT_mdPhiHatg4",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_jet2pt_0Lbg1Y1LlmT_mdPhiHatg4","h_jet2pt_0Lbg1Y1LlmT_mdPhiHatg4",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_jet3pt_0Lbg1Y1LlmT_mdPhiHatg4","h_jet3pt_0Lbg1Y1LlmT_mdPhiHatg4",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1LlmT_mdPhiHatg4 = new TH1D("h_leptonpt_0Lbg1Y1LlmT_mdPhiHatg4","h_leptonpt_0Lbg1Y1LlmT_mdPhiHatg4",25,0,500);

  // mdPhiHat4
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhiHat4 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhiHat4", "h_MR_0Lbg1Y1LlmT_mdPhiHat4", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhiHat4 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhiHat4", "h_R2_0Lbg1Y1LlmT_mdPhiHat4", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhiHat4 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhiHat4", "h_MR_R2_0Lbg1Y1LlmT_mdPhiHat4", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // mdPhig0p3 
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhig0p3", "h_MR_0Lbg1Y1LlmT_mdPhig0p3", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhig0p3", "h_R2_0Lbg1Y1LlmT_mdPhig0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhig0p3 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhig0p3", "h_MR_R2_0Lbg1Y1LlmT_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p3 = new TH2D("h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p3", "h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_njets_0Lbg1Y1LlmT_mdPhig0p3","h_njets_0Lbg1Y1LlmT_mdPhig0p3",15,0,15);
  TH1D * h_met_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_met_0Lbg1Y1LlmT_mdPhig0p3","h_met_0Lbg1Y1LlmT_mdPhig0p3",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_jet1pt_0Lbg1Y1LlmT_mdPhig0p3","h_jet1pt_0Lbg1Y1LlmT_mdPhig0p3",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_jet2pt_0Lbg1Y1LlmT_mdPhig0p3","h_jet2pt_0Lbg1Y1LlmT_mdPhig0p3",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_jet3pt_0Lbg1Y1LlmT_mdPhig0p3","h_jet3pt_0Lbg1Y1LlmT_mdPhig0p3",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1LlmT_mdPhig0p3 = new TH1D("h_leptonpt_0Lbg1Y1LlmT_mdPhig0p3","h_leptonpt_0Lbg1Y1LlmT_mdPhig0p3",25,0,500);

  // mdPhi0p3
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhi0p3 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhi0p3", "h_MR_0Lbg1Y1LlmT_mdPhi0p3", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhi0p3 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhi0p3", "h_R2_0Lbg1Y1LlmT_mdPhi0p3", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhi0p3 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhi0p3", "h_MR_R2_0Lbg1Y1LlmT_mdPhi0p3", nbins_MR, bins_MR, nbins_R2, bins_R2);

  // mdPhig0p5 
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhig0p5", "h_MR_0Lbg1Y1LlmT_mdPhig0p5", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhig0p5", "h_R2_0Lbg1Y1LlmT_mdPhig0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhig0p5 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhig0p5", "h_MR_R2_0Lbg1Y1LlmT_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH2D * h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p5 = new TH2D("h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p5", "h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_njets_0Lbg1Y1LlmT_mdPhig0p5","h_njets_0Lbg1Y1LlmT_mdPhig0p5",15,0,15);
  TH1D * h_met_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_met_0Lbg1Y1LlmT_mdPhig0p5","h_met_0Lbg1Y1LlmT_mdPhig0p5",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_jet1pt_0Lbg1Y1LlmT_mdPhig0p5","h_jet1pt_0Lbg1Y1LlmT_mdPhig0p5",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_jet2pt_0Lbg1Y1LlmT_mdPhig0p5","h_jet2pt_0Lbg1Y1LlmT_mdPhig0p5",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_jet3pt_0Lbg1Y1LlmT_mdPhig0p5","h_jet3pt_0Lbg1Y1LlmT_mdPhig0p5",20,0,1000);
  TH1D * h_leptonpt_0Lbg1Y1LlmT_mdPhig0p5 = new TH1D("h_leptonpt_0Lbg1Y1LlmT_mdPhig0p5","h_leptonpt_0Lbg1Y1LlmT_mdPhig0p5",25,0,500);

  // mdPhi0p5
  TH1D * h_MR_0Lbg1Y1LlmT_mdPhi0p5 = new TH1D("h_MR_0Lbg1Y1LlmT_mdPhi0p5", "h_MR_0Lbg1Y1LlmT_mdPhi0p5", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y1LlmT_mdPhi0p5 = new TH1D("h_R2_0Lbg1Y1LlmT_mdPhi0p5", "h_R2_0Lbg1Y1LlmT_mdPhi0p5", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y1LlmT_mdPhi0p5 = new TH2D("h_MR_R2_0Lbg1Y1LlmT_mdPhi0p5", "h_MR_R2_0Lbg1Y1LlmT_mdPhi0p5", nbins_MR, bins_MR, nbins_R2, bins_R2);


  // dimuon trajectory
  TH1D * h_Zmass_2mu = new TH1D("h_Zmass_2mu", "h_Zmass_2mu", 20, 50, 130);
  TH2D * h_R2_Zmass_2mu = new TH2D("h_R2_Zmass_2mu", "h_R2_Zmass_2mu", 20, 50, 130, nbins_MR, bins_MR);
  TH2D * h_MR_Zmass_2mu = new TH2D("h_MR_Zmass_2mu", "h_MR_Zmass_2mu", 20, 50, 130, 20, 0, 1.);

  TH1D * h_MR_2munoZmass = new TH1D("h_MR_2munoZmass", "h_MR_2munoZmass", nbins_MR, bins_MR);
  TH1D * h_R2_2munoZmass = new TH1D("h_R2_2munoZmass", "h_R2_2munoZmass", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2munoZmass = new TH2D("h_MR_R2_2munoZmass", "h_MR_R2_2munoZmass", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_2mu = new TH1D("h_MR_2mu", "h_MR_2mu", nbins_MR, bins_MR);
  TH1D * h_R2_2mu = new TH1D("h_R2_2mu", "h_R2_2mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2mu = new TH2D("h_MR_R2_2mu", "h_MR_R2_2mu", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D * h_lepton1pt_2mu = new TH1D("h_lepton1pt_2mu", "h_lepton1pt_2mu", 25, 0, 500);
  TH1D * h_lepton2pt_2mu = new TH1D("h_lepton2pt_2mu", "h_lepton2pt_2mu", 25, 0, 500);
  
  TH1D * h_MR_2mu0el = new TH1D("h_MR_2mu0el", "h_MR_2mu0el", nbins_MR, bins_MR);
  TH1D * h_R2_2mu0el = new TH1D("h_R2_2mu0el", "h_R2_2mu0el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2mu0el = new TH2D("h_MR_R2_2mu0el", "h_MR_R2_2mu0el", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D * h_PV_2mu0el = new TH1D("h_PV_2mu0el", "h_PV_2mu0el", 50, 0, 50);
  
  TH1D * h_njets_2mu0el = new TH1D("h_njets_2mu0el","h_njets_2mu0el",15,0,15);
  TH1D * h_nbjets_2mu0el = new TH1D("h_nbjets_2mu0el","h_nbjets_2mu0el",6,0,6);
  TH1D * h_met_2mu0el = new TH1D("h_met_2mu0el","h_met_2mu0el",20,0,1000);
  TH1D * h_jet1pt_2mu0el = new TH1D("h_jet1pt_2mu0el","h_jet1pt_2mu0el",20,0,1000);
  TH1D * h_jet2pt_2mu0el = new TH1D("h_jet2pt_2mu0el","h_jet2pt_2mu0el",20,0,1000);
  TH1D * h_jet3pt_2mu0el = new TH1D("h_jet3pt_2mu0el","h_jet3pt_2mu0el",20,0,1000);
  TH1D * h_lepton1pt_2mu0el = new TH1D("h_lepton1pt_2mu0el", "h_lepton1pt_2mu0el", 25, 0, 500);
  TH1D * h_lepton2pt_2mu0el = new TH1D("h_lepton2pt_2mu0el", "h_lepton2pt_2mu0el", 25, 0, 500);
  TH1D * h_HT_2mu0el = new TH1D("h_HT_2mu0el", "h_HT_2mu0el", 100, 0, 2000);

  TH1D * h_MR_0Lb2mu0el = new TH1D("h_MR_0Lb2mu0el", "h_MR_0Lb2mu0el", nbins_MR, bins_MR);
  TH1D * h_R2_0Lb2mu0el = new TH1D("h_R2_0Lb2mu0el", "h_R2_0Lb2mu0el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lb2mu0el = new TH2D("h_MR_R2_0Lb2mu0el", "h_MR_R2_0Lb2mu0el", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lb2mu0el = new TH1D("h_njets_0Lb2mu0el","h_njets_0Lb2mu0el",15,0,15);
  TH1D * h_met_0Lb2mu0el = new TH1D("h_met_0Lb2mu0el","h_met_0Lb2mu0el",20,0,1000);
  TH1D * h_jet1pt_0Lb2mu0el = new TH1D("h_jet1pt_0Lb2mu0el","h_jet1pt_0Lb2mu0el",20,0,1000);
  TH1D * h_jet2pt_0Lb2mu0el = new TH1D("h_jet2pt_0Lb2mu0el","h_jet2pt_0Lb2mu0el",20,0,1000);
  TH1D * h_jet3pt_0Lb2mu0el = new TH1D("h_jet3pt_0Lb2mu0el","h_jet3pt_0Lb2mu0el",20,0,1000);
  TH1D * h_lepton1pt_0Lb2mu0el = new TH1D("h_lepton1pt_0Lb2mu0el", "h_lepton1pt_0Lb2mu0el", 25, 0, 500);
  TH1D * h_lepton2pt_0Lb2mu0el = new TH1D("h_lepton2pt_0Lb2mu0el", "h_lepton2pt_0Lb2mu0el", 25, 0, 500);

  TH1D * h_MR_g1Mb2mu0el = new TH1D("h_MR_g1Mb2mu0el", "h_MR_g1Mb2mu0el", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb2mu0el = new TH1D("h_R2_g1Mb2mu0el", "h_R2_g1Mb2mu0el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb2mu0el = new TH2D("h_MR_R2_g1Mb2mu0el", "h_MR_R2_g1Mb2mu0el", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mb2mu0el = new TH1D("h_njets_g1Mb2mu0el","h_njets_g1Mb2mu0el",15,0,15);
  TH1D * h_nbjets_g1Mb2mu0el = new TH1D("h_nbjets_g1Mb2mu0el","h_nbjets_g1Mb2mu0el",6,0,6);
  TH1D * h_met_g1Mb2mu0el = new TH1D("h_met_g1Mb2mu0el","h_met_g1Mb2mu0el",20,0,1000);
  TH1D * h_jet1pt_g1Mb2mu0el = new TH1D("h_jet1pt_g1Mb2mu0el","h_jet1pt_g1Mb2mu0el",20,0,1000);
  TH1D * h_jet2pt_g1Mb2mu0el = new TH1D("h_jet2pt_g1Mb2mu0el","h_jet2pt_g1Mb2mu0el",20,0,1000);
  TH1D * h_jet3pt_g1Mb2mu0el = new TH1D("h_jet3pt_g1Mb2mu0el","h_jet3pt_g1Mb2mu0el",20,0,1000);
  TH1D * h_lepton1pt_g1Mb2mu0el = new TH1D("h_lepton1pt_g1Mb2mu0el", "h_lepton1pt_g1Mb2mu0el", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mb2mu0el = new TH1D("h_lepton2pt_g1Mb2mu0el", "h_lepton2pt_g1Mb2mu0el", 25, 0, 500);

  // Z no b Control region  
  TH1D * h_MR_0Lbg1Y2mu0el = new TH1D("h_MR_0Lbg1Y2mu0el", "h_MR_0Lbg1Y2mu0el", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y2mu0el = new TH1D("h_R2_0Lbg1Y2mu0el", "h_R2_0Lbg1Y2mu0el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y2mu0el = new TH2D("h_MR_R2_0Lbg1Y2mu0el", "h_MR_R2_0Lbg1Y2mu0el", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lbg1Y2mu0el = new TH1D("h_njets_0Lbg1Y2mu0el","h_njets_0Lbg1Y2mu0el",15,0,15);
  TH1D * h_met_0Lbg1Y2mu0el = new TH1D("h_met_0Lbg1Y2mu0el","h_met_0Lbg1Y2mu0el",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y2mu0el = new TH1D("h_jet1pt_0Lbg1Y2mu0el","h_jet1pt_0Lbg1Y2mu0el",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y2mu0el = new TH1D("h_jet2pt_0Lbg1Y2mu0el","h_jet2pt_0Lbg1Y2mu0el",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y2mu0el = new TH1D("h_jet3pt_0Lbg1Y2mu0el","h_jet3pt_0Lbg1Y2mu0el",20,0,1000);
  TH1D * h_lepton1pt_0Lbg1Y2mu0el = new TH1D("h_lepton1pt_0Lbg1Y2mu0el", "h_lepton1pt_0Lbg1Y2mu0el", 25, 0, 500);
  TH1D * h_lepton2pt_0Lbg1Y2mu0el = new TH1D("h_lepton2pt_0Lbg1Y2mu0el", "h_lepton2pt_0Lbg1Y2mu0el", 25, 0, 500);

  // Z with b Control region
  TH1D * h_MR_g1Mbg1Y2mu0el = new TH1D("h_MR_g1Mbg1Y2mu0el", "h_MR_g1Mbg1Y2mu0el", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1Y2mu0el = new TH1D("h_R2_g1Mbg1Y2mu0el", "h_R2_g1Mbg1Y2mu0el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1Y2mu0el = new TH2D("h_MR_R2_g1Mbg1Y2mu0el", "h_MR_R2_g1Mbg1Y2mu0el", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1Y2mu0el = new TH1D("h_njets_g1Mbg1Y2mu0el","h_njets_g1Mbg1Y2mu0el",15,0,15);
  TH1D * h_nbjets_g1Mbg1Y2mu0el = new TH1D("h_nbjets_g1Mbg1Y2mu0el","h_nbjets_g1Mbg1Y2mu0el",6,0,6);
  TH1D * h_met_g1Mbg1Y2mu0el = new TH1D("h_met_g1Mbg1Y2mu0el","h_met_g1Mbg1Y2mu0el",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1Y2mu0el = new TH1D("h_jet1pt_g1Mbg1Y2mu0el","h_jet1pt_g1Mbg1Y2mu0el",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1Y2mu0el = new TH1D("h_jet2pt_g1Mbg1Y2mu0el","h_jet2pt_g1Mbg1Y2mu0el",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1Y2mu0el = new TH1D("h_jet3pt_g1Mbg1Y2mu0el","h_jet3pt_g1Mbg1Y2mu0el",20,0,1000);
  TH1D * h_lepton1pt_g1Mbg1Y2mu0el = new TH1D("h_lepton1pt_g1Mbg1Y2mu0el", "h_lepton1pt_g1Mbg1Y2mu0el", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mbg1Y2mu0el = new TH1D("h_lepton2pt_g1Mbg1Y2mu0el", "h_lepton2pt_g1Mbg1Y2mu0el", 25, 0, 500);


  // dielectron trajectory
  TH1D * h_Zmass_2el = new TH1D("h_Zmass_2el", "h_Zmass_2el", 20, 50, 130);
  TH2D * h_R2_Zmass_2el = new TH2D("h_R2_Zmass_2el", "h_R2_Zmass_2el", 20, 50, 130, nbins_MR, bins_MR);
  TH2D * h_MR_Zmass_2el = new TH2D("h_MR_Zmass_2el", "h_MR_Zmass_2el", 20, 50, 130, 20, 0, 1.);

  TH1D * h_MR_2elnoZmass = new TH1D("h_MR_2elnoZmass", "h_MR_2elnoZmass", nbins_MR, bins_MR);
  TH1D * h_R2_2elnoZmass = new TH1D("h_R2_2elnoZmass", "h_R2_2elnoZmass", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2elnoZmass = new TH2D("h_MR_R2_2elnoZmass", "h_MR_R2_2elnoZmass", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_2el = new TH1D("h_MR_2el", "h_MR_2el", nbins_MR, bins_MR);
  TH1D * h_R2_2el = new TH1D("h_R2_2el", "h_R2_2el", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2el = new TH2D("h_MR_R2_2el", "h_MR_R2_2el", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D * h_lepton1pt_2el = new TH1D("h_lepton1pt_2el", "h_lepton1pt_2el", 25, 0, 500);
  TH1D * h_lepton2pt_2el = new TH1D("h_lepton2pt_2el", "h_lepton2pt_2el", 25, 0, 500);
  
  TH1D * h_MR_2el0mu = new TH1D("h_MR_2el0mu", "h_MR_2el0mu", nbins_MR, bins_MR);
  TH1D * h_R2_2el0mu = new TH1D("h_R2_2el0mu", "h_R2_2el0mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2el0mu = new TH2D("h_MR_R2_2el0mu", "h_MR_R2_2el0mu", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D * h_PV_2el0mu = new TH1D("h_PV_2el0mu", "h_PV_2el0mu", 50, 0, 50);
  
  TH1D * h_njets_2el0mu = new TH1D("h_njets_2el0mu","h_njets_2el0mu",15,0,15);
  TH1D * h_nbjets_2el0mu = new TH1D("h_nbjets_2el0mu","h_nbjets_2el0mu",6,0,6);
  TH1D * h_met_2el0mu = new TH1D("h_met_2el0mu","h_met_2el0mu",20,0,1000);
  TH1D * h_jet1pt_2el0mu = new TH1D("h_jet1pt_2el0mu","h_jet1pt_2el0mu",20,0,1000);
  TH1D * h_jet2pt_2el0mu = new TH1D("h_jet2pt_2el0mu","h_jet2pt_2el0mu",20,0,1000);
  TH1D * h_jet3pt_2el0mu = new TH1D("h_jet3pt_2el0mu","h_jet3pt_2el0mu",20,0,1000);
  TH1D * h_lepton1pt_2el0mu = new TH1D("h_lepton1pt_2el0mu", "h_lepton1pt_2el0mu", 25, 0, 500);
  TH1D * h_lepton2pt_2el0mu = new TH1D("h_lepton2pt_2el0mu", "h_lepton2pt_2el0mu", 25, 0, 500);
  TH1D * h_HT_2el0mu = new TH1D("h_HT_2el0mu", "h_HT_2el0mu", 100, 0, 2000);

  TH1D * h_MR_0Lb2el0mu = new TH1D("h_MR_0Lb2el0mu", "h_MR_0Lb2el0mu", nbins_MR, bins_MR);
  TH1D * h_R2_0Lb2el0mu = new TH1D("h_R2_0Lb2el0mu", "h_R2_0Lb2el0mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lb2el0mu = new TH2D("h_MR_R2_0Lb2el0mu", "h_MR_R2_0Lb2el0mu", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lb2el0mu = new TH1D("h_njets_0Lb2el0mu","h_njets_0Lb2el0mu",15,0,15);
  TH1D * h_met_0Lb2el0mu = new TH1D("h_met_0Lb2el0mu","h_met_0Lb2el0mu",20,0,1000);
  TH1D * h_jet1pt_0Lb2el0mu = new TH1D("h_jet1pt_0Lb2el0mu","h_jet1pt_0Lb2el0mu",20,0,1000);
  TH1D * h_jet2pt_0Lb2el0mu = new TH1D("h_jet2pt_0Lb2el0mu","h_jet2pt_0Lb2el0mu",20,0,1000);
  TH1D * h_jet3pt_0Lb2el0mu = new TH1D("h_jet3pt_0Lb2el0mu","h_jet3pt_0Lb2el0mu",20,0,1000);
  TH1D * h_lepton1pt_0Lb2el0mu = new TH1D("h_lepton1pt_0Lb2el0mu", "h_lepton1pt_0Lb2el0mu", 25, 0, 500);
  TH1D * h_lepton2pt_0Lb2el0mu = new TH1D("h_lepton2pt_0Lb2el0mu", "h_lepton2pt_0Lb2el0mu", 25, 0, 500);

  TH1D * h_MR_g1Mb2el0mu = new TH1D("h_MR_g1Mb2el0mu", "h_MR_g1Mb2el0mu", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb2el0mu = new TH1D("h_R2_g1Mb2el0mu", "h_R2_g1Mb2el0mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb2el0mu = new TH2D("h_MR_R2_g1Mb2el0mu", "h_MR_R2_g1Mb2el0mu", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mb2el0mu = new TH1D("h_njets_g1Mb2el0mu","h_njets_g1Mb2el0mu",15,0,15);
  TH1D * h_nbjets_g1Mb2el0mu = new TH1D("h_nbjets_g1Mb2el0mu","h_nbjets_g1Mb2el0mu",6,0,6);
  TH1D * h_met_g1Mb2el0mu = new TH1D("h_met_g1Mb2el0mu","h_met_g1Mb2el0mu",20,0,1000);
  TH1D * h_jet1pt_g1Mb2el0mu = new TH1D("h_jet1pt_g1Mb2el0mu","h_jet1pt_g1Mb2el0mu",20,0,1000);
  TH1D * h_jet2pt_g1Mb2el0mu = new TH1D("h_jet2pt_g1Mb2el0mu","h_jet2pt_g1Mb2el0mu",20,0,1000);
  TH1D * h_jet3pt_g1Mb2el0mu = new TH1D("h_jet3pt_g1Mb2el0mu","h_jet3pt_g1Mb2el0mu",20,0,1000);
  TH1D * h_lepton1pt_g1Mb2el0mu = new TH1D("h_lepton1pt_g1Mb2el0mu", "h_lepton1pt_g1Mb2el0mu", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mb2el0mu = new TH1D("h_lepton2pt_g1Mb2el0mu", "h_lepton2pt_g1Mb2el0mu", 25, 0, 500);

  // Z no b Control region  
  TH1D * h_MR_0Lbg1Y2el0mu = new TH1D("h_MR_0Lbg1Y2el0mu", "h_MR_0Lbg1Y2el0mu", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y2el0mu = new TH1D("h_R2_0Lbg1Y2el0mu", "h_R2_0Lbg1Y2el0mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y2el0mu = new TH2D("h_MR_R2_0Lbg1Y2el0mu", "h_MR_R2_0Lbg1Y2el0mu", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lbg1Y2el0mu = new TH1D("h_njets_0Lbg1Y2el0mu","h_njets_0Lbg1Y2el0mu",15,0,15);
  TH1D * h_met_0Lbg1Y2el0mu = new TH1D("h_met_0Lbg1Y2el0mu","h_met_0Lbg1Y2el0mu",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y2el0mu = new TH1D("h_jet1pt_0Lbg1Y2el0mu","h_jet1pt_0Lbg1Y2el0mu",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y2el0mu = new TH1D("h_jet2pt_0Lbg1Y2el0mu","h_jet2pt_0Lbg1Y2el0mu",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y2el0mu = new TH1D("h_jet3pt_0Lbg1Y2el0mu","h_jet3pt_0Lbg1Y2el0mu",20,0,1000);
  TH1D * h_lepton1pt_0Lbg1Y2el0mu = new TH1D("h_lepton1pt_0Lbg1Y2el0mu", "h_lepton1pt_0Lbg1Y2el0mu", 25, 0, 500);
  TH1D * h_lepton2pt_0Lbg1Y2el0mu = new TH1D("h_lepton2pt_0Lbg1Y2el0mu", "h_lepton2pt_0Lbg1Y2el0mu", 25, 0, 500);

  // Z with b Control Region
  TH1D * h_MR_g1Mbg1Y2el0mu = new TH1D("h_MR_g1Mbg1Y2el0mu", "h_MR_g1Mbg1Y2el0mu", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1Y2el0mu = new TH1D("h_R2_g1Mbg1Y2el0mu", "h_R2_g1Mbg1Y2el0mu", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1Y2el0mu = new TH2D("h_MR_R2_g1Mbg1Y2el0mu", "h_MR_R2_g1Mbg1Y2el0mu", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1Y2el0mu = new TH1D("h_njets_g1Mbg1Y2el0mu","h_njets_g1Mbg1Y2el0mu",15,0,15);
  TH1D * h_nbjets_g1Mbg1Y2el0mu = new TH1D("h_nbjets_g1Mbg1Y2el0mu","h_nbjets_g1Mbg1Y2el0mu",6,0,6);
  TH1D * h_met_g1Mbg1Y2el0mu = new TH1D("h_met_g1Mbg1Y2el0mu","h_met_g1Mbg1Y2el0mu",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1Y2el0mu = new TH1D("h_jet1pt_g1Mbg1Y2el0mu","h_jet1pt_g1Mbg1Y2el0mu",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1Y2el0mu = new TH1D("h_jet2pt_g1Mbg1Y2el0mu","h_jet2pt_g1Mbg1Y2el0mu",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1Y2el0mu = new TH1D("h_jet3pt_g1Mbg1Y2el0mu","h_jet3pt_g1Mbg1Y2el0mu",20,0,1000);
  TH1D * h_lepton1pt_g1Mbg1Y2el0mu = new TH1D("h_lepton1pt_g1Mbg1Y2el0mu", "h_lepton1pt_g1Mbg1Y2el0mu", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mbg1Y2el0mu = new TH1D("h_lepton2pt_g1Mbg1Y2el0mu", "h_lepton2pt_g1Mbg1Y2el0mu", 25, 0, 500);

  // dilepton trajectory
  TH1D * h_Zmass_2l = new TH1D("h_Zmass_2l", "h_Zmass_2l", 20, 50, 130);
  TH2D * h_R2_Zmass_2l = new TH2D("h_R2_Zmass_2l", "h_R2_Zmass_2l", 20, 50, 130, nbins_MR, bins_MR);
  TH2D * h_MR_Zmass_2l = new TH2D("h_MR_Zmass_2l", "h_MR_Zmass_2l", 20, 50, 130, 20, 0, 1.);

  TH1D * h_MR_2lnoZmass = new TH1D("h_MR_2lnoZmass", "h_MR_2lnoZmass", nbins_MR, bins_MR);
  TH1D * h_R2_2lnoZmass = new TH1D("h_R2_2lnoZmass", "h_R2_2lnoZmass", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2lnoZmass = new TH2D("h_MR_R2_2lnoZmass", "h_MR_R2_2lnoZmass", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_MR_2l = new TH1D("h_MR_2l", "h_MR_2l", nbins_MR, bins_MR);
  TH1D * h_R2_2l = new TH1D("h_R2_2l", "h_R2_2l", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2l = new TH2D("h_MR_R2_2l", "h_MR_R2_2l", nbins_MR, bins_MR, nbins_R2, bins_R2);
  TH1D * h_lepton1pt_2l = new TH1D("h_lepton1pt_2l", "h_lepton1pt_2l", 25, 0, 500);
  TH1D * h_lepton2pt_2l = new TH1D("h_lepton2pt_2l", "h_lepton2pt_2l", 25, 0, 500);
  
  TH1D * h_MR_2l0ol = new TH1D("h_MR_2l0ol", "h_MR_2l0ol", nbins_MR, bins_MR);
  TH1D * h_R2_2l0ol = new TH1D("h_R2_2l0ol", "h_R2_2l0ol", nbins_R2, bins_R2);
  TH2D * h_MR_R2_2l0ol = new TH2D("h_MR_R2_2l0ol", "h_MR_R2_2l0ol", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_2l0ol = new TH1D("h_njets_2l0ol","h_njets_2l0ol",15,0,15);
  TH1D * h_nbjets_2l0ol = new TH1D("h_nbjets_2l0ol","h_nbjets_2l0ol",6,0,6);
  TH1D * h_met_2l0ol = new TH1D("h_met_2l0ol","h_met_2l0ol",20,0,1000);
  TH1D * h_jet1pt_2l0ol = new TH1D("h_jet1pt_2l0ol","h_jet1pt_2l0ol",20,0,1000);
  TH1D * h_jet2pt_2l0ol = new TH1D("h_jet2pt_2l0ol","h_jet2pt_2l0ol",20,0,1000);
  TH1D * h_jet3pt_2l0ol = new TH1D("h_jet3pt_2l0ol","h_jet3pt_2l0ol",20,0,1000);
  TH1D * h_lepton1pt_2l0ol = new TH1D("h_lepton1pt_2l0ol", "h_lepton1pt_2l0ol", 25, 0, 500);
  TH1D * h_lepton2pt_2l0ol = new TH1D("h_lepton2pt_2l0ol", "h_lepton2pt_2l0ol", 25, 0, 500);
  TH1D * h_HT_2l0ol = new TH1D("h_HT_2l0ol", "h_HT_2l0ol", 100, 0, 2000);

  TH1D * h_MR_0Lb2l0ol = new TH1D("h_MR_0Lb2l0ol", "h_MR_0Lb2l0ol", nbins_MR, bins_MR);
  TH1D * h_R2_0Lb2l0ol = new TH1D("h_R2_0Lb2l0ol", "h_R2_0Lb2l0ol", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lb2l0ol = new TH2D("h_MR_R2_0Lb2l0ol", "h_MR_R2_0Lb2l0ol", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lb2l0ol = new TH1D("h_njets_0Lb2l0ol","h_njets_0Lb2l0ol",15,0,15);
  TH1D * h_met_0Lb2l0ol = new TH1D("h_met_0Lb2l0ol","h_met_0Lb2l0ol",20,0,1000);
  TH1D * h_jet1pt_0Lb2l0ol = new TH1D("h_jet1pt_0Lb2l0ol","h_jet1pt_0Lb2l0ol",20,0,1000);
  TH1D * h_jet2pt_0Lb2l0ol = new TH1D("h_jet2pt_0Lb2l0ol","h_jet2pt_0Lb2l0ol",20,0,1000);
  TH1D * h_jet3pt_0Lb2l0ol = new TH1D("h_jet3pt_0Lb2l0ol","h_jet3pt_0Lb2l0ol",20,0,1000);
  TH1D * h_lepton1pt_0Lb2l0ol = new TH1D("h_lepton1pt_0Lb2l0ol", "h_lepton1pt_0Lb2l0ol", 25, 0, 500);
  TH1D * h_lepton2pt_0Lb2l0ol = new TH1D("h_lepton2pt_0Lb2l0ol", "h_lepton2pt_0Lb2l0ol", 25, 0, 500);

  TH1D * h_MR_g1Mb2l0ol = new TH1D("h_MR_g1Mb2l0ol", "h_MR_g1Mb2l0ol", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mb2l0ol = new TH1D("h_R2_g1Mb2l0ol", "h_R2_g1Mb2l0ol", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mb2l0ol = new TH2D("h_MR_R2_g1Mb2l0ol", "h_MR_R2_g1Mb2l0ol", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mb2l0ol = new TH1D("h_njets_g1Mb2l0ol","h_njets_g1Mb2l0ol",15,0,15);
  TH1D * h_nbjets_g1Mb2l0ol = new TH1D("h_nbjets_g1Mb2l0ol","h_nbjets_g1Mb2l0ol",6,0,6);
  TH1D * h_met_g1Mb2l0ol = new TH1D("h_met_g1Mb2l0ol","h_met_g1Mb2l0ol",20,0,1000);
  TH1D * h_jet1pt_g1Mb2l0ol = new TH1D("h_jet1pt_g1Mb2l0ol","h_jet1pt_g1Mb2l0ol",20,0,1000);
  TH1D * h_jet2pt_g1Mb2l0ol = new TH1D("h_jet2pt_g1Mb2l0ol","h_jet2pt_g1Mb2l0ol",20,0,1000);
  TH1D * h_jet3pt_g1Mb2l0ol = new TH1D("h_jet3pt_g1Mb2l0ol","h_jet3pt_g1Mb2l0ol",20,0,1000);
  TH1D * h_lepton1pt_g1Mb2l0ol = new TH1D("h_lepton1pt_g1Mb2l0ol", "h_lepton1pt_g1Mb2l0ol", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mb2l0ol = new TH1D("h_lepton2pt_g1Mb2l0ol", "h_lepton2pt_g1Mb2l0ol", 25, 0, 500);

  // Z no b Control Region  
  TH1D * h_MR_0Lbg1Y2l0ol = new TH1D("h_MR_0Lbg1Y2l0ol", "h_MR_0Lbg1Y2l0ol", nbins_MR, bins_MR);
  TH1D * h_R2_0Lbg1Y2l0ol = new TH1D("h_R2_0Lbg1Y2l0ol", "h_R2_0Lbg1Y2l0ol", nbins_R2, bins_R2);
  TH2D * h_MR_R2_0Lbg1Y2l0ol = new TH2D("h_MR_R2_0Lbg1Y2l0ol", "h_MR_R2_0Lbg1Y2l0ol", nbins_MR, bins_MR, nbins_R2, bins_R2);
  
  TH1D * h_njets_0Lbg1Y2l0ol = new TH1D("h_njets_0Lbg1Y2l0ol","h_njets_0Lbg1Y2l0ol",15,0,15);
  TH1D * h_met_0Lbg1Y2l0ol = new TH1D("h_met_0Lbg1Y2l0ol","h_met_0Lbg1Y2l0ol",20,0,1000);
  TH1D * h_jet1pt_0Lbg1Y2l0ol = new TH1D("h_jet1pt_0Lbg1Y2l0ol","h_jet1pt_0Lbg1Y2l0ol",20,0,1000);
  TH1D * h_jet2pt_0Lbg1Y2l0ol = new TH1D("h_jet2pt_0Lbg1Y2l0ol","h_jet2pt_0Lbg1Y2l0ol",20,0,1000);
  TH1D * h_jet3pt_0Lbg1Y2l0ol = new TH1D("h_jet3pt_0Lbg1Y2l0ol","h_jet3pt_0Lbg1Y2l0ol",20,0,1000);
  TH1D * h_lepton1pt_0Lbg1Y2l0ol = new TH1D("h_lepton1pt_0Lbg1Y2l0ol", "h_lepton1pt_0Lbg1Y2l0ol", 25, 0, 500);
  TH1D * h_lepton2pt_0Lbg1Y2l0ol = new TH1D("h_lepton2pt_0Lbg1Y2l0ol", "h_lepton2pt_0Lbg1Y2l0ol", 25, 0, 500);

  // Z with b Control region
  TH1D * h_MR_g1Mbg1Y2l0ol = new TH1D("h_MR_g1Mbg1Y2l0ol", "h_MR_g1Mbg1Y2l0ol", nbins_MR, bins_MR);
  TH1D * h_R2_g1Mbg1Y2l0ol = new TH1D("h_R2_g1Mbg1Y2l0ol", "h_R2_g1Mbg1Y2l0ol", nbins_R2, bins_R2);
  TH2D * h_MR_R2_g1Mbg1Y2l0ol = new TH2D("h_MR_R2_g1Mbg1Y2l0ol", "h_MR_R2_g1Mbg1Y2l0ol", nbins_MR, bins_MR, nbins_R2, bins_R2);

  TH1D * h_njets_g1Mbg1Y2l0ol = new TH1D("h_njets_g1Mbg1Y2l0ol","h_njets_g1Mbg1Y2l0ol",15,0,15);
  TH1D * h_nbjets_g1Mbg1Y2l0ol = new TH1D("h_nbjets_g1Mbg1Y2l0ol","h_nbjets_g1Mbg1Y2l0ol",6,0,6);
  TH1D * h_met_g1Mbg1Y2l0ol = new TH1D("h_met_g1Mbg1Y2l0ol","h_met_g1Mbg1Y2l0ol",20,0,1000);
  TH1D * h_jet1pt_g1Mbg1Y2l0ol = new TH1D("h_jet1pt_g1Mbg1Y2l0ol","h_jet1pt_g1Mbg1Y2l0ol",20,0,1000);
  TH1D * h_jet2pt_g1Mbg1Y2l0ol = new TH1D("h_jet2pt_g1Mbg1Y2l0ol","h_jet2pt_g1Mbg1Y2l0ol",20,0,1000);
  TH1D * h_jet3pt_g1Mbg1Y2l0ol = new TH1D("h_jet3pt_g1Mbg1Y2l0ol","h_jet3pt_g1Mbg1Y2l0ol",20,0,1000);
  TH1D * h_lepton1pt_g1Mbg1Y2l0ol = new TH1D("h_lepton1pt_g1Mbg1Y2l0ol", "h_lepton1pt_g1Mbg1Y2l0ol", 25, 0, 500);
  TH1D * h_lepton2pt_g1Mbg1Y2l0ol = new TH1D("h_lepton2pt_g1Mbg1Y2l0ol", "h_lepton2pt_g1Mbg1Y2l0ol", 25, 0, 500);

  // Gen level plots
  TH1D* h_gen_toppt = new TH1D("h_gen_toppt", "h_gen_toppt", 50, 0, 1000);
  TH1D* h_gen_dRWb = new TH1D("h_gen_dRWb", "h_gen_dRWb", 200, 0, 5);
  TH2D* h_gen_toppt_dRWb = new TH2D("h_gen_toppt_dRWb", "h_gen_toppt_dRWb", 200, 0, 1000, 200, 0, 5);
  TH1D* h_gen_Wpt = new TH1D("h_gen_Wpt", "h_gen_Wpt", 50, 0, 1000);
  TH1D* h_gen_dRqq = new TH1D("h_gen_dRqq", "h_gen_dRqq", 200, 0, 5);
  TH2D* h_gen_Wpt_dRqq = new TH2D("h_gen_Wpt_dRqq", "h_gen_Wpt_dRqq", 200, 0, 1000, 200, 0, 5);

  TH1D* h_gen_top1pt_g1Mb0Ll = new TH1D("h_gen_top1pt_g1Mb0Ll", "h_gen_top1pt_g1Mb0Ll", 50, 0, 1000);
  TH1D* h_gen_dRWb_g1Mb0Ll = new TH1D("h_gen_dRWb_g1Mb0Ll", "h_gen_dRWb_g1Mb0Ll", 200, 0, 5);
  TH2D* h_gen_top1pt_dRWb_g1Mb0Ll = new TH2D("h_gen_top1pt_dRWb_g1Mb0Ll", "h_gen_top1pt_dRWb_g1Mb0Ll", 200, 0, 1000, 200, 0, 5);
  TH1D* h_gen_W1pt_g1Mb0Ll = new TH1D("h_gen_W1pt_g1Mb0Ll", "h_gen_W1pt_g1Mb0Ll", 50, 0, 1000);
  TH1D* h_gen_dRqq_g1Mb0Ll = new TH1D("h_gen_dRqq_g1Mb0Ll", "h_gen_dRqq_g1Mb0Ll", 200, 0, 5);
  TH2D* h_gen_W1pt_dRqq_g1Mb0Ll = new TH2D("h_gen_W1pt_dRqq_g1Mb0Ll", "h_gen_W1pt_dRqq_g1Mb0Ll", 200, 0, 1000, 200, 0, 5);

  // ------------------------------------------------------
  // -- Define the order of bins in the counts histogram --
  // ------------------------------------------------------

  ofile.count("NoCuts2genallmu", 0.0);
  ofile.count("NoCuts2genmu", 0.0);
  ofile.count("NoCuts2mu", 0.0);
  
  ofile.count("NoCuts", 0.0);
  ofile.count("Cleaning", 0.0);
  ofile.count("Pileup", 0.0);
  ofile.count("ISR", 0.0);
  ofile.count("TopPt", 0.0);
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
  ofile.count("1Mbg1W0Ll", 0.0);
  ofile.count("g2Mbg1W0Ll", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhiHat4", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhi0p3", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhig0p3", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhi0p5", 0.0);
  ofile.count("g1Mbg1W0Ll_mdPhig0p5", 0.0);
  ofile.count("g1Mb0Wg1uW0Ll", 0.0);
  ofile.count("0Lb0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhig0p3", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhi0p5", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhig0p5", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhiHatg4", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhiHat5", 0.0);
  ofile.count("0Lbg1W0Ll", 0.0);

  ofile.count("1Ll", 0.0);
  ofile.count("g1Mb1Ll", 0.0);
  ofile.count("g1Mbg1W1Ll", 0.0);
  ofile.count("g1Mbg1W1LlmT100", 0.0);
  ofile.count("1Mbg1W1LlmT100", 0.0);
  ofile.count("g2Mbg1W1LlmT100", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhiHat4", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhiHatg4", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhi0p3", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhig0p3", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhi0p5", 0.0);
  ofile.count("g1Mbg1W1LlmT100_mdPhig0p5", 0.0);
  ofile.count("g1Mbg1W1LlmT", 0.0);
  ofile.count("0Lb1Ll", 0.0);
  ofile.count("0Lbg1Y1Ll", 0.0);
  ofile.count("0Lbg1Y1LlmT100", 0.0);
  ofile.count("0Lbg1Y1LlmT", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhiHat4", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhiHatg4", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhi0p3", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhig0p3", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhi0p5", 0.0);
  ofile.count("0Lbg1Y1LlmT_mdPhig0p5", 0.0);

  ofile.count("presel2genallmu", 0.0);
  ofile.count("presel2genmu", 0.0);
  ofile.count("presel2mu", 0.0);
  ofile.count("2genallmu", 0.0);
  ofile.count("2genmu", 0.0);
  ofile.count("2munoZmass", 0.0);
  ofile.count("2mu", 0.0);
  ofile.count("2mu0el", 0.0);
  ofile.count("0Lb2mu0el", 0.0);
  ofile.count("g1Mb2mu0el", 0.0);
  ofile.count("0Lbg1Y2mu0el", 0.0);
  ofile.count("g1Mbg1Y2mu0el", 0.0);
  
  ofile.count("2elnoZmass", 0.0);
  ofile.count("2el", 0.0);
  ofile.count("2el0mu", 0.0);
  ofile.count("0Lb2el0mu", 0.0);
  ofile.count("g1Mb2el0mu", 0.0);
  ofile.count("0Lbg1Y2el0mu", 0.0);
  ofile.count("g1Mbg1Y2el0mu", 0.0);
  
  ofile.count("2lnoZmass", 0.0);
  ofile.count("2l", 0.0);
  ofile.count("2l0ol", 0.0);
  ofile.count("0Lb2l0ol", 0.0);
  ofile.count("g1Mb2l0ol", 0.0);
  ofile.count("0Lbg1Y2l0ol", 0.0);
  ofile.count("g1Mbg1Y2l0ol", 0.0);

  // ------------------------------------------------
  // -- Define order of cuts for ttbar composition --
  // ------------------------------------------------
  
  TH1D* TTallhad = new TH1D("counts_TTallhad","",1,0,1);
  TTallhad->SetBit(TH1::kCanRebin);
  TH1D* TTsemilep = new TH1D("counts_TTsemilep","",1,0,1);
  TTsemilep->SetBit(TH1::kCanRebin);
  TH1D* TTdilep = new TH1D("counts_TTdilep","",1,0,1);
  TTdilep->SetBit(TH1::kCanRebin);

  TTallhad->Fill("NoCuts2genallmu", 0.0);
  TTallhad->Fill("NoCuts2genmu", 0.0);
  TTallhad->Fill("NoCuts2mu", 0.0);
  TTallhad->Fill("NoCuts", 0.0);
  TTallhad->Fill("Cleaning", 0.0);
  TTallhad->Fill("Pileup", 0.0);
  TTallhad->Fill("ISR", 0.0);
  TTallhad->Fill("TopPt", 0.0);
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
  TTallhad->Fill("1Mbg1W0Ll", 0.0);
  TTallhad->Fill("g2Mbg1W0Ll", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhiHat4", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhi0p3", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhig0p3", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhi0p5", 0.0);
  TTallhad->Fill("g1Mbg1W0Ll_mdPhig0p5", 0.0);
  TTallhad->Fill("g1Mb0Wg1uW0Ll", 0.0);
  TTallhad->Fill("0Lb0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhig0p3", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p5", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhig0p5", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhiHatg4", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat5", 0.0);
  TTallhad->Fill("0Lbg1W0Ll", 0.0);
  TTallhad->Fill("1Ll", 0.0);
  TTallhad->Fill("g1Mb1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100", 0.0);
  TTallhad->Fill("1Mbg1W1LlmT100", 0.0);
  TTallhad->Fill("g2Mbg1W1LlmT100", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhiHat4", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhi0p3", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhig0p3", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhi0p5", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT100_mdPhig0p5", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT", 0.0);
  TTallhad->Fill("0Lb1Ll", 0.0);
  TTallhad->Fill("0Lbg1Y1Ll", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT100", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhiHat4", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhiHatg4", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhi0p3", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhig0p3", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhi0p5", 0.0);
  TTallhad->Fill("0Lbg1Y1LlmT_mdPhig0p5", 0.0);
  TTallhad->Fill("presel2genallmu", 0.0);
  TTallhad->Fill("presel2genmu", 0.0);
  TTallhad->Fill("presel2mu", 0.0);
  TTallhad->Fill("2genallmu", 0.0);
  TTallhad->Fill("2genmu", 0.0);
  TTallhad->Fill("2munoZmass", 0.0);
  TTallhad->Fill("2mu", 0.0);
  TTallhad->Fill("2mu0el", 0.0);
  TTallhad->Fill("0Lb2mu0el", 0.0);
  TTallhad->Fill("g1Mb2mu0el", 0.0);
  TTallhad->Fill("0Lbg1Y2mu0el", 0.0);
  TTallhad->Fill("g1Mbg1Y2mu0el", 0.0);
  TTallhad->Fill("2elnoZmass", 0.0);
  TTallhad->Fill("2el", 0.0);
  TTallhad->Fill("2el0mu", 0.0);
  TTallhad->Fill("0Lb2el0mu", 0.0);
  TTallhad->Fill("g1Mb2el0mu", 0.0);
  TTallhad->Fill("0Lbg1Y2el0mu", 0.0);
  TTallhad->Fill("g1Mbg1Y2el0mu", 0.0);
  TTallhad->Fill("2lnoZmass", 0.0);
  TTallhad->Fill("2l", 0.0);
  TTallhad->Fill("2l0ol", 0.0);
  TTallhad->Fill("0Lb2l0ol", 0.0);
  TTallhad->Fill("g1Mb2l0ol", 0.0);
  TTallhad->Fill("0Lbg1Y2l0ol", 0.0);
  TTallhad->Fill("g1Mbg1Y2l0ol", 0.0);



  TTsemilep->Fill("NoCuts2genallmu", 0.0);
  TTsemilep->Fill("NoCuts2genmu", 0.0);
  TTsemilep->Fill("NoCuts2mu", 0.0);
  TTsemilep->Fill("NoCuts", 0.0);
  TTsemilep->Fill("Cleaning", 0.0);
  TTsemilep->Fill("Pileup", 0.0);
  TTsemilep->Fill("ISR", 0.0);
  TTsemilep->Fill("TopPt", 0.0);
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
  TTsemilep->Fill("1Mbg1W0Ll", 0.0);
  TTsemilep->Fill("g2Mbg1W0Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHat4", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhi0p3", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhig0p3", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhi0p5", 0.0);
  TTsemilep->Fill("g1Mbg1W0Ll_mdPhig0p5", 0.0);
  TTsemilep->Fill("g1Mb0Wg1uW0Ll", 0.0);
  TTsemilep->Fill("0Lb0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhig0p3", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p5", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhig0p5", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHatg4", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat5", 0.0);
  TTsemilep->Fill("0Lbg1W0Ll", 0.0);
  TTsemilep->Fill("1Ll", 0.0);
  TTsemilep->Fill("g1Mb1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100", 0.0);
  TTsemilep->Fill("1Mbg1W1LlmT100", 0.0);
  TTsemilep->Fill("g2Mbg1W1LlmT100", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhiHat4", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhi0p3", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhig0p3", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhi0p5", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhig0p5", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT", 0.0);
  TTsemilep->Fill("0Lb1Ll", 0.0);
  TTsemilep->Fill("0Lbg1Y1Ll", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT100", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhiHat4", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhiHatg4", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhi0p3", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhig0p3", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhi0p5", 0.0);
  TTsemilep->Fill("0Lbg1Y1LlmT_mdPhig0p5", 0.0);
  TTsemilep->Fill("presel2genallmu", 0.0);
  TTsemilep->Fill("presel2genmu", 0.0);
  TTsemilep->Fill("presel2mu", 0.0);
  TTsemilep->Fill("2genallmu", 0.0);
  TTsemilep->Fill("2genmu", 0.0);
  TTsemilep->Fill("2munoZmass", 0.0);
  TTsemilep->Fill("2mu", 0.0);
  TTsemilep->Fill("2mu0el", 0.0);
  TTsemilep->Fill("0Lb2mu0el", 0.0);
  TTsemilep->Fill("g1Mb2mu0el", 0.0);
  TTsemilep->Fill("0Lbg1Y2mu0el", 0.0);
  TTsemilep->Fill("g1Mbg1Y2mu0el", 0.0);
  TTsemilep->Fill("2elnoZmass", 0.0);
  TTsemilep->Fill("2el", 0.0);
  TTsemilep->Fill("2el0mu", 0.0);
  TTsemilep->Fill("0Lb2el0mu", 0.0);
  TTsemilep->Fill("g1Mb2el0mu", 0.0);
  TTsemilep->Fill("0Lbg1Y2el0mu", 0.0);
  TTsemilep->Fill("g1Mbg1Y2el0mu", 0.0);
  TTsemilep->Fill("2lnoZmass", 0.0);
  TTsemilep->Fill("2l", 0.0);
  TTsemilep->Fill("2l0ol", 0.0);
  TTsemilep->Fill("0Lb2l0ol", 0.0);
  TTsemilep->Fill("g1Mb2l0ol", 0.0);
  TTsemilep->Fill("0Lbg1Y2l0ol", 0.0);
  TTsemilep->Fill("g1Mbg1Y2l0ol", 0.0);


  TTdilep->Fill("NoCuts2genallmu", 0.0);
  TTdilep->Fill("NoCuts2genmu", 0.0);
  TTdilep->Fill("NoCuts2mu", 0.0);
  TTdilep->Fill("NoCuts", 0.0);
  TTdilep->Fill("Cleaning", 0.0);
  TTdilep->Fill("Pileup", 0.0);
  TTdilep->Fill("ISR", 0.0);
  TTdilep->Fill("TopPt", 0.0);
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
  TTdilep->Fill("1Mbg1W0Ll", 0.0);
  TTdilep->Fill("g2Mbg1W0Ll", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhiHat4", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhi0p3", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhig0p3", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhi0p5", 0.0);
  TTdilep->Fill("g1Mbg1W0Ll_mdPhig0p5", 0.0);
  TTdilep->Fill("g1Mb0Wg1uW0Ll", 0.0);
  TTdilep->Fill("0Lb0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhig0p3", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p5", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhig0p5", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat4", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhiHatg4", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat5", 0.0);
  TTdilep->Fill("0Lbg1W0Ll", 0.0);
  TTdilep->Fill("1Ll", 0.0);
  TTdilep->Fill("g1Mb1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100", 0.0);
  TTdilep->Fill("1Mbg1W1LlmT100", 0.0);
  TTdilep->Fill("g2Mbg1W1LlmT100", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhiHat4", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhi0p3", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhig0p3", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhi0p5", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT100_mdPhig0p5", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT", 0.0);
  TTdilep->Fill("0Lb1Ll", 0.0);
  TTdilep->Fill("0Lbg1Y1Ll", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT100", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhiHat4", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhiHatg4", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhi0p3", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhig0p3", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhi0p5", 0.0);
  TTdilep->Fill("0Lbg1Y1LlmT_mdPhig0p5", 0.0);
  TTdilep->Fill("presel2genallmu", 0.0);
  TTdilep->Fill("presel2genmu", 0.0);
  TTdilep->Fill("presel2mu", 0.0);
  TTdilep->Fill("2genallmu", 0.0);
  TTdilep->Fill("2genmu", 0.0);
  TTdilep->Fill("2munoZmass", 0.0);
  TTdilep->Fill("2mu", 0.0);
  TTdilep->Fill("2mu0el", 0.0);
  TTdilep->Fill("0Lb2mu0el", 0.0);
  TTdilep->Fill("g1Mb2mu0el", 0.0);
  TTdilep->Fill("0Lbg1Y2mu0el", 0.0);
  TTdilep->Fill("g1Mbg1Y2mu0el", 0.0);
  TTdilep->Fill("2elnoZmass", 0.0);
  TTdilep->Fill("2el", 0.0);
  TTdilep->Fill("2el0mu", 0.0);
  TTdilep->Fill("0Lb2el0mu", 0.0);
  TTdilep->Fill("g1Mb2el0mu", 0.0);
  TTdilep->Fill("0Lbg1Y2el0mu", 0.0);
  TTdilep->Fill("g1Mbg1Y2el0mu", 0.0);
  TTdilep->Fill("2lnoZmass", 0.0);
  TTdilep->Fill("2l", 0.0);
  TTdilep->Fill("2l0ol", 0.0);
  TTdilep->Fill("0Lb2l0ol", 0.0);
  TTdilep->Fill("g1Mb2l0ol", 0.0);
  TTdilep->Fill("0Lbg1Y2l0ol", 0.0);
  TTdilep->Fill("g1Mbg1Y2l0ol", 0.0);


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry)
    {
      // Read event into memory
      stream.read(entry);
      
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

      ofile.count("NoCuts", w);

      // -------------------------------------------------------------------------
      // -- Get rid of the noise in data before you start filling ANY histogram --
      // -- by applying the filters:                                            --
      // -------------------------------------------------------------------------

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

      // ---------------------------------
      // -- Trigger requirement in Data --
      // ---------------------------------

      if (eventhelper_isRealData==1) { 
	// Run2012A:
	if ( eventhelper_run >= 190456 && eventhelper_run < 190762 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 1 //||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v5 == 1 ||
		 //triggerresultshelper_HLT_HT750_v1 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v3 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 0 &&
		 (////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v5 == 1) //&&
		  //triggerresultshelper_HLT_HT750_v1 == 1 &&
		 ////triggerresultshelper_HLT_DiPFJetAve320_v3 == 0 
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v3 == 1 ||
		 ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v5 == 1 //||
		 //triggerresultshelper_HLT_HT750_v1 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v3 == 1
		 ) ) continue;
	  }
	}

	if ( eventhelper_run >= 190762 && eventhelper_run < 191512 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1 //||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v6 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v4 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 0 &&
		 (////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v6 == 1) //&&
		 //triggerresultshelper_HLT_HT750_v2 == 1 &&
		 ////triggerresultshelper_HLT_DiPFJetAve320_v4 == 0
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1 ||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v6 == 1 //||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v4 == 1
		 ) ) continue;
	  }
	}

	if ( eventhelper_run >= 191512 && eventhelper_run < 193746 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 //||
		 //triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		 //triggerresultshelper_HLT_PFHT650_v7 == 1 ||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v5 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 0 &&
		 (////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		  triggerresultshelper_HLT_PFHT650_v7 == 1) //&&
		 //triggerresultshelper_HLT_HT750_v2 == 1 &&
		 ////triggerresultshelper_HLT_DiPFJetAve320_v5 == 0
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 ||
		 ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4 == 1 ||
		 triggerresultshelper_HLT_PFHT650_v7 == 1 //||
		 //triggerresultshelper_HLT_HT750_v2 == 1 ||
		 ////triggerresultshelper_HLT_DiPFJetAve320_v5 == 1
		 ) ) continue;
	  }
	}

	// 2012B
	if ( eventhelper_run >= 193746 && eventhelper_run < 196039) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5 == 1 ||
	       triggerresultshelper_HLT_PFHT650_v8 == 1 //||
	       //triggerresultshelper_HLT_HT750_v3 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v6 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 196039 && eventhelper_run < 197770 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6 == 1 ||
	       triggerresultshelper_HLT_PFHT650_v9 == 1 //||
	       //triggerresultshelper_HLT_HT750_v4 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v6 == 1
	       ) ) continue;
	}

	// Run2012C:
	if ( eventhelper_run >= 197770 && eventhelper_run < 199648 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet400_v6 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v1 == 1 //||
	       //triggerresultshelper_HLT_HT750_v5 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v7 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 199648 && eventhelper_run < 202807 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v8 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v3 == 1 //||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v9 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 202807 && eventhelper_run < 203734 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v4 == 1 //||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v10 == 1
	       ) ) continue;
	}

	// Run2012D:
	if ( eventhelper_run >= 203734 && eventhelper_run < 208940 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1 ||
	       ////triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10 == 1 ||
	       triggerresultshelper_HLT_PFNoPUHT650_v4 == 1 //||
	       //triggerresultshelper_HLT_HT750_v7 == 1 ||
	       ////triggerresultshelper_HLT_DiPFJetAve400_v10 == 1
	       ) ) continue;
	}
      }

      // ------------------------
      // -- Pileup reweighting --
      // ------------------------

      double num_vertices = pileupsummaryinfo[0].getTrueNumInteractions;
      h_TrueNumVertices->Fill(num_vertices,w);
      // get bin number in pileup histogram
      int pileup_bin = (int)ceil(num_vertices);
      double w_pileup = 1.;
      if(doPileupReweighting)
	w_pileup = h_pileup->GetBinContent(pileup_bin);      
      h_TrueNumVertices_reweighted->Fill(num_vertices,w*w_pileup);

      double w_old = w;
      w = w*w_pileup;

      ofile.count("Pileup",w);

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
        btagCSVMEEFull(partonFlavour, pt, fabs(eta), SFCSVMFl, dSFCSVMFl);
        btagCSVMEEFast(partonFlavour, pt, fabs(eta), SFCSVMFs, dSFCSVMFs);

        btagCSVLEEFull(partonFlavour, pt, fabs(eta), SFCSVLFl, dSFCSVLFl);
        btagCSVLEEFast(partonFlavour, pt, fabs(eta), SFCSVLFs, dSFCSVLFs);

	double eCSVM = 0;
	double eCSVL = 0;
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
	  if (fabs(partonFlavour) != 4 && fabs(partonFlavour) != 5) {
	    eCSVM = geteff2D(h_pt_eta_l_CSVMeff, pt, fabs(eta));
	    eCSVL = geteff2D(h_pt_eta_l_CSVLeff, pt, fabs(eta));
	  }
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
        h_jmass_jpt->Fill(jethelper4[i].mass, jethelper4[i].pt);
        h_d1pt_d2pt->Fill(jethelper4[i].daughter_0_pt, jethelper4[i].daughter_1_pt);
        h_d1m_d2m->Fill(jethelper4[i].daughter_0_mass, jethelper4[i].daughter_1_mass);
        h_jmass->Fill(jethelper4[i].mass);
	// New Andreas cuts:
        //if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
        if (!(jethelper4[i].mass > 65 && jethelper4[i].mass < 105)) continue;
	sY.push_back(jethelper4[i]);
        h_d1ptsel_d2ptsel->Fill(jethelper4[i].daughter_0_pt, jethelper4[i].daughter_1_pt);
        h_d1msel_d2msel->Fill(jethelper4[i].daughter_0_mass, jethelper4[i].daughter_1_mass);
        sjet2.push_back(jethelper4[i]);
        double massdrop = 1;
        double daughmass = -9;
        if (jethelper4[i].daughter_0_mass > jethelper4[i].daughter_1_mass) {
          daughmass = jethelper4[i].daughter_0_mass;
        } else {
          daughmass = jethelper4[i].daughter_1_mass;
        };
        massdrop = daughmass / jethelper4[i].mass;
        h_mdrp->Fill(massdrop, w);
        h_mdrp_jpt->Fill(massdrop, jethelper4[i].pt);
        double massdrop2 = 1;
        double daughmass2 = -9;
        if (jethelper4[i].daughter_0_pt > jethelper4[i].daughter_1_pt) {
          daughmass2 = jethelper4[i].daughter_0_mass;
        } else {
          daughmass2 = jethelper4[i].daughter_1_mass;
        };
        massdrop2 = daughmass2 / jethelper4[i].mass;
        h_mdrp2->Fill(massdrop2, w);
        h_mdrp2_jpt->Fill(massdrop2, jethelper4[i].pt);
        double dRd1d2 = fdeltaR(jethelper4[i].daughter_0_eta,
                                jethelper4[i].daughter_0_phi,
                                jethelper4[i].daughter_1_eta,
                                jethelper4[i].daughter_1_phi
                                );
        double yasym = (TMath::Min(pow(jethelper4[i].daughter_0_pt, 2),
                                  pow(jethelper4[i].daughter_1_pt, 2))*
                       pow(dRd1d2,2))/
          pow(jethelper4[i].mass,2);
        h_yasym->Fill(yasym, w);
        h_mdrp_yasym->Fill(massdrop, yasym);
        h_mdrp2_yasym->Fill(massdrop2, yasym);

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

      // W selection:
      std::vector<jethelper4_s> sWAK5;
      for (unsigned int i=0; i<jethelper4.size(); i++) {
        if (!(jethelper4[i].pt > 30) ) continue;
        if (!(fabs(jethelper4[i].eta) < 3) ) continue;

	// New Andreas cuts:
        //if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
        if (!(jethelper4[i].mass > 65 && jethelper4[i].mass < 105)) continue;
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
      // From supposed to be https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation#TauID
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


      // ---------------------------------------------------------------
      // -- Now define some more variables, and fill a few histograms --
      // ---------------------------------------------------------------

      // look at number of Ws and bs before any selection (besides cleaning and trigger)
      double nWs = sW.size();
      double nYs = sY.size();
      double nbs = sbjet.size();

      h_nW_Pileup->Fill(nWs, w);
      h_nb_Pileup->Fill(nbs, w);
      h_nW_nb_Pileup->Fill(nWs, nbs, w);
      h_nWAK5_Pileup->Fill(sWAK5.size(), w);

      h_met_Pileup->Fill(cmgbasemet2[0].et, w);
      if (sjet.size() > 0)
        h_j1pt_Pileup->Fill(sjet[0].pt, w);
      h_nj_Pileup->Fill(sjet.size(), w);

      
      // -------------------------------
      // -- Calculate Razor variables --
      // -------------------------------

      // Calculate MR and R2 ignoring leptons -- DEFAULT
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
	h_MR_Pileup->Fill(MR, w);
	MTR = CalcMTR(LVhemis[0], LVhemis[1], V3met);
	R2 = pow((MTR / MR),2);
	h_R2_Pileup->Fill(R2, w);
	h_MR_R2_Pileup->Fill(MR, R2, w);
	h_met_R2_Pileup->Fill(cmgbasemet2[0].et, R2, w);
      }
      
      // Calculate MR and R2 adding mus to MET
      TVector3 V3metmu;
      V3metmu.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      for (unsigned int i=0; i<V3mu.size(); i++) {
        V3metmu += V3mu[i];
      }
      h_metmu_Pileup->Fill(V3metmu.Pt(), w);

      double MTRmetmu = -9999;
      double R2metmu = -9999;
      if (LVhemis.size() == 2) {
        MTRmetmu = CalcMTR(LVhemis[0], LVhemis[1], V3metmu);
        R2metmu = pow((MTRmetmu / MR),2);
        h_R2metmu_Pileup->Fill(R2metmu, w);
        h_MR_R2metmu_Pileup->Fill(MR, R2, w);
        h_metmu_R2metmu_Pileup->Fill(V3metmu.Pt(), R2, w);
      }

      // Calculate MR and R2 adding electrons to MET
      TVector3 V3metel;
      V3metel.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      for (unsigned int i=0; i<V3el.size(); i++) {
        V3metel += V3el[i];
      }
      h_metel_Pileup->Fill(V3metel.Pt(), w);

      double MTRmetel = -9999;
      double R2metel = -9999;
      if (LVhemis.size() == 2) {
        MTRmetel = CalcMTR(LVhemis[0], LVhemis[1], V3metel);
        R2metel = pow((MTRmetel / MR),2);
        h_R2metel_Pileup->Fill(R2metel, w);
        h_MR_R2metel_Pileup->Fill(MR, R2, w);
        h_metel_R2metel_Pileup->Fill(V3metel.Pt(), R2, w);
      }
     

      // --------------------------------------------------------
      // -- Everything computed from generator level particles --
      // --------------------------------------------------------

      // find the tops and Ws and plot some basic quantities
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
	    h_gen_toppt->Fill(genparticlehelper[i].pt);
            double dRWb = fdeltaR(genparticlehelper[id1].eta, genparticlehelper[id1].phi,
				  genparticlehelper[id2].eta, genparticlehelper[id2].phi);
	    h_gen_dRWb->Fill(dRWb);
            h_gen_toppt_dRWb->Fill(genparticlehelper[i].pt, dRWb);
            for (unsigned int j=0; j<topdaughters.size(); j++) {
              if (fabs(topdaughters[j].pdgId)==24) {
                Ws.push_back(topdaughters[j]);
                h_gen_Wpt->Fill(topdaughters[j].pt);
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
	    double dRqq = fdeltaR(genparticlehelper[iWd1].eta, genparticlehelper[iWd1].phi,
				  genparticlehelper[iWd2].eta, genparticlehelper[iWd2].phi);
	    h_gen_dRqq->Fill(dRqq);
            h_gen_Wpt_dRqq->Fill(Ws[i].pt, dRqq);
	    
	    // make some plots for the Wtag efficiency, given a hadronic W
	    // first match gen W with CA8 jet
	    double dRgenCA8 = 999.;
	    int matched_index = -1;
	    for (int j=0; j<(int)jethelper4.size(); ++j){
	      if (!(jethelper4[i].pt > 30) ) continue;
	      if (!(fabs(jethelper4[i].eta) < 3) ) continue;
	      double dRgenCA8_temp = fdeltaR(Ws[i].eta, Ws[i].phi,
					     jethelper4[j].eta,jethelper4[j].phi);
	      if (dRgenCA8_temp < dRgenCA8 && dRgenCA8_temp < 0.8){
		dRgenCA8 = dRgenCA8_temp;
		matched_index = j;
	      }
	    }
	    if (matched_index != -1){ // we found a match
	      // now check whether there is a b from the top in the cone
	      bool found_bjet_in_cone = false;
	      for(int bjet=0; bjet<genparticlehelper.size(); bjet++){
		if(fabs(genparticlehelper[bjet].pdgId) == 5 && genparticlehelper[bjet].status == 3 && fabs(genparticlehelper[genparticlehelper[bjet].firstMother].pdgId) == 6){
		  double dR_CA8_b = fdeltaR(genparticlehelper[bjet].eta, genparticlehelper[bjet].phi,jethelper4[matched_index].eta,jethelper4[matched_index].phi);
		  if (dR_CA8_b < 0.8)
		    found_bjet_in_cone = true;
		}
	      }
	      if (!found_bjet_in_cone){
		h_pt_closestCA8jet->Fill(jethelper4[matched_index].pt,w);
		// now check whether the matched jet has mass in proper window
		if (jethelper4[matched_index].mass > 65 && jethelper4[matched_index].mass < 105){
		  h_pt_closestCA8jet_Wmass->Fill(jethelper4[matched_index].pt,w);
		  
		  // now check whether it is also tagged
		  // Match with the unpruned:
		  double prjmatch = 0;
		  int jpr = -1;
		  double dRmn = 100;
		  for (unsigned int j=0; j<jethelper5.size(); j++) {
		    double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, jethelper4[matched_index].eta, jethelper4[matched_index].phi);
		    if (dR < 0.7 && dR < dRmn) {
		      dRmn = dR;
		      prjmatch = 1;
		      jpr = j;
		      break;
		    }
		  }
		  if (prjmatch==1){
		    double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
		    double tau3 = jethelper5[jpr].tau3;
		    if (tau21 < 0.50) {
		    //if (tau21 < 0.46) {
		      //if (tau3 < 0.135){
			h_pt_closestCA8jet_Wmass_tagged->Fill(jethelper4[matched_index].pt,w);
		      //}
		    }
		  }
		}
	      }
	    }
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
	TTallhad->Fill("Pileup",w);
      else if(isTTsemilep)
	TTsemilep->Fill("Pileup",w);
      else if(isTTdilep)
	TTdilep->Fill("Pileup",w);
      h_PV_Pileup->Fill(eventhelperextra_numberOfPrimaryVertices,w);

      // ------------------------------------
      // -- Pick up trigger weights for MC --
      // ------------------------------------

      double w_trigger = 1.;
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
		//cout << xmin << " " << MR << " " << xmax << " " << ymin << " " << R2 << " " << ymax << " " << whlt << endl;
		break;
	      }
	    }
	  }
	  h_PV_1j->Fill(eventhelperextra_numberOfPrimaryVertices,w_old*w_trigger); // Included Trigger, But No PU weights
	  h_PV_reweighted_1j->Fill(eventhelperextra_numberOfPrimaryVertices,w*w_trigger); // current weight, includes Trigger, PU
	  h_HT_jet1pt_nosel->Fill(HT,sjet[0].pt,w*w_trigger);
	  h_HT_jet1pt_nosel_notrigw->Fill(HT,sjet[0].pt,w);
	  h_TrueNumVertices_reweighted_withtrigger->Fill(num_vertices,w*w_trigger);
	}
      }

      if (eventhelper_isRealData!=0) {
	if (sjet.size() > 0){
	  h_PV_1j->Fill(eventhelperextra_numberOfPrimaryVertices,w); 
	  h_PV_reweighted_1j->Fill(eventhelperextra_numberOfPrimaryVertices,w); 
	}
      }
      h_PV_nosel->Fill(eventhelperextra_numberOfPrimaryVertices,w_old); // No PU or trigger weights
      h_PV_reweighted_nosel->Fill(eventhelperextra_numberOfPrimaryVertices,w); // current weight, includes PU
      
      h_PV_HT_nosel->Fill(eventhelperextra_numberOfPrimaryVertices,HT,w_old); // No PU or trigger weights
      h_PV_HT_reweighted_nosel->Fill(eventhelperextra_numberOfPrimaryVertices,HT,w); // current weight, includes PU
      
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
	  if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80"
	      || sample == "T1ttcc_old" || sample == "T1t1t" )
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
      h_totalISRweight_nominal->Fill(1,w_ISR_nominal);
      h_totalISRweight_up->Fill(1,w_ISR_up);
      h_totalISRweight_down->Fill(1,w_ISR_down);
      
      if (doISRreweighting)
	w = w*w_ISR_nominal;
     
      ofile.count("ISR",w);
      if(isTTallhad)
	TTallhad->Fill("ISR",w);
      else if(isTTsemilep)
	TTsemilep->Fill("ISR",w);
      else if(isTTdilep)
	TTdilep->Fill("ISR",w);
      
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
      h_totalTopPTweight_nominal->Fill(1,w_TopPt_nominal);
      h_totalTopPTweight_up->Fill(1,w_TopPt_up);
      h_totalTopPTweight_down->Fill(1,w_TopPt_down);
      
      if(doTopPtreweighting)
	w = w*w_TopPt_nominal;

      ofile.count("TopPt",w);
      if(isTTallhad)
	TTallhad->Fill("TopPt",w);
      else if(isTTsemilep)
	TTsemilep->Fill("TopPt",w);
      else if(isTTdilep)
	TTdilep->Fill("TopPt",w);


      // -------------------------------
      // -- Stuff for Zinv estimation --
      // -------------------------------

      // GEN Z selection
      std::vector<genparticlehelper_s> genZ;      
      for (unsigned int i=0; i<genparticlehelper.size(); i++) {
        if (genparticlehelper[i].status != 3) continue;
        if (!(fabs(genparticlehelper[i].pdgId) == 23) ) continue;
	genZ.push_back(genparticlehelper[i]);
      }

      // GEN muon selection (to be used for dimuon acceptance)
      // We select status 1 particles, because in our DY sample, there seem to be
      // no extra muons...
      std::vector<genparticlehelper_s> genmu;
      std::vector<genparticlehelper_s> sgenmu;
      std::vector<TLorentzVector> LVsgenmu;
      for (unsigned int i=0; i<genparticlehelper.size(); i++) {
        if (!(fabs(genparticlehelper[i].pdgId) == 13) ) continue;
        if (!(genparticlehelper[i].status != 1) ) continue;
	genmu.push_back(genparticlehelper[i]);
	if (!(genparticlehelper[i].pt > 10.) ) continue;
	if (!(fabs(genparticlehelper[i].eta) < 2.4) ) continue;
	//cout << genparticlehelper[i].status << " " << genparticlehelper[i].pt << " " 
	//     << genparticlehelper[i].eta << endl;
	sgenmu.push_back(genparticlehelper[i]);
	TLorentzVector lvmu;
	lvmu.SetPtEtaPhiM(genparticlehelper[i].pt, genparticlehelper[i].eta,
			  genparticlehelper[i].phi, genparticlehelper[i].mass);
        LVsgenmu.push_back(lvmu);
      }


      // ---------------------
      // -- event selection --
      // ---------------------

      // dimuon selections for the muon acceptance and efficiency - for no cuts
      if (genmu.size() == 2 && genZ.size() == 1) {
	if (genmu[0].charge == -genmu[1].charge){
	  ofile.count("NoCuts2genallmu", w);
	  if(isTTallhad)
	    TTallhad->Fill("NoCuts2genallmu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("NoCuts2genallmu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("NoCuts2genallmu", w);
	}
      }
      if (sgenmu.size() == 2) {
	if (genmu[0].charge == -genmu[1].charge){
	  TLorentzVector LVZcand;
	  for (unsigned int m=0; m<LVsgenmu.size(); m++) {
	    LVZcand += LVsgenmu[m];
	  }
	  double Zmass = LVZcand.M();
	  if (Zmass >= 60 && Zmass <= 120){
	    ofile.count("NoCuts2genmu", w);
	    if(isTTallhad)
	      TTallhad->Fill("NoCuts2genmu", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("NoCuts2genmu", w);
	    else if(isTTdilep)
	      TTdilep->Fill("NoCuts2genmu", w);
	  }
	}
      }
      if (smuon.size() == 2) {
	TLorentzVector LVZcand;
	for (unsigned int m=0; m<LVmu.size(); m++) {
	  LVZcand += LVmu[m];
	}
	double Zmass = LVZcand.M();
	if (Zmass >= 60 && Zmass <= 120){
	  ofile.count("NoCuts2mu", w);
	  if(isTTallhad)
	    TTallhad->Fill("NoCuts2mu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("NoCuts2mu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("NoCuts2mu", w);
	}
      }

      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;
      ofile.count("HCAL_noise", w);
      h_MR_HCAL_noise->Fill(MR, w);
      h_R2_HCAL_noise->Fill(R2, w);
      h_MR_R2_HCAL_noise->Fill(MR, R2, w);
      h_PV_HCAL_noise->Fill(eventhelperextra_numberOfPrimaryVertices,w);
      if(isTTallhad)
	TTallhad->Fill("HCAL_noise", w);
      else if(isTTsemilep)
	TTsemilep->Fill("HCAL_noise", w);
      else if(isTTdilep)
	TTdilep->Fill("HCAL noise", w);

      // at least one good primary vertex
      if (!(svertex.size() > 0)) continue;
      ofile.count("vertexg0", w);
      h_MR_vertexg0->Fill(MR, w);
      h_R2_vertexg0->Fill(R2, w);
      h_MR_R2_vertexg0->Fill(MR, R2, w);
      h_PV_vertexg0->Fill(eventhelperextra_numberOfPrimaryVertices,w);
      if(isTTallhad)
	TTallhad->Fill("vertexg0", w);
      else if(isTTsemilep)
	TTsemilep->Fill("vertexg0", w);
      else if(isTTdilep)
	TTdilep->Fill("vertexg0", w);

      // at least three jets
      if (!(sjet.size() >= 3)) continue;
      ofile.count("njetge3", w);
      h_MR_njetge3->Fill(MR, w);
      h_R2_njetge3->Fill(R2, w);
      h_MR_R2_njetge3->Fill(MR, R2, w);
      h_PV_njetge3->Fill(eventhelperextra_numberOfPrimaryVertices,w);
      if(isTTallhad)
	TTallhad->Fill("njetge3", w);
      else if(isTTsemilep)
	TTsemilep->Fill("njetge3", w);
      else if(isTTdilep)
	TTdilep->Fill("njetge3", w);
      
      // Apply HLT weight and include it in the total weight:
      double w_nohlt = w;
      w = w*w_trigger;
 
      ofile.count("HLT", w);
      h_MR_HLT->Fill(MR, w);
      h_R2_HLT->Fill(R2, w);
      h_MR_R2_HLT->Fill(MR, R2, w);
      if(isTTallhad)
	TTallhad->Fill("HLT", w);
      else if(isTTsemilep)
	TTsemilep->Fill("HLT", w);
      else if(isTTdilep)
	TTdilep->Fill("HLT", w);

      h_PV_HLT->Fill(eventhelperextra_numberOfPrimaryVertices,w_old*w_trigger*w_TopPt_nominal*w_ISR_nominal);
      h_PV_reweighted_HLT->Fill(eventhelperextra_numberOfPrimaryVertices,w);
      
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
      
      h_minDeltaPhi_HLT->Fill(minDeltaPhi, w);
      h_MR_minDeltaPhi_HLT->Fill(MR, minDeltaPhi, w);
      h_R2_minDeltaPhi_HLT->Fill(R2, minDeltaPhi, w);
	
      h_minDeltaPhiHat_HLT->Fill(minDeltaPhiHat, w);
      h_MR_minDeltaPhiHat_HLT->Fill(MR, minDeltaPhiHat, w);
      h_R2_minDeltaPhiHat_HLT->Fill(R2, minDeltaPhiHat, w);
      
      // pt of first jet greater than 200 GeV
      if (!(sjet[0].pt > 200)) continue;
      ofile.count("jet1ptg200", w);
      h_MR_jet1ptg200->Fill(MR, w);
      h_R2_jet1ptg200->Fill(R2, w);
      h_MR_R2_jet1ptg200->Fill(MR, R2, w);
      h_PV_jet1ptg200->Fill(eventhelperextra_numberOfPrimaryVertices,w);

      h_njets_jet1ptg200->Fill(sjet.size(),w);
      h_nbjets_jet1ptg200->Fill(sbjet.size(),w);
      h_met_jet1ptg200->Fill(met.Pt(),w);
      h_jet1pt_jet1ptg200->Fill(sjet[0].pt,w);
      h_jet2pt_jet1ptg200->Fill(sjet[1].pt,w);
      h_jet3pt_jet1ptg200->Fill(sjet[2].pt,w);
      
      if(isTTallhad)
	TTallhad->Fill("jet1ptg200", w);
      else if(isTTsemilep)
	TTsemilep->Fill("jet1ptg200", w);
      else if(isTTdilep)
	TTdilep->Fill("jet1ptg200", w);

      h_nW_jet1ptg200->Fill(nWs, w);
      h_nb_jet1ptg200->Fill(nbs, w);
      h_nW_nb_jet1ptg200->Fill(nWs, nbs, w);
      h_nWAK5_jet1ptg200->Fill(sWAK5.size(), w);

      h_minDeltaPhi_jet1ptg200->Fill(minDeltaPhi, w);
      h_MR_minDeltaPhi_jet1ptg200->Fill(MR, minDeltaPhi, w);
      h_R2_minDeltaPhi_jet1ptg200->Fill(R2, minDeltaPhi, w);
      
      h_minDeltaPhiHat_jet1ptg200->Fill(minDeltaPhiHat, w);
      h_MR_minDeltaPhiHat_jet1ptg200->Fill(MR, minDeltaPhiHat, w);
      h_R2_minDeltaPhiHat_jet1ptg200->Fill(R2, minDeltaPhiHat, w);

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
	h_MR_SIG->Fill(MR, w);
	h_R2_SIG->Fill(R2, w);
	h_MR_R2_SIG->Fill(MR, R2, w);

	h_minDeltaPhi_SIG->Fill(minDeltaPhi, w);
	h_MR_minDeltaPhi_SIG->Fill(MR, minDeltaPhi, w);
	h_R2_minDeltaPhi_SIG->Fill(R2, minDeltaPhi, w);
	
	h_minDeltaPhiHat_SIG->Fill(minDeltaPhiHat, w);
	h_MR_minDeltaPhiHat_SIG->Fill(MR, minDeltaPhiHat, w);
	h_R2_minDeltaPhiHat_SIG->Fill(R2, minDeltaPhiHat, w);

	h_HT_SIG->Fill(HT, w);
	h_HT_jet1pt_SIG->Fill(HT,sjet[0].pt,w);
	h_HT_jet1pt_SIG_notrigw->Fill(HT,sjet[0].pt,w_nohlt);

	h_PV_SIG->Fill(eventhelperextra_numberOfPrimaryVertices,w_old*w_trigger*w_TopPt_nominal*w_ISR_nominal);
	h_PV_reweighted_SIG->Fill(eventhelperextra_numberOfPrimaryVertices,w);
	
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
	  h_MR_neleeq0->Fill(MR, w);
	  h_R2_neleeq0->Fill(R2, w);
	  h_MR_R2_neleeq0->Fill(MR, R2, w);
	  if(isTTallhad)
	    TTallhad->Fill("neleeq0", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("neleeq0", w);
	  else if(isTTdilep)
	    TTdilep->Fill("neleeq0", w);
	  
	  if (nloosemuons == 0) {
	    ofile.count("nmueq0", w);
	    h_MR_nmueq0->Fill(MR, w);
	    h_R2_nmueq0->Fill(R2, w);
	    h_MR_R2_nmueq0->Fill(MR, R2, w);
	    if(isTTallhad)
	      TTallhad->Fill("nmueq0", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("nmueq0", w);
	    else if(isTTdilep)
	      TTdilep->Fill("nmueq0", w);
	    
	    if (eventhelperextra_trackIso == 0){
	      ofile.count("trackIso", w);
	      h_MR_trackIso->Fill(MR, w);
	      h_R2_trackIso->Fill(R2, w);
	      h_MR_R2_trackIso->Fill(MR, R2, w);
	      if(isTTallhad)
		TTallhad->Fill("trackIso", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("trackIso", w);
	      else if(isTTdilep)
		TTdilep->Fill("trackIso", w);
	      
	      h_minDeltaPhi_0Ll->Fill(minDeltaPhi, w);
	      h_MR_minDeltaPhi_0Ll->Fill(MR, minDeltaPhi, w);
	      h_R2_minDeltaPhi_0Ll->Fill(R2, minDeltaPhi, w);
	      
	      h_minDeltaPhiHat_0Ll->Fill(minDeltaPhiHat, w);
	      h_MR_minDeltaPhiHat_0Ll->Fill(MR, minDeltaPhiHat, w);
	      h_R2_minDeltaPhiHat_0Ll->Fill(R2, minDeltaPhiHat, w);

	      h_PV_0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);

	      if (nmediumbs > 0){
		// apply btag weight if not data
		if (eventhelper_isRealData!=1) {
		  w = w*wCSVM;
		}
		ofile.count("g1Mb0Ll", w);
		h_MR_g1Mb0Ll->Fill(MR, w);
		h_R2_g1Mb0Ll->Fill(R2, w);
		h_MR_R2_g1Mb0Ll->Fill(MR, R2, w);
		h_PV_g1Mb0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mb0Ll", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mb0Ll", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mb0Ll", w);


		double gen_toppt = 0;
		double dRWb = 0;
		double gen_Wpt = 0;
		double dRqq = 0;
		for (unsigned int i=0; i<genparticlehelper.size(); i++) {
		  if (genparticlehelper[i].status != 3) continue;
		  if (fabs(genparticlehelper[i].pdgId) == 6) {
		    std::vector<genparticlehelper_s> topdaughters;
		    int id1 = genparticlehelper[i].firstDaughter;
		    int id2 = genparticlehelper[i].lastDaughter;
		    topdaughters.push_back(genparticlehelper[id1]);
		    topdaughters.push_back(genparticlehelper[id2]);

		    if ((fabs(topdaughters[0].pdgId) == 5 || fabs(topdaughters[0].pdgId) == 24) 
			&& (fabs(topdaughters[1].pdgId) == 5 || fabs(topdaughters[1].pdgId) == 24)) {
		      double toppt = genparticlehelper[i].pt;
		      double dRWb_ = fdeltaR(genparticlehelper[id1].eta, genparticlehelper[id1].phi,
					    genparticlehelper[id2].eta, genparticlehelper[id2].phi);
		      if (toppt > gen_toppt){
			gen_toppt = toppt;
			dRWb = dRWb_;
		      }
		    }
		  }
		  if (fabs(genparticlehelper[i].pdgId) == 24) {
		    std::vector<genparticlehelper_s> Wdaughters;
		    int id1 = genparticlehelper[i].firstDaughter;
		    int id2 = genparticlehelper[i].lastDaughter;
		    Wdaughters.push_back(genparticlehelper[id1]);
		    Wdaughters.push_back(genparticlehelper[id2]);

		    if ((fabs(Wdaughters[0].pdgId) <= 5 ) && (fabs(Wdaughters[1].pdgId) <= 5)) {
		      double Wpt = genparticlehelper[i].pt;
		      double dRqq_ = fdeltaR(genparticlehelper[id1].eta, genparticlehelper[id1].phi,
					    genparticlehelper[id2].eta, genparticlehelper[id2].phi);
		      if (Wpt > gen_Wpt){
			gen_Wpt = Wpt;
			dRqq = dRqq_;
		      }
		    }
		  }
		  
		}
		h_gen_top1pt_g1Mb0Ll->Fill(gen_toppt,w);
		h_gen_dRWb_g1Mb0Ll->Fill(dRWb,w);
		h_gen_top1pt_dRWb_g1Mb0Ll->Fill(gen_toppt, dRWb,w);
		h_gen_W1pt_g1Mb0Ll->Fill(gen_Wpt,w);
		h_gen_dRqq_g1Mb0Ll->Fill(dRqq,w);
		h_gen_W1pt_dRqq_g1Mb0Ll->Fill(gen_Wpt,dRqq,w);
		
		
		// g1Mb g1W 0Ll -- SIGNAL region
		if( sW.size() > 0){
		  ofile.count("g1Mbg1W0Ll",w);
		  h_MR_g1Mbg1W0Ll->Fill(MR, w);
		  h_R2_g1Mbg1W0Ll->Fill(R2, w);
		  h_MR_R2_g1Mbg1W0Ll->Fill(MR, R2, w);
		  h_njets_g1Mbg1W0Ll->Fill(sjet.size(),w);
		  h_nbjets_g1Mbg1W0Ll->Fill(nmediumbs,w);
		  h_met_g1Mbg1W0Ll->Fill(met.Pt(),w);
		  h_jet1pt_g1Mbg1W0Ll->Fill(sjet[0].pt,w);
		  h_jet2pt_g1Mbg1W0Ll->Fill(sjet[1].pt,w);
		  h_jet3pt_g1Mbg1W0Ll->Fill(sjet[2].pt,w);
		  h_dphimegajets_g1Mbg1W0Ll->Fill(fdeltaPhi(LVhemis[0].Phi(),LVhemis[1].Phi()),w);
		  h_dphijet1jet2_g1Mbg1W0Ll->Fill(fdeltaPhi(sjet[0].phi,sjet[1].phi),w);
		  h_HT_g1Mbg1W0Ll->Fill(HT,w);
		  h_PV_g1Mbg1W0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
		  h_minDeltaPhiHat_g1Mbg1W0Ll->Fill(minDeltaPhiHat, w);
		  h_MR_minDeltaPhiHat_g1Mbg1W0Ll->Fill(MR, minDeltaPhiHat, w);
		  h_R2_minDeltaPhiHat_g1Mbg1W0Ll->Fill(R2, minDeltaPhiHat, w);
		  h_minDeltaPhi_g1Mbg1W0Ll->Fill(minDeltaPhi, w);
		  h_MR_minDeltaPhi_g1Mbg1W0Ll->Fill(MR, minDeltaPhi, w);
		  h_R2_minDeltaPhi_g1Mbg1W0Ll->Fill(R2, minDeltaPhi, w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W0Ll", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W0Ll", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W0Ll", w);
		  
		  // separate ==1b and >=2b plots
		  if (nmediumbs == 1){
		    ofile.count("1Mbg1W0Ll",w);
		    h_MR_1Mbg1W0Ll->Fill(MR, w);
		    h_R2_1Mbg1W0Ll->Fill(R2, w);
		    h_MR_R2_1Mbg1W0Ll->Fill(MR, R2, w);
		    if(isTTallhad)
		      TTallhad->Fill("1Mbg1W0Ll", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("1Mbg1W0Ll", w);
		    else if(isTTdilep)
		      TTdilep->Fill("1Mbg1W0Ll", w);
		  } else { 
		    ofile.count("g2Mbg1W0Ll",w);
		    h_MR_g2Mbg1W0Ll->Fill(MR, w);
		    h_R2_g2Mbg1W0Ll->Fill(R2, w);
		    h_MR_R2_g2Mbg1W0Ll->Fill(MR, R2, w);
		    if(isTTallhad)
		      TTallhad->Fill("g2Mbg1W0Ll", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g2Mbg1W0Ll", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g2Mbg1W0Ll", w);
		  }

		  // add cut on mindeltaphihat
		  if (minDeltaPhiHat > 4){
		    ofile.count("g1Mbg1W0Ll_mdPhiHatg4",w);
		    h_MR_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, R2, w);
		    h_uw_MR_R2_g1Mbg1W0Ll_mdPhiHatg4->Fill(MR, R2, 1.);

		    h_njets_g1Mbg1W0Ll_mdPhiHatg4->Fill(sjet.size(),w);
		    h_nbjets_g1Mbg1W0Ll_mdPhiHatg4->Fill(nmediumbs,w);
		    h_met_g1Mbg1W0Ll_mdPhiHatg4->Fill(met.Pt(),w);
		    h_jet1pt_g1Mbg1W0Ll_mdPhiHatg4->Fill(sjet[0].pt,w);
		    h_jet2pt_g1Mbg1W0Ll_mdPhiHatg4->Fill(sjet[1].pt,w);
		    h_jet3pt_g1Mbg1W0Ll_mdPhiHatg4->Fill(sjet[2].pt,w);
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhiHatg4", w);
		  } else {
		    ofile.count("g1Mbg1W0Ll_mdPhiHat4",w);
		    h_MR_g1Mbg1W0Ll_mdPhiHat4->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhiHat4->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhiHat4->Fill(MR, R2, w);
		    
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhiHat4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhiHat4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhiHat4", w);
		  }

		  // add cut on mindeltaphi
		  if (minDeltaPhi > 0.3){
		    ofile.count("g1Mbg1W0Ll_mdPhig0p3",w);
		    h_MR_g1Mbg1W0Ll_mdPhig0p3->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhig0p3->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhig0p3->Fill(MR, R2, w);
		    h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p3->Fill(MR, R2, 1.);

		    h_njets_g1Mbg1W0Ll_mdPhig0p3->Fill(sjet.size(),w);
		    h_nbjets_g1Mbg1W0Ll_mdPhig0p3->Fill(nmediumbs,w);
		    h_met_g1Mbg1W0Ll_mdPhig0p3->Fill(met.Pt(),w);
		    h_jet1pt_g1Mbg1W0Ll_mdPhig0p3->Fill(sjet[0].pt,w);
		    h_jet2pt_g1Mbg1W0Ll_mdPhig0p3->Fill(sjet[1].pt,w);
		    h_jet3pt_g1Mbg1W0Ll_mdPhig0p3->Fill(sjet[2].pt,w);
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhig0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhig0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhig0p3", w);
		  } else {
		    ofile.count("g1Mbg1W0Ll_mdPhi0p3",w);
		    h_MR_g1Mbg1W0Ll_mdPhi0p3->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhi0p3->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhi0p3->Fill(MR, R2, w);
		    
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhi0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhi0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhi0p3", w);
		  }

		  // add cut on mindeltaphi
		  if (minDeltaPhi > 0.5){
		    ofile.count("g1Mbg1W0Ll_mdPhig0p5",w);
		    h_MR_g1Mbg1W0Ll_mdPhig0p5->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhig0p5->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhig0p5->Fill(MR, R2, w);
		    h_uw_MR_R2_g1Mbg1W0Ll_mdPhig0p5->Fill(MR, R2, 1.);

		    h_njets_g1Mbg1W0Ll_mdPhig0p5->Fill(sjet.size(),w);
		    h_nbjets_g1Mbg1W0Ll_mdPhig0p5->Fill(nmediumbs,w);
		    h_met_g1Mbg1W0Ll_mdPhig0p5->Fill(met.Pt(),w);
		    h_jet1pt_g1Mbg1W0Ll_mdPhig0p5->Fill(sjet[0].pt,w);
		    h_jet2pt_g1Mbg1W0Ll_mdPhig0p5->Fill(sjet[1].pt,w);
		    h_jet3pt_g1Mbg1W0Ll_mdPhig0p5->Fill(sjet[2].pt,w);
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhig0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhig0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhig0p5", w);
		  } else {
		    ofile.count("g1Mbg1W0Ll_mdPhi0p5",w);
		    h_MR_g1Mbg1W0Ll_mdPhi0p5->Fill(MR, w);
		    h_R2_g1Mbg1W0Ll_mdPhi0p5->Fill(R2, w);
		    h_MR_R2_g1Mbg1W0Ll_mdPhi0p5->Fill(MR, R2, w);
		    
		    if(isTTallhad)
		      TTallhad->Fill("g1Mbg1W0Ll_mdPhi0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mbg1W0Ll_mdPhi0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mbg1W0Ll_mdPhi0p5", w);
		  }

		} // end of sW.size() > 0
		else { // sW.size() == 0
		  if( aW.size() > 0) {
		    ofile.count("g1Mb0Wg1uW0Ll",w);
		    h_MR_g1Mb0Wg1uW0Ll->Fill(MR, w);
		    h_R2_g1Mb0Wg1uW0Ll->Fill(R2, w);
		    h_MR_R2_g1Mb0Wg1uW0Ll->Fill(MR, R2, w);
		    h_njets_g1Mb0Wg1uW0Ll->Fill(sjet.size(),w);
		    h_nbjets_g1Mb0Wg1uW0Ll->Fill(nmediumbs,w);
		    h_met_g1Mb0Wg1uW0Ll->Fill(met.Pt(),w);
		    h_jet1pt_g1Mb0Wg1uW0Ll->Fill(sjet[0].pt,w);
		    h_jet2pt_g1Mb0Wg1uW0Ll->Fill(sjet[1].pt,w);
		    h_jet3pt_g1Mb0Wg1uW0Ll->Fill(sjet[2].pt,w);
		    h_PV_g1Mb0Wg1uW0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
		    
		    if(isTTallhad)
		      TTallhad->Fill("g1Mb0Wg1uW0Ll", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("g1Mb0Wg1uW0Ll", w);
		    else if(isTTdilep)
		      TTdilep->Fill("g1Mb0Wg1uW0Ll", w);
		  } // end of aW.size() > 0
		} // end of sW.size() == 0
	      } // end of nmediumbs > 0
	      
	      if (nloosebs == 0){
		if (eventhelper_isRealData!=1) {
		  w = w*wCSVL;
		}
		ofile.count("0Lb0Ll", w);
		h_MR_0Lb0Ll->Fill(MR, w);
		h_R2_0Lb0Ll->Fill(R2, w);
		h_MR_R2_0Lb0Ll->Fill(MR, R2, w);

		h_minDeltaPhi_0Lb0Ll->Fill(minDeltaPhi, w);
		h_MR_minDeltaPhi_0Lb0Ll->Fill(MR, minDeltaPhi, w);
		h_R2_minDeltaPhi_0Lb0Ll->Fill(R2, minDeltaPhi, w);
		
		h_minDeltaPhiHat_0Lb0Ll->Fill(minDeltaPhiHat, w);
		h_MR_minDeltaPhiHat_0Lb0Ll->Fill(MR, minDeltaPhiHat, w);
		h_R2_minDeltaPhiHat_0Lb0Ll->Fill(R2, minDeltaPhiHat, w);

		h_PV_0Lb0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
		if(isTTallhad)
		  TTallhad->Fill("0Lb0Ll", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lb0Ll", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lb0Ll", w);
		
		// 0Lbg1uW0Ll -- QCD control region
		if( aW.size() > 0){
		  ofile.count("0Lbg1uW0Ll",w);
		  h_MR_0Lbg1uW0Ll->Fill(MR, w);
		  h_R2_0Lbg1uW0Ll->Fill(R2, w);
		  h_MR_R2_0Lbg1uW0Ll->Fill(MR, R2, w);
		  h_PV_0Lbg1uW0Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);

		  h_minDeltaPhi_0Lbg1uW0Ll->Fill(minDeltaPhi, w);
		  h_MR_minDeltaPhi_0Lbg1uW0Ll->Fill(MR, minDeltaPhi, w);
		  h_R2_minDeltaPhi_0Lbg1uW0Ll->Fill(R2, minDeltaPhi, w);

		  h_minDeltaPhiHat_0Lbg1uW0Ll->Fill(minDeltaPhiHat, w);
		  h_MR_minDeltaPhiHat_0Lbg1uW0Ll->Fill(MR, minDeltaPhiHat, w);
		  h_R2_minDeltaPhiHat_0Lbg1uW0Ll->Fill(R2, minDeltaPhiHat, w);

		  h_njets_0Lbg1uW0Ll->Fill(sjet.size(),w);
		  h_met_0Lbg1uW0Ll->Fill(met.Pt(),w);
		  h_jet1pt_0Lbg1uW0Ll->Fill(sjet[0].pt,w);
		  h_jet2pt_0Lbg1uW0Ll->Fill(sjet[1].pt,w);
		  h_jet3pt_0Lbg1uW0Ll->Fill(sjet[2].pt,w);
		  h_HT_0Lbg1uW0Ll->Fill(HT,w);

		  if(isTTallhad)
		    TTallhad->Fill("0Lbg1uW0Ll", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("0Lbg1uW0Ll", w);
		  else if(isTTdilep)
		    TTdilep->Fill("0Lbg1uW0Ll", w);
		  
		  // cut on mindDeltaPhi
		  if (minDeltaPhi < 0.3){
		    ofile.count("0Lbg1uW0Ll_mdPhi0p3",w);
		    h_MR_0Lbg1uW0Ll_mdPhi0p3->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhi0p3->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhi0p3->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p3->Fill(MR, R2, 1.);

		    h_njets_0Lbg1uW0Ll_mdPhi0p3->Fill(sjet.size(),w);
		    h_met_0Lbg1uW0Ll_mdPhi0p3->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1uW0Ll_mdPhi0p3->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1uW0Ll_mdPhi0p3->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1uW0Ll_mdPhi0p3->Fill(sjet[2].pt,w);

 		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		  } else {
		    h_MR_0Lbg1uW0Ll_mdPhig0p3->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhig0p3->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhig0p3->Fill(MR, R2, w);

 		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhig0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhig0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhig0p3", w);
		  } // end of minDeltaPhi < 0.3

		  // cut on mindDeltaPhi
		  if (minDeltaPhi < 0.5){
		    ofile.count("0Lbg1uW0Ll_mdPhi0p5",w);
		    h_MR_0Lbg1uW0Ll_mdPhi0p5->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhi0p5->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhi0p5->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1uW0Ll_mdPhi0p5->Fill(MR, R2, 1.);

		    h_njets_0Lbg1uW0Ll_mdPhi0p5->Fill(sjet.size(),w);
		    h_met_0Lbg1uW0Ll_mdPhi0p5->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1uW0Ll_mdPhi0p5->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1uW0Ll_mdPhi0p5->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1uW0Ll_mdPhi0p5->Fill(sjet[2].pt,w);

 		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p5", w);
		  } else {
		    h_MR_0Lbg1uW0Ll_mdPhig0p5->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhig0p5->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhig0p5->Fill(MR, R2, w);

 		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhig0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhig0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhig0p5", w);
		  } // end of minDeltaPhi < 0.5

		  if (minDeltaPhiHat < 4){
		    ofile.count("0Lbg1uW0Ll_mdPhiHat4",w);
		    h_MR_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1uW0Ll_mdPhiHat4->Fill(MR, R2, 1.);

		    h_njets_0Lbg1uW0Ll_mdPhiHat4->Fill(sjet.size(),w);
		    h_met_0Lbg1uW0Ll_mdPhiHat4->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1uW0Ll_mdPhiHat4->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1uW0Ll_mdPhiHat4->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1uW0Ll_mdPhiHat4->Fill(sjet[2].pt,w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat4", w);
		  } // end of minDeltaPhiHat < 4
		  else {
		    h_MR_0Lbg1uW0Ll_mdPhiHatg4->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhiHatg4->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhiHatg4->Fill(MR, R2, w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhiHatg4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHatg4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhiHatg4", w);
		  } // end of minDeltaPhiHat > 4

		  if (minDeltaPhiHat < 5){
		    ofile.count("0Lbg1uW0Ll_mdPhiHat5",w);
		    h_MR_0Lbg1uW0Ll_mdPhiHat5->Fill(MR, w);
		    h_R2_0Lbg1uW0Ll_mdPhiHat5->Fill(R2, w);
		    h_MR_R2_0Lbg1uW0Ll_mdPhiHat5->Fill(MR, R2, w);

		    h_njets_0Lbg1uW0Ll_mdPhiHat5->Fill(sjet.size(),w);
		    h_met_0Lbg1uW0Ll_mdPhiHat5->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1uW0Ll_mdPhiHat5->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1uW0Ll_mdPhiHat5->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1uW0Ll_mdPhiHat5->Fill(sjet[2].pt,w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1uW0Ll_mdPhiHat5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1uW0Ll_mdPhiHat5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1uW0Ll_mdPhiHat5", w);
		  } // end of minDeltaPhiHat < 5
		} // end of aW.size() > 0
		
		// 0Lbg1W0Ll
		if( sW.size() > 0){
		  ofile.count("0Lbg1W0Ll",w);
		  h_MR_0Lbg1W0Ll->Fill(MR, w);
		  h_R2_0Lbg1W0Ll->Fill(R2, w);
		  h_MR_R2_0Lbg1W0Ll->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("0Lbg1W0Ll", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("0Lbg1W0Ll", w);
		  else if(isTTdilep)
		    TTdilep->Fill("0Lbg1W0Ll", w);
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
	  
	  ofile.count("1Ll",w);
	  h_MR_1Ll->Fill(MR, w);
	  h_R2_1Ll->Fill(R2, w);
	  h_MR_R2_1Ll->Fill(MR, R2, w);
	  h_HT_1Ll->Fill(HT,w);
	  h_PV_1Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
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
	    h_MR_g1Mb1Ll->Fill(MR, w);
	    h_R2_g1Mb1Ll->Fill(R2, w);
	    h_MR_R2_g1Mb1Ll->Fill(MR, R2, w);
	    h_PV_g1Mb1Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);

	    h_mT_g1Mb1Ll->Fill(mT, w);
	    h_njets_g1Mb1Ll->Fill(sjet.size(),w);
	    h_nbjets_g1Mb1Ll->Fill(nmediumbs,w);
	    h_met_g1Mb1Ll->Fill(met.Pt(),w);
	    h_jet1pt_g1Mb1Ll->Fill(sjet[0].pt,w);
	    h_jet2pt_g1Mb1Ll->Fill(sjet[1].pt,w);
	    h_jet3pt_g1Mb1Ll->Fill(sjet[2].pt,w);
	    h_leptonpt_g1Mb1Ll->Fill(lepton.Pt(),w);
	    
	    if(isTTallhad)
	      TTallhad->Fill("g1Mb1Ll", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("g1Mb1Ll", w);
	    else if(isTTdilep)
	      TTdilep->Fill("g1Mb1Ll", w);
	    
	    if( sW.size() > 0 ){
	      ofile.count("g1Mbg1W1Ll",w);
	      h_MR_g1Mbg1W1Ll->Fill(MR, w);
	      h_R2_g1Mbg1W1Ll->Fill(R2, w);
	      h_MR_R2_g1Mbg1W1Ll->Fill(MR, R2, w);

	      h_mT_g1Mbg1W1Ll->Fill(mT, w);
	      h_MR_mT_g1Mbg1W1Ll->Fill(MR, mT, w);
	      h_R2_mT_g1Mbg1W1Ll->Fill(R2, mT, w);

	      h_njets_g1Mbg1W1Ll->Fill(sjet.size(),w);
	      h_nbjets_g1Mbg1W1Ll->Fill(nmediumbs,w);
	      h_met_g1Mbg1W1Ll->Fill(met.Pt(),w);
	      h_jet1pt_g1Mbg1W1Ll->Fill(sjet[0].pt,w);
	      h_jet2pt_g1Mbg1W1Ll->Fill(sjet[1].pt,w);
	      h_jet3pt_g1Mbg1W1Ll->Fill(sjet[2].pt,w);
	      h_leptonpt_g1Mbg1W1Ll->Fill(lepton.Pt(),w);
	      
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

		h_njets_g1Mbg1W1LlmT100->Fill(sjet.size(),w);
		h_nbjets_g1Mbg1W1LlmT100->Fill(nmediumbs,w);
		h_met_g1Mbg1W1LlmT100->Fill(met.Pt(),w);
		h_jet1pt_g1Mbg1W1LlmT100->Fill(sjet[0].pt,w);
		h_jet2pt_g1Mbg1W1LlmT100->Fill(sjet[1].pt,w);
		h_jet3pt_g1Mbg1W1LlmT100->Fill(sjet[2].pt,w);
		h_leptonpt_g1Mbg1W1LlmT100->Fill(lepton.Pt(),w);
		h_HT_g1Mbg1W1LlmT100->Fill(HT,w);
		h_minDeltaPhiHat_g1Mbg1W1LlmT100->Fill(minDeltaPhiHat, w);
		h_MR_minDeltaPhiHat_g1Mbg1W1LlmT100->Fill(MR, minDeltaPhiHat, w);
		h_R2_minDeltaPhiHat_g1Mbg1W1LlmT100->Fill(R2, minDeltaPhiHat, w);
		h_minDeltaPhi_g1Mbg1W1LlmT100->Fill(minDeltaPhi, w);
		h_MR_minDeltaPhi_g1Mbg1W1LlmT100->Fill(MR, minDeltaPhi, w);
		h_R2_minDeltaPhi_g1Mbg1W1LlmT100->Fill(R2, minDeltaPhi, w);
		h_PV_g1Mbg1W1LlmT100->Fill(eventhelperextra_numberOfPrimaryVertices,w_old*w_trigger*w_TopPt_nominal*w_ISR_nominal);
		h_PV_reweighted_g1Mbg1W1LlmT100->Fill(eventhelperextra_numberOfPrimaryVertices,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1W1LlmT100", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1W1LlmT100", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1W1LlmT100", w);

		// separate 1b and >=2b regions
		if (nmediumbs == 1){
		  ofile.count("1Mbg1W1LlmT100",w);
		  h_MR_1Mbg1W1LlmT100->Fill(MR, w);
		  h_R2_1Mbg1W1LlmT100->Fill(R2, w);
		  h_MR_R2_1Mbg1W1LlmT100->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("1Mbg1W1LlmT100", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("1Mbg1W1LlmT100", w);
		  else if(isTTdilep)
		    TTdilep->Fill("1Mbg1W1LlmT100", w);
		} else {
		  ofile.count("g2Mbg1W1LlmT100",w);
		  h_MR_g2Mbg1W1LlmT100->Fill(MR, w);
		  h_R2_g2Mbg1W1LlmT100->Fill(R2, w);
		  h_MR_R2_g2Mbg1W1LlmT100->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("g2Mbg1W1LlmT100", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g2Mbg1W1LlmT100", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g2Mbg1W1LlmT100", w);
		}
		
		// add cut on mindeltaphihat
		if (minDeltaPhiHat > 4){
		  ofile.count("g1Mbg1W1LlmT100_mdPhiHatg4",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(MR, R2, w);
		  h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(MR, R2, 1.);

		  h_njets_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(sjet.size(),w);
		  h_nbjets_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(nmediumbs,w);
		  h_met_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(met.Pt(),w);
		  h_jet1pt_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(sjet[0].pt,w);
		  h_jet2pt_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(sjet[1].pt,w);
		  h_jet3pt_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(sjet[2].pt,w);
		  h_leptonpt_g1Mbg1W1LlmT100_mdPhiHatg4->Fill(lepton.Pt(),w);

		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhiHatg4", w);
		} else {
		  ofile.count("g1Mbg1W1LlmT100_mdPhiHat4",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhiHat4->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhiHat4->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhiHat4->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhiHat4", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhiHat4", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhiHat4", w);
		}

		// add cut on mindeltaphi
		if (minDeltaPhi > 0.3){
		  ofile.count("g1Mbg1W1LlmT100_mdPhig0p3",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhig0p3->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhig0p3->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3->Fill(MR, R2, w);
		  h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p3->Fill(MR, R2, 1.);

		  h_njets_g1Mbg1W1LlmT100_mdPhig0p3->Fill(sjet.size(),w);
		  h_nbjets_g1Mbg1W1LlmT100_mdPhig0p3->Fill(nmediumbs,w);
		  h_met_g1Mbg1W1LlmT100_mdPhig0p3->Fill(met.Pt(),w);
		  h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p3->Fill(sjet[0].pt,w);
		  h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p3->Fill(sjet[1].pt,w);
		  h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p3->Fill(sjet[2].pt,w);
		  h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p3->Fill(lepton.Pt(),w);

		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhig0p3", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhig0p3", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhig0p3", w);
		} else {
		  ofile.count("g1Mbg1W1LlmT100_mdPhi0p3",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhi0p3->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhi0p3->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p3->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhi0p3", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhi0p3", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhi0p3", w);
		}

		// add cut on mindeltaphi
		if (minDeltaPhi > 0.5){
		  ofile.count("g1Mbg1W1LlmT100_mdPhig0p5",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhig0p5->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhig0p5->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5->Fill(MR, R2, w);
		  h_uw_MR_R2_g1Mbg1W1LlmT100_mdPhig0p5->Fill(MR, R2, 1.);

		  h_njets_g1Mbg1W1LlmT100_mdPhig0p5->Fill(sjet.size(),w);
		  h_nbjets_g1Mbg1W1LlmT100_mdPhig0p5->Fill(nmediumbs,w);
		  h_met_g1Mbg1W1LlmT100_mdPhig0p5->Fill(met.Pt(),w);
		  h_jet1pt_g1Mbg1W1LlmT100_mdPhig0p5->Fill(sjet[0].pt,w);
		  h_jet2pt_g1Mbg1W1LlmT100_mdPhig0p5->Fill(sjet[1].pt,w);
		  h_jet3pt_g1Mbg1W1LlmT100_mdPhig0p5->Fill(sjet[2].pt,w);
		  h_leptonpt_g1Mbg1W1LlmT100_mdPhig0p5->Fill(lepton.Pt(),w);

		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhig0p5", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhig0p5", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhig0p5", w);
		} else {
		  ofile.count("g1Mbg1W1LlmT100_mdPhi0p5",w);
		  h_MR_g1Mbg1W1LlmT100_mdPhi0p5->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT100_mdPhi0p5->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT100_mdPhi0p5->Fill(MR, R2, w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT100_mdPhi0p5", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT100_mdPhi0p5", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT100_mdPhi0p5", w);
		}

		// mT window
		if (mT > 30){
		  ofile.count("g1Mbg1W1LlmT",w);
		  h_MR_g1Mbg1W1LlmT->Fill(MR, w);
		  h_R2_g1Mbg1W1LlmT->Fill(R2, w);
		  h_MR_R2_g1Mbg1W1LlmT->Fill(MR, R2, w);
		  h_mT_g1Mbg1W1LlmT->Fill(mT, w);
		  h_MR_mT_g1Mbg1W1LlmT->Fill(MR, mT, w);
		  h_R2_mT_g1Mbg1W1LlmT->Fill(R2, mT, w);
		  if(isTTallhad)
		    TTallhad->Fill("g1Mbg1W1LlmT", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("g1Mbg1W1LlmT", w);
		  else if(isTTdilep)
		    TTdilep->Fill("g1Mbg1W1LlmT", w);
		} // end mT > 30
	      } // end mT < 100
	    } // end sW.size()
	  } // end nmediumbs > 0


	  if (nloosebs == 0){
	    if (eventhelper_isRealData!=1) {
	      w = w*wCSVL;
	    }
	    ofile.count("0Lb1Ll",w);
	    h_MR_0Lb1Ll->Fill(MR, w);
	    h_R2_0Lb1Ll->Fill(R2, w);
	    h_MR_R2_0Lb1Ll->Fill(MR, R2, w);	
	    h_PV_0Lb1Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
    
	    h_mT_0Lb1Ll->Fill(mT, w);
	    h_MR_mT_0Lb1Ll->Fill(MR, mT, w);
	    h_R2_mT_0Lb1Ll->Fill(R2, mT, w);

	    h_met_0Lb1Ll->Fill(met.Pt(),w);
	    h_njets_0Lb1Ll->Fill(sjet.size(),w);
	    h_jet1pt_0Lb1Ll->Fill(sjet[0].pt,w);
	    h_jet2pt_0Lb1Ll->Fill(sjet[1].pt,w);
	    h_jet3pt_0Lb1Ll->Fill(sjet[2].pt,w);
	    h_leptonpt_0Lb1Ll->Fill(lepton.Pt(),w);

	    if(isTTallhad)
	      TTallhad->Fill("0Lb1Ll", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("0Lb1Ll", w);
	    else if(isTTdilep)
	      TTdilep->Fill("0Lb1Ll", w);
	    
	    if( sY.size() > 0 ){
	      ofile.count("0Lbg1Y1Ll",w);
	      h_MR_0Lbg1Y1Ll->Fill(MR, w);
	      h_R2_0Lbg1Y1Ll->Fill(R2, w);
	      h_MR_R2_0Lbg1Y1Ll->Fill(MR, R2, w);
	      h_PV_0Lbg1Y1Ll->Fill(eventhelperextra_numberOfPrimaryVertices,w);
	      
	      h_mT_0Lbg1Y1Ll->Fill(mT, w);
	      h_MR_mT_0Lbg1Y1Ll->Fill(MR, mT, w);
	      h_R2_mT_0Lbg1Y1Ll->Fill(R2, mT, w);
	      
	      h_met_0Lbg1Y1Ll->Fill(met.Pt(),w);
	      h_njets_0Lbg1Y1Ll->Fill(sjet.size(),w);
	      h_jet1pt_0Lbg1Y1Ll->Fill(sjet[0].pt,w);
	      h_jet2pt_0Lbg1Y1Ll->Fill(sjet[1].pt,w);
	      h_jet3pt_0Lbg1Y1Ll->Fill(sjet[2].pt,w);
	      h_leptonpt_0Lbg1Y1Ll->Fill(lepton.Pt(),w);

	      if(isTTallhad)
		TTallhad->Fill("0Lbg1Y1Ll", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lbg1Y1Ll", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lbg1Y1Ll", w);
	      
	      // WJets Control Region
	      if (mT < 100){ 
		ofile.count("0Lbg1Y1LlmT100",w);
		h_MR_0Lbg1Y1LlmT100->Fill(MR, w);
		h_R2_0Lbg1Y1LlmT100->Fill(R2, w);
		h_MR_R2_0Lbg1Y1LlmT100->Fill(MR, R2, w);

		h_njets_0Lbg1Y1LlmT100->Fill(sjet.size(),w);
		h_met_0Lbg1Y1LlmT100->Fill(met.Pt(),w);
		h_jet1pt_0Lbg1Y1LlmT100->Fill(sjet[0].pt,w);
		h_jet2pt_0Lbg1Y1LlmT100->Fill(sjet[1].pt,w);
		h_jet3pt_0Lbg1Y1LlmT100->Fill(sjet[2].pt,w);
		h_leptonpt_0Lbg1Y1LlmT100->Fill(lepton.Pt(),w);

		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y1LlmT100", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y1LlmT100", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y1LlmT100", w);
		
		// mT window
		if (mT > 30){
		  ofile.count("0Lbg1Y1LlmT",w);
		  h_MR_0Lbg1Y1LlmT->Fill(MR, w);
		  h_R2_0Lbg1Y1LlmT->Fill(R2, w);
		  h_MR_R2_0Lbg1Y1LlmT->Fill(MR, R2, w);
		  h_uw_MR_R2_0Lbg1Y1LlmT->Fill(MR, R2, 1.);

		  h_mT_0Lbg1Y1LlmT->Fill(mT, w);
		  h_MR_mT_0Lbg1Y1LlmT->Fill(MR, mT, w);
		  h_R2_mT_0Lbg1Y1LlmT->Fill(R2, mT, w);

		  h_njets_0Lbg1Y1LlmT->Fill(sjet.size(),w);
		  h_met_0Lbg1Y1LlmT->Fill(met.Pt(),w);
		  h_jet1pt_0Lbg1Y1LlmT->Fill(sjet[0].pt,w);
		  h_jet2pt_0Lbg1Y1LlmT->Fill(sjet[1].pt,w);
		  h_jet3pt_0Lbg1Y1LlmT->Fill(sjet[2].pt,w);
		  h_leptonpt_0Lbg1Y1LlmT->Fill(lepton.Pt(),w);
		  h_HT_0Lbg1Y1LlmT->Fill(HT,w);
		  h_minDeltaPhiHat_0Lbg1Y1LlmT->Fill(minDeltaPhiHat, w);
		  h_MR_minDeltaPhiHat_0Lbg1Y1LlmT->Fill(MR, minDeltaPhiHat, w);
		  h_R2_minDeltaPhiHat_0Lbg1Y1LlmT->Fill(R2, minDeltaPhiHat, w);
		  h_minDeltaPhi_0Lbg1Y1LlmT->Fill(minDeltaPhi, w);
		  h_MR_minDeltaPhi_0Lbg1Y1LlmT->Fill(MR, minDeltaPhi, w);
		  h_R2_minDeltaPhi_0Lbg1Y1LlmT->Fill(R2, minDeltaPhi, w);
		  if(isTTallhad)
		    TTallhad->Fill("0Lbg1Y1LlmT", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("0Lbg1Y1LlmT", w);
		  else if(isTTdilep)
		    TTdilep->Fill("0Lbg1Y1LlmT", w);

		  // add cut on mindeltaphihat
		  if (minDeltaPhiHat > 4){
		    ofile.count("0Lbg1Y1LlmT_mdPhiHatg4",w);
		    h_MR_0Lbg1Y1LlmT_mdPhiHatg4->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhiHatg4->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1Y1LlmT_mdPhiHatg4->Fill(MR, R2, 1.);

		    h_njets_0Lbg1Y1LlmT_mdPhiHatg4->Fill(sjet.size(),w);
		    h_met_0Lbg1Y1LlmT_mdPhiHatg4->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1Y1LlmT_mdPhiHatg4->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1Y1LlmT_mdPhiHatg4->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1Y1LlmT_mdPhiHatg4->Fill(sjet[2].pt,w);
		    h_leptonpt_0Lbg1Y1LlmT_mdPhiHatg4->Fill(lepton.Pt(),w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhiHatg4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhiHatg4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhiHatg4", w);
		  } else {
		    ofile.count("0Lbg1Y1LlmT_mdPhiHat4",w);
		    h_MR_0Lbg1Y1LlmT_mdPhiHat4->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhiHat4->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhiHat4->Fill(MR, R2, w);
		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhiHat4", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhiHat4", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhiHat4", w);
		  }

		  // add cut on mindeltaphi
		  if (minDeltaPhi > 0.3){
		    ofile.count("0Lbg1Y1LlmT_mdPhig0p3",w);
		    h_MR_0Lbg1Y1LlmT_mdPhig0p3->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhig0p3->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhig0p3->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p3->Fill(MR, R2, 1.);

		    h_njets_0Lbg1Y1LlmT_mdPhig0p3->Fill(sjet.size(),w);
		    h_met_0Lbg1Y1LlmT_mdPhig0p3->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1Y1LlmT_mdPhig0p3->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1Y1LlmT_mdPhig0p3->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1Y1LlmT_mdPhig0p3->Fill(sjet[2].pt,w);
		    h_leptonpt_0Lbg1Y1LlmT_mdPhig0p3->Fill(lepton.Pt(),w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhig0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhig0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhig0p3", w);
		  } else {
		    ofile.count("0Lbg1Y1LlmT_mdPhi0p3",w);
		    h_MR_0Lbg1Y1LlmT_mdPhi0p3->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhi0p3->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhi0p3->Fill(MR, R2, w);
		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhi0p3", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhi0p3", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhi0p3", w);
		  }

		  // add cut on mindeltaphi
		  if (minDeltaPhi > 0.5){
		    ofile.count("0Lbg1Y1LlmT_mdPhig0p5",w);
		    h_MR_0Lbg1Y1LlmT_mdPhig0p5->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhig0p5->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhig0p5->Fill(MR, R2, w);
		    h_uw_MR_R2_0Lbg1Y1LlmT_mdPhig0p5->Fill(MR, R2, 1.);

		    h_njets_0Lbg1Y1LlmT_mdPhig0p5->Fill(sjet.size(),w);
		    h_met_0Lbg1Y1LlmT_mdPhig0p5->Fill(met.Pt(),w);
		    h_jet1pt_0Lbg1Y1LlmT_mdPhig0p5->Fill(sjet[0].pt,w);
		    h_jet2pt_0Lbg1Y1LlmT_mdPhig0p5->Fill(sjet[1].pt,w);
		    h_jet3pt_0Lbg1Y1LlmT_mdPhig0p5->Fill(sjet[2].pt,w);
		    h_leptonpt_0Lbg1Y1LlmT_mdPhig0p5->Fill(lepton.Pt(),w);

		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhig0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhig0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhig0p5", w);
		  } else {
		    ofile.count("0Lbg1Y1LlmT_mdPhi0p5",w);
		    h_MR_0Lbg1Y1LlmT_mdPhi0p5->Fill(MR, w);
		    h_R2_0Lbg1Y1LlmT_mdPhi0p5->Fill(R2, w);
		    h_MR_R2_0Lbg1Y1LlmT_mdPhi0p5->Fill(MR, R2, w);
		    if(isTTallhad)
		      TTallhad->Fill("0Lbg1Y1LlmT_mdPhi0p5", w);
		    else if(isTTsemilep)
		      TTsemilep->Fill("0Lbg1Y1LlmT_mdPhi0p5", w);
		    else if(isTTdilep)
		      TTdilep->Fill("0Lbg1Y1LlmT_mdPhi0p5", w);
		  }

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

      // dimuon selections for the muon acceptance and efficiency - after preselection
      if (genmu.size() == 2 && genZ.size() == 1) {
	if (genmu[0].charge == -genmu[1].charge){
	  ofile.count("presel2genallmu", w);
	  if(isTTallhad)
	    TTallhad->Fill("presel2genallmu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("presel2genallmu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("presel2genallmu", w);
	}
      }
      if (sgenmu.size() == 2) {
	if (genmu[0].charge == -genmu[1].charge){
	  TLorentzVector LVZcand;
	  for (unsigned int m=0; m<LVsgenmu.size(); m++) {
	    LVZcand += LVsgenmu[m];
	  }
	  double Zmass = LVZcand.M();
	  if (!(Zmass >= 60 && Zmass <= 120)) continue;
	  ofile.count("presel2genmu", w);
	  if(isTTallhad)
	    TTallhad->Fill("presel2genmu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("presel2genmu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("presel2genmu", w);
	}
      }
      if (smuon.size() == 2) {
	TLorentzVector LVZcand;
	for (unsigned int m=0; m<LVmu.size(); m++) {
	  LVZcand += LVmu[m];
	}
	double Zmass = LVZcand.M();
	if (Zmass >= 60 && Zmass <= 120){
	  ofile.count("presel2mu", w);
	  if(isTTallhad)
	    TTallhad->Fill("presel2mu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("presel2mu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("presel2mu", w);
	}
      }


      // Need to use R2metmu to define the signal region
      if (MR > 800 && R2metmu > 0.08){
	// GEN
	if (genmu.size() == 2 && genZ.size() == 1) {
	  if (genmu[0].charge == -genmu[1].charge){
	    ofile.count("2genallmu", w);
	    if(isTTallhad)
	      TTallhad->Fill("2genallmu", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2genallmu", w);
	    else if(isTTdilep)
	      TTdilep->Fill("2genallmu", w);
	  }
	}
	if (sgenmu.size() == 2) {
	  if (genmu[0].charge == -genmu[1].charge) {
	    TLorentzVector LVZcand;
	    for (unsigned int m=0; m<LVsgenmu.size(); m++) {
	      LVZcand += LVsgenmu[m];
	    }
	    double Zmass = LVZcand.M();
	    if (!(Zmass >= 60 && Zmass <= 120)) continue;
	    ofile.count("2genmu", w);
	    if(isTTallhad)
	      TTallhad->Fill("2genmu", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2genmu", w);
	    else if(isTTdilep)
	      TTdilep->Fill("genmu", w);
	  }
	}

	// reco
	TLorentzVector LVZcand;
	if (ntightmuons == 2 && nloosemuons == 2) {
	  // Make sure that the muons are opposite-signed:
	  if (!(smuon[0].charge == -smuon[1].charge)) continue;
	  // See if the 2mu makes a Z:
	  for (unsigned int m=0; m<LVmu.size(); m++) {
	    LVZcand += LVmu[m];
	  }
	  double Zmass = LVZcand.M();
	  h_Zmass_2mu->Fill(Zmass, w);
	  h_MR_Zmass_2mu->Fill(MR, Zmass, w);
	  h_R2_Zmass_2mu->Fill(R2metmu, Zmass, w);
	  h_Zmass_2l->Fill(Zmass, w);
	  h_MR_Zmass_2l->Fill(MR, Zmass, w);
	  h_R2_Zmass_2l->Fill(R2metmu, Zmass, w);

	  ofile.count("2munoZmass",w);
	  h_MR_2munoZmass->Fill(MR, w);
	  h_R2_2munoZmass->Fill(R2metmu, w);
	  h_MR_R2_2munoZmass->Fill(MR, R2metmu, w);
	  if(isTTallhad)
	    TTallhad->Fill("2munoZmass", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2munoZmass", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2munoZmass", w);

	  ofile.count("2lnoZmass",w);
	  h_MR_2lnoZmass->Fill(MR, w);
	  h_R2_2lnoZmass->Fill(R2metmu, w);
	  h_MR_R2_2lnoZmass->Fill(MR, R2metmu, w);
	  if(isTTallhad)
	    TTallhad->Fill("2lnoZmass", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2lnoZmass", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2lnoZmass", w);
	  
	  if (!(Zmass >= 60 && Zmass <= 120)) continue;
	  ofile.count("2mu",w);
	  h_MR_2mu->Fill(MR, w);
	  h_R2_2mu->Fill(R2metmu, w);
	  h_MR_R2_2mu->Fill(MR, R2metmu, w);
	  h_lepton1pt_2mu->Fill(smuon[0].pt,w);
	  h_lepton2pt_2mu->Fill(smuon[1].pt,w);
	  if(isTTallhad)
	    TTallhad->Fill("2mu", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2mu", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2mu", w);

	  ofile.count("2l",w);
	  h_MR_2l->Fill(MR, w);
	  h_R2_2l->Fill(R2metmu, w);
	  h_MR_R2_2l->Fill(MR, R2metmu, w);
	  h_lepton1pt_2l->Fill(smuon[0].pt,w);
	  h_lepton2pt_2l->Fill(smuon[1].pt,w);
	  if(isTTallhad)
	    TTallhad->Fill("2l", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2l", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2l", w);
	  
	  if (nlooseelectrons == 0){
	    ofile.count("2mu0el",w);
	    h_MR_2mu0el->Fill(MR, w);
	    h_R2_2mu0el->Fill(R2metmu, w);
	    h_MR_R2_2mu0el->Fill(MR, R2metmu, w);
	    h_njets_2mu0el->Fill(sjet.size(),w);
	    h_nbjets_2mu0el->Fill(nmediumbs,w);
	    h_met_2mu0el->Fill(V3metmu.Pt(),w);
	    h_jet1pt_2mu0el->Fill(sjet[0].pt,w);
	    h_jet2pt_2mu0el->Fill(sjet[1].pt,w);
	    h_jet3pt_2mu0el->Fill(sjet[2].pt,w);
	    h_lepton1pt_2mu0el->Fill(smuon[0].pt,w);
	    h_lepton2pt_2mu0el->Fill(smuon[1].pt,w);
	    h_HT_2mu0el->Fill(HT,w);
	    h_PV_2mu0el->Fill(eventhelperextra_numberOfPrimaryVertices,w);
	    if(isTTallhad)
	      TTallhad->Fill("2mu0el", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2mu0el", w);
	    else if(isTTdilep)
	      TTdilep->Fill("2mu0el", w);

	    ofile.count("2l0ol",w);
	    h_MR_2l0ol->Fill(MR, w);
	    h_R2_2l0ol->Fill(R2metmu, w);
	    h_MR_R2_2l0ol->Fill(MR, R2metmu, w);
	    h_njets_2l0ol->Fill(sjet.size(),w);
	    h_nbjets_2l0ol->Fill(nmediumbs,w);
	    h_met_2l0ol->Fill(V3metmu.Pt(),w);
	    h_jet1pt_2l0ol->Fill(sjet[0].pt,w);
	    h_jet2pt_2l0ol->Fill(sjet[1].pt,w);
	    h_jet3pt_2l0ol->Fill(sjet[2].pt,w);
	    h_lepton1pt_2l0ol->Fill(smuon[0].pt,w);
	    h_lepton2pt_2l0ol->Fill(smuon[1].pt,w);
	    h_HT_2l0ol->Fill(HT,w);
	    if(isTTallhad)
	      TTallhad->Fill("2l0ol", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2l0ol", w);
	    else if(isTTdilep)
	      TTdilep->Fill("2l0ol", w);
	    
	    if (nloosebs == 0){
	      if (eventhelper_isRealData!=1) {
		w = w*wCSVL;
	      }
	      ofile.count("0Lb2mu0el",w);
	      h_MR_0Lb2mu0el->Fill(MR, w);
	      h_R2_0Lb2mu0el->Fill(R2metmu, w);
	      h_MR_R2_0Lb2mu0el->Fill(MR, R2metmu, w);
	      h_njets_0Lb2mu0el->Fill(sjet.size(),w);
	      h_met_0Lb2mu0el->Fill(V3metmu.Pt(),w);
	      h_jet1pt_0Lb2mu0el->Fill(sjet[0].pt,w);
	      h_jet2pt_0Lb2mu0el->Fill(sjet[1].pt,w);
	      h_jet3pt_0Lb2mu0el->Fill(sjet[2].pt,w);
	      h_lepton1pt_0Lb2mu0el->Fill(smuon[0].pt,w);
	      h_lepton2pt_0Lb2mu0el->Fill(smuon[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("0Lb2mu0el", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lb2mu0el", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lb2mu0el", w);

	      ofile.count("0Lb2l0ol",w);
	      h_MR_0Lb2l0ol->Fill(MR, w);
	      h_R2_0Lb2l0ol->Fill(R2metmu, w);
	      h_MR_R2_0Lb2l0ol->Fill(MR, R2metmu, w);
	      h_njets_0Lb2l0ol->Fill(sjet.size(),w);
	      h_met_0Lb2l0ol->Fill(V3metmu.Pt(),w);
	      h_jet1pt_0Lb2l0ol->Fill(sjet[0].pt,w);
	      h_jet2pt_0Lb2l0ol->Fill(sjet[1].pt,w);
	      h_jet3pt_0Lb2l0ol->Fill(sjet[2].pt,w);
	      h_lepton1pt_0Lb2l0ol->Fill(smuon[0].pt,w);
	      h_lepton2pt_0Lb2l0ol->Fill(smuon[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("0Lb2l0ol", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lb2l0ol", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lb2l0ol", w);

	      if (nYs >= 1){ // Z no b, mu CR 
		ofile.count("0Lbg1Y2mu0el",w);
		h_MR_0Lbg1Y2mu0el->Fill(MR, w);
		h_R2_0Lbg1Y2mu0el->Fill(R2metmu, w);
		h_MR_R2_0Lbg1Y2mu0el->Fill(MR, R2metmu, w);
		h_njets_0Lbg1Y2mu0el->Fill(sjet.size(),w);
		h_met_0Lbg1Y2mu0el->Fill(V3metmu.Pt(),w);
		h_jet1pt_0Lbg1Y2mu0el->Fill(sjet[0].pt,w);
		h_jet2pt_0Lbg1Y2mu0el->Fill(sjet[1].pt,w);
		h_jet3pt_0Lbg1Y2mu0el->Fill(sjet[2].pt,w);
		h_lepton1pt_0Lbg1Y2mu0el->Fill(smuon[0].pt,w);
		h_lepton2pt_0Lbg1Y2mu0el->Fill(smuon[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y2mu0el", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y2mu0el", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y2mu0el", w);

		ofile.count("0Lbg1Y2l0ol",w);
		h_MR_0Lbg1Y2l0ol->Fill(MR, w);
		h_R2_0Lbg1Y2l0ol->Fill(R2metmu, w);
		h_MR_R2_0Lbg1Y2l0ol->Fill(MR, R2metmu, w);
		h_njets_0Lbg1Y2l0ol->Fill(sjet.size(),w);
		h_met_0Lbg1Y2l0ol->Fill(V3metmu.Pt(),w);
		h_jet1pt_0Lbg1Y2l0ol->Fill(sjet[0].pt,w);
		h_jet2pt_0Lbg1Y2l0ol->Fill(sjet[1].pt,w);
		h_jet3pt_0Lbg1Y2l0ol->Fill(sjet[2].pt,w);
		h_lepton1pt_0Lbg1Y2l0ol->Fill(smuon[0].pt,w);
		h_lepton2pt_0Lbg1Y2l0ol->Fill(smuon[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y2l0ol", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y2l0ol", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y2l0ol", w);
	      }// nYs >= 1
	    }// end nloosebs == 0

	    if (nmediumbs >= 1){
	      if (eventhelper_isRealData!=1) {
		w = w*wCSVM;
	      }
	      ofile.count("g1Mb2mu0el",w);
	      h_MR_g1Mb2mu0el->Fill(MR, w);
	      h_R2_g1Mb2mu0el->Fill(R2metmu, w);
	      h_MR_R2_g1Mb2mu0el->Fill(MR, R2metmu, w);
	      h_njets_g1Mb2mu0el->Fill(sjet.size(),w);
	      h_nbjets_g1Mb2mu0el->Fill(nmediumbs,w);
	      h_met_g1Mb2mu0el->Fill(V3metmu.Pt(),w);
	      h_jet1pt_g1Mb2mu0el->Fill(sjet[0].pt,w);
	      h_jet2pt_g1Mb2mu0el->Fill(sjet[1].pt,w);
	      h_jet3pt_g1Mb2mu0el->Fill(sjet[2].pt,w);
	      h_lepton1pt_g1Mb2mu0el->Fill(smuon[0].pt,w);
	      h_lepton2pt_g1Mb2mu0el->Fill(smuon[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mb2mu0el", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mb2mu0el", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mb2mu0el", w);

	      ofile.count("g1Mb2l0ol",w);
	      h_MR_g1Mb2l0ol->Fill(MR, w);
	      h_R2_g1Mb2l0ol->Fill(R2metmu, w);
	      h_MR_R2_g1Mb2l0ol->Fill(MR, R2metmu, w);
	      h_njets_g1Mb2l0ol->Fill(sjet.size(),w);
	      h_nbjets_g1Mb2l0ol->Fill(nmediumbs,w);
	      h_met_g1Mb2l0ol->Fill(V3metmu.Pt(),w);
	      h_jet1pt_g1Mb2l0ol->Fill(sjet[0].pt,w);
	      h_jet2pt_g1Mb2l0ol->Fill(sjet[1].pt,w);
	      h_jet3pt_g1Mb2l0ol->Fill(sjet[2].pt,w);
	      h_lepton1pt_g1Mb2l0ol->Fill(smuon[0].pt,w);
	      h_lepton2pt_g1Mb2l0ol->Fill(smuon[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mb2l0ol", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mb2l0ol", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mb2l0ol", w);
	    
	      if (nYs >= 1){ // Z with b, mu CR
		ofile.count("g1Mbg1Y2mu0el",w);
		h_MR_g1Mbg1Y2mu0el->Fill(MR, w);
		h_R2_g1Mbg1Y2mu0el->Fill(R2metmu, w);
		h_MR_R2_g1Mbg1Y2mu0el->Fill(MR, R2metmu, w);
		h_njets_g1Mbg1Y2mu0el->Fill(sjet.size(),w);
		h_nbjets_g1Mbg1Y2mu0el->Fill(nmediumbs,w);
		h_met_g1Mbg1Y2mu0el->Fill(V3metmu.Pt(),w);
		h_jet1pt_g1Mbg1Y2mu0el->Fill(sjet[0].pt,w);
		h_jet2pt_g1Mbg1Y2mu0el->Fill(sjet[1].pt,w);
		h_jet3pt_g1Mbg1Y2mu0el->Fill(sjet[2].pt,w);
		h_lepton1pt_g1Mbg1Y2mu0el->Fill(smuon[0].pt,w);
		h_lepton2pt_g1Mbg1Y2mu0el->Fill(smuon[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1Y2mu0el", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1Y2mu0el", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1Y2mu0el", w);

		ofile.count("g1Mbg1Y2l0ol",w);
		h_MR_g1Mbg1Y2l0ol->Fill(MR, w);
		h_R2_g1Mbg1Y2l0ol->Fill(R2metmu, w);
		h_MR_R2_g1Mbg1Y2l0ol->Fill(MR, R2metmu, w);
		h_njets_g1Mbg1Y2l0ol->Fill(sjet.size(),w);
		h_nbjets_g1Mbg1Y2l0ol->Fill(nmediumbs,w);
		h_met_g1Mbg1Y2l0ol->Fill(V3metmu.Pt(),w);
		h_jet1pt_g1Mbg1Y2l0ol->Fill(sjet[0].pt,w);
		h_jet2pt_g1Mbg1Y2l0ol->Fill(sjet[1].pt,w);
		h_jet3pt_g1Mbg1Y2l0ol->Fill(sjet[2].pt,w);
		h_lepton1pt_g1Mbg1Y2l0ol->Fill(smuon[0].pt,w);
		h_lepton2pt_g1Mbg1Y2l0ol->Fill(smuon[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1Y2l0ol", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1Y2l0ol", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1Y2l0ol", w);
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
	  h_Zmass_2el->Fill(Zmass, w);
	  h_MR_Zmass_2el->Fill(MR, Zmass, w);
	  h_R2_Zmass_2el->Fill(R2metel, Zmass, w);
	  h_Zmass_2l->Fill(Zmass, w);
	  h_MR_Zmass_2l->Fill(MR, Zmass, w);
	  h_R2_Zmass_2l->Fill(R2metel, Zmass, w);

	  ofile.count("2elnoZmass",w);
	  h_MR_2elnoZmass->Fill(MR, w);
	  h_R2_2elnoZmass->Fill(R2metel, w);
	  h_MR_R2_2elnoZmass->Fill(MR, R2metel, w);
	  if(isTTallhad)
	    TTallhad->Fill("2elnoZmass", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2elnoZmass", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2elnoZmass", w);
	  ofile.count("2lnoZmass",w);
	  h_MR_2lnoZmass->Fill(MR, w);
	  h_R2_2lnoZmass->Fill(R2metel, w);
	  h_MR_R2_2lnoZmass->Fill(MR, R2metel, w);
	  if(isTTallhad)
	    TTallhad->Fill("2lnoZmass", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2lnoZmass", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2lnoZmass", w);

     	  if (!(Zmass >= 60 && Zmass <= 120)) continue;
	  ofile.count("2el",w);
	  h_MR_2el->Fill(MR, w);
	  h_R2_2el->Fill(R2metel, w);
	  h_MR_R2_2el->Fill(MR, R2metel, w);
	  h_lepton1pt_2el->Fill(selectron[0].pt,w);
	  h_lepton2pt_2el->Fill(selectron[1].pt,w);
	  if(isTTallhad)
	    TTallhad->Fill("2el", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2el", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2el", w);
	  ofile.count("2l",w);
	  h_MR_2l->Fill(MR, w);
	  h_R2_2l->Fill(R2metel, w);
	  h_MR_R2_2l->Fill(MR, R2metel, w);
	  if(isTTallhad)
	    TTallhad->Fill("2l", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2l", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2l", w);
	  
	  if (nloosemuons == 0){
	    ofile.count("2el0mu",w);
	    h_MR_2el0mu->Fill(MR, w);
	    h_R2_2el0mu->Fill(R2metel, w);
	    h_MR_R2_2el0mu->Fill(MR, R2metel, w);
	    h_njets_2el0mu->Fill(sjet.size(),w);
	    h_nbjets_2el0mu->Fill(nmediumbs,w);
	    h_met_2el0mu->Fill(V3metel.Pt(),w);
	    h_jet1pt_2el0mu->Fill(sjet[0].pt,w);
	    h_jet2pt_2el0mu->Fill(sjet[1].pt,w);
	    h_jet3pt_2el0mu->Fill(sjet[2].pt,w);
	    h_lepton1pt_2el0mu->Fill(selectron[0].pt,w);
	    h_lepton2pt_2el0mu->Fill(selectron[1].pt,w);
	    h_HT_2el0mu->Fill(HT,w);
	    h_PV_2el0mu->Fill(eventhelperextra_numberOfPrimaryVertices,w);
	    if(isTTallhad)
	      TTallhad->Fill("2el0mu", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2el0mu", w);
	    else if(isTTdilep)
	      TTdilep->Fill("2el0mu", w);
	    ofile.count("2l0ol",w);
	    h_MR_2l0ol->Fill(MR, w);
	    h_R2_2l0ol->Fill(R2metel, w);
	    h_MR_R2_2l0ol->Fill(MR, R2metel, w);
	    h_njets_2l0ol->Fill(sjet.size(),w);
	    h_nbjets_2l0ol->Fill(nmediumbs,w);
	    h_met_2l0ol->Fill(V3metel.Pt(),w);
	    h_jet1pt_2l0ol->Fill(sjet[0].pt,w);
	    h_jet2pt_2l0ol->Fill(sjet[1].pt,w);
	    h_jet3pt_2l0ol->Fill(sjet[2].pt,w);
	    h_lepton1pt_2l0ol->Fill(selectron[0].pt,w);
	    h_lepton2pt_2l0ol->Fill(selectron[1].pt,w);
	    h_HT_2l0ol->Fill(HT,w);
	    if(isTTallhad)
	      TTallhad->Fill("2l0ol", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("2l0ol", w);
	    else if(isTTdilep)
	      TTdilep->Fill("2l0ol", w);
	    
	    if (nloosebs == 0){
	      if (eventhelper_isRealData!=1) {
		w = w*wCSVL;
	      }
	      ofile.count("0Lb2el0mu",w);
	      h_MR_0Lb2el0mu->Fill(MR, w);
	      h_R2_0Lb2el0mu->Fill(R2metel, w);
	      h_MR_R2_0Lb2el0mu->Fill(MR, R2metel, w);
	      h_njets_0Lb2el0mu->Fill(sjet.size(),w);
	      h_met_0Lb2el0mu->Fill(V3metel.Pt(),w);
	      h_jet1pt_0Lb2el0mu->Fill(sjet[0].pt,w);
	      h_jet2pt_0Lb2el0mu->Fill(sjet[1].pt,w);
	      h_jet3pt_0Lb2el0mu->Fill(sjet[2].pt,w);
	      h_lepton1pt_0Lb2el0mu->Fill(selectron[0].pt,w);
	      h_lepton2pt_0Lb2el0mu->Fill(selectron[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("0Lb2el0mu", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lb2el0mu", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lb2el0mu", w);
	      ofile.count("0Lb2l0ol",w);
	      h_MR_0Lb2l0ol->Fill(MR, w);
	      h_R2_0Lb2l0ol->Fill(R2metel, w);
	      h_MR_R2_0Lb2l0ol->Fill(MR, R2metel, w);
	      h_njets_0Lb2l0ol->Fill(sjet.size(),w);
	      h_met_0Lb2l0ol->Fill(V3metel.Pt(),w);
	      h_jet1pt_0Lb2l0ol->Fill(sjet[0].pt,w);
	      h_jet2pt_0Lb2l0ol->Fill(sjet[1].pt,w);
	      h_jet3pt_0Lb2l0ol->Fill(sjet[2].pt,w);
	      h_lepton1pt_0Lb2l0ol->Fill(selectron[0].pt,w);
	      h_lepton2pt_0Lb2l0ol->Fill(selectron[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("0Lb2l0ol", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("0Lb2l0ol", w);
	      else if(isTTdilep)
		TTdilep->Fill("0Lb2l0ol", w);

	      if (nYs >= 1){
		ofile.count("0Lbg1Y2el0mu",w);
		h_MR_0Lbg1Y2el0mu->Fill(MR, w);
		h_R2_0Lbg1Y2el0mu->Fill(R2metel, w);
		h_MR_R2_0Lbg1Y2el0mu->Fill(MR, R2metel, w);
		h_njets_0Lbg1Y2el0mu->Fill(sjet.size(),w);
		h_met_0Lbg1Y2el0mu->Fill(V3metel.Pt(),w);
		h_jet1pt_0Lbg1Y2el0mu->Fill(sjet[0].pt,w);
		h_jet2pt_0Lbg1Y2el0mu->Fill(sjet[1].pt,w);
		h_jet3pt_0Lbg1Y2el0mu->Fill(sjet[2].pt,w);
		h_lepton1pt_0Lbg1Y2el0mu->Fill(selectron[0].pt,w);
		h_lepton2pt_0Lbg1Y2el0mu->Fill(selectron[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y2el0mu", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y2el0mu", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y2el0mu", w);

		ofile.count("0Lbg1Y2l0ol",w);
		h_MR_0Lbg1Y2l0ol->Fill(MR, w);
		h_R2_0Lbg1Y2l0ol->Fill(R2metel, w);
		h_MR_R2_0Lbg1Y2l0ol->Fill(MR, R2metel, w);
		h_njets_0Lbg1Y2l0ol->Fill(sjet.size(),w);
		h_met_0Lbg1Y2l0ol->Fill(V3metel.Pt(),w);
		h_jet1pt_0Lbg1Y2l0ol->Fill(sjet[0].pt,w);
		h_jet2pt_0Lbg1Y2l0ol->Fill(sjet[1].pt,w);
		h_jet3pt_0Lbg1Y2l0ol->Fill(sjet[2].pt,w);
		h_lepton1pt_0Lbg1Y2l0ol->Fill(selectron[0].pt,w);
		h_lepton2pt_0Lbg1Y2l0ol->Fill(selectron[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("0Lbg1Y2l0ol", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("0Lbg1Y2l0ol", w);
		else if(isTTdilep)
		  TTdilep->Fill("0Lbg1Y2l0ol", w);
	      }// nYs >= 1
	    }// end nloosebs == 0
	    
	    if (nmediumbs >= 1){
	      if (eventhelper_isRealData!=1) {
		w = w*wCSVM;
	      }
	      ofile.count("g1Mb2el0mu",w);
	      h_MR_g1Mb2el0mu->Fill(MR, w);
	      h_R2_g1Mb2el0mu->Fill(R2metel, w);
	      h_MR_R2_g1Mb2el0mu->Fill(MR, R2metel, w);
	      h_njets_g1Mb2el0mu->Fill(sjet.size(),w);
	      h_nbjets_g1Mb2el0mu->Fill(nmediumbs,w);
	      h_met_g1Mb2el0mu->Fill(V3metel.Pt(),w);
	      h_jet1pt_g1Mb2el0mu->Fill(sjet[0].pt,w);
	      h_jet2pt_g1Mb2el0mu->Fill(sjet[1].pt,w);
	      h_jet3pt_g1Mb2el0mu->Fill(sjet[2].pt,w);
	      h_lepton1pt_g1Mb2el0mu->Fill(selectron[0].pt,w);
	      h_lepton2pt_g1Mb2el0mu->Fill(selectron[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mb2el0mu", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mb2el0mu", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mb2el0mu", w);

	      ofile.count("g1Mb2l0ol",w);
	      h_MR_g1Mb2l0ol->Fill(MR, w);
	      h_R2_g1Mb2l0ol->Fill(R2metel, w);
	      h_MR_R2_g1Mb2l0ol->Fill(MR, R2metel, w);
	      h_njets_g1Mb2l0ol->Fill(sjet.size(),w);
	      h_nbjets_g1Mb2l0ol->Fill(nmediumbs,w);
	      h_met_g1Mb2l0ol->Fill(V3metel.Pt(),w);
	      h_jet1pt_g1Mb2l0ol->Fill(sjet[0].pt,w);
	      h_jet2pt_g1Mb2l0ol->Fill(sjet[1].pt,w);
	      h_jet3pt_g1Mb2l0ol->Fill(sjet[2].pt,w);
	      h_lepton1pt_g1Mb2l0ol->Fill(selectron[0].pt,w);
	      h_lepton2pt_g1Mb2l0ol->Fill(selectron[1].pt,w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mb2l0ol", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mb2l0ol", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mb2l0ol", w);
	    
	      if (nYs >= 1){
		ofile.count("g1Mbg1Y2el0mu",w);
		h_MR_g1Mbg1Y2el0mu->Fill(MR, w);
		h_R2_g1Mbg1Y2el0mu->Fill(R2metel, w);
		h_MR_R2_g1Mbg1Y2el0mu->Fill(MR, R2metel, w);
		h_njets_g1Mbg1Y2el0mu->Fill(sjet.size(),w);
		h_nbjets_g1Mbg1Y2el0mu->Fill(nmediumbs,w);
		h_met_g1Mbg1Y2el0mu->Fill(V3metel.Pt(),w);
		h_jet1pt_g1Mbg1Y2el0mu->Fill(sjet[0].pt,w);
		h_jet2pt_g1Mbg1Y2el0mu->Fill(sjet[1].pt,w);
		h_jet3pt_g1Mbg1Y2el0mu->Fill(sjet[2].pt,w);
		h_lepton1pt_g1Mbg1Y2el0mu->Fill(selectron[0].pt,w);
		h_lepton2pt_g1Mbg1Y2el0mu->Fill(selectron[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1Y2el0mu", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1Y2el0mu", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1Y2el0mu", w);
		
		ofile.count("g1Mbg1Y2l0ol",w);
		h_MR_g1Mbg1Y2l0ol->Fill(MR, w);
		h_R2_g1Mbg1Y2l0ol->Fill(R2metel, w);
		h_MR_R2_g1Mbg1Y2l0ol->Fill(MR, R2metel, w);
		h_njets_g1Mbg1Y2l0ol->Fill(sjet.size(),w);
		h_nbjets_g1Mbg1Y2l0ol->Fill(nmediumbs,w);
		h_met_g1Mbg1Y2l0ol->Fill(V3metel.Pt(),w);
		h_jet1pt_g1Mbg1Y2l0ol->Fill(sjet[0].pt,w);
		h_jet2pt_g1Mbg1Y2l0ol->Fill(sjet[1].pt,w);
		h_jet3pt_g1Mbg1Y2l0ol->Fill(sjet[2].pt,w);
		h_lepton1pt_g1Mbg1Y2l0ol->Fill(selectron[0].pt,w);
		h_lepton2pt_g1Mbg1Y2l0ol->Fill(selectron[1].pt,w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1Y2l0ol", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1Y2l0ol", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1Y2l0ol", w);
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
