//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmd.h"
#include "utils.h"
#include "systutils.h"
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
  if (!fhlt){
    fhlt = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/hlteff/hlteff_HT_jpt_singlel.root");
  }
  if (!fhlt){
    cout << "Could not find trigger efficiency root file... Where did you put it??" << endl;
    return 1;
  }
  TH2D* h_hlteff = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph");
  TH2D* h_hlteff_up = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph_up");
  TH2D* h_hlteff_low = (TH2D*)fhlt->Get("h_HT_j1pt_0_effph_low");
    
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


  double sigmaJEC = -1.;
  double sigmaJECCA8 = -1.;
  double sigmaHLT = -1.;
  double sigmabtagFl = 0;
  double sigmabtagFs = 0;

  bool doPileupReweighting = false;
  if (sample != "Data"){
    doPileupReweighting = true;
  }
  // ---------------------------------------
  // --- Get the correct pileup histogram --
  // ---------------------------------------

  TString pileupname = "pileup_weights.root"; // default
  TFile* fpileup;
  fpileup = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/pileup/"+pileupname);
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
    //fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_TTJets.root");
  } else if (sample == "T1ttcc_DM10" || sample == "T1ttcc_DM25" || sample == "T1ttcc_DM80" 
	     || sample == "T1ttcc_old" || sample == "T2tt" || sample == "T1t1t") {
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_T1ttcc.root");
    //fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_T1ttcc.root");
  } else {
    fbeff = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/btageff/btageff_QCD.root");
    //fbeff = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/btageff/btageff_QCD.root");
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


  // Save the total weight for cross check purposes
  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);

  // main plots

  int nbn2 = 12;
  Double_t bn_tmp2[] = {200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,500.,1000.};
  Double_t* bn2 = getVariableBinEdges(nbn2+1,bn_tmp2);

  TH1D* h_pt_CA8jet = new TH1D("h_pt_CA8jet","h_pt_CA8jet",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass = new TH1D("h_pt_CA8jet_Wmass","h_pt_CA8jet_Wmass",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged = new TH1D("h_pt_CA8jet_Wmass_antitagged","h_pt_CA8jet_Wmass_antitagged",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged = new TH1D("h_pt_CA8jet_Wmass_tagged","h_pt_CA8jet_Wmass_tagged",nbn2,bn2);

  TH1D* h_pt_CA8jet_nomdPhi = new TH1D("h_pt_CA8jet_nomdPhi","h_pt_CA8jet_nomdPhi",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_nomdPhi","h_pt_CA8jet_Wmass_nomdPhi",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_antitagged_nomdPhi","h_pt_CA8jet_Wmass_antitagged_nomdPhi",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_tagged_nomdPhi","h_pt_CA8jet_Wmass_tagged_nomdPhi",nbn2,bn2);

  TH1D* h_pt_CA8jet_mdPhi0p5 = new TH1D("h_pt_CA8jet_mdPhi0p5","h_pt_CA8jet_mdPhi0p5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_mdPhi0p5","h_pt_CA8jet_Wmass_mdPhi0p5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_antitagged_mdPhi0p5","h_pt_CA8jet_Wmass_antitagged_mdPhi0p5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_tagged_mdPhi0p5","h_pt_CA8jet_Wmass_tagged_mdPhi0p5",nbn2,bn2);

  // Additional checks
  TH1D* h_pt_CA8jet_full = new TH1D("h_pt_CA8jet_full","h_pt_CA8jet_full",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full = new TH1D("h_pt_CA8jet_Wmass_full","h_pt_CA8jet_Wmass_full",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full = new TH1D("h_pt_CA8jet_Wmass_antitagged_full","h_pt_CA8jet_Wmass_antitagged_full",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full = new TH1D("h_pt_CA8jet_Wmass_tagged_full","h_pt_CA8jet_Wmass_tagged_full",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full_nomatch = new TH1D("h_pt_CA8jet_Wmass_tagged_full_nomatch","h_pt_CA8jet_Wmass_tagged_full_nomatch",50,0,1000);

  TH1D* h_pt_CA8jet_full_nomdPhi = new TH1D("h_pt_CA8jet_full_nomdPhi","h_pt_CA8jet_full_nomdPhi",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_full_nomdPhi","h_pt_CA8jet_Wmass_full_nomdPhi",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_antitagged_full_nomdPhi","h_pt_CA8jet_Wmass_antitagged_full_nomdPhi",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_tagged_full_nomdPhi","h_pt_CA8jet_Wmass_tagged_full_nomdPhi",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full_nomatch_nomdPhi = new TH1D("h_pt_CA8jet_Wmass_tagged_full_nomatch_nomdPhi","h_pt_CA8jet_Wmass_tagged_full_nomatch_nomdPhi",50,0,1000);

  TH1D* h_pt_CA8jet_full_mdPhi0p5 = new TH1D("h_pt_CA8jet_full_mdPhi0p5","h_pt_CA8jet_full_mdPhi0p5",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_full_mdPhi0p5","h_pt_CA8jet_Wmass_full_mdPhi0p5",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_antitagged_full_mdPhi0p5","h_pt_CA8jet_Wmass_antitagged_full_mdPhi0p5",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_tagged_full_mdPhi0p5","h_pt_CA8jet_Wmass_tagged_full_mdPhi0p5",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_tagged_full_nomatch_mdPhi0p5 = new TH1D("h_pt_CA8jet_Wmass_tagged_full_nomatch_mdPhi0p5","h_pt_CA8jet_Wmass_tagged_full_nomatch_mdPhi0p5",50,0,1000);

  TH1D* h_pt_CA8jet_full_eta1 = new TH1D("h_pt_CA8jet_full_eta1","h_pt_CA8jet_full_eta1",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full_eta1 = new TH1D("h_pt_CA8jet_Wmass_full_eta1","h_pt_CA8jet_Wmass_full_eta1",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full_eta1 = new TH1D("h_pt_CA8jet_Wmass_antitagged_full_eta1","h_pt_CA8jet_Wmass_antitagged_full_eta1",50,0,1000);

  TH1D* h_eta_CA8jet = new TH1D("h_eta_CA8jet","h_eta_CA8jet",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass = new TH1D("h_eta_CA8jet_Wmass","h_eta_CA8jet_Wmass",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass_antitagged = new TH1D("h_eta_CA8jet_Wmass_antitagged","h_eta_CA8jet_Wmass_antitagged",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass_tagged = new TH1D("h_eta_CA8jet_Wmass_tagged","h_eta_CA8jet_Wmass_tagged",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass_tagged_nomatch = new TH1D("h_eta_CA8jet_Wmass_tagged_nomatch","h_eta_CA8jet_Wmass_tagged_nomatch",50,-2.5,2.5);

  TH1D* h_mass_CA8jet = new TH1D("h_mass_CA8jet","h_mass_CA8jet",50,0,500);
  TH1D* h_mass_CA8jet_Wmass = new TH1D("h_mass_CA8jet_Wmass","h_mass_CA8jet_Wmass",50,0,500);
  TH1D* h_mass_CA8jet_Wmass_antitagged = new TH1D("h_mass_CA8jet_Wmass_antitagged","h_mass_CA8jet_Wmass_antitagged",50,0,500);
  TH1D* h_mass_CA8jet_Wmass_tagged = new TH1D("h_mass_CA8jet_Wmass_tagged","h_mass_CA8jet_Wmass_tagged",50,0,500);
  TH1D* h_mass_CA8jet_Wmass_tagged_nomatch = new TH1D("h_mass_CA8jet_Wmass_tagged_nomatch","h_mass_CA8jet_Wmass_tagged_nomatch",50,0,500);

  TH2D* h_mass_eta_CA8jet = new TH2D("h_mass_eta_CA8jet","h_mass_eta_CA8jet",50,0,500,50,-2.5,2.5);
  TH2D* h_mass_eta_CA8jet_Wmass = new TH2D("h_mass_eta_CA8jet_Wmass","h_mass_eta_CA8jet_Wmass",50,0,500,50,-2.5,2.5);
  TH2D* h_mass_eta_CA8jet_Wmass_antitagged = new TH2D("h_mass_eta_CA8jet_Wmass_antitagged","h_mass_eta_CA8jet_Wmass_antitagged",50,0,500,50,-2.5,2.5);

  TH2D* h_mass_pt_CA8jet = new TH2D("h_mass_pt_CA8jet","h_mass_pt_CA8jet",50,0,500,50,0,1000);
  TH2D* h_mass_pt_CA8jet_Wmass = new TH2D("h_mass_pt_CA8jet_Wmass","h_mass_pt_CA8jet_Wmass",50,0,500,50,0,1000);
  TH2D* h_mass_pt_CA8jet_Wmass_antitagged = new TH2D("h_mass_pt_CA8jet_Wmass_antitagged","h_mass_pt_CA8jet_Wmass_antitagged",50,0,500,50,0,1000);
  
  TH2D* h_pt_eta_CA8jet = new TH2D("h_pt_eta_CA8jet","h_pt_eta_CA8jet",50,0,1000,50,-2.5,2.5);
  TH2D* h_pt_eta_CA8jet_Wmass = new TH2D("h_pt_eta_CA8jet_Wmass","h_pt_eta_CA8jet_Wmass",50,0,1000,50,-2.5,2.5);
  TH2D* h_pt_eta_CA8jet_Wmass_antitagged = new TH2D("h_pt_eta_CA8jet_Wmass_antitagged","h_pt_eta_CA8jet_Wmass_antitagged",50,0,1000,50,-2.5,2.5);
  
  
  // relax all jet cuts
  TH1D* h_pt_CA8jet_full_nojet = new TH1D("h_pt_CA8jet_full_nojet","h_pt_CA8jet_full_nojet",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full_nojet = new TH1D("h_pt_CA8jet_Wmass_full_nojet","h_pt_CA8jet_Wmass_full_nojet",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full_nojet = new TH1D("h_pt_CA8jet_Wmass_antitagged_full_nojet","h_pt_CA8jet_Wmass_antitagged_full_nojet",50,0,1000);

  TH1D* h_pt_CA8jet_full_nojet_eta1 = new TH1D("h_pt_CA8jet_full_nojet_eta1","h_pt_CA8jet_full_nojet_eta1",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_full_nojet_eta1 = new TH1D("h_pt_CA8jet_Wmass_full_nojet_eta1","h_pt_CA8jet_Wmass_full_nojet_eta1",50,0,1000);
  TH1D* h_pt_CA8jet_Wmass_antitagged_full_nojet_eta1 = new TH1D("h_pt_CA8jet_Wmass_antitagged_full_nojet_eta1","h_pt_CA8jet_Wmass_antitagged_full_nojet_eta1",50,0,1000);

  TH1D* h_eta_CA8jet_nojet = new TH1D("h_eta_CA8jet_nojet","h_eta_CA8jet_nojet",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass_nojet = new TH1D("h_eta_CA8jet_Wmass_nojet","h_eta_CA8jet_Wmass_nojet",50,-2.5,2.5);
  TH1D* h_eta_CA8jet_Wmass_antitagged_nojet = new TH1D("h_eta_CA8jet_Wmass_antitagged_nojet","h_eta_CA8jet_Wmass_antitagged_nojet",50,-2.5,2.5);

  TH1D* h_mass_CA8jet_nojet = new TH1D("h_mass_CA8jet_nojet","h_mass_CA8jet_nojet",50,0,500);
  TH1D* h_mass_CA8jet_Wmass_nojet = new TH1D("h_mass_CA8jet_Wmass_nojet","h_mass_CA8jet_Wmass_nojet",50,0,500);
  TH1D* h_mass_CA8jet_Wmass_antitagged_nojet = new TH1D("h_mass_CA8jet_Wmass_antitagged_nojet","h_mass_CA8jet_Wmass_antitagged_nojet",50,0,500);

  TH2D* h_mass_eta_CA8jet_nojet = new TH2D("h_mass_eta_CA8jet_nojet","h_mass_eta_CA8jet_nojet",50,0,500,50,-2.5,2.5);
  TH2D* h_mass_eta_CA8jet_Wmass_nojet = new TH2D("h_mass_eta_CA8jet_Wmass_nojet","h_mass_eta_CA8jet_Wmass_nojet",50,0,500,50,-2.5,2.5);
  TH2D* h_mass_eta_CA8jet_Wmass_antitagged_nojet = new TH2D("h_mass_eta_CA8jet_Wmass_antitagged_nojet","h_mass_eta_CA8jet_Wmass_antitagged_nojet",50,0,500,50,-2.5,2.5);

  TH2D* h_mass_pt_CA8jet_nojet = new TH2D("h_mass_pt_CA8jet_nojet","h_mass_pt_CA8jet_nojet",50,0,500,50,0,1000);
  TH2D* h_mass_pt_CA8jet_Wmass_nojet = new TH2D("h_mass_pt_CA8jet_Wmass_nojet","h_mass_pt_CA8jet_Wmass_nojet",50,0,500,50,0,1000);
  TH2D* h_mass_pt_CA8jet_Wmass_antitagged_nojet = new TH2D("h_mass_pt_CA8jet_Wmass_antitagged_nojet","h_mass_pt_CA8jet_Wmass_antitagged_nojet",50,0,500,50,0,1000);
  
  TH2D* h_pt_eta_CA8jet_nojet = new TH2D("h_pt_eta_CA8jet_nojet","h_pt_eta_CA8jet_nojet",50,0,1000,50,-2.5,2.5);
  TH2D* h_pt_eta_CA8jet_Wmass_nojet = new TH2D("h_pt_eta_CA8jet_Wmass_nojet","h_pt_eta_CA8jet_Wmass_nojet",50,0,1000,50,-2.5,2.5);
  TH2D* h_pt_eta_CA8jet_Wmass_antitagged_nojet = new TH2D("h_pt_eta_CA8jet_Wmass_antitagged_nojet","h_pt_eta_CA8jet_Wmass_antitagged_nojet",50,0,1000,50,-2.5,2.5);



  // in bins of njets

  TH1D* h_pt_CA8jet_njet3 = new TH1D("h_pt_CA8jet_njet3","h_pt_CA8jet_njet3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_njet3 = new TH1D("h_pt_CA8jet_Wmass_njet3","h_pt_CA8jet_Wmass_njet3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_njet3 = new TH1D("h_pt_CA8jet_Wmass_antitagged_njet3","h_pt_CA8jet_Wmass_antitagged_njet3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_njet3 = new TH1D("h_pt_CA8jet_Wmass_tagged_njet3","h_pt_CA8jet_Wmass_tagged_njet3",nbn2,bn2);

  TH1D* h_pt_CA8jet_njet4 = new TH1D("h_pt_CA8jet_njet4","h_pt_CA8jet_njet4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_njet4 = new TH1D("h_pt_CA8jet_Wmass_njet4","h_pt_CA8jet_Wmass_njet4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_njet4 = new TH1D("h_pt_CA8jet_Wmass_antitagged_njet4","h_pt_CA8jet_Wmass_antitagged_njet4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_njet4 = new TH1D("h_pt_CA8jet_Wmass_tagged_njet4","h_pt_CA8jet_Wmass_tagged_njet4",nbn2,bn2);

  TH1D* h_pt_CA8jet_njet5 = new TH1D("h_pt_CA8jet_njet5","h_pt_CA8jet_njet5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_njet5 = new TH1D("h_pt_CA8jet_Wmass_njet5","h_pt_CA8jet_Wmass_njet5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_njet5 = new TH1D("h_pt_CA8jet_Wmass_antitagged_njet5","h_pt_CA8jet_Wmass_antitagged_njet5",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_njet5 = new TH1D("h_pt_CA8jet_Wmass_tagged_njet5","h_pt_CA8jet_Wmass_tagged_njet5",nbn2,bn2);

  // in bins of HT: 1 = [-500]; 2 = [500-600]; 3 = [600-700]; 4 = [700-]

  TH1D* h_pt_CA8jet_HT1 = new TH1D("h_pt_CA8jet_HT1","h_pt_CA8jet_HT1",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_HT1 = new TH1D("h_pt_CA8jet_Wmass_HT1","h_pt_CA8jet_Wmass_HT1",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_HT1 = new TH1D("h_pt_CA8jet_Wmass_antitagged_HT1","h_pt_CA8jet_Wmass_antitagged_HT1",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_HT1 = new TH1D("h_pt_CA8jet_Wmass_tagged_HT1","h_pt_CA8jet_Wmass_tagged_HT1",nbn2,bn2);

  TH1D* h_pt_CA8jet_HT2 = new TH1D("h_pt_CA8jet_HT2","h_pt_CA8jet_HT2",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_HT2 = new TH1D("h_pt_CA8jet_Wmass_HT2","h_pt_CA8jet_Wmass_HT2",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_HT2 = new TH1D("h_pt_CA8jet_Wmass_antitagged_HT2","h_pt_CA8jet_Wmass_antitagged_HT2",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_HT2 = new TH1D("h_pt_CA8jet_Wmass_tagged_HT2","h_pt_CA8jet_Wmass_tagged_HT2",nbn2,bn2);

  TH1D* h_pt_CA8jet_HT3 = new TH1D("h_pt_CA8jet_HT3","h_pt_CA8jet_HT3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_HT3 = new TH1D("h_pt_CA8jet_Wmass_HT3","h_pt_CA8jet_Wmass_HT3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_HT3 = new TH1D("h_pt_CA8jet_Wmass_antitagged_HT3","h_pt_CA8jet_Wmass_antitagged_HT3",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_HT3 = new TH1D("h_pt_CA8jet_Wmass_tagged_HT3","h_pt_CA8jet_Wmass_tagged_HT3",nbn2,bn2);

  TH1D* h_pt_CA8jet_HT4 = new TH1D("h_pt_CA8jet_HT4","h_pt_CA8jet_HT4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_HT4 = new TH1D("h_pt_CA8jet_Wmass_HT4","h_pt_CA8jet_Wmass_HT4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_antitagged_HT4 = new TH1D("h_pt_CA8jet_Wmass_antitagged_HT4","h_pt_CA8jet_Wmass_antitagged_HT4",nbn2,bn2);
  TH1D* h_pt_CA8jet_Wmass_tagged_HT4 = new TH1D("h_pt_CA8jet_Wmass_tagged_HT4","h_pt_CA8jet_Wmass_tagged_HT4",nbn2,bn2);


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

      // ---------------------------------
      // -- Trigger requirement in Data --
      // ---------------------------------

      if (eventhelper_isRealData==1) { 
	// Run2012A:
	if ( eventhelper_run >= 190456 && eventhelper_run < 190762 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (! (triggerresultshelper_HLT_PFJet320_v3 == 1)) 
	      continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (! (triggerresultshelper_HLT_PFJet320_v3 == 0 
		   &&
		   triggerresultshelper_HLT_PFHT650_v5 == 1 
		   )) continue;
	  } else {
	    if (! (triggerresultshelper_HLT_PFJet320_v3 == 1
		   || triggerresultshelper_HLT_PFHT650_v5 == 1 
		   ) ) continue;
	  }
	}
	
	if ( eventhelper_run >= 190762 && eventhelper_run < 191512 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 0 
		 &&
		 triggerresultshelper_HLT_PFHT650_v6 == 1
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v4 == 1 
		 || triggerresultshelper_HLT_PFHT650_v6 == 1 
		 ) ) continue;
	  }
	}

	if ( eventhelper_run >= 191512 && eventhelper_run < 193746 ) {
	  if (fsample.find("_Jet_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 
		 )) continue;
	  } else if (fsample.find("_HT_") < fsample.length()) {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 0 
		 &&
		 triggerresultshelper_HLT_PFHT650_v7 == 1
		 )) continue;
	  } else {
	    if (!
		(triggerresultshelper_HLT_PFJet320_v5 == 1 
		 || triggerresultshelper_HLT_PFHT650_v7 == 1 
		 ) ) continue;
	  }
	}

	// 2012B
	if ( eventhelper_run >= 193746 && eventhelper_run < 196039) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 
	       || triggerresultshelper_HLT_PFHT650_v8 == 1 
	       ) ) continue;
	}

	if ( eventhelper_run >= 196039 && eventhelper_run < 197770 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v5 == 1 
	       || triggerresultshelper_HLT_PFHT650_v9 == 1 
	       ) ) continue;
	}

	// Run2012C:
	if ( eventhelper_run >= 197770 && eventhelper_run < 199648 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet400_v6 == 1 
	       || triggerresultshelper_HLT_PFNoPUHT650_v1 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 199648 && eventhelper_run < 202807 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v8 == 1 
	       || triggerresultshelper_HLT_PFNoPUHT650_v3 == 1
	       ) ) continue;
	}

	if ( eventhelper_run >= 202807 && eventhelper_run < 203734 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1
	       || triggerresultshelper_HLT_PFNoPUHT650_v4 == 1
	       ) ) continue;
	}

	// Run2012D:
	if ( eventhelper_run >= 203734 && eventhelper_run < 208940 ) {
	  if (!
              (triggerresultshelper_HLT_PFJet320_v9 == 1 
	       || triggerresultshelper_HLT_PFNoPUHT650_v4 == 1 
	       ) ) continue;
	}
      }

      // ------------------------
      // -- Pileup reweighting --
      // ------------------------

      double num_vertices = pileupsummaryinfo[0].getTrueNumInteractions;
      // get bin number in pileup histogram
      int pileup_bin = (int)ceil(num_vertices);
      double w_pileup = 1.;
      if(doPileupReweighting)
	w_pileup = h_pileup->GetBinContent(pileup_bin);      

      w = w*w_pileup;

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
	// Put the jet in a TLorentzVector and scale it with JEC SF
	TLorentzVector jlnojecSF;
	jlnojecSF.SetPtEtaPhiE(cmgpfjet[i].pt, cmgpfjet[i].eta,
			cmgpfjet[i].phi, cmgpfjet[i].energy);
	TLorentzVector jl;
	jl.SetPxPyPzE(jlnojecSF.Px()*jecSF, jlnojecSF.Py()*jecSF, jlnojecSF.Pz()*jecSF, jlnojecSF.E()*jecSF);
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
	  SFCSVL = (SFCSVLFl + sigmabtagFl*dSFCSVLFl)*(SFCSVLFs + sigmabtagFs*dSFCSVLFs);
	  SFCSVM = (SFCSVMFl + sigmabtagFl*dSFCSVMFl)*(SFCSVMFs + sigmabtagFs*dSFCSVMFs);
	} else { // FullSim
	  if (fabs(partonFlavour) == 5) {
	    eCSVM = geteff1D(h_pt_b_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_b_CSVLeff, pt);
	  }
	  if (fabs(partonFlavour) != 5) {
	    eCSVM = geteff1D(h_pt_lc_CSVMeff, pt);
	    eCSVL = geteff1D(h_pt_lc_CSVLeff, pt);
	  }
	  SFCSVL = (SFCSVLFl + sigmabtagFl*dSFCSVLFl);
	  SFCSVM = (SFCSVMFl + sigmabtagFl*dSFCSVMFl);
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
	// begin JEC SF
	double jecUnc = 0.;
	AK5PFCHSJECunc(jethelper4[i].pt, jethelper4[i].eta, jecUnc);
	double jecSFCA8 = 1. + (sigmaJECCA8 * fabs(jecUnc));
	// Put the jet in a TLorentzVector and scale it with JEC SF
	TLorentzVector jlCA8nojecSF;
	jlCA8nojecSF.SetPtEtaPhiE(jethelper4[i].pt, jethelper4[i].eta,
				  jethelper4[i].phi, jethelper4[i].energy);
	TLorentzVector jlCA8;
	jlCA8.SetPxPyPzE(jlCA8nojecSF.Px()*jecSFCA8, jlCA8nojecSF.Py()*jecSFCA8, 
			 jlCA8nojecSF.Pz()*jecSFCA8, jlCA8nojecSF.E()*jecSFCA8);
	jethelper4[i].pt = jlCA8.Pt();
	jethelper4[i].mass = jlCA8.M();
	if (!(jethelper4[i].pt > 200) ) continue;
        if (!(fabs(jethelper4[i].eta) < 2.4) ) continue;
        sjet2.push_back(jethelper4[i]);
	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
	sY.push_back(jethelper4[i]);
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
	if (tau21 >= 0.50) {
          aW.push_back(jethelper4[i]);
        }
	if (!(tau21 < 0.50) ) continue;
        sW.push_back(jethelper4[i]);
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


      // Propagate JEC to MET
      TVector3 V3metnojecSF;
      V3metnojecSF.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      TLorentzVector metnojecSF;
      metnojecSF.SetPtEtaPhiE(cmgbasemet2[0].pt, 0, cmgbasemet2[0].phi, cmgbasemet2[0].energy);
      // MET with JEC SF
      TVector3 V3met;
      V3met.SetXYZ(V3metnojecSF.Px()-METcorrfromJEC_px, V3metnojecSF.Py()-METcorrfromJEC_py, 0.0);
      TLorentzVector met;
      met.SetPxPyPzE(V3met.Px(), V3met.Py(), 0.0, V3met.Pt());


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
		//cout << xmin << " " << MR << " " << xmax << " " << ymin << " " << R2 << " " << ymax << " " << whlt << endl;
		break;
	      }
	    }
	  }
	}
      }

      
      // ---------------------
      // -- event selection --
      // ---------------------

      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;

      // at least one good primary vertex
      if (!(svertex.size() > 0)) continue;

      // at least three jets
      if (!(sjet.size() >= 1)) continue;
      
      // Apply HLT weight and include it in the total weight:
      w = w*w_trigger;
    
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
      
	
      // ----------------------------------------------------------------------------------------------------
      // 0 Lepton trajectory
      // ----------------------------------------------------------------------------------------------------
      if (nlooseelectrons == 0){
	
	if (nloosemuons == 0) {
	  
	  if (eventhelperextra_trackIso == 0){
	    
	      if (nloosebs == 0){
		if (eventhelper_isRealData!=1) {
		  w = w*wCSVL;
		}

		if (sjet2.size() >0){
		  h_pt_CA8jet_full_nojet->Fill(sjet2[0].pt,w);
		  if(fabs(sjet2[0].eta)<1)
		    h_pt_CA8jet_full_nojet_eta1->Fill(sjet2[0].pt,w);
		  h_eta_CA8jet_nojet->Fill(sjet2[0].eta,w);
		  h_mass_CA8jet_nojet->Fill(sjet2[0].mass,w);
		  h_mass_eta_CA8jet_nojet->Fill(sjet2[0].mass,sjet2[0].eta,w);
		  h_pt_eta_CA8jet_nojet->Fill(sjet2[0].pt,sjet2[0].eta,w);
 		  h_mass_pt_CA8jet_nojet->Fill(sjet2[0].mass,sjet2[0].pt,w);

		  if (sjet2[0].mass > 70 && sjet2[0].mass < 100){
		    h_pt_CA8jet_Wmass_full_nojet->Fill(sjet2[0].pt,w);
		    if(fabs(sjet2[0].eta)<1)
		      h_pt_CA8jet_Wmass_full_nojet_eta1->Fill(sjet2[0].pt,w);
		    h_eta_CA8jet_Wmass_nojet->Fill(sjet2[0].eta,w);
		    h_mass_CA8jet_Wmass_nojet->Fill(sjet2[0].mass,w);
		    h_mass_eta_CA8jet_Wmass_nojet->Fill(sjet2[0].mass,sjet2[0].eta,w);
		    h_pt_eta_CA8jet_Wmass_nojet->Fill(sjet2[0].pt,sjet2[0].eta,w);
		    h_mass_pt_CA8jet_Wmass_nojet->Fill(sjet2[0].mass,sjet2[0].pt,w);

		    // Match with the unpruned:
		    double prjmatch = 0;
		    int jpr = -1;
		    double dRmn = 100;
		    for (unsigned int j=0; j<jethelper5.size(); j++) {
		      double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, sjet2[0].eta, sjet2[0].phi);
		      if (dR < 0.7 && dR < dRmn) {
			dRmn = dR;
			prjmatch = 1;
			jpr = j;
			break;
		      }
		    }
		    if (prjmatch==1){
		      double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
		      if (tau21 >= 0.50) {
			h_pt_CA8jet_Wmass_antitagged_full_nojet->Fill(sjet2[0].pt,w);
			if(fabs(sjet2[0].eta) < 1)
			  h_pt_CA8jet_Wmass_antitagged_full_nojet_eta1->Fill(sjet2[0].pt,w);
			h_eta_CA8jet_Wmass_antitagged_nojet->Fill(sjet2[0].eta,w);
			h_mass_CA8jet_Wmass_antitagged_nojet->Fill(sjet2[0].mass,w);
			h_mass_eta_CA8jet_Wmass_antitagged_nojet->Fill(sjet2[0].mass,sjet2[0].eta,w);
			h_pt_eta_CA8jet_Wmass_antitagged_nojet->Fill(sjet2[0].pt,sjet2[0].eta,w);
			h_mass_pt_CA8jet_Wmass_antitagged_nojet->Fill(sjet2[0].mass,sjet2[0].pt,w);
		      }
		      
		    }
		  }
		}

		// ---------------------------------------------------------------------------------------
		// -- Here start the regions we want to compare
		// ---------------------------------------------------------------------------------------
		if (!(sjet.size() >= 3)) continue;
		
		// Compute the minDeltaPhi variable, taking the first three jets into account
		double minDeltaPhi = 99.;
		for (int jet=0; jet<3; ++jet){
		  double mdphi = fdeltaPhi(sjet[jet].phi,V3met.Phi());
		  if (mdphi < minDeltaPhi)
		    minDeltaPhi = mdphi;
		}
		
		// pt of first jet greater than 200 GeV
		if (!(sjet[0].pt > 200)) continue;
		if (sjet2.size() == 0) continue;

		// No mindeltaphi cut
		if(sjet2[0].pt>=1000){
		  h_pt_CA8jet_nomdPhi->Fill(999,w);
		} else { 
		  h_pt_CA8jet_nomdPhi->Fill(sjet2[0].pt,w);
		}
		h_pt_CA8jet_full_nomdPhi->Fill(sjet2[0].pt,w);
		    
		if (sjet2[0].mass > 70 && sjet2[0].mass < 100){
		  if(sjet2[0].pt>=1000){
		    h_pt_CA8jet_Wmass_nomdPhi->Fill(999,w);
		  } else {
		    h_pt_CA8jet_Wmass_nomdPhi->Fill(sjet2[0].pt,w);
		  }
		  h_pt_CA8jet_Wmass_full_nomdPhi->Fill(sjet2[0].pt,w);

		  // Match with the unpruned:
		  double prjmatch = 0;
		  int jpr = -1;
		  double dRmn = 100;
		  for (unsigned int j=0; j<jethelper5.size(); j++) {
		    double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, sjet2[0].eta, sjet2[0].phi);
		    if (dR < 0.7 && dR < dRmn) {
		      dRmn = dR;
		      prjmatch = 1;
		      jpr = j;
		      break;
		    }
		  }
		  if (prjmatch==1){
		    double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
		    if (tau21 >= 0.50) {
		      if(sjet2[0].pt>=1000){
			h_pt_CA8jet_Wmass_antitagged_nomdPhi->Fill(999,w);
		      } else {
			h_pt_CA8jet_Wmass_antitagged_nomdPhi->Fill(sjet2[0].pt,w);
		      }
		      h_pt_CA8jet_Wmass_antitagged_full_nomdPhi->Fill(sjet2[0].pt,w);
		    } else {
		      if(sjet2[0].pt>=1000){
			h_pt_CA8jet_Wmass_tagged_nomdPhi->Fill(999,w);
		      } else {
			h_pt_CA8jet_Wmass_tagged_nomdPhi->Fill(sjet2[0].pt,w);
		      }
		      h_pt_CA8jet_Wmass_tagged_full_nomdPhi->Fill(sjet2[0].pt,w);
		      // Check whether the jet matches with a W
		      bool matched = false;
		      for (unsigned int i=0; i<genparticlehelper.size() && !matched; i++) {
			if (genparticlehelper[i].status != 3) continue;
			if (fabs(genparticlehelper[i].pdgId) != 24) continue;
			double dr = fdeltaR(sjet2[0].eta,sjet2[0].phi, genparticlehelper[i].eta,genparticlehelper[i].phi);
			if (dr < 0.8){
			  matched = true;
			}
		      }
		      if(!matched){
			h_pt_CA8jet_Wmass_tagged_full_nomatch_nomdPhi->Fill(sjet2[0].pt,w);
		      }
		    }
		  }
		}


		// mindeltaphi > 0.5
		if(minDeltaPhi > 0.5){
		  if(sjet2[0].pt>=1000){
		    h_pt_CA8jet_mdPhi0p5->Fill(999,w);
		  } else { 
		    h_pt_CA8jet_mdPhi0p5->Fill(sjet2[0].pt,w);
		  }
		  h_pt_CA8jet_full_mdPhi0p5->Fill(sjet2[0].pt,w);
		  
		  if (sjet2[0].mass > 70 && sjet2[0].mass < 100){
		    if(sjet2[0].pt>=1000){
		      h_pt_CA8jet_Wmass_mdPhi0p5->Fill(999,w);
		    } else {
		      h_pt_CA8jet_Wmass_mdPhi0p5->Fill(sjet2[0].pt,w);
		    }
		    h_pt_CA8jet_Wmass_full_mdPhi0p5->Fill(sjet2[0].pt,w);
		    
		    // Match with the unpruned:
		    double prjmatch = 0;
		    int jpr = -1;
		    double dRmn = 100;
		    for (unsigned int j=0; j<jethelper5.size(); j++) {
		      double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, sjet2[0].eta, sjet2[0].phi);
		      if (dR < 0.7 && dR < dRmn) {
			dRmn = dR;
			prjmatch = 1;
			jpr = j;
			break;
		      }
		    }
		    if (prjmatch==1){
		      double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
		      if (tau21 >= 0.50) {
			if(sjet2[0].pt>=1000){
			  h_pt_CA8jet_Wmass_antitagged_mdPhi0p5->Fill(999,w);
			} else {
			  h_pt_CA8jet_Wmass_antitagged_mdPhi0p5->Fill(sjet2[0].pt,w);
			}
			h_pt_CA8jet_Wmass_antitagged_full_mdPhi0p5->Fill(sjet2[0].pt,w);
		      } else {
			if(sjet2[0].pt>=1000){
			  h_pt_CA8jet_Wmass_tagged_mdPhi0p5->Fill(999,w);
			} else {
			  h_pt_CA8jet_Wmass_tagged_mdPhi0p5->Fill(sjet2[0].pt,w);
			}
			h_pt_CA8jet_Wmass_tagged_full_mdPhi0p5->Fill(sjet2[0].pt,w);
			// Check whether the jet matches with a W
			bool matched = false;
			for (unsigned int i=0; i<genparticlehelper.size() && !matched; i++) {
			  if (genparticlehelper[i].status != 3) continue;
			  if (fabs(genparticlehelper[i].pdgId) != 24) continue;
			  double dr = fdeltaR(sjet2[0].eta,sjet2[0].phi, genparticlehelper[i].eta,genparticlehelper[i].phi);
			  if (dr < 0.8){
			    matched = true;
			  }
			}
			if(!matched){
			  h_pt_CA8jet_Wmass_tagged_full_nomatch_mdPhi0p5->Fill(sjet2[0].pt,w);
			}
		      }
		    }
		  }
		}


		if (minDeltaPhi < 0.3){
		  if(sjet2.size() > 0){
		    if(sjet2[0].pt>=1000){
		      h_pt_CA8jet->Fill(999,w);
		    } else { 
		      h_pt_CA8jet->Fill(sjet2[0].pt,w);
		    }
		    if(sjet.size() == 3){
		      h_pt_CA8jet_njet3->Fill(sjet2[0].pt,w);
		    } else if(sjet.size() == 4){
		      h_pt_CA8jet_njet4->Fill(sjet2[0].pt,w);
		    } else {
		      h_pt_CA8jet_njet5->Fill(sjet2[0].pt,w);
		    }
		    if(HT < 500){
		      h_pt_CA8jet_HT1->Fill(sjet2[0].pt,w);
		    } else if(HT < 600){
		      h_pt_CA8jet_HT2->Fill(sjet2[0].pt,w);
		    } else if(HT < 700){
		      h_pt_CA8jet_HT3->Fill(sjet2[0].pt,w);
		    } else {
		      h_pt_CA8jet_HT4->Fill(sjet2[0].pt,w);
		    }

		    h_pt_CA8jet_full->Fill(sjet2[0].pt,w);
		    if(fabs(sjet2[0].eta)<1)
		      h_pt_CA8jet_full_eta1->Fill(sjet2[0].pt,w);
		    h_eta_CA8jet->Fill(sjet2[0].eta,w);
		    h_mass_CA8jet->Fill(sjet2[0].mass,w);
		    h_mass_eta_CA8jet->Fill(sjet2[0].mass,sjet2[0].eta,w);
		    h_pt_eta_CA8jet->Fill(sjet2[0].pt,sjet2[0].eta,w);
		    h_mass_pt_CA8jet->Fill(sjet2[0].mass,sjet2[0].pt,w);
		    
		    if (sjet2[0].mass > 70 && sjet2[0].mass < 100){
		      if(sjet2[0].pt>=1000){
			h_pt_CA8jet_Wmass->Fill(999,w);
		      } else {
			h_pt_CA8jet_Wmass->Fill(sjet2[0].pt,w);
		      }
		      if(sjet.size() == 3){
			h_pt_CA8jet_Wmass_njet3->Fill(sjet2[0].pt,w);
		      } else if(sjet.size() == 4){
			h_pt_CA8jet_Wmass_njet4->Fill(sjet2[0].pt,w);
		      } else {
			h_pt_CA8jet_Wmass_njet5->Fill(sjet2[0].pt,w);
		      }
		      if(HT < 500){
			h_pt_CA8jet_Wmass_HT1->Fill(sjet2[0].pt,w);
		      } else if(HT < 600){
			h_pt_CA8jet_Wmass_HT2->Fill(sjet2[0].pt,w);
		      } else if(HT < 700){
			h_pt_CA8jet_Wmass_HT3->Fill(sjet2[0].pt,w);
		      } else {
			h_pt_CA8jet_Wmass_HT4->Fill(sjet2[0].pt,w);
		      }

		      h_pt_CA8jet_Wmass_full->Fill(sjet2[0].pt,w);
		      if(fabs(sjet2[0].eta)<1)
			h_pt_CA8jet_Wmass_full_eta1->Fill(sjet2[0].pt,w);
		      h_eta_CA8jet_Wmass->Fill(sjet2[0].eta,w);
		      h_mass_CA8jet_Wmass->Fill(sjet2[0].mass,w);
		      h_mass_eta_CA8jet_Wmass->Fill(sjet2[0].mass,sjet2[0].eta,w);
		      h_pt_eta_CA8jet_Wmass->Fill(sjet2[0].pt,sjet2[0].eta,w);
		      h_mass_pt_CA8jet_Wmass->Fill(sjet2[0].mass,sjet2[0].pt,w);
		      
		      // Match with the unpruned:
		      double prjmatch = 0;
		      int jpr = -1;
		      double dRmn = 100;
		      for (unsigned int j=0; j<jethelper5.size(); j++) {
			double dR = fdeltaR(jethelper5[j].eta, jethelper5[j].phi, sjet2[0].eta, sjet2[0].phi);
			if (dR < 0.7 && dR < dRmn) {
			  dRmn = dR;
			  prjmatch = 1;
			  jpr = j;
			  break;
			}
		      }
		      if (prjmatch==1){
			double tau21 = jethelper5[jpr].tau2 / jethelper5[jpr].tau1;
			if (tau21 >= 0.50) {
			  if(sjet2[0].pt>=1000){
			    h_pt_CA8jet_Wmass_antitagged->Fill(999,w);
			  } else {
			    h_pt_CA8jet_Wmass_antitagged->Fill(sjet2[0].pt,w);
			  }
			  if(sjet.size() == 3){
			    h_pt_CA8jet_Wmass_antitagged_njet3->Fill(sjet2[0].pt,w);
			  } else if(sjet.size() == 4){
			    h_pt_CA8jet_Wmass_antitagged_njet4->Fill(sjet2[0].pt,w);
			  } else {
			    h_pt_CA8jet_Wmass_antitagged_njet5->Fill(sjet2[0].pt,w);
			  }
			  if(HT < 500){
			    h_pt_CA8jet_Wmass_antitagged_HT1->Fill(sjet2[0].pt,w);
			  } else if(HT < 600){
			    h_pt_CA8jet_Wmass_antitagged_HT2->Fill(sjet2[0].pt,w);
			  } else if(HT < 700){
			    h_pt_CA8jet_Wmass_antitagged_HT3->Fill(sjet2[0].pt,w);
			  } else {
			    h_pt_CA8jet_Wmass_antitagged_HT4->Fill(sjet2[0].pt,w);
			  }
			  h_pt_CA8jet_Wmass_antitagged_full->Fill(sjet2[0].pt,w);
			  if(fabs(sjet2[0].eta)<1)
			    h_pt_CA8jet_Wmass_antitagged_full_eta1->Fill(sjet2[0].pt,w);
			  h_eta_CA8jet_Wmass_antitagged->Fill(sjet2[0].eta,w);
			  h_mass_CA8jet_Wmass_antitagged->Fill(sjet2[0].mass,w);
			  h_mass_eta_CA8jet_Wmass_antitagged->Fill(sjet2[0].mass,sjet2[0].eta,w);
			  h_pt_eta_CA8jet_Wmass_antitagged->Fill(sjet2[0].pt,sjet2[0].eta,w);
			  h_mass_pt_CA8jet_Wmass_antitagged->Fill(sjet2[0].mass,sjet2[0].pt,w);
			} else {
			  if(sjet2[0].pt>=1000){
			    h_pt_CA8jet_Wmass_tagged->Fill(999,w);
			  } else {
			    h_pt_CA8jet_Wmass_tagged->Fill(sjet2[0].pt,w);
			  }
			  if(sjet.size() == 3){
			    h_pt_CA8jet_Wmass_tagged_njet3->Fill(sjet2[0].pt,w);
			  } else if(sjet.size() == 4){
			    h_pt_CA8jet_Wmass_tagged_njet4->Fill(sjet2[0].pt,w);
			  } else {
			    h_pt_CA8jet_Wmass_tagged_njet5->Fill(sjet2[0].pt,w);
			  }
			  if(HT < 500){
			    h_pt_CA8jet_Wmass_tagged_HT1->Fill(sjet2[0].pt,w);
			  } else if(HT < 600){
			    h_pt_CA8jet_Wmass_tagged_HT2->Fill(sjet2[0].pt,w);
			  } else if(HT < 700){
			    h_pt_CA8jet_Wmass_tagged_HT3->Fill(sjet2[0].pt,w);
			  } else {
			    h_pt_CA8jet_Wmass_tagged_HT4->Fill(sjet2[0].pt,w);
			  }
			  h_pt_CA8jet_Wmass_tagged_full->Fill(sjet2[0].pt,w);
			  h_eta_CA8jet_Wmass_tagged->Fill(sjet2[0].eta,w);
			  h_mass_CA8jet_Wmass_tagged->Fill(sjet2[0].mass,w);
			  // Check whether the jet matches with a W
			  bool matched = false;
			  for (unsigned int i=0; i<genparticlehelper.size() && !matched; i++) {
			    if (genparticlehelper[i].status != 3) continue;
			    if (fabs(genparticlehelper[i].pdgId) != 24) continue;
			    double dr = fdeltaR(sjet2[0].eta,sjet2[0].phi, genparticlehelper[i].eta,genparticlehelper[i].phi);
			    if (dr < 0.8){
			      matched = true;
			    }
			  }
			  if(!matched){
			  h_pt_CA8jet_Wmass_tagged_full_nomatch->Fill(sjet2[0].pt,w);
			  h_eta_CA8jet_Wmass_tagged_nomatch->Fill(sjet2[0].eta,w);
			  h_mass_CA8jet_Wmass_tagged_nomatch->Fill(sjet2[0].mass,w);
			  }
			}
		      }
		    }
		  }
		} // end of mdphi		
	      } // end of nloosebs == 0
	  } // end veto iso track
	} // end veto muon
      }  // end veto electron
      
    } // end event loop
  
  fhlt->Close();
  fpileup->Close();
  fbeff->Close();
  stream.close();
  ofile.close();
  return 0;
}
