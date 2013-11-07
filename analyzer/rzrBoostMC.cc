//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmd.h"
#include "utils.h"

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
  TFile* fhlt = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/hlteff/extr_eff0_sm2.root");
  TH2D* h_hlteff = (TH2D*)fhlt->Get("hBinValues");
  
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
  
  outputFile ofile(cmdline.outputfilename);

  string fsample = cmdline.filelist;
  double xsect = cmdline.xsect;
  double totweight = cmdline.totweight;
  double lumi = cmdline.lumi;

  // Calculate the normalization factor for the event weights
  // The original MC weight will be divided by this quantity
  double weightnorm = 1.;
  if (xsect != -1 && totweight != -1 && lumi != -1) {
    weightnorm = (xsect*lumi)/totweight;
  }

  cout << "lumi: " << lumi << endl;
  cout << "xsect: " << xsect << endl;
  cout << "totweight: " << totweight << endl;
  cout << "weightnorm: " << weightnorm << endl;

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------

  double MRmx = 4000;

  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);

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

  // Plots of relevant quantities with tag for the cut flow

  TH1D* h_nW_Cleaning = new TH1D("h_nW_Cleaning", "h_nW_Cleaning", 5, 0, 5);
  TH1D* h_nb_Cleaning = new TH1D("h_nb_Cleaning", "h_nb_Cleaning", 5, 0, 5);
  TH2D* h_nW_nb_Cleaning = new TH2D("h_nW_nb_Cleaning", "h_nW_nb_Cleaning", 5, 0, 5, 5, 0, 5);
  TH1D* h_nWAK5_Cleaning = new TH1D("h_nWAK5_Cleaning", "h_nWAK5_Cleaning", 5, 0, 5);
  TH2D* h_met_R2_Cleaning = new TH2D("h_met_R2_Cleaning", "h_met_R2_Cleaning", 50, 0, 1000, 25, 0, 1);
  TH1D* h_met_Cleaning = new TH1D("h_met_Cleaning", "h_met_Cleaning", 50, 0, 1000);
  TH1D* h_j1pt_Cleaning = new TH1D("h_j1pt_Cleaning", "h_j1pt_Cleaning", 50, 0, 1000);
  TH1D* h_nj_Cleaning = new TH1D("h_nj_Cleaning", "h_nj_Cleaning", 20, 0, 20);

  TH1D* h_nW_jet1ptg200 = new TH1D("h_nW_jet1ptg200", "h_nW_jet1ptg200", 5, 0, 5);
  TH1D* h_nb_jet1ptg200 = new TH1D("h_nb_jet1ptg200", "h_nb_jet1ptg200", 5, 0, 5);
  TH2D* h_nW_nb_jet1ptg200 = new TH2D("h_nW_nb_jet1ptg200", "h_nW_nb_jet1ptg200", 5, 0, 5, 5, 0, 5);
  TH1D* h_nWAK5_jet1ptg200 = new TH1D("h_nWAK5_jet1ptg200", "h_nWAK5_jet1ptg200", 5, 0, 5);

  // MR, R2 plots for the different steps in the selection
  // need at least two jets to be able to compute MR and R2

  TH1D* h_MR_Cleaning = new TH1D("h_MR_Cleaning", "h_MR_Cleaning", 20, 0, MRmx);
  TH1D* h_R2_Cleaning = new TH1D("h_R2_Cleaning", "h_R2_Cleaning", 25, 0, 1);
  TH2D* h_MR_R2_Cleaning = new TH2D("h_MR_R2_Cleaning", "h_MR_R2_Cleaning", 20, 0, MRmx, 25, 0, 1);

  TH1D* h_MR_HCAL_noise = new TH1D("h_MR_HCAL_noise", "h_MR_HCAL_noise", 20, 0, MRmx);
  TH1D* h_R2_HCAL_noise = new TH1D("h_R2_HCAL_noise", "h_R2_HCAL_noise", 25, 0, 1);
  TH2D* h_MR_R2_HCAL_noise = new TH2D("h_MR_R2_HCAL_noise", "h_MR_R2_HCAL_noise", 20, 0, MRmx, 25, 0, 1);

  TH1D* h_MR_vertexg0 = new TH1D("h_MR_vertexg0", "h_MR_vertexg0", 20, 0, MRmx);
  TH1D* h_R2_vertexg0 = new TH1D("h_R2_vertexg0", "h_R2_vertexg0", 25, 0, 1);
  TH2D* h_MR_R2_vertexg0 = new TH2D("h_MR_R2_vertexg0", "h_MR_R2_vertexg0", 20, 0, MRmx, 25, 0, 1);

  TH1D* h_MR_njetge3 = new TH1D("h_MR_njetge3", "h_MR_njetge3", 20, 0, MRmx);
  TH1D* h_R2_njetge3 = new TH1D("h_R2_njetge3", "h_R2_njetge3", 25, 0, 1);
  TH2D* h_MR_R2_njetge3 = new TH2D("h_MR_R2_njetge3", "h_MR_R2_njetge3", 20, 0, MRmx, 25, 0, 1);

  TH1D* h_MR_HLT = new TH1D("h_MR_HLT", "h_MR_HLT", 20, 0, MRmx);
  TH1D* h_R2_HLT = new TH1D("h_R2_HLT", "h_R2_HLT", 25, 0, 1);
  TH2D* h_MR_R2_HLT = new TH2D("h_MR_R2_HLT", "h_MR_R2_HLT", 20, 0, MRmx, 25, 0, 1);

  TH1D* h_MR_jet1ptg200 = new TH1D("h_MR_jet1ptg200", "h_MR_jet1ptg200", 20, 0, MRmx);
  TH1D* h_R2_jet1ptg200 = new TH1D("h_R2_jet1ptg200", "h_R2_jet1ptg200", 25, 0, 1);
  TH2D* h_MR_R2_jet1ptg200 = new TH2D("h_MR_R2_jet1ptg200", "h_MR_R2_jet1ptg200", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_MR_SIG = new TH1D("h_MR_SIG", "h_MR_SIG", 20, 0, MRmx);
  TH1D * h_R2_SIG = new TH1D("h_R2_SIG", "h_R2_SIG", 25, 0, 1);
  TH2D * h_MR_R2_SIG = new TH2D("h_MR_R2_SIG", "h_MR_R2_SIG", 20, 0, MRmx, 25, 0, 1);

  // 0 lepton trajectory
  TH1D * h_MR_neleeq0 = new TH1D("h_MR_neleeq0", "h_MR_neleeq0", 20, 0, MRmx);
  TH1D * h_R2_neleeq0 = new TH1D("h_R2_neleeq0", "h_R2_neleeq0", 25, 0, 1);
  TH2D * h_MR_R2_neleeq0 = new TH2D("h_MR_R2_neleeq0", "h_MR_R2_neleeq0", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_MR_nmueq0 = new TH1D("h_MR_nmueq0", "h_MR_nmueq0", 20, 0, MRmx);
  TH1D * h_R2_nmueq0 = new TH1D("h_R2_nmueq0", "h_R2_nmueq0", 25, 0, 1);
  TH2D * h_MR_R2_nmueq0 = new TH2D("h_MR_R2_nmueq0", "h_MR_R2_nmueq0", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_MR_trackIso = new TH1D("h_MR_trackIso", "h_MR_trackIso", 20, 0, MRmx);
  TH1D * h_R2_trackIso = new TH1D("h_R2_trackIso", "h_R2_trackIso", 25, 0, 1);
  TH2D * h_MR_R2_trackIso = new TH2D("h_MR_R2_trackIso", "h_MR_R2_trackIso", 20, 0, MRmx, 25, 0, 1);

  // g1Mb 0Ll
  TH1D * h_MR_g1Mb0Ll = new TH1D("h_MR_g1Mb0Ll", "h_MR_g1Mb0Ll", 20, 0, MRmx);
  TH1D * h_R2_g1Mb0Ll = new TH1D("h_R2_g1Mb0Ll", "h_R2_g1Mb0Ll", 25, 0, 1);
  TH2D * h_MR_R2_g1Mb0Ll = new TH2D("h_MR_R2_g1Mb0Ll", "h_MR_R2_g1Mb0Ll", 20, 0, MRmx, 25, 0, 1);

  // g1Mb g1W 0Ll ; Signal box: >= 1 Mb; >= 1 W; 0 Ll
  TH1D * h_MR_g1Mbg1W0Ll = new TH1D("h_MR_g1Mbg1W0Ll", "h_MR_g1Mbg1W0Ll", 20, 0, MRmx);
  TH1D * h_R2_g1Mbg1W0Ll = new TH1D("h_R2_g1Mbg1W0Ll", "h_R2_g1Mbg1W0Ll", 25, 0, 1);
  TH2D * h_MR_R2_g1Mbg1W0Ll = new TH2D("h_MR_R2_g1Mbg1W0Ll", "h_MR_R2_g1Mbg1W0Ll", 20, 0, MRmx, 25, 0, 1);
  TH1D * h_njets_g1Mbg1W0Ll = new TH1D("h_njets_g1Mbg1W0Ll","h_njets_g1Mbg1W0Ll",10,0,10);

  // 0Lb 0Ll
  TH1D * h_MR_0Lb0Ll = new TH1D("h_MR_0Lb0Ll", "h_MR_0Lb0Ll", 20, 0, MRmx);
  TH1D * h_R2_0Lb0Ll = new TH1D("h_R2_0Lb0Ll", "h_R2_0Lb0Ll", 25, 0, 1);
  TH2D * h_MR_R2_0Lb0Ll = new TH2D("h_MR_R2_0Lb0Ll", "h_MR_R2_0Lb0Ll", 20, 0, MRmx, 25, 0, 1);

  // QCD control region: 0 Lb; >= 1 uW; 0 Ll
  TH1D * h_MR_0Lbg1uW0Ll = new TH1D("h_MR_0Lbg1uW0Ll", "h_MR_0Lbg1uW0Ll", 20, 0, MRmx);
  TH1D * h_R2_0Lbg1uW0Ll = new TH1D("h_R2_0Lbg1uW0Ll", "h_R2_0Lbg1uW0Ll", 25, 0, 1);
  TH2D * h_MR_R2_0Lbg1uW0Ll = new TH2D("h_MR_R2_0Lbg1uW0Ll", "h_MR_R2_0Lbg1uW0Ll", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_minDeltaPhi_0Lbg1uW0Ll = new TH1D("h_minDeltaPhi_0Lbg1uW0Ll", "h_minDeltaPhi_0Lbg1uW0Ll", 50, 0, 5);
  TH2D * h_MR_minDeltaPhi_0Lbg1uW0Ll = new TH2D("h_MR_minDeltaPhi_0Lbg1uW0Ll", "h_MR_minDeltaPhi_0Lbg1uW0Ll", 20, 0, MRmx, 50, 0, 5);
  TH2D * h_R2_minDeltaPhi_0Lbg1uW0Ll = new TH2D("h_R2_minDeltaPhi_0Lbg1uW0Ll", "h_R2_minDeltaPhi_0Lbg1uW0Ll", 25, 0, 1, 50, 0, 5);
  
  // QCD control region: 0 Lb; >= 1 uW; 0 Ll + minDeltaPhi < 0.3
  TH1D * h_MR_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_MR_0Lbg1uW0Ll_mdPhi0p3", "h_MR_0Lbg1uW0Ll_mdPhi0p3", 20, 0, MRmx);
  TH1D * h_R2_0Lbg1uW0Ll_mdPhi0p3 = new TH1D("h_R2_0Lbg1uW0Ll_mdPhi0p3", "h_R2_0Lbg1uW0Ll_mdPhi0p3", 25, 0, 1);
  TH2D * h_MR_R2_0Lbg1uW0Ll_mdPhi0p3 = new TH2D("h_MR_R2_0Lbg1uW0Ll_mdPhi0p3", "h_MR_R2_0Lbg1uW0Ll_mdPhi0p3", 20, 0, MRmx, 25, 0, 1);

  // QCD control region: 0 Lb; >= 1 W; 0 Ll
  TH1D * h_MR_0Lbg1W0Ll = new TH1D("h_MR_0Lbg1W0Ll", "h_MR_0Lbg1W0Ll", 20, 0, MRmx);
  TH1D * h_R2_0Lbg1W0Ll = new TH1D("h_R2_0Lbg1W0Ll", "h_R2_0Lbg1W0Ll", 25, 0, 1);
  TH2D * h_MR_R2_0Lbg1W0Ll = new TH2D("h_MR_R2_0Lbg1W0Ll", "h_MR_R2_0Lbg1W0Ll", 20, 0, MRmx, 25, 0, 1);


  // 1 loose lepton trajectory
  TH1D * h_MR_1Ll = new TH1D("h_MR_1Ll", "h_MR_1Ll", 20, 0, MRmx);
  TH1D * h_R2_1Ll = new TH1D("h_R2_1Ll", "h_R2_1Ll", 25, 0, 1);
  TH2D * h_MR_R2_1Ll = new TH2D("h_MR_R2_1Ll", "h_MR_R2_1Ll", 20, 0, MRmx, 25, 0, 1);

  // g1Mb1Ll
  TH1D * h_MR_g1Mb1Ll = new TH1D("h_MR_g1Mb1Ll", "h_MR_g1Mb1Ll", 20, 0, MRmx);
  TH1D * h_R2_g1Mb1Ll = new TH1D("h_R2_g1Mb1Ll", "h_R2_g1Mb1Ll", 25, 0, 1);
  TH2D * h_MR_R2_g1Mb1Ll = new TH2D("h_MR_R2_g1Mb1Ll", "h_MR_R2_g1Mb1Ll", 20, 0, MRmx, 25, 0, 1);

  // g1Mbg1W1Ll ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll
  TH1D * h_MR_g1Mbg1W1Ll = new TH1D("h_MR_g1Mbg1W1Ll", "h_MR_g1Mbg1W1Ll", 20, 0, MRmx);
  TH1D * h_R2_g1Mbg1W1Ll = new TH1D("h_R2_g1Mbg1W1Ll", "h_R2_g1Mbg1W1Ll", 25, 0, 1);
  TH2D * h_MR_R2_g1Mbg1W1Ll = new TH2D("h_MR_R2_g1Mbg1W1Ll", "h_MR_R2_g1Mbg1W1Ll", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_mT_g1Mbg1W1Ll = new TH1D("h_mT_g1Mbg1W1Ll", "h_mT_g1Mbg1W1Ll", 50, 0, 500);
  TH2D * h_MR_mT_g1Mbg1W1Ll = new TH2D("h_MR_mT_g1Mbg1W1Ll", "h_MR_mT_g1Mbg1W1Ll", 20, 0, MRmx, 50, 0, 500);
  TH2D * h_R2_mT_g1Mbg1W1Ll = new TH2D("h_R2_mT_g1Mbg1W1Ll", "h_R2_mT_g1Mbg1W1Ll", 25, 0, 1, 50, 0, 500);

  // g1Mbg1W1LlmT ; TTj control region: >= 1 Mb; >= 1 W; 1 Ll; 30<mT<100
  TH1D * h_MR_g1Mbg1W1LlmT = new TH1D("h_MR_g1Mbg1W1LlmT", "h_MR_g1Mbg1W1LlmT", 20, 0, MRmx);
  TH1D * h_R2_g1Mbg1W1LlmT = new TH1D("h_R2_g1Mbg1W1LlmT", "h_R2_g1Mbg1W1LlmT", 25, 0, 1);
  TH2D * h_MR_R2_g1Mbg1W1LlmT = new TH2D("h_MR_R2_g1Mbg1W1LlmT", "h_MR_R2_g1Mbg1W1LlmT", 20, 0, MRmx, 25, 0, 1);

  TH1D * h_mT_g1Mbg1W1LlmT = new TH1D("h_mT_g1Mbg1W1LlmT", "h_mT_g1Mbg1W1LlmT", 50, 0, 500);
  TH2D * h_MR_mT_g1Mbg1W1LlmT = new TH2D("h_MR_mT_g1Mbg1W1LlmT", "h_MR_mT_g1Mbg1W1LlmT", 20, 0, MRmx, 50, 0, 500);
  TH2D * h_R2_mT_g1Mbg1W1LlmT = new TH2D("h_R2_mT_g1Mbg1W1LlmT", "h_R2_mT_g1Mbg1W1LlmT", 25, 0, 1, 50, 0, 500);


  // dimuon trajectory
  TH1D * h_MR_2mu = new TH1D("h_MR_2mu", "h_MR_2mu", 20, 0, MRmx);
  TH1D * h_R2_2mu = new TH1D("h_R2_2mu", "h_R2_2mu", 25, 0, 1);
  TH2D * h_MR_R2_2mu = new TH2D("h_MR_R2_2mu", "h_MR_R2_2mu", 20, 0, MRmx, 25, 0, 1);
  
  TH1D * h_MR_2mu0el = new TH1D("h_MR_2mu0el", "h_MR_2mu0el", 20, 0, MRmx);
  TH1D * h_R2_2mu0el = new TH1D("h_R2_2mu0el", "h_R2_2mu0el", 25, 0, 1);
  TH2D * h_MR_R2_2mu0el = new TH2D("h_MR_R2_2mu0el", "h_MR_R2_2mu0el", 20, 0, MRmx, 25, 0, 1);
  
  TH1D * h_MR_0Lb2mu0el = new TH1D("h_MR_0Lb2mu0el", "h_MR_0Lb2mu0el", 20, 0, MRmx);
  TH1D * h_R2_0Lb2mu0el = new TH1D("h_R2_0Lb2mu0el", "h_R2_0Lb2mu0el", 25, 0, 1);
  TH2D * h_MR_R2_0Lb2mu0el = new TH2D("h_MR_R2_0Lb2mu0el", "h_MR_R2_0Lb2mu0el", 20, 0, MRmx, 25, 0, 1);
  
  TH1D * h_MR_g1Mb2mu0el = new TH1D("h_MR_g1Mb2mu0el", "h_MR_g1Mb2mu0el", 20, 0, MRmx);
  TH1D * h_R2_g1Mb2mu0el = new TH1D("h_R2_g1Mb2mu0el", "h_R2_g1Mb2mu0el", 25, 0, 1);
  TH2D * h_MR_R2_g1Mb2mu0el = new TH2D("h_MR_R2_g1Mb2mu0el", "h_MR_R2_g1Mb2mu0el", 20, 0, MRmx, 25, 0, 1);
  

  // Gen level plots
  TH1D* h_gen_toppt = new TH1D("h_gen_toppt", "h_gen_toppt", 50, 0, 1000);
  TH1D* h_gen_dRWb = new TH1D("h_gen_dRWb", "h_gen_dRWb", 200, 0, 5);
  TH2D* h_gen_toppt_dRWb = new TH2D("h_gen_toppt_dRWb", "h_gen_toppt_dRWb", 200, 0, 1000, 200, 0, 5);
  TH1D* h_gen_Wpt = new TH1D("h_gen_Wpt", "h_gen_Wpt", 50, 0, 1000);
  TH1D* h_gen_dRqq = new TH1D("h_gen_dRqq", "h_gen_dRqq", 200, 0, 5);
  TH2D* h_gen_Wpt_dRqq = new TH2D("h_gen_Wpt_dRqq", "h_gen_Wpt_dRqq", 200, 0, 1000, 200, 0, 5);

  // Define the order of bins in the counts histogram:
  
  ofile.count("NoCuts", 0.0);
  ofile.count("Cleaning", 0.0);
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
  ofile.count("0Lb0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll", 0.0);
  ofile.count("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  ofile.count("0Lbg1W0Ll", 0.0);

  ofile.count("1Ll", 0.0);
  ofile.count("g1Mb1Ll", 0.0);
  ofile.count("g1Mbg1W1Ll", 0.0);
  ofile.count("g1Mbg1W1LlmT", 0.0);

  ofile.count("2mu", 0.0);
  ofile.count("2mu0el", 0.0);
  ofile.count("0Lb2mu0el", 0.0);
  ofile.count("g1Mb2mu0el", 0.0);
  
  TH1D* TTallhad = new TH1D("counts_TTallhad","",1,0,1);
  TTallhad->SetBit(TH1::kCanRebin);
  TH1D* TTsemilep = new TH1D("counts_TTsemilep","",1,0,1);
  TTsemilep->SetBit(TH1::kCanRebin);
  TH1D* TTdilep = new TH1D("counts_TTdilep","",1,0,1);
  TTdilep->SetBit(TH1::kCanRebin);
 
  TTallhad->Fill("NoCuts", 0.0);
  TTallhad->Fill("Cleaning", 0.0);
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
  TTallhad->Fill("0Lb0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll", 0.0);
  TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTallhad->Fill("0Lbg1W0Ll", 0.0);
  TTallhad->Fill("1Ll", 0.0);
  TTallhad->Fill("g1Mb1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1Ll", 0.0);
  TTallhad->Fill("g1Mbg1W1LlmT", 0.0);
  TTallhad->Fill("2mu", 0.0);
  TTallhad->Fill("2mu0el", 0.0);
  TTallhad->Fill("0Lb2mu0el", 0.0);
  TTallhad->Fill("g1Mb2mu0el", 0.0);

  TTsemilep->Fill("NoCuts", 0.0);
  TTsemilep->Fill("Cleaning", 0.0);
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
  TTsemilep->Fill("0Lb0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll", 0.0);
  TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTsemilep->Fill("0Lbg1W0Ll", 0.0);
  TTsemilep->Fill("1Ll", 0.0);
  TTsemilep->Fill("g1Mb1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1Ll", 0.0);
  TTsemilep->Fill("g1Mbg1W1LlmT", 0.0);
  TTsemilep->Fill("2mu", 0.0);
  TTsemilep->Fill("2mu0el", 0.0);
  TTsemilep->Fill("0Lb2mu0el", 0.0);
  TTsemilep->Fill("g1Mb2mu0el", 0.0);

  TTdilep->Fill("NoCuts", 0.0);
  TTdilep->Fill("Cleaning", 0.0);
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
  TTdilep->Fill("0Lb0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll", 0.0);
  TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p3", 0.0);
  TTdilep->Fill("0Lbg1W0Ll", 0.0);
  TTdilep->Fill("1Ll", 0.0);
  TTdilep->Fill("g1Mb1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1Ll", 0.0);
  TTdilep->Fill("g1Mbg1W1LlmT", 0.0);
  TTdilep->Fill("2mu", 0.0);
  TTdilep->Fill("2mu0el", 0.0);
  TTdilep->Fill("0Lb2mu0el", 0.0);
  TTdilep->Fill("g1Mb2mu0el", 0.0);



  TH1::SetDefaultSumw2();



  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  //nevents = 10000;
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

      // ---------------------------------------------
      // -- First look at generator level particles --
      // ---------------------------------------------

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
      std::vector<TLorentzVector> sjetl;
      std::vector<cmgpfjet_s> sbjet;
      std::vector<cmgpfjet_s> slbjet;
      for (unsigned int i=0; i<cmgpfjet.size(); i++) {
	if (!(cmgpfjet[i].pt > 30) ) continue;
	if (!(fabs(cmgpfjet[i].eta) < 3) ) continue;
	//if (!(cmgpfjet[i].neutralHadronEnergyFraction < 0.99) ) continue;
	if (!(cmgpfjet[i].component_5_fraction + cmgpfjet[i].component_6_fraction < 0.99) ) continue;
	//if (!(cmgpfjet[i].neutralEmEnergyFraction < 0.99) ) continue;
	if (!(cmgpfjet[i].component_4_fraction < 0.99) ) continue;
	if (!(cmgpfjet[i].nConstituents > 1) ) continue;
	if (fabs(cmgpfjet[i].eta) < 2.4) {
	  if (!(cmgpfjet[i].component_1_fraction > 0) ) continue;
	  if (!(cmgpfjet[i].component_1_number > 0) ) continue;
	  if (!(cmgpfjet[i].component_2_fraction < 0.99) ) continue;
	}
	sjet.push_back(cmgpfjet[i]);
	if (cmgpfjet[i].combinedSecondaryVertexBJetTags > 0.679) {
	  sbjet.push_back(cmgpfjet[i]);
	}
	if (cmgpfjet[i].combinedSecondaryVertexBJetTags > 0.244) {
	  slbjet.push_back(cmgpfjet[i]);
	}
	TLorentzVector jl;
	jl.SetPtEtaPhiE(cmgpfjet[i].pt, cmgpfjet[i].eta,
			cmgpfjet[i].phi, cmgpfjet[i].energy);
	sjetl.push_back(jl);
      }


      // CA8
      // W selection:
      std::vector<jethelper4_s> sjet2;
      std::vector<jethelper4_s> sW;
      std::vector<jethelper4_s> aW;
      for (unsigned int i=0; i<jethelper4.size(); i++) {
        if (!(jethelper4[i].pt > 30) ) continue;
        if (!(fabs(jethelper4[i].eta) < 3) ) continue;

        h_jmass_jpt->Fill(jethelper4[i].mass, jethelper4[i].pt);
        h_d1pt_d2pt->Fill(jethelper4[i].daughter_0_pt, jethelper4[i].daughter_1_pt);
        h_d1m_d2m->Fill(jethelper4[i].daughter_0_mass, jethelper4[i].daughter_1_mass);
        h_jmass->Fill(jethelper4[i].mass);
	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
        h_d1ptsel_d2ptsel->Fill(jethelper4[i].daughter_0_pt, jethelper4[i].daughter_1_pt);
        h_d1msel_d2msel->Fill(jethelper4[i].daughter_0_mass, jethelper4[i].daughter_1_mass);
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
      for (unsigned int i=0; i<cmgmuon1.size(); i++) {
	if (!(cmgmuon1[i].pt > 25) ) continue;
	if (!(fabs(cmgmuon1[i].eta) < 2.4) ) continue;
	smuon.push_back(cmgmuon1[i]);
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
      // Electrons - tight
      std::vector<cmgelectron1_s> selectron;
      for (unsigned int i=0; i<cmgelectron1.size(); i++) {
        if (!(cmgelectron1[i].pt > 30) ) continue;
        if (!(fabs(cmgelectron1[i].eta) < 2.5) ) continue;
        // veto 1.442 < |eta| < 1.556? --> is already done in the collection ??
        if (!(fabs(cmgelectron1[i].eta) < 1.442 && fabs(cmgelectron1[i].eta) < 1.556) ) continue;
        selectron.push_back(cmgelectron1[i]);
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
      double nbs = sbjet.size();

      h_nW_Cleaning->Fill(nWs, w);
      h_nb_Cleaning->Fill(nbs, w);
      h_nW_nb_Cleaning->Fill(nWs, nbs, w);

      h_nWAK5_Cleaning->Fill(sWAK5.size(), w);
      
      // ---------------------
      // -- Razor variables --
      // ---------------------


      TVector3 metl;
      metl.SetPtEtaPhi(cmgbasemet2[0].et, 0, cmgbasemet2[0].phi);
      TLorentzVector met;
      met.SetPtEtaPhiE(cmgbasemet2[0].pt, 0, cmgbasemet2[0].phi, cmgbasemet2[0].energy);
      //cout << "MET info: " << cmgbasemet2[0].pt << " " << cmgbasemet2[0].phi << " " << cmgbasemet2[0].energy << endl;
      std::vector<TLorentzVector> hemis = CombineJets(sjetl);

      double MR = -9999;
      double MTR = -9999;
      double R2 = -9999;
      if (hemis.size() == 2) {
	MR = CalcMR(hemis[0], hemis[1]);
	if (MR == MR) {
	  h_MR_Cleaning->Fill(MR, w);
	  MTR = CalcMTR(hemis[0], hemis[1], metl);
	  R2 = pow((MTR / MR),2);
	  h_R2_Cleaning->Fill(R2, w);
	  h_MR_R2_Cleaning->Fill(MR, R2, w);
	  h_met_R2_Cleaning->Fill(cmgbasemet2[0].et, R2, w);
	}
      }
      //cout << MR << " " << MTR << " " << R2 << endl;


      // ---------------------
      // -- fill histograms --
      // ---------------------


      h_met_Cleaning->Fill(cmgbasemet2[0].et, w);
      if (sjet.size() > 0) {
        h_j1pt_Cleaning->Fill(sjet[0].pt, w);
      }
      h_nj_Cleaning->Fill(sjet.size(), w);

      // ---------------------
      // -- event selection --
      // ---------------------

      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;
      ofile.count("HCAL_noise", w);
      h_MR_HCAL_noise->Fill(MR, w);
      h_R2_HCAL_noise->Fill(R2, w);
      h_MR_R2_HCAL_noise->Fill(MR, R2, w);
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
      if(isTTallhad)
	TTallhad->Fill("njetge3", w);
      else if(isTTsemilep)
	TTsemilep->Fill("njetge3", w);
      else if(isTTdilep)
	TTdilep->Fill("njetge3", w);
      
      // Calculate the HLT weight and include it in the total weight:
      double whlt = 1;
      if (eventhelper_isRealData==0) {
	for (int i=1; i<h_hlteff->GetNbinsX()+1; i++) {
	  double xmin = h_hlteff->GetXaxis()->GetBinLowEdge(i);
	  double xmax = h_hlteff->GetXaxis()->GetBinUpEdge(i);
	  if (!(MR >= xmin && MR < xmax)) continue;
	  for (int j=1; j<h_hlteff->GetNbinsY()+1; j++) {
	    double ymin = h_hlteff->GetYaxis()->GetBinLowEdge(j);
	    double ymax = h_hlteff->GetYaxis()->GetBinUpEdge(j);
	    if (R2 >= ymin && R2 < ymax) {
	      whlt = h_hlteff->GetBinContent(i, j);
	      //cout << xmin << " " << MR << " " << xmax << " " << ymin << " " << R2 << " " << ymax << " " << whlt << endl;
	      break;
	    }
	  }
	}
	
	w = w*whlt;
      }
      //w = w*whlt;
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

      // pt of first jet greater than 200 GeV
      if (!(sjet[0].pt > 200)) continue;
      ofile.count("jet1ptg200", w);
      h_MR_jet1ptg200->Fill(MR, w);
      h_R2_jet1ptg200->Fill(R2, w);
      h_MR_R2_jet1ptg200->Fill(MR, R2, w);
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


      // Only select events in MR-R2 SIG region 
      if (!(MR > 800 && R2 > 0.08)) continue;
      ofile.count("SIG", w);
      h_MR_SIG->Fill(MR, w);
      h_R2_SIG->Fill(R2, w);
      h_MR_R2_SIG->Fill(MR, R2, w);
      if(isTTallhad)
	TTallhad->Fill("SIG", w);
      else if(isTTsemilep)
	TTsemilep->Fill("SIG", w);
      else if(isTTdilep)
	TTdilep->Fill("SIG", w);


      // Compute the minDeltaPhi variable, taking the first three jets into account
      double minDeltaPhi = 99.;
      for (int jet=0; jet<3; ++jet){
	double mdphi = fdeltaPhi(sjet[jet].phi,metl.Phi());
	if (mdphi < minDeltaPhi)
	  minDeltaPhi = mdphi;
      }


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

      // 0 Lepton
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

	    if (nmediumbs > 0){
	      ofile.count("g1Mb0Ll", w);
	      h_MR_g1Mb0Ll->Fill(MR, w);
	      h_R2_g1Mb0Ll->Fill(R2, w);
	      h_MR_R2_g1Mb0Ll->Fill(MR, R2, w);
	      if(isTTallhad)
		TTallhad->Fill("g1Mb0Ll", w);
	      else if(isTTsemilep)
		TTsemilep->Fill("g1Mb0Ll", w);
	      else if(isTTdilep)
		TTdilep->Fill("g1Mb0Ll", w);
	      
	      // g1Mb g1W 0Ll
	      if( sW.size() > 0){
		ofile.count("g1Mbg1W0Ll",w);
		h_MR_g1Mbg1W0Ll->Fill(MR, w);
		h_R2_g1Mbg1W0Ll->Fill(R2, w);
		h_MR_R2_g1Mbg1W0Ll->Fill(MR, R2, w);
		h_njets_g1Mbg1W0Ll->Fill(sjet.size(),w);
		if(isTTallhad)
		  TTallhad->Fill("g1Mbg1W0Ll", w);
		else if(isTTsemilep)
		  TTsemilep->Fill("g1Mbg1W0Ll", w);
		else if(isTTdilep)
		  TTdilep->Fill("g1Mbg1W0Ll", w);
		
	      } // end of sW.size() > 0
	    } // end of nmediumbs > 0

	    if (nloosebs == 0){
	      ofile.count("0Lb0Ll", w);
	      h_MR_0Lb0Ll->Fill(MR, w);
	      h_R2_0Lb0Ll->Fill(R2, w);
	      h_MR_R2_0Lb0Ll->Fill(MR, R2, w);
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
		h_minDeltaPhi_0Lbg1uW0Ll->Fill(minDeltaPhi, w);
		h_MR_minDeltaPhi_0Lbg1uW0Ll->Fill(MR, minDeltaPhi, w);
		h_R2_minDeltaPhi_0Lbg1uW0Ll->Fill(R2, minDeltaPhi, w);
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
		  if(isTTallhad)
		    TTallhad->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		  else if(isTTsemilep)
		    TTsemilep->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		  else if(isTTdilep)
		    TTdilep->Fill("0Lbg1uW0Ll_mdPhi0p3", w);
		} // end of minDeltaPhi < 0.3
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
      
      
      // 1 Loose, not tight lepton 
      if(nlooseleptons == 1 && ntightleptons == 0){
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
	if(isTTallhad)
	  TTallhad->Fill("1Ll", w);
	else if(isTTsemilep)
	  TTsemilep->Fill("1Ll", w);
	else if(isTTdilep)
	  TTdilep->Fill("1Ll", w);
	
	if (nmediumbs > 0){
	  ofile.count("g1Mb1Ll",w);
	  h_MR_g1Mb1Ll->Fill(MR, w);
	  h_R2_g1Mb1Ll->Fill(R2, w);
	  h_MR_R2_g1Mb1Ll->Fill(MR, R2, w);
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
	    if(isTTallhad)
	      TTallhad->Fill("g1Mbg1W1Ll", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("g1Mbg1W1Ll", w);
	    else if(isTTdilep)
	      TTdilep->Fill("g1Mbg1W1Ll", w);

	    // mT window
	    if (mT > 30 && mT < 100){
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

	    }
	  } // end sW.size()
	} // end nmediumbs > 0
      } // end nlooseleptons == 1

      if (ntightmuons == 2 && nloosemuons == 2) {
	ofile.count("2mu",w);
	h_MR_2mu->Fill(MR, w);
	h_R2_2mu->Fill(R2, w);
	h_MR_R2_2mu->Fill(MR, R2, w);
	if(isTTallhad)
	  TTallhad->Fill("2mu", w);
	else if(isTTsemilep)
	  TTsemilep->Fill("2mu", w);
	else if(isTTdilep)
	  TTdilep->Fill("2mu", w);

	if (nlooseelectrons == 0){
	  ofile.count("2mu0el",w);
	  h_MR_2mu0el->Fill(MR, w);
	  h_R2_2mu0el->Fill(R2, w);
	  h_MR_R2_2mu0el->Fill(MR, R2, w);
	  if(isTTallhad)
	    TTallhad->Fill("2mu0el", w);
	  else if(isTTsemilep)
	    TTsemilep->Fill("2mu0el", w);
	  else if(isTTdilep)
	    TTdilep->Fill("2mu0el", w);

	  if (nloosebs == 0){
	    ofile.count("0Lb2mu0el",w);
	    h_MR_0Lb2mu0el->Fill(MR, w);
	    h_R2_0Lb2mu0el->Fill(R2, w);
	    h_MR_R2_0Lb2mu0el->Fill(MR, R2, w);
	    if(isTTallhad)
	      TTallhad->Fill("0Lb2mu0el", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("0Lb2mu0el", w);
	    else if(isTTdilep)
	      TTdilep->Fill("0Lb2mu0el", w);
	  }// end nloosebs == 0

	  if (nmediumbs >= 1){
	    ofile.count("g1Mb2mu0el",w);
	    h_MR_g1Mb2mu0el->Fill(MR, w);
	    h_R2_g1Mb2mu0el->Fill(R2, w);
	    h_MR_R2_g1Mb2mu0el->Fill(MR, R2, w);
	    if(isTTallhad)
	      TTallhad->Fill("g1Mb2mu0el", w);
	    else if(isTTsemilep)
	      TTsemilep->Fill("g1Mb2mu0el", w);
	    else if(isTTdilep)
	      TTdilep->Fill("g1Mb2mu0el", w);
	  } // end nmediumbs >= 1
	} // end nlooseelectrons == 0
      } // end ntightmuons == 2
      
    } // end for

  fhlt->Close();
  stream.close();
  ofile.close();
  return 0;
}
