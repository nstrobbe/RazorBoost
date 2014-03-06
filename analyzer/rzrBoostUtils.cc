//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmd.h"
#include "utils.h"
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
  TFile* fhlt = TFile::Open("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/hlteff/extr_eff0_sm2.root");
  if (!fhlt){
    fhlt = TFile::Open("/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/hlteff/extr_eff0_sm2.root");
  }
  if (!fhlt){
    cout << "Could not find trigger efficiency root file... Where did you put it??" << endl;
    return 1;
  }
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
    cout << "Will do top pt reweighting" << endl;
  }

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

  // Total weight
  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);
  TH1D* h_totweight_filter = new TH1D("h_totweight_filter", "h_totweight_filter", 1, 1, 2);

  // histograms for total ISR weight: Sum(all events) w_ISR
  // needed to do the reweighting properly without changing the overall cross section. 
  TH1D* h_totalISRweight_nominal = new TH1D("h_totalISRweight_nominal", "h_totalISRweight_nominal", 1, 1, 2); // nominal ISR weight
  TH1D* h_totalISRweight_up = new TH1D("h_totalISRweight_up", "h_totalISRweight_up", 1, 1, 2); //  ISR up weight
  TH1D* h_totalISRweight_down = new TH1D("h_totalISRweight_down", "h_totalISRweight_down", 1, 1, 2); // ISR down weight

  // Total weight for the TopPt reweighting
  TH1D* h_totalTopPTweight_nominal = new TH1D("h_totalTopPTweight_nominal", "h_totalTopPTweight_nominal", 1, 1, 2); // nominal TopPt weight
  TH1D* h_totalTopPTweight_up = new TH1D("h_totalTopPTweight_up", "h_totalTopPTweight_up", 1, 1, 2); //  TopPT up weight
  TH1D* h_totalTopPTweight_down = new TH1D("h_totalTopPTweight_down", "h_totalTopPTweight_down", 1, 1, 2); // TopPt down weight

  // Total weight for the TopPt reweighting on top of the nominal ISR reweighting
  TH1D* h_total_ISRnominal_TopPTweight_nominal = new TH1D("h_total_ISRnominal_TopPTweight_nominal", "h_total_ISRnominal_TopPTweight_nominal", 1, 1, 2); // nominal ISR weight + TopPt nominal
  TH1D* h_total_ISRnominal_TopPTweight_up = new TH1D("h_total_ISRnominal_TopPTweight_up", "h_total_ISRnominal_TopPTweight_up", 1, 1, 2); //  ISR nominal weight + TopPt up
  TH1D* h_total_ISRnominal_TopPTweight_down = new TH1D("h_total_ISRnominal_TopPTweight_down", "h_total_ISRnominal_TopPTweight_down", 1, 1, 2); // ISR weight + TopPt down

  TH1D* h_pileup = new TH1D("h_pileup","h_pileup",1000,0,100);

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
        w = geneventinfoproduct_weight;
      }
      
      // Write every ith event:
      if (entry % 10000 == 0) cout << entry << endl;
      
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file rzrBTanalyzer.h to find out what structs
      // are available. Each struct contains the field "selected", which
      // can be set as needed. Call saveSelectedObjects() before a call to
      // addEvent if you wish to save only the selected objects.
      
      fillObjects();

      // Get rid of the events with wrong kinematics in the MadGraph samples
      //if (eventhelper_isRealData!=1) {
      //  if (triggerresultshelper2_totalKinematicsFilterPath==0) continue;
      //}
      h_totweight_filter->Fill(1,w);

      // Fill number of pileup interactions
      h_pileup->Fill(pileupsummaryinfo[0].getTrueNumInteractions);

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
      h_totalISRweight_nominal->Fill(1,w*w_ISR_nominal);
      h_totalISRweight_up->Fill(1,w*w_ISR_up);
      h_totalISRweight_down->Fill(1,w*w_ISR_down);

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
      h_totalTopPTweight_nominal->Fill(1,w*w_TopPt_nominal);
      h_totalTopPTweight_up->Fill(1,w*w_TopPt_up);
      h_totalTopPTweight_down->Fill(1,w*w_TopPt_down);

      h_total_ISRnominal_TopPTweight_nominal->Fill(1,w*w_TopPt_nominal*w_ISR_nominal);
      h_total_ISRnominal_TopPTweight_up->Fill(1,w*w_TopPt_up*w_ISR_nominal);
      h_total_ISRnominal_TopPTweight_down->Fill(1,w*w_TopPt_down*w_ISR_nominal);






    } // end event loop
  
  fhlt->Close();
  stream.close();
  ofile.close();
  return 0;
}
