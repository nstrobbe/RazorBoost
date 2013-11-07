//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Jun 12 20:22:53 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
#include "rzrBTanalyzercmd.h"
#include "utils.h"
#include <iostream>

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

  // If running on SMS, will have to do some extra stuff
  string SMS = "";
  if (argc > 6){
    SMS = argv[6];
  }
  
  TFile* fstopxsect = new TFile("/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/smsinput/stop.root");
  TH1D* hstopxsect = (TH1D*)fstopxsect->Get("stop");

  string fcountname = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/smsinput/T2tt_counts.root";
  if (SMS == "T1ttcc")
    fcountname = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/smsinput/T1ttcc_counts.root";

  TFile* fsmscount = new TFile(fcountname.c_str());
  TH1D* h_mstop_mLSP_nevents = (TH1D*)fsmscount->Get("h_mstop_mLSP_nevents");



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
  
  ofstream ofile;
  ofile.open(cmdline.outputfilename.c_str());
  ofile << "Weight \t MR \t R2" << endl;

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

  // For SMSs we need to get the cross section and total number of events in the event loop (it depends on the mass point considered)
  if(SMS != "")
    weightnorm = lumi;

  cout << "lumi: " << lumi << endl;
  cout << "xsect: " << xsect << endl;
  cout << "totweight: " << totweight << endl;
  cout << "weightnorm: " << weightnorm << endl;

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------

  

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  //nevents = 10000;
  for(int entry=0; entry < nevents; ++entry)
    {
      // Read event into memory
      stream.read(entry);
      
      
      // Count events and get the total weight contibuted by the event
      double w = 1.;
      if (geneventinfoproduct_weight != 0) {
        w = geneventinfoproduct_weight*weightnorm;
      }

      // Set SMS weights
      double mt1 = lheeventproducthelper_mt1;
      double mz1 = lheeventproducthelper_mz1;
      
      // Do not run on all T2tt points
      if (SMS == "T2tt") {
	int mt1_int = static_cast<int>(mt1);
	int mz1_int = static_cast<int>(mz1);
	if (mz1_int == 0) continue;
	if (mt1_int < 500) continue;
	if (mt1_int % 50 != 0) continue;
	if (mz1_int % 50 != 0  && mz1_int != 1) continue;
      }
      if (mt1>0) {
	if(SMS == "T2tt"){
	  for (int b=1; b<hstopxsect->GetNbinsX()+1; b++) {
	    if (mt1 == hstopxsect->GetBinCenter(b)) {
	      xsect = hstopxsect->GetBinContent(b);
	      break;
	    }
	  }
	} else if (SMS == "T1ttcc"){
	  xsect = 0.0243547; // all points have mgluino = 1 TeV
	}
	// get the event count from the histogram
        double it = 1;
        double iz = 1;
	if (SMS == "T2tt"){
	  for (int b=1; b<h_mstop_mLSP_nevents->GetNbinsX()+1; ++b){
	    double xbinedge = h_mstop_mLSP_nevents->GetXaxis()->GetBinLowEdge(b);
	    if (mt1 == xbinedge) {
	      it = b;
	      break;
	    }
	  }
	} 
	if (SMS == "T1ttcc") {
	  // Something got screwed up in the histogram, need to work around it
	  for (int b=1; b<h_mstop_mLSP_nevents->GetNbinsX()+1; ++b){
	    double xbinedge = h_mstop_mLSP_nevents->GetXaxis()->GetBinLowEdge(b);
	    double xbinedgehigh = h_mstop_mLSP_nevents->GetXaxis()->GetBinUpEdge(b);
	    if (mt1 >= xbinedge && mt1 < xbinedgehigh) {
	      it = b;
	      break;
	    }
	  }
	}
	for (int b=1; b<h_mstop_mLSP_nevents->GetNbinsY()+1; ++b){
	  double ybinedge = h_mstop_mLSP_nevents->GetYaxis()->GetBinLowEdge(b);
	  if (mz1 == ybinedge || (mz1 == 1 && ybinedge == 0) ) {
	    iz = b;
	    break;
	  }
	}
        totweight = h_mstop_mLSP_nevents->GetBinContent(it, iz);
        w = (w*xsect)/totweight;
      }
      

      // Write every ith event:
      if (entry % 10000 == 0) cout << entry << endl;
      
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file rzrBTanalyzer.h to find out what structs
      // are available. Each struct contains the field "selected", which
      // can be set as needed. Call saveSelectedObjects() before a call to
      // addEvent if you wish to save only the selected objects.
      
      fillObjects();

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
	// New Andreas cuts:
        if (!(jethelper4[i].mass > 70 && jethelper4[i].mass < 100)) continue;
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
	  MTR = CalcMTR(hemis[0], hemis[1], metl);
	  R2 = pow((MTR / MR),2);
	}
      }
      //cout << MR << " " << MTR << " " << R2 << endl;


      // ---------------------
      // -- fill histograms --
      // ---------------------


      // ---------------------
      // -- event selection --
      // ---------------------

      // Additional HCAL noise cleaning
      double dphi_PF_CALO_met = fdeltaPhi(cmgbasemet2[0].phi,calomet[0].phi);
      if (fabs(dphi_PF_CALO_met - TMath::Pi()) < 1 ) continue;

      // at least one good primary vertex
      if (!(svertex.size() > 0)) continue;

      // at least three jets
      if (!(sjet.size() >= 3)) continue;
      
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

      // pt of first jet greater than 200 GeV
      if (!(sjet[0].pt > 200)) continue;

      // Only select events in MR-R2 FULL region 
      if (!(MR > 600 && R2 > 0.04)) continue;



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
		
	if (nloosemuons == 0) {
	  
	  if (eventhelperextra_trackIso == 0){

	    if (nmediumbs > 0){
	      
	      // g1Mb g1W 0Ll
	      if( sW.size() > 0){
		ofile << w << "\t" << MR << "\t" << R2 << endl;
	      } // end of sW.size() > 0
	    } // end of nmediumbs > 0

	  } // end veto iso track
	} // end veto muon
      }  // end veto electron
      
      

      
    } // end for

  fhlt->Close();
  stream.close();
  ofile.close();
  return 0;
}
