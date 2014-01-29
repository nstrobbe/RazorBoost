//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Sun Oct 28 00:36:05 2012 by mkntanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
//#include "rzrBTsmsskim.h"
#include "rzrBTanalyzercmd.h"

#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/pdg.h"
#else
#include "pdg.h"
#endif

using namespace std;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line

  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  double _mSTOP = atof(argv[6]);
  double _mLSP = atof(argv[7]);
  double _mG = atof(argv[8]);

  cout << "mstop, mLSP, mG: " << _mSTOP << ", " << _mLSP << ", " << _mG << endl;
  // Get names of ntuple files to be processed and open chain of ntuples

  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

  // Select variables to be read

  selectVariables(stream);

  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the following
  // line

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
  
  outputFile ofile(cmdline.outputfilename, stream);

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------

  TH1D* h_mstop = new TH1D("h_mstop","h_mstop",300,0,1500);
  TH1D* h_mLSP = new TH1D("h_mLSP","h_mLSP",300,0,1500);
  TH1D* h_mg = new TH1D("h_mg","h_mg",300,0,1500);
  TH2D* h_mstop_mLSP = new TH2D("h_mstop_mLSP","h_mstop_mLSP",300,0,1500,300,0,1500);
  TH2D* h_mg_mLSP = new TH2D("h_mg_mLSP","h_mg_mLSP",300,0,1500,300,0,1500);

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry)
  //for(int entry=0; entry<100; ++entry)
    {
      // Read event into memory
      stream.read(entry);

      if (entry % 100000 == 0) 
	cout << entry << endl;
      
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file rzrBTanalyzer.h to find out what structs
      // are available. Each struct contains the field "selected", which
      // can be set as needed. Call saveSelectedObjects() before a call to
      // addEvent if you wish to save only the selected objects.
      
      //fillObjects();

      double mt1 = lheeventproducthelper_mt1;
      double mz1 = lheeventproducthelper_mz1;
      double mg = lheeventproducthelper_mg;

      // pick some points, range available: Mstop-225to1200_mLSP-0to1000
      // chosen points: 600-100 ; 750-100 ; 750-300 ; 800-500
      //cout << "stop, lsp: " << mt1 << " " << mz1 << endl; 
      if (!( mt1 == _mSTOP && mz1 ==_mLSP && mg == _mG)) continue;      
      //cout << "Event passed" << endl;
      ofile.addEvent(1.0);
      h_mstop->Fill(mt1);
      h_mLSP->Fill(mz1);
      h_mg->Fill(mg);
      h_mstop_mLSP->Fill(mt1,mz1);
      h_mg_mLSP->Fill(mg,mz1);
      
      // ---------------------
      // -- event selection --
      // ---------------------
      
      
      // ---------------------
      // -- fill histograms --
      // ---------------------	  
      
    }
  
  ofile.close();
  stream.close();
  return 0;
}
