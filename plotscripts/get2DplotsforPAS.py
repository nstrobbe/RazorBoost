import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':


    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140917/forPAS"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140917/summary/"
    inputfile_bg = basedir + "../rzrBoostMC_bg.root"
    inputfile_sig = basedir + "../rzrBoostMC_T1ttcc_1000_325_300.root"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/2DplotsPAS.root","RECREATE")
    infile_bg = TFile.Open(inputfile_bg)
    infile_sig = TFile.Open(inputfile_sig)

    # Integrated luminosity in fb-1s
    intlumi = 19.7 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    #############################################################
    ## compare bg with signal
    #############################################################

    # total BG
    hdict_bg = plotTools.ConstructHDict(infile_bg.Get("h_MR_R2_jet1ptg200"), 
                                        name="Background", 
                                        title="",
                                        xtitle="M_{R} (GeV)", ytitle="R^{2}",
                                        drawoption="colz", palette="SMS") 
    canvasname = "MR_R2_jet1ptg200_bg"
    plotTools.Plot2DPAS(hdict_bg,outputdir,outfile,cname=canvasname,scale="Yes",logscale=True,lumitext="Simulation")
            
    # total BG
    hdict_sig = plotTools.ConstructHDict(infile_sig.Get("h_MR_R2_jet1ptg200"), 
                                        name="Signal", 
                                        title="",
                                        xtitle="M_{R} (GeV)", ytitle="R^{2}",
                                        drawoption="colz", palette="SMS") 
    canvasname = "MR_R2_jet1ptg200_sig"
    plotTools.Plot2DPAS(hdict_sig,outputdir,outfile,cname=canvasname,scale="Yes",logscale=True,lumitext="Simulation")
            
    
    outfile.Close()
    infile_bg.Close()
