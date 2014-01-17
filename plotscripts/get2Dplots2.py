import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    #if len(sys.argv) < 3:
    #    print "Run as: python %s <outputdir> <inputfile1> " % (sys.argv[0])
    #    sys.exit()
    #outputdir = sys.argv[1]
    #inputfile = sys.argv[2] # total bg histograms

    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140114"
    outputdir = "./test_2D/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140114/summary/"
    inputfile_TTJ = basedir + "rzrBoostMC_TTJets.root"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/2Dplots2.root","RECREATE")
    infile_TTJ = TFile.Open(inputfile_TTJ)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR_R2"]

    # comparison plots to be made:
    # 1. comparison between process in signal region and relevant control region
    #    * QCD
    #    * TTJets
    #    * WJets
    #    * ZJets
    # 2. comparison between different control regions
    #    * different QCD regions
    #    * different TTjets regions

    # build hdictlist for TTJ:
    for var in vars:
        hdict_TTJ = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1Ll"),
                                             name="TTJets CR", title="R2 vs MR for TTJets in TTJets Contol region",
                                             xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                             drawoption="colz", palette="SMS")
        
        canvasname = var+"_TTJ"
        plotTools.Plot2D(hdict_TTJ,outputdir,outfile,cname=canvasname,scale="No",logscale=True)




    outfile.Close()
    infile_TTJ.Close()
