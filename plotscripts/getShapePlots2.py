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

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140101"
    inputfile_TTJ = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140101/summary/rzrBoostMC_TTJets.root"
    inputfile_QCD = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140101/summary/rzrBoostMC_QCD.root"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile_TTJ = TFile.Open(inputfile_TTJ)
    infile_QCD = TFile.Open(inputfile_QCD)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR","R2"]

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
                                             name="TTJets CR, no mT cut", color=rt.kCyan+3,
                                             title="Shape comparison for TTJets in the various TTJets Contol regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmtwindow = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT"),
                                             name="TTJets CR, 30 < mT < 100", color=rt.kCyan,
                                             title="Shape comparison for TTJets in the various TTJets Contol regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_TTJmt = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                             name="TTJets CR, mT < 100", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets in the various TTJets Contol regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_TTJ,hdict_TTJmt,hdict_TTJmtwindow]
        canvasname = var+"_comparison_TTJ_mT"
        rtitle = "#frac{mT cut}{no mT cut}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for TTJets in Signal and TTJ Control region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmt = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                             name="TTJets CR, mT < 100", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets in Signal and TTJ control region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_TTJmt]
        canvasname = var+"_comparison_TTJ_SIG"
        rtitle = "#frac{TTJ}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # build hdictlist for QCD:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for QCD in the Signal and various QCD Control regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_QCD_mdphi = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),
                                             name="QCD CR, minDeltaPhi < 0.3", color=rt.kMagenta,
                                             title="Shape comparison for QCD in the Signal and various QCD Control regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_QCD_mdphihat4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhiHat4"),
                                             name="QCD CR, minDeltaPhiHat < 4", color=rt.kMagenta+2,
                                             title="Shape comparison for QCD in the Signal and various QCD Control regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_QCD_mdphihat5 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhiHat5"),
                                             name="QCD CR, minDeltaPhiHat < 5", color=rt.kMagenta-9,
                                             title="Shape comparison for QCD in the Signal and various QCD Control regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_QCD_mdphi,hdict_QCD_mdphihat4,hdict_QCD_mdphihat5]
        canvasname = var+"_comparison_QCD"
        rtitle = "#frac{QCD}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")


    outfile.Close()
    infile_TTJ.Close()
