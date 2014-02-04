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

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140131_ttjsemilep/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140131_ttjsemilep/"
    inputfile_TTJ = basedir + "rzrBoostMC_ttj_semilep.root"

    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile_TTJ = TFile.Open(inputfile_TTJ)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR","R2"]

    # build hdictlist for TTJ:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for TTJets semileptonic in Signal and TTJ Control region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmt = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                             name="TTJets CR, mT < 100", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets semileptonic in Signal and TTJ control region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_TTJmt]
        canvasname = var+"_comparison_TTJsemilep_SIG"
        rtitle = "#frac{TTJ}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")


    for var in vars:
        hdict_matched = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1uW_matched"),
                                                 name="untagged W matched to gen W", color=rt.kBlack,
                                                 title="Comparison between untagged W's that are matched or unmatched to gen W's",
                                                 appear_in_ratio="Ref",
                                                 xtitle=var)
        
        hdict_unmatched = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1uW_unmatched"),
                                                   name="untagged W not matched to gen W", color=rt.kRed,
                                                   title="Comparison between untagged W's that are matched or unmatched to gen W's",
                                                   appear_in_ratio="Yes",
                                                   xtitle=var)
        hdictlist = [hdict_matched,hdict_unmatched]
        canvasname = var+"_comparison_matched_unmatched_W"
        rtitle = "#frac{not matched}{matched}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")


    outfile.Close()
    infile_TTJ.Close()

