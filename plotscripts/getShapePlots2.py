import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "Run as: python %s <outputdir> <inputfile1> " % (sys.argv[0])
        sys.exit()
    outputdir = sys.argv[1]
    inputfile = sys.argv[2] # total bg histograms

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile = TFile.Open(inputfile)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR","R2"]

    # build hdictlist for TTJ control regions:
    for var in vars:
        hdict_TTJ = plotTools.ConstructHDict(infile.Get("h_"+var+"_g1Mbg1W1Ll"),
                                             name="TTJets CR", color=rt.kCyan+3,
                                             title="Shape comparison for TTJets Control Regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmtwindow = plotTools.ConstructHDict(infile.Get("h_"+var+"_g1Mbg1W1LlmT"),
                                             name="TTJets CR, 30 < mT < 100", color=rt.kCyan,
                                             title="Shape comparison for TTJets Control Regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_TTJ = plotTools.ConstructHDict(infile.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                             name="TTJets CR, mT < 100", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets Control Regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_TTJ,hdict_TTJmtwindow,hdict_TTJ]
        canvasname = var+"_comparison_TTJ_mT"
        rtitle = "#frac{mT cut}{no mT cut}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle)


    outfile.Close()
