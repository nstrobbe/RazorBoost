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

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    inputfile_1000_325_300_old = basedir + "rzrBoostMC_T1ttcc_325_300.root"
    inputfile_1000_325_300_new = basedir + "rzrBoostMC_T1ttcc_1000_325_300.root"

    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile_1000_325_300_old = TFile.Open(inputfile_1000_325_300_old)
    infile_1000_325_300_new = TFile.Open(inputfile_1000_325_300_new)



    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR","R2"]
    yt = ["Events/(100 GeV)","Events/(0.01)"]
    sf = [100,0.01]
    cuts = ["g1Mbg1W0Ll","g1Mbg1W1LlmT100"]
    # leg = plotTools.ConstructLDict(0.7,0.88,0.7,0.88,"")
    
    for cut in cuts:
        for var in vars:
            hdict_old = plotTools.ConstructHDict(infile_1000_325_300_old.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_325_300 OLD", color=rt.kRed,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Ref", xtitle=var, ytitle=yt[vars.index(var)])
            hdict_new = plotTools.ConstructHDict(infile_1000_325_300_new.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_325_300 NEW", color=rt.kBlue,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Yes", xtitle=var, ytitle=yt[vars.index(var)])
            
            hdictlist=[hdict_old,hdict_new]
            canvasname = "Signal_comparison_oldvsnew_1000_325_300__"+var+"__"+cut
            rtitle = "#frac{NEW}{OLD}"
            # plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="No",legdict=leg)
            plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Width",scalefactor=sf[vars.index(var)])

    vars = ["met","njets","nbjets","jet1pt","jet2pt","jet3pt"]
    cuts = ["g1Mbg1W0Ll","g1Mbg1W1LlmT100"]
    # leg = plotTools.ConstructLDict(0.7,0.88,0.7,0.88,"")
    
    for cut in cuts:
        for var in vars:
            hdict_old = plotTools.ConstructHDict(infile_1000_325_300_old.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_325_300 OLD", color=rt.kRed,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Ref", xtitle=var)
            hdict_new = plotTools.ConstructHDict(infile_1000_325_300_new.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_325_300 NEW", color=rt.kBlue,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Yes", xtitle=var)
            
            hdictlist=[hdict_old,hdict_new]
            canvasname = "Signal_comparison_oldvsnew_1000_325_300__"+var+"__"+cut
            rtitle = "#frac{NEW}{OLD}"
            # plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="No",legdict=leg)
            plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="No")

    outfile.Close()
    infile_1000_325_300_old.Close()
    infile_1000_325_300_new.Close()

