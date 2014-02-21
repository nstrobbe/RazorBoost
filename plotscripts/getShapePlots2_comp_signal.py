import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

def runA():
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    inputfile_1000_325_300_old = basedir + "rzrBoostMC_T1ttcc_325_300.root"
    inputfile_1000_325_300_new = basedir + "rzrBoostMC_T1ttcc_1000_325_300.root"
    inputfile_1000_510_500_old = basedir + "rzrBoostMC_T1ttcc_510_500.root"
    inputfile_1000_510_500_new = basedir + "rzrBoostMC_T1ttcc_1000_510_500.root"

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile_1000_325_300_old = TFile.Open(inputfile_1000_325_300_old)
    infile_1000_325_300_new = TFile.Open(inputfile_1000_325_300_new)
    infile_1000_510_500_old = TFile.Open(inputfile_1000_510_500_old)
    infile_1000_510_500_new = TFile.Open(inputfile_1000_510_500_new)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    
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

    for cut in cuts:
        for var in vars:
            hdict_old = plotTools.ConstructHDict(infile_1000_510_500_old.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_510_500 OLD", color=rt.kRed,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Ref", xtitle=var, ytitle=yt[vars.index(var)])
            hdict_new = plotTools.ConstructHDict(infile_1000_510_500_new.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_510_500 NEW", color=rt.kBlue,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Yes", xtitle=var, ytitle=yt[vars.index(var)])
            
            hdictlist=[hdict_old,hdict_new]
            canvasname = "Signal_comparison_oldvsnew_1000_510_500__"+var+"__"+cut
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

    for cut in cuts:
        for var in vars:
            hdict_old = plotTools.ConstructHDict(infile_1000_510_500_old.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_510_500 OLD", color=rt.kRed,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Ref", xtitle=var)
            hdict_new = plotTools.ConstructHDict(infile_1000_510_500_new.Get("h_"+var+"_"+cut),
                                                 name="T1ttcc_1000_510_500 NEW", color=rt.kBlue,
                                                 title="Comparison between old and new signal scan, for selection "+cut,
                                                 appear_in_ratio="Yes", xtitle=var)
            
            hdictlist=[hdict_old,hdict_new]
            canvasname = "Signal_comparison_oldvsnew_1000_510_500__"+var+"__"+cut
            rtitle = "#frac{NEW}{OLD}"
            # plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="No",legdict=leg)
            plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="No")

    outfile.Close()
    infile_1000_325_300_old.Close()
    infile_1000_325_300_new.Close()
    infile_1000_510_500_old.Close()
    infile_1000_510_500_new.Close()

def runB():
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140214_signal_comparison/"
    inputfile_1200_110_100 = basedir + "rzrBoostMC_T1ttcc_1200_110_100.root"
    inputfile_1200_125_100 = basedir + "rzrBoostMC_T1ttcc_1200_125_100.root"
    inputfile_1200_180_100 = basedir + "rzrBoostMC_T1ttcc_1200_180_100.root"

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2B.root","RECREATE")
    infile_1200_110_100 = TFile.Open(inputfile_1200_110_100)
    infile_1200_125_100 = TFile.Open(inputfile_1200_125_100)
    infile_1200_180_100 = TFile.Open(inputfile_1200_180_100)


    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    
    vars = ["MR","R2"]
    yt = ["Events/(100 GeV)","Events/(0.01)"]
    sf = [100,0.01]
    cuts = ["SIG","g1Mbg1W0Ll","g1Mbg1W1LlmT100"]
    # leg = plotTools.ConstructLDict(0.7,0.88,0.7,0.88,"")
    
    for cut in cuts:
        for var in vars:
            hdict_10 = plotTools.ConstructHDict(infile_1200_110_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_110_100", color=rt.kCyan,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Ref", xtitle=var, ytitle=yt[vars.index(var)])
            hdict_25 = plotTools.ConstructHDict(infile_1200_125_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_125_100", color=rt.kCyan+2,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Yes", xtitle=var, ytitle=yt[vars.index(var)])
            hdict_80 = plotTools.ConstructHDict(infile_1200_180_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_180_100", color=rt.kCyan+4,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Yes", xtitle=var, ytitle=yt[vars.index(var)])
            
            
            hdictlist=[hdict_10,hdict_25,hdict_80]
            canvasname = "Signal_comparison_DMlines_1200_1xx_100__"+var+"__"+cut
            rtitle = "#frac{DM}{DM10}"
            # plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="No",legdict=leg)
            plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Width",scalefactor=sf[vars.index(var)])



    vars = ["met","njets","nbjets","jet1pt","jet2pt","jet3pt"]
    cuts = ["g1Mbg1W0Ll","g1Mbg1W1LlmT100"]
    leg = plotTools.ConstructLDict(0.65,0.88,0.7,0.88,"")
    
    for cut in cuts:
        for var in vars:
            hdict_10 = plotTools.ConstructHDict(infile_1200_110_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_110_100", color=rt.kCyan,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Ref", xtitle=var)
            hdict_25 = plotTools.ConstructHDict(infile_1200_125_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_125_100", color=rt.kCyan+2,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Yes", xtitle=var)
            hdict_80 = plotTools.ConstructHDict(infile_1200_180_100.Get("h_"+var+"_"+cut),
                                                name="T1ttcc_1200_180_100", color=rt.kCyan+4,
                                                title="Comparison between different DM lines for new signal scan, for selection "+cut,
                                                appear_in_ratio="Yes", xtitle=var)

            hdictlist=[hdict_10,hdict_25,hdict_80]
            canvasname = "Signal_comparison_DMlines_1200_1xx_100__"+var+"__"+cut
            rtitle = "#frac{DM}{DM10}"
            # plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="No",legdict=leg)
            plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,legdict=leg,cname=canvasname,ratiotitle=rtitle,scale="No")


    outfile.Close()
    infile_1200_110_100.Close()
    infile_1200_125_100.Close()
    infile_1200_180_100.Close()


if __name__ == '__main__':

    # set root styles
    plotTools.SetBoostStyle()

    runA()
    runB()
