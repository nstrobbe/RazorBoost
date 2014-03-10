import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools


def makeplot(sigregion, cregion, basedir, estdir, outputdir, MConly=True):
    inputfile_data = basedir + "/summary/rzrBoostMC_data.root" # data histograms
    inputfile_totalbg = basedir + "/rzrBoostMC_bg.root" # data histograms
    inputfile_estimate = estdir + "/BGestimate_simple_"+sigregion+"_from_"+cregion+".root" # bg estimate histograms

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/BGplotsSimple_"+sigregion+"_from_"+cregion+".root","RECREATE")
    infile_bg = TFile.Open(inputfile_totalbg)
    infile_data = TFile.Open(inputfile_data)
    infile_estimate = TFile.Open(inputfile_estimate)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    
    vars = ["MR","R2"]

    # get the 1D projections for the estimate
    hname = "h_MR_R2_"+cregion
    h_2D = infile_estimate.Get(hname)
    h_1D = {}
    h_1D["MR"] = h_2D.ProjectionX("h_1D_MR")
    h_1D["R2"] = h_2D.ProjectionY("h_1D_R2")

    # build legend dictionary
    leg = plotTools.ConstructLDict(0.5,0.87,0.6,0.87)
    
    # build hdictlists
    for var in vars:
        htitle = "Comparison Data vs BG estimate for region "+sigregion+" from region "+cregion
        hdict_data = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+sigregion),
                                              name="Data", color=rt.kBlack,
                                              title=htitle, drawoption="E1X0 P", markerstyle=20,
                                              appear_in_ratio="Ref", xtitle=var)
        hdict_estimate = plotTools.ConstructHDict(h_1D[var],
                                                  name="Data-driven BG estimate", color=rt.kCyan+2,
                                                  title=htitle, drawoption="E2", fillstyle=3002,
                                                  appear_in_ratio="Yes", xtitle=var)
        hdict_bg = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+sigregion),
                                            name="Full MC BG estimate", color=rt.kRed+2,
                                            title=htitle, drawoption="E0", fillstyle=3002,
                                            appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_estimate,hdict_data]
        if MConly:
            hdictlist=[hdict_estimate,hdict_bg,hdict_data]
        canvasname = var+"_comparison_data_estimate_"+sigregion+"_from_"+cregion
        rtitle = "#frac{BG estimate}{Data}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,
                                  ratiotitle=rtitle,scale=False,legdict=leg,
                                  cdim=[696,550],ratiodim=0.3)
        
    
    outfile.Close()


if __name__ == '__main__':

    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140306_noISR_btag_TopPt/"
    estdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140306_noISR_btag_TopPt/"
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140306_noISR_btag_TopPt/"

    # set root styles
    plotTools.SetBoostStyle()

    makeplot("g1Mbg1W1LlmT100", "0Lbg1Y1LlmT", basedir, estdir, outputdir)
    makeplot("1Mbg1W1LlmT100", "0Lbg1Y1LlmT", basedir, estdir, outputdir)
    makeplot("g2Mbg1W1LlmT100", "1Mbg1W1LlmT100", basedir, estdir, outputdir)
    makeplot("0Lbg1uW0Ll_mdPhiHatg4","0Lbg1uW0Ll_mdPhiHat4", basedir, estdir, outputdir)
    makeplot("g1Mbg1W1LlmT100_mdPhiHatg4","g1Mbg1W1LlmT100_mdPhiHat4", basedir, estdir, outputdir)
    makeplot("0Lbg1Y1LlmT_mdPhiHatg4","0Lbg1Y1LlmT_mdPhiHat4", basedir, estdir, outputdir)
