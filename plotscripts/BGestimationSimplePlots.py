import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/estimate_20140204/"

    sigregion = "g2Mbg1W1LlmT100"
    cregion = "1Mbg1W1LlmT100"

    inputfile_data = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140204/summary/rzrBoostMC_data.root" # data histograms
    inputfile_totalbg = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140204/rzrBoostMC_bg.root" # data histograms
    inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/estimate_20140204/BGestimate_simple_"+sigregion+"_from_"+cregion+".root" # bg estimate histograms

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/BGplotsSimple_"+sigregion+"_from_"+cregion+".root","RECREATE")
    infile_bg = TFile.Open(inputfile_totalbg)
    infile_data = TFile.Open(inputfile_data)
    infile_estimate = TFile.Open(inputfile_estimate)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
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
        #hdictlist=[hdict_estimate,hdict_bg,hdict_data]
        hdictlist=[hdict_estimate,hdict_data]
        canvasname = var+"_comparison_data_estimate_"+sigregion+"_from_"+cregion
        rtitle = "#frac{BG estimate}{Data}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,
                                  ratiotitle=rtitle,scale=False,legdict=leg)
        
    
    outfile.Close()
