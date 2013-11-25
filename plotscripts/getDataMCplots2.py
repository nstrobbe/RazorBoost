import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20131125_varbin"
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20131122/summary/"
    analyzer ="rzrBoostMC"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/plots.root","RECREATE")

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    plotTools.SetBoostStyle()
    
    print "Will make plots for integrated luminosity of %.3f fb-1" % (intlumi)

    # define all the datasets we want to plot, and their colors
    # backgrounds
    mc_datasets = ["QCD","TTJets","WJetsToLNu","Wbb","Top","TTX","ZJetsToNuNu","DYJetsToLL","DYToBB","DYToCC","VV","VVV"]
    mc_colors   = [rt.kMagenta,rt.kRed,rt.kGreen+1,rt.kGreen+3,rt.kCyan,rt.kCyan+2,rt.kOrange,rt.kOrange+2,rt.kOrange+7,rt.kOrange+9,rt.kBlue+1,rt.kBlue-3]
    flist = []
    for d in mc_datasets:
        f = TFile.Open(inputdir+analyzer+"_"+d+".root")
        flist.append(f)
    # signal
    sig_datasets = ["T1ttcc_325_300"]
    sig_colors   = [rt.kGray]
    fsiglist = []
    for d in sig_datasets:
        f = TFile.Open(inputdir+analyzer+"_"+d+".root")
        fsiglist.append(f)
    # data
    fdata = TFile.Open(inputdir+analyzer+"_data.root")
        
    # make the dictionaries to pass to the plot routine
    vars = ["MR","R2"]
    cuts = ["Cleaning","HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
            "SIG","neleeq0","nmueq0","trackIso",
            "g1Mb0Ll","g1Mbg1W0Ll","0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100","g1Mbg1W1LlmT",
            "2munoZmass","2mu","2mu0el","0Lb2mu0el","0Lbg1Y2mu0el","g1Mb2mu0el","g1Mbg1Y2mu0el"
            ]

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = "Data/MC comparison plot"
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle)
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle)
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events",markerstyle=20)

            # now make the actual plot
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="No")
            if var == "MR":
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events/(100 GeV)",markerstyle=20)
                plotTools.PlotDataMC(hlist,hdict_data2,hsiglist,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo="Selection %s"%(cut),
                                     ratiotitle="Data/MC", logscale=True, scale="Width", scalefactor=100)
            else:
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events/(0.01)",markerstyle=20)
                
                plotTools.PlotDataMC(hlist,hdict_data2,hsiglist,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo="Selection %s"%(cut),
                                     ratiotitle="Data/MC", logscale=True, scale="Width", scalefactor=0.01)

            htitle = "Data/MC comparison plot with total MC integral scaled to match Data"
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s_scaled"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="Yes")

            htitle = "Data/MC comparison plot with QCD scaled up to match Data in first bin"
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s_QCDscaled"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="QCD")

#             htitle = "Data/MC comparison plot with QCD scaled up by factor of 1.98"
#             plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                  cname="DataMC_%s_%s_QCDscaled_factor1p98"%(var,cut), plotinfo="Selection %s"%(cut),
#                                  ratiotitle="Data/MC", logscale=True, scale="QCD", scalefactor=1.98)

#             htitle = "Data/MC comparison plot with QCD scaled up by factor of 2.09"
#             plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                  cname="DataMC_%s_%s_QCDscaled_factor2p09"%(var,cut), plotinfo="Selection %s"%(cut),
#                                  ratiotitle="Data/MC", logscale=True, scale="QCD", scalefactor=2.09)

    # make plot of minDeltaPhiHat
    hname = "h_minDeltaPhiHat_0Lbg1uW0Ll"
    htitle = "Data/MC comparison plot"
    hlist = []
    for i in range(len(mc_datasets)):
        if not flist[i]: continue
        hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle)
        hlist.append(hdict)
        
    hsiglist = []
    for i in range(len(sig_datasets)):
        if not fsiglist[i]: continue
        hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle)
        hsiglist.append(hdict)

    hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle="MR",ytitle="Events",markerstyle=20)

    plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
                         cname="DataMC_minDeltaPhiHat_0Lbg1uW0Ll", plotinfo="Selection 0Lbg1uW0Ll",
                         ratiotitle="Data/MC", logscale=True, scale="No")
    

    outfile.Close()

