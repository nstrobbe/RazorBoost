import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140917/forPAS/"
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140917/summary/"
    analyzer = "rzrBoostMC"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/plots.root","RECREATE")

    # Integrated luminosity in fb-1s
    intlumi = 19.7 # ABCD

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
    sig_datasets = ["T1ttcc_1000_325_300"]
    sig_colors   = [rt.kGray]
    fsiglist = []
    for d in sig_datasets:
        f = TFile.Open(inputdir+analyzer+"_"+d+".root")
        fsiglist.append(f)
    # data
    fdata = TFile.Open(inputdir+analyzer+"_data.root")
        
    # make the dictionaries to pass to the plot routine
    vars = ["MR","R2"]
    vartitles = ["M_{R} (GeV)","R^{2}"]
    cuts = ["g1Mbg1W0Ll_mdPhig0p5",
            "g1Mbg1W1LlmT100_mdPhig0p5",
            "0Lbg1Y1LlmT_mdPhig0p5",
            "0Lbg1uW0Ll_mdPhi0p3"
            ]
    plotinfos = ["Signal region","T region", "W region", "Q region"]
    legd = plotTools.ConstructLDict(0.65,0.92,0.5,0.8,ncolumns=2)

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = ""
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=vartitles[vars.index(var)],ytitle="Events",markerstyle=20)
            
            # now make the actual plot
            plotTools.PlotDataMCPAS(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                 ratiotitle="Data/MC ", logscale=True, scale="No", style="CMS")

            # scale according to bin width; need to adjust y axis title according to variable
            if var == "MR":
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,
                                                       xtitle=vartitles[vars.index(var)],ytitle="Events/(100 GeV)",markerstyle=20)
                #hdict_data2=0
                plotTools.PlotDataMCPAS(hlist,hdict_data2,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                     ratiotitle="Data/MC ", logscale=True, scale="Width", scalefactor=100, style="CMS")
            else:
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,
                                                       xtitle=vartitles[vars.index(var)],ytitle="Events/(0.01)",markerstyle=20)
                #hdict_data2=0
                plotTools.PlotDataMCPAS(hlist,hdict_data2,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                     ratiotitle="Data/MC ", logscale=True, scale="Width", scalefactor=0.01, style="CMS")

    ##############################################################################################################

    # Make plot of signal region without data
    vars = ["MR","R2"]
    vartitles = ["M_{R} (GeV)","R^{2}"]
    cuts = ["g1Mbg1W0Ll_mdPhig0p5"]
    plotinfos = ["Signal region"]
    sf = [100,0.01]
    for cut in cuts:
        for vi,var in enumerate(vars):
            hname = "h_%s_%s" % (var,cut)
            htitle = ""
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,
                                                 xtitle=vartitles[vars.index(var)],ytitle="Events/(%s GeV)"%(sf[vi]))
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,
                                                 xtitle=vartitles[vars.index(var)],ytitle="Events/(%s GeV)"%(sf[vi]))
                hsiglist.append(hdict)

            plotTools.PlotDataMCPAS(hlist,0,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s_width_nodata"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                 ratiotitle="Data/MC ", logscale=True, scale="Width", scalefactor=sf[vi], style="CMS")


    ##############################################################################################################

    # make mT plots
    vars = ["mT"]
    vartitles = ["m_{T} (GeV)"]
    cuts = ["g1Mbg1W1Ll_rebin",
            "0Lbg1Y1Ll_rebin",
            ]
    plotinfos = ["T region, no selection on m_{T} and #Delta#phi_{min}", "W region, no selection on m_{T} and #Delta#phi_{min}"]
    legd = plotTools.ConstructLDict(0.65,0.92,0.5,0.8,ncolumns=2)

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = ""
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=vartitles[vars.index(var)],ytitle="Events",markerstyle=20)
            
            # now make the actual plot
            plotTools.PlotDataMCPAS(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                 ratiotitle="Data/MC ", logscale=True, scale="No", style="CMS")

    ##############################################################################################################

    # make mindeltaphi plots
    vars = ["minDeltaPhi"]
    vartitles = ["#Delta#phi_{min}"]
    cuts = ["0Lbg1uW0Ll_rebin",
            ]
    plotinfos = ["Q region, no selection on #Delta#phi_{min}"]
    legd = plotTools.ConstructLDict(0.65,0.92,0.5,0.8,ncolumns=2)

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = ""
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=vartitles[vars.index(var)])
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=vartitles[vars.index(var)],ytitle="Events",markerstyle=20)
            
            # now make the actual plot
            plotTools.PlotDataMCPAS(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo=plotinfos[cuts.index(cut)],
                                 ratiotitle="Data/MC ", logscale=True, scale="No", style="CMS")

    ##############################################################################################################


    outfile.Close()

