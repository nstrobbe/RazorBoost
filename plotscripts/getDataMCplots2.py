import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass"
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass/summary/"
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
    #mc_datasets = ["QCD","TTJets","WJetsToLNu","Top","TTX","ZJetsToNuNu","DYJetsToLL_PtZ","VV","VVV"]
    #mc_colors   = [rt.kMagenta,rt.kRed,rt.kGreen+1,rt.kCyan,rt.kCyan+2,rt.kOrange,rt.kOrange+2,rt.kBlue+1,rt.kBlue-3]
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
            "g1Mb0Ll","g1Mbg1W0Ll","1Mbg1W0Ll","g2Mbg1W0Ll","g1Mbg1W0Ll_mdPhiHatg4","g1Mbg1W0Ll_mdPhiHat4","g1Mb0Wg1uW0Ll",
            "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhi0p5","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100",
            "1Mbg1W1LlmT100","g2Mbg1W1LlmT100","g1Mbg1W1LlmT","g1Mbg1W1LlmT100_mdPhiHatg4","g1Mbg1W1LlmT100_mdPhiHat4",
            "0Lb1Ll","0Lbg1Y1Ll","0Lbg1Y1LlmT100","0Lbg1Y1LlmT",
            "0Lbg1Y1LlmT_mdPhiHatg4","0Lbg1Y1LlmT_mdPhiHat4",
            "2munoZmass","2mu","2mu0el","0Lb2mu0el","0Lbg1Y2mu0el","g1Mb2mu0el","g1Mbg1Y2mu0el",
            "2elnoZmass","2el","2el0mu","0Lb2el0mu","0Lbg1Y2el0mu","g1Mb2el0mu","g1Mbg1Y2el0mu",
            "2lnoZmass","2l","2l0ol","0Lb2l0ol","0Lbg1Y2l0ol","g1Mb2l0ol","g1Mbg1Y2l0ol",
            "g1Mbg1W0Ll_mdPhig0p3","g1Mbg1W0Ll_mdPhi0p3","g1Mbg1W1LlmT100_mdPhig0p3","g1Mbg1W1LlmT100_mdPhi0p3",
            "0Lbg1Y1LlmT_mdPhig0p3","0Lbg1Y1LlmT_mdPhi0p3",
            "g1Mbg1W0Ll_mdPhig0p5","g1Mbg1W0Ll_mdPhi0p5","g1Mbg1W1LlmT100_mdPhig0p5","g1Mbg1W1LlmT100_mdPhi0p5",
            "0Lbg1Y1LlmT_mdPhig0p5","0Lbg1Y1LlmT_mdPhi0p5",
            ]

    legd = plotTools.ConstructLDict(0.6,0.87,0.5,0.8,ncolumns=2)

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = "Data/MC comparison plot"
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=var)
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=var)
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events",markerstyle=20)
            #hdict_data = 0
            
            # now make the actual plot
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="No")

            # scale according to bin width; need to adjust y axis title according to variable
            if var == "MR":
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events/(100 GeV)",markerstyle=20)
                #hdict_data2=0
                plotTools.PlotDataMC(hlist,hdict_data2,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo="Selection %s"%(cut),
                                     ratiotitle="Data/MC", logscale=True, scale="Width", scalefactor=100)
            else:
                hdict_data2 = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events/(0.01)",markerstyle=20)
                #hdict_data2=0
                plotTools.PlotDataMC(hlist,hdict_data2,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                     cname="DataMC_%s_%s_width"%(var,cut), plotinfo="Selection %s"%(cut),
                                     ratiotitle="Data/MC", logscale=True, scale="Width", scalefactor=0.01)

#            htitle = "Data/MC comparison plot with total MC integral scaled to match Data"
#            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                 cname="DataMC_%s_%s_scaled"%(var,cut), plotinfo="Selection %s"%(cut),
#                                 ratiotitle="Data/MC", logscale=True, scale="Yes")

#            htitle = "Data/MC comparison plot with QCD scaled up to match Data in first bin"
#            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                 cname="DataMC_%s_%s_QCDscaled"%(var,cut), plotinfo="Selection %s"%(cut),
#                                 ratiotitle="Data/MC", logscale=True, scale="QCD")

#             htitle = "Data/MC comparison plot with QCD scaled up by factor of 1.98"
#             plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                  cname="DataMC_%s_%s_QCDscaled_factor1p98"%(var,cut), plotinfo="Selection %s"%(cut),
#                                  ratiotitle="Data/MC", logscale=True, scale="QCD", scalefactor=1.98)

#             htitle = "Data/MC comparison plot with QCD scaled up by factor of 2.09"
#             plotTools.PlotDataMC(hlist,hdict_data,hsiglist,outputdir=outputdir, outfile=outfile,
#                                  cname="DataMC_%s_%s_QCDscaled_factor2p09"%(var,cut), plotinfo="Selection %s"%(cut),
#                                  ratiotitle="Data/MC", logscale=True, scale="QCD", scalefactor=2.09)

    # make plot of other variables
    vars = ["minDeltaPhiHat","minDeltaPhi"]
    cuts = ["HLT","jet1ptg200","SIG","0Ll","0Lb0Ll","0Lbg1uW0Ll"]
    htitle = "Data/MC comparison plot"
    for var in vars:
        for cut in cuts:
            hname = "h_"+var+"_"+cut
            print "checking stuff: ", hname
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

            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_"+hname.replace("h_",""), plotinfo="Selection "+cut,
                                 ratiotitle="Data/MC", logscale=True, scale="No")

    hnames = ["h_mT_g1Mbg1W1Ll","h_mT_0Lbg1Y1Ll","h_mT_0Lb1Ll"]
    htitle = "Data/MC comparison plot"
    for hname in hnames:
        var = hname.split("_")[1]
        
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

        plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                             cname="DataMC_"+hname.replace("h_",""), plotinfo="Selection "+hname.replace("h_mT_",""),
                             ratiotitle="Data/MC", logscale=True, scale="No")

    #hnames = ["h_TrueNumVertices","h_TrueNumVertices_reweighted"]
    hnames = ["h_PV","h_PV_reweighted"
              , "h_PV_nosel", "h_PV_reweighted_nosel"
              , "h_PV_1j", "h_PV_reweighted_1j"
              , "h_PV_HLT", "h_PV_reweighted_HLT"
              , "h_PV_SIG", "h_PV_reweighted_SIG" 
              , "h_PV_g1Mbg1W1LlmT100", "h_PV_reweighted_g1Mbg1W1LlmT100"
              ]
    htitle = "Data/MC comparison plot"
    sels = ["","",
            "No Selection", "No Selection",
            "njets >= 1", "njets >= 1",
            "njets >= 3", "njets >= 3",
            "SIG", "SIG",
            "g1Mbg1W1LlmT100", "g1Mbg1W1LlmT100"]
    for hname in hnames:
        var = hname.split("_")[1]
        
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

        plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                             cname="DataMC_"+hname.replace("h_",""), plotinfo=sels[hnames.index(hname)],
                             ratiotitle="Data/MC", logscale=True, scale="No")


    vars = ["njets","nbjets","met","jet1pt","jet2pt","jet3pt",
            "leptonpt","lepton1pt","lepton2pt",
            "HT","PV"]

    cuts = ["SIG","g1Mbg1W0Ll","g1Mbg1W0Ll_mdPhiHatg4","g1Mb0Wg1uW0Ll",
            "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhi0p5","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5",
            "g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100",
            "0Lb1Ll","0Lbg1Y1Ll","0Lbg1Y1LlmT","0Lbg1Y1LlmT100",
            "2mu","2mu0el","0Lb2mu0el","g1Mb2mu0el","0Lbg1Y2mu0el","g1Mbg1Y2mu0el",
            "2el","2el0mu","0Lb2el0mu","g1Mb2el0mu","0Lbg1Y2el0mu","g1Mbg1Y2el0mu",
            "2l","2l0ol","0Lb2l0ol","g1Mb2l0ol","0Lbg1Y2l0ol","g1Mbg1Y2l0ol"
            "g1Mbg1W0Ll_mdPhig0p3","g1Mbg1W1LlmT100_mdPhig0p3","0Lbg1Y1LlmT_mdPhig0p3",
            "g1Mbg1W0Ll_mdPhig0p5","g1Mbg1W1LlmT100_mdPhig0p5","0Lbg1Y1LlmT_mdPhig0p5"
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
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="No")


    vars = ["dphimegajets","dphijet1jet2"]
    cuts = ["g1Mbg1W0Ll"]
    legd2 = plotTools.ConstructLDict(0.15,0.37,0.55,0.85,ncolumns=2)
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
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd2,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="No")


    vars = ["minDeltaPhiHat","minDeltaPhi"]
    cuts = ["g1Mbg1W0Ll","g1Mbg1W1LlmT100","0Lbg1Y1LlmT"]
    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = "Data/MC comparison plot"
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=var,ytitle="Events")
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=var,ytitle="Events")
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events",markerstyle=20)

            # now make the actual plot
            plotTools.PlotDataMC(hlist,hdict_data,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="No")

    for cut in cuts:
        for var in vars:
            hname = "h_%s_%s" % (var,cut)
            htitle = "MC composition plot"
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=var,ytitle="Events")
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=var,ytitle="Events")
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events",markerstyle=20)
            plotTools.PlotDataMC(hlist,0,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s_nodata"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=False, scale="No")


    # Make plot of signal region without data
    vars = ["MR","R2"]
    cuts = ["g1Mbg1W0Ll_mdPhiHatg4","g1Mbg1W0Ll_mdPhig0p3"]
    sf = [100,0.01]
    for cut in cuts:
        for vi,var in enumerate(vars):
            hname = "h_%s_%s" % (var,cut)
            htitle = "MC composition plot"
            hlist = []
            for i in range(len(mc_datasets)):
                if not flist[i]: continue
                hdict = plotTools.ConstructHDict(flist[i].Get(hname),name=mc_datasets[i],color=mc_colors[i],title=htitle,xtitle=var,ytitle="Events/(%s GeV)"%(sf[vi]))
                hlist.append(hdict)
        
            hsiglist = []
            for i in range(len(sig_datasets)):
                if not fsiglist[i]: continue
                hdict = plotTools.ConstructHDict(fsiglist[i].Get(hname),name=sig_datasets[i],color=sig_colors[i],title=htitle,xtitle=var,ytitle="Events/(%s GeV)"%(sf[vi]))
                hsiglist.append(hdict)

            hdict_data = plotTools.ConstructHDict(fdata.Get(hname),name="data",color=rt.kBlack,title=htitle,xtitle=var,ytitle="Events",markerstyle=20)
            plotTools.PlotDataMC(hlist,0,hsiglist,legdict=legd,outputdir=outputdir, outfile=outfile,
                                 cname="DataMC_%s_%s_width_nodata"%(var,cut), plotinfo="Selection %s"%(cut),
                                 ratiotitle="Data/MC", logscale=True, scale="Width", scalefactor=sf[vi])


    outfile.Close()

