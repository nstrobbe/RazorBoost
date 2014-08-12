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

    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass"
    #basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass/summary/"
    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140610_FullStatusReport"
    #basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140610_FullStatusReport/summary/"
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140730_preApp_comments"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140730_preApp_comments/summary/"
    inputfile_TTJ = basedir + "rzrBoostMC_TTJets.root"
    inputfile_Top = basedir + "rzrBoostMC_Top.root"
    inputfile_QCD = basedir + "rzrBoostMC_QCD.root"
    inputfile_WJets = basedir + "rzrBoostMC_WJetsToLNu.root"
    inputfile_DYJets = basedir + "rzrBoostMC_DYJetsToLL.root"
    inputfile_ZJetsToNuNu = basedir + "rzrBoostMC_ZJetsToNuNu.root"
    inputfile_data = basedir + "rzrBoostMC_data.root"
    inputfile_bg = basedir + "../rzrBoostMC_bg.root"
    inputfile_1000_325_300 = basedir + "rzrBoostMC_T1ttcc_1000_325_300.root"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots2.root","RECREATE")
    infile_TTJ = TFile.Open(inputfile_TTJ)
    infile_Top = TFile.Open(inputfile_Top)
    infile_QCD = TFile.Open(inputfile_QCD)
    infile_WJets = TFile.Open(inputfile_WJets)
    infile_DYJets = TFile.Open(inputfile_DYJets)
    infile_Zinv = TFile.Open(inputfile_ZJetsToNuNu)
    infile_data = TFile.Open(inputfile_data)
    infile_bg = TFile.Open(inputfile_bg)
    infile_1000_325_300 = TFile.Open(inputfile_1000_325_300)

    # Integrated luminosity in fb-1s
    intlumi = 19.712 # ABCD

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

    leg = plotTools.ConstructLDict(0.7,0.85,0.65,0.87)

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

    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll_mdPhiHatg4"),
                                             name=" S region", color=rt.kBlack,
                                             title="Shape comparison for TTJets in S and T region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmt = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhiHatg4"),
                                             name=" T region", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets in Sand T region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_TTJmt]
        canvasname = var+"_comparison_TTJ_TvsS"
        rtitle = "#frac{T}{S}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll_mdPhig0p5"),
                                             name=" S region", color=rt.kBlack,
                                             title="Shape comparison for TTJets in S and T region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_TTJmt = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhig0p5"),
                                             name=" T region", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets in Sand T region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_TTJmt]
        canvasname = var+"_comparison_TTJ_TvsS_mdphi"
        rtitle = "#frac{T}{S}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    # compare Wjets and TTjets in signal region and TTJets control region
    for var in vars:
        hdict_TTJ_S = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll"),
                                               name="TTJets in SIG region", color=rt.kBlack,
                                               title="Shape comparison for TTJets and WJets in the Signal and TTJets Control region",
                                               appear_in_ratio="Ref", xtitle=var)
        hdict_TTJ_T = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                               name="TTJets in TTJets CR", color=rt.kCyan+2,
                                               title="Shape comparison for TTJets and WJets in the Signal and TTJets Control region",
                                               appear_in_ratio="Yes", xtitle=var)
        hdict_WJ_S = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W0Ll"),
                                              name="WJets in SIG region", color=rt.kGreen+2,
                                              title="Shape comparison for TTJets and WJets in the Signal and TTJets Control region",
                                              appear_in_ratio="Yes", xtitle=var)
        hdict_WJ_T = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W1LlmT100"),
                                              name="WJets in TTJets CR", color=rt.kSpring+9,
                                              title="Shape comparison for TTJets and WJets in the Signal and TTJets Control region",
                                              appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_TTJ_S, hdict_TTJ_T,hdict_WJ_S,hdict_WJ_T]
        canvasname = var+"_comparison_TTJ_WJ"
        rtitle = "#frac{other}{TTJets SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # compare TOP and TTjets in signal
    for var in vars:
        hdict_TTJ_S = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll_mdPhig0p5"),
                                               name="TTJets in S region", color=rt.kMagenta+3,
                                               title="Shape comparison for TTJets and Single Top in the Signal region",
                                               appear_in_ratio="Ref", xtitle=var, drawoption="EP")
        hdict_Top_S = plotTools.ConstructHDict(infile_Top.Get("h_"+var+"_g1Mbg1W0Ll_mdPhig0p5"),
                                              name="Single top in S region", color=rt.kMagenta+1,
                                              title="Shape comparison for TTJets and Single Top in the Signal region",
                                              appear_in_ratio="Yes", xtitle=var, drawoption="EP")
        hdictlist=[hdict_TTJ_S,hdict_Top_S]
        canvasname = var+"_comparison_TTJ_Top_S"
        rtitle = "#frac{Single top}{TTJ}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # compare TOP and TTjets in TTJets control region
    for var in vars:
        hdict_TTJ_T = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhig0p5"),
                                               name="TTJets in T CR", color=rt.kMagenta+3,
                                               title="Shape comparison for TTJets and Single Top in the TTJets Control region",
                                               appear_in_ratio="Ref", xtitle=var,drawoption="EP")
        hdict_Top_T = plotTools.ConstructHDict(infile_Top.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhig0p5"),
                                              name="Single top in T CR", color=rt.kMagenta+1,
                                              title="Shape comparison for TTJets and Single Top in the TTJets Control region",
                                              appear_in_ratio="Yes", xtitle=var,drawoption="EP")
        hdictlist=[hdict_TTJ_T,hdict_Top_T]
        canvasname = var+"_comparison_TTJ_Top_T"
        rtitle = "#frac{Single Top}{TTJ}"
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

    # build hdictlist for QCD:
    for var in vars:
        hdict_QCD_mdphihat4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhiHat4"),
                                             name="QCD CR", color=rt.kMagenta,
                                             title="Shape comparison for QCD MC in the QCD CR and in the SIG region with mdPhiHat > or < 4",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIG_mdphihat4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll_mdPhiHat4"),
                                             name="SIGNAL region, mdPhiHat < 4", color=rt.kOrange-3,
                                             title="Shape comparison for QCD MC in the QCD CR and in the SIG region with mdPhiHat > or < 4",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_SIG_mdphihatg4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll_mdPhiHatg4"),
                                             name="SIGNAL region, mdPhiHat > 4", color=rt.kBlue-3,
                                             title="Shape comparison for QCD MC in the QCD CR and in the SIG region with mdPhiHat > or < 4",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_QCD_mdphihat4,hdict_SIG_mdphihat4,hdict_SIG_mdphihatg4]
        canvasname = var+"_comparison_QCD_SIG_mdphihat"
        rtitle = "#frac{SIG}{QCD}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # build hdictlist for QCD:
    for var in vars:
        hdict_QCD_mdphihat4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhiHat4"),
                                             name=" Q region", color=rt.kMagenta,
                                             title="Shape comparison for QCD MC in the Q region and the S region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIG_mdphihatg4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll_mdPhiHatg4"),
                                             name=" S region", color=rt.kBlue-3,
                                             title="Shape comparison for QCD MC in the Q region and in the S region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_QCD_mdphihat4,hdict_SIG_mdphihatg4]
        canvasname = var+"_comparison_QCD_QvsS"
        rtitle = "#frac{S}{Q}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    # build hdictlist for QCD with mdphi:
    for var in vars:
        hdict_QCD_mdphihat4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),
                                             name=" Q region", color=rt.kMagenta,
                                             title="Shape comparison for QCD MC in the Q region and the S region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIG_mdphihatg4 = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll_mdPhig0p5"),
                                             name=" S region", color=rt.kBlue-3,
                                             title="Shape comparison for QCD MC in the Q region and in the S region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_QCD_mdphihat4,hdict_SIG_mdphihatg4]
        canvasname = var+"_comparison_QCD_QvsS_mdphi"
        rtitle = "#frac{S}{Q}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    # build hdictlist for Wjets:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for WJets in the Signal and various WJets Control regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_WJets_mt100 = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_0Lbg1Y1LlmT100"),
                                             name="WJets CR, mT < 100", color=rt.kGreen+1,
                                             title="Shape comparison for WJets in the Signal and various WJets Control regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdict_WJets_mt = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_0Lbg1Y1LlmT"),
                                             name="WJets CR, 30 < mT < 100", color=rt.kGreen+2,
                                             title="Shape comparison for WJets in the Signal and various WJets Control regions",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_WJets_mt100,hdict_WJets_mt]
        canvasname = var+"_comparison_WJets"
        rtitle = "#frac{WJets}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # build hdictlist for Wjets:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W0Ll_mdPhiHatg4"),
                                             name="S region", color=rt.kBlack,
                                             title="Shape comparison for WJets in the S and W region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_WJets_mt = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhiHatg4"),
                                             name="W region", color=rt.kGreen+2,
                                             title="Shape comparison for WJets in the S and W region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_WJets_mt]
        canvasname = var+"_comparison_WJ_WvsS"
        rtitle = "#frac{W}{S}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    # build hdictlist for Wjets:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W0Ll_mdPhig0p5"),
                                             name="S region", color=rt.kBlack,
                                             title="Shape comparison for WJets in the S and W region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_WJets_mt = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhig0p5"),
                                             name="W region", color=rt.kGreen+2,
                                             title="Shape comparison for WJets in the S and W region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_WJets_mt]
        canvasname = var+"_comparison_WJ_WvsS_mdphi"
        rtitle = "#frac{W}{S}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes",legdict=leg)

    leps = ["2el0mu","2mu0el","2l0ol"]
    # build hdictlist for Zjets 0b:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_Zinv.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="ZJetsToNuNu in SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for ZJets in the Signal and various ZJets Control regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdictlist=[hdict_SIG]
        for i,lep in enumerate(leps):
            hdict_ZJets = plotTools.ConstructHDict(infile_DYJets.Get("h_"+var+"_0Lbg1Y"+lep),
                                                   name="ZJets CR, 0b, "+lep, color=rt.kOrange+2*i,
                                                   title="Shape comparison for ZJets in the Signal and various ZJets Control regions",
                                                   appear_in_ratio="Yes", xtitle=var)
            hdictlist.append(hdict_ZJets)
        canvasname = var+"_comparison_ZJets0b"
        rtitle = "#frac{ZJets}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # build hdictlist for Zjets g1b:
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_Zinv.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="ZJetsToNuNu in SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for ZJets in the Signal and various ZJets Control regions",
                                             appear_in_ratio="Ref", xtitle=var)
        hdictlist=[hdict_SIG]
        for i,lep in enumerate(leps):
            hdict_ZJets = plotTools.ConstructHDict(infile_DYJets.Get("h_"+var+"_g1Mbg1Y"+lep),
                                                   name="ZJets CR, >= 1b, "+lep, color=rt.kOrange+2*i,
                                                   title="Shape comparison for ZJets in the Signal and various ZJets Control regions",
                                                   appear_in_ratio="Yes", xtitle=var)
            hdictlist.append(hdict_ZJets)
        canvasname = var+"_comparison_ZJets1b"
        rtitle = "#frac{ZJets}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    # build hdictlist for Wlnu vs Zinv in QCD CR:
    for var in vars:
        hdict_W = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),
                                           name="WJetsToLNu", color=rt.kGreen+2,
                                           title="Shape comparison in the QCD control region",
                                           appear_in_ratio="Ref", xtitle=var)
        hdict_Zinv = plotTools.ConstructHDict(infile_Zinv.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),
                                               name="ZJetsToNuNu", color=rt.kOrange,
                                               title="Shape comparison in the QCD Control region",
                                               appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_Zinv,hdict_W]
        canvasname = var+"_comparison_W_Zinv_QCD"
        rtitle = "#frac{Znunu}{Wlnu}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    ###################################
    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for TTJets in Signal and Validation region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIGlike = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mb0Wg1uW0Ll"),
                                             name="SIGNAL like region", color=rt.kCyan+2,
                                             title="Shape comparison for TTJets in Signal and Validation region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_SIGlike]
        canvasname = var+"_comparison_SIG_SIGlike_TTJ"
        rtitle = "#frac{SIGlike}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for QCD in Signal and Validation region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIGlike = plotTools.ConstructHDict(infile_QCD.Get("h_"+var+"_g1Mb0Wg1uW0Ll"),
                                             name="SIGNAL like region", color=rt.kMagenta+2,
                                             title="Shape comparison for QCD in Signal and Validation region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_SIGlike]
        canvasname = var+"_comparison_SIG_SIGlike_QCD"
        rtitle = "#frac{SIGlike}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_SIG = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mbg1W0Ll"),
                                             name="SIGNAL region", color=rt.kBlack,
                                             title="Shape comparison for WJets in Signal and Validation region",
                                             appear_in_ratio="Ref", xtitle=var)
        hdict_SIGlike = plotTools.ConstructHDict(infile_WJets.Get("h_"+var+"_g1Mb0Wg1uW0Ll"),
                                             name="SIGNAL like region", color=rt.kGreen+2,
                                             title="Shape comparison for WJets in Signal and Validation region",
                                             appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_SIG,hdict_SIGlike]
        canvasname = var+"_comparison_SIG_SIGlike_WJets"
        rtitle = "#frac{SIGlike}{SIG}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")


    #############################################################3######
    # build hdictlist to compare dphimegajets
    dphis = ["dphimegajets","dphijet1jet2"]
    legd2 = plotTools.ConstructLDict(0.2,0.4,0.65,0.85)
    for dphi in dphis:
        hdict_QCD = plotTools.ConstructHDict(infile_QCD.Get("h_"+dphi+"_g1Mbg1W0Ll"),
                                             name="QCD", color=rt.kMagenta,
                                             title="",
                                             appear_in_ratio="Ref", xtitle=dphi)
        hdict_TTJ = plotTools.ConstructHDict(infile_TTJ.Get("h_"+dphi+"_g1Mbg1W0Ll"),
                                             name="TTJets", color=rt.kRed,
                                             title="",
                                             appear_in_ratio="Yes", xtitle=dphi)

        hdictlist=[hdict_QCD,hdict_TTJ]
        canvasname = dphi+"_comparison"
        rtitle = "#frac{TTJ}{QCD}"
        plotTools.Plot1D(hdictlist,outputdir,outfile,legdict=legd2,cname=canvasname,scale="Yes")

    # build hdictlist to compare HT
    legd3 = plotTools.ConstructLDict(0.6,0.87,0.65,0.85)
    hdict_data = plotTools.ConstructHDict(infile_data.Get("h_HT_SIG"),
                                          name="Data", color=rt.kBlack,
                                          title="HT after baseline selection",
                                          xtitle="HT")
    hdict_QCD = plotTools.ConstructHDict(infile_QCD.Get("h_HT_SIG"),
                                         name="QCD", color=rt.kMagenta,
                                         title="HT after baseline selection",
                                         xtitle="HT")
    hdict_WJets = plotTools.ConstructHDict(infile_WJets.Get("h_HT_SIG"),
                                         name="WJets", color=rt.kGreen+2,
                                         title="HT after baseline selection",
                                         xtitle="HT")
    hdict_DYJets = plotTools.ConstructHDict(infile_DYJets.Get("h_HT_SIG"),
                                            name="DYJets", color=rt.kOrange+7,
                                            title="HT after baseline selection",
                                            xtitle="HT")
    hdict_Zinv = plotTools.ConstructHDict(infile_Zinv.Get("h_HT_SIG"),
                                         name="ZJetsToNuNu", color=rt.kOrange,
                                         title="HT after baseline selection",
                                         xtitle="HT")
    
    hdictlist=[hdict_data,hdict_QCD,hdict_WJets,hdict_DYJets,hdict_Zinv]
    canvasname = "HT_comparison"
    plotTools.Plot1D(hdictlist,outputdir,outfile,legdict=legd3,logscale=True,cname=canvasname,scale="Yes")

    #############################################

    for var in vars:
        hdict_T_mdphihatg4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhiHatg4"),
                                                      name="T region, minDeltaPhiHat > 4", color=rt.kRed,
                                                      title="Shape comparison for total background in T region",
                                                      appear_in_ratio="Ref", xtitle=var)
        hdict_T_mdphihat4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhiHat4"),
                                                     name="T region, minDeltaPhiHat < 4", color=rt.kRed+2,
                                                     title="Shape comparison for total background in T region",
                                                     appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_T_mdphihatg4,hdict_T_mdphihat4]
        canvasname = var+"_comparison_T_mdphihat"
        rtitle = "#frac{min#Delta#Phi > 4}{min#Delta#Phi < 4}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_W_mdphihatg4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhiHatg4"),
                                                      name="W region, minDeltaPhiHat > 4", color=rt.kRed,
                                                      title="Shape comparison for total background in W region",
                                                      appear_in_ratio="Ref", xtitle=var)
        hdict_W_mdphihat4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhiHat4"),
                                                     name="W region, minDeltaPhiHat < 4", color=rt.kRed+2,
                                                     title="Shape comparison for total background in W region",
                                                     appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_W_mdphihatg4,hdict_W_mdphihat4]
        canvasname = var+"_comparison_W_mdphihat"
        rtitle = "#frac{min#Delta#Phi > 4}{min#Delta#Phi < 4}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_T_mdphihatg4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhiHatg4"),
                                                      name="T region, minDeltaPhiHat > 4", color=rt.kRed,
                                                      title="Shape comparison for data in T region",
                                                      appear_in_ratio="Ref", xtitle=var)
        hdict_T_mdphihat4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_g1Mbg1W1LlmT100_mdPhiHat4"),
                                                     name="T region, minDeltaPhiHat < 4", color=rt.kRed+2,
                                                     title="Shape comparison for data in T region",
                                                     appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_T_mdphihatg4,hdict_T_mdphihat4]
        canvasname = var+"_comparison_T_mdphihat_data"
        rtitle = "#frac{min#Delta#Phi > 4}{min#Delta#Phi < 4}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    for var in vars:
        hdict_W_mdphihatg4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhiHatg4"),
                                                      name="W region, minDeltaPhiHat > 4", color=rt.kRed,
                                                      title="Shape comparison for data in W region",
                                                      appear_in_ratio="Ref", xtitle=var)
        hdict_W_mdphihat4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_0Lbg1Y1LlmT_mdPhiHat4"),
                                                     name="W region, minDeltaPhiHat < 4", color=rt.kRed+2,
                                                     title="Shape comparison for data in W region",
                                                     appear_in_ratio="Yes", xtitle=var)
        hdictlist=[hdict_W_mdphihatg4,hdict_W_mdphihat4]
        canvasname = var+"_comparison_W_mdphihat_data"
        rtitle = "#frac{min#Delta#Phi > 4}{min#Delta#Phi < 4}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,ratiotitle=rtitle,scale="Yes")

    ########################################################
    # Make comparison of gen level W pt
    hname = "h_gen_Wpt"
    htitle = ""
    xt = "generator W pT"
    legd4 = plotTools.ConstructLDict(0.5,0.87,0.7,0.85)
    hdict_TTJ = plotTools.ConstructHDict(infile_TTJ.Get(hname),
                                         name="TTJets", color=rt.kRed, 
                                         title=htitle,
                                         xtitle=xt)
    hdict_1000_325_300 = plotTools.ConstructHDict(infile_1000_325_300.Get(hname),
                                         name="T1ttcc_1000_325_300", color=rt.kCyan+2, 
                                         title=htitle,
                                         xtitle=xt)
    #hdict_1200_125_100 = plotTools.ConstructHDict(infile_1000_125_100.Get(hname),
    #                                     name="T1ttcc_1000_125_100", color=rt.kTeal, 
    #                                     title=htitle,
    #                                     xtitle=xt)
    hdictlist=[hdict_TTJ,hdict_1000_325_300]
    canvasname = "comparison_genWpt"
    plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="Yes",legdict=legd4)

    ########################################################
    # Make comparison of W pt in signal region
    hname = "h_Wpt_g1Mbg1W0Ll_mdPhig0p5"
    htitle = ""
    xt = "W pT"
    legd4 = plotTools.ConstructLDict(0.5,0.87,0.7,0.85)
    hdict_bg = plotTools.ConstructHDict(infile_bg.Get(hname),
                                         name="Total background (MC)", color=rt.kRed, 
                                         title=htitle,
                                         xtitle=xt,drawoption="histE2")
    hdict_1000_325_300 = plotTools.ConstructHDict(infile_1000_325_300.Get(hname),
                                         name="T1ttcc_1000_325_300", color=rt.kCyan+2, 
                                         title=htitle,
                                         xtitle=xt,drawoption="histE2")
    hdictlist=[hdict_bg,hdict_1000_325_300]
    canvasname = "comparison_Wpt"
    plotTools.Plot1D(hdictlist,outputdir,outfile,cname=canvasname,scale="Yes",legdict=legd4)


    outfile.Close()
    infile_TTJ.Close()
    infile_QCD.Close()
    infile_WJets.Close()
    infile_DYJets.Close()
    infile_Zinv.Close()
    infile_data.Close()
    infile_bg.Close()
    infile_1000_325_300.Close()
