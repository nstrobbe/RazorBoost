import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

import plotTools

if __name__ == '__main__':

    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140522_newWtagger/"
    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140528_newWtagger_global/"
    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140606_TopinTTJets/"
    #outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140620_tests/"
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140814/"

    #region = "g1Mb0Wg1uW0Ll"
    #region = "g1Mbg1W0Ll_mdPhiHat4"
    #region = "g1Mbg1W0Ll_mdPhiHatg4"
    #region = "g1Mbg1W0Ll_mdPhi0p3"
    #region = "g1Mbg1W0Ll_mdPhig0p3"
    region = "g1Mbg1W0Ll_mdPhi0p5"
    #region = "g1Mbg1W0Ll_mdPhig0p5"
    
    postfix = ""

    #inputfile_data = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass/summary/rzrBoostMC_data.root" # data histograms
    #inputfile_totalbg = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass/rzrBoostMC_bg.root" # data histograms
    #inputfile_data = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140610_FullStatusReport/summary/rzrBoostMC_data.root" # data histograms
    #inputfile_totalbg = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140610_FullStatusReport/rzrBoostMC_bg.root" # data histograms
    inputfile_data = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140814/summary/rzrBoostMC_data.root" # data histograms
    inputfile_totalbg = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140814/rzrBoostMC_bg.root" # data histograms
    
    #inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140522_newWtagger/BGestimate_"+region+"_QCDWJTTJ_2.root" # bg estimate histograms
    #inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140528_newWtagger_global/BGestimate_"+region+"_QCDWJTTJ_2.root" # bg estimate histograms
    #inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140610_FullStatusReport/BGestimate_"+region+"_QCDWJTTJ.root" # bg estimate histograms
    inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140814/BGestimate_"+region+"_QCDWJTTJ.root" # bg estimate histograms
    #inputfile_estimate = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/closuretest_20140619_test_QCDbinbybin/BGestimate_"+region+"_QCDWJTTJ.root" # bg estimate histograms

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/BGplots_"+region+postfix+".root","RECREATE")
    infile_bg = TFile.Open(inputfile_totalbg)
    infile_data = TFile.Open(inputfile_data)
    infile_estimate = TFile.Open(inputfile_estimate)

    # Integrated luminosity in fb-1s
    intlumi = 19.712 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR","R2"]

    # Get the separate components
    hname_WJ = "h_WJetsToLNu_in_"+region+"_from_0Lbg1Y1LlmT_mdPhig0p5"
    h_2D_WJ = infile_estimate.Get(hname_WJ)
    #h_2D_WJ.Scale(1.4)
    h_1D_WJ = {}
    h_1D_WJ["MR"] = h_2D_WJ.ProjectionX("h_1D_MR_WJ")
    h_1D_WJ["R2"] = h_2D_WJ.ProjectionY("h_1D_R2_WJ")
    hname_TTJ = "h_TTJets_in_"+region+"_from_g1Mbg1W1LlmT100_mdPhig0p5"
    h_2D_TTJ = infile_estimate.Get(hname_TTJ)
    h_1D_TTJ = {}
    h_1D_TTJ["MR"] = h_2D_TTJ.ProjectionX("h_1D_MR_TTJ")
    h_1D_TTJ["R2"] = h_2D_TTJ.ProjectionY("h_1D_R2_TTJ")
    hname_QCD = "h_QCD_in_"+region+"_from_0Lbg1uW0Ll_mdPhi0p3"
    h_2D_QCD = infile_estimate.Get(hname_QCD)
    #h_2D_QCD_toadd = h_2D_QCD.Clone("to_add")
    #h_2D_QCD_toadd.Scale(0.4)
    #h_2D_QCD.Scale(1.4)
    h_1D_QCD = {}
    h_1D_QCD["MR"] = h_2D_QCD.ProjectionX("h_1D_MR_QCD")
    h_1D_QCD["R2"] = h_2D_QCD.ProjectionY("h_1D_R2_QCD")

    # instead of using a stack (not supported yet in plotTools), make sum of histograms
    h_1D_TTJ_WJ = {}
    h_1D_TTJ_WJ["MR"] = h_1D_TTJ["MR"].Clone("h_1D_MR_TTJ_WJ")
    h_1D_TTJ_WJ["MR"].Add(h_1D_WJ["MR"])
    h_1D_TTJ_WJ["R2"] = h_1D_TTJ["R2"].Clone("h_1D_R2_TTJ_WJ")
    h_1D_TTJ_WJ["R2"].Add(h_1D_WJ["R2"])

    h_1D_QCD_TTJ_WJ = {}
    print h_1D_TTJ_WJ["MR"]
    h_1D_QCD_TTJ_WJ["MR"] = h_1D_TTJ_WJ["MR"].Clone("h_1D_MR_QCD_TTJ_WJ")
    h_1D_QCD_TTJ_WJ["MR"].Add(h_1D_QCD["MR"])
    h_1D_QCD_TTJ_WJ["R2"] = h_1D_TTJ_WJ["R2"].Clone("h_1D_R2_QCD_TTJ_WJ")
    h_1D_QCD_TTJ_WJ["R2"].Add(h_1D_QCD["R2"])

    # get the 1D projections for the total estimate
    hname = "h_MR_R2_"+region
    h_2D = infile_estimate.Get(hname)
    #h_2D.Add(h_2D_QCD_toadd)
    h_1D = {}
    h_1D["MR"] = h_2D.ProjectionX("h_1D_MR")
    h_1D["R2"] = h_2D.ProjectionY("h_1D_R2")

    # build legend dictionary
    leg = plotTools.ConstructLDict(0.5,0.87,0.6,0.87)
    
    # build hdictlists
    for var in vars:
        htitle = "Comparison Data vs BG estimate for region "+region
        hdict_data = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+region),
                                              name="Data", color=rt.kBlack,
                                              title=htitle, drawoption="E1X0 P", markerstyle=20,
                                              appear_in_ratio="Ref", xtitle=var)
        hdict_estimate = plotTools.ConstructHDict(h_1D[var],
                                                  name="Data-driven BG estimate", color=rt.kCyan+2,
                                                  title=htitle, drawoption="E2", fillstyle=3002,
                                                  appear_in_ratio="Yes", xtitle=var)
        hdict_estimate_WJ = plotTools.ConstructHDict(h_1D_WJ[var],
                                                  name="WJets estimate", color=rt.kGreen+2,
                                                  title=htitle, drawoption="hist", 
                                                  appear_in_ratio="No", xtitle=var)
        hdict_estimate_TTJ = plotTools.ConstructHDict(h_1D_TTJ_WJ[var],
                                                  name="TTJ estimate", color=rt.kRed+2,
                                                  title=htitle, drawoption="hist", 
                                                  appear_in_ratio="No", xtitle=var)
        hdict_estimate_QCD = plotTools.ConstructHDict(h_1D_QCD_TTJ_WJ[var],
                                                  name="QCD estimate", color=rt.kMagenta+1,
                                                  title=htitle, drawoption="hist", 
                                                  appear_in_ratio="No", xtitle=var)
        hdict_bg = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+region),
                                            name="Full MC BG estimate", color=rt.kRed+2,
                                            title=htitle, drawoption="E0", fillstyle=3002,
                                            appear_in_ratio="Yes", xtitle=var)
        #hdictlist=[hdict_estimate,hdict_bg,hdict_data]
        hdictlist=[hdict_estimate,hdict_estimate_QCD,hdict_estimate_TTJ,hdict_estimate_WJ,hdict_data]
        canvasname = var+"_comparison_data_estimate_"+region+postfix#+"_log"
        rtitle = "#frac{BG estimate}{Data}"
        plotTools.Plot1DWithRatio(hdictlist,outputdir,outfile,cname=canvasname,
                                  ratiotitle=rtitle,scale=False,legdict=leg,
                                  cdim=[696,550],ratiodim=0.3,logscale=False)
        
    
    outfile.Close()
