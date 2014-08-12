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
    outputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/plots_20140730_preApp_comments"
    # outputdir = "./test_2D/"
    basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140730_preApp_comments/summary/"
    #basedir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass/summary/"
    #basedir_old = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140331_noISR_btag_TopPt/summary/"
    # inputfile_TTJ = basedir + "rzrBoostMC_TTJets.root"
    inputfile_data = basedir + "rzrBoostMC_data.root"
    #inputfile_data_old = basedir_old + "rzrBoostMC_data.root"
    inputfile_bg = basedir + "../rzrBoostMC_bg.root"
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/2Dplots2.root","RECREATE")
    # infile_TTJ = TFile.Open(inputfile_TTJ)
    infile_data = TFile.Open(inputfile_data)
    #infile_data_old = TFile.Open(inputfile_data_old)
    infile_bg = TFile.Open(inputfile_bg)

    # Integrated luminosity in fb-1s
    intlumi = 19.712 # ABCD

    # set root styles
    plotTools.SetBoostStyle()
    
    vars = ["MR_R2"]

    # comparison plots to be made:
    # 1. comparison between process in signal region and relevant control region
    #    * QCD
    #    * TTJets
    #    * WJets
    #    * ZJets
    # 2. comparison between different control regions
    #    * different QCD regions
    #    * different TTjets regions

    # build hdictlist for TTJ:
    #for var in vars:
    #    hdict_TTJ = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W1Ll"),
    #                                         name="TTJets CR", title="R2 vs MR for TTJets in TTJets Control region",
    #                                         xtitle=var.split("_")[0], ytitle=var.split("_")[1],
    #                                         drawoption="colz", palette="SMS")
    #    
    #    hdict_SIG = plotTools.ConstructHDict(infile_TTJ.Get("h_"+var+"_g1Mbg1W0Ll"),
    #                                         name="SIG", title="R2 vs MR for TTJets in SIG Control region",
    #                                         xtitle=var.split("_")[0], ytitle=var.split("_")[1],
    #                                         drawoption="colz", palette="SMS")
    #    
    #    canvasname = var+"_TTJ"
    #    plotTools.Plot2D(hdict_TTJ,outputdir,outfile,cname=canvasname,scale="No",logscale=True)

    #    canvasname = var+"_ratio_TTJ_SIG"
    #    plotTools.Plot2DRatio(hdict_TTJ,hdict_SIG,outputdir,outfile,cname=canvasname,scale="Yes",ctitle="Ratio TTJ CR / SIG region for TTJets MC",ztitle="TTJ/SIG")

    cuts = ["g1Mbg1W1LlmT100_mdPhig0p5","0Lbg1Y1LlmT_mdPhig0p5",
            "g1Mbg1W1LlmT100_mdPhi0p5","0Lbg1Y1LlmT_mdPhi0p5"]
    names = ["T mdPhi > 0.5", "W mdPhi > 0.5","T mdPhi < 0.5","W mdPhi < 0.5"]
    for var in vars:
        for i,cut in enumerate(cuts):
            # total BG
            hdict_bg = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+cut), 
                                                name=names[i]+" CR", title="R2 vs MR distribution for the total background in the %s region"%(names[i]),
                                                xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                drawoption="colztext", palette="SMS") 
            canvasname = var+"_"+cut+"_bg"
            plotTools.Plot2D(hdict_bg,outputdir,outfile,cname=canvasname,scale="No",logscale=True)
            
            # data
            hdict_data = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+cut), 
                                                  name=names[i]+" CR", title="R2 vs MR distribution for the data in the %s region"%(names[i]),
                                                  xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                  drawoption="colztext", palette="SMS") 
            canvasname = var+"_"+cut+"_data"
            plotTools.Plot2D(hdict_data,outputdir,outfile,cname=canvasname,scale="No",logscale=True)
                        
            # ratio
            canvasname = var+"_"+cut+"_ratio_data_bg"
            plotTools.Plot2DRatio(hdict_data,hdict_bg,outputdir,outfile,cname=canvasname,scale="No",
                                  ctitle="Ratio Data / MC for %s region"%(names[i]),ztitle="Data/MC")


    #############################################################
    ## normal 2D plots
    #############################################################
    vars = ["MR_minDeltaPhiHat","R2_minDeltaPhiHat",
            "MR_minDeltaPhi","R2_minDeltaPhi"]
    cuts = ["0Lbg1Y1LlmT","g1Mbg1W1LlmT100"]
    names = ["W","T"]
    for var in vars:
        for i,cut in enumerate(cuts):
            # total BG
            hdict_bg = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+cut), 
                                                name=names[i]+" CR", 
                                                title="%s vs %s distribution for the total background in the %s region"%(var.split("_")[1],var.split("_")[0],names[i]),
                                                xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                drawoption="colztext", palette="SMS") 
            canvasname = var+"_"+cut+"_bg"
            plotTools.Plot2D(hdict_bg,outputdir,outfile,cname=canvasname,scale="No",logscale=True)
            
            # data
            hdict_data = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+cut), 
                                                  name=names[i]+" CR", 
                                                  title="%s vs %s distribution for the data in the %s region"%(var.split("_")[1],var.split(\
"_")[0],names[i]),
                                                  xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                  drawoption="colztext", palette="SMS") 
            canvasname = var+"_"+cut+"_data"
            plotTools.Plot2D(hdict_data,outputdir,outfile,cname=canvasname,scale="No",logscale=True)
                        
            # ratio
            canvasname = var+"_"+cut+"_ratio_data_bg"
            plotTools.Plot2DRatio(hdict_data,hdict_bg,outputdir,outfile,cname=canvasname,scale="No",
                                  ctitle="Ratio Data / MC for %s region"%(names[i]),ztitle="Data/MC")

    #############################################################
    ## Shape comparisons
    #############################################################
    vars = ["MR_R2"]
    cuts = ["g1Mbg1W1LlmT100","0Lbg1Y1LlmT"]
    for var in vars:
        for cut in cuts:
            # mdphi > 0.5
            hdict_mdphihatg4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+cut+"_mdPhig0p5"), 
                                                        name=cut+" CR", title="R2 vs MR distribution for the total background in the %s region"%(cut),
                                                        xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                        drawoption="colztext", palette="SMS") 

            # mdphi < 0.5
            hdict_mdphihat4 = plotTools.ConstructHDict(infile_bg.Get("h_"+var+"_"+cut+"_mdPhi0p5"), 
                                                       name=cut+" CR", title="R2 vs MR distribution for the total background in the %s region"%(cut),
                                                       xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                       drawoption="colztext", palette="SMS") 
        
            # ratio
            canvasname = var+"_"+cut+"_ratio_mdphi"
            plotTools.Plot2DRatio(hdict_mdphihatg4,hdict_mdphihat4,outputdir,outfile,cname=canvasname,scale="Yes",
                                  ctitle="Ratio minDeltaPhi > 0.5 / minDeltaPhi < 0.5 for %s region"%(cut),ztitle="Ratio")

    for var in vars:
        for cut in cuts:
            # mdphi > 0.5
            hdict_mdphihatg4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+cut+"_mdPhig0p5"), 
                                                        name=cut+" CR", title="R2 vs MR distribution for data in the %s region"%(cut),
                                                        xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                        drawoption="colztext", palette="SMS") 

            # mdphi < 0.5
            hdict_mdphihat4 = plotTools.ConstructHDict(infile_data.Get("h_"+var+"_"+cut+"_mdPhi0p5"), 
                                                       name=cut+" CR", title="R2 vs MR distribution for data in the %s region"%(cut),
                                                       xtitle=var.split("_")[0], ytitle=var.split("_")[1],
                                                       drawoption="colztext", palette="SMS") 
        
            # ratio
            canvasname = var+"_"+cut+"_ratio_mdphi_data"
            plotTools.Plot2DRatio(hdict_mdphihatg4,hdict_mdphihat4,outputdir,outfile,cname=canvasname,scale="Yes",
                                  ctitle="Ratio minDeltaPhi > 0.5 / minDeltaPhi < 0.5 for %s region"%(cut),ztitle="Ratio")

    #############################
    # compare old vs new tagger #
    #############################

    # old tagger
    #hdict_old = plotTools.ConstructHDict(infile_data_old.Get("h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4"), 
    #                                     name="Old tagger", title="R2 vs MR distribution for data in the Signal region",
    #                                     xtitle="MR", ytitle="R2",
    #                                     drawoption="colztext", palette="SMS") 

    # new tagger
    #hdict_new = plotTools.ConstructHDict(infile_data.Get("h_MR_R2_g1Mbg1W0Ll_mdPhiHatg4"), 
    #                                     name="New tagger", title="R2 vs MR distribution for data in the Signal region",
    #                                     xtitle="MR", ytitle="R2",
    #                                     drawoption="colztext", palette="SMS") 
        
    # ratio
    #canvasname = "Wtagger_old_new_ratio_data"
    #plotTools.Plot2DRatio(hdict_new,hdict_old,outputdir,outfile,cname=canvasname,scale="No",
    #                      ctitle="Ratio New tagger / Old tagger for signal region",ztitle="New/Old")



    outfile.Close()
    # infile_TTJ.Close()
    infile_data.Close()
    #infile_data_old.Close()
    infile_bg.Close()
