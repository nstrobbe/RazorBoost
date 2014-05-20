from ROOT import *
import os
import plotTools

def getGluinoCrossSection(mgluino,f):
    h = f.Get("gluino8TeV_NLONLL")
    binnr = h.GetXaxis().FindBin(mgluino)
    xsec = h.GetBinContent(binnr)
    return xsec

def getStopCrossSection(mstop,f):
    h = f.Get("stop8TeV_NLONLL")
    binnr = h.GetXaxis().FindBin(mstop)
    xsec = h.GetBinContent(binnr)
    return xsec

if __name__ == "__main__":
    gStyle.SetOptStat(0)
    gStyle.SetLabelSize(0.045,"xyz")
    gStyle.SetTitleSize(0.06,"xyz")
    gROOT.ForceStyle()
    plotTools.SetColorPaletteSMS()

    dirbase = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results"
    indir = dirbase + "/results_20140325_T2tt/"
    bfile = TFile.Open(dirbase+"/results_20140325_noISR_btag_TopPt/rzrBoostMC_bg.root")
    counts = bfile.Get("counts")

    outdir = dirbase + "/results_20140325_T2tt/"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
            
    crosssectionfilename = "referenceXSecs.root"
    cf = TFile.Open(crosssectionfilename)
    
    #DMs = ["T1ttcc_DM-10","T1ttcc_DM-25","T1ttcc_DM-80"]
    DMs = ["T2tt"]
    #DMs = ["T1ttcc_old"]

    for DM in DMs:
        infile = TFile.Open(indir+"/rzrBoostMC_SMS_"+DM+".root")
        
        outfile = TFile.Open(outdir+"/"+DM+".root","RECREATE")
        
        #h_nevents = infile.Get("h_mstop_mLSP_Pileup")
        h_nevents = infile.Get("h_mstop_mLSP_ISR")
        
        basename = "h_mstop_mLSP"
        selections = ["HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
                      "SIG","neleeq0","nmueq0","trackIso",
                      "g1Mb0Ll","g1Mbg1W0Ll","1Mbg1W0Ll","g2Mbg1W0Ll",
                      "g1Mbg1W0Ll_mdPhiHat4","g1Mbg1W0Ll_mdPhiHatg4",
                      "g1Mb0Wg1uW0Ll",
                      "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5","0Lbg1W0Ll",
                      "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100","1Mbg1W1LlmT100","g2Mbg1W1LlmT100",
                      "g1Mbg1W1LlmT",
                      "0Lb1Ll","0Lbg1Y1Ll","0Lbg1Y1LlmT100","0Lbg1Y1LlmT",
                      "2munoZmass","2mu","2mu0el","0Lb2mu0el","0Lbg1Y2mu0el","g1Mb2mu0el","g1Mbg1Y2mu0el",
                      "2elnoZmass","2el","2el0mu","0Lb2el0mu","0Lbg1Y2el0mu","g1Mb2el0mu","g1Mbg1Y2el0mu",
                      "2lnoZmass","2l","2l0ol","0Lb2l0ol","0Lbg1Y2l0ol","g1Mb2l0ol","g1Mbg1Y2l0ol",
                      ]

        for selection in selections:
            # make plots
            hname = basename + "_" + selection
            print hname
            h = infile.Get(hname)
            print h
            # efficiency
            h_eff = h.Clone("eff")
            h_eff.Divide(h_nevents)
            if "T1ttcc" in DM:
                h_eff.GetXaxis().SetTitle("mGluino")
            if "T1ttcc_old" in DM:
                h_eff.GetXaxis().SetTitle("mStop")
            if "T2tt" in DM:
                h_eff.GetXaxis().SetTitle("mStop")
            h_eff.GetYaxis().SetTitle("mLSP")
            h_eff.GetZaxis().SetTitle("Efficiency")
            h_eff.SetTitle("Efficiency for "+DM)
            c_eff = TCanvas("efficiency_"+DM+"_"+selection)
            c_eff.cd()
            c_eff.SetRightMargin(0.21)
            c_eff.SetLeftMargin(0.14)
            c_eff.SetBottomMargin(0.14)
            h_eff.DrawCopy("colz")
            c_eff.SaveAs(outdir+"/efficiency_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_eff.Write()
            c_eff.Close()
            
            # #events
            h_events = h_eff.Clone("h_events")
            for x in range(1,h_eff.GetNbinsX()+1):
                for y in range(1,h_eff.GetNbinsY()+1):
                    eff = h_eff.GetBinContent(x,y)
                    if "T1ttcc_old" in DM:
                        sf = 0.0243547
                    elif "T1ttcc" in DM:
                        sf = getGluinoCrossSection(h_eff.GetXaxis().GetBinLowEdge(x),cf)
                    elif "T2tt" in DM:
                        sf = getStopCrossSection(h_eff.GetXaxis().GetBinLowEdge(x),cf)
                    #print h_eff.GetXaxis().GetBinLowEdge(x), sf
                    h_events.SetBinContent(x,y,eff*sf)
                    
            h_events.Scale(19789)
            if "T1ttcc" in DM:
                h_events.GetXaxis().SetTitle("mGluino")
            if "T1ttcc_old" in DM:
                h_eff.GetXaxis().SetTitle("mStop")
            if "T2tt" in DM:
                h_events.GetXaxis().SetTitle("mStop")
            h_events.GetYaxis().SetTitle("mLSP")
            h_events.GetZaxis().SetTitle("Events")
            h_events.SetTitle("Number of events for "+DM+" (19.789 fb-1)")
            c_events = TCanvas("events_"+DM+"_"+selection)
            c_events.cd()
            c_events.SetRightMargin(0.21)
            c_events.SetLeftMargin(0.14)
            c_events.SetBottomMargin(0.14)
            c_events.SetLogz()
            h_events.DrawCopy("colz")
            c_events.SaveAs(outdir+"/events_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_events.Write()
            c_events.Close()
            
            # significance
            
            binnumber = counts.GetXaxis().FindBin(selection)
            B = counts.GetBinContent(binnumber)
            print B
            h_signif = h_events.Clone("signif")
            h_signif.Scale(1./TMath.Sqrt(B))
            if "T1ttcc" in DM:
                h_signif.GetXaxis().SetTitle("mGluino")
            if "T1ttcc_old" in DM:
                h_eff.GetXaxis().SetTitle("mStop")
            if "T2tt" in DM:
                h_signif.GetXaxis().SetTitle("mStop")
            h_signif.GetYaxis().SetTitle("mLSP")
            h_signif.GetZaxis().SetTitle("#frac{S}{#sqrt{B}}")
            h_signif.SetTitle("significance S/sqrt(B) for "+DM)
            c_signif = TCanvas("significance_"+DM+"_"+selection)
            c_signif.cd()
            c_signif.SetRightMargin(0.21)
            c_signif.SetLeftMargin(0.14)
            c_signif.SetBottomMargin(0.14)
            h_signif.DrawCopy("colz")
            c_signif.SaveAs(outdir+"/significance_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_signif.Write()
            c_signif.Close()
            
        infile.Close()
        outfile.Close()
    
        # print region, selection, h_eff.GetBinContent(1), h_events.GetBinContent(1), B, h_signif.GetBinContent(1)
    
    bfile.Close()
    cf.Close()
