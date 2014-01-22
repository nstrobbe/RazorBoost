from ROOT import *
import plotTools

def getGluinoCrossSection(mgluino,f):
    h = f.Get("gluino8TeV_NLONLL")
    binnr = h.GetXaxis().FindBin(mgluino)
    xsec = h.GetBinContent(binnr)
    return xsec

if __name__ == "__main__":
    gStyle.SetOptStat(0)

    plotTools.SetColorPaletteSMS()

    outdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results"
    bfile = TFile.Open(outdir+"/results_20140122/rzrBoostMC_bg.root")
    counts = bfile.Get("counts")

    crosssectionfilename = "referenceXSecs.root"
    cf = TFile.Open(crosssectionfilename)
    
    DMs = ["DM-10","DM-25","DM-80"]

    for DM in DMs:
        infile = TFile.Open(outdir+"/results_20140122_T1ttcc/rzrBoostMC_SMS_T1ttcc_"+DM+".root")
        
        outfile = TFile.Open(outdir+"/plots_20140122/T1ttcc_"+DM+".root","RECREATE")
        
        h_nevents = infile.Get("h_mstop_mLSP_Pileup")
        
        basename = "h_mstop_mLSP"
        selections = ["HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
                      "SIG","neleeq0","nmueq0","trackIso",
                      "g1Mb0Ll","g1Mbg1W0Ll","g1Mb0Wg1uW0Ll",
                      "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5","0Lbg1W0Ll",
                      "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100","g1Mbg1W1LlmT",
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
            h_eff = h.Clone()
            h_eff.Divide(h_nevents)
            h_eff.GetXaxis().SetTitle("mGluino")
            h_eff.GetYaxis().SetTitle("mLSP")
            h_eff.SetTitle("Efficiency for T1ttcc "+DM)
            c_eff = TCanvas("efficiency_T1ttcc_"+DM+"_"+selection)
            c_eff.cd()
            c_eff.SetRightMargin(0.15)
            h_eff.Draw("colz")
            c_eff.SaveAs(outdir+"/plots_20140122/efficiency_T1ttcc_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_eff.Write()
            
            # #events
            h_events = h_eff.Clone("h_events")
            for x in range(1,h_eff.GetNbinsX()+1):
                for y in range(1,h_eff.GetNbinsY()+1):
                    eff = h_eff.GetBinContent(x,y)
                    sf = getGluinoCrossSection(h_eff.GetXaxis().GetBinLowEdge(x),cf)
                    #print h_eff.GetXaxis().GetBinLowEdge(x), sf
                    h_events.SetBinContent(x,y,eff*sf)
                    
            h_events.Scale(19789)
            h_events.GetXaxis().SetTitle("mGluino")
            h_events.GetYaxis().SetTitle("mLSP")
            h_events.SetTitle("Number of events for T1ttcc "+DM+" (19.789 fb-1)")
            c_events = TCanvas("events_T1ttcc_"+DM+"_"+selection)
            c_events.cd()
            c_events.SetRightMargin(0.15)
            c_events.SetLogz()
            h_events.Draw("colz")
            c_events.SaveAs(outdir+"/plots_20140122/events_T1ttcc_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_events.Write()
        
            # significance
            
            binnumber = counts.GetXaxis().FindBin(selection)
            B = counts.GetBinContent(binnumber)
            print B
            h_signif = h_events.Clone()
            h_signif.Scale(1./TMath.Sqrt(B))
            h_signif.GetXaxis().SetTitle("mGluino")
            h_signif.GetYaxis().SetTitle("mLSP")
            h_signif.SetTitle("significance S/sqrt(B) for T1ttcc "+DM)
            c_signif = TCanvas("significance_T1ttcc_"+DM+"_"+selection)
            c_signif.cd()
            c_signif.SetRightMargin(0.15)
            h_signif.Draw("colz")
            c_signif.SaveAs(outdir+"/plots_20140122/significance_T1ttcc_"+DM+"_"+selection+".pdf")
            outfile.cd()
            c_signif.Write()
            
        infile.Close()
        outfile.Close()
    
        # print region, selection, h_eff.GetBinContent(1), h_events.GetBinContent(1), B, h_signif.GetBinContent(1)
    
    bfile.Close()
    cf.Close()
