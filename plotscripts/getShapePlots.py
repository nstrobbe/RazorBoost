import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

SIGNAME = "T1ttcc_325_300"
#SIGNAME = "T2tt_400_50"

def gethistos(a,fs,cut):
    histos = {}
    for sample,f in fs.iteritems():
        histos[sample] = f.Get("h_"+a+"_"+cut)
    return histos

def gethistfromdict(histos,name):
    hist = 0
    for sample in histos.keys():
        if name in sample:
            hist = histos[sample]
    return hist

def makeplots(var,fs,bgcut,sigcut,signame,outputdir,outputfile,ext):
    histosMR = gethistos(var,fs,bgcut)
    h_TTj = gethistfromdict(histosMR,"TTJets").Clone()
    h_TTj.Sumw2()
    h_TTj.SetFillColor(rt.kRed)
    h_TTj.SetLineColor(rt.kRed)
    h_Wj = gethistfromdict(histosMR,"WJetsToLNu").Clone()
    h_Wj.Sumw2()
    h_Wj.SetFillColor(rt.kGreen+1)
    h_Wj.SetLineColor(rt.kGreen+1)
    h_VV = gethistfromdict(histosMR,"VV").Clone()
    h_VV.Sumw2()
    h_VV.SetFillColor(rt.kBlue)
    h_VV.SetLineColor(rt.kBlue)
    h_QCD = gethistfromdict(histosMR,"QCD").Clone()
    h_QCD.Sumw2()
    h_QCD.SetFillColor(rt.kMagenta)
    h_QCD.SetLineColor(rt.kMagenta)
    h_T = gethistfromdict(histosMR,"Top").Clone()
    h_T.Sumw2()
    h_T.SetFillColor(rt.kCyan)
    h_T.SetLineColor(rt.kCyan)
    h_zjets = gethistfromdict(histosMR,"ZJetsToNuNu").Clone()
    h_zjets.Sumw2()
    h_zjets.SetFillColor(rt.kOrange)
    h_zjets.SetLineColor(rt.kOrange)
    
    print "got histograms"

    scalefactor = (h_VV.Integral() + h_Wj.Integral() + h_TTj.Integral() + h_QCD.Integral() + h_T.Integral() + h_zjets.Integral())
    print "got scalefactor ", scalefactor
    h_VV.Scale(1./scalefactor)
    h_TTj.Scale(1./scalefactor)
    h_Wj.Scale(1./scalefactor)
    h_QCD.Scale(1./scalefactor)
    h_T.Scale(1./scalefactor)
    h_zjets.Scale(1./scalefactor)
    mc = rt.THStack()
    mc.Add(h_T, "hist")
    mc.Add(h_zjets, "hist")
    mc.Add(h_VV, "hist")
    mc.Add(h_TTj, "hist")
    mc.Add(h_Wj, "hist")
    mc.Add(h_QCD, "hist")
    
    htotal = h_QCD.Clone()
    htotal.Add(h_zjets)
    htotal.Add(h_T)
    htotal.Add(h_VV)
    htotal.Add(h_TTj)
    htotal.Add(h_Wj)
    htotal.Scale(1./htotal.Integral())

    sigs = gethistos(var,fs,sigcut)
    sig = gethistfromdict(sigs,signame).Clone()
    sig.Sumw2()
    sig.Scale(1./sig.Integral())
    #sig.SetMarkerStyle(20)
    sig.SetLineColor(rt.kBlack)
    
    legend = rt.TLegend(0.5,0.5,0.87,0.87,"")
    legend.AddEntry(sig," "+signame+" "+sigcut,"l")
    legend.AddEntry(h_QCD," QCD "+bgcut,"f")
    legend.AddEntry(h_TTj," TTjets "+bgcut,"f")
    legend.AddEntry(h_Wj," Wjets "+bgcut,"f")
    legend.AddEntry(h_VV," VV "+bgcut,"f")
    legend.AddEntry(h_zjets," Zinv "+bgcut,"f")
    legend.AddEntry(h_T," T "+bgcut,"f")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    
    canvas = rt.TCanvas(var+"_"+bgcut+"_"+signame+ext,"")
    pad1 = rt.TPad("pad1","",0,0.25,1,1)
    pad1.SetBottomMargin(0)
    pad1.SetLogy(1)
    pad1.Draw()
    pad1.cd()
    mc.Draw()
    mc.GetYaxis().SetTitle("A.U.")
    mc.GetYaxis().SetTitleSize(0.055)
    mc.GetYaxis().SetTitleOffset(0.8)
    mc.GetYaxis().SetLabelSize(0.05)
    mc.SetTitle("Control region " + bgcut)
    #pad1.Update()
    #print mc.GetXaxis()
    #mc.GetXaxis().SetTitle(var)
    mc.SetMinimum(0.00005)
    sig.Draw("HISTsame")
    legend.Draw("same")
    canvas.cd()
    pad2 = rt.TPad("pad2","",0,0,1,0.25)
    pad2.SetBottomMargin(0.35)#0.25
    pad2.SetTopMargin(0)#0.05
    pad2.SetGridy(1)
    pad2.Draw()
    pad2.cd()
    ratio = sig.Clone()
    ratio.Divide(htotal)
    ratio.SetTitle("")
    ratio.SetName("histoRatio")
    ratio.GetYaxis().SetRangeUser(0,2.2)
    ratio.GetYaxis().SetNdivisions(4,8,0)
    ratio.GetYaxis().SetLabelSize(0.14)
    ratio.GetYaxis().SetTitleSize(0.15)
    ratio.GetYaxis().SetTitle("SR/CR ")
    ratio.GetYaxis().SetTitleOffset(0.22)
    ratio.GetXaxis().SetLabelSize(0.13)
    ratio.GetXaxis().SetTitleSize(0.15)
    ratio.GetXaxis().SetTickLength(0.1)
    ratio.GetXaxis().SetTitle(var)
    ratio.SetStats(0)
    ratio.Draw("pe")

    rt.SetOwnership(pad1, False) # to avoid seg fault
    rt.SetOwnership(pad2, False) # to avoid seg fault

    canvas.cd()
    canvas.SaveAs(outputdir+"/"+var+"_"+bgcut+"_"+signame+ext+".pdf")
    outfile.cd()
    canvas.Write()


    
if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "Run as: python %s <preselection> <outputdir>" % (sys.argv[0])
        sys.exit()
    #samplesfile = sys.argv[1]
    outputdir = sys.argv[2]
    presel = sys.argv[1]

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/plots.root","RECREATE")

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    #intlumi = 5.126 # AB

    datasets = ["QCD","TTJets","WJetsToLNu","VV","Top","ZJetsToNuNu","data",SIGNAME]
    fs = {}
    for d in datasets:
        f = TFile.Open(presel+"/summary/rzrBTevsel_plots_"+d+".root")
        fs[d]=f

    # set root styles
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetTextFont(42)
        
    #make plots
    varis = ["MR","R2"]
    regions = ["Full","SB","SIG"]
    cuts = ["g1Mbg1W0Ll","g1Mbg1W0Ll_Wbveto","g1Mbg1W0Ll_Wbveto_MinMassDiff100","g1Mbg1W1Ll","g1Mbg1W1Ll_Wbveto","g1Mbg1W1Ll_Wbveto_MinMassDiff100","0Lbg1W0Ll","g1uW0Ll","0Lb0Ll","g1Mbg1uW0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_MinMassDiff100"]
    sigcut = "g1Mbg1W0Ll"
    sigcut_Wbveto = "g1Mbg1W0Ll_Wbveto"
    sigcut_Wbveto_MinMassDiff100 = "g1Mbg1W0Ll_Wbveto_MinMassDiff100"
    for var in varis:
        for region in regions:
            for cut in cuts:
                nm = region + "_" + cut
                if cut is "":
                    nm = region
                #makeplots(var,fs,nm,region+"_"+sigcut,"TTJets",outputdir,outfile,"")
                #makeplots(var,fs,nm,region+"_"+sigcut,"QCD",outputdir,outfile,"")
                #makeplots(var,fs,nm,region+"_"+sigcut_Wbveto,"TTJets",outputdir,outfile,"_Wbveto")
                #makeplots(var,fs,nm,region+"_"+sigcut_Wbveto,"QCD",outputdir,outfile,"_Wbveto")
                #makeplots(var,fs,nm,region+"_"+sigcut_Wbveto_MinMassDiff100,"TTJets",outputdir,outfile,"_Wbveto_MinMassDiff100")
                #makeplots(var,fs,nm,region+"_"+sigcut_Wbveto_MinMassDiff100,"QCD",outputdir,outfile,"_Wbveto_MinMassDiff100")

    # data/MC comparison of control regions
    cuts2 = ["g1Mbg1W0Ll","g1Mbg1W0Ll_Wbveto","g1Mbg1W0Ll_Wbveto_MinMassDiff100","g1Mbg1W1Ll","g1Mbg1W1Ll_Wbveto","g1Mbg1W1Ll_Wbveto_MinMassDiff100","0Lbg1W0Ll","g1uW0Ll","0Lb0Ll","g1Mbg1uW0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_MinMassDiff100"]
    for var in varis:
        for region in regions:
            for cut in cuts2:
                print var, region, cut
                makeDataMCplots(var,fs,region+"_"+cut,outputdir,outfile,intlumi)
                makeDataMCplotsScaled(var,fs,region+"_"+cut,outputdir,outfile,intlumi)

    print "now make njets plots"

    for region in regions:
        makeDataMCplots("njets",fs,region+"_g1Mbg1W0Ll_Wbveto",outputdir,outfile,intlumi)
        makeDataMCplotsScaled("njets",fs,region+"_g1Mbg1W0Ll_Wbveto",outputdir,outfile,intlumi)
        makeDataMCplots("njets",fs,region+"_g1Mbg1W0Ll_Wbveto_MinMassDiff100",outputdir,outfile,intlumi)
        makeDataMCplotsScaled("njets",fs,region+"_g1Mbg1W0Ll_Wbveto_MinMassDiff100",outputdir,outfile,intlumi)

    outfile.Close()

