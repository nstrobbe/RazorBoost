import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

SIGNAME = "T1ttcc_325_300"
#SIGNAME = "T2tt_400_50"

def gethistos(a,fs,cut):
    histos = {}
    for sample,f in fs.iteritems():
        # print sample,f,"h_"+a+"_"+cut
        histos[sample] = f.Get("h_"+a+"_"+cut)
    return histos

def gethistfromdict(histos,name):
    hist = 0
    for sample in histos.keys():
        if name in sample:
            hist = histos[sample]
    return hist

def makeDataMCplot(var,fs,bgcut,outputdir,outputfile,intlumi,scaled=False):
    # get all the histograms out of the dictionary
    histosMR = gethistos(var,fs,bgcut)
    h_QCD = gethistfromdict(histosMR,"QCD").Clone()
    h_QCD.Sumw2()
    h_QCD.SetFillColor(rt.kMagenta)
    h_QCD.SetLineColor(rt.kMagenta)
    h_TTj = gethistfromdict(histosMR,"TTJets").Clone()
    h_TTj.Sumw2()
    h_TTj.SetFillColor(rt.kRed)
    h_TTj.SetLineColor(rt.kRed)
    h_T = gethistfromdict(histosMR,"Top").Clone()
    h_T.Sumw2()
    h_T.SetFillColor(rt.kCyan)
    h_T.SetLineColor(rt.kCyan)
    h_TTX = gethistfromdict(histosMR,"TTX").Clone()
    h_TTX.Sumw2()
    h_TTX.SetFillColor(rt.kCyan+2)
    h_TTX.SetLineColor(rt.kCyan+2)
    h_Wj = gethistfromdict(histosMR,"WJetsToLNu").Clone()
    h_Wj.Sumw2()
    h_Wj.SetFillColor(rt.kGreen+1)
    h_Wj.SetLineColor(rt.kGreen+1)
    h_VV = gethistfromdict(histosMR,"VV").Clone()
    h_VV.Sumw2()
    h_VV.SetFillColor(rt.kBlue+1)
    h_VV.SetLineColor(rt.kBlue+1)
    h_VVV = gethistfromdict(histosMR,"VVV").Clone()
    h_VVV.Sumw2()
    h_VVV.SetFillColor(rt.kBlue-3)
    h_VVV.SetLineColor(rt.kBlue-3)
    h_zjets = gethistfromdict(histosMR,"ZJetsToNuNu").Clone()
    h_zjets.Sumw2()
    h_zjets.SetFillColor(rt.kOrange)
    h_zjets.SetLineColor(rt.kOrange)
    h_zll = gethistfromdict(histosMR,"DYJetsToLL").Clone()
    h_zll.Sumw2()
    h_zll.SetFillColor(rt.kOrange+2)
    h_zll.SetLineColor(rt.kOrange+2)
    h_sig = gethistfromdict(histosMR,SIGNAME).Clone()
    h_sig.Sumw2()
    h_sig.SetFillColor(rt.kGray)
    h_sig.SetLineColor(rt.kGray)
    print "got histograms"

    # Add histograms to the stack, starting with the smallest one
    mc = rt.THStack()
    mc.Add(h_VVV,"hist")
    mc.Add(h_VV, "hist")
    mc.Add(h_zll,"hist")
    mc.Add(h_zjets, "hist")
    mc.Add(h_TTX,"hist")
    mc.Add(h_T, "hist")
    mc.Add(h_Wj, "hist")
    mc.Add(h_TTj, "hist")
    mc.Add(h_QCD, "hist")
    mc.Add(h_sig, "hist")
    
    # Make total histogram with only the backgrounds, 
    # will be used for the ratio plot
    htotal = h_QCD.Clone()
    htotal.Add(h_VV)
    htotal.Add(h_VVV)
    htotal.Add(h_T)
    htotal.Add(h_TTX)
    htotal.Add(h_zll)
    htotal.Add(h_zjets)
    htotal.Add(h_TTj)
    htotal.Add(h_Wj)

    # Get the data histogram
    data = gethistfromdict(histosMR,"data").Clone()
    data.SetMarkerStyle(20)
    data.SetLineColor(rt.kBlack)

    # scale mc to data if requested
    if scaled:
        # get the total mc count, including underflow and overflow bins
        mc_int = htotal.Integral(0,htotal.GetNbinsX()+1)
        if mc_int == 0:
            mc_int = 1.
        data_int = data.Integral()
        sf = data_int/mc_int
        h_QCD.Scale(sf)
        h_TTj.Scale(sf)
        h_T.Scale(sf)
        h_TTX.Scale(sf)
        h_Wj.Scale(sf)
        h_VV.Scale(sf)
        h_VVV.Scale(sf)
        h_zjets.Scale(sf)
        h_zll.Scale(sf)
        
    
    # Make the legend, ordering histograms according to the stack
    legend = rt.TLegend(0.63,0.4,0.87,0.8,"")
    legend.AddEntry(data," data ("+ str(intlumi) + "/fb)","epl")
    legend.AddEntry(h_QCD," QCD ","f")
    legend.AddEntry(h_TTj," TTjets ","f")
    legend.AddEntry(h_Wj," Wjets ","f")
    legend.AddEntry(h_T," T ","f")
    legend.AddEntry(h_TTX," TTX ","f")
    legend.AddEntry(h_zjets," Zinv ","f")
    legend.AddEntry(h_zll," DYToLL ","f")
    legend.AddEntry(h_VV," VV ","f")
    legend.AddEntry(h_VVV," VVV ","f")
    legend.AddEntry(h_sig," %s "%(SIGNAME),"f")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    
    # Get the maximum of the histograms, so that we can set the Y-axis range
    maxi = data.GetMaximum()
    if htotal.GetMaximum() > maxi:
        maxi = htotal.GetMaximum()

    # String to distinguish output files
    sca = ""
    if scaled:
        sca = "Scaled_"

    # Make the  canvas
    canvas = rt.TCanvas("DataMC_"+sca+var+"_"+bgcut+"_withSignal","")
    pad1 = rt.TPad("pad1","",0,0.25,1,1)
    pad1.SetBottomMargin(0)
    pad1.SetLogy(1)
    pad1.Draw()
    pad1.cd()
    data.GetYaxis().SetTitle("Events")
    data.GetYaxis().SetTitleSize(0.055)
    data.GetYaxis().SetTitleOffset(0.8)
    data.GetYaxis().SetLabelSize(0.05)
    data.SetTitle("Data-MC Comparison ")
    if scaled:
        data.SetTitle("Data-MC Comparison with MC scaled to match data ")
    data.SetMaximum(2.*maxi)
    data.SetMinimum(0.007)
    data.Draw("EP")
    mc.Draw("same")
    data.Draw("EPsame")
    legend.Draw("same")
    t = '#scale[1.1]{Selection ' + bgcut + '}'
    tex = rt.TLatex(0.52,0.82,t)
    tex.SetNDC();
    tex.Draw("same");

    # Make the ratioplot in pad2
    canvas.cd()
    pad2 = rt.TPad("pad2","",0,0,1,0.25)
    pad2.SetBottomMargin(0.35)
    pad2.SetTopMargin(0)
    pad2.SetGridy(1)
    pad2.Draw()
    pad2.cd()
    ratio = data.Clone()
    ratio.Divide(htotal)
    ratio.SetTitle("")
    ratio.SetName("histoRatio")
    ratio.GetYaxis().SetRangeUser(0,2.2)
    ratio.GetYaxis().SetNdivisions(4,8,0)
    ratio.GetYaxis().SetLabelSize(0.14)
    ratio.GetYaxis().SetTitleSize(0.15)
    ratio.GetYaxis().SetTitle("Data/MC ")
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
    cname = "DataMC_"+sca+var+"_"+bgcut+"_withSignal_"+SIGNAME+".pdf"
    canvas.SaveAs(outputdir+"/"+cname)
    outfile.cd()
    canvas.Write()

    
if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Run as: python %s  <outputdir>" % (sys.argv[0])
        sys.exit()
    outputdir = sys.argv[1]

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/plots.root","RECREATE")

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    #intlumi = 5.126 # AB

    print "Will make plots for integrated luminosity of %.3f fb-1" % (intlumi)

    datasets = ["QCD","TTJets","WJetsToLNu","VV","Top","ZJetsToNuNu","DYJetsToLL","TTX","VVV","dataSingleMu",SIGNAME]
    fs = {}
    for d in datasets:
        f = TFile.Open("results_BoostMC_1/summary/rzrBoostMC_"+d+".root")
        fs[d]=f

    # set root styles
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetTextFont(42)
        
    # make plots
    varis = ["MR","R2"]
    cuts = ["Cleaning","HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
            "SIG","neleeq0","nmueq0","trackIso",
            "g1Mb0Ll","g1Mbg1W0Ll","0Lb0Ll","0Lbg1uW0Ll","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll", 
            "2mu","2mu0el","0Lb2mu0el","g1Mb2mu0el" 
            ]

    for var in varis:
        for cut in cuts:
            print var, cut
            makeDataMCplot(var,fs,cut,outputdir,outfile,intlumi,False) # Scaled to cross section
            makeDataMCplot(var,fs,cut,outputdir,outfile,intlumi,True) # Scaled to match data

    outfile.Close()

