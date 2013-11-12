import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

# Do not actually need this function anymore...
def CloneHist1D(hold,name=""):
    hname = name
    if hname == "":
        hname = hold.GetName()+"_new"
    nx = hold.GetNbinsX()
    xmin = hold.GetXaxis().GetXmin()
    xmax = hold.GetXaxis().GetXmax()

    newh = rt.TH1D(hname,"",nx,xmin,xmax)
    for bin in range(hold.GetNbinsX()):
        newh.SetBinContent(bin,hold.GetBinContent(bin))
        newh.SetBinError(bin,hold.GetBinError(bin))
    
    return newh

def makeplots(hdict,outputdir,outfile,cname,ratiotitle):
    # placeholder for draw objects
    rootEvil = []

    # Place to store histogram that will be used as reference
    h_ref = None

    # Make the canvas, which will have two pads
    canvas = rt.TCanvas(cname,"")
    pad1 = rt.TPad("pad1","",0,0.25,1,1)
    pad1.SetBottomMargin(0)
    pad1.SetLogy(1)
    pad1.Draw()
    pad1.cd()

    # Make the legend
    legend = rt.TLegend(0.5,0.5,0.87,0.87,"")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)

    # Get histograms from a dictionary, and plot them 
    print "Getting all histograms"
    first = 0
    for name,stuff in hdict.iteritems():
        h = stuff[0] # Get the histogram
        sf = h.Integral(0,h.GetNbinsX()+1)
        h.Sumw2() # need to put this otherwise errors in ratio plot are wrong
        h.Scale(1./sf)
                
        if stuff[3]:
            print "Found ref histo"
            h_ref=h.Clone("h_ref")
        h.GetYaxis().SetTitle("A.U.")
        h.GetYaxis().SetTitleSize(0.055)
        h.GetYaxis().SetTitleOffset(0.8)
        h.GetYaxis().SetLabelSize(0.05)
        h.SetLineColor(stuff[1])
        h.SetLineWidth(2)
        h.SetTitle(stuff[2])
        legend.AddEntry(h,name,"l")
        h.SetMinimum(0.00005)
        drawoption = "HIST"
        if first > 0:
            drawoption = "HIST same"
        rootEvil.append(h.DrawClone(drawoption))
        first = first+1
        
    print "Drew all histograms"
    legend.Draw("same")
    
    # Make the second canvas, with the ratios
    canvas.cd()
    pad2 = rt.TPad("pad2","",0,0,1,0.25)
    pad2.SetBottomMargin(0.35)#0.25
    pad2.SetTopMargin(0)#0.05
    pad2.SetGridy(1)
    pad2.Draw()
    pad2.cd()

    # Get all the ratio plots
    print "Getting all ratios"
    firstB = 0
    for name,stuff in hdict.iteritems():
        if stuff[3]: continue
        h = stuff[0]
        ratio = h.Clone()
        #ratio.Sumw2()
        #ratio.Scale(1./ratio.Integral(0,ratio.GetNbinsX()+1))
        h_ref_ratio = h_ref.Clone()
        #h_ref_ratio.Sumw2()
        #h_ref_ratio.Scale(1./h_ref_ratio.Integral(0,h_ref_ratio.GetNbinsX()+1))
        ratio.Divide(h_ref_ratio)
        ratio.SetMarkerColor(stuff[1])
        ratio.SetMarkerStyle(7)
        ratio.SetTitle("")
        ratio.SetName("ratio"+str(firstB))
        ratio.GetYaxis().SetRangeUser(0,2.2)
        ratio.GetYaxis().SetNdivisions(4,8,0)
        ratio.GetYaxis().SetLabelSize(0.14)
        ratio.GetYaxis().SetTitleSize(0.1)
        ratio.GetYaxis().SetTitle(ratiotitle)
        ratio.GetYaxis().SetTitleOffset(0.22)
        ratio.GetXaxis().SetLabelSize(0.13)
        ratio.GetXaxis().SetTitleSize(0.15)
        ratio.GetXaxis().SetTickLength(0.1)
        ratio.GetXaxis().SetTitle(var)
        ratio.SetStats(0)
        drawoption="E P"
        if firstB > 0:
            drawoption=drawoption + " SAME"
        rootEvil.append(ratio.DrawClone(drawoption))
        firstB = firstB+1
        
    print "Drew the ratio plot"

    rt.SetOwnership(pad1, False) # to avoid seg fault
    rt.SetOwnership(pad2, False) # to avoid seg fault
    
    canvas.cd()
    canvas.SaveAs(outputdir+"/"+cname+".pdf")
    outfile.cd()
    canvas.Write()
    canvas.Close()

    
if __name__ == '__main__':

    if len(sys.argv) < 5:
        print "Run as: python %s <outputdir> <inputfile1> <inputfile2> <inputfile3>" % (sys.argv[0])
        sys.exit()
    outputdir = sys.argv[1]
    inputfile = sys.argv[2] # total bg histograms
    inputfile2 = sys.argv[3] # QCD only histograms
    inputfile3 = sys.argv[4] # TTJets only histograms

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/shapeplots.root","RECREATE")
    infile = TFile.Open(inputfile)
    infile2 = TFile.Open(inputfile2)
    infile3 = TFile.Open(inputfile3)

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD

    # set root styles
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetTextFont(42)
    rt.TH1.SetDefaultSumw2()
    
    vars = ["MR","R2"]

    # build hdict for TTJ control regions:
    for var in vars:
        hdict = {}
        hdict["TTJets CR"]                = [infile.Get("h_"+var+"_g1Mbg1W1Ll")     ,rt.kCyan+3,"Shape comparison for TTJets Control Regions",True]
        hdict["TTJets CR, 30 < mT < 100"] = [infile.Get("h_"+var+"_g1Mbg1W1LlmT")   ,rt.kCyan  ,"Shape comparison for TTJets Control Regions",False]
        hdict["TTJets CR, mT < 100"]      = [infile.Get("h_"+var+"_g1Mbg1W1LlmT100"),rt.kCyan+2,"Shape comparison for TTJets Control Regions",False]
        cname = var+"_comparison_TTJ_mT"
        ratiotitle = "#frac{mT cut}{no mT cut}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    # QCD control regions
    for var in vars:
        hdict = {}
        hdict["QCD CR"]                    = [infile.Get("h_"+var+"_0Lbg1uW0Ll")         ,rt.kMagenta  ,"Shape comparison for QCD Control Regions",True]
        hdict["QCD CR, minDeltaPhi < 0.3"] = [infile.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),rt.kMagenta+2,"Shape comparison for QCD Control Regions",False]
        cname = var+"_comparison_QCD_minDeltaPhi"
        ratiotitle = "#frac{minDeltaPhi cut}{no minDeltaPhi cut}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    # shape of QCD in Signal region vs shape of total QCD CR region
    for var in vars:
        hdict = {}
        hdict["QCD in SIG"]   = [infile2.Get("h_"+var+"_g1Mbg1W0Ll")        ,rt.kBlack    ,"Shape comparison for QCD in SIG vs total QCD CR",True]
        hdict["Total QCD CR"] = [infile.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),rt.kMagenta+2,"Shape comparison for QCD in SIG vs total QCD CR",False]
        cname = var+"_comparison_QCD"
        ratiotitle = "#frac{Total QCD CR}{QCD in SIG}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    # shape of TTJ in Signal region vs shape of total TTJ CR region
    for var in vars:
        hdict = {}
        hdict["TTJets in SIG"]   = [infile3.Get("h_"+var+"_g1Mbg1W0Ll")    ,rt.kBlack ,"Shape comparison for TTJets in SIG vs total TTJ CR",True]
        hdict["Total TTJets CR"] = [infile.Get("h_"+var+"_g1Mbg1W1LlmT100"),rt.kCyan+2,"Shape comparison for TTJets in SIG vs total TTJ CR",False]
        cname = var+"_comparison_TTJets"
        ratiotitle = "#frac{Total TTJets CR}{TTJets in SIG}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    # shape of QCD in Signal region vs shape of QCD in QCD CR region
    for var in vars:
        hdict = {}
        hdict["QCD in SIG"] = [infile2.Get("h_"+var+"_g1Mbg1W0Ll")         ,rt.kBlack    ,"Shape comparison for QCD, SIG vs CR",True]
        hdict["QCD in CR"]  = [infile2.Get("h_"+var+"_0Lbg1uW0Ll_mdPhi0p3"),rt.kMagenta+2,"Shape comparison for QCD, SIG vs CR",False]
        cname = var+"_comparison_QCD_SIGvsCR"
        ratiotitle = "#frac{QCD in CR}{QCD in SIG}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    # shape of TTJ in Signal region vs shape of TTJ in TTJ CR region
    for var in vars:
        hdict = {}
        hdict["TTJets in SIG"] = [infile3.Get("h_"+var+"_g1Mbg1W0Ll")     ,rt.kBlack ,"Shape comparison for TTJets, SIG vs CR",True]
        hdict["TTJets in CR"]  = [infile3.Get("h_"+var+"_g1Mbg1W1LlmT100"),rt.kCyan+2,"Shape comparison for TTJets, SIG vs CR",False]
        cname = var+"_comparison_TTJets_SIGvsCR"
        ratiotitle = "#frac{TTJets in CR}{TTJets in SIG}"
        makeplots(hdict,outputdir,outfile,cname,ratiotitle)

    outfile.Close()
