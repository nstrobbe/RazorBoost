import sys, os
from string import *
import ROOT as rt
from ROOT import TFile

def make2Dplot(varx,vary,cut,infile,outfile,process):
    hname = "h_%s_%s_%s" % (varx,vary,cut)
    h = infile.Get(hname)
    h.GetYaxis().SetTitle(vary)
    h.GetXaxis().SetTitle(varx)
    h.GetZaxis().SetTitle("Events")
    h.GetZaxis().SetTitleOffset(1.5)
    h.SetMaximum(10000)
    h.SetTitle("2D plot of %s vs %s for %s"%(vary,varx,process))
    
    cname = "%s_%s_%s_%s"%(varx,vary,cut,process)
    draw2Dplot(h,outfile,cut,cname)

def make2Dratioplot(varx,vary,cut,infile1,infile2,outfile):
    #""" Will divide histogram from infile1 by histogram from infile2 """
    hname = "h_%s_%s_%s" % (varx,vary,cut)
    h1 = infile1.Get(hname)
    h1.GetZaxis().SetTitle("Data/MC")
    h2 = infile2.Get(hname)
    h2.GetZaxis().SetTitle("Data/MC")
    # ratio = h1.Clone("new")
    nx = h1.GetNbinsX()
    xmin = h1.GetXaxis().GetXmin()
    xmax = h1.GetXaxis().GetXmax()
    ny = h1.GetNbinsY()
    ymin = h1.GetYaxis().GetXmin()
    ymax = h1.GetYaxis().GetXmax()
    
    # Do our own clone, as otherwise Z axis is messed up when drawing with 'colz'
    # Again some random ROOT evil...
    ratio = rt.TH2D("ratio","",nx,xmin,xmax,ny,ymin,ymax)
    for bin in range(h1.GetNbinsX()*h1.GetNbinsY()):
        ratio.SetBinContent(bin,h1.GetBinContent(bin))
    ratio.Divide(h2)
    ratio.GetYaxis().SetTitle(vary)
    ratio.GetXaxis().SetTitle(varx)
    ratio.GetZaxis().SetTitle("Data/MC")
    ratio.GetZaxis().SetTitleOffset(1.5)
    ratio.SetMaximum(3)
    ratio.SetTitle("Data/MC ratio plot in 2D ")
    
    cname = "%s_%s_%s_%s"%(varx,vary,cut,"DataMC")
    draw2Dplot(ratio,outfile,cut,cname)

def draw2Dplot(h,outfile,cut,cname):
    # Make the  canvas
    canvas = rt.TCanvas(cname,"")
    canvas.cd()
    canvas.SetRightMargin(0.17)
    h.Draw("colz")
    t = '#scale[1.1]{Selection ' + cut + '}'
    tex = rt.TLatex(0.52,0.82,t)
    tex.SetNDC();
    tex.Draw("same");

    canvas.cd()
    canvas.SaveAs(outputdir+"/"+cname+".pdf")
    outfile.cd()
    canvas.Write()
    canvas.Close() # close histogram, better for speed when making lots of plots

if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print "Run as: python %s  <outputdir> <bg histogram dir> (<data histogram dir>)" % (sys.argv[0])
        sys.exit()
    outputdir = sys.argv[1]
    histdir = sys.argv[2]
    datadir = histdir
    if len(sys.argv) < 4: 
        print "Will take %s as directory for data histogram as well" % (histdir)
    else:
        datadir = sys.argv[3]

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    outfile = TFile.Open(outputdir+"/2Dplots.root","RECREATE")

    # Integrated luminosity in fb-1s
    intlumi = 19.789 # ABCD
    analyzer = "rzrBoostMC"
    analyzer = "rzrBTevsel_plots"

    print "Will make plots for integrated luminosity of %.3f fb-1" % (intlumi)

    # For the 2D plots we only need the total background, 
    # not the separate components

    f_bg = TFile.Open(histdir+"/"+analyzer+"_bg.root")
    f_data = TFile.Open(datadir+"/"+analyzer+"_data.root")

    # set root styles
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetTextFont(42)

    cuts = ["Cleaning","HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
            "SIG","neleeq0","nmueq0","trackIso",
            "g1Mb0Ll","g1Mbg1W0Ll","0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT","g1Mbg1W1LlmT100", 
            "2mu","2mu0el","0Lb2mu0el","g1Mb2mu0el" 
            ]
    cuts = ["SIG" ]

    for cut in cuts:
        print "Making plot for", cut
        make2Dplot("MR","R2",cut,f_bg,outfile,"bg")        
        make2Dplot("MR","R2",cut,f_data,outfile,"data")
        make2Dratioplot("MR","R2",cut,f_data,f_bg,outfile)

    # Make 2D plots of mindeltaphi
    cut = "0Lbg1uW0Ll"
    make2Dplot("MR","minDeltaPhi",cut,f_bg,outfile,"bg")        
    make2Dplot("MR","minDeltaPhi",cut,f_data,outfile,"data")        
    make2Dratioplot("MR","minDeltaPhi",cut,f_data,f_bg,outfile)        
    make2Dplot("R2","minDeltaPhi",cut,f_bg,outfile,"bg")        
    make2Dplot("R2","minDeltaPhi",cut,f_data,outfile,"data")        
    make2Dratioplot("R2","minDeltaPhi",cut,f_data,f_bg,outfile)        



        
    outfile.Close()
