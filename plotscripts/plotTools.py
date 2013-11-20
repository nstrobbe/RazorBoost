##################################################################
### This file will contain usefull plotting routines that can  ###
### be easily included in other scripts.                       ###
##################################################################
import sys, os
import ROOT as rt

def SetBoostStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetTextFont(42)
    rt.TH1.SetDefaultSumw2()

def ConstructHDict(h, name="name", color=rt.kBlue, title="title", appear_in_ratio="Yes"
                   , linestyle=1, linewidth=2, markerstyle=7, markersize=1, fillstyle=0
                   , xtitle="", ytitle=""):
    hdict = {}
    hdict["name"] = name                       # what will appear in the legend
    hdict["histogram"] = h                     # the actual histogram
    hdict["color"] = color                     # the color for the histogram
    hdict["title"] = title                     # title of the histogram
    hdict["appear in ratio"] = appear_in_ratio # can be "Yes", "No", "Ref"
    hdict["linestyle"] = linestyle             # line style
    hdict["linewidth"] = linewidth             # line width
    hdict["markerstyle"] = markerstyle         # marker style
    hdict["markersize"] = markersize           # marker size
    hdict["fillstyle"] = fillstyle             # fill style
    hdict["xtitle"] = xtitle                   # x axis title
    hdict["ytitle"] = ytitle                   # y axis title
    return hdict

def ConstructLDict(xmin,xmax,ymin,ymax,title):
    ldict = {}
    legdict["xmin"] = xmin
    legdict["ymin"] = ymin
    legdict["xmax"] = xmax
    legdict["ymax"] = ymax
    legdict["title"] = title
    return ldict

def Plot2D():
    pass

def Plot1DWithRatio(hdictlist=0,outputdir="plots",outfile=0,legdict=0,cname="canvas"
                    ,ratiotitle="ratio",logscale=False,scale=True):
    # First do some checks on the input
    if hdictlist == 0:
        print "You did not give a list of histogram dictionaries, will quit"
        sys.exit(0)
    if outfile == 0:
        print "You did not pass me a root file to store the plots. I will only produce pdf files."
    if not os.path.isdir(outputdir):
        print "Output directory doesn't exist yet"
        print "Will create directory %s" % (outputdir)
        os.makedirs(outputdir)

    # placeholder for draw objects
    rootEvil = []

    # Place to store histogram that will be used as reference in ratio plot
    h_ref = None

    # Make the canvas, which will have two pads
    canvas = rt.TCanvas(cname,"")
    pad1 = rt.TPad("pad1","",0,0.25,1,1)
    pad1.SetBottomMargin(0)
    if logscale:
        pad1.SetLogy(1)
    pad1.Draw()
    pad1.cd()

    # Make the legend
    legend = rt.TLegend(0.5,0.5,0.87,0.87,"")
    if legdict != 0:
        legend = rt.TLegend(legdict["xmin"],legdict["ymin"],legdict["xmax"],legdict["ymax"],legdict["title"])
    legend.SetFillColor(0)
    legend.SetBorderSize(0)

    # Get histograms from a list of dictionaries, and plot them 
    print "Getting all histograms"

    first = 0 # to keep track whether this is the first histogram we're plotting
    for hdict in hdictlist:
        h = hdict["histogram"] # Get the histogram
        h.Sumw2() # need to put this otherwise errors in ratio plot are wrong
        if scale:
            sf = h.Integral(0,h.GetNbinsX()+1)
            h.Scale(1./sf)
                
        if hdict["appear in ratio"] == "Ref":
            print "Found ref histo"
            h_ref=h.Clone("h_ref")

        if scale:
            h.GetYaxis().SetTitle("A.U.")
        else:
            h.GetYaxis().SetTitle("Events")

        h.GetYaxis().SetTitleSize(0.055)
        h.GetYaxis().SetTitleOffset(0.8)
        h.GetYaxis().SetLabelSize(0.05)
        h.SetLineColor(hdict["color"])
        h.SetLineWidth(hdict["linewidth"])
        h.SetTitle(hdict["title"])
        legend.AddEntry(h,hdict["name"],"l")
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
    for hdict in hdictlist:
        if hdict["appear in ratio"]!= "Yes": continue

        h = hdict["histogram"]
        ratio = h.Clone()
        h_ref_ratio = h_ref.Clone()
        ratio.Divide(h_ref_ratio)

        ratio.SetMarkerColor(hdict["color"])
        ratio.SetMarkerStyle(hdict["markerstyle"])
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
        ratio.GetXaxis().SetTitle(hdict["xtitle"])
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
    if outfile != 0:
        outfile.cd()
        canvas.Write()
    canvas.Close()

if __name__ == "__main__":
    print "This is the tools script, which isn't supposed to be run on it's own"
