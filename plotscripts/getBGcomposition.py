import sys
from string import *
from ROOT import TFile
from math import sqrt 

def printCut(cut,info,signame="T1ttcc",option="Counts"):
    ttj = 0
    wj = 0
    qcd = 0
    diboson = 0
    top = 0
    zjets = 0
    zll = 0
    triboson = 0
    ttX = 0
    Wbb = 0
    DYhad = 0
    sig = 0
    for sample,t in info.iteritems():
        #print sample, t, cut
        c = t[0][cut]
        w = t[1]
        effxs = c 
        if "TTJets" in sample:
            ttj = ttj + effxs
        if "WJets" in sample:
            wj = wj + effxs
        if "QCD" in sample:
            qcd = qcd + effxs
        if "WW" in sample:
            diboson = diboson + effxs
        if "WZ" in sample:
            diboson = diboson + effxs
        if "ZZ" in sample:
            diboson = diboson + effxs
        if "_T_" in sample:
            top = top + effxs
        if "_Tbar_" in sample:
            top = top + effxs
        if "_ZJetsToNuNu_" in sample:
            zjets = zjets + effxs
        if "_DYJetsToLL" in sample: 
            zll = zll + effxs
        if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
            triboson = triboson + effxs
        if "_TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
            ttX = ttX + effxs
        if "_DYToCC_" in sample or "_DYToBB_" in sample:
            DYhad = DYhad + effxs
        if "_Wbb_" in sample:
            Wbb = Wbb + effxs
        if signame in sample:
            sig = sig + effxs
    total = ttj + wj + qcd + diboson + top + zjets + zll + triboson + ttX + Wbb + DYhad
    if option == "Counts":
        print "Cut %s \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d" % (cut,qcd,ttj,wj,diboson,top,zjets,zll,triboson,ttX,Wbb,DYhad,sig)
        row = "%s & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %d & %.3g \\\\ \n" % (cut,qcd,ttj,wj,diboson,top,zjets,zll,triboson,ttX,Wbb,DYhad,total,sig) 
        return row.replace("%","\%")
    
    elif option == "Percentage":
        if total == 0:
            total = 1
        print "Cut %s \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f" % (cut,qcd*100./total,ttj*100./total,wj*100./total,diboson*100./total,top*100./total,zjets*100./total,zll*100./total,triboson*100./total,ttX*100./total,Wbb*100./total,DYhad*100./total)
        row = " & {0:.1%} & {1:.1%} & {2:.1%} & {3:.1%} & {4:.1%} & {5:.1%} & {6:.1%} & {7:.1%} & {8:.1%} & {9:.1%} & {10:.1%} \\\\ \n".format(qcd/total,ttj/total,wj/total,diboson/total,top/total,zjets/total,zll/total,triboson/total,ttX/total,Wbb/total,DYhad/total )
        return row.replace("%","\%")

    
def getTTjetsComposition(cut,counts_allhad,counts_semilep,counts_dilep):
    allhad = counts_allhad[cut]
    semilep = counts_semilep[cut]
    dilep = counts_dilep[cut]
    total = allhad + semilep + dilep
    total = float(total)
    if total <= 0:
        print "Cut %s \t / \t / \t / " % cut
    else: 
        print "Cut %s \t %.2f \t %.2f \t %.2f" % (cut,allhad/total,semilep/total,dilep/total)


def makeCutflowTable(info,intlumi,cuts,signame):
    table = open("cutflow.tex",'w')
    # get the table header stuff
    header = """\\begin{sidewaystable}[p]
\\centering
\\caption{Cutflow table, event counts are normalized to $%(lumi)s\\textrm{fb}^{-1}$. }
\\begin{tabular}{| l || c | c | c | c | c | c | c | c | c | c | c || c || c |}
\\hline
Cut & QCD & TTJets & WJets & VV & T & Zinv & Zll & VVV & TTX & Wbb & DYhad & Total & Signal \\\\ \\hline
""" % {"lumi":intlumi}
    table.write(header)

    # put the body of the table
    for cut in cuts:
        table.write( printCut(cut,info,signame,"Counts") )
        table.write( printCut(cut,info,signame,"Percentage") )
        table.write("\\hline\n")
        
    # put the table closing stuff
    footer = """\\end{tabular}
\\label{tab:cutflow}
\\end{sidewaystable}
"""
    table.write(footer)
    
    table.close()

def allcuts():
    cuts = ["NoCuts","Cleaning","HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
            "SIG","neleeq0","nmueq0","trackIso",
            "g1Mb0Ll","g1Mbg1W0Ll","1Mbg1W0Ll","g2Mbg1W0Ll","g1Mbg1W0Ll_mdPhiHat4","g1Mbg1W0Ll_mdPhiHatg4","g1Mb0Wg1uW0Ll",
            "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll","g1Mbg1W1LlmT100","1Mbg1W1LlmT100","g2Mbg1W1LlmT100","g1Mbg1W1LlmT",
            "0Lb1Ll","0Lbg1Y1Ll","0Lbg1Y1LlmT100","0Lbg1Y1LlmT",
            "2munoZmass","2mu","2mu0el","0Lb2mu0el","0Lbg1Y2mu0el","g1Mb2mu0el","g1Mbg1Y2mu0el",
            "2elnoZmass","2el","2el0mu","0Lb2el0mu","0Lbg1Y2el0mu","g1Mb2el0mu","g1Mbg1Y2el0mu",
            "2lnoZmass","2l","2l0ol","0Lb2l0ol","0Lbg1Y2l0ol","g1Mb2l0ol","g1Mbg1Y2l0ol",
            ]
    return cuts

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Run as: python %s <samples file> " % (sys.argv[0])
        sys.exit()
    samplesfile = sys.argv[1]
    
    datasets = open(samplesfile).readlines()

    names = {}
    info = {}
    counts_allhad = {}
    counts_semilep = {}
    counts_dilep = {}
        
    # Integrated luminosity in fb-1s
    intlumi = 19.789

    # input information
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140204"
    analyzer = "rzrBoostMC"
    signame = "T1ttcc_325_300"
    sigxs = 0.0243547
    
    # get signal info
    fs = TFile.Open(inputdir + "/summary/" + analyzer + "_" + signame + ".root")
    histos = fs.Get("counts")
    countss = {}
    for bin in range(histos.GetNbinsX()):
        countss[histos.GetXaxis().GetBinLabel(bin+1)] = histos.GetBinContent(bin+1)
    weights = sigxs/countss["NoCuts"]
    info[signame] = (countss,weights)
    
    
    for d in datasets:
        xsect = -1
        totweight = -1
        lumi = -1
        d = strip(d)
        d = split(d)
        if len(d) >= 2:
            sampletype = d[1]
            if sampletype == "mc":
                lumi = intlumi
        if len(d) >= 3:
            xsect = d[2]
        if len(d) == 4:
            totweight = d[3]
        names['xsect'] = xsect
        names['totweight'] = totweight 
        names['lumi'] = lumi
        d = split(d[0], '%')
        # print d
        user = d[0]
        sample = d[1]
        samplew_ = sample.replace('/', '_')
        names['samplew_'] = samplew_

        # print samplew_
        # get the count:
        f = TFile.Open(inputdir+"/"+analyzer+samplew_+"histo.root")
        histo = f.Get("counts")
        counts = {}
        for bin in range(histo.GetNbinsX()):
            counts[histo.GetXaxis().GetBinLabel(bin+1)] = histo.GetBinContent(bin+1)
        #print "sample, xsec, counts", samplew_, xsect, counts["NoCuts"]
        if counts["NoCuts"] == 0: continue
        weight = float(xsect)/counts["NoCuts"]
        info[samplew_] = (counts,weight)

        # now also get TTbar composition
        if "_TTJets_" in samplew_:
            h_allhad = f.Get("counts_TTallhad")
            h_semilep = f.Get("counts_TTsemilep")
            h_dilep = f.Get("counts_TTdilep")
            for bin in range(h_allhad.GetNbinsX()):
                counts_allhad[h_allhad.GetXaxis().GetBinLabel(bin+1)] = h_allhad.GetBinContent(bin+1)
            for bin in range(h_semilep.GetNbinsX()):
                counts_semilep[h_semilep.GetXaxis().GetBinLabel(bin+1)] = h_semilep.GetBinContent(bin+1)
            for bin in range(h_dilep.GetNbinsX()):
                counts_dilep[h_dilep.GetXaxis().GetBinLabel(bin+1)] = h_dilep.GetBinContent(bin+1)

        f.Close()


    cuts = allcuts()
    print "Make cutflow table"

    print "Cut \t \t qcd \t ttj \t wj \t diboson \t single top \t Zinv \t Zll \t VVV \t TTX \t Wbb \t DYhad \t signal \n"
    print "Counts and percentage for %d fb-1 of data" % (intlumi)
    makeCutflowTable(info,intlumi,cuts,signame)

    print "\n"
    
    print "TTbar composition"
    print "Cut \t allhad \t semilep \t dilep "
    for cut in cuts:
        getTTjetsComposition(cut,counts_allhad,counts_semilep,counts_dilep)


