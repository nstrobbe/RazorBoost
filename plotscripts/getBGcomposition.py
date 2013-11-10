import sys
from string import *
from ROOT import TFile
from math import sqrt 

def printCutPercentage(cut,info):
    ttj = 0
    wj = 0
    qcd =0
    diboson= 0
    top = 0
    zjets = 0
    zll = 0
    triboson = 0
    ttX = 0
    for sample,t in info.iteritems():
        c = t[0][cut]
        w = t[1]
        effxs = c #*w*intlumi
        if "_TTJets_" in sample:
            ttj = ttj + effxs
        if "_WJetsToLNu_" in sample:
            wj = wj + effxs
        if "_QCD_" in sample:
            qcd = qcd + effxs
        if "_WW_" in sample:
            diboson = diboson + effxs
        if "_WZ_" in sample:
            diboson = diboson + effxs
        if "_ZZ_" in sample:
            diboson = diboson + effxs
        if "_T_" in sample:
            top = top + effxs
        if "_Tbar_" in sample:
            top = top + effxs
        if "_ZJetsToNuNu_" in sample:
            zjets = zjets + effxs
        if "_DYJets" in sample: 
            zll = zll + effxs
        if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
            triboson = triboson + effxs
        if "_TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
            ttX = ttX + effxs
    total = ttj + wj + qcd + diboson + top + zjets + zll + triboson + ttX
    print "Cut %s \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f" % (cut,qcd*100./total,ttj*100./total,wj*100./total,diboson*100./total,top*100./total,zjets*100./total,zll*100./total,triboson*100./total,ttX*100./total) 

def printCut(cut,info,intlumi):
    ttj = 0
    wj = 0
    qcd = 0
    diboson = 0
    top = 0
    zjets = 0
    zll = 0
    triboson = 0
    ttX = 0
    #sig = 0
    for sample,t in info.iteritems():
        c = t[0][cut]
        w = t[1]
        effxs = c #*w*intlumi
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
        if "_DYJets" in sample: 
            zll = zll + effxs
        if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
            triboson = triboson + effxs
        if "_TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
            ttX = ttX + effxs
        # if "T1ttcc" in sample:
        #    sig = sig + effxs
    
    print "Cut %s \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d" % (cut,qcd,ttj,wj,diboson,top,zjets,zll,triboson,ttX) 
    
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


def getrow(cut,info,intlumi):
    ttj = 0
    wj = 0
    qcd = 0
    diboson = 0
    top = 0
    zjets = 0
    zll = 0
    triboson = 0
    ttX = 0
    #sig = 0
    for sample,t in info.iteritems():
        c = t[0][cut]
        w = t[1]
        effxs = c #*w*intlumi
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
        if "_DYJets" in sample: 
            zll = zll + effxs
        if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
            triboson = triboson + effxs
        if "_TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
            ttX = ttX + effxs
        # if "T1ttcc" in sample:
        #    sig = sig + effxs
    
    totalbg = ttj+wj+qcd+diboson+top+zjets+zll+triboson+ttX
    #sign = sig/sqrt(totalbg)
    #row = "%s & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.2g \\\\ \n" % (cut,qcd,ttj,wj,diboson,top,zjets,totalbg,sig,sign) 
    row = "%s & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g \\\\ \n" % (cut,qcd,ttj,wj,diboson,top,zjets,zll,triboson,ttX,totalbg) 
    print row
    return row.replace("%","\%")


def getrowpercentage(cut,info,intlumi):
    ttj = 0
    wj = 0
    qcd = 0
    diboson = 0
    top = 0
    zjets = 0
    zll = 0
    triboson = 0
    ttX = 0
    for sample,t in info.iteritems():
        c = t[0][cut]
        w = t[1]
        effxs = c #*w*intlumi
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
        if "_DYJets" in sample: 
            zll = zll + effxs
        if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
            triboson = triboson + effxs
        if "_TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
            ttX = ttX + effxs
    
    totalbg = ttj+wj+qcd+diboson+top+zjets+zll+triboson+ttX
    row = " & {0:.1%} & {1:.1%} & {2:.1%} & {3:.1%} & {4:.1%} & {5:.1%} & {6:.1%} & {7:.1%} & {8:.1%} \\\\ \n".format(qcd/totalbg,ttj/totalbg,wj/totalbg,diboson/totalbg,top/totalbg,zjets/totalbg,zll/totalbg,triboson/totalbg,ttX/totalbg)
    print row
    return row.replace("%","\%")

def makeCutflowTable(info,intlumi):
    table = open("cutflow.tex",'w')
    # get the table header stuff
    header = """\\begin{sidewaystable}[p]
\\centering
\\caption{Cutflow table, event counts are normalized to $19.789\\textrm{fb}^{-1}$. }
\\begin{tabular}{| l || c | c | c | c | c | c | c | c | c || c |}
\\hline
Cut & QCD & TTJets & WJets & VV & T & Zinv & Zll & VVV & TTX & Total \\\\ \\hline
"""
    table.write(header)

    # put the body of the table
    table.write(getrow("NoCuts",info,intlumi))
    table.write(getrow("Cleaning",info,intlumi))
    table.write(getrow("vertexg0",info,intlumi))
    table.write(getrow("njetge3",info,intlumi))
    table.write(getrow("HLT",info,intlumi))
    table.write(getrow("jet1ptg200",info,intlumi))
    table.write(getrow("SIG",info,intlumi))
    table.write(getrowpercentage("SIG",info,intlumi))

    table.write("\\hline \\hline \n")

    table.write(getrow("neleeq0",info,intlumi))
    table.write(getrow("nmueq0",info,intlumi))
    table.write(getrow("trackIso",info,intlumi))

    table.write("\\hline \n")

    table.write(getrow("g1Mb0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("g1Mb0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("g1Mbg1W0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("g1Mbg1W0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("0Lb0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("0Lb0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("0Lbg1uW0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("0Lbg1uW0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("0Lbg1W0Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("0Lbg1W0Ll",info,intlumi).replace("_","\\_"))

    table.write("\\hline \n")

    table.write(getrow("1Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("1Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("g1Mb1Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("g1Mb1Ll",info,intlumi).replace("_","\\_"))
    table.write(getrow("g1Mbg1W1Ll",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("g1Mbg1W1Ll",info,intlumi).replace("_","\\_"))

    table.write("\\hline \n")

    table.write(getrow("2mu",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("2mu",info,intlumi).replace("_","\\_"))
    table.write(getrow("2mu0el",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("2mu0el",info,intlumi).replace("_","\\_"))
    table.write(getrow("0Lb2mu0el",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("0Lb2mu0el",info,intlumi).replace("_","\\_"))
    table.write(getrow("g1Mb2mu0el",info,intlumi).replace("_","\\_"))
    table.write(getrowpercentage("g1Mb2mu0el",info,intlumi).replace("_","\\_"))

    # put the table closing stuff
    footer = """\\hline
\\end{tabular}
\\label{tab:cutflow}
\\end{sidewaystable}
"""
    table.write(footer)
    
    table.close()

def allcuts():
    cuts = ["NoCuts","Cleaning","HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
            "SIG","neleeq0","nmueq0","trackIso",
            "g1Mb0Ll","g1Mbg1W0Ll","0Lb0Ll","0Lbg1uW0Ll","0Lbg1W0Ll",
            "1Ll","g1Mb1Ll","g1Mbg1W1Ll", 
            "2mu","2mu0el","0Lb2mu0el","g1Mb2mu0el" 
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
        
    # Integrated luminosity in pb-1s
    intlumi = 19789.
    #intlumi = 5126.

    # get signal info
    #fs = TFile.Open("selge3j/summary/rzrBTevsel_plots_T1ttcc_325_300.root")
    #histos = fs.Get("counts")
    #countss = {}
    #for bin in range(histos.GetNbinsX()):
    #    countss[histos.GetXaxis().GetBinLabel(bin+1)] = histos.GetBinContent(bin+1)
    # print "xsec, counts", xsect, counts["NoCuts"]
    #weights = 0.0243547/countss["NoCuts"]
    #info["T1ttcc"] = (countss,weights)
    
    
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
        f = TFile.Open("results_BoostMC_1/rzrBoostMC"+samplew_+"histo.root")
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
    print "Cut \t \t qcd \t ttj \t wj \t diboson \t single top \t Zinv \t Zll \t VVV \t TTX \n"

    print "Counts for %d pb-1 of data" % (intlumi,)
    for cut in cuts:
        printCut(cut,info,intlumi)
        
    print "\n"
    
    print "Background composition in percent"
    for cut in cuts:
        printCutPercentage(cut,info)

    print "\n"
    
    print "TTbar composition"
    print "Cut \t allhad \t semilep \t dilep "
    for cut in cuts:
        getTTjetsComposition(cut,counts_allhad,counts_semilep,counts_dilep)


    print "Make cutflow table"
    makeCutflowTable(info,intlumi)
