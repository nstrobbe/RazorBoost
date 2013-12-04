from ROOT import *

def make_info_dict(base,
                   QCD_fromMC=True,TTJets_fromMC=True,WJetsToLNu_fromMC=True,Top_fromMC=True,
                   ZJetsToNuNu_fromMC=True,DYJetsToLL_fromMC=True,VV_fromMC=True,
                   Vhad_fromMC=True,VVV_fromMC=True,TTX_fromMC=True,
                   QCD_CR="",TTJets_CR="",WJetsToLNu_CR="",Top_CR="",
                   ZJetsToNuNu_CR="",DYJetsToLL_CR="",VV_CR="",
                   Vhad_CR="",VVV_CR="",TTX_CR="",
                   ):
    # this will contain all the info pertaining to all background processes.
    info = dict()
    # for each process, form a list with following information:
    # [filename, fromMC, CR if required]
    info["QCD"] = [base+"QCD.root", QCD_fromMC, QCD_CR]
    info["TTJets"] = [base+"TTJets.root", TTJets_fromMC, TTJets_CR] 
    info["WJetsToLNu"] = [base+"WJetsToLNu.root", WJetsToLNu_fromMC, WJetsToLNu_CR]
    info["Top"] = [base+"Top.root", Top_fromMC, Top_CR]
    info["ZJetsToNuNu"] = [base+"ZJetsToNuNu.root", ZJetsToNuNu_fromMC, ZJetsToNuNu_CR]
    info["DYJetsToLL"] = [base+"DYJetsToLL.root", DYJetsToLL_fromMC, DYJetsToLL_CR]
    info["VV"] = [base+"VV.root", VV_fromMC, VV_CR]
    info["VVV"] = [base+"VVV.root", VVV_fromMC, VVV_CR]
    info["TTX"] = [base+"TTX.root", TTX_fromMC, TTX_CR]
    info["Vhad"] = [base+"Vhad.root", Vhad_fromMC, Vhad_CR]

    return info

def calcBGfromCR(process, SIG, CR, fdataname, infodict, outfile, bin_by_bin=True):
    hCRname = "h_MR_R2_" + CR
    hSIGname = "h_MR_R2_" + SIG

    # get the data histogram
    fdata = TFile.Open(fdataname)
    outfile.cd()
    hdata = fdata.Get(hCRname)
    hdata.Sumw2()

    # get the histogram with other MC
    # for now assume all but the process considered have to be taken from MC
    otherBG = 0
    for name,info in infodict.iteritems():
        # skip the process we want to determine:
        if name == process: continue
        # open the file and get the histogram
        f = TFile.Open(info[0])
        outfile.cd()
        h = f.Get(hCRname)
        h.Sumw2()
        # if this is the first histogram we look at:
        if otherBG == 0:
            otherBG = h.Clone()
            otherBG.Sumw2()
        else:
            otherBG.Add(h)
        f.Close()

    # get the MC ratio
    fprocess = TFile.Open(infodict[process][0])
    outfile.cd()
    h_SIG = fprocess.Get(hSIGname)
    h_SIG.Sumw2()
    h_CR = fprocess.Get(hCRname)
    h_CR.Sumw2()
    h_ratio = h_SIG.Clone()
    h_ratio.Sumw2()
    h_ratio.Divide(h_CR)
    ratio = h_SIG.Integral()/h_CR.Integral() if h_CR.Integral()>0 else 1
    
    # now do the actual estimate: (data - others) * ratio
    BGestimate = hdata.Clone("h_"+process+"_in_"+SIG+"_from_"+CR)
    BGestimate.Sumw2()
    BGestimate.Add(otherBG,-1)
    if bin_by_bin:
        BGestimate.Multiply(h_ratio)
    else:
        BGestimate.Scale(ratio)
        
    outfile.cd()
    BGestimate.Write()
    fdata.Close()
    fprocess.Close()

    return BGestimate

    
def doBGestimate(region,infodict,fdata,extra_info="",bin_by_bin=True):
    # We need to sum the contributions of all the BG components
    # Some will be taken straight from MC, others from control region
    # all this information is provided in the "info" dictionary

    # file to store the estimate
    outfile = TFile.Open("BGestimate_"+region+"_"+extra_info+".root","RECREATE")
    BGestimate = 0
    for name,info in infodict.iteritems():
        print "Determine background", name
        # if we want to take it from MC:
        if info[1]:
            # open the file
            f = TFile.Open(info[0])
            outfile.cd()
            # get the histogram
            h = f.Get("h_MR_R2_"+region)
            h.Sumw2()
            # if this is the first histogram we look at:
            if BGestimate == 0:
                BGestimate = h.Clone()
                BGestimate.Sumw2()
            else:
                htemp = h.Clone()
                htemp.Sumw2()
                BGestimate.Add(htemp)
            f.Close()
        else:
            # we need to take it from a control region
            h = calcBGfromCR(name, region, info[2], fdata, infodict, outfile, bin_by_bin)
            print h
            BGestimate.Add(h)
            
    outfile.cd()
    BGestimate.Write()
    outfile.Close()

    
if __name__ == "__main__":

    print "Will perform the background estimate"

    # directory containing all the raw histograms
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20131125/summary/" 

    # analyzer name
    analyzer = "rzrBoostMC"
    
    # make the dictionary with process information
    SIG_info = make_info_dict(inputdir+analyzer+"_",
                              QCD_fromMC=False,QCD_CR="0Lbg1uW0Ll_mdPhi0p3",
                              TTJets_fromMC=False,TTJets_CR="g1Mbg1W1LlmT100")

    #doBGestimate(region,infodict,fdata,extra_info,bin_by_bin)
    doBGestimate("g1Mbg1W0Ll",SIG_info,inputdir+analyzer+"_data.root","",True) # BG estimate for signal region, using bin-by-bin ratio

    doBGestimate("g1Mbg1W0Ll",SIG_info,inputdir+analyzer+"_data.root","global",False) # BG estimate for signal region, using global MC ratio

    # need to add functionality to do for example: estimate TTJets from TTJets CR ( taking QCD from QCD CR and rest from MC)

