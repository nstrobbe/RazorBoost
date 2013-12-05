from ROOT import *

def make_info_dict(base,
                   QCD_CR="MC",TTJets_CR="MC",WJetsToLNu_CR="MC",Top_CR="MC",
                   ZJetsToNuNu_CR="MC",DYJetsToLL_CR="MC",VV_CR="MC",
                   Vhad_CR="MC",VVV_CR="MC",TTX_CR="MC",
                   QCD_binbybin=True,TTJets_binbybin=True,WJetsToLNu_binbybin=True,Top_binbybin=True,
                   ZJetsToNuNu_binbybin=True,DYJetsToLL_binbybin=True,VV_binbybin=True,
                   Vhad_binbybin=True,VVV_binbybin=True,TTX_binbybin=True
                   ):
    # this will contain all the info pertaining to all background processes.
    info = dict()
    # for each process, form a list with following information:
    # [filename, fromMC, CR if required]
    info["QCD"] = [base+"QCD.root", QCD_CR, QCD_binbybin]
    info["TTJets"] = [base+"TTJets.root", TTJets_CR, TTJets_binbybin] 
    info["WJetsToLNu"] = [base+"WJetsToLNu.root", WJetsToLNu_CR, WJetsToLNu_binbybin]
    info["Top"] = [base+"Top.root", Top_CR, Top_binbybin]
    info["ZJetsToNuNu"] = [base+"ZJetsToNuNu.root", ZJetsToNuNu_CR, ZJetsToNuNu_binbybin]
    info["DYJetsToLL"] = [base+"DYJetsToLL.root", DYJetsToLL_CR, DYJetsToLL_binbybin]
    info["VV"] = [base+"VV.root", VV_CR, VV_binbybin]
    info["VVV"] = [base+"VVV.root", VVV_CR, VVV_binbybin]
    info["TTX"] = [base+"TTX.root", TTX_CR, TTX_binbybin]
    info["Vhad"] = [base+"Vhad.root", Vhad_CR, Vhad_binbybin]

    return info

def calcBGfromCR(process, SIG, CR, fdataname, infodict, outfile):
    hCRname = "h_MR_R2_" + CR
    hSIGname = "h_MR_R2_" + SIG

    # get the data histogram
    fdata = TFile.Open(fdataname)
    outfile.cd()
    hdata = fdata.Get(hCRname)
    hdata.Sumw2()

    # get the histogram with other MC
    # Some processes are to be taken from MC, others from a different control region
    otherBG = 0
    info_CR = infodict[CR]
    for name,info in info_CR.iteritems():
        # skip the process we want to determine:
        if name == process: continue
        filename = info[0]
        poss_CR = info[1]               
        if poss_CR == "MC":
            # take process straight from MC
            # open the file and get the histogram
            f = TFile.Open(filename)
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
        else:
            # take process from other control region
            h = calcBGfromCR(name,CR,poss_CR,fdataname, infodict, outfile)
            if otherBG == 0:
                otherBG = h.Clone()
                otherBG.Sumw2()
            else:
                otherBG.Add(h)

    # get the MC ratio
    fprocess = TFile.Open(infodict[CR][process][0])
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
    bin_by_bin = infodict[SIG][process][2]
    if bin_by_bin:
        BGestimate.Multiply(h_ratio)
    else:
        BGestimate.Scale(ratio)
        
    outfile.cd()
    BGestimate.Write()
    fdata.Close()
    fprocess.Close()

    return BGestimate

    
def doBGestimate(region,infodict,fdata,extra_info=""):
    # We need to sum the contributions of all the BG components
    # Some will be taken straight from MC, others from control region
    # all this information is provided in the "info" dictionary

    # file to store the estimate
    outfile = TFile.Open("BGestimate_"+region+extra_info+".root","RECREATE")
    BGestimate = 0
    region_info = infodict[region]
    for name,info in region_info.iteritems():
        print "Determine background", name
        # unpack info
        filename = info[0]
        poss_CR = info[1]
        # if we want to take it from MC:
        if poss_CR == "MC":
            # open the file
            f = TFile.Open(filename)
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
            h = calcBGfromCR(name, region, poss_CR, fdata, infodict, outfile)
            if BGestimate == 0:
                BGestimate = h.Clone()
                BGestimate.Sumw2()
            else:
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
    
    # make the dictionary with process information, one per region
    SIG_info = make_info_dict(inputdir+analyzer+"_",
                              QCD_CR="0Lbg1uW0Ll_mdPhi0p3",
                              TTJets_CR="g1Mbg1W1LlmT100")

    QCD_info = make_info_dict(inputdir+analyzer+"_",
                              ZJetsToNuNu_CR="0Lbg1Y2mu0el",ZJetsToNuNu_binbybin=False)
    #QCD_info2 = make_info_dict(inputdir+analyzer+"_")

    TTJ_info = make_info_dict(inputdir+analyzer+"_",
                              QCD_CR="0Lbg1uW0Ll_mdPhi0p3")
    #TTJ_info2 = make_info_dict(inputdir+analyzer+"_")             

    Zll_info = make_info_dict(inputdir+analyzer+"_")

    # put all this in one big dictionary:
    info = {"g1Mbg1W0Ll":SIG_info,
            "0Lbg1uW0Ll_mdPhi0p3":QCD_info,
            "g1Mbg1W1LlmT100":TTJ_info,
            "0Lbg1Y2mu0el": Zll_info}


    #doBGestimate(region,infodict,fdata,extra_info)
    doBGestimate("g1Mbg1W0Ll",info,inputdir+analyzer+"_data.root","_test") # BG estimate for signal region, using bin-by-bin ratio

    #doBGestimate("g1Mbg1W0Ll",info,inputdir+analyzer+"_data.root","global") # BG estimate for signal region, using global MC ratio

