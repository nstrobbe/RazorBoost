from ROOT import *

### for Zinv estimation
def get_Zll_acceptance(rootfilename,selection):
    f = TFile.Open(rootfilename)
    h = f.Get("counts")
    num = h.GetBinContent(h.GetXaxis().FindBin(selection+"2genmu"))
    den = h.GetBinContent(h.GetXaxis().FindBin(selection+"2genallmu"))
    acc = 0
    if den != 0:
        acc = num/den
    f.Close()
    return acc

def get_Zll_efficiency(rootfilename,selection):
    f = TFile.Open(rootfilename)
    h = f.Get("counts")
    num = h.GetBinContent(h.GetXaxis().FindBin(selection+"2mu"))
    den = h.GetBinContent(h.GetXaxis().FindBin(selection+"2genmu"))
    eff = 0
    if den != 0:
        eff = num/den
    f.Close()
    return eff

def get_Zll_scalefactor(rootfilename,selection):
    BR_nunu = 20.0
    BR_mumu = 3.366
    acc = get_Zll_acceptance(rootfilename,selection)
    eff = get_Zll_efficiency(rootfilename,selection)
    sf = BR_nunu / (BR_mumu * acc * eff)
    return sf

def get_Zvv_estimation(CR, selection_for_sf, fdataname, infodict):
    outfile = TFile.Open("BGestimate_Zinv_from_"+CR+"_sf_from_"+selection_for_sf+".root","RECREATE")
    
    # first get Zmumu estimation in the control region. For this get data and subtract rest
    hCR = "h_MR_R2_" + CR

    f_data = TFile.Open(fdataname)
    h_data = f_data.Get(hCR)

    otherBG = 0
    process_to_skip = ["ZJetsToNuNu","DYJetsToLNu","DYJetsToLNu_PtZ"]
    for name,info in infodict.iteritems():
        # skip the process we want to determine:
        if name in process_to_skip: continue
        filename = info[0]
        # open the file and get the histogram
        f = TFile.Open(filename)
        outfile.cd()
        h = f.Get(hCR)
        h.Sumw2()
        # if this is the first histogram we look at:
        if otherBG == 0:
            otherBG = h.Clone()
            otherBG.Sumw2()
        else:
            otherBG.Add(h)
        f.Close()

    h_est = h_data.Clone("h_MR_R2_Zinv_from_"+CR)
    h_est.Sumw2()
    h_est.Add(otherBG,-1)
    
    # scale factor containing branching ratios and acceptance + efficiency
    sf = get_Zll_scalefactor(rootfilename,selection_for_sf)
    h_est.Scale(sf)
    
    # TODO: add a factor for Wtag efficiency

    outfile.cd()
    h_est.Write()
    outfile.Close()


### For usual simple estimation
def make_info_dict(base,
                   QCD_CR="MC",TTJets_CR="MC",WJetsToLNu_CR="MC",Top_CR="MC",
                   ZJetsToNuNu_CR="MC",DYJetsToLL_CR="MC",VV_CR="MC",
                   Vhad_CR="MC",VVV_CR="MC",TTX_CR="MC",
                   QCD_binbybin=True,TTJets_binbybin=True,WJetsToLNu_binbybin=True,Top_binbybin=True,
                   ZJetsToNuNu_binbybin=True,DYJetsToLL_binbybin=True,VV_binbybin=True,
                   Vhad_binbybin=True,VVV_binbybin=True,TTX_binbybin=True,
                   QCD_include=True,TTJets_include=True,WJetsToLNu_include=True,Top_include=True,
                   ZJetsToNuNu_include=True,DYJetsToLL_include=True,VV_include=True,
                   Vhad_include=True,VVV_include=True,TTX_include=True,
                   ):
    # this will contain all the info pertaining to all background processes.
    info = dict()
    # for each process, form a list with following information:
    # [filename, fromMC, CR if required]
    info["QCD"] = [base+"QCD.root", QCD_CR, QCD_binbybin, QCD_include]
    info["TTJets"] = [base+"TTJets.root", TTJets_CR, TTJets_binbybin, TTJets_include]
    info["WJetsToLNu"] = [base+"WJetsToLNu.root", WJetsToLNu_CR, WJetsToLNu_binbybin, WJetsToLNu_include]
    info["Top"] = [base+"Top.root", Top_CR, Top_binbybin, Top_include]
    info["ZJetsToNuNu"] = [base+"ZJetsToNuNu.root", ZJetsToNuNu_CR, ZJetsToNuNu_binbybin, ZJetsToNuNu_include]
    info["DYJetsToLL"] = [base+"DYJetsToLL.root", DYJetsToLL_CR, DYJetsToLL_binbybin, DYJetsToLL_include]
    info["VV"] = [base+"VV.root", VV_CR, VV_binbybin, VV_include]
    info["VVV"] = [base+"VVV.root", VVV_CR, VVV_binbybin, VVV_include]
    info["TTX"] = [base+"TTX.root", TTX_CR, TTX_binbybin, TTX_include]
    info["Vhad"] = [base+"Vhad.root", Vhad_CR, Vhad_binbybin, Vhad_include]

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
        if not info[3]:
            print "skipping process", name
            continue
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
        if not info[3]:
            print "Skipping background", name, "Should be obtained in a different way"
            continue
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

### Use only one (total) bg histogram instead of summing over different components
def simple_estimate(sigregion,cregion,fdataname,fbgname):
    outfile = TFile.Open("BGestimate_simple_"+sigregion+"_from_"+cregion+".root","RECREATE")

    fdata = TFile.Open(fdataname)
    fbg = TFile.Open(fbgname)
    
    # estimate_sig = observed_CR * (N_sig_all/N_CR_all)_MC

    obs = fdata.Get("h_MR_R2_"+cregion)
    obs.Sumw2()
    N_sig = fbg.Get("h_MR_R2_"+sigregion)
    N_sig.Sumw2()
    N_cr = fbg.Get("h_MR_R2_"+cregion)
    N_cr.Sumw2()

    pred = obs.Clone()
    #pred.Sumw2()
    pred.Multiply(N_sig)
    pred.Divide(N_cr)

    outfile.cd()
    pred.Write()
    outfile.Close()
    
if __name__ == "__main__":

    print "Will perform the background estimate"

    # directory containing all the raw histograms
    inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140331_noISR_btag_TopPt/summary/" 
    inputdirbase = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140331_noISR_btag_TopPt/" 

    # analyzer name
    analyzer = "rzrBoostMC"
    
    # make the dictionary with process information, one per region
    SIG_info = make_info_dict(inputdir+analyzer+"_",
                              QCD_CR="0Lbg1uW0Ll_mdPhiHat4",
                              TTJets_CR="g1Mbg1W1LlmT100",
                              WJetsToLNu_CR="0Lbg1Y1LlmT"
                              )

    SIGlike_info = make_info_dict(inputdir+analyzer+"_",
                                  QCD_CR="0Lbg1uW0Ll_mdPhiHat4",
                                  TTJets_CR="g1Mbg1W1LlmT100",
                                  WJetsToLNu_CR="0Lbg1Y1LlmT"
                                  )
    
    # signal region with mindeltaphihat > 4
    SIG_dphig4_info = make_info_dict(inputdir+analyzer+"_",
                                    QCD_CR="0Lbg1uW0Ll_mdPhiHat4",QCD_binbybin=True,
                                    TTJets_CR="g1Mbg1W1LlmT100",
                                    WJetsToLNu_CR="0Lbg1Y1LlmT"
                                    )
    SIG_dphig4_info2 = make_info_dict(inputdir+analyzer+"_",
                                      QCD_CR="0Lbg1uW0Ll_mdPhiHat4",QCD_binbybin=True,
                                      TTJets_CR="g1Mbg1W1LlmT100_mdPhiHatg4",
                                      WJetsToLNu_CR="0Lbg1Y1LlmT_mdPhiHatg4"
                                      )
    
    # signal region with mindeltaphihat < 4
    SIG_dphi4_info = make_info_dict(inputdir+analyzer+"_",
                                    QCD_CR="0Lbg1uW0Ll_mdPhiHat4",QCD_binbybin=True,
                                    TTJets_CR="g1Mbg1W1LlmT100",
                                    WJetsToLNu_CR="0Lbg1Y1LlmT"
                                    )
    SIG_dphi4_info2 = make_info_dict(inputdir+analyzer+"_",
                                     QCD_CR="0Lbg1uW0Ll_mdPhiHat4",QCD_binbybin=True,
                                     TTJets_CR="g1Mbg1W1LlmT100_mdPhiHatg4",
                                     WJetsToLNu_CR="0Lbg1Y1LlmT_mdPhiHatg4"
                                     )
    
    #QCD_info = make_info_dict(inputdir+analyzer+"_",
    #                          ZJetsToNuNu_CR="0Lbg1Y2mu0el",ZJetsToNuNu_binbybin=False)
    QCD_info = make_info_dict(inputdir+analyzer+"_")

    WJ_info = make_info_dict(inputdir+analyzer+"_"
                             #,QCD_CR="0Lbg1uW0Ll_mdPhiHat4"
                             )

    TTJ_info = make_info_dict(inputdir+analyzer+"_",
                              QCD_CR="0Lbg1uW0Ll_mdPhiHat4"
                              )
    #TTJ_info2 = make_info_dict(inputdir+analyzer+"_")             

    #Zll_info = make_info_dict(inputdir+analyzer+"_")

                                
    # put all this in one big dictionary:
    # for signal region estimation
    info = {"g1Mbg1W0Ll":SIG_info,
            "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
            "g1Mbg1W1LlmT100":TTJ_info,
            "0Lbg1Y1LlmT": WJ_info
            }

    # for signal like region estimation
    info2 = {"g1Mb0Wg1uW0Ll":SIGlike_info,
             "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
             "g1Mbg1W1LlmT100":TTJ_info,
             "0Lbg1Y1LlmT": WJ_info
             }

    info_dphig4 = {"g1Mbg1W0Ll_mdPhiHatg4":SIG_dphig4_info,
                   "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
                   "g1Mbg1W1LlmT100":TTJ_info,
                   "0Lbg1Y1LlmT": WJ_info
                   }
    info2_dphig4 = {"g1Mbg1W0Ll_mdPhiHatg4":SIG_dphig4_info,
                    "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
                    "g1Mbg1W1LlmT100_mdPhiHatg4":TTJ_info,
                    "0Lbg1Y1LlmT_mdPhiHatg4": WJ_info
                   }

    info_dphi4 = {"g1Mbg1W0Ll_mdPhiHat4":SIG_dphi4_info,
                  "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
                  "g1Mbg1W1LlmT100":TTJ_info,
                  "0Lbg1Y1LlmT": WJ_info
                  }
    info2_dphi4 = {"g1Mbg1W0Ll_mdPhiHat4":SIG_dphi4_info2,
                   "0Lbg1uW0Ll_mdPhiHat4":QCD_info,
                   "g1Mbg1W1LlmT100_mdPhiHatg4":TTJ_info,
                   "0Lbg1Y1LlmT_mdPhiHatg4": WJ_info
                   }

    #doBGestimate(region,infodict,fdata,extra_info)

    #doBGestimate("g1Mbg1W0Ll",info,inputdir+analyzer+"_data.root","_QCDWJTTJ_Feb5") # BG estimate for signal region, using bin-by-bin ratio
    #doBGestimate("g1Mbg1W0Ll_mdPhiHatg4",info_dphig4,inputdir+analyzer+"_data.root","_QCDWJTTJ_Feb5") # BG estimate for signal region, using bin-by-bin ratio

    #doBGestimate("g1Mbg1W0Ll",info,inputdir+analyzer+"_data.root","global") # BG estimate for signal region, using global MC ratio

    #regions_Zll_est = ["NoCuts","presel",""]
    #for region in regions_Zll_est:
    #    print "Acceptance", region, ",", get_Zll_acceptance(inputdir+"rzrBoostMC_DYJetsToLL.root",region)
    #    print "Efficiency", region, ",", get_Zll_efficiency(inputdir+"rzrBoostMC_DYJetsToLL.root",region)
    #    print "Scale factor", region, ",", get_Zll_scalefactor(inputdir+"rzrBoostMC_DYJetsToLL.root",region)

    #Zinv_info = make_info_dict(inputdir+analyzer+"_")
    #get_Zvv_estimation("g1Mbg1Y2mu0el", "", inputdir+analyzer+"_data.root", Zinv_info)



    ###############################################################################
    ##  Perform some closure tests doing the full estimation on several regions  ##
    ###############################################################################

    # signal-like region
    doBGestimate("g1Mb0Wg1uW0Ll",info2,inputdir+analyzer+"_data.root","_QCDWJTTJ") 

    # signal region mdphihat<4
    doBGestimate("g1Mbg1W0Ll_mdPhiHat4",info_dphi4,inputdir+analyzer+"_data.root","_QCDWJTTJ") 
 
    # signal region mdphihat<4, using T and W with mdphihat>4
    doBGestimate("g1Mbg1W0Ll_mdPhiHat4",info2_dphi4,inputdir+analyzer+"_data.root","_QCDWJTTJ_2") 
 



    ############################################################
    ##  Perform some closure tests using the simple estimate  ##
    ############################################################

    # TTJ and WJ composition
    simple_estimate("g1Mbg1W1LlmT100","0Lbg1Y1LlmT",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")
    # TTJ and WJ composition, using W and T with mdphihat>4
    simple_estimate("g1Mbg1W1LlmT100_mdPhiHatg4","0Lbg1Y1LlmT_mdPhiHatg4",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

    # Btag modeling
    simple_estimate("1Mbg1W1LlmT100","0Lbg1Y1LlmT",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")
    simple_estimate("g2Mbg1W1LlmT100","1Mbg1W1LlmT100",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

    # Modeling of minDeltaPhiHat for QCD
    simple_estimate("0Lbg1uW0Ll_mdPhiHatg4","0Lbg1uW0Ll_mdPhiHat4",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

    # Modeling of minDeltaPhiHat for TTJ
    simple_estimate("g1Mbg1W1LlmT100_mdPhiHatg4","g1Mbg1W1LlmT100_mdPhiHat4",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

    # Modeling of minDeltaPhiHat for WJ
    simple_estimate("0Lbg1Y1LlmT_mdPhiHatg4","0Lbg1Y1LlmT_mdPhiHat4",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

    # lepton modeling
    simple_estimate("g1Mbg1W0Ll_mdPhiHat4","g1Mbg1W1LlmT100_mdPhiHat4",inputdir+analyzer+"_data.root",inputdirbase+analyzer+"_bg.root")

