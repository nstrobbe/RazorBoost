#!/usr/bin/env python

import os,sys
from string import *

anl = 'rzrBoostMC'

dwork = '/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/'
dflist = dwork+'filelists/'
dflisttmp = dwork+'fileliststmp/'
druns = dwork+'runs/'
dscr = dwork+'scripts/'
djob = dwork+'jobs/'
dtemplates = dwork+'templates/'
#nmscrtmp = dtemplates+'scr_'+anl+'.py'
nmscrtmp = dtemplates+'scr_'+anl+'.py'
cscrtmp = open(nmscrtmp).read()
nmjobtmp = dtemplates+'job.job'
cjobtmp = open(nmjobtmp).read()
dres = dwork+'results/'
drestmp = dwork+'resultstmp/'
exe = dwork+anl


# Get the samples list file:
if len(sys.argv) < 2:
    print 'Run as python '+sys.argv[0]+' <sample list filename>'
    sys.exit()

samplesfile = sys.argv[1]

datasets = open(samplesfile).readlines()

# Loop over the datasets
for d in datasets:
    d = strip(d)
    d = split(d)
    d = split(d[0], '%')
    print d
    user = d[0]
    sample = d[1]
    samplew_ = sample.replace('/', '_')
    nmhisto = dres+anl+samplew_+'histo.root'
    nmhistostmp = drestmp+anl+samplew_+'*'
    os.system('hadd -f '+nmhisto+' '+nmhistostmp)
             
os.system('hadd -f '+dres+'/'+anl+'_data.root '+dres+'/'+anl+'_Jet*.root '+dres+'/'+anl+'_HT*.root')
#os.system('hadd -f '+dres+'/'+anl+'_dataSingleMu.root '+dres+'/'+anl+'_SingleMu*.root ')

os.system('hadd -f '+dres+'/'+anl+'_QCD.root '+dres+'/'+anl+'_QCD*.root')
os.system('hadd -f '+dres+'/'+anl+'_VV.root '+dres+'/'+anl+'_WW*.root '+dres+'/'+anl+'_ZZ*.root '+dres+'/'+anl+'_WZ*.root')
os.system('hadd -f '+dres+'/'+anl+'_TTJets.root '+dres+'/'+anl+'_TTJets*.root ')
os.system('hadd -f '+dres+'/'+anl+'_Top.root '+dres+'/'+anl+'_Tbar*.root '+dres+'/'+anl+'_T_*.root')
os.system('hadd -f '+dres+'/'+anl+'_WJetsToLNu.root '+dres+'/'+anl+'_WJetsToLNu*.root')
os.system('hadd -f '+dres+'/'+anl+'_ZJetsToNuNu.root '+dres+'/'+anl+'_ZJetsToNuNu*.root')
os.system('hadd -f '+dres+'/'+anl+'_DYToCC.root '+dres+'/'+anl+'_DYToCC*.root')
os.system('hadd -f '+dres+'/'+anl+'_DYToBB.root '+dres+'/'+anl+'_DYToBB*.root')
os.system('hadd -f '+dres+'/'+anl+'_Wbb.root '+dres+'/'+anl+'_Wbb*.root')
os.system('hadd -f '+dres+'/'+anl+'_Vhad.root '+dres+'/'+anl+'_DYToCC.root '+dres+'/'+anl+'_DYToBB.root '+dres+'/'+anl+'_Wbb.root')
os.system('hadd -f '+dres+'/'+anl+'_DYJetsToLL.root '+dres+'/'+anl+'_DYJetsToLL_M*.root')
os.system('hadd -f '+dres+'/'+anl+'_DYJetsToLL_PtZ.root '+dres+'/'+anl+'_DYJetsToLL_PtZ*.root')
os.system('hadd -f '+dres+'/'+anl+'_VVV.root '+dres+'/'+anl+'_WZZ*aMCatNLO*.root '+dres+'/'+anl+'_WWW*aMCatNLO*.root '+dres+'/'+anl+'_WWZ*.root '  +dres+'/'+anl+'_WWG*.root '+dres+'/'+anl+'_ZZZ*.root')
os.system('hadd -f '+dres+'/'+anl+'_TTX.root '+dres+'/'+anl+'_ttbarZ*.root '+dres+'/'+anl+'_TTbarW*.root '+dres+'/'+anl+'_TTG*.root '+dres+'/'+anl+'_TTWW*.root')

os.system('hadd -f '+dres+'/'+anl+'_bg_oldDYsamples.root '+dres+'/'+anl+'_QCD.root '+dres+'/'+anl+'_VV.root '+dres+'/'+anl+'_TTJets.root '+dres+'/'+anl+'_WJetsToLNu.root '+dres+'/'+anl+'_Top.root '+dres+'/'+anl+'_ZJetsToNuNu.root '+dres+'/'+anl+'_DYJetsToLL.root '+dres+'/'+anl+'_VVV.root '+dres+'/'+anl+'_TTX.root '+dres+'/'+anl+'_Wbb.root '+dres+'/'+anl+'_DYToCC.root '+dres+'/'+anl+'_DYToBB.root')
os.system('hadd -f '+dres+'/'+anl+'_bg.root '+dres+'/'+anl+'_QCD.root '+dres+'/'+anl+'_VV.root '+dres+'/'+anl+'_TTJets.root '+dres+'/'+anl+'_WJetsToLNu.root '+dres+'/'+anl+'_Top.root '+dres+'/'+anl+'_ZJetsToNuNu.root '+dres+'/'+anl+'_DYJetsToLL_PtZ.root '+dres+'/'+anl+'_VVV.root '+dres+'/'+anl+'_TTX.root '+dres+'/'+anl+'_Wbb.root '+dres+'/'+anl+'_DYToCC.root '+dres+'/'+anl+'_DYToBB.root')

