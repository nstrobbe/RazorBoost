#!/usr/bin/env python

import os,sys
from string import *

anl = 'rzrBoostsys'

dwork = '/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/RazorBoost/analyzer/'
#dwork = '/afs/cern.ch/work/s/ssekmen/RazorBoost/analyzer/'
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
dres = dwork+'results_sys/'
drestmp = dwork+'resultstmp_sys/'
exe = dwork+anl

# Get the samples list file:
if len(sys.argv) < 2:
    print 'Run as python '+sys.argv[0]+' <sample list filename>'
    sys.exit()

samplesfile = sys.argv[1]

datasets = open(samplesfile).readlines()

# Find the systematic directories in resultstmp and make the same under results:
sysdirs = []
systmpdirs = os.popen('ls '+drestmp).readlines()
for s in systmpdirs:
    s = split(strip(s), '/')[-1]
    print s
    if 'sys' in s:
        sysdirs.append(s)

print sysdirs

# Loop over the systematic directories:
for s in sysdirs:
    dresnsys = dres+s+'/'
    if os.path.exists(dresnsys):
        os.system('rm -rf '+dresnsys)
    os.mkdir(dresnsys)
    drestmpnsys = drestmp+s+'/'
    # Loop over the datasets
    for d in datasets:
        d = strip(d)
        d = split(d)
        d = split(d[0], '%')
        print d
        user = d[0]
        sample = d[1]
        samplew_ = sample.replace('/', '_')
        nmhisto = dresnsys+anl+samplew_+'histo.root'
        nmhistostmp = drestmpnsys+anl+samplew_+'*'
        print ''
        print 'hadd -f '+nmhisto+' '+nmhistostmp
        os.system('hadd -f '+nmhisto+' '+nmhistostmp)

    os.system('hadd -f '+dresnsys+'/'+anl+'_QCD.root '+dresnsys+'/'+anl+'_QCD*.root')
    os.system('hadd -f '+dresnsys+'/'+anl+'_TTJets.root '+dresnsys+'/'+anl+'_TTJets*.root '+dresnsys+'/'+anl+'_Tbar*.root '+dresnsys+'/'+anl+'_T_*.root')
    os.system('hadd -f '+dresnsys+'/'+anl+'_WJetsToLNu.root '+dresnsys+'/'+anl+'_WJetsToLNu*.root')

    os.system('hadd -f '+dresnsys+'/'+anl+'_BGoth.root '+dresnsys+'/'+anl+'_WW_*.root '+dresnsys+'/'+anl+'_ZZ_*.root '+dresnsys+'/'+anl+'_WZ_*.root '+dresnsys+'/'+anl+'_ZJetsToNuNu*.root '+dresnsys+'/'+anl+'_DYJetsToLL_HT*.root '+dresnsys+'/'+anl+'_WZZ*aMCatNLO*.root '+dresnsys+'/'+anl+'_WWW*aMCatNLO*.root '+dresnsys+'/'+anl+'_WWZ*.root '  +dresnsys+'/'+anl+'_WWG*.root '+dresnsys+'/'+anl+'_ZZZ*.root '+dresnsys+'/'+anl+'_ttbarZ*.root '+dresnsys+'/'+anl+'_TTbarW*.root '+dresnsys+'/'+anl+'_TTG*.root '+dresnsys+'/'+anl+'_TTWW*.root '+dresnsys+'/'+anl+'_Wbb*.root '+dresnsys+'/'+anl+'_DYToCC*.root '+dresnsys+'/'+anl+'_DYToBB*.root')



