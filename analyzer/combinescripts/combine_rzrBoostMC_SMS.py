#!/usr/bin/env python

import os,sys
from string import *

anl = 'rzrBoostMC_SMS'

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
             

os.system('hadd -f '+dres+'/'+anl+'_T1ttcc_DM-10.root '+dres+'/'+anl+'*T1ttcc*DM-10*.root')
os.system('hadd -f '+dres+'/'+anl+'_T1ttcc_DM-25.root '+dres+'/'+anl+'*T1ttcc*DM-25*.root')
os.system('hadd -f '+dres+'/'+anl+'_T1ttcc_DM-80.root '+dres+'/'+anl+'*T1ttcc*DM-80*.root')

