#!/usr/bin/env python

import os,sys
from string import *

anl = 'rzrBTsmsskim'

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

# Do some cleaning:
#os.system('rm -f '+drestmp+'*')
#os.system('rm -f '+dflisttmp+'*')

# Get the samples list file:
if len(sys.argv) < 5:
    print 'Run as python '+sys.argv[0]+' <sample list filename> <mstop> <mLSP> <mGluino>'
    sys.exit()

samplesfile = sys.argv[1]

datasets = open(samplesfile).readlines()

mstop = sys.argv[2]
mLSP = sys.argv[3]
mG = sys.argv[4]

names = {}
names['exe'] = exe
names['drestmp'] = drestmp
names['anl'] = anl
names['mstop'] = mstop
names['mLSP'] = mLSP
names['mG'] = mG


# Total number of runs
nrunstot = 0
# Number of input files per run
nf = 10
# Integrated luminosity in pb-1s
intlumi = 19789
# Loop over the datasets
for d in datasets:
    # Read the sample name, and when relevant,
    #    - sample type: data or mc
    #    - cross section in pb-1
    #    - totalweight to be used for normalization (equal to
    #      total number of events when the events are unweighted)
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
    print d
    user = d[0]
    sample = d[1]
    samplew_ = sample.replace('/', '_')
    names['samplew_'] = samplew_
    # Read the filelist
    nmflist = dflist+'filelist'+samplew_+'.txt'
    lflist = open(nmflist).readlines()
    lenflist = len(lflist)
    # Number of runs
    nruns = (lenflist / nf) + 1
    # Delete any existing directories and files
    nmflisttmps = dflisttmp+'filelist'+samplew_+'*.txt'
    os.system('rm '+nmflisttmps)
    nmscrs = dscr+'scr_'+anl+samplew_+'*.py'
    os.system('rm '+nmscrs)
    print samplew_, lenflist, nruns
    # Partition the files in filelist into files for each run
    for f in range(0, lenflist):
        nrun = 1
        for i in range(0, lenflist, nf):
            nmflisttmp = dflisttmp+'filelist'+samplew_+str(nrun)+'.txt'
            if (f>=i and f<i+nf):
                open(nmflisttmp,'a').write(lflist[f])
            nrun = nrun+1
    # Make the scripts, rundirs, jobs, and submit
    for nrun in range(1, nruns+1):
        names['nrun'] = nrun
        nmflisttmp = dflisttmp+'filelist'+samplew_+str(nrun)+'.txt'
        names['nmflisttmp'] = nmflisttmp
        nmhistotmp = drestmp+anl+samplew_+str(nrun)+'.root'
        names['nmhistotmp'] = nmhistotmp
        # Make the run directory
        nmdrun = 'run_'+anl+samplew_+str(nrun)
        names['nmdrun'] = nmdrun
        # Make the job script
        nmscr = dscr+'scr_'+anl+samplew_+str(nrun)+'.py'
        names['nmscr'] = nmscr
        cscr = cscrtmp % names
        open(nmscr,'w').write(cscr)
        os.system('chmod +x '+nmscr)
        os.chdir('/tmp/nstrobbe/')
        #os.chdir("/tmp/ssekmen/")
        os.system('bsub -q 1nh '+nmscr)
        #print names
        nrunstot = nrunstot + 1

print '\n*** TOTAL NUMBER OF RUNS:', nrunstot, '\n' 
