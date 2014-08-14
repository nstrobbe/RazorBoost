#!/usr/bin/env python

import os,sys
from string import *

# Get the samples list file:
if len(sys.argv) < 3:
    print 'Run as python '+sys.argv[0]+' <analyzername> ><sample list filename>'
    sys.exit()

anl = sys.argv[1]
samplesfile = sys.argv[2]


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
dres = dwork+'results_sys/'
drestmp = dwork+'resultstmp_sys/'
exe = dwork+anl


datasets = open(samplesfile).readlines()


names = {}
names['exe'] = exe
names['drestmp'] = drestmp
names['anl'] = anl

# Find the systematic directories in resultstmp and make the same under results:
sysdirs = []
systmpdirs = os.popen('ls '+drestmp).readlines()
for s in systmpdirs:
    s = split(strip(s), '/')[-1]
    print s
    if 'sys' in s:
        #if float(s.replace("sys","")) >= 448:
        sysdirs.append(s)

print sysdirs

# Total number of runs
nrunstot = 0
# Number of input files per run
nf = 8
# Integrated luminosity in pb-1s
intlumi = 19712
# Loop over the datasets

print 'MISSING FILES: \n'

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
    #print d
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

    # Loop over the systematic directories:
    for s in sysdirs:
        drestmpnsys = drestmp+s+'/'
        dscrnsys = dscr+s+'/'
        # Partition the files in filelist into files for each run
        for nrun in range(1, nruns+1):
            nmflisttmp = dflisttmp+'filelist'+samplew_+str(nrun)+'.txt'
            nmhistotmp = drestmpnsys+anl+samplew_+str(nrun)+'.root'
            nmscr = dscrnsys+'scr_'+anl+samplew_+str(nrun)+'.py'
            if not os.path.exists(nmhistotmp):
                if os.path.exists(nmflisttmp):
                    print nmhistotmp
                    os.chdir("/tmp/nstrobbe/")
                    os.system('bsub -q 1nh '+nmscr)
                    nrunstot = nrunstot + 1


print '\n'+str(nrunstot)+' JOBS RESUBMITTED\n'
