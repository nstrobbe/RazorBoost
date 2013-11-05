#!/usr/bin/env python

import os,sys
from string import *

# First setup root
os.environ['SCRAM_ARCH'] = 'slc5_amd64_gcc462'
os.environ['ROOTSYS'] = '/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms6/'
path = os.environ['ROOTSYS']+'/bin:'+os.environ['PATH']
os.environ['PATH'] = path
os.environ['LD_LIBRARY_PATH'] = '/afs/cern.ch/work/n/nstrobbe/CMSSW_5_2_5/lib/slc5_amd64_gcc462:/afs/cern.ch/work/n/nstrobbe/CMSSW_5_2_5/external/slc5_amd64_gcc462/lib:/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_2_5/lib/slc5_amd64_gcc462:/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_2_5/external/slc5_amd64_gcc462/lib:/afs/cern.ch/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib64:/afs/cern.ch/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/d-cache/dcap/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/d-cache/dcap/lib64:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/glite/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/glite/lib64:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/globus/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/lcg/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/lcg/lib64:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/external/usr/lib64:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/external/usr/lib:/afs/cern.ch/project/gd/LCG-share/3.2.11-1/classads/lib64/:/afs/cern.ch/cms/ccs/wm/scripts/Crab/CRAB_2_8_3/external/sqlite//3.3.5/lib'

dwork = '/afs/cern.ch/work/n/nstrobbe/RazorBoost/rzrBoostMCanalyzer/'
exe = '%(exe)s'
anl = '%(anl)s'
xsect = '%(xsect)s'
totweight = '%(totweight)s'
lumi = '%(lumi)s'
nmflisttmp = '%(nmflisttmp)s'
nmhistotmp = '%(nmhistotmp)s'
histo = anl+'_histograms.root'

# Args:
args = nmflisttmp+' '+xsect+' '+totweight+' '+lumi+' '+histo

# Run
os.system(exe+' '+args)

# Copy whatever results in the temporary results directory
os.system('cp '+histo+' '+nmhistotmp)

