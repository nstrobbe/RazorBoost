# check the file size of all root files
# if too small, something went wrong, so delete file

import os,sys

dir_to_check = "results_tmp"
if len(sys.argv) > 1: 
    dir_to_check = sys.argv[1]

ndeleted = 0
systmpdirs = os.popen('ls '+ dir_to_check).readlines()
for s in systmpdirs:
    fullpath = dir_to_check + "/" + s.strip()
    print fullpath
    # now check all files
    fnames = os.popen('ls '+ fullpath).readlines()
    for fname in fnames:
        fullname = fullpath+"/"+fname.strip()
        stats = os.stat(fullname)
        if stats.st_size < 1000: 
            print "Deleting ", fullname, " with size ", stats.st_size
            os.remove(fullname)
            ndeleted = ndeleted + 1

print "DELETED %s FILES" % (ndeleted)
