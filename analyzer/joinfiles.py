import glob

if __name__ == "__main__":

    indir = "resultstmp"
    outfilename = "BGnonQCD.dat"
    files_extension = ".dat"

    things_to_veto = ["QCD"]
    things_to_keep = []

    files = glob.glob(indir+"/*"+files_extension)
    count = 0
    with open(outfilename,'w') as outfile:
        for f in files:
            stop = False
            if len(things_to_veto) > 0:
                for veto in things_to_veto:
                    if veto in f:
                        stop = True
                        break
            if stop: continue

            stop2 = False
            if len(things_to_keep) > 0:
                anything_to_keep = False
                for keep in things_to_keep:
                    if keep in f:
                        anything_to_keep = True
                        break
                if not anything_to_keep:
                    stop2 = True
            if stop2: continue
            
            with open(f) as infile:
                for line in infile:
                    if "MR" in line and count == 1: continue
                    outfile.write(line)
            count = 1
    print "Concatenated all files into one big file"
