import glob

if __name__ == "__main__":

    indir = "resultstmp"
    files = glob.glob(indir+"/*T2tt*.dat")
    count = 0
    with open("T2tt.dat",'w') as outfile:
        for f in files:
            with open(f) as infile:
                for line in infile:
                    if "MR" in line and count == 1: continue
                    outfile.write(line)
            count = 1
    print "Concatenated all files into one big file"
