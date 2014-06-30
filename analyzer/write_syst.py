import random as rd
import sys
import datetime

def generate_line():
    # Here we generate one line in the systematics file
    # it consists of 15 gaussians, and then a uniform number for the pdfs

    sysline = ""
    for i in range(15):
        number = rd.gauss(0,1)
        sysline = sysline + str(number) + " "

    # randrange returns integers
    pdfnumber = rd.randrange(1,301)
    sysline = sysline + str(pdfnumber)
    return sysline

def generate_lineOnlyPDF():
    # Here we generate one line in the systematics file
    # it consists of 15 zeros, and then a uniform number for the pdfs

    sysline = ""
    for i in range(15):
        number = 0.
        sysline = sysline + str(number) + " "

    # randrange returns integers
    pdfnumber = rd.randrange(1,301)
    sysline = sysline + str(pdfnumber)
    return sysline

def generate_lineNoPDF():
    # Here we generate one line in the systematics file
    # it consists of 15 gaussians

    sysline = ""
    for i in range(15):
        number = rd.gauss(0,1)
        sysline = sysline + str(number) + " "

    return sysline


if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print "Usage: python %s nsyst" % (sys.argv[0])
        exit()

    nsyst = int(sys.argv[1])
    today = datetime.date.today()
    fname = "systematics_%s.txt" % (str(today).replace("-",""))

    systfile = open(fname,'w')

    for i in range(nsyst):
        l = generate_line()
        systfile.write(l+"\n")

    systfile.close()
