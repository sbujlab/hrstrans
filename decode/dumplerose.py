#!/usr/bin/python

import re, sys

#filename = 'Fwd_r12p5_STD.f'
#planes = [ 'q1ex', 'dent', 'dext', 'q3en', 'q3ex', 'fp' ]
#prefix = 'r12p5'

filename = 'prex_forward.f'
planes = [ 'sen', 'sex','col', 'q1ex', 'q2ex', 'den', 'q3en', 'q3ex', 'fp' ]
prefix = 'sp'




vars = [ 'x', 't', 'y', 'p' ]

for plane in planes:
    matrix = []
    for var in vars:
        #  start
        fdef = re.compile("function " + var + "_" + prefix + "_" + plane )
        matwork = re.compile(var + "_" +prefix+"_" + plane )

        minre = re.compile("xmin")
        maxre = re.compile("xmax")
        avre  = re.compile('avdat.*/\s+([-+]?(\d+(\.\d*)?|\.\d+)(E[-+]?\d+)?)/')
        coeffre = re.compile("coeff")
        retre = re.compile("return")
        endarr = re.compile("/")

        singlere = re.compile("coeff\(\s+(\d+)\)\s*\*x([1-5])1\s*$")

#        scalere  = re.compile('(([-+]?(\d+(\.\d*)?|\.\d+)(E[-+]?\d+)?)\s*,)+')
        
#        print "function " + var + "_sp_" + plane

        f = open(filename, 'r');

        line = f.readline()
        while not fdef.search(line):
            line = f.readline()
                #  Get avdat
        while not avre.search(line):
             line = f.readline()
        m = avre.search(line)
        avdat = float(m.group(1))

                #  Get Min, max data
        while not minre.search(line):
            line = f.readline()
        line = f.readline()
        m = line[6:].split(",")
        xmin = []
        for i in range(5):
            xmin.append(float(m[i]))
#        print xmin
        while not maxre.search(line):
            line = f.readline()
        line = f.readline()
        m = line[6:].split(",")
        xmax = []
        for i in range(5):
            xmax.append(float(m[i]))
#        print xmax

        # parse coefficients
        while not coeffre.search(line):
            line = f.readline()

        line = f.readline()
        coeffs = []
        notdone = True;
        while notdone:
            #print "Parsing ", line
            tokens = line[6:].split(",")
            for tok in tokens:
                if endarr.search(tok):
                    notdone = False
                else:
                   # print tok
                    try: 
                        coeffs.append(float(tok))
                    except ValueError:
                        continue
            line = f.readline()
#        print coeffs

        vmap = 5*[-1]

        while not matwork.search(line):
            line = f.readline()
        while not retre.search(line):
            line = f.readline()
            if singlere.search(line):
                m = singlere.search(line)
#                print line, m.groups(1), m.groups(2)
#                print line, m.groups()
#                print m.group(1), m.group(2)
                varia = int(m.group(2))
                index = int(m.group(1))
                vmap[varia-1] = index-1
                
        thesecoeff = []
        for i in range(5):
            if vmap[i] != -1:
                thesecoeff.append( coeffs[vmap[i]]*2/(xmax[i]-xmin[i]) )
            else :
                thesecoeff.append(0.0)


        matrix.append(thesecoeff)

        f.close()


    thesecoeff = []
    for i in range(5):
        if i == 4:
            thesecoeff.append(1.0)
        else:
            thesecoeff.append(0.0)
    matrix.append(thesecoeff)

    print plane
    for i in range(5):
        for j in range(5):
            sys.stdout.write("%7.3f " % (matrix[i][j]))
        print ""











