#!/usr/bin/python

# Author: Frank Gelhausen

import sys, getopt, subprocess
import matplotlib.pyplot as plt
sys.path.append("/usr/local/lib/python2.7/site-packages/")
import RNA
import re
import math
intaRNAPath = "../../../IntaRNA_basepairprobs/bin/IntaRNA"

# VR1 straight mRNA: AGCUUGGUACAAGCGCAUCUUCUACUUCAACUUCUAGAA
# VsiRNA: UUCGCGUAGAAGAUGAAGUUG

def main(argv):
    seqT = ""
    seqQ = ""
    method = ""
    title = ""
    cTinit = 1e-5
    cQinit = 1e-5
    paramFile = None
    plotFile = None
    temperature = 37

    try:
        opts, args = getopt.getopt(argv, "hq:t:", ["qConc=", "tConc=", "paramFile=", "method=", "plot="])
    except getopt.GetoptError:
        printHelp()
    for opt, arg in opts:
        if opt == "-h":
            printHelp()
        elif opt == "-t":
            seqT = arg
        elif opt == "-q":
            seqQ = arg
        elif opt == "--tConc":
            cTinit = float(arg)
        elif opt == "--qConc":
            cQinit = float(arg)
        elif opt == "--paramFile":
            paramFile = arg
        elif opt == "--plot":
            plotFile = arg
        elif opt == "--method":
            method = arg

    if seqT == "" or seqQ == "":
        printHelp()

    if method != "I" and method != "C":
        printHelp()

    # parse paramFile
    if paramFile != None:
        f = open(paramFile, "r")
        for line in f:
            split = line.split("=")
            param = map(str.strip, split)
            if param[0] == "temperature":
                temperature = float(param[1])
            elif param[0] == "energy":
                if param[1] != "V":
                    print("--energy != V not supported")
                    sys.exit(2)
            elif param[0] == "energyVrna":
                print("--energyVrna not supported")
                sys.exit(2)
            elif param[0] == "tAccConstr":
                print("--tAccConstr not supported")
                sys.exit(2)
            elif param[0] == "qAccConstr":
                print("--qAccConstr not supported")
                sys.exit(2)

    # TODO: required because of a bug in ViennaRNA (calling deprecated method)
    # (see: https://github.com/ViennaRNA/ViennaRNA/issues/69)
    RNA.co_pf_fold("C", None)

    # if IntaRNA ...
    if method == "I":
        print("Using IntaRNA ...")
        title = "IntaRNA"
        fcTT = getInteractionFeature(seqT, seqT, temperature, "EallTotal")
        fcQQ = getInteractionFeature(seqQ, seqQ, temperature, "EallTotal")
        fcTQ = getInteractionFeature(seqT, seqQ, temperature, "EallTotal")
        Tc = getInteractionFeature(seqT, seqQ, temperature, "Eall1")
        Qc = getInteractionFeature(seqT, seqQ, temperature, "Eall2")
    else:
        print("Using rnaCofold ...")
        title = "RNAcofold"
        (fcTQ, fcTT, fcQQ, Tc, Qc) = rnaCofold("sequences.fa")

    print("free energies: ", fcTQ, fcTT, fcQQ, Tc, Qc)

    # binding energy:
    print("binding energy: ", fcTQ - Tc - Qc)

    # plot concentration
    if plotFile != None:
        xList = []
        cTList = []
        cTTList = []
        cTQList = []
        cQQList = []

        for i in range(1, 50):
            xList.append(i)

            (cTQ, cTT, cQQ, cT, cQ) = RNA.get_concentrations(fcTQ, fcTT, fcQQ, Tc, Qc, cTinit, i)
            cTList.append(cT)
            cTTList.append(cTT)
            cTQList.append(cTQ)
            cQQList.append(cQQ)

        plt.xscale('log')
        plt.plot(xList, cTQList, color='black', label="A.si", linewidth=3)
        plt.plot(xList, cTTList, color='red', label="A.A")
        plt.plot(xList, cTList, color='blue', label="A", linewidth=3)
        plt.plot(xList, cQQList, color='pink', label="si")

        plt.xlim(1, 50)
        plt.ylim(0, 20)
        plt.xlabel("total siRNA concentration [nmol/L]", fontsize=16)
        plt.ylabel("concentration [nmol/L]", fontsize=16)
        plt.legend(loc='upper left')
        plt.suptitle(title, fontsize=20)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)

        plt.text(5, -7.5, "Binding energies: ΔF(A) = %.2f kcal/mol"%(fcTQ - Tc - Qc), ha='center', fontsize=16)
        plt.subplots_adjust(bottom=0.3)

        plt.savefig(plotFile)

    # output concentration
    (cTQ, cTT, cQQ, cT, cQ) = RNA.get_concentrations(fcTQ, fcTT, fcQQ, Tc, Qc, cTinit, cQinit)
    print("concentrations: ", cTQ, cTT, cQQ, cT, cQ)

# returns given feature for RNAFold
def rnaFold(inputFile):
    res = subprocess.check_output(["RNAFold", "-P", "./rna_turner1999.par", "--noPS", inputFile]).decode('ascii', 'ignore').splitlines()
    split = res[2].split(" ", 1)
    print("mfe: ", split[1])
    return split[0]

# returns given feature for RNACofold
def rnaCofold(inputFile):
    res = subprocess.check_output(["RNAcofold", "-P", "./rna_turner1999.par", "--noPS", "--all_pf=2", inputFile]).decode('ascii', 'ignore').splitlines()
    fe = res[7].split("\t")
    return (float(fe[0]), float(fe[1]), float(fe[2]), float(fe[3]), float(fe[4]))

# returns given feature for IntaRNA interaction
def getInteractionFeature(seqT, seqQ, temperature, feature):
    res = subprocess.check_output([intaRNAPath, "-q", seqQ, "-t", seqT, "--mode=M", "--energyVRNA=Turner99", "--outMode=E", "--noseed", "--temperature=%d"%(temperature)]).decode('ascii', 'ignore')
    for line in res.splitlines():
        if feature in line:
            return float(line.split(" ")[1])
    return 0

# prints usage information
def printHelp():
    print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>] [--method <I|C>] [--paramFile <FileName>] [--paramFile <FileName.png>]")
    sys.exit(2)

if __name__ == "__main__":
    main(sys.argv[1:])
