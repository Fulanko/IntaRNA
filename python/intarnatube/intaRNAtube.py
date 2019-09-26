#!/usr/bin/python

# Author: Frank Gelhausen

import sys, getopt, subprocess
sys.path.append("/usr/local/lib/python2.7/site-packages/")
import RNA
intaRNAPath = "../../../IntaRNA_basepairprobs/bin/IntaRNA"

def main(argv):
    seqT = ""
    seqQ = ""
    cTinit = 1e-5
    cQinit = 1e-5
    paramFile = None
    temperature = 37

    try:
        opts, args = getopt.getopt(argv, "hq:t:", ["qConc=", "tConc=", "paramFile="])
    except getopt.GetoptError:
        print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>] [--paramFile <FileName>]")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>] [--paramFile <FileName>]")
            sys.exit()
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

    if seqT == "" or seqQ == "":
        print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>] [--paramFile <FileName>]")
        sys.exit(2)

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

    (ss, pfT) = RNA.fold_compound(seqT).pf()
    (ss, pfQ) = RNA.fold_compound(seqQ).pf()

    # TODO: required because of a bug in ViennaRNA (calling deprecated method)
    # (see: https://github.com/ViennaRNA/ViennaRNA/issues/69)
    RNA.co_pf_fold("C", None)

    pfTT = getInteractionFeature(seqT, seqT, temperature, "Eall")
    pfQQ = getInteractionFeature(seqQ, seqQ, temperature, "Eall")
    pfTQ = getInteractionFeature(seqT, seqQ, temperature, "Eall")

    (cTQ, cTT, cQQ, cT, cQ) = RNA.get_concentrations(pfTQ, pfTT, pfQQ, pfT, pfQ, cTinit, cQinit)
    cTQ /= (cTinit + cQinit)
    cTT /= (cTinit + cQinit)
    cQQ /= (cTinit + cQinit)
    cT /= (cTinit + cQinit)
    cQ /= (cTinit + cQinit)

    print(cTQ, cTT, cQQ, cT, cQ)

# returns given feature for IntaRNA interaction
def getInteractionFeature(seqT, seqQ, temperature, feature):
    res = subprocess.check_output([intaRNAPath, "-q", seqQ, "-t", seqT, "--mode=M", "--outMode=E", "--noseed", "--temperature=%d"%(temperature)])
    for line in res.splitlines():
        if feature in line:
            return float(line.split(" ")[1])
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])
