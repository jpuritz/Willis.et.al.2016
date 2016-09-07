##################################################
##  Script to simulate RADseq-like data         ##
##  author:  Deren Eaton                        ##
##  contact: deren.eaton@yale.edu               ##
##  date:    5/8/14                             ##
##  version: 1.03                               ##
##################################################

##################################################
##  change log
##  1.03: - paired data file names have _R1_ & _R2_
##  1.02: - dropout is only affected by arg 2 not data type
##  1.01: - barcodes differ by 2 bp by default
##################################################


## load modules                                        
import egglib
import sys
import numpy
import gzip
import os

## load args
indel = float(sys.argv[1])
sitemut = int(sys.argv[2])
nrepets = int(sys.argv[3])   ## must be multiple of 100 or 1K
Ninds = int(sys.argv[4])     ## per tiptax
frags = str(sys.argv[5])
datatype = str(sys.argv[6])
outhandle = str(sys.argv[7])

## check args
if datatype in ['rad','gbs','pairgbs','ddrad','pairddrad']:
    print '\tsimulating '+datatype+" data"
else:
    print 'datatype not recognized'
    sys.exit()

if nrepets < 1000:
    divis = 100
else:
    divis = 1000

if nrepets % divis:
    if nrepets < 1000:
        print "nloci not multiple of 100 " ; sys.exit()
    else:
        print "nloci not multiple of 1K " ; sys.exit()


print "\tsimulating %i loci at an average of 20X coverage across 2 tip taxa with %i samples per taxon" % (nrepets, Ninds)
print "\tindels arise at frequency of %f per mutation" % indel
if sitemut:
    print "\tmutations in restriction site cause locus dropout"
else:
    print "\tmutations in restriction site =",bool(sitemut)
print "\tsequencing error rate = 0.001"

## fixed args    
copies  = 20
locuslength = 1000
N = 70000
u = (7.0e-9)*locuslength
theta = 4*N*u
print '\ttheta=4Nu=',theta/locuslength

if int(frags.split(",")[0]) < 200:
    print "\tmin fragment length allows read overlaps/adapter sequences "


""" cuts at Pst1 (and) at EcoR1, attaches barcode to Pst1 """
CUT1 = "AATT"    ## C|TGCAG
CUT2 = "CCG"     ## G|AATT


tiptax = list("ABCDX")
tiptax = [i+j for i,j in zip(list("11114"),list("ABCDX"))]
outgroup = "4X"

"random number seeds used in test example"
#numpy.random.seed(998877)

###Edits by JBP
#Change tree structure to two pairs of populations 

"simulate species with divergence times"
A=B=C=D=X=Ninds   ## number of diploid individuals to sample

###Edits by JBP
#Change tree structure to two pairs of populations 
"clade 1 = ((A,B),(C,D))"
"tree = ((A,B),(C,D))"
newick = "(((A:0.2,B:0.2):0.2,(C:0.2,D:0.2):0.2):0.2,X:0.6):0.2;"
###
divscale=1.0


###Edits by JBP 
#Changed branch lengths
nodeABX          = 20.0*divscale
nodeAB           = 4.0*divscale
nodeABC          = 1.0*divscale

# sets the two parameter classes
###Edits by JBP
#Changed simulation to include recombination, migration within population pairs, and to include proper fusion events to model tree above
paramSet = egglib.simul.CoalesceParamSet(singleSamples=None, doubleSamples=[A,B,C,D,X],M=0.0, rho=10.0, nsites=100)
"clade 1"
paramSet.populationFusion(nodeAB, 0,1)       ## B into A
paramSet.populationFusion(nodeAB, 2,3) 	     ## C into A
paramSet.populationFusion(nodeABX, 0,2)       ## D into A
paramSet.changePairwiseMigrationRate(0.1, 0, 1, 4*N*0.01)
paramSet.changePairwiseMigrationRate(0.1, 3, 2, 4*N*0.01)
paramSet.changePairwiseMigrationRate(0.1, 2, 3, 4*N*0.01)
paramSet.changePairwiseMigrationRate(0.1, 1, 0, 4*N*0.01)
"together and outgroup"
paramSet.populationFusion(nodeABX, 0,4) ## X into A

mutator = egglib.simul.CoalesceFiniteAlleleMutator(theta=theta,
                                                   alleles= 4,
                                                   randomAncestralState=True)
mutator.setSites(locuslength)


def twodiffs(a,b):
    "requires two base differences between barcodes"
    if len(a) == len(b):
        t = [a[i]==b[i] for i in range(len(a))]
        if t.count(False) > 1:
            return True


def barcoder(aligns, Ninds, tiptax, outgroup, outhandle, Barcodes):
    """ takes an align object and puts names to sequences and
    makes a barcode map """

    " append names"
    names = []
    for name in tiptax:
        for i in range(Ninds):
            names.append(name+str(i))
            names.append(name+str(i))

    " appends names for allele 1 vs. allele 2 for each diploid sample "
    l = 1
    for i in aligns:
        for s,j in zip(i,names):
            s.name = j
    
    " make dictionary with list of loci for each sample "
    D = {}
    for i in set(names):
        D[i] = []
    for samp in D:
        for i in aligns:
            l = []
            for j in i:
                if j.name == samp:
                    l.append(j.sequence)
            D[samp].append(l) 

    " barcodes to append to start of loci "
    bases = list("ATGC")
    if not Barcodes:
        print "creating new barcode map"
        Barcodes[names[0]] = "CATCAT"
        for name in names:
            while name not in Barcodes.keys():
                bbs = numpy.random.randint(0,3,6)
                bar = "".join([bases[i] for i in bbs])
                if all([twodiffs(bar,bb) for bb in Barcodes.values()]):
                    if not any([i in bar for i in [CUT1,CUT2]]):
                        Barcodes[name] = bar
        " creates random barcodes and writes map to file "
        with open(outhandle+".barcodes",'w') as barout:
            bnames = list(Barcodes.keys())
            bnames.sort()
            for bar in bnames:
                if outgroup not in bar:
                    print >>barout, "\t".join([bar,Barcodes[bar]])
    return D,Barcodes



def mutate(base):
    "introduce point sequencing errors"
    nbase = list("ATGC")[numpy.random.randint(0,4)]
    while nbase == base:
        nbase = list("ATGC")[numpy.random.randint(0,4)]
    return nbase


def revcomp(seq):
    S = seq.replace('A','t').replace("T",'a').replace("C",'g').replace("G",'c')
    return S[::-1].upper()


def dressitup(D, Barcodes, sitemut, indel, locuslength, divis, addon, frags):
    """ return a list of sequences all dressed up """
    SEQ1s     = []
    SEQ2s     = []
    Illumina_P1_adapter = "ACGACGCTCTTCCGATCT"
    Illumina_P2_adapter = "AGATCGGAAGAGCTCGTATG"
    for loc in range(divis):
        f1,f2 = frags.split(',')
        frag = numpy.random.randint(int(f1),int(f2))

        " can change to random sequencing depth. Even number 2-20 "
	#copies  = numpy.random.randint(2,42,1)
	#while copies % 2:
        #     copies  = numpy.random.randint(2,42,1)

        if indel:
            skip = 5
            "indel = probability site is checked for indel-mutation, skip first 'skip' bp"
            inds = numpy.random.binomial(locuslength-skip,indel)
            "if indel, for each indel site randomly choose location on locus"
            wh = [numpy.random.randint(skip,locuslength) for i in range(inds)]
            "randomly select size of indel (default size 1)"
            le = [numpy.random.randint(1,2) for i in wh]
            where = [wh,le]

        for samp in D:
            """ exclude the outgroup taxon used for detecting restriction site mutations and indels
            use both chromosomes """
            ###Edits by JBP
            #Changed the random read distribution to a gamma distribution instead of normal
            copies  = max(int(numpy.random.gamma(1.6,20)-2),0)
            while copies % 2:
                copies = max(int(numpy.random.gamma(1.6,20)-2),0)
                #copies  = numpy.random.randint(2,42,1)

            if sitemut:
                " remove loci if mutation arises in restriction site"
                #dropout = len(CUT+1)
                dropout = sitemut
                L1 = [list(D[samp][loc][0]) for copy in range(copies/2) if D[samp][loc][0][0:dropout] == D['4X0'][loc][0][0:dropout]] 
                L2 = [list(D[samp][loc][1]) for copy in range(copies/2) if D[samp][loc][0][0:dropout] == D['4X0'][loc][0][0:dropout]] 
                if datatype not in ['rad']:
                    "if double digested (ddrad/gbs) check other end for locus dropout "
                    #dropout2 = -len(CUT2)-1
                    dropout2 = -1*sitemut
                    L1 = [list(D[samp][loc][0]) for copy in range(copies/2) if D[samp][loc][0][dropout2:] == D['4X0'][loc][0][dropout2:]] 
                    L2 = [list(D[samp][loc][1]) for copy in range(copies/2) if D[samp][loc][0][dropout2:] == D['4X0'][loc][0][dropout2:]]
            else:
                L1 = [list(D[samp][loc][0]) for copy in range(copies/2)] 
                L2 = [list(D[samp][loc][1]) for copy in range(copies/2)]
            L = L1+L2

            for copy in range(len(L)):
                ll = locuslength
                """ introduce indel if mutation present at site
                relative to the outgroup """
                if indel:
                    for ww,ee in zip(where[0],where[1]):
                        df = locuslength-ll
                        ww -= df
                        if L[copy][ww] != D['4X0'][loc][0][ww]:
                            L[copy] = L[copy][:ww] + L[copy][ww+ee:]
                            ll -= ee

                " introduce sequencing errors"
                Es = numpy.random.binomial(ll,0.001)
                Elocs = numpy.random.randint(0,ll,Es)
                eL = L[copy]
                for err in Elocs:
                    eL[err] = mutate(eL[err])

                " shorten fragment to potentially overlapping length"
                eL = eL[:frag]

                " make fastQ formatted"
                if datatype in ['ddrad','pairddrad']:
                    eL1 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+CUT1+"".join(eL)+revcomp(CUT2)+Illumina_P2_adapter+"A"*30
                else:
                    eL1 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+CUT1+"".join(eL)+revcomp(CUT1)+Illumina_P2_adapter+"A"*30

                if datatype in ['gbs']:
                    eL2 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+CUT1+revcomp("".join(eL))+revcomp(CUT1)+Illumina_P2_adapter+"A"*30
                    
                " sequence read1 from 5' end "
                startL = eL1.index("CGATCT")
                "---->"
                ii = startL+len("CGATCT")
                sss = eL1[ii:ii+100]
                if "X" not in samp:
                    SEQ1s.append("@lane1_fakedata"+str(loc+addon)+"_"+str(copy)+" 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

                if datatype in ['gbs']:
                    " sequence read1 from 3' end for GBS"
                    startL = eL2.index("CGATCT")
                    "---->"
                    ii = startL+len("CGATCT")
                    sss = eL2[ii:ii+100]
                    if "X" not in samp:
                        SEQ1s.append("@lane1_fakedata"+str(loc+addon)+"_"+str(copy)+" 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

                " sequence read2 from 3' end for paired read "
                startR = eL1.rindex("AGATCG")
                "<----"
                ii = len(eL1)-startR
                sss = eL1[-1*(ii+100):-ii]
                if "X" not in samp:
                    SEQ2s.append("@lane1_fakedata"+str(loc+addon)+"_"+str(copy)+" 1:N:0:"+"\n"+revcomp(sss)+"\n+\n"+("B"*len(sss))+"\n")

                if datatype in ['pairgbs']:
                    ' double check this '
                    " sequence read2 from 5' end for GBS"
                    startL = eL2.index("CGATCT")
                    "---->"
                    ii = startL+len("CGATCT")
                    sss = eL2[ii:ii+100]
                    if "X" not in samp:
                        SEQ2s.append("@lane1_fakedata"+str(loc+addon)+"_"+str(copy)+" 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")
                    

    return SEQ1s,SEQ2s




print '\t',
D = []
Barcodes = {}

out1 = gzip.open(outhandle+"_R1_.fastq.gz",'wb')
if datatype in ['pairddrad','pairgbs']:
    out2 = gzip.open(outhandle+"_R2_.fastq.gz",'wb')

if datatype in ['gbs']:
    nrepets = nrepets/2  # because forward and reverse sequenced

for i in range(0, nrepets, divis):
    " performs simulation of 1K loci at a time"
    aligns = egglib.simul.coalesce(paramSet, mutator, divis)
    " creates barcodes if not made yet and fill D dict."
    D,Barcodes = barcoder(aligns, Ninds, tiptax, outgroup, outhandle, Barcodes)
    " dresses up data to be fastq and puts in errors, indels, etc."
    SEQ1s,SEQ2s = dressitup(D, Barcodes, sitemut, indel, locuslength, divis, i, frags)
    out1.write("".join(SEQ1s))
    if datatype in ['pairddrad','pairgbs']:
        out2.write("".join(SEQ2s))
    pol_data = map(egglib.Align.polymorphism, aligns)
    Fst = [i['Fst'] for i in pol_data]
    print 'average Fst:', sum(Fst)/nrepets
    print '.',
out1.close()
if datatype in ['pairddrad','pairgbs']:
    out2.close()
