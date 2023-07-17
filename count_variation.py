'''
only snps are considered !
Negative vlues are upstream on the +ve strand
alleles on the -ve strands are flipped !

'''

import os
import numpy as np
from collections import defaultdict

def reverse_complement(seq):
    comp={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", ",":","}
    return ("".join ([comp.get (el, "N") for el in seq]))
    ##return "".join([comp.get (base,"N") for base in seq[::-1]])



position_file=snakemake.input.position_file
frq_file=snakemake.input.frq_file
count_out_file=snakemake.output.count_file
type_out_file=snakemake.output.type_file
myrange=snakemake.params.myrange



# position_file="/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN/CHR1/feature.txt"
# frq_file="/cluster/work/pausch/naveen/CNS/GENOME2/FRQ/NEW/CHR1/my.frq"
# count_out_file="test_mutation.txt"
# type_out_file="test_mutation_type.txt"
# myrange=[-51,21]


##frequency and allele information
print ('reading the maf file')
mafs=dict ()
alleles=defaultdict(list)
with open (frq_file, "rt") as inf:
    for lnum,line in enumerate(inf):
        if lnum >0:
            mychr,pos,vartype,reflen,altlen,nal,ref,alt,maf=line.rstrip().split()
            #mychr, snp, alt, ref, af, nobs=line.rstrip().split()
            if vartype == "SNP" and "*" not in ref and "*" not in alt:
                mypos=int(pos)
                maf=float (maf)
                mafs[mypos]=maf
                alleles [mypos]=[ref, alt]
count_out=open (count_out_file, "w")
type_out=open (type_out_file, "w")


btypes =dict ()
can_mut_type = dict ()
noncan_mut_type = dict ()
n=0
print ('reading the feature file')
with open (position_file, "rt") as inf:
    for line in inf:
        mychr, mypos, strand,myid,bps,score=line.rstrip().split()
        mycord = mychr + "_" + mypos
        mutations = []
        if bps [3] == "T" and bps [5] == "A":
            btype='canonical'
        else:
            btype='noncanonical'
        btypes [mycord] = btype
        for dist in range(myrange[0],myrange [1]):
            if strand == "+":
                mybase=int(mypos) + dist
                if mybase in alleles :
                    myref=alleles [mybase] [0]
                    myalt=alleles [mybase] [1]
            else:
                mybase=int(mypos) - dist
                if mybase in alleles:
                    myref=reverse_complement(  alleles[mybase] [0]   )
                    myalt=reverse_complement(  alleles[mybase] [1]   )

            if mafs.get (mybase,1) != 1.0 :
                mutations.append (1)
                my_type="\t".join([myref, myalt])
                tw = "\t".join ([mycord, myid, bps, btype, str(dist), my_type])
                type_out.write (f"{tw}\n")
            else:
                mutations.append (0)
        mutations=[str(el) for el in mutations]
        tw=[mycord, myid, bps, btype] + mutations
        tw ="\t".join (tw)
        count_out.write (f"{tw}\n")

type_out.close ()
count_out.close()
