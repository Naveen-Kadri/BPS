'''
some additional checks --
writes a file with anchor at one bp upstream of the bps [7bps]
this eases visulisation of mutation load  [and others]
T at pos 4 and A at pos 6


'''
# res_file = "/cluster/work/pausch/naveen/CNS/BPS/BPSII/CHR1/bpp_res.txt"
# seq_file="/cluster/work/pausch/naveen/CNS/BPS/BPSII/CHR1/seq.txt"
# fasta_file="/cluster/work/pausch/naveen/CNS/REFFASTA/chr1.fa"
# out_file="todel.txt"
# out =open (out_file, "w")


res_file = snakemake.input.res_file
seq_file= snakemake.input.seq_file
fasta_file=snakemake.input.fasta_file
out_file=snakemake.output.out_file



#print ("extract contig from the  fasta file")
def reverse_complement(seq):
    comp={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    return "".join([comp.get (base) for base in seq[::-1]])
contig=""    
print ("reading the fasta file")
with open (fasta_file, "rt") as inf:
    for line in inf:
        if line[0]!=">":
            contig+=line.rstrip()
#print (len (contig))
contig="0"+contig  #to avoid confusion of zero-indexing



print ("hash all intronic sequences")
seq = dict ()
with open (seq_file) as inf:
    for line in inf:
        if line[0]==">":
            myid=line.rstrip()
            transcript, intron_nr, mychr, strand, mystart, myend=myid.split("_")
        else:
            seq [myid]= line
            myid=""

print ("read the result file")
out =open (out_file, "w")


with open (res_file) as inf:
    for lnum,line in enumerate(inf):
        if lnum>0:
            myid, bps, bp_pos, sc_bps, sc_ppt, sc, zsc, zsc_ppt, zsc=line.rstrip().split ()
            transcript, intron_nr, mychr, strand, mystart, myend=myid.split("_")
            bp_pos= int (bp_pos)
            myend=int (myend)
            myseq=seq[myid]
            if strand == "+":
                st=myend - (bp_pos-1+5)
                en=st+ 6
                checkseq=contig [st: (en +1)]
                feature_pos=st-1
            elif strand == "-":
                en=myend +  (bp_pos -1 + 5)  #end in python should be +1 [not included]

                st=en-6
                checkseq=reverse_complement(contig [st:(en+1)])
                feature_pos=en+1
            if bps == checkseq:
                out.write (f"{mychr}\t{feature_pos}\t{strand}\t{myid}\t{bps}\t{zsc}\n")
            else:
                print ("problem!")
                exit ()
                



                

