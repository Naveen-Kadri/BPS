"""
all exons are considered [non-coding exons included]
only unique sequneces are retained....  as there are ~10k introns that are duplicated w.r.t intron ends
>ENSBTAT00000048616_11_16_+_44406237_44406327 and
>ENSBTAT00000009606_11_14_+_82386589_82386679 have the same intronic sequneces they are on bta 14 and 16
so the uniqueness is checked only within a chromosome
there are 33 duplicated sequences even after removing duplicates within chromosomes
"""
import gzip
from collections import defaultdict
gtf_file=snakemake.input.gtf_file
fasta_file=snakemake.input.fasta_file

mychr=str(snakemake.params.mychr)
min_intronic_bases=snakemake.params.min_intronic_bases
max_intronic_bases=snakemake.params.max_intronic_bases


out_all_info_file=snakemake.output.allinfo
out_seq_file=snakemake.output.intronic_seq
out_seq_info_file=snakemake.output.seq_info_file


out_all_info = open (out_all_info_file, "w")
header="\t".join (["gene_id", "gene_name", "transcript", "strand", "intron", "start", "end", "SS5", "SS3"])
out_all_info.write (f"{header}\n")

out_seq_info = open (out_seq_info_file, "w")
out_seq_info.write (f"{header}\n")

out_seq=open (out_seq_file, "w")



gene_info=dict()
strands=dict()
transcripts = dict ()
exons=defaultdict(dict)
seen=dict()


def reverse_complement(seq):
    comp={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    return "".join([comp.get (base) for base in seq[::-1]])


contig=""    
print ("reading the fasta file")
with open (fasta_file, "rt") as inf:
    for line in inf:
        if line[0]!=">":
            contig+=line.rstrip()
print (len (contig))
contig="0"+contig  #to avoid confusion of zero-indexing

#function to collect all info from the gtf
def collect_info(attributes):
    attris=attributes.split(";")
    myattris=dict()
    for el in attris:
        splits=el.lstrip().split (" ")
        if len(splits)==2:
            myattris[splits[0]]=splits[1].strip('"')
    return (myattris)

infos_to_print =['transcript', 'exon']


with gzip.open(gtf_file, "rt") as inf:
    for line in inf:
        seqid,source,ftype,start, end,score,strand, phase,attributes=line.rstrip().split("\t")
        if seqid != mychr:
            continue
        start=int (start)
        end = int (end)
        if ftype in infos_to_print:
            myinfo=collect_info(attributes)
            mytr = myinfo.get("transcript_id")
            gene_id=myinfo.get ("gene_id")
            gene_name=myinfo.get ("gene_name", "NA")
            gene_info[mytr] = "\t".join ([gene_id, gene_name])
            
            if strand == "-":
                start, end = end, start
            if ftype == "transcript":
                transcripts [mytr] = [start, end]
                strands [mytr] = strand
            elif ftype=="exon" and myinfo.get ("gene_biotype") == "protein_coding":
                exon_nr = int (myinfo.get ("exon_number"))
                exons  [mytr]  [exon_nr] = [start, end]

#gene_biotype "protein_coding"
for tr in exons:
    myexons=exons[tr]
    mystrand=strands [tr]
    my_gene_info=gene_info[tr]
    myexons = sorted  (myexons.items()) 
    myexons=list (map (lambda tup : tup [1], myexons))
    nex = len (myexons) -1
    myintrons = []
    
    if mystrand == "+":
        for i in range (0, nex):
            myintrons.append ([   myexons [i] [1] +1, myexons [i+1] [0] - 1   ])
    else:
        for i in range (0, nex):
            myintrons.append ([   myexons [i] [1] -1, myexons [i+1] [0] + 1   ])
                
    for k, exon in enumerate(myexons):
        mystart=exon[0]
        myend=exon [1]
        
        if mystrand == "+":
            ss5=contig [mystart:mystart+2]
            ss3=contig [myend-1:myend+1]
        else:
            ss5=reverse_complement(contig [mystart-1:mystart+1])
            ss3=reverse_complement(contig [myend:myend+2])
                                    

        tw="\t".join ([my_gene_info  ,tr, mystrand, "exon" + str (k+1), str(exon[0]), str(exon [1]), ss5, ss3 ])
        out_all_info.write (f"{tw}\n")
        
    for k, intron in enumerate(myintrons):
        mystart=intron[0]
        myend=intron [1]
        #mycord=str(mystart) + "_"  + str (myend)

        intron_length = abs(myend - mystart) + 1
        if mystrand == "+":
            ss5=contig [mystart:mystart+2]
            ss3=contig [myend-1:myend+1]
            intronic_seq = contig [mystart : myend+1]
        else:
            ss5=reverse_complement(contig [mystart-1:mystart+1])
            ss3=reverse_complement(contig [myend:myend+2])
            intronic_seq = reverse_complement(contig [myend : mystart+1])

        tw="\t".join ([my_gene_info,tr, mystrand, "intron" + str (k+1), str(mystart), str(myend), ss5, ss3 ])
        out_all_info.write (f"{tw}\n")
        
        
        nbases=min (intron_length,  max_intronic_bases)
        if nbases >= min_intronic_bases:
            intronic_seq_sub=intronic_seq[-1*nbases:]
            if 'NN' not in intronic_seq_sub:
                if intronic_seq_sub not in seen:
                    out_seq.write (f">{tr}_{k+1}_{mychr}_{mystrand}_{mystart}_{myend}\n")
                    out_seq.write (f"{intronic_seq_sub}\n")
                    out_seq_info.write (f"{tw}\n")
                    seen[intronic_seq_sub]=1
                
                        
out_all_info.close ()
out_seq.close ()
out_seq_info.close ()
