from collections import defaultdict

res_file = snakemake.input[0]
prefix=snakemake.params [0]

RES = defaultdict(list)
header="\t".join   ( ["myid", "bps", "bps_pos", "sc_bps", "sc_ppt", "sc", "zsc", "zsc_ppt", "zsc"])

with open (res_file) as inf:
    for line in inf:
        if line[0]!="#":
            myid, bps, bp_pos, sc_bps, sc_ppt, sc, zsc, zsc_ppt, zsc=line.rstrip().split ()
            transcript, intron_nr, mychr, strand, mystart, myend=myid.split("_")
            RES [mychr].append  ( line)


for mychr in RES:
    print (mychr)
    myres = RES [mychr]
    with open (prefix + mychr + "/bpp_res.txt", "w") as out:
        out.write (f"{header}\n")
        for el in myres:
            out.write (f"{el}")
        
