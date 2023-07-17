in_files = snakemake.input.in_files
out_file=snakemake.output.out_file

out=open (out_file, 'w')
for myfile in in_files:
    with open (myfile) as inf:
        for i, line in enumerate (inf):
            if i>0:
                out.write (f"{line}")

out.close()




