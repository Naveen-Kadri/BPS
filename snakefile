localrules:cat, cat_seq_info, cat_mutation
OUT_DIR ="/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN5/"


DATA=["NEW"]
chromosomes=range (1,30)
myrange = [-20, 29]

#without Ns
genome_mut= 100 * 29_418_137/2_489_368_272



gtf_file="/cluster/work/pausch/inputs/ref/BTA/UCD1.2/ensembl_104/Bos_taurus.ARS-UCD1.2.104.chr.gtf3.gz"
fasta_file="/cluster/work/pausch/naveen/CNS/REFFASTA/chr{chr}.fa"
frq_file="/cluster/work/pausch/naveen/CNS/GENOME2/FRQ/NEW/CHR{chr}/my.frq"
ac_file="/cluster/work/pausch/naveen/CNS/GENOME2/FRQ/NEW/CHR{chr}/alt.frq"


#files for bpp
BPP="/cluster/work/pausch/group_bin/BPP/BP_PPT.py"

##bpp training files
scPPT_file="/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/BTA/scPPT_bta.txt"
pwmBP_file="/cluster/work/pausch/group_bin/BPP/demo/pwmBP_human.txt"


rule all:
    input:
        OUT_DIR + "seq_info.txt",
        OUT_DIR + "seq.txt",
        OUT_DIR + "BP_logo.pdf",
        expand(OUT_DIR + "MUTATION/{data}/mutation_count{chr}.txt",data=DATA,chr=chromosomes),
        expand(OUT_DIR + "MUTATION/{data}/mutation_type{chr}.txt", data=DATA,chr=chromosomes),
        expand(OUT_DIR + "MUTATION/{data}/bps_mutation_{data}.pdf",data=DATA),
        expand(OUT_DIR + "{data}/distance_to_3ss.pdf",data=DATA),
        expand(OUT_DIR + "{data}/intron_bp_types.pdf",data=DATA)
        

        
rule extract_seq:
    input:
        fasta_file=fasta_file,
        gtf_file=gtf_file
    params:
        mychr="{chr}",
        min_intronic_bases=20,
        max_intronic_bases=150

    output:
        allinfo=OUT_DIR + "CHR{chr}/all_info.txt",
        intronic_seq=OUT_DIR + "CHR{chr}/seq.txt",
        seq_info_file=OUT_DIR + "CHR{chr}/seq_info.txt"
    script:
        "extract_seq.py"
        
rule cat:
    input:
        infiles=expand (OUT_DIR + "CHR{chr}/seq.txt", chr=chromosomes )
    output:
        outfile=OUT_DIR + "seq.txt"
    shell:
        "cat {input.infiles} > {output.outfile}"

rule cat_seq_info:
    input:
        infiles=expand (OUT_DIR + "CHR{chr}/seq_info.txt", chr=chromosomes)
    output:
        outfile=OUT_DIR + "seq_info.txt"
    script:
        "cat_info.R"

rule bpp:
    input:
        seq=OUT_DIR + "seq.txt",
        scPPT=scPPT_file,
        pwmBP=pwmBP_file
    output:
        OUT_DIR + "result_bpp.txt"
    shell:
        BPP + " -i {input.seq} -p {input.scPPT} -b {input.pwmBP} -r 1 > {output}"

rule split_res:
    input:
        OUT_DIR + "result_bpp.txt"
    params:
        OUT_DIR + "CHR"
    output:
        expand(OUT_DIR + "CHR{chr}/bpp_res.txt",chr=chromosomes)
    script:
        "split_result.py"

rule format_for_mutation:
    input:
        res_file=OUT_DIR + "CHR{chr}/bpp_res.txt",
        seq_file=OUT_DIR + "CHR{chr}/seq.txt",
        fasta_file=fasta_file,
    output:
        out_file=OUT_DIR + "CHR{chr}/feature.txt"
    script:
        "format_for_mutation.py"

rule cat_bpp_res:
    input:
        in_files=expand (OUT_DIR + "CHR{chr}/bpp_res.txt",chr=chromosomes)
    output:
        out_file=OUT_DIR + "formatted_bpp_res.txt"
    script:
        "cat_bpp_res.py"
        
rule plot_logo:
    input:
        in_file=OUT_DIR + "formatted_bpp_res.txt"
    output:
        plot_file=OUT_DIR + "BP_logo.pdf"
    script:
        "logo.R"
        
        
rule count_variation:
    input:
        position_file=OUT_DIR + "CHR{chr}/feature.txt",
        frq_file=frq_file
    output:
        count_file=OUT_DIR + "MUTATION/{data}/mutation_count{chr}.txt",
        type_file=OUT_DIR+"MUTATION/{data}/mutation_type{chr}.txt",
    params:
        myrange=myrange
    script:
        "count_variation.py"

# rule count_ac:
#     input:
#         position_file=OUT_DIR + "CHR{chr}/feature.txt",
#         frq_file=ac_file
#     output:
#         count_file=OUT_DIR + "AC/{data}/ac_count{chr}.txt"
#     params:
#         myrange=myrange
#     script:
#         "count_ac.py"
        
# rule cat_AF:
#     input:
#         count_files=expand (OUT_DIR + "AC/{{data}}/ac_count{chr}.txt", chr=chromosomes  ),
#     output:
#         count_file=OUT_DIR + "AC/{data}/ac_count.txt"
#     run:
#         shell ("cat {input.count_files} > {output.count_file}")

        
rule cat_mutation:
    input:
        count_files=expand (OUT_DIR + "MUTATION/{{data}}/mutation_count{chr}.txt", chr=chromosomes  ),
        type_files=expand (OUT_DIR+"MUTATION/{{data}}/mutation_type{chr}.txt",chr=chromosomes)
    output:
        count_file=OUT_DIR + "MUTATION/{data}/mutation_count.txt",
        type_file=OUT_DIR + "MUTATION/{data}/mutation_type.txt"
    run:
        shell ("cat {input.count_files} > {output.count_file}"),
        shell ("cat {input.type_files} >  {output.type_file}")
        
rule combine_info:
    input:
        info_file=OUT_DIR + "seq_info.txt",
        mutation_file=OUT_DIR + "MUTATION/{data}/mutation_count.txt",
        bpp_res=OUT_DIR + "formatted_bpp_res.txt"
    output:
        out_file=OUT_DIR + "MUTATION/{data}/allinfo_mutation_count.txt"
    params:
        myrange=myrange
    script:
        "combine_all_info_mutation_count.R"


rule combine_info_AC:
    input:
        info_file=OUT_DIR + "seq_info.txt",
        mutation_file=OUT_DIR + "AC/{data}/ac_count.txt",
        bpp_res=OUT_DIR + "formatted_bpp_res.txt"
    output:
        out_file=OUT_DIR + "AC/{data}/allinfo_mutation_count.txt"
    params:
        myrange=myrange
    script:
        "combine_all_info_mutation_count.R"


        
rule plot_mut_rate:
    input:
        in_file=rules.combine_info.output.out_file
    output:
        plot_file=OUT_DIR + "MUTATION/{data}/bps_mutation_{data}.pdf"
    params:
        myrange=myrange,
        genome_mut=genome_mut
        
    script:
        "plot_mutation_rate.R"
    
rule plot_distance:
    input:
        in_file=rules.combine_info.output.out_file
    output:
        plot_file=OUT_DIR + "{data}/distance_to_3ss.pdf"
    script:
        "plot_distance_from_3ss.R"


rule plot_intron_types:
    input:
        in_file=rules.combine_info.output.out_file
    output:
        plot_file=OUT_DIR + "{data}/intron_bp_types.pdf"
    script:
        "plot_intron_bp_types.R"
        
        
