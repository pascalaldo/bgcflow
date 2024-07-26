#%
# final_output: data/processed/{name}/fastani/df_fastani.csv
# description: Do pairwise Average Nucleotide Identity (ANI) calculation across all
#   samples.
# category: QC and Data Selection
# link:
# - https://github.com/ParBLiSS/FastANI
# references:
# - 'Jain, C., Rodriguez-R, L.M., Phillippy, A.M. et al. High throughput ANI analysis
#   of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun 9, 5114
#   (2018). [https://doi.org/10.1038/s41467-018-07641-9](https://doi.org/10.1038/s41467-018-07641-9)'
#%
rule fastani:
    input:
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fastani_infile="data/interim/fastani/{name}/fastani_in.txt",
        fastani_out="data/interim/fastani/{name}/fastani_out.tsv",
        fastani_matrix="data/interim/fastani/{name}/fastani_out.tsv.matrix",
    conda:
        "../envs/fastani.yaml"
    threads: 32
    log:
        "logs/fastani/fastani-{name}.log",
    shell:
        """
        for fna in {input.fna}
        do
            echo $fna >> {output.fastani_infile}
        done
        fastANI --ql {output.fastani_infile} --rl {output.fastani_infile} -t {threads} --matrix -o {output.fastani_out} 2>> {log}
        """


rule fastani_convert:
    input:
        fastani_matrix="data/interim/fastani/{name}/fastani_out.tsv.matrix",
    output:
        df_fastani="data/processed/{name}/fastani/df_fastani.csv",
    log:
        "logs/fastani/fastani-convert-{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/convert_triangular_matrix.py {input.fastani_matrix} {output.df_fastani} 2>> {log}
        """
