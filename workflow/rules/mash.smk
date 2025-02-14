#%
# final_output: data/processed/{name}/mash/df_mash.csv
# description: Calculate distance estimation for all samples using MinHash.
# category: QC and Data Selection
# link:
# - https://github.com/marbl/Mash
# references:
# - 'Mash: fast genome and metagenome distance estimation using MinHash. Ondov BD,
#     Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. [Genome
#     Biol. 2016 Jun 20;17(1):132. doi: 10.1186/s13059-016-0997-x.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)'
# - 'Mash Screen: high-throughput sequence containment estimation for genome discovery.
#     Ondov BD, Starrett GJ, Sappington A, Kostic A, Koren S, Buck CB, Phillippy AM.
#     [Genome Biol. 2019 Nov 5;20(1):232. doi: 10.1186/s13059-019-1841-x.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1841-x)'
#%
rule mash:
    input:
        fna=fexpand("data/interim/all/fasta/{accession}.fna", accession=RULE_FUNCTIONS["mash"]["accessions"]),
        # fna=lambda wildcards: get_fasta_inputs(wildcards.name, get_samples_df()),
    output:
        mash_infile="data/interim/{stage}/mash/{name}/mash_in.txt",
        triangle_dist="data/interim/{stage}/mash/{name}/triangle_distance_matrix.tsv",
    conda:
        "../envs/mash.yaml"
    threads: 32
    params:
        sketch_size=1000,
    log:
        "logs/{stage}/mash/mash-triangle-{name}.log",
    shell:
        """
        for fna in {input.fna}
        do
            echo $fna >> {output.mash_infile}
        done
        (mash triangle -s {params.sketch_size} -p {threads} -l {output.mash_infile} >> {output.triangle_dist}) 2>> {log}
        """


rule mash_convert:
    input:
        mash_matrix="data/interim/{stage}/mash/{name}/triangle_distance_matrix.tsv",
    output:
        df_mash="data/processed/{stage}/{name}/mash/df_mash.csv",
    log:
        "logs/{stage}/mash/mash-convert-{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/convert_triangular_matrix.py {input.mash_matrix} {output.df_mash} 2>> {log}
        """
