#%
# final_output: data/processed/{name}/tables/df_checkm_stats.csv
# description: Assess genome quality with CheckM.
# category: QC and Data Selection
# link:
# - https://github.com/Ecogenomics/CheckM
# references:
# - 'Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2014. Assessing
#   the quality of microbial genomes recovered from isolates, single cells, and metagenomes.
#   [Genome Research, 25: 1043-1055.](https://genome.cshlp.org/content/25/7/1043.long)'
#%
rule install_checkm:
    output:
        checkm_db=directory("resources/checkm/"),
    conda:
        "../envs/checkm.yaml"
    log:
        "logs/checkm/checkm-install_checkm.log",
    shell:
        """
        (cd resources && wget -O checkm_data_2015_01_16.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz) &>> {log}
        (cd resources && mkdir -p checkm && tar -xvf checkm_data_2015_01_16.tar.gz -C checkm && rm checkm_data_2015_01_16.tar.gz) &>> {log}
        checkm data setRoot resources/checkm &>> {log}
        """


rule checkm:
    input:
        fna=lambda wildcards: expand("data/interim/all/fasta/{accession}.fna", accession=RULE_FUNCTIONS["checkm"][wildcards.stage]["accessions"](wildcards.name)),
        checkm_db="resources/checkm/",
    output:
        fna=temp(directory("data/interim/{stage}/checkm/{name}_fna")),
        stat="data/interim/{stage}/checkm/{name}/storage/bin_stats_ext.tsv",
        checkm_dir=directory("data/interim/{stage}/checkm/{name}"),
    conda:
        "../envs/checkm.yaml"
    log:
        "logs/{stage}/checkm/checkm_{name}.log",
    params:
        checkm_log="data/interim/{stage}/checkm/{name}/checkm_{name}.log",
    threads: 16
    shell:
        """
        mkdir -p {output.fna}
        for f in {input.fna}; do cp $f {output.fna}/.; done
        checkm lineage_wf -t {threads} --reduced_tree -x fna {output.fna} {output.checkm_dir} &>> {log}
        """


rule checkm_out:
    input:
        stat="data/interim/{stage}/checkm/{name}/storage/bin_stats_ext.tsv",
    output:
        stat_processed=report(
            "data/processed/{stage}/{name}/tables/df_checkm_stats.csv",
            caption="../report/table-checkm.rst",
            category="Quality Control",
        ),
    log:
        "logs/{stage}/checkm/checkm_out_{name}.log",
    params:
        checkm_json=directory("data/interim/{stage}/checkm/{name}/json/"),
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        mkdir -p {params.checkm_json}
        python workflow/bgcflow/bgcflow/data/get_checkm_data.py {input.stat} {params.checkm_json} {output.stat_processed} 2>> {log}
        """
