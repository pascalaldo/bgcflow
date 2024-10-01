#%
# final_output: data/processed/{name}/tables/df_seqfu_stats.csv
# description: Calculate sequence statistics using SeqFu.
# category: QC and Data Selection
# link:
# - https://github.com/telatin/seqfu2
# references:
# - 'Telatin, A., Birolo, G., & Fariselli, P. SeqFu [Computer software]. GITHUB: [https://github.com/telatin/seqfu2](https://github.com/telatin/seqfu2)'
#%
rule seqfu_stats:
    input:
        fna="data/interim/all/fasta/{strains}.fna",
    output:
        json="data/interim/{stage}/seqfu/{strains}.json",
    conda:
        "../envs/seqfu.yaml"
    log:
        "logs/{stage}/seqfu/seqfu/seqfu-{strains}.log",
    params:
        precision=3,
    shell:
        """
        seqfu stats {input.fna} --json -b --gc --precision 3 > {output.json} 2>> {log}
        """


rule seqfu_combine:
    input:
        json=fexpand(
            "data/interim/{{stage}}/seqfu/{strains}.json",
            strains=RULE_FUNCTIONS["seqfu"]["strains"],
        ),
    output:
        all_csv=report(
            "data/processed/{stage}/{name}/tables/df_seqfu_stats.csv",
            caption="../report/table-seqfu.rst",
            category="{name}",
            subcategory="Quality Control",
        ),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/{stage}/seqfu/seqfu-{name}.log",
    shell:
        """
        TMPDIR="data/interim/{wildcards.stage}/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_seqfu.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/make_seqfu_table.py $INPUT_JSON {output.all_csv} &>> {log}
        rm $INPUT_JSON
        """
