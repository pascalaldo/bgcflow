def get_genome_ids(name, df_samples):
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    return selection

rule pankb_isosource:
    input:
        qc.output
    params:
        genome_ids=lambda wildcards: get_genome_ids(wildcards.name, filter_samples_qc(wildcards, DF_SAMPLES)),
    output:
        isosource="data/processed/{name}/pankb/source_info/df_ncbi_isolation_src.csv"
    log:
        "logs/pankb_data_prep/pankb_isosource_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep isosource {params.genome_ids} \
            -o {output.isosource} > {log} 2>&1
        """
