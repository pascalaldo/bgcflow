#%
# final_output: data/processed/{name}/automlst_wrapper/final.newick
# description: Simplified Tree building using [autoMLST](https://github.com/NBChub/automlst-simplified-wrapper)
# category: Phylogenomic Placement
# link:
# - https://github.com/KatSteinke/automlst-simplified-wrapper
# references:
# - 'Mohammad Alanjary, Katharina Steinke, Nadine Ziemert, AutoMLST: an automated
#   web server for generating multi-locus species trees highlighting natural product
#   potential,[Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W276–W282](https://doi.org/10.1093/nar/gkz282)'
#%
rule install_automlst_wrapper:
    output:
        folder=directory("resources/automlst-simplified-wrapper-main"),
        reduced_core="resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    conda:
        "../envs/automlst_wrapper.yaml"
    log:
        "logs/automlst_wrapper/install_automlst_wrapper.log",
    params:
        source="https://github.com/NBChub/automlst-simplified-wrapper",
        version="0.1.2"
    shell:
        """
        set -e
        mkdir -p resources 2>> {log}
        wget {params.source}/archive/refs/tags/v{params.version}.zip -O resources/automlst-simplified-wrapper-v{params.version}.zip 2>> {log}
        (cd resources && unzip -o automlst-simplified-wrapper-v{params.version}.zip && rm automlst-simplified-wrapper-v{params.version}.zip) &>> {log}
        (cd resources/automlst-simplified-wrapper-{params.version} && unzip -o reducedcore.zip) &>> {log}
        cp -r resources/automlst-simplified-wrapper-{params.version}/* {output.folder}/. 2>> {log}
        rm -rf resources/automlst-simplified-wrapper-{params.version} 2>> {log}
        """


rule prep_automlst_gbk:
    input:
        gbk="data/interim/{stage}/processed-genbank/{strains}.gbk",
    output:
        auto_gbk=temp("data/interim/{stage}/automlst_wrapper/{name}/{strains}.gbk"),
    conda:
        "../envs/automlst_wrapper.yaml"
    log:
        "logs/{stage}/automlst_wrapper/prep_automlst_gbk/prep_automlst_gbk-{name}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/features/prep_automlst.py {input.gbk} {output.auto_gbk} {wildcards.strains} 2>> {log}
        """


rule automlst_wrapper:
    input:
        # gbk=lambda wildcards: get_automlst_inputs(wildcards.name, filter_samples_qc(wildcards, get_samples_df())),
        gbk=fexpand("data/interim/{{stage}}/automlst_wrapper/{{name}}/{accession}.gbk", accession=RULE_FUNCTIONS["automlst_wrapper"]["accessions"]),
        reduced_core="resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    output:
        tree="data/interim/{stage}/automlst_wrapper/{name}/raxmlpart.txt.treefile",
        gbk_to_sql=temp("data/interim/{stage}/automlst_wrapper/{name}/genbank_to_sql"),
    log:
        "logs/{stage}/automlst_wrapper/automlst_wrapper/automlst_wrapper-{name}.log",
    conda:
        "../envs/automlst_wrapper.yaml"
    threads: 8
    resources:
        #TODO Fix
        tmpdir="data/interim/{stage}/automlst_wrapper/tmpdir/",
    shell:
        """
        mkdir -p "data/interim/{wildcards.stage}/automlst_wrapper/{wildcards.name}/singles" 2>> {log}
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py data/interim/{wildcards.stage}/automlst_wrapper/{wildcards.name} {threads} &>> {log}
        """


rule automlst_wrapper_out:
    input:
        tree="data/interim/{stage}/automlst_wrapper/{name}/raxmlpart.txt.treefile",
        # organism_info=lambda wildcards: expand("data/interim/prokka/{strains}/organism_info.txt",
        #             strains=[s for s in list(filter_samples_qc(wildcards, PEP_PROJECTS[wildcards.name].sample_table).index)],
        # ),
        # organism_info=lambda wildcards: expand("data/interim/prokka/{strains}/organism_info.txt",
        #             strains=list(get_samples_for_project_from_df(filter_samples_qc(wildcards, get_samples_df()), wildcards.name).index)),
        organism_info=fexpand("data/interim/all/prokka/{accession}/organism_info.txt", accession=RULE_FUNCTIONS["automlst_wrapper"]["accessions"]),
        # gtdb=lambda wildcards: expand("data/interim/gtdb/{strains}.json",
        #     name=wildcards.name,
        #     strains=[s for s in list(filter_samples_qc(wildcards, PEP_PROJECTS[wildcards.name].sample_table).index)])
        # gtdb=lambda wildcards: expand("data/interim/gtdb/{strains}.json", strains=list(get_samples_for_project_from_df(filter_samples_qc(wildcards, get_samples_df()), wildcards.name).index)),
        # gtdb=fexpand("data/interim/{stage}/gtdb/{accession}.json", stage=RULE_FUNCTIONS["automlst_wrapper"]["accessions_stage"], accession=RULE_FUNCTIONS["automlst_wrapper"]["accessions"]),
        gtdb="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
    output:
        genomes_tree="data/processed/{stage}/{name}/automlst_wrapper/df_genomes_tree.csv",
        mlst_genes="data/processed/{stage}/{name}/automlst_wrapper/df_mlst_genes.csv",
        final_tree="data/processed/{stage}/{name}/automlst_wrapper/final.newick",
        raxmlpart="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt",
        raxmlpart_bionj="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt.bionj",
        raxmlpart_contree="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt.contree",
        raxmlpart_iqtree="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt.iqtree",
        raxmlpart_log="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt.log",
        raxmlpart_mldist="data/processed/{stage}/{name}/automlst_wrapper/raxmlpart.txt.mldist",
        # raxmlpart_treefile="data/processed/{name}/automlst_wrapper/raxmlpart.txt.treefile",
    log:
        "logs/{stage}/automlst_wrapper/automlst_wrapper/automlst_wrapper_out-{name}.log",
    priority: 100
    params:
        automlst_interim="data/interim/{stage}/automlst_wrapper/{name}/",
        automlst_processed="data/processed/{stage}/{name}/automlst_wrapper/",
        prokka_interim="data/interim/{stage}/prokka",
        organism_info="data/interim/all/prokka/",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_phylo_tree.py {params.automlst_interim} {params.automlst_processed} {params.prokka_interim} {input.gtdb} {params.organism_info} 2>> {log}
        """
