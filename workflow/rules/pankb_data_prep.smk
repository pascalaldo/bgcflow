def get_genome_ids(name, df_samples):
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    return selection

def get_family_list():
    with open(checkpoints.pankb_family_list.get().output[0], 'r') as f:
        family_list = [family.strip() for family in f.readlines()]
    return family_list

def genome_list_dummy_input(wildcards):
    RULE_FUNCTIONS["pankb_data_prep"]["genomes"](wildcards)
    return []

rule pankb_genome_list:
    input:
        dummy=genome_list_dummy_input
    output:
        genomes="data/processed/{stage}/{name}/pankb/genomes.txt"
    log:
        "logs/{stage}/pankb_data_prep/pankb_genome_list_{name}.log"
    run:
        # genomes = get_genome_ids(wildcards.name, filter_samples_qc(wildcards, get_samples_df()))
        genomes = RULE_FUNCTIONS["pankb_data_prep"]["genomes"](wildcards)
        with open(output.genomes, "w") as f:
            f.writelines(f"{genome}\n" for genome in genomes)

rule pankb_select_and_merge:
    input:
        genomes=fexpand("data/processed/{{stage}}/{name}/pankb/genomes.txt", name=RULE_FUNCTIONS["pankb_data_prep"]["projects"]),
        csv=fexpand("data/processed/{{stage}}/{name}/{{directory}}/{{filename}}.csv", name=RULE_FUNCTIONS["pankb_data_prep"]["projects"]),
    output: "data/processed/{stage}/pankb/merged/{directory}/{filename}.csv",
    run:
        import pandas as pd

        df = None
        for genomes_file, csv_file in zip(input.genomes, input.csv):
            add_df = pd.read_csv(csv_file, index_col=0, low_memory=False)
            if add_df.index.name == "genome_id":
                with open(genomes_file, 'r') as f:
                    genome_list = [genome.strip() for genome in f.readlines()]
                add_df = add_df.loc[genome_list, :]
            df = pd.concat([df, add_df])
        df.to_csv(output[0])

checkpoint pankb_family_list:
    input:
        gtdb_merged="data/processed/{stage}/pankb/merged/tables/df_gtdb_meta.csv"
    output:
        families="data/processed/{stage}/pankb/families.txt"
    run:
        import pandas as pd

        df = pd.read_csv(input.gtdb_merged, index_col=0, low_memory=False)
        families = [family.replace('f__', '') for family in list(set(df["Family"].tolist()))]
        with open(output.families, "w") as f:
            f.writelines(f"{family}\n" for family in families)

rule pankb_isosource:
    input:
        genomes="data/processed/{stage}/{name}/pankb/genomes.txt"
    output:
        isosource="data/processed/{stage}/{name}/pankb/source_info/df_ncbi_isolation_src.csv"
    log:
        "logs/{stage}/pankb_data_prep/pankb_isosource_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep isosource \
            -g {input.genomes} \
            -o {output.isosource} > {log} 2>&1
        """

rule pankb_species_summary:
    input:
        gtdb_meta="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
    output:
        csv="data/processed/{stage}/{name}/pankb/summary.csv",
        json="data/processed/{stage}/pankb/pankb/web_data/species/{name}/info_panel.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_species_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep species {wildcards.name} \
            --gtdb_meta {input.gtdb_meta} \
            --summary {input.summary_v2} \
            --gp_binary {input.gp_binary} \
            --output_json {output.json} \
            -o {output.csv} > {log} 2>&1
        """

rule pankb_family_summary:
    input:
        gtdb_merged="data/processed/{stage}/pankb/merged/tables/df_gtdb_meta.csv",
        seqfu_stats="data/processed/{stage}/pankb/merged/tables/df_seqfu_stats.csv",
    output:
        summary="data/processed/{stage}/pankb/family/{family}/summary.csv"
    log:
        "logs/{stage}/pankb_data_prep/pankb_family_{family}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep family {wildcards.family} \
            --gtdb_meta {input.gtdb_merged} \
            --seqfu {input.seqfu_stats} \
            -o {output.summary} > {log} 2>&1
        """

rule pankb_full_summary:
    input:
        gtdb_merged="data/processed/{stage}/pankb/merged/tables/df_gtdb_meta.csv",
        seqfu_stats="data/processed/{stage}/pankb/merged/tables/df_seqfu_stats.csv",
    output:
        summary="data/processed/{stage}/pankb/pankb/full_summary.csv"
    log:
        "logs/{stage}/pankb_data_prep/pankb_full_summary.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep full \
            --gtdb_meta {input.gtdb_merged} \
            --seqfu {input.seqfu_stats} \
            -o {output.summary} > {log} 2>&1
        """

rule pankb_pangenome_summary:
    input:
        species_summary="data/processed/{stage}/pankb/pankb/summary.csv",
    output:
        species_list="data/processed/{stage}/pankb/pankb/web_data/species_list.json",
        genome_count="data/processed/{stage}/pankb/pankb/organism_genome_count.json",
        gene_count="data/processed/{stage}/pankb/pankb/organism_gene_count.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_all.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep all \
            --species_summary {input.species_summary} \
            --output_genome_count {output.genome_count} \
            --output_gene_count {output.gene_count} \
            --output_json {output.species_list} > {log} 2>&1
        """

rule pankb_mash:
    input:
        genomes="data/processed/{stage}/{name}/pankb/genomes.txt",
        gtdb_meta="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
        mash_file="data/processed/{stage}/{name}/mash/df_mash.csv",
    output:
        mash_list="data/processed/{stage}/{name}/pankb/mash_list.csv",
    log:
        "logs/{stage}/pankb_data_prep/pankb_mash_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep mash \
            --genomes {input.genomes} \
            --gtdb_meta {input.gtdb_meta} \
            --mash_file {input.mash_file} \
            -o {output.mash_list} > {log} 2>&1
        """

rule pankb_eggnog:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        eggnog_table="data/processed/{stage}/{name}/eggnog_roary/emapper.annotations",
        reference="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
    output:
        eggnog_summary="data/processed/{stage}/{name}/pankb/df_pangene_eggnog_summary.csv"
    log:
        "logs/{stage}/pankb_data_prep/pankb_eggnog_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep eggnog \
            --summary {input.summary_v2} \
            --gp_binary {input.gp_binary} \
            --eggnog_table {input.eggnog_table} \
            --reference {input.reference} \
            -o {output.eggnog_summary} > {log} 2>&1
        """
rule pankb_heatmap:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        gp_locustag="data/interim/{stage}/roary/{name}/df_gene_presence_locustag.csv",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{stage}/{name}/pankb/df_pangene_eggnog_summary.csv",
        mash_list="data/processed/{stage}/{name}/pankb/mash_list.csv",
        isosource="data/processed/{stage}/{name}/pankb/source_info/df_ncbi_isolation_src.csv",
        # species_info="data/processed/{name}/tables/df_ncbi_meta.csv",
    output:
        heatmap_target="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_target.json.gz",
        heatmap_core="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_core.json.gz",
        heatmap_accessory="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_accessory.json.gz",
        heatmap_1_15="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_1_15.json.gz",
        heatmap_above_1="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_above_1.json.gz",
        heatmap_only_1="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_only_1.json.gz",
        source_info_core="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_core.json.gz",
        source_info_accessory="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_accessory.json.gz",
        source_info_1_15="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_1_15.json.gz",
        source_info_above_1="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_above_1.json.gz",
        source_info_only_1="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_only_1.json.gz",
    params:
        heatmap_base="data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_",
        source_info_base="data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_",
    log:
        "logs/{stage}/pankb_data_prep/pankb_heatmap_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep heatmap \
            --summary {input.summary_v2} \
            --gp_binary {input.gp_binary} \
            --gp_locustag {input.gp_locustag} \
            --eggnog_summary {input.eggnog_summary} \
            --mash_list {input.mash_list} \
            --isosource {input.isosource} \
            --sourceinfo_base {params.source_info_base} \
            --heatmap_base {params.heatmap_base} > {log} 2>&1
        """

rule pankb_cog:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        eggnog_summary="data/processed/{stage}/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        all_cog="data/processed/{stage}/pankb/pankb/web_data/species/{name}/All.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_cog_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep cog \
            --gp_binary {input.gp_binary} \
            --eggnog_summary {input.eggnog_summary} \
            --output_json {output.all_cog} > {log} 2>&1
        """

rule pankb_cog_distribution:
    input:
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{stage}/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        cog_distribution="data/processed/{stage}/pankb/pankb/web_data/species/{name}/COG_distribution.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_distribution_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep distribution \
            --summary {input.summary_v2} \
            --eggnog_summary {input.eggnog_summary} \
            --output_json {output.cog_distribution} > {log} 2>&1
        """

rule pankb_heaps:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
    output:
        gene_freq="data/processed/{stage}/pankb/pankb/web_data/species/{name}/gene_freq.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_heaps_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep heaps \
            --gp_binary {input.gp_binary} \
            --output_json {output.gene_freq} > {log} 2>&1
        """

rule pankb_locustag:
    input:
        gp_locustag="data/interim/{stage}/roary/{name}/df_gene_presence_locustag.csv",
        all_locustag="data/interim/{stage}/alleleome/{name}/all_locustag.csv",
    output:
        locustag=directory("data/processed/{stage}/pankb/pankb/web_data/species/{name}/gene_locustag/"),
    log:
        "logs/{stage}/pankb_data_prep/pankb_locustag_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep locustag \
            --gp_locustag {input.gp_locustag} \
            --all_locustag {input.all_locustag} \
            -o {output.locustag} > {log} 2>&1
        """

rule pankb_genome_page:
    input:
        full_summary="data/processed/{stage}/pankb/pankb/full_summary.csv",
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        isosource="data/processed/{stage}/{name}/pankb/source_info/df_ncbi_isolation_src.csv",
        # species_info="data/processed/{stage}/{name}/tables/df_ncbi_meta.csv",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{stage}/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        genome_page_dir=directory("data/processed/{stage}/pankb/pankb/web_data/species/{name}/genome_page/")
    log:
        "logs/{stage}/pankb_data_prep/pankb_genome_page_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep genome \
            --summary {input.summary_v2} \
            --gp_binary {input.gp_binary} \
            --species_summary {input.full_summary} \
            --isosource {input.isosource} \
            --eggnog_summary {input.eggnog_summary} \
            -o {output.genome_page_dir} > {log} 2>&1
        """

rule pankb_landing_page:
    input:
        species_summary="data/processed/{stage}/pankb/pankb/summary.csv",
    output:
        pankb_dimension="data/processed/{stage}/pankb/pankb/pankb_dimension.json",
        species_genome_gene="data/processed/{stage}/pankb/pankb/species_genome_gene.json",
    log:
        "logs/{stage}/pankb_data_prep/pankb_landing_page.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep landing \
            --species_summary {input.species_summary} \
            --output_pankb_dimension {output.pankb_dimension} \
            --output_json {output.species_genome_gene} > {log} 2>&1
        """

rule pankb_copy_phylogenetic_tree:
    input:
        "data/processed/{stage}/{name}/automlst_wrapper/final.newick"
    output:
        "data/processed/{stage}/pankb/pankb/web_data/species/{name}/phylogenetic_tree.newick"
    shell:
        """
        cp {input} {output}
        """

rule pankb_gzip:
    input: "data/processed/{stage}/{name}/pankb/{filename}",
    output: "data/processed/{stage}/{name}/pankb/{filename}.gz",
    shell:
        """
        gzip -c {input} > {output}
        """

rule pankb_copy_panalleleome:
    input:
        step_line="data/processed/{stage}/{name}/alleleome/Pan/step_line.json",
        dn_ds="data/processed/{stage}/{name}/alleleome/Pan/dn_ds.json",
    output:
        step_line="data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/step_line.json",
        dn_ds="data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/dn_ds.json",
        gene_data=directory("data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/gene_data/"),
    params:
        gene_data_in="data/processed/{stage}/{name}/alleleome/gene_data/",
    shell:
        """
        cp {input.step_line} {output.step_line}
        cp {input.dn_ds} {output.dn_ds}
        ln -sr {params.gene_data_in} {output.gene_data}
        """

rule pankb_all:
    input:
        fexpand(
            [
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/genome_page/",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/gene_locustag/",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/gene_freq.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/COG_distribution.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/All.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/heatmap_target.json.gz",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/source_info_core.json.gz",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/info_panel.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/step_line.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/dn_ds.json",
                "data/processed/{stage}/pankb/pankb/web_data/species/{name}/panalleleome/gene_data/",
                "data/processed/{stage}/pankb/pankb/web_data/species_list.json",
            ],
            stage=RULE_FUNCTIONS["pankb_data_prep"]["stages"],
            name=RULE_FUNCTIONS["pankb_data_prep"]["projects"],
        ),
