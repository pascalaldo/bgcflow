def get_genome_ids(name, df_samples):
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    return selection

def get_family_list():
    with open(checkpoints.pankb_family_list.get().output[0], 'r') as f:
        family_list = [family.strip() for family in f.readlines()]
    return family_list

rule pankb_genome_list:
    input:
        rules.qc.output
    output:
        genomes="data/processed/{name}/pankb/genomes.txt"
    log:
        "logs/pankb_data_prep/pankb_genome_list_{name}.log"
    run:
        genomes = get_genome_ids(wildcards.name, filter_samples_qc(wildcards, DF_SAMPLES))
        with open(output.genomes, "w") as f:
            f.writelines(f"{genome}\n" for genome in genomes)

rule pankb_select_and_merge:
    input:
        genomes=expand("data/processed/{name}/pankb/genomes.txt", name=PROJECT_IDS),
        csv=expand("data/processed/{name}/{{directory}}/{{filename}}.csv", name=PROJECT_IDS),
    output: "data/processed/pankb/{directory}/{filename}.csv",
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
        gtdb_merged="data/processed/pankb/tables/df_gtdb_meta.csv"
    output:
        families="data/processed/pankb/families.txt"
    run:
        import pandas as pd

        df = pd.read_csv(input.gtdb_merged, index_col=0, low_memory=False)
        families = [family.replace('f__', '') for family in list(set(df["Family"].tolist()))]
        with open(output.families, "w") as f:
            f.writelines(f"{family}\n" for family in families)

rule pankb_isosource:
    input:
        genomes="data/processed/{name}/pankb/genomes.txt"
    output:
        isosource="data/processed/{name}/pankb/source_info/df_ncbi_isolation_src.csv"
    log:
        "logs/pankb_data_prep/pankb_isosource_{name}.log"
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
        gtdb_meta="data/processed/{name}/tables/df_gtdb_meta.csv",
        summary_v2="data/interim/alleleome/{name}/pangene_v2.csv",
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
    output:
        csv="data/processed/{name}/pankb/summary.csv",
        json="data/processed/{name}/pankb/summary.json",
    log:
        "logs/pankb_data_prep/pankb_species_{name}.log"
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
        gtdb_merged="data/processed/pankb/tables/df_gtdb_meta.csv",
    output:
        summary="data/processed/pankb/family/{family}/summary.csv"
    log:
        "logs/pankb_data_prep/pankb_family_{family}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep family {wildcards.family} \
            --gtdb_meta {input.gtdb_merged} \
            -o {output.summary} > {log} 2>&1
        """

rule pankb_full_summary:
    input:
        gtdb_merged="data/processed/pankb/tables/df_gtdb_meta.csv",
    output:
        summary="data/processed/pankb/pankb/full_summary.csv"
    log:
        "logs/pankb_data_prep/pankb_full_summary.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep full \
            --gtdb_meta {input.gtdb_merged} \
            -o {output.summary} > {log} 2>&1
        """

rule pankb_pangenome_summary:
    input:
        species_summary="data/processed/pankb/pankb/summary.csv",
    output:
        species_list="data/processed/pankb/pankb/species_list.json",
        genome_count="data/processed/pankb/pankb/organism_genome_count.json",
        gene_count="data/processed/pankb/pankb/organism_gene_count.json",
    log:
        "logs/pankb_data_prep/pankb_all.log"
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
        genomes="data/processed/{name}/pankb/genomes.txt",
        gtdb_meta="data/processed/{name}/tables/df_gtdb_meta.csv",
        mash_file="data/processed/{name}/mash/df_mash.csv",
    output:
        mash_list="data/processed/{name}/pankb/mash_list.csv",
    log:
        "logs/pankb_data_prep/pankb_mash_{name}.log"
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
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        summary_v2="data/interim/alleleome/{name}/pangene_v2.csv",
        eggnog_table="data/processed/{name}/eggnog_roary/emapper.annotations",
        reference="data/interim/roary/{name}/pan_genome_reference.fa",
    output:
        eggnog_summary="data/processed/{name}/pankb/df_pangene_eggnog_summary.csv"
    log:
        "logs/pankb_data_prep/pankb_eggnog_{name}.log"
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
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gp_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        summary_v2="data/interim/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{name}/pankb/df_pangene_eggnog_summary.csv",
        mash_list="data/processed/{name}/pankb/mash_list.csv",
        isosource="data/processed/{name}/pankb/source_info/df_ncbi_isolation_src.csv",
        species_info="data/processed/{name}/tables/df_ncbi_meta.csv",
    output:
        heatmap="data/processed/{name}/pankb/heatmap_target.json",
    log:
        "logs/pankb_data_prep/pankb_heatmap_{name}.log"
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
            --species_info {input.species_info} \
            --output_json {output.heatmap} > {log} 2>&1
        """

rule pankb_cog:
    input:
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        eggnog_summary="data/processed/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        all_cog="data/processed/{name}/pankb/All.json",
    log:
        "logs/pankb_data_prep/pankb_cog_{name}.log"
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
        summary_v2="data/interim/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        cog_distribution="data/processed/{name}/pankb/COG_distribution.json",
    log:
        "logs/pankb_data_prep/pankb_distribution_{name}.log"
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
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
    output:
        gene_freq="data/processed/{name}/pankb/gene_freq.json",
    log:
        "logs/pankb_data_prep/pankb_heaps_{name}.log"
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
        gp_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        all_locustag="data/interim/alleleome/{name}/all_locustag.csv",
    output:
        locustag=directory("data/processed/{name}/pankb/gene_locustag/"),
    log:
        "logs/pankb_data_prep/pankb_locustag_{name}.log"
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
        full_summary="data/processed/pankb/pankb/full_summary.csv",
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        isosource="data/processed/{name}/pankb/source_info/df_ncbi_isolation_src.csv",
        species_info="data/processed/{name}/tables/df_ncbi_meta.csv",
        summary_v2="data/interim/alleleome/{name}/pangene_v2.csv",
        eggnog_summary="data/processed/{name}/pankb/df_pangene_eggnog_summary.csv",
    output:
        genome_page_dir=directory("data/processed/{name}/pankb/genome_page/")
    log:
        "logs/pankb_data_prep/pankb_genome_page_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep genome \
            --summary {input.summary_v2} \
            --gp_binary {input.gp_binary} \
            --species_summary {input.full_summary} \
            --isosource {input.isosource} \
            --species_info {input.species_info} \
            --eggnog_summary {input.eggnog_summary} \
            -o {output.genome_page_dir} > {log} 2>&1
        """

rule pankb_landing_page:
    input:
        species_summary="data/processed/pankb/pankb/summary.csv",
    output:
        pankb_dimension="data/processed/pankb/pankb/pankb_dimension.json",
        species_genome_gene="data/processed/pankb/pankb/species_genome_gene.json",
    log:
        "logs/pankb_data_prep/pankb_landing_page.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_data_prep landing \
            --species_summary {input.species_summary} \
            --output_pankb_dimension {output.pankb_dimension} \
            --output_json {output.species_genome_gene} > {log} 2>&1
        """

rule pankb_gzip:
    input: "data/processed/{name}/pankb/{filename}",
    output: "data/processed/{name}/pankb/{filename}.gz",
    shell:
        """
        gzip -c {input} > {output}
        """