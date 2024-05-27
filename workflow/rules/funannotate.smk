if len(CUSTOM_FNA) > 0:
    rule copy_custom_fasta:
        input:
            fna = lambda wildcards: DF_SAMPLES.loc[wildcards, "input_file"]
        output:
            fna = "data/interim/fasta/{custom_fna}.fna"
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/prokka/copy_custom_fasta/copy_custom_fasta-{custom_fna}.log"
        shell:
            """
            if [[ {input} == *.fna || {input} == *.fasta || {input} == *.fa ]]
            then
                cp {input} {output} 2>> {log}
            else
                echo "ERROR: Wrong Extension:" {input} >> {log}
                exit 1
            fi
            """

if len(STRAINS_FNA) > 0:
    rule extract_meta_prokka:
        input:
            fna = "data/interim/fasta/{strains_fna}.fna",
        output:
            org_info = "data/interim/prokka/{strains_fna}/organism_info.txt",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/prokka/extract_meta_prokka/extract_meta_prokka-{strains_fna}.log"
        params:
            samples_path = bgcflow_util_dir / "samples.csv",
        shell:
            """
            python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.strains_fna} \
                "{params.samples_path}" data/interim/assembly_report/ data/interim/prokka/ 2>> {log}
            """
    try:
        if os.path.isfile(config["resources_path"]["pfam_for_prokka"]):
            prokka_use_pfam = f'--hmms {config["resources_path"]["pfam_for_prokka"]}'
        else:
            prokka_use_pfam = ""
    except KeyError:
        prokka_use_pfam = ""

    rule funannotate_setup:
        output:
            funannotate_dbs=directory("resources/funannotate/"),
        conda:
            "../envs/funannotate.yaml"
        log: "logs/funannotate/funannotate_setup/funannotate-setup.log"
        shell:
            """
            funannotate setup --install all -d {output.funannotate_dbs}
            """

    rule funannotate_clean:
        input:
            fna = "data/interim/fasta/{strains_fna}.fna",
        output:
            clean_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.fa",
        conda:
            "../envs/funannotate.yaml"
        log: "logs/funannotate/funannotate_clean/funannotate-{strains_fna}.log"
        shell:
            """
            cp {input.fna} {output.clean_fa}
            """
        # shell:
        #     """
        #     funannotate clean -i {input.fna} --minlen 1000 -o {output.clean_fa} 2> {log}
        #     """

    rule funannotate_sort:
        input:
            clean_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.fa",
        output:
            sorted_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.sorted.fa",
        conda:
            "../envs/funannotate.yaml"
        log: "logs/funannotate/funannotate_sort/funannotate-{strains_fna}.log"
        shell:
            """
            funannotate sort -i {input.clean_fa} -o {output.sorted_fa} 2> {log}
            """
    
    rule funannotate_mask:
        input:
            sorted_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.sorted.fa",
        output:
            masked_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.sorted.masked.fa",
        conda:
            "../envs/funannotate.yaml"
        log: "logs/funannotate/funannotate_mask/funannotate-{strains_fna}.log"
        threads: 12
        shell:
            """
            funannotate mask -i {input.sorted_fa} --cpus {threads} -o {output.masked_fa} 2> {log}
            """

    rule funannotate_predict:
        input:
            funannotate_dbs="resources/funannotate/",
            masked_fa = "data/interim/funannotate/{strains_fna}/{strains_fna}.clean.sorted.masked.fa",
            org_info = "data/interim/prokka/{strains_fna}/organism_info.txt",
        output:
            gff = "data/interim/prokka/{strains_fna}/{strains_fna}.gff",
            faa = "data/interim/prokka/{strains_fna}/{strains_fna}.faa",
            gbk = "data/interim/prokka/{strains_fna}/{strains_fna}.gbk",
            txt = "data/interim/prokka/{strains_fna}/{strains_fna}.txt",
            tsv = "data/interim/prokka/{strains_fna}/{strains_fna}.tsv",
            fna = temp("data/interim/prokka/{strains_fna}/{strains_fna}.fna"),
            sqn = temp("data/interim/prokka/{strains_fna}/{strains_fna}.sqn"),
            fsa = temp("data/interim/prokka/{strains_fna}/{strains_fna}.fsa"),
            tbl = temp("data/interim/prokka/{strains_fna}/{strains_fna}.tbl"),
        conda:
            "../envs/funannotate.yaml"
        log: "logs/funannotate/funannotate_predict/funannotate-{strains_fna}.log"
        threads: 4
        shell:
            """
            FUNANNOTATE_DB="{input.funannotate_dbs}" funannotate predict \
                -i {input.masked_fa} \
                --species "`cut -d "," -f 2 {input.org_info}`" \
                --strain "`cut -d "," -f 3 {input.org_info}`" \
                --cpus {threads} \
                -o "data/interim/prokka/{wildcards.strains_fna}/" 2> {log}
            """

    rule format_gbk:
        input:
            gbk_prokka = "data/interim/prokka/{strains_fna}/{strains_fna}.gbk",
            gtdb_json = "data/interim/gtdb/{strains_fna}.json",
        output:
            gbk_processed = "data/interim/processed-genbank/{strains_fna}.gbk"
        conda:
            "../envs/bgc_analytics.yaml"
        params:
            version = __version__,
        log: "logs/prokka/format_gbk/format_gbk-{strains_fna}.log"
        shell:
            """
            python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} \
                {params.version} {input.gtdb_json} {wildcards.strains_fna} {output.gbk_processed} 2> {log}
            """

    rule copy_prokka_gbk:
        input:
            gbk = "data/interim/processed-genbank/{strains_fna}.gbk",
            summary = "data/interim/prokka/{strains_fna}/{strains_fna}.txt",
            tsv = "data/interim/prokka/{strains_fna}/{strains_fna}.tsv"
        output:
            gbk = report("data/processed/{name}/genbank/{strains_fna}.gbk", \
                caption="../report/file-genbank.rst", category="{name}", subcategory="Annotated Genbanks"),
            summary = "data/processed/{name}/genbank/{strains_fna}.txt",
            tsv = "data/processed/{name}/genbank/{strains_fna}.tsv",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/prokka/copy_prokka_gbk/copy_prokka_gbk_-{strains_fna}-{name}.log"
        shell:
            """
            cp {input.gbk} {output.gbk} 2>> {log}
            cp {input.summary} {output.summary} 2>> {log}
            cp {input.tsv} {output.tsv} 2>> {log}
            """
