#%
# description: Annotate using Prokka.
# category: Functional Annotation
# link:
# - https://github.com/tseemann/prokka
# references:
# - 'Seemann T. Prokka: rapid prokaryotic genome annotation. [Bioinformatics 2014 Jul
#   15;30(14):2068-9. PMID:24642063](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517?login=false)'
#%
try:
    if Path(config["resources_path"]["RNAmmer"]).is_file():
        prokka_params_rna = "--rnammer"
        rule rnammer_setup:
            output:
                "resources/rnammer_test.txt"
            priority: 50
            conda:
                "../envs/prokka.yaml"
            log: "logs/prokka/rnammer_setup.log"
            shell:
                """
                ln -s $PWD/resources/RNAmmer/rnammer $CONDA_PREFIX/bin/rnammer 2>> {log}
                rnammer -S bac -m lsu,ssu,tsu -gff - resources/RNAmmer/example/ecoli.fsa >> {output}
                """
    else:
        prokka_params_rna = ""
        pass
except KeyError:
    prokka_params_rna = ""
    pass

rule prokka_db_setup:
    input:
        table = "resources/prokka_db/{prokka_db}.json"
    output:
        refgbff = "resources/prokka_db/{prokka_db}.gbff",
        ncbi_tempdir = temp(directory("resources/prokka_db/{prokka_db}")),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/prokka/prokka_db_setup/prokka_db_setup-{prokka_db}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_prokka_db.py {wildcards.prokka_db} \
            {input.table} {output.refgbff} 2>> {log}
        """
rule extract_meta_prokka:
    input:
        fna = "data/interim/all/fasta/{strains_fna}.fna",
        samples_path = bgcflow_util_dir / "samples.csv",
        assembly_report= "data/interim/all/assembly_report/{strains_fna}.json",
    output:
        org_info = "data/interim/all/prokka/{strains_fna}/organism_info.txt",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/all/prokka/extract_meta_prokka/extract_meta_prokka-{strains_fna}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.strains_fna} \
            "{input.samples_path}" {input.assembly_report} data/interim/all/prokka/ 2>> {log}
        """
try:
    if os.path.isfile(config["resources_path"]["pfam_for_prokka"]):
        prokka_use_pfam = f'--hmms {config["resources_path"]["pfam_for_prokka"]}'
    else:
        prokka_use_pfam = ""
except KeyError:
    prokka_use_pfam = ""

rule prokka:
    input:
        fna = "data/interim/all/fasta/{strains_fna}.fna",
        org_info = "data/interim/all/prokka/{strains_fna}/organism_info.txt",
        # refgbff = lambda wildcards: get_prokka_refdb(wildcards, "file", RULE_FUNCTIONS["prokka"]["samples"](wildcards), PROKKA_DB_MAP)
        refgbff=[]
    output:
        gff = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.gff",
        faa = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.faa",
        gbk = temp("data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.gbk"), # Temporary file now, use processed-genbank files in other rules
        txt = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.txt",
        tsv = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.tsv",
        fna = temp("data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.fna"),
        sqn = temp("data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.sqn"),
        fsa = temp("data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.fsa"),
        tbl = temp("data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.tbl"),
    conda:
        "../envs/prokka.yaml"
    log: "logs/{stage}/prokka/prokka/prokka-{strains_fna}.log"
    params:
        increment = 10,
        evalue = "1e-05",
        rna_detection = prokka_params_rna,
        # refgbff = lambda wildcards: get_prokka_refdb(wildcards, "params", RULE_FUNCTIONS["prokka"][wildcards.stage]["samples"](), PROKKA_DB_MAP),
        refgbff = "",
        use_pfam = prokka_use_pfam,
        kingdom = KINGDOM.capitalize(),
    threads: 4
    shell:
        """
        prokka --outdir data/interim/{wildcards.stage}/prokka/{wildcards.strains_fna} --force \
            --kingdom {params.kingdom} \
            {params.refgbff} --prefix {wildcards.strains_fna} --genus "`cut -d "," -f 1 {input.org_info}`" \
            --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" \
            --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} \
            {params.use_pfam} --evalue {params.evalue} {input.fna} &> {log}
        """

rule format_gbk:
    input:
        gbk_prokka = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.gbk",
        gtdb_json = "data/interim/{stage}/gtdb/{strains_fna}.json",
    output:
        gbk_processed = "data/interim/{stage}/processed-genbank/{strains_fna}.gbk"
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        version = __version__,
    log: "logs/{stage}/prokka/format_gbk/format_gbk-{strains_fna}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} \
            {params.version} {input.gtdb_json} {wildcards.strains_fna} {output.gbk_processed} 2> {log}
        """

rule copy_prokka_gbk:
    input:
        gbk = "data/interim/{stage}/processed-genbank/{strains_fna}.gbk",
        summary = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.txt",
        tsv = "data/interim/{stage}/prokka/{strains_fna}/{strains_fna}.tsv",
    output:
        gbk = report("data/processed/{stage}/{name}/genbank/{strains_fna}.gbk", \
            caption="../report/file-genbank.rst", category="{name}", subcategory="Annotated Genbanks"),
        summary = "data/processed/{stage}/{name}/genbank/{strains_fna}.txt",
        tsv = "data/processed/{stage}/{name}/genbank/{strains_fna}.tsv",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/{stage}/prokka/copy_prokka_gbk/copy_prokka_gbk_-{strains_fna}-{name}.log"
    shell:
        """
        cp {input.gbk} {output.gbk} 2>> {log}
        cp {input.summary} {output.summary} 2>> {log}
        cp {input.tsv} {output.tsv} 2>> {log}
        """
