#%
# final_output: data/processed/{name}/cblaster/genomes/
# description: Build diamond database of genomes for cblaster search.
# category: Genome Mining
# link:
# - https://github.com/gamcil/cblaster
# references:
# - 'Gilchrist, C., Booth, T. J., van Wersch, B., van Grieken, L., Medema, M. H.,
#   & Chooi, Y. (2021). cblaster: a remote search tool for rapid identification and
#   visualisation of homologous gene clusters (Version 1.3.9) [Computer software].
#   [https://doi.org/10.1101/2020.11.08.370601](https://doi.org/10.1101/2020.11.08.370601)'
# - 'Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND.
#   [Nat. Methods 12, 59–60 (2015)](https://www.nature.com/articles/nmeth.3176).'
#%
rule cblaster_genome_db:
    input:
        gbk=lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES, ext="gbk", path="processed-genbank"),
    output:
        folder_interim=directory("data/interim/cblaster/{name}/genomes/"),
        folder_processed=directory("data/processed/{name}/cblaster/genomes/"),
        sql="data/interim/cblaster/{name}/genomes/cblaster_genome_db.sqlite3",
        dmnd="data/interim/cblaster/{name}/genomes/cblaster_genome_db.dmnd",
        fasta="data/interim/cblaster/{name}/genomes/cblaster_genome_db.fasta",
        version="data/interim/cblaster/{name}/genomes/cblaster_genome_db.log",
    conda:
        "../envs/cblaster.yaml"
    log:
        "logs/cblaster/cblaster_db_genomes_{name}.log",
    threads: 8
    params:
        db_prefix="data/interim/cblaster/{name}/genomes/cblaster_genome_db",
        batch_size=50,
    shell:
        """
        cblaster --version >> {output.version}
        diamond --version >> {output.version}
        cat {output.version} >> {log}
        cblaster config --email dummy@cblaster.com 2>> {log}
        cblaster makedb --cpus {threads} -b {params.batch_size} -n {params.db_prefix} {input.gbk} 2>> {log}
        cp -r {output.folder_interim} {output.folder_processed} 2>> {log}
        """


rule cblaster_bgc_db:
    input:
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
    output:
        folder_interim=directory("data/interim/cblaster/{name}/bgcs/{version}/"),
        folder_processed=directory("data/processed/{name}/cblaster/bgcs/{version}/"),
        sql="data/interim/cblaster/{name}/bgcs/{version}/cblaster_bgc_db.sqlite3",
        dmnd="data/interim/cblaster/{name}/bgcs/{version}/cblaster_bgc_db.dmnd",
        fasta="data/interim/cblaster/{name}/bgcs/{version}/cblaster_bgc_db.fasta",
    conda:
        "../envs/cblaster.yaml"
    params:
        db_prefix="data/interim/cblaster/{name}/bgcs/{version}/cblaster_bgc_db",
        antismash_dir="data/interim/bgcs/{name}/{version}/*/*region*.gbk",
        batch_size=50,
    log:
        "logs/cblaster/cblaster_db_bgc_{name}_{version}.log",
    threads: 8
    shell:
        """
        cblaster config --email dummy@cblaster.com
        cblaster makedb --cpus {threads} -b {params.batch_size} -n {params.db_prefix} {params.antismash_dir} 2>> {log}
        cp -r {output.folder_interim} {output.folder_processed} 2>> {log}
        """
