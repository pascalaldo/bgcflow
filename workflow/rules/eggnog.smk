#%
# final_output: data/interim/eggnog/{strains}/
# description: Annotate samples with eggNOG database (http://eggnog5.embl.de)
# category: Functional Annotation
# link:
# - https://github.com/eggnogdb/eggnog-mapper
# references:
# - 'eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction
#     at the metagenomic scale. Carlos P. Cantalapiedra, Ana Hernandez-Plaza, Ivica
#     Letunic, Peer Bork, Jaime Huerta-Cepas. 2021. [Molecular Biology and Evolution,
#     msab293](https://doi.org/10.1093/molbev/msab293)'
# - 'eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology
#     resource based on 5090 organisms and 2502 viruses. Jaime Huerta-Cepas, Damian
#     Szklarczyk, Davide Heller, Ana Hernández-Plaza, Sofia K Forslund, Helen Cook,
#     Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars J Jensen, Christian von Mering,
#     Peer Bork [Nucleic Acids Res. 2019 Jan 8; 47(Database issue): D309–D314. doi: 10.1093/nar/gky1085](https://academic.oup.com/nar/article/47/D1/D309/5173662)'
#%

rule install_eggnog:
    output:
        eggnog_db=directory("resources/eggnog_db"),
        dmnd="resources/eggnog_db/bacteria.dmnd",
    conda:
        "../envs/eggnog.yaml"
    log:
        "logs/eggnog/install_eggnog.log",
    shell:
        """
        download_eggnog_data.py --data_dir {output.eggnog_db} -y &>> {log}
        create_dbs.py -m diamond --dbname bacteria --taxa Bacteria --data_dir {output.eggnog_db} -y &>> {log}
        """


rule eggnog:
    input:
        faa="data/interim/prokka/{strains}/{strains}.faa",
        eggnog_db="resources/eggnog_db",
        dmnd="resources/eggnog_db/bacteria.dmnd",
    output:
        eggnog_dir=directory("data/interim/eggnog/{strains}/"),
        tempdir=temp(directory("data/interim/eggnog/tmp/{strains}"))
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log:
        "logs/eggnog/eggnog/eggnog-{strains}.log",
    shell:
        """
        mkdir -p {output.eggnog_dir}
        mkdir -p {output.tempdir}
        emapper.py -i {input.faa} --decorate_gff "yes" --excel --cpu {threads} -o {wildcards.strains} --output_dir {output.eggnog_dir} --data_dir {input.eggnog_db} --temp_dir {output.tempdir} &>> {log}
        """
