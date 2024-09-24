#%
# final_output: data/processed/{name}/tables/gtdbtk.bac120.summary.tsv
# description: Taxonomic placement with GTDB-Tk
# category: Phylogenomic Placement
# link:
# - https://github.com/Ecogenomics/GTDBTk
# references:
# - 'Chaumeil PA, et al. 2019. GTDB-Tk: A toolkit to classify genomes with the Genome
#   Taxonomy Database. Bioinformatics, btz848.'
# - Parks DH, et al. 2020. A complete domain-to-species taxonomy for Bacteria and
#   Archaea. Nature Biotechnology, [https://doi.org/10.1038/s41587-020-0501-8]([https://doi.org/10.1038/s41587-020-0501-8).
# - Parks DH, et al. 2018. A standardized bacterial taxonomy based on genome phylogeny
#   substantially revises the tree of life. Nature Biotechnology, [http://dx.doi.org/10.1038/nbt.4229](http://dx.doi.org/10.1038/nbt.4229).
#%
# Read release version from config
try:
    gtdb_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdb_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdb_release = "220.0"
    gtdb_release_version = "r220"

if "." in str(gtdb_release):
    gtdb_release_major, gtdb_release_minor = str(gtdb_release).split(".")
else:
    gtdb_release_major = gtdb_release
    gtdb_release_minor = "0"
gtdb_extra_path = ""
if int(gtdb_release_major) >= 220:
    gtdb_extra_path = "gtdbtk_package/full_package/"

# Decide to use ani screen or not
try:
    if config["rule_parameters"]["gtdbtk"]["ani_screen"]:
        ani_screen = f"--mash_db resources/gtdb-tk-{gtdb_release_version}.msh"
    else:
        ani_screen = "--skip_ani_screen"
except KeyError:
    ani_screen = "--skip_ani_screen"

rule install_gtdbtk:
    output:
        gtdbtk=directory("resources/gtdbtk/"),
    priority: 40
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk/gtdbtk-install_gtdbtk.log",
    params:
        release_major=gtdb_release_major,
        release_minor=gtdb_release_minor,
        release_version=gtdb_release_version,
        extra_path=gtdb_extra_path,
    shell:
        """
        (cd resources && wget https://data.gtdb.ecogenomic.org/releases/release{params.release_major}/{params.release_major}.{params.release_minor}/auxillary_files/{params.extra_path}gtdbtk_{params.release_version}_data.tar.gz -nc) 2>> {log}
        (cd resources && mkdir -p gtdbtk && tar -xvzf gtdbtk_{params.release_version}_data.tar.gz -C "gtdbtk" --strip 1 && rm gtdbtk_{params.release_version}_data.tar.gz) &>> {log}
        """

checkpoint prepare_gtdbtk_input:
    input:
        gtdb_meta="data/interim/{stage}/gtdb/{name}/tables/df_gtdb_meta.csv",
        json_list=lambda wildcards: expand("data/interim/all/assembly_report/{accession}.json", accession=RULE_FUNCTIONS["gtdbtk_simple"][wildcards.stage]["accessions"](wildcards.name)),
        fna=lambda wildcards: expand("data/interim/all/fasta/{accession}.fna", accession=RULE_FUNCTIONS["gtdbtk_simple"][wildcards.stage]["accessions"](wildcards.name)),
    output:
        fnadir=directory("data/interim/{stage}/gtdbtk/{name}/fasta/"),
        fnalist="data/interim/{stage}/gtdbtk/{name}/fasta_list.txt",
    log:
        "logs/{stage}/gtdbtk/prepare_gtdbtk_input/{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        TMPDIR="data/interim/{wildcards.stage}/tmp/gtdbtk/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_FNA="$TMPDIR/df_fna_gtdbtk.txt"
        INPUT_JSON="$TMPDIR/df_json_gtdbtk.txt"
        echo '{input.fna}' > $INPUT_FNA
        echo '{input.json_list}' > $INPUT_JSON
        if [ -s diff.txt ]; then
            python workflow/bgcflow/bgcflow/data/gtdbtk_prep.py $INPUT_FNA $INPUT_JSON {output.fnadir} {output.fnalist} 2>> {log}
        else
            mkdir -p {output.fnadir}
            touch {output.fnalist}
        fi
        rm $INPUT_FNA
        rm $INPUT_JSON
        rm -r $TMPDIR
        """

rule gtdbtk:
    input:
        gtdbtk="resources/gtdbtk/",
        fnadir="data/interim/{stage}/gtdbtk/{name}/fasta/",
    output:
        fnalist="data/interim/{stage}/gtdbtk/{name}/fasta_list_success.txt",
        gtdbtk_dir=directory("data/interim/{stage}/gtdbtk/{name}/result/"),
        tmpdir=temp(directory("data/interim/{stage}/gtdbtk/{name}_tmp/")),
        summary_interim="data/interim/{stage}/gtdbtk/{name}/result/classify/gtdbtk.bac120.summary.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/{stage}/gtdbtk/gtdbtk/gtdbtk_{name}.log",
    threads: 32
    params:
        ani_screen=ani_screen,
    shell:
        """
        set -e
        mkdir -p {output.tmpdir}
        gtdbtk classify_wf --genome_dir {input.fnadir} --out_dir {output.gtdbtk_dir} --cpus {threads} --pplacer_cpus 1 --tmpdir {output.tmpdir} {params.ani_screen} &>> {log}
        tail -n +2 {output.summary_interim} | cut -f1 > {output.fnalist}
        """

rule gtdbtk_fna_fail:
    output:
        "data/interim/{stage}/gtdbtk/{name}/fasta_list_fail.txt",
    log:
        "logs/gtdbtk/{stage}/gtdbtk/gtdbtk_{name}.log",
    shell:
        """
        echo "WARNING: No genomes are eligible for GTDB-Tk classification. Please check if the genome ids already exists in GTDB. Returning empty ouput." > {log}
        echo -e "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" > {output}
        """

def evaluate_gtdbtk_input(wildcards):
    """
    Evaluate whether there is a valid genome to use for GTDB-Tk or not based on the content of the output file.
    The function reads the content of the output file using the method open() of the returned file.
    This way, Snakemake is able to automatically download the file if it is generated in a cloud environment
    without a shared filesystem. If the file has content, the function returns the path to the
    `fasta_list_success.txt` file, otherwise it returns the path to the `fasta_list_fail.txt` file.

    Args:
        wildcards (object): A Snakemake wildcard object.

    Returns:
        str: The path to the `fasta_list_success.txt` file if the file has content, otherwise the path to the
        `fasta_list_fail.txt` file.
    """
    gtdtbk_input = checkpoints.prepare_gtdbtk_input.get(stage=wildcards.stage, name=wildcards.name).output["fnalist"]
    sys.stderr.write(f"GTDB-Tk checkpoint - Reading input file: {gtdtbk_input}\n")
    with gtdtbk_input.open() as f:
        textfile = f.readlines()
        sys.stderr.write(f"GTDB-Tk checkpoint - Found {len(textfile)} genomes to process...\n")
        if len(textfile) > 0:
            sys.stderr.write("GTDB-Tk checkpoint - Running GTDB-Tk classify_wf...\n")
            return "data/interim/{stage}/gtdbtk/{name}/fasta_list_success.txt",
        else:
            sys.stderr.write("GTDB-Tk checkpoint - No input passed. Skipping GTDB-tk run...\n")
            return "data/interim/{stage}/gtdbtk/{name}/fasta_list_fail.txt",

rule evaluate_gtdbtk_input:
    input:
        evaluate_gtdbtk_input
    output:
        summary_processed="data/processed/{stage}/{name}/tables/gtdbtk.bac120.summary.tsv"
    log:
        "logs/{stage}/gtdbtk/gtdbtk/evaluate_gtdbtk_{name}.log",
    shell:
        "cp {input} {output.summary_processed} 2>> {log}"
