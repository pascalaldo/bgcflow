if NCBI == []:
    pass
else:

    rule ncbi_genome_download:
        output:
            fna="data/interim/fasta/{ncbi}.fna",
            assembly_report="data/interim/assembly_report/{ncbi}.txt",
            json_report="data/interim/assembly_report/{ncbi}.json",
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/ncbi/ncbi_genome_download/ncbi_genome_download_{ncbi}.log",
        params:
            kingdom=KINGDOM
        shell:
            """
            if [[ {wildcards.ncbi} == GCF* ]]
            then
                source="refseq"
            elif [[ {wildcards.ncbi} == GCA* ]]
            then
                source="genbank"
            else
                echo "accession must start with GCA or GCF" >> {log}
            fi
            ncbi-genome-download -s $source -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download -P -N --verbose {params.kingdom} 2>> {log}
            gunzip -c data/raw/ncbi/download/$source/{params.kingdom}/{wildcards.ncbi}/*.fna.gz > {output.fna}
            cp data/raw/ncbi/download/$source/{params.kingdom}/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/$source/{params.kingdom}/{wildcards.ncbi}
            python workflow/bgcflow/bgcflow/data/get_assembly_information.py {output.assembly_report} {output.json_report} {wildcards.ncbi} 2>> {log}
            """

    rule extract_ncbi_information:
        input:
            all_json=lambda wildcards: get_ncbi_assembly_inputs(
                wildcards.name, DF_SAMPLES
            ),
        output:
            ncbi_meta_path=report(
                "data/processed/{name}/tables/df_ncbi_meta.csv",
                caption="../report/table-ncbi_meta.rst",
                category="{name}",
                subcategory="NCBI Genome Overview",
            ),
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/ncbi/extract_ncbi_information/extract_ncbi_information-{name}.log",
        shell:
            """
            TMPDIR="data/interim/tmp/{wildcards.name}"
            mkdir -p $TMPDIR
            INPUT_JSON="$TMPDIR/df_ncbi_meta_input.txt"
            echo '{input.all_json}' > $INPUT_JSON
            python workflow/bgcflow/bgcflow/data/extract_ncbi_information.py \
                $INPUT_JSON {output.ncbi_meta_path} 2>> {log}
            rm $INPUT_JSON
            """

    rule download_patric_tables:
        output:
            patric_genome_summary="resources/patric_meta/genome_summary",
            patric_genome_metadata="resources/patric_meta/genome_metadata",
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/patric/download_patric_tables.log",
        shell:
            """
            wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary -O {output.patric_genome_summary} 2>> {log}
            wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata -O {output.patric_genome_metadata} 2>> {log}
            """

    rule extract_patric_meta:
        input:
            ncbi_meta_path="data/processed/{name}/tables/df_ncbi_meta.csv",
            patric_genome_summary="resources/patric_meta/genome_summary",
            patric_genome_metadata="resources/patric_meta/genome_metadata",
        output:
            patric_meta_path=report(
                "data/processed/{name}/tables/df_patric_meta.csv",
                caption="../report/table-patric_meta.rst",
                category="{name}",
                subcategory="Patric Genome Overview",
            ),
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/patric/extract_patric_information/extract_patric_information-{name}.log",
        shell:
            """
            python workflow/bgcflow/bgcflow/data/extract_patric_meta.py \
                {input.ncbi_meta_path} {input.patric_genome_summary} {input.patric_genome_metadata} {output.patric_meta_path} 2>> {log}
            """

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
