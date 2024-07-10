#datasets summary genome taxon Lactobacillaceae  --assembly-source refseq --as-json-lines | dataformat tsv genome > Lactobacillaceae_refseq.tsv

rule ncbi_dataset_for_taxon:
    output:
        tsv="data/interim/ncbi_datasets/taxon/{taxon}_original.tsv",
    params:
        taxon=lambda wildcards: wildcards.taxon.replace('_', ' '),
        fields=",".join([
            "accession",
            "organism-name",
            "organism-tax-id",
            "source_database",
            "assminfo-name",
            "assminfo-level",
            "ani-best-ani-match-assembly",
            "checkm-completeness",
            "checkm-contamination",
            "assmstats-gc-percent",
            "assmstats-number-of-scaffolds",
            "assmstats-number-of-contigs",
            "assmstats-total-sequence-len",
            "assmstats-contig-n50",
            "assmstats-scaffold-n50"]),
        reference=lambda wildcards: ("--reference" if TAXONS.loc[wildcards.taxon, "reference_only"] else ""),
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_for_taxon_{taxon}.log",
    shell:
        """
            datasets summary genome taxon {params.taxon} --assembly-version latest --assembly-source RefSeq {params.reference} --as-json-lines | dataformat tsv genome --fields {params.fields} > {output.tsv} 2> {log}
        """

checkpoint ncbi_dataset_tsv_to_samples_csv:
    input:
        tsv="data/interim/ncbi_datasets/taxon/{taxon}.tsv",
    output:
        csv="data/interim/ncbi_datasets/taxon/{taxon}.csv",
    log:
        "logs/ncbi_datasets/ncbi_dataset_tsv_to_csv_{taxon}.log",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", low_memory=False, header=0, index_col=0)
        df["source"] = "ncbi_taxon"
        try:
            CUSTOM_SAMPLES
        except NameError:
            pass
        else:
            if not CUSTOM_SAMPLES is None:
                df = df.join(CUSTOM_SAMPLES)
        df["genome_id"] = df.index
        df["closest_placement_reference"] = df["ani-best-ani-match-assembly"] if "ani-best-ani-match-assembly" in df.columns else ""
        df.to_csv(output.csv)

def get_accessions_for_taxon(taxon):
    genome_list = checkpoints.ncbi_dataset_tsv_to_samples_csv.get(taxon=taxon).output.csv
    df = pd.read_csv(genome_list)
    return df.index.to_list()
def get_all_accessions():
    accessions = []
    for taxon in TAXONS.index.to_list():
        accessions.extend(get_accessions_for_taxon(taxon))
    return accessions
def get_taxon_for_accession(accession):
    for taxon in TAXONS.index.to_list():
        genome_list = checkpoints.ncbi_dataset_tsv_to_samples_csv.get(taxon=taxon).output.csv
        df = pd.read_csv(genome_list)
        if accession in df.index:
            return taxon
    raise ValueError(f"No taxon found for {accession}")


rule ncbi_dataset_download_genome_for_taxon_dehydrated:
    input:
        genome_list="data/interim/ncbi_datasets/taxon/{taxon}.genome_list",
    output:
        dehydrated_dataset="data/interim/ncbi_datasets/taxon/{taxon}.dehydrated.zip",
        readme="data/interim/ncbi_datasets/datasets/{taxon}/README.md",
        fetch="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl",
        dataset_catalog="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/dataset_catalog.json",
    params:
        taxon_out_dir="data/interim/ncbi_datasets/datasets/{taxon}/"
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_download_genome_for_taxon_dehydrated_{taxon}.log",
    shell:
        """
            datasets download genome accession --inputfile {input.genome_list} --dehydrated --filename {output.dehydrated_dataset} > {log} 2>&1
            unzip -o -d {params.taxon_out_dir} {output.dehydrated_dataset}
        """

rule ncbi_dataset_rehydrate:
    input:
        genome_list="data/interim/ncbi_datasets/taxon/{taxon}.genome_list",
        fetch="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl",
    output:
        dummy="data/interim/ncbi_datasets/taxon/{taxon}.dummy",
    params:
        taxon_out_dir="data/interim/ncbi_datasets/datasets/{taxon}/"
        # dataset_dir=lambda wildcards: f"data/interim/ncbi_datasets/datasets/{wildcards.taxon}/",
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_redydrate_{taxon}.log",
    shell:
        """
            datasets rehydrate --directory {params.taxon_out_dir} > {log} 2>&1
            while IFS="" read -r p || [ -n "$p" ]
            do
                ls -w 1 {params.taxon_out_dir}/ncbi_dataset/data/$p/${{p}}_*_genomic.fna | wc -l | grep "^1$" > /dev/null
                mv {params.taxon_out_dir}/ncbi_dataset/data/$p/${{p}}_*_genomic.fna {params.taxon_out_dir}/ncbi_dataset/data/$p/$p.fna
            done < {input.genome_list}
            touch {output.dummy}
        """

rule ncbi_insert_custom_genomes:
    input:
        dummy="data/interim/ncbi_datasets/taxon/{taxon}.dummy",
        samples_file="data/interim/custom_genomes/samples.csv",
    output:
        dummy="data/interim/ncbi_datasets/taxon/{taxon}-custom.dummy",
        tmp=temp("data/interim/ncbi_datasets/temp/{taxon}/"),
    params:
        taxon_out_dir="data/interim/ncbi_datasets/datasets/{taxon}/"
        # dataset_dir=lambda wildcards: f"data/interim/ncbi_datasets/datasets/{wildcards.taxon}/",
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        """
            while IFS="" read -r p || [ -n "$p" ]
            do
                GENOME_ID=`echo "$p" | cut --delimiter=',' --fields=1`
                GENOME_SOURCE=`echo "$p" | cut --delimiter=',' --fields=2`
                GENOME_PATH=`echo "$p" | cut --delimiter=',' --fields=3`
                mkdir -p {params.taxon_out_dir}/ncbi_dataset/data/$GENOME_ID/
                if [[ "$GENOME_SOURCE" == "local" ]]; then
                    cp "$GENOME_PATH" {params.taxon_out_dir}/ncbi_dataset/data/$GENOME_ID/$GENOME_ID.fna
                elif [[ "$GENOME_SOURCE" == "local" ]]; then
                    datasets download genome accession $GENOME_ID
                    unzip {output.tmp}/ncbi_dataset.zip -d {output.tmp} -x README.md
                    rm {output.tmp}/ncbi_dataset.zip
                    mv {output.tmp}/ncbi_dataset/data/$GENOME_ID/*_genomic.fna {params.taxon_out_dir}/ncbi_dataset/data/$GENOME_ID/$GENOME_ID.fna
                    cat {output.tmp}/ncbi_dataset/data/assembly_data_report.jsonl >> {params.taxon_out_dir}/ncbi_dataset/data/assembly_data_report.jsonl
                    rm -r {output.tmp}/ncbi_dataset/
                else
                    echo "Incorrect source '$GENOME_SOURCE' for '$GENOME_ID'"
                    exit 1
                fi
            done < {input.samples_file}
            touch {output.dummy}
        """

rule ncbi_dataset_collect:
    input:
        dummy=lambda wildcards: expand("data/interim/ncbi_datasets/taxon/{taxon}.dummy", taxon=get_taxon_for_accession(wildcards.accession)),
        dummy_custom=lambda wildcards: expand("data/interim/ncbi_datasets/taxon/{taxon}-custom.dummy", taxon=get_taxon_for_accession(wildcards.accession)),
        jsonl_report=lambda wildcards: expand("data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl", taxon=get_taxon_for_accession(wildcards.accession)),
    output:
        fna="data/interim/fasta/{accession}.fna",
        json_report="data/interim/assembly_report/{accession}.json",
    params:
        fna=lambda wildcards: expand("data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/{accession}/{accession}.fna", taxon=get_taxon_for_accession(wildcards.accession), accession=wildcards.accession),
    conda:
        "../envs/data_processing.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_collect_{accession}.log",
    shell:
        """
            if [ ! -f {output.fna} ]
            then
                mv {params.fna} {output.fna}
            fi
            grep '^{{"accession":"{wildcards.accession}"' {input.jsonl_report} > {output.json_report}
            python workflow/scripts/add_info_to_assembly_report.py {output.json_report}
        """