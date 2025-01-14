rule original_locus_tags:
    input:
        ncbi_gff=fexpand("data/interim/all/gff/{accession}.gff", accession=RULE_FUNCTIONS["original_locus_tags"]["accession"]),
        prokka_gbk=fexpand("data/interim/{{stage}}/processed-genbank/{accession}.gbk", accession=RULE_FUNCTIONS["original_locus_tags"]["accession"]),
        samples="data/processed/species/samples/{name}.csv",
    output:
        csv="data/processed/{stage}/{name}/tables/df_locus_tag_mapping.csv",
    params:
        ncbi_gff="data/interim/all/gff/",
        prokka_gbk="data/interim/{stage}/processed-genbank/",
    conda:
        "../envs/locus_tag_mapping.yaml"
    log:
        "logs/{stage}/original_locus_tags/original_locus_tags-{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/map_locus_tags.py \
            --samples {input.samples} \
            --ncbi_gff {params.ncbi_gff} \
            --prokka_gbk {params.prokka_gbk} \
            -o {output.csv} > {log} 2>&1
        """