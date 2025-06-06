#%
# final_output: data/processed/{name}/bigscape/result_as{version}/index.html
# description: Cluster BGCs using BiG-SCAPE
# category: Genome Mining
# link:
# - https://github.com/medema-group/BiG-SCAPE
# references:
# - Navarro-Muñoz, J.C., Selem-Mojica, N., Mullowney, M.W. et al. A computational
#   framework to explore large-scale biosynthetic diversity. [Nat Chem Biol 16, 60–68
#   (2020)](https://doi.org/10.1038/s41589-019-0400-9)
#%
rule install_bigscape:
    output:
        bigscape=directory("resources/BiG-SCAPE/"),
    conda:
        "../envs/bigscape.yaml"
    log:
        "logs/bigscape/bigscape-install_bigscape.log",
    params:
        release="1.1.9",
    shell:
        """
        (cd resources && wget https://github.com/medema-group/BiG-SCAPE/archive/refs/tags/v{params.release}.zip) &>> {log}
        (cd resources && unzip -o v{params.release}.zip && mv BiG-SCAPE-{params.release}/ BiG-SCAPE/ && rm v{params.release}.zip) &>> {log}
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz) &>> {log}
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm) &>> {log}
        """


rule bigscape:
    input:
        bigscape="resources/BiG-SCAPE",
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        antismash_dir="data/interim/bgcs/{name}/{version}/",
    output:
        index="data/interim/bigscape/{name}_antismash_{version}/index.html",
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir="data/interim/bigscape/{name}_antismash_{version}/",
        label="{name}_antismash_{version}",
    log:
        "logs/bigscape/{name}_antismash_{version}/bigscape.log",
    threads: 32
    shell:
        """
        python {input.bigscape}/bigscape.py -i {input.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose &>> {log}
        """


## Unused rule, might be deleted
rule bigscape_no_mibig:
    input:
        bigscape="resources/BiG-SCAPE",
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        antismash_dir="data/interim/bgcs/{name}/{version}/",
    output:
        index="data/interim/bigscape/no_mibig_{name}_antismash_{version}/index.html",
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir="data/interim/bigscape/no_mibig_{name}_antismash_{version}/",
        label="{name}_antismash_{version}",
    log:
        "logs/bigscape/{name}_antismash_{version}/bigscape.log",
    threads: 32
    shell:
        """
        python {input.bigscape}/bigscape.py -i {input.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --verbose --no_classify --mix &>> {log}
        """

rule copy_bigscape_zip:
    input:
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        index="data/interim/bigscape/no_mibig_{name}_antismash_{version}/index.html",
    output:
        zip="data/processed/{name}/bigscape/no_mibig_bigscape_as{version}.zip",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/copy_bigscape_zip/copy_bigscape_zip-{name}-{version}.log",
    shell:
        """
        topdir=$PWD
        (cd data/interim/bigscape && zip -r $topdir/{output.zip} no_mibig_{wildcards.name}_antismash_{wildcards.version} \
            -x {wildcards.name}_antismash_{wildcards.version}/cache/**\* &>> $topdir/{log})
        """


rule bigscape_to_cytoscape:
    input:
        index="data/interim/bigscape/{name}_antismash_{version}/index.html",
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        mibig_bgc_table="resources/mibig/df_mibig_bgcs.csv/",
        as_dir="data/interim/bgcs/{name}/{version}",
        df_genomes_path="data/processed/{name}/tables/df_antismash_{version}_summary.csv",
    output:
        output_dir=directory(
            "data/processed/{name}/bigscape/for_cytoscape_antismash_{version}"
        ),
        bgc_mapping="data/processed/{name}/bigscape/{name}_bigscape_as_{version}_mapping.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/bigscape_to_cytoscape/bigscape_to_cytoscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_bigscape_to_cytoscape.py data/interim/bigscape/{wildcards.name}_antismash_{wildcards.version}/ {input.as_dir} {input.df_genomes_path} {input.mibig_bgc_table} {output.output_dir} 2>> {log}
        cp {input.bgc_mapping} {output.bgc_mapping}
        """


rule copy_bigscape:
    input:
        index="data/interim/bigscape/{name}_antismash_{version}/index.html",
        cytoscape="data/processed/{name}/bigscape/for_cytoscape_antismash_{version}",
    output:
        main="data/processed/{name}/bigscape/result_as{version}/index.html",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/copy_bigscape_zip/copy_bigscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/bigscape_copy.py {input.index} data/processed/{wildcards.name}/bigscape/result_as{wildcards.version} 2>> {log}
        """
