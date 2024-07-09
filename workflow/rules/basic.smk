__version__ = "0.10.0"


# roary.smk #
def get_prokka_outputs(name, df_samples, ext="gff", path="prokka"):
    """
    Given a project name, find the corresponding sample file to use

    Arguments:
        name {str} -- project name
        df_samples {pd.DataFrame} -- sample table
        ext {str} -- file extension

    Returns:
        output {list} -- list of prokka outputs
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    assert path in ["prokka", "processed-genbank"]
    if path == "prokka":
        output = [f"data/interim/{path}/{s}/{s}.{ext}" for s in selection]
    elif path == "processed-genbank":
        output = [f"data/interim/{path}/{s}.{ext}" for s in selection]
    return output


# prokka.smk #
def get_prokka_refdb(genome_id, params, df_samples, mapping_file, config=config):
    """
    Given a genome id, find which prokka-db input to use.

    params:
        - "table" - will return the corresponding prokka-db table to use
        - "file" - will return the corresponding reference gbks
        - "params" - will return prokka protein params and the corresponding file

    Arguments:
        genome_id {str} -- genome id
        params {str} -- "table" or "file" or "params"
        df_samples {pd.DataFrame} -- sample table
        mapping_file {dict} -- mapping file between prokka-db and reference gbks
        config {dict} -- config file

    Returns:
        output {str} -- prokka-db table or reference gbks or prokka protein params
    """

    if not "ref_annotation" in df_samples.columns:
        if params == "file":
            return []
        else:
            return ""
    prokka_db = df_samples.loc[genome_id, "ref_annotation"].iloc[0]
    name = df_samples.loc[genome_id, "name"].iloc[0]

    if not os.path.isfile(str(prokka_db)):
        if params == "file":
            output = []
        else:
            output = ""
    elif params == "table":
        output = prokka_db
    elif params == "file":
        output = f"resources/prokka_db/{mapping_file[prokka_db]}.gbff"
    elif params == "params":
        output = f"--proteins resources/prokka_db/{mapping_file[prokka_db]}.gbff"
    else:
        sys.stderr.write(f"Second argument should be: table, file, or params.\n")
        raise
    return output

# automlst_wrapper.smk #
def get_automlst_inputs(name, df_samples):
    """
    Given a project name, find the corresponding sample file to use

    Arguments:
        name {str} -- project name
        df_samples {pd.DataFrame} -- sample table

    Returns:
        output {list} -- list of automlst gbk files
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/automlst_wrapper/{name}/{s}.gbk" for s in selection]
    return output

def get_samples_for_project_from_df(df_samples, name):
    cols = list({"sample_paths", "prokka-db", "gtdb_paths", "name"} & set(df_samples.columns.to_list()))
    df_exploded = df_samples.explode(cols)
    return df_exploded.loc[df_exploded["name"] == name, :]

def get_unclassified_accessions(taxon):
    import pandas as pd
    df = pd.read_csv(checkpoints.fix_gtdb_taxonomy.get(taxon=taxon).output.meta, low_memory=False, header=0, index_col=0)
    unclassified = df.loc[df["Species"] == "s__", :]
    return unclassified.index.tolist()
