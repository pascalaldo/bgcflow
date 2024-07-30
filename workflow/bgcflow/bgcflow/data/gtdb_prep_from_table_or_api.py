import json
import logging
import sys
from pathlib import Path

import pandas as pd
import requests

from random import randint
from time import sleep


log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

metadata_keywords = {
    "metadata_nucleotide": [
        "trna_aa_count",
        "contig_count",
        "n50_contigs",
        "longest_contig",
        "scaffold_count",
        "n50_scaffolds",
        "longest_scaffold",
        "genome_size",
        "gc_percentage",
        "ambiguous_bases",
        "gc_count",
        "l50_contigs",
        "l50_scaffolds",
        "mean_contig_length",
        "mean_scaffold_length",
    ],
    "metadata_gene": [
        "checkm_completeness",
        "checkm_contamination",
        "checkm_strain_heterogeneity",
        "lsu_5s_count",
        "ssu_count",
        "lsu_23s_count",
        "protein_count",
        "coding_density",
        "lsu_23s_contig_len",
        "lsu_23s_length",
        "lsu_23s_query_id",
        "lsu_5s_contig_len",
        "lsu_5s_length",
        "lsu_5s_query_id",
        "lsu_silva_23s_blast_align_len",
        "lsu_silva_23s_blast_bitscore",
        "lsu_silva_23s_blast_evalue",
        "lsu_silva_23s_blast_perc_identity",
        "lsu_silva_23s_blast_subject_id",
        "lsu_silva_23s_taxonomy",
        "checkm_marker_count",
        "checkm_marker_lineage",
        "checkm_marker_set_count",
        "coding_bases",
        "mimag_high_quality",
        "mimag_low_quality",
        "mimag_medium_quality",
        "ssu_contig_len",
        "ssu_gg_blast_align_len",
        "ssu_gg_blast_bitscore",
        "ssu_gg_blast_evalue",
        "ssu_gg_blast_perc_identity",
        "ssu_gg_blast_subject_id",
        "ssu_gg_taxonomy",
        "ssu_length",
        "ssu_query_id",
        "ssu_silva_blast_align_len",
        "ssu_silva_blast_bitscore",
        "ssu_silva_blast_evalue",
        "ssu_silva_blast_perc_identity",
        "ssu_silva_blast_subject_id",
        "ssu_silva_taxonomy",
        "total_gap_length",
        "trna_count",
        "trna_selenocysteine_count",
    ],
    "metadata_ncbi": [
        "ncbi_genbank_assembly_accession",
        "ncbi_strain_identifiers",
        "ncbi_assembly_level",
        "ncbi_assembly_name",
        "ncbi_assembly_type",
        "ncbi_bioproject",
        "ncbi_biosample",
        "ncbi_country",
        "ncbi_date",
        "ncbi_genome_category",
        "ncbi_genome_representation",
        "ncbi_isolate",
        "ncbi_isolation_source",
        "ncbi_lat_lon",
        "ncbi_molecule_count",
        "ncbi_protein_count",
        "ncbi_refseq_category",
        "ncbi_seq_rel_date",
        "ncbi_spanned_gaps",
        "ncbi_species_taxid",
        "ncbi_ssu_count",
        "ncbi_submitter",
        "ncbi_taxid",
        "ncbi_total_gap_length",
        "ncbi_translation_table",
        "ncbi_trna_count",
        "ncbi_unspanned_gaps",
        "ncbi_contig_count",
        "ncbi_contig_n50",
        "ncbi_ncrna_count",
        "ncbi_organism_name",
        "ncbi_rrna_count",
        "ncbi_scaffold_count",
        "ncbi_scaffold_l50",
        "ncbi_scaffold_n50",
        "ncbi_scaffold_n75",
        "ncbi_scaffold_n90",
        "ncbi_total_length",
        "ncbi_ungapped_length",
        "ncbi_wgs_master",
    ],
    "metadata_type_material": [
        "gtdb_type_designation_ncbi_taxa",
        "gtdb_type_designation_ncbi_taxa_sources",
        "gtdb_type_species_of_genus",
    ],
    "metadataTaxonomy": [
        "ncbi_taxonomy",
        "ncbi_taxonomy_unfiltered",
        "gtdb_representative",
        "gtdb_genome_representative",
        "ncbi_type_material_designation",
    ],
}


def empty_result(genome_id):
    empty_result = {
        "genome_id": genome_id,
        "gtdb_url": None,
        "gtdb_release": None,
        "gtdb_taxonomy": {
            "domain": "d__",
            "phylum": "p__",
            "class": "c__",
            "order": "o__",
            "family": "f__",
            "genus": "g__",
            "species": "s__",
        },
    }
    return empty_result


def create_gtdb_metadata_from_table_or_api(
    genomes_csv,
    gtdb_metadata_path,
    gtdb_release,
    out_path,
    tax_path,
    offline_mode=False,
    api_base_url="https://gtdb-api.ecogenomic.org",
    metadata_keywords=metadata_keywords,
):
    """
    This function creates GTDB metadata from a given table and writes it to a JSON file.

    Parameters:
    genome_id (str): The genome ID to be processed.
    gtdb_metadata_path (str): The path to the GTDB metadata file.
    gtdb_release (str): The GTDB release version.
    out_path (str): The path to the output file.
    metadata_keywords (dict): The metadata keywords to be processed.

    Returns:
    None
    """
    gtdb_mapping = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    logging.info("Reading GTDB metadata from %s", gtdb_metadata_path)
    df_gtdb_metadata = pd.read_csv(gtdb_metadata_path, sep="\t", low_memory=False)
    df_gtdb_metadata = df_gtdb_metadata.set_index("accession")

    df_genomes = pd.read_csv(genomes_csv, low_memory=False, index_col=0, header=0)
    genome_ids = df_genomes.index.to_list()

    with open(out_path, "w") as file:
        for genome_id in genome_ids:
            if genome_id.startswith("GCF_"):
                logging.debug("Detected Refseq entry")
                accession = f"RS_{genome_id}"
            elif genome_id.startswith("GCA_"):
                logging.debug("Detected Refseq entry")
                accession = f"GB_{genome_id}"

            if accession in df_gtdb_metadata.index:
                logging.info("Accession %s found in metadata", accession)
                metadata = df_gtdb_metadata.loc[accession, :].to_dict()
                output = populate_from_gtdb_table(
                    genome_id, metadata, gtdb_release, metadata_keywords, gtdb_mapping
                )
            else:
                logging.warning(
                    "Accession %s not found in table. Trying API.", accession
                )
                output = gtdb_prep_from_api(genome_id, df_genomes, tax_path, release=gtdb_release, offline_mode=offline_mode, api_base_url=api_base_url)

            logging.info("Metadata creation completed for genome_id %s", genome_id)

            Path(out_path).parent.mkdir(parents=True, exist_ok=True)
            logging.debug("Writing metadata to %s", out_path)
            json.dump(output, file)
            file.write("\n")


def populate_from_gtdb_table(
    genome_id, metadata, gtdb_release, metadata_keywords, gtdb_mapping
):
    output = {}
    output["genome_id"] = genome_id
    output[
        "gtdb_url"
    ] = f"https://api.gtdb.ecogenomic.org/genome/{genome_id}/taxon-history"
    output["gtdb_release"] = gtdb_release
    output["gtdb_taxonomy"] = {
        gtdb_mapping[v.split("__")[0]]: v for v in metadata["gtdb_taxonomy"].split(";")
    }
    output["metadata_url"] = f"https://api.gtdb.ecogenomic.org/genome/{genome_id}/card"
    output["metadata"] = {}
    output["metadata"]["genome"] = {"accession": genome_id, "name": genome_id}
    for keywords, values in metadata_keywords.items():
        output["metadata"][keywords] = {k: metadata[k] for k in values}
    output["metadata"]["detail"] = "genome_found"
    return output

def gtdb_prep_from_api(
    genome_id,
    df_samples,
    tax_path,
    release="R207",
    offline_mode=False,
    api_base_url="https://gtdb-api.ecogenomic.org",
):  # what happen if it does not find?
    """
    Given a genome id and the samples Pandas dataframe, write  a JSON file containing taxonomic information from GTDB API. The script will first search using the closest taxonomic placement (NCBI accession id), then using the genus information provided by user. If no information is provided, return an empty taxonomic information.
    """

    class EmptyTaxError(Exception):
        """Raised when this script returns empty dict"""

        pass

    class PlacementError(Exception):
        """Raised when this script returns empty dict"""

        pass

    def find_taxonomy(query, genome_id, gtdb_tax, api_base_url, offline_mode=False):
        """
        Helper script to decide taxonomic placement for a given query
        """
        if type(offline_mode) == str:
            offline_mode = True if offline_mode.lower() == "true" else False
        assert type(offline_mode) == bool, "offline_mode must be a boolean"
        logging.info(f"Offline mode: {offline_mode}")
        if offline_mode is True:
            logging.info(
                "WARNING: Running in offline mode. Returning empty dataframe..."
            )
            gtdb_tax = empty_result(genome_id)
            return gtdb_tax
        else:
            # If closest placement reference is provided, try finding taxonomic information from GTDB API
            if query.closest_placement_reference.values[0] != "":
                try:
                    logging.info(
                        "Inferring taxonomic placement from provided closest reference...."
                    )
                    gtdb_tax = get_ncbi_taxon_GTDB(
                        query.closest_placement_reference.values[0],
                        api_base_url,
                        release=release,
                    )
                    gtdb_tax["genome_id"] = genome_id
                except KeyError as e:
                    raise PlacementError(
                        f"Cannot infer taxonomic placement from provided closest reference. Make sure the accession id: {query.closest_placement_reference.values[0]} is part of GTDB release: {e}"
                    )

            # If NCBI accession is provided, try to find taxonomic information from GTDB API
            elif query.source.values[0] == "ncbi" and ("genus" in query.columns):
                try:
                    logging.info(
                        "Inferring taxonomic placement from NCBI accession...."
                    )
                    gtdb_tax = get_ncbi_taxon_GTDB(
                        query.genome_id.values[0], api_base_url, release=release
                    )
                except KeyError:
                    if query.genus.values[0] != "":
                        logging.info(
                            "Inferring taxonomic placement from provided genus information...."
                        )
                        gtdb_tax["genome_id"] = genome_id
                        gtdb_tax.update(
                            get_parent_taxon_GTDB(
                                query.genus.values[0], "genus", release
                            )
                        )
                        gtdb_tax["gtdb_taxonomy"][
                            "species"
                        ] = f"s__{gtdb_tax['gtdb_taxonomy']['genus'].split('__')[-1]} sp."
                    else:
                        gtdb_tax = empty_result(genome_id)

            # Try to get taxonomic information from genus information
            elif ("genus" in query.columns) and (query.genus.values[0] != ""):
                logging.info(
                    "Inferring taxonomic placement from provided genus information...."
                )
                gtdb_tax["genome_id"] = genome_id
                gtdb_tax.update(
                    get_parent_taxon_GTDB(query.genus.values[0], "genus", release)
                )
                try:
                    gtdb_tax["gtdb_taxonomy"][
                        "species"
                    ] = f"s__{gtdb_tax['gtdb_taxonomy']['genus'].split('__')[-1]} sp."
                except KeyError:
                    gtdb_tax = empty_result(genome_id)

            # If no information is found, return an empty dict
            else:
                gtdb_tax = empty_result(genome_id)

            return gtdb_tax

    # get query by subsetting samples df with genome id
    query = df_samples[df_samples.loc[:, "genome_id"] == genome_id].fillna("")

    # create empty container
    gtdb_tax = {}

    # Starting process
    logging.info(f"Fetching taxonomic information for {genome_id}...")
    # Go through user provided taxonomic placement
    if any(os.path.isfile(t) for t in tax_path.split()):
        try:
            gtdb_tax = get_user_defined_classification(genome_id, tax_path)
        except KeyError:
            logging.warning(
                f"{genome_id}: Not found in user provided taxonomic placement..."
            )
            gtdb_tax = find_taxonomy(
                query, genome_id, gtdb_tax, api_base_url, offline_mode=offline_mode
            )
    else:
        gtdb_tax = find_taxonomy(
            query, genome_id, gtdb_tax, api_base_url, offline_mode=offline_mode
        )

    if gtdb_tax == {}:
        raise EmptyTaxError(
            "Oops, this shouldn't happen. It returns an empty dict. Something is wrong with the script."
        )

    return gtdb_tax


def get_user_defined_classification(genome_id, tax_path):
    """
    Get taxonomic information from user provided GTDB-like output
    """
    shell_input = tax_path.split()
    dfList = [
        pd.read_csv(s, sep="\t").set_index("user_genome", drop=False)
        for s in shell_input
    ]
    df_tax_raw = pd.concat(dfList, axis=0)

    # drop duplicates! causes error
    logging.debug(f"Checking user provided taxonomy table from: {shell_input}")
    logging.info(
        "Checking if user provided taxonomy table contains duplicated genome ids.."
    )
    df_tax_dup = df_tax_raw["user_genome"].duplicated()
    if df_tax_dup.any():
        logging.warning(
            f"Found duplicated genome ids: {list(df_tax_raw[df_tax_dup]['user_genome'].unique())}"
        )
        duplicates_mask = df_tax_raw.duplicated(keep="first")
        df_tax = df_tax_raw[~duplicates_mask].copy()
        logging.debug("Making sure duplicated genome ids values are identical...")
        assert (
            not df_tax.duplicated().any()
        ), "Two or more genome ids have more than one taxonomic placement! Please check your taxonomy files!"
    else:
        df_tax = df_tax_raw.copy()

    level_dict = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    query = df_tax.loc[genome_id, "classification"].split(";")

    result = {}
    result["genome_id"] = genome_id
    result["gtdb_url"] = "user provided classification"
    result["gtdb_release"] = "unknown"
    result["gtdb_taxonomy"] = {level_dict[q.split("__")[0]]: q for q in query}
    logging.info("Using user provided GTDB classification.")
    return result


def get_ncbi_taxon_GTDB(accession, api_base_url, release="R207"):
    """
    Given an NCBI accession, return a json object of taxonomic information from GTDB API
    """

    def gtdb_api_request(accession, api_type):
        if api_type == "taxonomy":
            api_url = f"{api_base_url}/genome/{accession}/taxon-history"
        elif api_type == "summary":
            api_url = f"{api_base_url}/genome/{accession}/card"
        logging.debug(f"Requesting GTDB API: {api_url}")

        js = None
        max_retries = 3
        retries_left = max_retries
        while retries_left > 0:
            if retries_left != max_retries:
                sleep(randint(5, 20))
            response = requests.get(api_url)
            try:
                js = response.json()
                retries_left = 0
            except json.JSONDecodeError:
                logging.critical(
                    f"Cannot decode response from GTDB API. Make sure this is a valid url: {api_url}"
                )
                retries_left = retries_left - 1
                # raise
        if js is None:
            raise

        return js, api_url

    # Mapping to bgcflow format
    level_dict = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    # get taxonomic information
    js_tax, api_url_tax = gtdb_api_request(accession, "taxonomy")
    result = {}
    result["genome_id"] = accession
    result["gtdb_url"] = api_url_tax
    result["gtdb_release"] = release
    card_detail = "Genome found"
    try:
        if js_tax == []:
            logging.warning(
                f"Genome id: {accession} is in GTDB but has no taxonomic placement. Returning empty values."
            )
            result["gtdb_taxonomy"] = {
                level_dict[k]: f"{k}__" for k in level_dict.keys()
            }
            card_detail = (
                f"Genome not found - no taxonomic assignment in release {release}"
            )
        else:
            logging.info(js_tax)
            tax = [tax for tax in js_tax if tax["release"] == release]
            if len(tax) == 1:
                result["gtdb_taxonomy"] = tax[0]
            elif len(tax) == 0:
                logging.warning(
                    f"Genome id: {accession} is in GTDB but has no taxonomic placement in release {release}. Returning empty values."
                )
                result["gtdb_taxonomy"] = {
                    level_dict[k]: f"{k}__" for k in level_dict.keys()
                }
                card_detail = (
                    f"Genome not found - no taxonomic assignment in release {release}"
                )
            else:
                raise
            result["gtdb_taxonomy"].pop("release", None)
            result["gtdb_taxonomy"] = {
                level_dict[k]: result["gtdb_taxonomy"][k]
                for k in result["gtdb_taxonomy"].keys()
            }
    except KeyError as err:
        if err.args[0] == "gtdb_taxonomy":
            logging.critical(
                f"Malformed genome id: {accession}. Make sure to use the right NCBI genome accession format."
            )
            raise
        elif err.args[0] == release:
            logging.critical(f"Cannot find genome id: {accession} in GTDB API.")
            raise

    # get other metadata from genome summary
    js_sum, api_url_sum = gtdb_api_request(accession, "summary")
    result["metadata_url"] = api_url_sum
    result["metadata"] = js_sum

    if "detail" in result["metadata"].keys():
        pass
    else:
        result["metadata"]["detail"] = card_detail

    return result


def get_parent_taxon_GTDB(taxon, level, api_base_url, release="R207"):
    """
    Given a taxon and its level, return a json object of parent taxons from GTDB API
    """
    level_dict = {
        "domain": "d",
        "phylum": "p",
        "class": "c",
        "order": "o",
        "family": "f",
        "genus": "g",
        "species": "s",
    }

    try:
        query = f"{level_dict[level]}__{taxon}"
    except KeyError:
        logging.critical(
            f"Incorrect taxon level format. Please choose from available format: {list(level_dict.keys())}."
        )
        raise

    api_url = f"{api_base_url}/taxonomy/partial/{query}/all-releases"
    response = requests.get(api_url)
    js = response.json()
    result = {}
    result["gtdb_url"] = api_url
    result["gtdb_release"] = release
    result["gtdb_taxonomy"] = {}

    level_dict_rev = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    try:
        for item in js:
            if release in item["release"]:
                t = item["taxonomy"]
                taxonomy = {level_dict_rev[k]: t[k] for k in t.keys()}

                result["gtdb_taxonomy"] = taxonomy
    except TypeError:
        logging.critical(f"{js['message']}")
        raise

    return result


if __name__ == "__main__":
    create_gtdb_metadata_from_table_or_api(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], offline_mode=sys.argv[6], api_base_url=sys.argv[7])
