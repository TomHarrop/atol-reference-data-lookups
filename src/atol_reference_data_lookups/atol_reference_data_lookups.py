from .taxdump_tree import TaxdumpTree
from atol_reference_data_lookups import logger
from argparse import ArgumentParser, Namespace
from pathlib import Path
import importlib.resources as pkg_resources
import os
import sys
import json


def parse_args() -> Namespace:
    # Note to self, this returns a Path
    package_files_path = pkg_resources.files(__package__)

    parser = ArgumentParser()

    input_group = parser.add_argument_group("Input")
    ref_group = parser.add_argument_group("Reference data")
    options_group = parser.add_argument_group("General options")

    taxid_group = input_group.add_mutually_exclusive_group()

    taxid_group.add_argument("--taxid", help="A single NCBI TaxId to look up", type=int)
    taxid_group.add_argument(
        "--taxid-list",
        help=(
            """
            A file containing a list NCBI TaxIds to look up, one per line.
            """
        ),
        type=Path,
    )

    ref_group.add_argument(
        "--nodes", required=True, help="NCBI nodes.dmp file from taxdump", type=Path
    )

    ref_group.add_argument(
        "--names", required=True, help="NCBI names.dmp file from taxdump", type=Path
    )

    ref_group.add_argument(
        "--taxids_to_busco_dataset_mapping",
        required=True,
        help=(
            """
              BUSCO placement file from
              https://busco-data.ezlab.org/v5/data/placement_files/
            """
        ),
        type=Path,
    )

    ref_group.add_argument(
        "--taxids_to_augustus_dataset_mapping",
        help=(
            """
            File that maps Augustus datasets to NCBI TaxIDs. See
            config/taxid_to_augustus_dataset.tsv
            """
        ),
        default=Path(package_files_path, "config", "taxid_to_augustus_dataset.tsv"),
        type=Path,
    )

    options_group.add_argument(
        "--cache_dir",
        help=(
            """
            Directory to cache the NCBI taxonomy after processing
            """
        ),
        default=Path(
            os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache")),
            "atol_reference_data_lookups",
        ),
    )

    return parser.parse_args()


def read_taxid_list(taxid_list_file: Path) -> list[int]:
    with open(taxid_list_file, "rt") as f:
        taxid_list = [int(x) for x in f.read().splitlines()]
    return taxid_list


def write_json_output(data_dict: dict) -> None:
    json.dump(data_dict, sys.stdout)


def main() -> None:
    args = parse_args()

    if args.taxid is not None:
        query_taxids = [args.taxid]
    else:
        query_taxids = read_taxid_list(taxid_list_file=args.taxid_list)

    taxdump_tree = TaxdumpTree(
        args.nodes,
        args.names,
        args.taxids_to_busco_dataset_mapping,
        args.taxids_to_augustus_dataset_mapping,
        args.cache_dir,
    )

    logger.info(f"Looking up {len(query_taxids)} query_taxids")

    taxonomy_reference_data = {}
    for query_taxid in query_taxids:
        ancestor_taxids = taxdump_tree.get_ancestor_taxids(query_taxid)
        genetic_code_id, mitochondrial_genetic_code_id = taxdump_tree.get_genetic_codes(
            query_taxid
        )

        taxonomy_reference_data[query_taxid] = {
            "busco_dataset_name": taxdump_tree.get_busco_lineage(
                query_taxid, ancestor_taxids
            ),
            "augustus_dataset_name": taxdump_tree.get_augustus_lineage(
                query_taxid, ancestor_taxids
            ),
            "genetic_code_id": int(genetic_code_id),
            "mitochondrial_genetic_code_id": int(mitochondrial_genetic_code_id),
        }

    logger.info("Finished lookups")

    write_json_output(taxonomy_reference_data)
