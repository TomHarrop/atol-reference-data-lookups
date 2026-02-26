from .taxdump_tree import TaxdumpTree
from atol_reference_data_lookups import logger
from argparse import ArgumentParser
from pathlib import Path
import importlib.resources as pkg_resources
import os


def parse_args():
    # Note to self, this returns a Path
    package_files_path = pkg_resources.files(__package__)

    parser = ArgumentParser()

    input_group = parser.add_argument_group("Input")
    options_group = parser.add_argument_group("General options")

    input_group.add_argument(
        "--nodes", required=True, help="NCBI nodes.dmp file from taxdump", type=Path
    )

    input_group.add_argument(
        "--names", required=True, help="NCBI names.dmp file from taxdump", type=Path
    )

    input_group.add_argument(
        "--taxids_to_busco_dataset_mapping",
        required=True,
        help="BUSCO placement file from https://busco-data.ezlab.org/v5/data/placement_files/",
        type=Path,
    )

    input_group.add_argument(
        "--taxids_to_augustus_dataset_mapping",
        help="File that maps Augustus datasets to NCBI TaxIDs. See config/taxid_to_augustus_dataset.tsv",
        default=Path(package_files_path, "config", "taxid_to_augustus_dataset.tsv"),
        type=Path,
    )

    options_group.add_argument(
        "--cache_dir",
        help="Directory to cache the NCBI taxonomy after processing",
        default=Path(
            os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache")),
            "atol_reference_data_lookups",
        ),
    )

    return parser.parse_args()


def main():

    args = parse_args()

    taxump_tree = TaxdumpTree(
        args.nodes,
        args.names,
        args.taxids_to_busco_dataset_mapping,
        args.taxids_to_augustus_dataset_mapping,
        args.cache_dir,
    )

    raise ValueError(taxump_tree)
