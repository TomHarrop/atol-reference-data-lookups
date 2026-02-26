from pathlib import Path
from argparse import ArgumentParser
import importlib.resources as pkg_resources


def parse_args():
    # Note to self, this returns a Path
    package_files_path = pkg_resources.files(__package__)

    parser = ArgumentParser()

    input_group = parser.add_argument_group("Input")

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

    return parser.parse_args()


def main():

    parse_args()

    raise ValueError("Here we are")
