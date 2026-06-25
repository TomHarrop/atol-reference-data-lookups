#!/usr/bin/env python3

import argparse
import json
from sys import stdout
import urllib.parse

import requests


def parse_arguments() -> argparse.Namespace:

    parser = argparse.ArgumentParser()

    _ = parser.add_argument("taxon_id", type=int)

    return parser.parse_args()


# FIXME. Hard coded defaults for now.
_api_url = "https://goat.genomehubs.org/api/v2/"

_endpoints = {
    "record": "record",
}

_headers = {"Accept": "application/json"}

# only asking for the first result
_params = {"result": "taxon", "size": 1, "taxonomy": "ncbi"}


def main():
    args = parse_arguments()

    records_url = urllib.parse.urljoin(_api_url, _endpoints.get("record", ""))

    taxid_results = requests.get(
        records_url, headers=_headers, params={**_params, "recordId": args.taxon_id}
    )

    json_results = taxid_results.json()

    first_hit_attributes = (
        json_results.get("records", {})[0].get("record", {}).get("attributes", {})
    )

    # get e.g. the genome_size attribute from this result
    # jq '.genome_size.value'
    json.dump(first_hit_attributes, stdout)


if __name__ == "__main__":
    main()
