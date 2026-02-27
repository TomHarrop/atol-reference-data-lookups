#!/usr/bin/env python3

import skbio.tree._exception

from atol_reference_data_lookups import logger
from atol_reference_data_lookups.io import read_busco_mapping
from atol_reference_data_lookups.tree import (
    generate_augustus_tree,
    generate_taxonomy_tree,
    get_node,
    read_taxdump_file,
    read_taxdump_nodes,
)


class TaxdumpTree:

    def __init__(
        self,
        nodes_file,
        names_file,
        taxids_to_busco_dataset_mapping,
        taxids_to_augustus_dataset_mapping,
        cache_dir,
    ):
        logger.info(f"Reading NCBI taxonomy from {nodes_file}")
        self.nodes, nodes_changed = read_taxdump_file(
            nodes_file, cache_dir, "nodes_slim"
        )

        logger.info(f"Reading full NCBI nodes from {nodes_file}")
        self.nodes_full, nodes_full_changed = read_taxdump_nodes(
            nodes_file, cache_dir
        )

        logger.info(f"Reading NCBI taxon names from {names_file}")
        names, names_changed = read_taxdump_file(names_file, cache_dir, "names")

        update_tree = any([nodes_changed, names_changed])

        self.tree = generate_taxonomy_tree(
            names, self.nodes, cache_dir, update_tree=update_tree
        )

        logger.info(
            f"Reading BUSCO dataset mapping from {taxids_to_busco_dataset_mapping}"
        )
        self.busco_mapping = read_busco_mapping(taxids_to_busco_dataset_mapping)
        logger.info(
            f"    ... found {len(self.busco_mapping)} datasets in BUSCO mapping file"
        )

        self.augustus_mapping, self.augustus_tree, self.augustus_tip_names = (
            generate_augustus_tree(
                self.tree, taxids_to_augustus_dataset_mapping, cache_dir, update_tree
            )
        )

    def get_ancestor_taxids(self, taxid):
        logger.debug(f"Looking up ancestors for taxid {taxid}")
        node = get_node(self.tree, taxid)
        ancestor_taxids = [x.name for x in node.ancestors()]
        logger.debug(f"ancestor_taxids: {ancestor_taxids}")
        return ancestor_taxids

    def get_busco_lineage(self, taxid, ancestor_taxids):
        """
        Find the closest ancestor that is in the BUSCO taxid map and return the
        lineage name.
        """
        logger.debug(f"Looking up BUSCO dataset name for taxid {taxid}")
        logger.debug(f"Checking ancestor_taxids {ancestor_taxids}")

        for ancestor_taxid in ancestor_taxids:
            if int(ancestor_taxid) in self.busco_mapping:
                return self.busco_mapping[int(ancestor_taxid)]
            if ancestor_taxid in self.busco_mapping:
                return self.busco_mapping[ancestor_taxid]

        return None

    def get_augustus_lineage(self, taxid, ancestor_taxids):
        logger.debug(f"Looking up Augustus dataset name for taxid {taxid}")

        # Include taxid in the search, in case it has been trained
        search_taxids = [taxid] + list(ancestor_taxids)

        closest_taxid_in_augustus_tree = None
        while not closest_taxid_in_augustus_tree and len(search_taxids) > 0:
            current = search_taxids.pop(0)
            logger.debug(f"Checking Augustus tree for {current}")
            try:
                closest_taxid_in_augustus_tree = get_node(
                    self.augustus_tree, current
                )
            except skbio.tree._exception.MissingNodeError:
                logger.debug(f"Node {current} not in Augustus tree")

        if not closest_taxid_in_augustus_tree:
            return None

        logger.debug(
            f"Found closest_taxid_in_augustus_tree {closest_taxid_in_augustus_tree.name}"
        )

        dist_to_dataset = {}
        for augustus_taxid in self.augustus_mapping:
            logger.debug(f"Calculating distance to {augustus_taxid}")
            try:
                dest_node = get_node(self.augustus_tree, augustus_taxid)
                dist_to_dataset[augustus_taxid] = (
                    closest_taxid_in_augustus_tree.distance(
                        dest_node, use_length=False
                    )
                )
                logger.debug(f"    ...{dist_to_dataset[augustus_taxid]}")
            except skbio.tree._exception.MissingNodeError:
                logger.debug(f"{augustus_taxid} not in Augustus tree")

        logger.debug(f"Distances: {dist_to_dataset}")
        closest_dataset = min(dist_to_dataset, key=dist_to_dataset.get)

        logger.debug(
            f"Closest dataset is {closest_dataset} with dist {dist_to_dataset[closest_dataset]}"
        )

        return self.augustus_mapping[closest_dataset]

    def get_genetic_codes(self, taxid):
        """
        Look up the genetic_code_id and mitochondrial_genetic_code_id for a
        given taxid.

        Returns a tuple of (genetic_code_id, mitochondrial_genetic_code_id).
        """
        logger.debug(f"Looking up genetic codes for taxid {taxid}")
        row = self.nodes_full.loc[taxid]
        return (row["genetic_code_id"], row["mitochondrial_genetic_code_id"])
