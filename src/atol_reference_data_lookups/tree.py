import shelve
from pathlib import Path

import pandas as pd
import skbio.io
import skbio.tree._exception
import skbio.tree._tree
from skbio.tree import TreeNode

from atol_reference_data_lookups import logger
from atol_reference_data_lookups.cache import compute_sha256, open_cache
from atol_reference_data_lookups.io import read_augustus_mapping


def read_taxdump_file(file_path, cache_dir, scheme):
    """
    Reads the taxdump file and caches it in a shelve file.

    Return a tuple of the data and a boolean indicating whether the cache was
    updated.
    """
    cache_file = Path(cache_dir, f"{Path(file_path).stem}_{scheme}.db")
    Path.mkdir(cache_file.parent, exist_ok=True, parents=True)
    current_checksum = compute_sha256(file_path)

    with shelve.open(cache_file) as cache:
        if (
            "data" in cache
            and "checksum" in cache
            and cache["checksum"] == current_checksum
        ):
            logger.info(f"Reading {scheme} from cache {cache_file}")
            return (cache["data"], False)
        else:
            data = skbio.io.read(
                file_path, "taxdump", into=pd.DataFrame, scheme=scheme
            )
            logger.info(f"Writing {scheme} to cache {cache_file}")
            cache["data"] = data
            cache["checksum"] = current_checksum
            return (data, True)


def read_taxdump_nodes(file_path, cache_dir):
    """
    Read the full nodes.dmp file into a DataFrame indexed by tax_id,
    with caching.
    """
    cache_file = Path(cache_dir, "nodes_full_df.db")
    nodes_full_df_checksum = compute_sha256(file_path)

    with shelve.open(cache_file) as cache:
        if (
            "nodes_full_df" in cache
            and "nodes_full_df_checksum" in cache
            and cache["nodes_full_df_checksum"] == nodes_full_df_checksum
        ):
            logger.info(f"Reading full nodes from cache {cache_file}")
            return (cache["nodes_full_df"], False)
        else:
            data = skbio.io.read(
                file_path, "taxdump", into=pd.DataFrame, scheme="nodes"
            )
            logger.info(f"Writing full nodes to cache {cache_file}")
            cache["nodes_full_df"] = data
            cache["nodes_full_df_checksum"] = nodes_full_df_checksum
            return (data, True)


def generate_taxonomy_tree(names, nodes, cache_dir, update_tree=False):
    cache_file = Path(cache_dir, "taxonomy_tree.db")
    with shelve.open(cache_file) as cache:
        if "tree" in cache and not update_tree:
            logger.info(f"Reading taxonomy tree from {cache_file}")
            return cache["tree"]
        else:
            logger.info("Generating taxonomy tree")
            tree = TreeNode.from_taxdump(nodes)
            logger.info("Indexing tree")
            tree.index_tree()
            cache["tree"] = tree
            return tree


def get_node(tree, taxid):
    try:
        node = tree.find(int(taxid))
    except skbio.tree._exception.MissingNodeError:
        logger.debug(f"Node {taxid} not found, trying a string search")
        node = tree.find(str(taxid))
        logger.debug(f"    ... found {node}")

    if isinstance(node, skbio.tree._tree.TreeNode):
        return node
    else:
        logger.warning(f"Node for taxid {taxid} not found in tree.")


def find_lower_ranks(tree, top_rank="species", excluded_ranks=None):
    if excluded_ranks is None:
        excluded_ranks = ["no rank"]
    rank_list = _recursive_find_lower_ranks(tree, top_rank)
    return [rank for rank in sorted(set(rank_list)) if rank not in excluded_ranks]


def _recursive_find_lower_ranks(
    tree, top_rank="species", rank_list=None, top_rank_or_lower=False
):
    if rank_list is None:
        rank_list = []
    if tree.rank == top_rank:
        top_rank_or_lower = True
    if top_rank_or_lower:
        for node in tree.traverse():
            rank_list.append(node.rank)
    else:
        for node in tree.children:
            _recursive_find_lower_ranks(node, top_rank, rank_list, top_rank_or_lower)
    return rank_list


def generate_augustus_tree(
    tree, taxids_to_augustus_dataset_mapping, cache_dir, update_tree=False
):
    cache_file = Path(cache_dir, "augustus_tree.db")
    taxids_to_augustus_dataset_mapping_checksum = compute_sha256(
        taxids_to_augustus_dataset_mapping
    )

    logger.info(
        f"Reading Augustus dataset mapping from {taxids_to_augustus_dataset_mapping}"
    )
    augustus_mapping = read_augustus_mapping(taxids_to_augustus_dataset_mapping)
    logger.info(
        f"    ... found {len(augustus_mapping.keys())} datasets in Augustus mapping file"
    )

    with shelve.open(cache_file) as cache:
        if (
            "augustus_tree" in cache
            and "augustus_tip_names" in cache
            and "taxids_to_augustus_dataset_mapping_checksum" in cache
            and cache["taxids_to_augustus_dataset_mapping_checksum"]
            == taxids_to_augustus_dataset_mapping_checksum
            and not update_tree
        ):
            logger.info(f"Reading Augustus tree from {cache_file}")
            return augustus_mapping, cache["augustus_tree"], cache["augustus_tip_names"]

        else:
            logger.info("Pruning tree for Augustus datasets")
            augustus_tree = tree.copy(deep=True)
            initial_node_count = int(augustus_tree.count())

            augustus_nodes = [
                get_node(augustus_tree, x) for x in augustus_mapping.keys()
            ]
            augustus_node_names = [x.name for x in augustus_nodes]
            logger.debug(
                f"Found {len(augustus_node_names)} Augustus taxids in tree:\n{augustus_node_names}"
            )

            for node in augustus_nodes:
                if node.has_children():
                    node_children = node.children
                    logger.debug(
                        f"Node {node.name} has children {[x.name for x in node_children]}"
                    )
                    for child in node_children:
                        logger.debug(f"Removing node {child.name}")
                        removed = node.remove(child)
                        logger.debug(f"Removed? {removed}")

                    if not node.is_tip():
                        node.pop()

                logger.debug(f"Is node {node.name} a tip? {node.is_tip()}")

            sheared_augustus_tree = augustus_tree.shear(
                names=augustus_node_names,
                prune=False,
                inplace=False,
                strict=False,
            )

            final_node_count = int(sheared_augustus_tree.count())
            nodes_removed = initial_node_count - final_node_count

            logger.debug(f"    ... NCBI tree had {initial_node_count} nodes.")
            logger.debug(f"    ... Removed {nodes_removed} nodes.")
            logger.debug(f"    ... Augustus tree has {final_node_count} nodes.")

            augustus_tip_names = [x.name for x in sheared_augustus_tree.tips()]

            logger.info("Caching the pruned Augustus tree")
            cache["augustus_tree"] = sheared_augustus_tree
            cache["augustus_tip_names"] = augustus_tip_names
            cache["taxids_to_augustus_dataset_mapping_checksum"] = (
                taxids_to_augustus_dataset_mapping_checksum
            )

            return augustus_mapping, sheared_augustus_tree, augustus_tip_names
