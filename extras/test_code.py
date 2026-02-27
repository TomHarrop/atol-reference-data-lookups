from atol_reference_data_lookups.taxdump_tree import TaxdumpTree
from atol_reference_data_lookups.tree import get_node
from pathlib import Path

import skbio.io
import pandas as pd

nodes_path = Path("resources/new_taxdump/nodes.dmp")

taxdump_tree = TaxdumpTree(
    nodes_path,
    Path("resources/new_taxdump/names.dmp"),
    Path(
        "resources/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz"
    ),
    Path("src/atol_reference_data_lookups/config/taxid_to_augustus_dataset.tsv"),
    Path("test-output/cache"),
)

my_node = get_node(taxdump_tree.tree, 2792576)

genetic_code_id, mitochondrial_genetic_code_id = taxdump_tree.get_genetic_codes(2792576)


