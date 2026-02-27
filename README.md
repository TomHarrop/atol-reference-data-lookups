## atol-reference-data-lookups

Use the NCBI Taxdump files to look up reference databases for the AToL genome
assembly process.

### Installation

The
[BioContainer](https://quay.io/repository/biocontainers/atol-reference-data-lookups?tab=tags)
is the only supported method of running `atol-reference-data-lookups`.

*e.g.* with Apptainer/Singularity:

```bash
apptainer exec \
  docker://quay.io/biocontainers/atol-reference-data-lookups:0.0.1--pyhdfd78af_0 \
  atol-reference-data-lookups
```

### Usage

`atol-reference-data-lookups` takes a single NCBI TaxId, or a path to a plain
text file containing a list of NCBI TaxIds (one per line).

It prints the results of the lookups to Standard Output in JSON format , *e.g.*

```bash
 $ atol-reference-data-lookups \
        --taxid 172942 \
        --nodes resources/new_taxdump/nodes.dmp \
        --names resources/new_taxdump/names.dmp \
        --taxids_to_busco_dataset_mapping resources/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz

{"172942": {"busco_dataset_name": "sauropsida", "augustus_dataset_name": "Xenopus_tropicalis", "genetic_code_id": 1, "mitochondrial_genetic_code_id": 2}}
```

You also need to provide some reference data. 

> [!TIP] A convenience script to download the reference data is
> [included](README.md#reference-data).


- The **nodes** and **names** files are "nodes.dmp" and "names.dmp" from NCBI's
  [new_taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/). 
- **taxids_to_busco_dataset_mapping** is the "mapping_taxids-busco_dataset_name"
  file from
  [busco-data.ezlab.org/v5/data/placement_files](https://busco-data.ezlab.org/v5/data/placement_files/).
- **taxids_to_augustus_dataset_mapping** is a mapping of Augustus training
  datasets to NCBI TaxID, [shipped with the package](src/atol_reference_data_lookups/config/taxid_to_augustus_dataset.tsv).


```
usage: atol-reference-data-lookups [-h] [--taxid TAXID | --taxid-list TAXID_LIST] --nodes NODES --names NAMES
                                   --taxids_to_busco_dataset_mapping TAXIDS_TO_BUSCO_DATASET_MAPPING
                                   [--taxids_to_augustus_dataset_mapping TAXIDS_TO_AUGUSTUS_DATASET_MAPPING] [--cache_dir CACHE_DIR]

options:
  -h, --help            show this help message and exit

Input:
  --taxid TAXID         A single NCBI TaxId to look up
  --taxid-list TAXID_LIST
                        A file containing a list NCBI TaxIds to look up, one per line.

Reference data:
  --nodes NODES         NCBI nodes.dmp file from taxdump
  --names NAMES         NCBI names.dmp file from taxdump
  --taxids_to_busco_dataset_mapping TAXIDS_TO_BUSCO_DATASET_MAPPING
                        BUSCO placement file from https://busco-data.ezlab.org/v5/data/placement_files/
  --taxids_to_augustus_dataset_mapping TAXIDS_TO_AUGUSTUS_DATASET_MAPPING
                        File that maps Augustus datasets to NCBI TaxIDs. See config/taxid_to_augustus_dataset.tsv

General options:
  --cache_dir CACHE_DIR
                        Directory to cache the NCBI taxonomy after processing
```

### Performance notes

`atol-reference-data-lookups` uses
[`skbio.tree`](https://scikit.bio/docs/latest/tree.html) to import the NCBI
Taxdump and build the query trees. This takes a while, but the results are
automatically cached on the first run to speed up following runs. After the
trees are loaded into memory, searches are fast.

Override the default cache directory with the `--cache_dir` argument.

The cache is automatically invalidated if any of the reference data files
change.

### Reference data

Download the reference data by running `get-remote-files`. Files will be
downloaded to the `./resources` directory. This is hard-coded.

```
atol-reference-data-lookups version 0.1.dev13+gcbdbebe90.d20260226
usage: get-remote-files [-h] [-n]
                        [--parallel_downloads PARALLEL_DOWNLOADS]

options:
  -h, --help            show this help message and exit
  -n                    Dry run
  --parallel_downloads PARALLEL_DOWNLOADS
                        Number of parallel downloads
```

### TODO:

- [ ] Lookup OATK database