## atol-reference-data-lookups

Use the NCBI Taxdump files to look up reference databases for the AToL genome
assembly process.

### Reference data

Download the reference data by running `get-remote-files`. Files will be
downloaded to the `resources` directory. This is hard-coded.

### Plan

 - download the taxdump if it isn't cached
 - load the tree into memory
 - run through the taxids once and do all the lookups