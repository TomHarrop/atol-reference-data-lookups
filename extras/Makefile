get_remote_files:
	get-remote-files --help

test: test_taxid test_tax_list

test_taxid:
	atol-reference-data-lookups \
		--taxid 172942 \
		--nodes resources/new_taxdump/nodes.dmp \
		--names resources/new_taxdump/names.dmp \
		--taxids_to_busco_dataset_mapping resources/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz \
		> test-output/test_taxid.json

test_tax_list:
	atol-reference-data-lookups \
		--taxid-list test-data/taxid_list.txt \
		--nodes resources/new_taxdump/nodes.dmp \
		--names resources/new_taxdump/names.dmp \
		--taxids_to_busco_dataset_mapping resources/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz \
		> test-output/test_tax_list.json