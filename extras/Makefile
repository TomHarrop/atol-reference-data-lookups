shell = /bin/bash

local_version := $(shell python3 -m setuptools_scm)
tmpdir := $(shell mktemp -d)
branch := $(shell git rev-parse --abbrev-ref HEAD)

changelog: CHANGELOG.md

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

CHANGELOG.md: .git/index
	gitchangelog > $(tmpdir)/CHANGELOG.md
	apptainer exec docker://davidanson/markdownlint-cli2:v0.13.0 \
		markdownlint-cli2 \
		--fix \
		:$(tmpdir)/CHANGELOG.md || true
	mv $(tmpdir)/CHANGELOG.md ./CHANGELOG.md
	git add CHANGELOG.md
	git commit -m 'changelog !cosmetic'
	git push origin $(branch)