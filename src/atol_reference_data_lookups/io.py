import tarfile

from atol_reference_data_lookups import logger


def _extract_tarfile(file_path):
    with tarfile.open(file_path, "r:gz") as tar:
        for member in tar.getmembers():
            if member.isfile() and not member.name.startswith("."):
                for line in tar.extractfile(member).read().decode().splitlines():
                    yield line


def read_gzip_textfile(file_path):
    file_string = file_path.as_posix()
    if file_string.endswith(".tar.gz") or file_string.endswith(".tgz"):
        f = _extract_tarfile(file_path)
    else:
        import gzip

        f = gzip.open(file_path, "rt")

    for i, line in enumerate(f, 1):
        if "\x00" in line:
            raise ValueError(f"Null bytes at line {i} of {file_path}")
        yield line


def read_busco_mapping(taxids_to_busco_dataset_mapping):
    dataset_mapping = read_gzip_textfile(taxids_to_busco_dataset_mapping)
    next(dataset_mapping)  # skip the header
    taxid_to_dataset = {}
    for mapping in dataset_mapping:
        splits = mapping.strip().split(maxsplit=1)
        taxid_to_dataset[int(splits[0])] = str(splits[1])
    logger.debug(taxid_to_dataset)
    return taxid_to_dataset


def read_augustus_mapping(taxids_to_augustus_dataset_mapping):
    taxid_to_dataset = {}
    with open(taxids_to_augustus_dataset_mapping, "rt") as f:
        for line in f:
            splits = line.strip().split(maxsplit=1)
            taxid_to_dataset[int(splits[0])] = str(splits[1])
    logger.debug(taxid_to_dataset)
    return taxid_to_dataset


def read_oatk_mapping(oatk_taxid_file):
    taxid_to_hmm = {}
    with open(oatk_taxid_file, "rt") as f:
        for line in f:
            splits = line.strip().split(maxsplit=2)
            taxid_to_hmm[int(splits[1])] = str(splits[0])
    logger.debug(taxid_to_hmm)
    return taxid_to_hmm
