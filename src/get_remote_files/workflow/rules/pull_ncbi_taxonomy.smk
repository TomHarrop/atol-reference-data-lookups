#!/usr/bin/env python3

taxdump_url = "http://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"

taxdump_files = [
    "citations.dmp",
    "delnodes.dmp",
    "division.dmp",
    "excludedfromtype.dmp",
    "fullnamelineage.dmp",
    "gencode.dmp",
    "host.dmp",
    "images.dmp",
    "merged.dmp",
    "names.dmp",
    "nodes.dmp",
    "rankedlineage.dmp",
    "taxidlineage.dmp",
    "typematerial.dmp",
    "typeoftype.dmp",
]


rule expand_taxdump:
    input:
        taxdump="resources/new_taxdump.tar.gz",
    output:
        [Path("resources/new_taxdump", x).as_posix() for x in taxdump_files],
        timestamp="resources/new_taxdump/TIMESTAMP",
    container:
        "docker://debian:stable-20250113"
    params:
        outdir=subpath(output[0], parent=True),
    shell:
        "mkdir -p {params.outdir} && "
        "tar -xzf {input.taxdump} -C {params.outdir} && "
        "printf $(date -Iseconds) > {output.timestamp}"


rule download_taxdump_file:
    output:
        taxdump=temp("resources/new_taxdump.tar.gz"),
    log:
        "resources/download_taxdump_file.log",
    shadow:
        "minimal"
    container:
        "docker://quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10"
    resources:
        runtime=60,
    params:
        url=taxdump_url,
        name=lambda wildcards: Path(taxdump_url).name,
    shell:
        "wget "
        "--retry-connrefused "
        "--waitretry=5 "
        "--read-timeout=20 "
        "--timeout=15 -t 5 "
        "{params.url} "
        "-O {params.name} "
        "&> {log} "
        "&& "
        "wget {params.url}.md5 -O {params.name}.md5 &>> {log} && "
        "md5sum -c {params.name}.md5 &>> {log} && "
        "mv {params.name} {output.taxdump}"
