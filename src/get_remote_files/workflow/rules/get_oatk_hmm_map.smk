hmm_map_url = (
    "https://raw.githubusercontent.com/c-zhou/OatkDB/refs/heads/main/v20230921/TAXID"
)


rule get_oatk_hmm_map:
    output:
        taxid_file=("resources/" "oatk.TAXID.tsv"),
        timestamp="resources/oatk.TAXID.TIMESTAMP",
    log:
        "resources/get_oatk_hmm_map.log",
    shadow:
        "minimal"
    container:
        "docker://quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10"
    params:
        url=hmm_map_url,
    shell:
        "wget {params.url} -O {output.taxid_file} &> {log}"
        "&& "
        "printf $(date -Iseconds) > {output.timestamp}"
