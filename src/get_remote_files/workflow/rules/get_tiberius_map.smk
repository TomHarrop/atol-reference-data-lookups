rule get_tiberius_map:
    output:
        tiberius_map="resources/tiberius_map.tsv.gz",
        timestamp="resources/tiberius_map.TIMESTAMP",
    log:
        "resources/get_tiberius_map.log",
    shell:
        "get-tiberius-map "
        "2> {log} "
        "| gzip -9 > {output.tiberius_map} "
        "&& "
        "printf $(date -Iseconds) > {output.timestamp}"
