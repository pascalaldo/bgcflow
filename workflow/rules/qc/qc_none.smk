checkpoint qc:
    output: "data/processed/qc/qc_passed.dummy",
    shell:
        """
            touch {output}
        """
