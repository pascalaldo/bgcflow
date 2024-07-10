rule custom_samples:
    output:
        samples_file="data/interim/custom_genomes/samples.csv",
    run:
        CUSTOM_SAMPLES[["source", "path"]].to_csv(output.samples_file)