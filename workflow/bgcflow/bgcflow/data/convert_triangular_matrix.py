"""
Convert a triangular matrix (missing everything above the diagonal, such as those output by mash) to a complete, 
symmetrical matrix
"""
import csv
import logging
from itertools import chain
from os.path import basename, splitext
from sys import argv

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s - %(levelname)s] - %(message)s")
    logger = logging.getLogger()

    triangle_matrix_fp, full_matrix_output_fp = argv[1], argv[2]

    genome_names = []

    with open(triangle_matrix_fp) as file:
        logger.info("Reading triangle matrix")
        for i, line in enumerate(file):
            if i == 0:
                n_genomes = line.strip()
                assert n_genomes.isnumeric(), f"First line should be an integer, got '{n_genomes}'"
                n_genomes = int(n_genomes)
                logger.info(f"Found {n_genomes} genomes. Intializing full matrix")
                full_matrix = [[0] * n_genomes for _ in range(n_genomes)]
                logger.info("Populating full matrix")
                continue

            row_values = line.split("\t")

            genome_path = row_values.pop(0)
            genome_file = basename(genome_path)
            genome_name, _ = splitext(genome_file)
            genome_names.append(genome_name)

            row_distances = [float(val) for val in row_values]

            for j, distance in enumerate(row_distances):
                full_matrix[i - 1][j] = distance
                full_matrix[j][i - 1] = distance
    
    logger.info("Finished populating full matrix. Writing to CSV file")

    with open(full_matrix_output_fp, "w") as outfile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(chain(["genome_id"], genome_names))
        for i, row in enumerate(full_matrix):
            genome_name = genome_names[i]
            csv_writer.writerow(chain([genome_name], row))
    
    logger.info(f"Finished writing full matrix to {full_matrix_output_fp}")
