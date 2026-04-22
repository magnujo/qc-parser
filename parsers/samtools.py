import csv
from glob import glob
import os


def parse_files(paths_to_library_root_folder, output_tsv_path):
    # First pass: collect all metrics to build column headers
    metrics = []
    full_paths = []

    
    for filepath in paths_to_library_root_folder:
        pattern = f"{filepath}/stats/aligns/samtools_stats/*txt"
        matches = glob(pattern)
        num_of_files = len(matches)
        expected_num_of_files = 1
        
        assert num_of_files == expected_num_of_files, f"Expected exactly {expected_num_of_files} files matching pattern, but found {num_of_files}: {matches}"
        
        filepath = matches[0]
        full_paths.append(filepath)
        
        with open(filepath) as f:
            for line in f:
                if line.startswith("SN\t"):
                    metric = line.strip().split("\t")[1].rstrip(":")
                    if metric not in metrics:
                        metrics.append(metric)

    # Second pass: extract values per file
    with open(output_tsv_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["filename"] + metrics)

        for filepath in full_paths:
            row = {m: "" for m in metrics}
            with open(filepath) as f:
                for line in f:
                    if line.startswith("SN\t"):
                        parts = line.strip().split("\t")
                        row[parts[1].rstrip(":")] = parts[2]
            writer.writerow([filepath] + [row[m] for m in metrics])
    
