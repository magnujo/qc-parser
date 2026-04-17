import csv
import glob
import os



def parse_files(file_paths, output_tsv_path):
    # First pass: collect all metrics to build column headers
    metrics = []
    print('collecting metrics')
    for filepath in file_paths:
        with open(filepath) as f:
            for line in f:
                if line.startswith("SN\t"):
                    metric = line.strip().split("\t")[1].rstrip(":")
                    if metric not in metrics:
                        metrics.append(metric)

    # Second pass: extract values per file
    print('writing summary')
    with open(output_tsv_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["filename"] + metrics)

        for filepath in file_paths:
            row = {m: "" for m in metrics}
            with open(filepath) as f:
                for line in f:
                    if line.startswith("SN\t"):
                        parts = line.strip().split("\t")
                        row[parts[1].rstrip(":")] = parts[2]
            writer.writerow([os.path.basename(filepath)] + [row[m] for m in metrics])
    print('done')