import csv
from glob import glob

def parse_files(paths_to_library_root_folder, output_tsv_path):
    metrics_seen = set()
    metrics_ordered = []  # preserves insertion order without O(n) `not in` checks
    all_rows = []         # stores (filepath, {metric: value}) — avoids re-reading files

    for filepath in paths_to_library_root_folder:
        pattern = f"{filepath}/stats/aligns/samtools_stats/*txt"
        matches = glob(pattern)
        num_of_files = len(matches)
        expected_num_of_files = 1

        assert num_of_files == expected_num_of_files, (
            f"Expected exactly {expected_num_of_files} files matching pattern, "
            f"but found {num_of_files}: {matches}"
        )

        filepath = matches[0]
        row = {}

        with open(filepath) as f:
            for line in f:
                if line.startswith("SN\t"):
                    parts = line.split("\t")  # no need to strip before splitting
                    metric = parts[1].rstrip(":")
                    value = parts[2].strip()
                    row[metric] = value
                    if metric not in metrics_seen:
                        metrics_seen.add(metric)
                        metrics_ordered.append(metric)

        all_rows.append((filepath, row))

    # Single write pass
    with open(output_tsv_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["filename"] + metrics_ordered)

        for filepath, row in all_rows:
            writer.writerow([filepath] + [row.get(m, "") for m in metrics_ordered])