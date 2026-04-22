from pathlib import Path


def parse_file_path(file_path):
    path = Path(file_path)
    return path.parts, path.stem, path.suffix

# Parses metadat from the current 20261704 standard path
def parse_path_metadata(file_path):
    parts = parse_file_path(file_path)[0]
    metadata = {}
    library_path = parts[7]
    metadata['date'] = parts[8].split("_")[0]
    flowcell_info = parts[8].split("_")[1]
    metadata['flowcell_id'] = flowcell_info[1:]
    metadata['flowcell_side'] = flowcell_info[0]
    metadata['pipeline_version'] = parts[9]
    metadata['pipeline_hash'] = parts[10]
    library_file = parts[15].split("_")[1]
    metadata['lane'] = parts[15].split("_")[2]
    
    assert library_path == library_file, "Library path and library file do not match"
    metadata['library_id'] = library_path
    
    return metadata


