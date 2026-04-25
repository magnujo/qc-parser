from email import utils
from zipfile import Path, ZipFile
import shutil


def check_zip_contains_file(input_file_path, expected_file):
  with  ZipFile(input_file_path, "r") as zf:
        assert expected_file in zf.namelist(), f"{expected_file} not found in ZIP"

def unzip_file(input_file_path, output_dir):
    with ZipFile(input_file_path, 'r') as zip_ref:
        zip_ref.extractall(output_dir)
    
def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()  

def remove_directory(dir_path):
    if dir_path.is_dir():
        shutil.rmtree(dir_path)
        
def is_numeric(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def extract_path_metadata(zip_path: str) -> dict:
    """
    Expected structure:
    /datasets/caeg_production/libraires/lv7/008/003/{libid}/{date}_{fc}/{version}/{hash}/
        stats/reads/fastqc/{data_type}/{filename}_fastqc.zip
    """
    parts = Path(zip_path).parts
    # parts[0] = '/'
    # parts[1] = 'datasets'
    # parts[2] = 'caeg_production'
    # parts[3] = 'libraires'
    # parts[4] = 'lv7'    \
    # parts[5] = '008'     } sharding dirs — ignored
    # parts[6] = '003'    /
    # parts[7] = libid
    # parts[8] = date_fc  e.g. '20231015_HV3TWDSX7'
    # parts[9] = version  e.g. 'v1.08'
    # parts[10]= hash
    # parts[11]= 'stats'
    # parts[12]= 'reads'
    # parts[13]= 'fastqc'
    # parts[14]= data_type  e.g. 'trim', 'raw'
    # parts[15]= filename   e.g. 'Lib_LV7001856478_L004_singleton_fastqc.zip'

    libid            = parts[7]
    date_fc          = parts[8]
    pipeline_version = parts[9]
    pipeline_hash    = parts[10]
    data_type        = parts[14]
    filename         = parts[15]

    # Split 'date_fc' → date and flowcell on the FIRST underscore only,
    # since flowcell IDs can also contain underscores
    date_str, _, flowcell = date_fc.partition("_")
    try:
        from datetime import datetime
        run_date = datetime.strptime(date_str, "%Y%m%d").date()
    except ValueError:
        run_date = None  # handle unexpected formats gracefully

    # Parse filename: Lib_{libid}_{lane}_{subtype}_fastqc.zip
    fname_match = re.match(
        r"Lib_(?P<lib>[^_]+)_(?P<lane>L\d+)_(?P<subtype>.+)_fastqc\.zip",
        filename
    )
    if fname_match:
        lane    = fname_match.group("lane")
        subtype = fname_match.group("subtype")   # R1, R2, singleton, etc.
    else:
        lane    = None
        subtype = None

    return {
        "libid":            libid,
        "run_date":         run_date,
        "flowcell":         flowcell,
        "pipeline_version": pipeline_version,
        "pipeline_hash":    pipeline_hash,
        "data_type":        data_type,
        "lane":             lane,
        "subtype":          subtype,
    }
    
