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
    
