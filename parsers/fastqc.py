from email import utils
from os import replace
from pathlib import Path
import re
from wsgiref import headers
import pandas as pd
from pyparsing import line
import utils

def is_numeric(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False

def parse_fastqc_file(filepath):
    """Parse a FastQC data file into a dict of modules."""
    current_module = None
    current_status = None
    current_headers = None
    previous_line = None
    rows = []
    current_header_value = None
    result = {}
    header_values = {}
    lines_done = set()
    module_lines_done = set()

    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()            
            line_parts = line.split("\t")

            assert i >= 0, f"Expected at least one line in the file, got empty file: {filepath}"
            
            # Skip empty lines
            assert line, "Unexpected empty line in FastQC data file"
            assert line != '', "Unexpected empty line in FastQC data file"

            if i == 0:
                # Parse program metadata from the first line
                assert line.startswith("##"), f"Expected first line to start with '##', got: {line}"
                assert 'metadata' not in lines_done, "Found multiple metadata lines in FastQC data file"
                assert len(line_parts) == 2, f"Expected 2 columns in the first line, got: {line_parts}"
                program_name = line_parts[0][2:]
                assert program_name == "FastQC", f"Expected program name 'FastQC', got: {program_name}"
                result["metadata"] = {f"{program_name}_version": line_parts[1]}
                result['modules'] = {}
                previous_line = 'metadata'
                lines_done.add(previous_line)

            if i == 1:
                # Special parsing for the first line after metadata, which should be a module header
                assert previous_line == 'metadata', f"Found module header without ending previous metadata on line {i}"
                assert line.startswith(">>") and not line.startswith(">>END_MODULE"), f"Expected module header in line 1, got: {line}"
                current_module, current_status, previous_line, module_lines_done = parse_module_header(line_parts, module_lines_done, i)
                
            if i > 1:
                assert not line.startswith("##"), f"Expected line {i} to not start with '##', got: {line}"
            
                # Parse module headers
                if line.startswith(">>") and not line.startswith(">>END_MODULE"):
                    assert previous_line == 'end_module', f"Found module header without ending previous module on line {i}"
                    current_module, current_status, previous_line, module_lines_done = parse_module_header(line_parts, module_lines_done, i)

                # Parse table headers
                elif line.startswith("#") and not line.startswith("##") and not is_numeric(line_parts[1]):
                    current_headers = [line_parts[0][1:]] + line_parts[1:]
                    assert current_module is not None, f"Found header without a preceding module header on line {i}"
                    assert len(current_headers) > 1, f"Expected at least one header column in line: {line}"
                    assert previous_line in ['module_header', 'header_value'], f"Found header without a preceding module header on line {i}"
                    previous_line = 'header'
                    assert previous_line not in module_lines_done, f"Found multiple header lines for the same module on line {i}"
                    module_lines_done.add(previous_line)
                
                # Parse special header values (lines that start with # and have a numeric value in the second column)
                elif line.startswith("#") and not line.startswith("##") and is_numeric(line_parts[1]):
                    current_header_value = [line_parts[0][1:]] + line_parts[1:]
                    assert current_module is not None, f"Found header without a preceding module header on line {i}"
                    assert len(current_header_value) > 1, f"Expected at least one header column in line: {line}"
                    assert previous_line == 'module_header', f"Found header without a preceding module header on line {i}"
                    assert len(current_header_value) == 2, f"Expected exactly 2 columns in header line with numeric value, got: {current_header_value}"
                    header_values[current_header_value[0]] = current_header_value[1]
                    previous_line = 'header_value'
                    module_lines_done.add(previous_line)

                # Parse end module
                elif line.startswith(">>END_MODULE"):
                    assert current_module is not None, f"Found '>>END_MODULE' without a preceding module header on line {i}"
                    assert current_status is not None, f"Found '>>END_MODULE' without a preceding module status on line {i}"
                    assert previous_line in ['data','module_header'], f"Found '>>END_MODULE' without preceding data on line {i}"
                    result['modules'][f"{current_module}"] = {"status": current_status, "header_values": header_values, "table": {"headers": current_headers, "rows": rows}}
                    current_module = None
                    current_status = None
                    current_headers = None
                    rows = []
                    header_values = {}
                    previous_line = 'end_module'
                    assert previous_line not in module_lines_done, f"Found multiple '>>END_MODULE' lines for the same module on line {i}"
                    module_lines_done = set() # Reset module lines done for the next module
                
                else:
                    assert previous_line in ['header', 'data'], f"Found data row without preceding header or data on line {i}"
                    assert current_module and current_headers and current_status, f"Found data row without a preceding module header and headers on line {i}"
                    assert current_headers is not None, f"Found '>>END_MODULE' without a preceding module headers on line {i}"
                    assert not line.startswith(">>"), f"Expected data line, got module header: {line}"
                    assert not line.startswith("#"), f"Expected data line, got header: {line}"
                    assert len(line_parts) == len(current_headers), f"Expected {len(current_headers)} columns in data row, got: {line_parts}"
                    values = line.split("\t")
                    rows.append(dict(zip(current_headers, values)))
                    previous_line = 'data' 
    
    assert len(result['modules']) == 11, f"Expected 12 modules in the result, got: {len(result['modules'])}"
    assert len(result) == 2, f"Expected result length: 2, got: {len(result)}"
    assert isinstance(result, dict), f"Expected result type: dict, got: {type(result)}"
                 
    return result

def parse_module_header(line_parts, module_lines_done, i):
    assert len(line_parts) == 2, f"Expected 2 columns in module header, got: {line_parts}"
    assert len(line_parts) == 2, f"Expected 2 columns in module header, got: {line_parts}"
    current_module = line_parts[0][2:]
    current_status = line_parts[1]
    assert current_status in ["pass", "warn", "fail"], f"Unexpected module status: {current_status} in line: {line}"
    assert current_module, f"Module name cannot be empty in line: {line}"
    previous_line = 'module_header'               
    assert previous_line not in module_lines_done, f"Found multiple module header lines for the same module on line {i}"
    module_lines_done.add(previous_line)

    return current_module, current_status, previous_line, module_lines_done


def parse_fastqc_zip(zip_path):
    try:
        output_dir = Path("tmp")
        extracted_path = output_dir / zip_path.with_suffix('').name
        expected_file_name = Path("fastqc_data.txt")
        data_path = extracted_path / expected_file_name
        utils.unzip_file(zip_path, 'tmp')
        assert data_path.is_file(), f"Expected file {expected_file_name} not found in {extracted_path}"
        result = parse_fastqc_file(data_path)
        return result
    finally:
        utils.remove_directory(extracted_path)

def parse_fastqc_zips(zip_paths):
    results = []
    for zip_path in zip_paths:
        result = parse_fastqc_zip(zip_path)
        results.append(result)
    return results

def get_basic_stats() -> pd.DataFrame:
    combined_df = pd.DataFrame()
    results = parse_fastqc_zips(Path("data").glob("*.zip"))
    for result in results:
        status_data = [{'Measure': f"{k} status", 'Value': v['status']} for k, v in result['modules'].items()]
        basic_stats = result['modules']['Basic Statistics']['table']['rows']
        header_data = [{'Measure': i, 'Value': j} for k, v in result['modules'].items() for i, j in v['header_values'].items()]
        df = pd.DataFrame(basic_stats + header_data + status_data)
        df = df.set_index('Measure').T.reset_index(drop=True)
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df
