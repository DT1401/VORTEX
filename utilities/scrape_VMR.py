#!/usr/bin/env python3

import os
import argparse
import requests
import pandas as pd
import tempfile
import sys
import re

def download_vmr(url="https://ictv.global/vmr/current"):
    print(f"Downloading latest ICTV VMR from {url}...")
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error downloading VMR: {e}")
        sys.exit(1)
        
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx')
    with open(temp_file.name, 'wb') as f:
        f.write(response.content)
    return temp_file.name

def process_vmr(xlsx_path, output_csv, host_filter=None):
    print("Reading VMR Excel file...")
    try:
        xls = pd.ExcelFile(xlsx_path)
        sheet_names = xls.sheet_names
        
        target_sheet = None
        for sheet in sheet_names:
            if 'VMR MSL' in sheet.upper():
                target_sheet = sheet
                break
        
        if target_sheet:
            print(f"Found target worksheet: '{target_sheet}'")
            df = pd.read_excel(xls, sheet_name=target_sheet)
        else:
            fallback_index = 1 if len(sheet_names) > 1 else 0
            target_sheet = sheet_names[fallback_index]
            print(f"Warning: 'VMR MSL' sheet not found. Falling back to sheet: '{target_sheet}'")
            df = pd.read_excel(xls, sheet_name=target_sheet)
            
    except Exception as e:
        print(f"Error reading Excel file. Ensure 'openpyxl' is installed. {e}")
        sys.exit(1)
        
    accession_col = next((col for col in df.columns if 'GENBANK' in str(col).upper() and 'ACCESSION' in str(col).upper()), None)
            
    if not accession_col:
        print("Error: Could not find a 'Virus GENBANK accession' column in the selected worksheet.")
        sys.exit(1)

    ranks_and_meta = [
        'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 
        'Family', 'Subfamily', 'Genus', 'Species', 
        'Virus name(s)', 'Host source'
    ]
    keep_cols = [c for c in ranks_and_meta if c in df.columns]
    
    # Drop rows where the entire accession cell is blank
    df = df.dropna(subset=[accession_col])
    
    # Optional: Filter by Host source
    if host_filter:
        if 'Host source' in df.columns:
            print(f"Filtering for Host source containing: '{host_filter}' (case-insensitive)...")
            df['Host source'] = df['Host source'].fillna("Unknown")
            df = df[df['Host source'].str.contains(host_filter, case=False, na=False)]
            if df.empty:
                print(f"Warning: No records found matching host '{host_filter}'. Exiting.")
                sys.exit(0)
        else:
            print("Warning: 'Host source' column not found in VMR. Cannot apply filter.")

    acc_pattern = re.compile(r'[A-Za-z]{1,4}_?\d{5,8}')
    df['Genome_List'] = df[accession_col].astype(str).apply(lambda x: acc_pattern.findall(x))
    
    df = df.explode('Genome_List')
    
    df['Genome'] = df['Genome_List'].str.upper()
    df = df.dropna(subset=['Genome'])
    
    final_cols = ['Genome'] + keep_cols
    
    final_df = df[final_cols].drop_duplicates(subset=['Genome'])
    final_df = final_df.fillna("Unknown")
    
    final_df.to_csv(output_csv, index=False)
    print(f"Successfully saved {len(final_df)} metadata records to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Scrape the latest ICTV VMR and generate a VORTEX-compatible metadata CSV.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path (e.g., ictv_metadata.csv)")
    parser.add_argument("--url", default="https://ictv.global/vmr/current", help="VMR Download URL (defaults to the current latest release)")
    parser.add_argument("--host", default=None, help="Optional: Filter by 'Host source' field using a case-insensitive substring match, e.g., 'bacteria' or 'archaea,bacteria'")
    args = parser.parse_args()
    
    temp_xlsx = None
    try:
        temp_xlsx = download_vmr(args.url)
        process_vmr(temp_xlsx, args.output, host_filter=args.host)
    finally:
        if temp_xlsx and os.path.exists(temp_xlsx):
            os.remove(temp_xlsx)

if __name__ == "__main__":
    main()
