import os
import pandas as pd
import requests
import zipfile

# Paths and constants
metadata_xlsx = "/usr/xtmp/nd141/projects/data/supplementary-data-4-210306.xlsx"
base_url = "https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/"
data_dir = "/usr/xtmp/nd141/projects/data"
output_dir = "./output"  # adjust as needed
os.makedirs(output_dir, exist_ok=True)

def get_replicates_for_tf(tf_name, metadata_file):
    """
    Read the metadata excel to find samples/replicates for a given TF.
    """
    df = pd.read_excel(metadata_file, sheet_name='GEO_GPL19756_GSE147927')
    # Filter to the TF
    df_tf = df[df['Yeast Target Common Name'].str.strip().str.lower() == tf_name.lower()]
    # Extract columns: Replicate, Sample ID
    replicates = df_tf[['Replicate', 'Sample ID']].drop_duplicates()
    return replicates

def download_and_extract_zip(sample_id, dest_folder):
    """
    Download zip file for sample_id if not exists, and extract its contents to dest_folder.
    """
    zip_filename = os.path.join(dest_folder, f"{sample_id}_YEP.zip")
    if not os.path.isfile(zip_filename):
        url = base_url + f"{sample_id}_YEP.zip"
        print(f"Downloading {url} ...")
        r = requests.get(url)
        r.raise_for_status()
        with open(zip_filename, 'wb') as f:
            f.write(r.content)
    else:
        print(f"{zip_filename} already downloaded.")
    
    # Extract
    extract_path = os.path.join(dest_folder, str(sample_id))
    if not os.path.isdir(extract_path):
        print(f"Extracting {zip_filename} to {extract_path}")
        with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
    else:
        print(f"{extract_path} already extracted.")
    
    return extract_path

def load_all_motif_bed_files(extract_path, sample_id):
    """
    Read all motif bed files for a sample from a folder until motif files are missing.
    Return a dict of motif number -> dataframe.
    Filename pattern: {sample_id}_Motif_{motif_num}_FourColor.bed
    """
    motif_dfs = {}
    motif_num = 1
    while True:
        motif_file = os.path.join(extract_path, f"{sample_id}_YEP/{sample_id}_Motif_{motif_num}_FourColor.bed")
        if not os.path.isfile(motif_file):
            # No more motifs
            break
        print(f"Reading motif file: {motif_file}")
        df = pd.read_csv(motif_file, delimiter='\t', header=None)
        df.columns = ['chrom', 'start', 'end', 'col3', 'col4', 'strand']
        motif_dfs[motif_num] = df
        motif_num += 1
    return motif_dfs

def download_tf_chexmix_bed(tf_name, dest_dir):
    """
    Downloads the ChExMix bed file for the TF from GitHub if not already downloaded.
    Returns the local file path.
    """
    tf_name_uc = tf_name.capitalize()  # 'Abf1' etc, assuming file uses capital first letter
    filename = f"{tf_name_uc}_CX.bed"
    url = f"https://github.com/CEGRcode/2021-Rossi_Nature/raw/master/04_ChExMix_Peaks/{filename}"
    out_path = os.path.join(dest_dir, filename)
    
    if not os.path.exists(out_path):
        print(f"Downloading {filename} from {url} ...")
        r = requests.get(url)
        r.raise_for_status()
        with open(out_path, 'wb') as f:
            f.write(r.content)
        print(f"Saved to {out_path}")
    else:
        print(f"{filename} already downloaded.")
    
    return out_path

def load_tf_chexmix_bed(tf_name, data_dir):
    """
    Download (if needed) and load the TF's ChExMix bed file.
    Returns a pandas dataframe with renamed columns for merging.
    """
    os.makedirs(data_dir, exist_ok=True)
    chexmix_bed_file = download_tf_chexmix_bed(tf_name, data_dir)
    
    # Read bed file
    df = pd.read_csv(chexmix_bed_file, sep='\t', header=None)
    # Assign columns (assuming bed format with >= 6 columns)
    df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    # Rename start/end for merging suffixes
    df.rename(columns={'start': 'start_abf1', 'end': 'end_abf1'}, inplace=True)
    return df

def merge_motif_with_chexmix(motif_dfs, chexmix_df):
    """
    For each motif df, merge with chexmix_df on 'chrom' and filter rows where chexmix 'start' lies within motif start and end.
    Returns a dataframe concatenating all motif results with motif numbers as an extra column.
    """
    all_results = []
    for motif_num, motif_df in motif_dfs.items():
        merged = pd.merge(chexmix_df, motif_df, on='chrom', suffixes=('_abf1', '_four'))
        print(motif_df)
        condition = (merged['start_abf1'] >= merged['start']) & (merged['start_abf1'] <= merged['end'])
        filtered = merged[condition]
        # select columns and add motif number
        filtered = filtered[['chrom', 'start_abf1', 'end_abf1', 'score', 'strand_four']]
        filtered['motif'] = motif_num
        all_results.append(filtered)
    if all_results:
        combined = pd.concat(all_results).drop_duplicates()
    else:
        combined = pd.DataFrame(columns=['chrom', 'start_abf1', 'end_abf1', 'score', 'strand_four', 'motif'])
    return combined

def main(tf_name):
    # Step 0: Read metadata and find out replicates
    replicates = get_replicates_for_tf(tf_name, metadata_xlsx)
    print(f"Found {len(replicates)} replicates for TF {tf_name}.")
    
    # Load chexmix for TF (download & read from GitHub)
    chexmix_df = load_tf_chexmix_bed(tf_name, data_dir)
    
    aggregated_results = []
    for idx, row in replicates.iterrows():
        sample_id = str(row['Sample ID'])
        replicate_name = row['Replicate']
        print(f"Processing replicate: {replicate_name} Sample ID: {sample_id}")
        # Step 1 & 2 : Download and extract
        extract_path = download_and_extract_zip(sample_id, data_dir)
        # Step 3: Load all motif files for this replicate
        motif_dfs = load_all_motif_bed_files(extract_path, sample_id)
        if not motif_dfs:
            print(f"No motif files found for sample {sample_id}, skipping.")
            continue
        # Step 4-6: Merge motifs and chexmix and concatenate per replicate
        combined = merge_motif_with_chexmix(motif_dfs, chexmix_df)
        # Add sample and replicate info
        combined['sample_id'] = sample_id
        combined['replicate'] = replicate_name
        aggregated_results.append(combined)
    
    if aggregated_results:
        final_df = pd.concat(aggregated_results).drop_duplicates(subset=['chrom', 'start_abf1', 'end_abf1', 'score'])
    else:
        print("No results found for any replicates.")
        return
    
    final_df.columns = ['chr', 'start', 'end', 'peakVal', 'strand', 'motif', 'sample_id', 'replicate']
    # Save results    
    output_bed_file = os.path.join(output_dir, f"{tf_name.lower()}_rossi_peak_w_strand.bed")
    print(f"Saving final results to {output_bed_file}")
    final_df.to_csv(output_bed_file, sep='\t', index=False)
    print("Pipeline finished.")

# Example to run:
# if __name__ == "__main__":
#     tf = "Abf1"
#     main(tf)


def run_all_tfs(metadata_file):
    # Read the metadata Excel sheet
    df = pd.read_excel(metadata_file, sheet_name='GEO_GPL19756_GSE147927')
    # Get unique TF names (dropna to avoid null entries)
    unique_tfs = df['Yeast Target Common Name'].dropna().unique()
    
    for tf_name in unique_tfs:
        print(f"\nRunning pipeline for TF: {tf_name}")
        try:
            main(tf_name)
        except Exception as e:
            print(f"Error processing {tf_name}: {e}")

# Then call
run_all_tfs(metadata_xlsx)



