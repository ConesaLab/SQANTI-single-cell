import os
import sys
import glob
import pandas as pd

def fill_design_table(args):
    # Read the design CSV file into a DataFrame
    df = pd.read_csv(args.inDESIGN, sep=",")

    # Check for required columns
    required_columns = {'sampleID', 'file_acc'}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"ERROR: Missing required columns: {missing_columns}")

    # Create the new columns based on the existing ones
    df['classification_file'] = df.apply(
        lambda row: os.path.join(
            args.out_dir, row['file_acc'], f"{row['sampleID']}_classification.txt"
        ),
        axis=1
    )
    df['junction_file'] = df.apply(
        lambda row: os.path.join(
            args.out_dir, row['file_acc'], f"{row['sampleID']}_junctions.txt"
        ),
        axis=1
    )

    if 'cell_association_file' not in df.columns:
        df['cell_association_file'] = ''

    for index, row in df.iterrows():
        if pd.isna(row['cell_association_file']) or not row['cell_association_file']:
            file_acc = row['file_acc']
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            assoc_file_path = os.path.join(
                args.input_dir, f"{file_acc}_cell_association.txt"
            )

            if bam_files:
                df.at[index, 'cell_association_file'] = os.path.abspath(bam_files[0])
            elif os.path.isfile(assoc_file_path):
                df.at[index, 'cell_association_file'] = os.path.abspath(
                    assoc_file_path
                )

    # Write the DataFrame back to the CSV file
    df.to_csv(args.inDESIGN, sep=',', index=False)

    return df


