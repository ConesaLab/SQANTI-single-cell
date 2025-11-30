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

    if 'cell_association' not in df.columns:
        df['cell_association'] = ''
    if 'abundance' not in df.columns:
        df['abundance'] = ''

    for index, row in df.iterrows():
        # Handle cell_association
        if pd.isna(row['cell_association']) or not row['cell_association']:
            file_acc = row['file_acc']
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            assoc_file_path = os.path.join(
                args.input_dir, f"{file_acc}_cell_association.txt"
            )

            if bam_files:
                df.at[index, 'cell_association'] = os.path.abspath(bam_files[0])
            elif os.path.isfile(assoc_file_path):
                df.at[index, 'cell_association'] = os.path.abspath(assoc_file_path)
        
        # Validation based on mode
        has_assoc = pd.notna(df.at[index, 'cell_association']) and df.at[index, 'cell_association']
        has_abundance = pd.notna(df.at[index, 'abundance']) and df.at[index, 'abundance']

        if args.mode == 'reads':
            if not has_assoc:
                raise ValueError(
                    f"ERROR: 'cell_association' is required for sample {row['sampleID']} in 'reads' mode."
                )
        elif args.mode == 'isoforms':
            if not has_assoc and not has_abundance:
                raise ValueError(
                    f"ERROR: Either 'cell_association' or 'abundance' is required for sample {row['sampleID']} in 'isoforms' mode."
                )

    # Write the DataFrame back to the CSV file
    df.to_csv(args.inDESIGN, sep=',', index=False)

    return df


