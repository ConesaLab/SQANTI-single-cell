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





    for index, row in df.iterrows():
        # Validate 'abundance' if present
        has_abundance = False
        if 'abundance' in df.columns:
            ab = row.get('abundance')
            if pd.notna(ab) and ab:
                if os.path.isdir(str(ab)):
                     df.loc[index, 'abundance'] = os.path.abspath(str(ab))
                     has_abundance = True
                else:
                     print(f"[WARNING] Provided abundance directory not found: {ab}", file=sys.stderr)

        # Initial check for 'cell_association'
        has_assoc = False
        if 'cell_association' in df.columns:
             ca = row.get('cell_association')
             if pd.notna(ca) and ca:
                 if os.path.isfile(str(ca)):
                      df.loc[index, 'cell_association'] = os.path.abspath(str(ca))
                      has_assoc = True
                 else:
                      print(f"[WARNING] Provided cell_association file not found: {ca}", file=sys.stderr)
        
        # Only auto-discover cell_association if missing AND it is strictly required
        needs_assoc = (getattr(args, 'mode', '') == 'reads') or (getattr(args, 'mode', '') == 'isoforms' and not has_abundance)
        
        if not has_assoc and needs_assoc:
            file_acc = row['file_acc']
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            assoc_file_path = os.path.join(args.input_dir, f"{file_acc}_cell_association.txt")

            if bam_files:
                df.loc[index, 'cell_association'] = os.path.abspath(bam_files[0])
                has_assoc = True
            elif os.path.isfile(assoc_file_path):
                df.loc[index, 'cell_association'] = os.path.abspath(assoc_file_path)
                has_assoc = True

        # Handle coverage
        if 'coverage' in df.columns:
            cov = row.get('coverage')
            if pd.notna(cov) and cov:
                if os.path.isfile(str(cov)):
                     df.loc[index, 'coverage'] = os.path.abspath(str(cov))

        # Handle SR_bam
        if 'SR_bam' in df.columns:
            srb = row.get('SR_bam')
            if pd.notna(srb) and srb:
                if os.path.isfile(str(srb)):
                     df.loc[index, 'SR_bam'] = os.path.abspath(str(srb))

        # Handle input_file
        if 'input_file' in df.columns:
            inf = row.get('input_file')
            if pd.notna(inf) and inf:
                if os.path.isfile(str(inf)):
                     df.loc[index, 'input_file'] = os.path.abspath(str(inf))
                else:
                     print(f"[WARNING] Provided input_file not found: {inf}", file=sys.stderr)
        
        # Validation based on mode
        final_has_assoc = 'cell_association' in df.columns and pd.notna(df.loc[index, 'cell_association']) and str(df.loc[index, 'cell_association']).strip() != ''
        final_has_abundance = 'abundance' in df.columns and pd.notna(df.loc[index, 'abundance']) and str(df.loc[index, 'abundance']).strip() != ''

        if args.mode == 'reads':
            if not final_has_assoc:
                raise ValueError(
                    f"ERROR: 'cell_association' is required for sample {row['sampleID']} in 'reads' mode."
                )
        elif args.mode == 'isoforms':
            if not final_has_assoc and not final_has_abundance:
                raise ValueError(
                    f"ERROR: Either 'cell_association' or 'abundance' is required for sample {row['sampleID']} in 'isoforms' mode."
                )

    # Write the DataFrame back to the CSV file
    df.to_csv(args.inDESIGN, sep=',', index=False)

    return df


