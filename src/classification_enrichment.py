import os
import sys
import re
import subprocess
import pandas as pd
import pysam
import hashlib


def annotate_with_ujc_hash(args, df):
    def format_chr(chr_val):
        return f"chr{chr_val}" if not str(chr_val).startswith("chr") else str(chr_val)

    def create_jxn_string(row):
        base = f"{format_chr(row['chr'])}_{row['strand']}_"
        if pd.isna(row["jxn_string"]):
            return f"{base}monoexon_{row['associated_transcript']}"
        junctions = row["jxn_string"].split("_")
        return f"{base}" + "_".join(junctions[2:])

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        print(f"**** Calculating UJCs for {file_acc}...", file=sys.stdout)

        gffread_tmp = f"{outputPathPrefix}_corrected.gtf.gffread_tmp"
        gffread_cmd = f"gffread {outputPathPrefix}_corrected.gtf -T -o {gffread_tmp} && mv {gffread_tmp} {outputPathPrefix}_corrected.gtf"

        try:
            subprocess.check_call(gffread_cmd, shell=True)
        except subprocess.CalledProcessError:
            print(f"ERROR running command: {gffread_cmd}\n Missing or failed gffread", file=sys.stderr)
            sys.exit(-1)

        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.gtf | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.gtf"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,"chr"$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""

        if subprocess.check_call(introns_cmd, shell=True) != 0:
            print(f"ERROR running command: {introns_cmd}\nMissing GTFTOOLS",
                  file=sys.stderr)
            sys.exit(-1)

        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")

        if subprocess.check_call(ujc_cmd, shell=True) != 0:
            print(f"ERROR running command: {ujc_cmd}\nMissing BEDTOOLS",
                  file=sys.stderr)
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}tmp_introns.bed")

        with open(f"{outputPathPrefix}tmp_UJC.txt", 'r') as f_in, \
             open(f"{outputPathPrefix}tmp_UJC_chr.txt", 'w') as f_out:
            for line in f_in:
                parts = line.strip().split('\t')
                chr_part = parts[1].split('_')[0]
                parts[1] = '_'.join([format_chr(chr_part)] + parts[1].split('_')[1:])
                f_out.write('\t'.join(parts) + '\n')
        os.rename(f"{outputPathPrefix}tmp_UJC_chr.txt", f"{outputPathPrefix}tmp_UJC.txt")

        classfile = f"{outputPathPrefix}_classification.txt"
        clas_df = pd.read_csv(
            classfile, sep="\t", usecols=[0, 1, 2, 7], dtype="str"
        )
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(
            f"{outputPathPrefix}tmp_UJC.txt", sep='\t',
            names=["isoform", "jxn_string"], dtype="str"
        )
        merged_df = pd.merge(clas_df, ujc_df, on="isoform", how="left")
        merged_df["jxn_string"] = merged_df.apply(create_jxn_string, axis=1)
        merged_df['jxnHash'] = merged_df['jxn_string'].apply(
            lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest()
        )

        original_class_df = pd.read_csv(classfile, sep='\t', dtype=str, low_memory=False)
        new_cols_df = merged_df[['isoform', 'jxn_string', 'jxnHash']]

        final_df = pd.merge(original_class_df, new_cols_df, on='isoform', how='left')
        final_df.to_csv(f"{outputPathPrefix}_classification_tmp.txt", index=False, sep='\t')
        os.rename(f"{outputPathPrefix}_classification_tmp.txt", classfile)

        print(f"**** UJCs successfully calculated and added to {file_acc} classification", file=sys.stdout)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")


def annotate_with_cell_metadata(args, df):
    """
    Extract isoform, UMI, and cell barcode information from BAM file,
    cell association file, or abundance folder, and add them to classification.
    """
    def _parse_abundance_folder(folder_path):
        cell_dict = {} # isoform -> { "CB": "bc1,bc2", "FL": "cnt1,cnt2" }
        
        matrix_file = os.path.join(folder_path, "matrix.mtx")
        features_file = os.path.join(folder_path, "features.tsv")
        barcodes_file = os.path.join(folder_path, "barcodes.tsv")
        
        if not all(os.path.exists(f) for f in [matrix_file, features_file, barcodes_file]):
             print(f"[ERROR] Abundance folder {folder_path} missing required files (matrix.mtx, features.tsv, barcodes.tsv).", file=sys.stderr)
             return cell_dict

        # Load features (1-based index)
        features = []
        try:
            with open(features_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if parts:
                        features.append(parts[0]) # Isoform ID is 1st col
        except Exception as e:
             print(f"[ERROR] Failed to read features.tsv: {e}", file=sys.stderr)
             return cell_dict

        # Load barcodes (1-based index)
        barcodes = []
        try:
            with open(barcodes_file, 'r') as f:
                for line in f:
                    barcodes.append(line.strip())
        except Exception as e:
             print(f"[ERROR] Failed to read barcodes.tsv: {e}", file=sys.stderr)
             return cell_dict

        # Parse Matrix
        try:
            with open(matrix_file, 'r') as f:
                header_passed = False
                for line in f:
                    if line.startswith('%'):
                        continue
                    if not header_passed:
                        # Dimensions line: rows cols entries
                        header_passed = True
                        continue
                    
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        row_idx = int(parts[0]) - 1 # 0-based
                        col_idx = int(parts[1]) - 1 # 0-based
                        count = parts[2] # Keep as string
                        
                        if row_idx < len(features) and col_idx < len(barcodes):
                            isoform = features[row_idx]
                            barcode = barcodes[col_idx]
                            
                            if isoform not in cell_dict:
                                cell_dict[isoform] = {"CB": [], "FL": []}
                            
                            cell_dict[isoform]["CB"].append(barcode)
                            cell_dict[isoform]["FL"].append(count)
                            
        except Exception as e:
             print(f"[ERROR] Failed to read matrix.mtx: {e}", file=sys.stderr)
             return cell_dict

        # Join lists
        for iso in cell_dict:
            cell_dict[iso]["CB"] = ",".join(cell_dict[iso]["CB"])
            cell_dict[iso]["FL"] = ",".join(cell_dict[iso]["FL"])
            
        return cell_dict

    def _parse_isoform_tsv_for_cell_association(tsv_file):
        cell_dict = {}
        try:
            # First try reading with header
            df_assoc = pd.read_csv(tsv_file, sep='\t', dtype=str, comment='#')
            
            # Identify ID column
            id_col = None
            possible_id_cols = ['pbid', 'isoform_id', 'isoform', 'id']
            for col in df_assoc.columns:
                if col.lower() in possible_id_cols:
                    id_col = col
                    break
            
            # Identify Barcode column
            cb_col = None
            possible_cb_cols = ['cell_barcodes', 'cbs', 'cb']
            for col in df_assoc.columns:
                if col.lower() in possible_cb_cols:
                    cb_col = col
                    break
            
            # Fallback to indices if columns not found by name
            if id_col is None or cb_col is None:
                # Re-read without header to check if it's a headerless 2-col file
                # Or just assume the dataframe we have is using first row as header incorrectly
                # But safer to check dimensions or just use iloc if we assume 1st and 2nd cols
                
                # If we couldn't find headers, let's try to assume it might be headerless
                # or the headers are just weird. 
                # Strategy: Check if we have at least 2 columns.
                if df_assoc.shape[1] >= 2:
                    # If we suspect headerless, the first row we read might be data.
                    # But pandas read_csv defaults to header=0. 
                    # Let's re-read with header=None to be safe if we failed to find headers
                    df_assoc = pd.read_csv(tsv_file, sep='\t', header=None, dtype=str, comment='#')
                    id_col = 0
                    cb_col = 1
                else:
                     print(f"[ERROR] {tsv_file} has fewer than 2 columns and no recognized headers.", file=sys.stderr)
                     return cell_dict

            # Iterate and populate
            # If we re-read with header=None, columns are integers 0, 1
            # If we used original df, columns are strings
            
            for _, row in df_assoc.iterrows():
                pbid = row[id_col]
                cbs = row[cb_col]
                if pd.notna(pbid) and pd.notna(cbs):
                    cell_dict[pbid] = {"CB": cbs}
                    
        except Exception as e:
            print(f"[ERROR] Failed to parse {tsv_file}: {e}", file=sys.stderr)
        return cell_dict

    def _parse_bam_for_cell_association(bam_file):
        cell_dict = {}
        print(f"Extracting CBs and UMIs from {bam_file} ...", file=sys.stderr)
        if not os.path.isfile(bam_file):
            print(f"[ERROR] BAM file {bam_file} not accessible. Skipping...",
                  file=sys.stderr)
            return cell_dict
        try:
            with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
                for read in bam:
                    if read.has_tag("XM") and read.has_tag("CB"):
                        umi = read.get_tag("XM")
                        cb = read.get_tag("CB")
                        isoform = read.query_name
                        cell_dict[isoform] = {"UMI": umi, "CB": cb}
        except Exception as e:
            print(f"[ERROR] Failed to parse BAM file {bam_file}: {e}",
                  file=sys.stderr)
        return cell_dict

    def _parse_tsv_for_cell_association(tsv_file):
        cell_dict = {}
        try:
            df_assoc = pd.read_csv(tsv_file, sep='\t', dtype=str, comment='#')
            required = ['id', 'cell_barcode']
            if not all(col in df_assoc.columns for col in required):
                print(f"[ERROR] {tsv_file} missing required columns.",
                      file=sys.stderr)
                return cell_dict
            for _, row in df_assoc.iterrows():
                iso, cb = row['id'], row['cell_barcode']
                umi = row.get('umi')
                cell_dict[iso] = {"CB": cb}
                if umi:
                    cell_dict[iso]["UMI"] = umi
        except Exception as e:
            print(f"[ERROR] Failed to parse {tsv_file}: {e}", file=sys.stderr)
        return cell_dict

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        
        # Get paths (renamed cell_association)
        assoc_path = row.get('cell_association')
        abundance_path = row.get('abundance')

        cell_data = {}
        
        if args.mode == "isoforms":
            # Priority: Abundance > Cell Association (or use abundance if present)
            if pd.notna(abundance_path) and abundance_path:
                print(f"INFO: Using abundance folder: {abundance_path}", file=sys.stdout)
                cell_data = _parse_abundance_folder(abundance_path)
            elif pd.notna(assoc_path) and assoc_path:
                print(f"INFO: Using cell association file: {assoc_path}", file=sys.stdout)
                cell_data = _parse_isoform_tsv_for_cell_association(assoc_path)
                
        elif args.mode == "reads":
            if pd.notna(assoc_path) and assoc_path:
                print(f"INFO: Using cell association file: {assoc_path}", file=sys.stdout)
                if assoc_path.endswith(".bam"):
                    cell_data = _parse_bam_for_cell_association(assoc_path)
                else:
                    cell_data = _parse_tsv_for_cell_association(assoc_path)

        # Determine classification file path
        if 'classification_file' in row and pd.notna(row['classification_file']) and row['classification_file'] != '':
            class_file = row['classification_file']
        else:
            class_file = f"{outputPathPrefix}_classification.txt"

        if not cell_data:
            print(f"[INFO] No cell data for {file_acc}. Skipping.",
                  file=sys.stdout)
            continue

        if os.path.isfile(class_file):
            try:
                class_df = pd.read_csv(class_file, sep='\t',
                                       low_memory=False, dtype=str)

                if args.mode == "reads":
                    class_df['UMI'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("UMI")
                    )
                    class_df['CB'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )

                    for i, r in class_df.iterrows():
                        match = re.match(r"(.+)_dup\d+$", r['isoform'])
                        if match:
                            primary = match.group(1)
                            if primary in cell_data:
                                class_df.at[i, 'UMI'] = cell_data[primary].get('UMI')
                                class_df.at[i, 'CB'] = cell_data[primary].get('CB')
                                
                elif args.mode == "isoforms":
                    class_df['CB'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )
                    # Add FL column if available (from abundance)
                    class_df['FL'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("FL")
                    )

                for col in class_df.select_dtypes(include=['bool']).columns:
                    class_df[col] = class_df[col].map({True: 'TRUE', False: 'FALSE'})

                class_df = class_df.fillna('NA')
                class_df.to_csv(class_file, index=False, sep="\t")
            except Exception as e:
                print(f"[ERROR] Could not merge cell data: {e}",
                      file=sys.stderr)
        else:
            print(f"[INFO] Classification file for {file_acc} not found at {class_file}.",
                  file=sys.stdout)

        if 'junction_file' in row and pd.notna(row['junction_file']) and row['junction_file'] != '':
            junctions_path = row['junction_file']
        else:
            junctions_path = f"{outputPathPrefix}_junctions.txt"
        if os.path.isfile(junctions_path):
            try:
                junc_df = pd.read_csv(junctions_path, sep='\t', dtype=str)
                if 'isoform' in junc_df.columns:
                    junc_df['CB'] = junc_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )
                    junc_df.to_csv(junctions_path, sep='\t', index=False)
                    print(f"**** CB column added to {file_acc} junctions.",
                          file=sys.stdout)
                else:
                    print(f"[ERROR] 'isoform' col not in {junctions_path}",
                          file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] Could not add CB to junctions: {e}",
                      file=sys.stderr)
        else:
            print(f"[INFO] Junctions file for {file_acc} not found.",
                  file=sys.stdout)

        print(f"**** Cell data added to {file_acc}", file=sys.stdout)



