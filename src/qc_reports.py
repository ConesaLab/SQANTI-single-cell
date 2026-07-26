import os
import sys
import subprocess
import pandas as pd
from paths import reportAssetsPath
    
def generate_report(args, df):
    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"
        print(f"**** Generating SQANTI3 report for {file_acc}...",
              file=sys.stdout)

        if os.path.isfile(class_file):
            try:
                flags = []
                if args.ignore_cell_summary:
                    flags.append("--ignore_cell_summary")
                if getattr(args, 'include_ORF', False):
                    flags.append("--include_ORF")
                if args.CAGE_peak:
                    flags.append("--CAGE_peak")
                if args.polyA_motif_list:
                    flags.append("--polyA_motif_list")
                
                cell_summary_file = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
                if os.path.isfile(cell_summary_file):
                    flags.extend(["--cell_summary", cell_summary_file])
                
                # Check for clustering results
                # Clustering is usually one level up from sampleID if run per file_acc
                clustering_file = os.path.join(os.path.dirname(outputPathPrefix), "clustering", "umap_results.csv")
                if os.path.isfile(clustering_file):
                    flags.extend(["--clustering", clustering_file])
                
                if hasattr(args, 'refGTF') and args.refGTF:
                    flags.extend(["--refGTF", f'"{args.refGTF}"'])

                cmd = (
                    f"Rscript {reportAssetsPath}/SQANTI-sc_report.R "
                    f"\"{class_file}\" \"{junc_file}\" {args.report} \"{outputPathPrefix}\" "
                    f"{args.mode} {' '.join(flags)}"
                )
                subprocess.run(cmd, shell=True, check=True)
                print(f"**** SQANTI3 report generated for {file_acc}",
                      file=sys.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error generating report for {class_file}: {e}",
                      file=sys.stdout)
        else:
            print(f"Classification file for {file_acc} not found.",
                  file=sys.stdout)


def generate_multisample_report(args, df):
    try:
        total_samples = df.shape[0]
    except Exception:
        total_samples = 0

    if total_samples < 2:
        print("[INFO] Design has fewer than 2 samples. Skipping multisample report.",
              file=sys.stdout)
        return

    cell_summaries = []
    class_files = []
    color_groups = []
    shape_groups = []
    shade_groups = []

    has_color_col = 'color_group' in df.columns
    has_shape_col = 'shape_group' in df.columns
    has_shade_col = 'shade_group' in df.columns

    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        cell_summary = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
        if os.path.isfile(cell_summary):
            cell_summaries.append(os.path.abspath(cell_summary))
            if has_color_col:
                color_groups.append(str(row.get('color_group', '')))
            if has_shape_col:
                shape_groups.append(str(row.get('shape_group', '')))
            if has_shade_col:
                shade_groups.append(str(row.get('shade_group', '')))
        else:
            print(f"[INFO] Cell summary not found for {file_acc} ({sampleID}). Skipping this sample.",
                  file=sys.stdout)

        class_file = f"{outputPathPrefix}_classification.txt"
        if os.path.isfile(class_file):
            class_files.append(os.path.abspath(class_file))

    if len(cell_summaries) < 2:
        print("[INFO] Fewer than 2 cell summaries found. Skipping multisample report.",
              file=sys.stdout)
        return

    prefix = getattr(args, 'multisample_report_prefix', 'SQANTI_sc_multisample_report')
    files_arg = ",".join(cell_summaries)
    out_dir = os.path.abspath(args.out_dir)
    report_fmt = args.report
    mode = args.mode

    class_files_flag = ""
    if len(class_files) >= 2:
        class_files_flag = f' --class_files "{",".join(class_files)}"'

    pca_features_flag = ""
    pca_features = getattr(args, 'pca_features', None)
    if pca_features:
        pca_features_flag = f' --pca_features "{os.path.abspath(pca_features)}"'

    group_flags = ""
    if color_groups and any(v for v in color_groups):
        group_flags += f' --color_group "{",".join(color_groups)}"'
    if shape_groups and any(v for v in shape_groups):
        group_flags += f' --shape_group "{",".join(shape_groups)}"'
    if shade_groups and any(v for v in shade_groups):
        group_flags += f' --shade_group "{",".join(shade_groups)}"'

    cmd = (
        f"Rscript {reportAssetsPath}/SQANTI-sc_multisample_report.R "
        f"--files \"{files_arg}\" --out_dir \"{out_dir}\" "
        f"--mode {mode} --report {report_fmt} --prefix \"{prefix}\""
        f"{class_files_flag}{pca_features_flag}{group_flags}"
    )

    print("**** Generating multisample SQANTI-sc report...", file=sys.stdout)
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("**** Multisample SQANTI-sc report generated.", file=sys.stdout)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Multisample report generation failed: {e}", file=sys.stderr)


