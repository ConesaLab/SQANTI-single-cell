import os
import sys
import subprocess
from paths import utilitiesPath
    
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
                if args.skipORF:
                    flags.append("--skipORF")
                if args.CAGE_peak:
                    flags.append("--CAGE_peak")
                if args.polyA_motif_list:
                    flags.append("--polyA_motif_list")
                cell_summary_file = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
                if os.path.isfile(cell_summary_file):
                    flags.extend(["--cell_summary", cell_summary_file])

                cmd = (
                    f"Rscript {utilitiesPath}/SQANTI-sc_reads.R "
                    f"{class_file} {junc_file} {args.report} {outputPathPrefix} "
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
    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        cell_summary = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
        if os.path.isfile(cell_summary):
            cell_summaries.append(os.path.abspath(cell_summary))
        else:
            print(f"[INFO] Cell summary not found for {file_acc} ({sampleID}). Skipping this sample.",
                  file=sys.stdout)

    if len(cell_summaries) < 2:
        print("[INFO] Fewer than 2 cell summaries found. Skipping multisample report.",
              file=sys.stdout)
        return

    prefix = getattr(args, 'multisample_report_prefix', 'SQANTI_sc_multisample_report')
    files_arg = ",".join(cell_summaries)
    out_dir = os.path.abspath(args.out_dir)
    report_fmt = args.report
    mode = args.mode

    cmd = (
        f"Rscript {utilitiesPath}/SQANTI-sc_multisample.R "
        f"--files \"{files_arg}\" --out_dir \"{out_dir}\" "
        f"--mode {mode} --report {report_fmt} --prefix \"{prefix}\""
    )

    print("**** Generating multisample SQANTI-sc report...", file=sys.stdout)
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("**** Multisample SQANTI-sc report generated.", file=sys.stdout)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Multisample report generation failed: {e}", file=sys.stderr)


