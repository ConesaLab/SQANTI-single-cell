import sys

import filter_io
from filter_pipeline import run_downstream
from rescue_args import build_rescue_parser
from transcript_rescue import run_transcript_rescue


def main():
    args = build_rescue_parser().parse_args()
    try:
        df = filter_io.read_sample_table(args.inDESIGN, args.qc_dir, labelled_summary_ok=True)
        mode, evidence = run_transcript_rescue(args, df)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
    run_downstream(args, df, mode, evidence)


if __name__ == "__main__":
    main()
