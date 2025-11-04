import os
import sys
import pandas as pd

def write_gene_counts_by_cell(args, df):
    structural_categories = [
        'antisense','full-splice_match','fusion','genic','genic_intron',
        'incomplete-splice_match','intergenic','novel_in_catalog','novel_not_in_catalog'
    ]

    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping gene counts.",
                  file=sys.stdout)
            continue

        try:
            needed_cols = ['associated_gene', 'structural_category', 'jxnHash', 'CB']
            class_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                   usecols=lambda c: c in needed_cols or c == 'isoform')
            for col in ['associated_gene', 'structural_category', 'CB']:
                if col in class_df.columns:
                    try:
                        class_df[col] = class_df[col].astype('category')
                    except Exception:
                        pass
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        for col in ['associated_gene', 'structural_category', 'CB']:
            if col not in class_df.columns:
                class_df[col] = ''
        if 'jxnHash' not in class_df.columns:
            class_df['jxnHash'] = class_df.get('isoform', pd.Series(index=class_df.index, dtype=str))

        filtered = class_df[(class_df['CB'].notna()) & (class_df['CB'] != '')].copy()
        if filtered.empty:
            out_file = f"{outputPathPrefix}_gene_counts.csv"
            empty_cols = (['associated_gene', 'CB'] + structural_categories +
                          ['total_read_count', 'unique_jxnHash_counts', 'flag_annotated_gene'])
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB rows for {file_acc}. Wrote empty gene counts.", file=sys.stdout)
            continue

        grouped_size = (
            filtered
            .groupby(['associated_gene', 'CB', 'structural_category'], observed=True)
            .size()
            .unstack(fill_value=0)
        )

        for cat in structural_categories:
            if cat not in grouped_size.columns:
                grouped_size[cat] = 0
        grouped_size = grouped_size[structural_categories]

        total_counts = (
            filtered
            .groupby(['associated_gene', 'CB'], observed=True)
            .size()
            .rename('total_read_count')
        )

        unique_jxn_counts = (
            filtered
            .groupby(['associated_gene', 'CB'], observed=True)['jxnHash']
            .nunique()
            .rename('unique_jxnHash_counts')
        )

        agg_df = grouped_size.merge(total_counts, left_index=True, right_index=True)
        agg_df = agg_df.merge(unique_jxn_counts, left_index=True, right_index=True)
        agg_df = agg_df.reset_index()

        agg_df['flag_annotated_gene'] = agg_df['associated_gene'].apply(
            lambda g: 0 if isinstance(g, str) and g.startswith('novel') else 1
        )

        final_cols = (['associated_gene', 'CB'] + structural_categories +
                      ['total_read_count', 'unique_jxnHash_counts', 'flag_annotated_gene'])
        agg_df = agg_df[final_cols]

        out_file = f"{outputPathPrefix}_gene_counts.csv"
        try:
            agg_df.to_csv(out_file, index=False)
            print(f"**** Gene counts by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


def write_ujc_counts_by_cell(args, df):
    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping UJC counts.",
                  file=sys.stdout)
            continue

        try:
            needed_cls_cols = ['isoform', 'associated_gene', 'structural_category', 'exons', 'jxnHash', 'CB']
            cls_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                 usecols=lambda c: c in needed_cls_cols)
            for col in ['associated_gene', 'structural_category', 'CB']:
                if col in cls_df.columns:
                    try:
                        cls_df[col] = cls_df[col].astype('category')
                    except Exception:
                        pass
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        for col in ['associated_gene', 'structural_category', 'CB', 'jxnHash']:
            if col not in cls_df.columns:
                cls_df[col] = ''
        if 'exons' not in cls_df.columns:
            cls_df['exons'] = ''
        cls_df['exons_num'] = pd.to_numeric(cls_df['exons'], errors='coerce')
        cls_df = cls_df[(cls_df['CB'].notna()) & (cls_df['CB'] != '') & (cls_df['jxnHash'].notna()) & (cls_df['jxnHash'] != '')].copy()

        if cls_df.empty:
            out_file = f"{outputPathPrefix}_ujc_counts.csv"
            empty_cols = ['jxnHash', 'CB', 'read_count', 'associated_gene',
                          'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical',
                          'structural_category', 'flag_MEI', 'flag_annotated_gene']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB/jxnHash rows for {file_acc}. Wrote empty UJC counts.", file=sys.stdout)
            continue

        def majority_value(series):
            counts = series.value_counts(dropna=True)
            return counts.idxmax() if len(counts) > 0 else ''

        def flag_mei(series):
            return 1 if (series.fillna(0).astype(float) > 1).any() else 0

        multi_df = cls_df[cls_df['exons_num'] > 1].copy()
        mono_df = cls_df[cls_df['exons_num'] <= 1].copy()

        base_agg_multi = multi_df.groupby(['jxnHash', 'CB'], observed=True).agg(
            read_count=('isoform', 'size'),
            associated_gene=('associated_gene', majority_value),
            structural_category=('structural_category', majority_value),
            flag_MEI=('exons_num', flag_mei)
        ).reset_index()

        if not mono_df.empty:
            base_agg_mono = mono_df.groupby(['jxnHash', 'CB', 'associated_gene', 'structural_category'], observed=True).agg(
                read_count=('isoform', 'size')
            ).reset_index()
            base_agg_mono['flag_MEI'] = 0
        else:
            base_agg_mono = pd.DataFrame(columns=['jxnHash','CB','associated_gene','structural_category','read_count','flag_MEI'])

        jxn_counts = None
        if os.path.isfile(junc_file):
            try:
                junc_needed = ['isoform', 'junction_category', 'canonical']
                jdf = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False,
                                  usecols=lambda c: c in junc_needed)
                iso_to_hash = cls_df[['isoform', 'jxnHash']].drop_duplicates()
                jdf = pd.merge(jdf, iso_to_hash, on='isoform', how='left')
                jdf = jdf[jdf['jxnHash'].notna() & (jdf['jxnHash'] != '')]

                def _type(row):
                    if row['junction_category'] == 'known' and row['canonical'] == 'canonical':
                        return 'known_canonical'
                    if row['junction_category'] == 'known' and row['canonical'] == 'non_canonical':
                        return 'known_non_canonical'
                    if row['junction_category'] == 'novel' and row['canonical'] == 'canonical':
                        return 'novel_canonical'
                    if row['junction_category'] == 'novel' and row['canonical'] == 'non_canonical':
                        return 'novel_non_canonical'
                    return None

                jdf['jxn_type'] = jdf.apply(_type, axis=1)
                jdf = jdf[jdf['jxn_type'].notna()]
                jxn_counts = jdf.groupby(['jxnHash', 'jxn_type']).size().unstack(fill_value=0)
                for col in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']:
                    if col not in jxn_counts.columns:
                        jxn_counts[col] = 0
                jxn_counts = jxn_counts[['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']]
                jxn_counts = jxn_counts.reset_index()
            except Exception as e:
                print(f"[WARNING] Failed to compute junction counts for {file_acc}: {e}", file=sys.stderr)

        if jxn_counts is not None and not base_agg_multi.empty:
            out_multi = pd.merge(base_agg_multi, jxn_counts, on='jxnHash', how='left')
            for col in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']:
                if col not in out_multi.columns:
                    out_multi[col] = 0
                out_multi[col] = out_multi[col].fillna(0).astype(int)
        else:
            out_multi = base_agg_multi.copy()
            out_multi['known_canonical'] = 0
            out_multi['known_non_canonical'] = 0
            out_multi['novel_canonical'] = 0
            out_multi['novel_non_canonical'] = 0

        if not base_agg_mono.empty:
            out_mono = base_agg_mono.copy()
            out_mono['known_canonical'] = 0
            out_mono['known_non_canonical'] = 0
            out_mono['novel_canonical'] = 0
            out_mono['novel_non_canonical'] = 0
        else:
            out_mono = base_agg_mono

        required_cols = ['jxnHash','CB','read_count','associated_gene','structural_category','flag_MEI',
                         'known_canonical','known_non_canonical','novel_canonical','novel_non_canonical']
        for df_norm in (out_multi, out_mono):
            for col in required_cols:
                if col not in df_norm.columns:
                    df_norm[col] = 0 if col in ['read_count','flag_MEI','known_canonical','known_non_canonical','novel_canonical','novel_non_canonical'] else ''
        out_df = pd.concat([out_multi[required_cols], out_mono[required_cols]], ignore_index=True)

        out_df['flag_MEI'] = out_df.groupby(['associated_gene','CB'])['read_count'] \
            .transform(lambda s: (s == s.max()).astype(int))

        out_df['flag_annotated_gene'] = out_df['associated_gene'].apply(
            lambda g: 0 if isinstance(g, str) and g.startswith('novel') else 1
        )

        final_cols = ['jxnHash', 'CB', 'read_count', 'associated_gene',
                      'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical',
                      'structural_category', 'flag_MEI', 'flag_annotated_gene']
        out_df = out_df[final_cols]

        out_file = f"{outputPathPrefix}_ujc_counts.csv"
        try:
            out_df.to_csv(out_file, index=False)
            print(f"**** UJC counts by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


def write_cv_by_cell(args, df):
    import numpy as np

    def safe_to_numeric(series):
        return pd.to_numeric(series, errors='coerce')

    def aggregate_cv(df_grouped, coord_col, abs_diff_col):
        agg = df_grouped.agg(
            mean_abs_diff=(abs_diff_col, 'mean'),
            std_abs_diff=(abs_diff_col, 'std'),
            count=(abs_diff_col, 'size'),
            associated_gene=('associated_gene', lambda s: '|'.join(sorted(set(s.dropna())))),
            gene_count=('associated_gene', lambda s: len(set([g for g in s.dropna()])))
        ).reset_index()

        group_keys = ['chrom', 'strand', coord_col, 'CB']
        cv_vals = df_grouped[abs_diff_col].apply(
            lambda s: np.nan if pd.isna(s.mean()) or s.mean() == 0 else float(np.std(s, ddof=0)) / float(s.mean())
        ).reset_index(name='cv')
        agg = pd.merge(agg, cv_vals, on=group_keys, how='left')

        agg['flag_mean_0'] = agg['mean_abs_diff'].apply(lambda x: 1 if pd.notna(x) and x == 0 else 0)
        agg['flag_std_0'] = agg['std_abs_diff'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
        agg['flag_ref_match'] = agg['flag_mean_0']
        agg['flag_cv_0'] = agg.apply(lambda r: 1 if pd.notna(r['cv']) and r['cv'] == 0 and r['flag_mean_0'] != 1 else 0, axis=1)
        agg['flag_cv_gt_0'] = agg['cv'].apply(lambda x: 1 if pd.notna(x) and x > 0 else 0)
        agg['flag_multi_gene'] = agg['gene_count'].apply(lambda x: 1 if x > 1 else 0)

        def is_annotated(genestr):
            genes = [g for g in str(genestr).split('|') if g]
            if len(genes) == 0:
                return 0
            return 1 if all(not str(g).startswith('novel') for g in genes) else 0

        agg['flag_annotated_gene'] = agg['associated_gene'].apply(is_annotated)

        agg = agg.rename(columns={coord_col: 'coord'})
        return agg

    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping CV.", file=sys.stdout)
            continue
        if not os.path.isfile(junc_file):
            print(f"[INFO] Junctions file not found for {file_acc}. Skipping CV.", file=sys.stdout)
            continue

        try:
            class_needed = ['isoform', 'associated_gene', 'CB']
            cls_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                 usecols=lambda c: c in class_needed)
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        try:
            j_needed = ['isoform', 'chrom', 'strand', 'genomic_start_coord', 'genomic_end_coord',
                        'diff_to_Ref_start_site', 'diff_to_Ref_end_site']
            jdf = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False,
                              usecols=lambda c: c in j_needed)
        except Exception as e:
            print(f"[ERROR] Failed reading {junc_file}: {e}", file=sys.stderr)
            continue

        jdf = pd.merge(jdf, cls_df, on='isoform', how='left')
        jdf = jdf[(jdf['CB'].notna()) & (jdf['CB'] != '')].copy()
        if jdf.empty:
            out_file = f"{outputPathPrefix}_cv.csv"
            empty_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                          'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                          'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                          'flag_cv_0', 'flag_cv_gt_0', 'CB']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB rows for {file_acc}. Wrote empty CV table.", file=sys.stdout)
            continue

        jdf['genomic_start_coord'] = safe_to_numeric(jdf['genomic_start_coord'])
        jdf['genomic_end_coord'] = safe_to_numeric(jdf['genomic_end_coord'])
        jdf['diff_to_Ref_start_site'] = safe_to_numeric(jdf['diff_to_Ref_start_site'])
        jdf['diff_to_Ref_end_site'] = safe_to_numeric(jdf['diff_to_Ref_end_site'])

        jdf = jdf.dropna(subset=['diff_to_Ref_start_site', 'diff_to_Ref_end_site'])
        if jdf.empty:
            out_file = f"{outputPathPrefix}_cv.csv"
            empty_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                          'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                          'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                          'flag_cv_0', 'flag_cv_gt_0', 'CB']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No diff-to-ref data for {file_acc}. Wrote empty CV table.", file=sys.stdout)
            continue

        jdf['ref_junction_start'] = jdf['genomic_start_coord'] + jdf['diff_to_Ref_start_site']
        jdf['ref_junction_end'] = jdf['genomic_end_coord'] + jdf['diff_to_Ref_end_site']
        jdf['abs_diff_start'] = jdf['diff_to_Ref_start_site'].abs()
        jdf['abs_diff_end'] = jdf['diff_to_Ref_end_site'].abs()

        start_group = jdf.groupby(['chrom', 'strand', 'ref_junction_start', 'CB'], dropna=False, observed=True)
        end_group = jdf.groupby(['chrom', 'strand', 'ref_junction_end', 'CB'], dropna=False, observed=True)

        cv_start = aggregate_cv(start_group, 'ref_junction_start', 'abs_diff_start')
        cv_end = aggregate_cv(end_group, 'ref_junction_end', 'abs_diff_end')

        cv_start['flag_donor'] = cv_start['strand'].apply(lambda s: 1 if s == '+' else 0)
        cv_start['flag_acceptor'] = cv_start['strand'].apply(lambda s: 1 if s == '-' else 0)
        cv_end['flag_donor'] = cv_end['strand'].apply(lambda s: 1 if s == '-' else 0)
        cv_end['flag_acceptor'] = cv_end['strand'].apply(lambda s: 1 if s == '+' else 0)

        out_df = pd.concat([cv_start, cv_end], ignore_index=True)

        final_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                      'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                      'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                      'flag_cv_0', 'flag_cv_gt_0', 'CB']
        for col in final_cols:
            if col not in out_df.columns:
                out_df[col] = np.nan if col in ['mean_abs_diff', 'std_abs_diff', 'cv'] else 0
        out_df = out_df[final_cols]

        out_df['coord'] = pd.to_numeric(out_df['coord'], errors='coerce').astype('Int64')

        out_file = f"{outputPathPrefix}_cv.csv"
        try:
            out_df.to_csv(out_file, index=False)
            print(f"**** CV by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


