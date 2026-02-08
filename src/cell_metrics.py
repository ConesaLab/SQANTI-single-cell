import os
import sys
import numpy as np
import pandas as pd
import warnings
from pandas.errors import PerformanceWarning

def calculate_metrics_per_cell(args, df):
    structural_categories = [
        'full-splice_match', 'incomplete-splice_match', 'novel_in_catalog',
        'novel_not_in_catalog', 'genic', 'antisense', 'fusion', 'intergenic',
        'genic_intron'
    ]
    cat_to_root = {
        'full-splice_match': 'FULLSPLICEMATCH',
        'incomplete-splice_match': 'INCOMPLETESPLICEMATCH',
        'novel_in_catalog': 'NOVELINCATALOG',
        'novel_not_in_catalog': 'NOVELNOTINCATALOG',
        'genic': 'GENIC', 'antisense': 'ANTISENSE', 'fusion': 'FUSION',
        'intergenic': 'INTERGENIC', 'genic_intron': 'GENICINTRON'
    }
    cat_to_tag = {
        'full-splice_match': 'FSM',
        'incomplete-splice_match': 'ISM',
        'novel_in_catalog': 'NIC',
        'novel_not_in_catalog': 'NNC',
        'genic': 'Genic', 'antisense': 'Antisense', 'fusion': 'Fusion',
        'intergenic': 'Intergenic', 'genic_intron': 'Genic_intron'
    }
    abbr_pairs = [
        ('FSM', 'full-splice_match'),
        ('ISM', 'incomplete-splice_match'),
        ('NIC', 'novel_in_catalog'),
        ('NNC', 'novel_not_in_catalog')
    ]

    def final_count_name(cat):
        tag = cat_to_tag[cat]
        return 'Genic_Genomic' if cat == 'genic' else tag

    def safe_prop(numer, denom):
        numer = numer.astype(float)
        denom = denom.astype(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            arr = np.where(denom > 0, (numer / denom) * 100.0, 0.0)
        return pd.Series(arr, index=numer.index)

    warnings.simplefilter('ignore', PerformanceWarning)

    for _, r in df.iterrows():
        file_acc = r['file_acc']
        sampleID = r['sampleID']
        prefix = os.path.join(args.out_dir, file_acc, sampleID)
        if 'classification_file' in r and pd.notna(r['classification_file']) and r['classification_file'] != '':
            class_file = r['classification_file']
        else:
            class_file = f"{prefix}_classification.txt"
            
        if 'junction_file' in r and pd.notna(r['junction_file']) and r['junction_file'] != '':
            junc_file = r['junction_file']
        else:
            junc_file = f"{prefix}_junctions.txt"

        out_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"

        try:
            cls = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False)
        except Exception:
            continue
        try:
            junc = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False)
        except Exception:
            junc = pd.DataFrame()

        for col in ['CB','isoform','associated_gene','structural_category','exons','length','ref_length',
                    'all_canonical','subcategory','chrom','UMI','jxn_string','associated_transcript',
                    'RTS_stage','predicted_NMD','within_CAGE_peak','polyA_motif_found','perc_A_downstream_TTS','diff_to_gene_TSS','coding','min_cov','ratio_TSS']:
            if col not in cls.columns:
                cls[col] = np.nan
        for c in ['exons','length','ref_length','perc_A_downstream_TTS','diff_to_gene_TSS','min_cov','ratio_TSS']:
            cls[c] = pd.to_numeric(cls[c], errors='coerce')

        if args.mode == 'isoforms':
            cls['CB'] = cls['CB'].fillna('')
            if 'FL' in cls.columns:
                cls['FL'] = cls['FL'].fillna('1')
                cls['CB'] = cls['CB'].str.split(',')
                cls['FL'] = cls['FL'].astype(str).str.split(',')
                cls = cls.explode(['CB', 'FL'])
                cls['_count'] = pd.to_numeric(cls['FL'], errors='coerce').fillna(1)
            else:
                cls = cls.assign(CB=cls['CB'].str.split(',')).explode('CB')
                cls['_count'] = 1
            cls['CB'] = cls['CB'].fillna('')
        else:
            cls['_count'] = 1

        cls_valid = cls[(cls['CB'].notna()) & (cls['CB'] != '')].copy()

        total_reads = cls_valid.groupby('CB')['_count'].sum().rename('total_reads')
        if args.mode == 'isoforms':
            total_umi = pd.Series(dtype=float) # No UMI
        else:
            total_umi = cls_valid.groupby('CB')['UMI'].nunique().rename('total_UMI')
        
        reads_no_mono = cls_valid[cls_valid['exons'] != 1].groupby('CB')['_count'].sum().rename('total_reads_no_monoexon')
        
        if args.mode == 'isoforms':
             summary = pd.DataFrame(total_reads).join(reads_no_mono, how='outer').fillna(0)
        else:
             summary = pd.DataFrame(total_reads).join(total_umi, how='outer').join(reads_no_mono, how='outer').fillna(0)

        cat_counts = cls_valid.groupby(['CB','structural_category'])['_count'].sum().unstack(fill_value=0)
        for c in structural_categories:
            if c not in cat_counts.columns:
                cat_counts[c] = 0
        cat_counts = cat_counts[structural_categories]
        cat_counts.columns = [final_count_name(c) for c in structural_categories]
        summary = summary.join(cat_counts, how='left').fillna(0)
        for c in structural_categories:
            tag = final_count_name(c)
            summary[f"{tag}_prop"] = safe_prop(summary[tag], summary['total_reads'])

        summary['Genes_in_cell'] = cls_valid.groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        if args.mode != 'isoforms':
            summary['UJCs_in_cell'] = cls_valid[cls_valid['exons'] > 1].groupby('CB')['jxn_string'].nunique().reindex(summary.index, fill_value=0)

        mt = cls_valid[cls_valid['chrom'] == 'MT'].groupby('CB')['_count'].sum().rename('MT_reads_count')
        summary = summary.join(mt, how='left').fillna({'MT_reads_count': 0})
        summary['MT_perc'] = safe_prop(summary['MT_reads_count'], summary['total_reads'])

        anno = (~cls_valid['associated_gene'].fillna('').str.startswith('novel'))
        summary['Annotated_genes'] = cls_valid[anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        summary['Novel_genes'] = cls_valid[~anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)

        if not junc.empty:
            if 'CB' not in junc.columns or (junc['CB'].fillna('') == '').all():
                iso_to_cb = cls_valid[['isoform','CB']].dropna().drop_duplicates()
                junc = pd.merge(junc, iso_to_cb, on='isoform', how='left')
            jv = junc[(junc['CB'].notna()) & (junc['CB'] != '')].copy()
            if not jv.empty:
                jv['junction_type'] = jv['junction_category'].astype(str) + '_' + jv['canonical'].astype(str)
                counts = jv.groupby(['CB','junction_type']).size().unstack(fill_value=0)
                for tp in ['known_canonical','known_non_canonical','novel_canonical','novel_non_canonical']:
                    if tp not in counts.columns:
                        counts[tp] = 0
                counts['total_junctions'] = counts.sum(axis=1)
                counts = counts.rename(columns={
                    'known_canonical':'Known_canonical_junctions',
                    'known_non_canonical':'Known_non_canonical_junctions',
                    'novel_canonical':'Novel_canonical_junctions',
                    'novel_non_canonical':'Novel_non_canonical_junctions'
                })
                for src, dst in [
                    ('Known_canonical_junctions','Known_canonical_junctions_prop'),
                    ('Known_non_canonical_junctions','Known_non_canonical_junctions_prop'),
                    ('Novel_canonical_junctions','Novel_canonical_junctions_prop'),
                    ('Novel_non_canonical_junctions','Novel_non_canonical_junctions_prop')]:
                    counts[dst] = safe_prop(counts[src].reindex(counts.index, fill_value=0), counts['total_junctions'])
                summary = summary.join(counts, how='left').fillna(0)
            else:
                summary[['Known_canonical_junctions','Known_non_canonical_junctions','Novel_canonical_junctions','Novel_non_canonical_junctions','total_junctions',
                         'Known_canonical_junctions_prop','Known_non_canonical_junctions_prop','Novel_canonical_junctions_prop','Novel_non_canonical_junctions_prop']] = 0
        else:
            summary[['Known_canonical_junctions','Known_non_canonical_junctions','Novel_canonical_junctions','Novel_non_canonical_junctions','total_junctions',
                     'Known_canonical_junctions_prop','Known_non_canonical_junctions_prop','Novel_canonical_junctions_prop','Novel_non_canonical_junctions_prop']] = 0

        sublevels = {
            'full-splice_match': ['alternative_3end','alternative_3end5end','alternative_5end','reference_match','mono-exon'],
            'incomplete-splice_match': ['3prime_fragment','internal_fragment','5prime_fragment','intron_retention','mono-exon'],
            'novel_in_catalog': ['combination_of_known_junctions','combination_of_known_splicesites','intron_retention','mono-exon_by_intron_retention','mono-exon'],
            'novel_not_in_catalog': ['at_least_one_novel_splicesite','intron_retention'],
            'genic': ['mono-exon','multi-exon'], 'antisense': ['mono-exon','multi-exon'],
            'fusion': ['intron_retention','multi-exon'], 'intergenic': ['mono-exon','multi-exon'],
            'genic_intron': ['mono-exon','multi-exon']
        }
        def _normalize_sub_lv(lv: str) -> str:
            return lv.replace('-', '_')
        def subkey(cat, lv):
            tag = cat_to_tag[cat]
            return f"{tag}_{_normalize_sub_lv(lv)}_prop"
        for cat in structural_categories:
            denom = summary[final_count_name(cat)]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            tbl = sub.groupby(['CB','subcategory'])['_count'].sum().unstack(fill_value=0) if not sub.empty else pd.DataFrame()
            for lv in sublevels[cat]:
                numer = tbl.get(lv, pd.Series(0, index=summary.index)).reindex(summary.index, fill_value=0)
                summary[subkey(cat, lv.replace('-', '_'))] = safe_prop(numer, denom).fillna(0)

        gene_counts = cls_valid.groupby(['CB','associated_gene'])['_count'].sum().rename('read_count').reset_index()
        gene_counts['gene_type'] = np.where(gene_counts['associated_gene'].fillna('').str.startswith('novel'), 'novel', 'annotated')
        bins = gene_counts.groupby(['CB','gene_type']).agg(
            bin1_count=('read_count', lambda s: (s == 1).sum()),
            bin2_4_count=('read_count', lambda s: ((s >= 2) & (s <= 4)).sum()),
            bin5_9_count=('read_count', lambda s: ((s >= 5) & (s <= 9)).sum()),
            bin10plus_count=('read_count', lambda s: (s >= 10).sum()),
            total_genes_in_type=('associated_gene','nunique')
        ).reset_index()
        def bin_props(df, gene_kind, out_prefix):
            out = pd.DataFrame(index=summary.index)
            keyed = df[df['gene_type'] == gene_kind].set_index('CB') if not df.empty else pd.DataFrame(index=summary.index)
            for label, src in [(f"{out_prefix}_bin1_perc", 'bin1_count'), (f"{out_prefix}_bin2_4_perc", 'bin2_4_count'), (f"{out_prefix}_bin5_9_perc", 'bin5_9_count'), (f"{out_prefix}_bin10plus_perc", 'bin10plus_count')]:
                if not keyed.empty and src in keyed.columns:
                    out[label] = safe_prop(keyed[src].reindex(summary.index, fill_value=0), keyed['total_genes_in_type'].reindex(summary.index, fill_value=0)).fillna(0)
                else:
                    out[label] = 0
            return out
        summary = summary.join(bin_props(bins, 'annotated', 'anno')).join(bin_props(bins, 'novel', 'novel'))

        if args.mode != 'isoforms':
            gene_ujc = cls_valid[cls_valid['exons'] > 1].groupby(['CB','associated_gene'])['jxn_string'].nunique().rename('ujc_count').reset_index()
            gene_ujc['gene_type'] = np.where(gene_ujc['associated_gene'].fillna('').str.startswith('novel'), 'novel', 'annotated')
            ujc_bins = gene_ujc.groupby(['CB','gene_type']).agg(
                ujc_bin1_count=('ujc_count', lambda s: (s == 1).sum()),
                ujc_bin2_3_count=('ujc_count', lambda s: ((s >= 2) & (s <= 3)).sum()),
                ujc_bin4_5_count=('ujc_count', lambda s: ((s >= 4) & (s <= 5)).sum()),
                ujc_bin6plus_count=('ujc_count', lambda s: (s >= 6).sum()),
                total_genes_in_type_ujc=('associated_gene','nunique')
            ).reset_index()
            def ujc_props(df, gene_kind, out_prefix):
                out = pd.DataFrame(index=summary.index)
                keyed = df[df['gene_type'] == gene_kind].set_index('CB') if not df.empty else pd.DataFrame(index=summary.index)
                for label, src in [(f"{out_prefix}_ujc_bin1_perc", 'ujc_bin1_count'), (f"{out_prefix}_ujc_bin2_3_perc", 'ujc_bin2_3_count'), (f"{out_prefix}_ujc_bin4_5_perc", 'ujc_bin4_5_count'), (f"{out_prefix}_ujc_bin6plus_perc", 'ujc_bin6plus_count')]:
                    if not keyed.empty and src in keyed.columns:
                        out[label] = safe_prop(keyed[src].reindex(summary.index, fill_value=0), keyed['total_genes_in_type_ujc'].reindex(summary.index, fill_value=0)).fillna(0)
                    else:
                        out[label] = 0
                return out
            summary = summary.join(ujc_props(ujc_bins, 'annotated', 'anno')).join(ujc_props(ujc_bins, 'novel', 'novel'))

        def compute_lenbins_by_cb(df_group):
            df_g = df_group.copy()
            conditions = [
                (df_g['length'] <= 250),
                (df_g['length'] > 250) & (df_g['length'] <= 500),
                (df_g['length'] > 500) & (df_g['length'] <= 1000),
                (df_g['length'] > 1000) & (df_g['length'] <= 2000),
                (df_g['length'] > 2000)
            ]
            choices = ['two_fifty', 'five_hund', 'short', 'mid', 'long']
            df_g['len_cat'] = np.select(conditions, choices, default='unknown')
            
            counts = df_g.groupby(['CB', 'len_cat'])['_count'].sum().unstack(fill_value=0)
            # Ensure all cols exist
            for c in choices:
                if c not in counts.columns:
                    counts[c] = 0
            return counts
        lo = compute_lenbins_by_cb(cls_valid)
        for dst, src in [('Total_250b_length_prop','two_fifty'),('Total_500b_length_prop','five_hund'),('Total_short_length_prop','short'),('Total_mid_length_prop','mid'),('Total_long_length_prop','long')]:
            summary[dst] = safe_prop(lo[src].reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        mono = cls_valid[cls_valid['exons'] == 1]
        if not mono.empty:
            lm = compute_lenbins_by_cb(mono)
            for dst, src in [('Total_250b_length_mono_prop','two_fifty'),('Total_500b_length_mono_prop','five_hund'),('Total_short_length_mono_prop','short'),('Total_mid_length_mono_prop','mid'),('Total_long_length_mono_prop','long')]:
                summary[dst] = safe_prop(lm[src].reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        else:
            for nm in ['Total_250b_length_mono_prop','Total_500b_length_mono_prop','Total_short_length_mono_prop','Total_mid_length_mono_prop','Total_long_length_mono_prop']:
                summary[nm] = 0
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = summary[final_count_name(cat)]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            if sub.empty:
                for nm in ['250b_length_prop','500b_length_prop','short_length_prop','mid_length_prop','long_length_prop']:
                    summary[f"{tag}_{nm}"] = 0
                for nm in ['250b_length_mono_prop','500b_length_mono_prop','short_length_mono_prop','mid_length_mono_prop','long_length_mono_prop']:
                    summary[f"{tag}_{nm}"] = 0
            else:
                lc = compute_lenbins_by_cb(sub)
                for dst, src in [('250b_length_prop','two_fifty'),('500b_length_prop','five_hund'),('short_length_prop','short'),('mid_length_prop','mid'),('long_length_prop','long')]:
                    summary[f"{tag}_{dst}"] = safe_prop(lc[src].reindex(summary.index, fill_value=0), denom).fillna(0)
                subm = sub[sub['exons'] == 1]
                if subm.empty:
                    for nm in ['250b_length_mono_prop','500b_length_mono_prop','short_length_mono_prop','mid_length_mono_prop','long_length_mono_prop']:
                        summary[f"{tag}_{nm}"] = 0
                else:
                    lcm = compute_lenbins_by_cb(subm)
                    for dst, src in [('250b_length_mono_prop','two_fifty'),('500b_length_mono_prop','five_hund'),('short_length_mono_prop','short'),('mid_length_mono_prop','mid'),('long_length_mono_prop','long')]:
                        summary[f"{tag}_{dst}"] = safe_prop(lcm[src].reindex(summary.index, fill_value=0), denom).fillna(0)

        # Reference body coverage: parameterized threshold and export cutoff for plotting
        ref_cov_min = float(getattr(args, 'ref_cov_min_pct', 45.0))
        cls_valid['ref_body_cov_flag'] = (cls_valid['length'] / cls_valid['ref_length'] * 100.0) >= ref_cov_min
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            denom = summary[final_count_name(cat)]
            cov = sub[sub['ref_body_cov_flag']].groupby('CB')['_count'].sum() if not sub.empty else pd.Series(dtype=float)
            summary[f"{tag}_ref_coverage_prop"] = safe_prop(cov.reindex(summary.index, fill_value=0), denom).fillna(0)
        # Store the cutoff used so downstream plotting can incorporate it in titles
        summary['ref_cov_min_pct'] = ref_cov_min

        rts = cls_valid[cls_valid['RTS_stage'].astype(str) == 'TRUE'].groupby('CB')['_count'].sum()
        summary['RTS_prop_in_cell'] = safe_prop(rts.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        # Add per-category RTS across all structural categories
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['RTS_stage'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{tag}_RTS_prop"] = safe_prop(numer, denom).fillna(0)
        for abbr, cat in abbr_pairs:
            root = cat_to_root[cat]
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['RTS_stage'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_RTS_prop"] = safe_prop(numer, denom).fillna(0)

        noncanon = cls_valid[cls_valid['all_canonical'] == 'non_canonical'].groupby('CB')['_count'].sum()
        summary['Non_canonical_prop_in_cell'] = safe_prop(noncanon.reindex(summary.index, fill_value=0), summary['total_reads_no_monoexon']).fillna(0)
        # Add per-category Non-canonical across all structural categories (exclude mono-exon)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'non_canonical')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{tag}_noncanon_prop"] = safe_prop(numer, denom).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'non_canonical')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_noncanon_prop"] = safe_prop(numer, denom).fillna(0)

        intr = cls_valid[cls_valid['perc_A_downstream_TTS'] >= 60].groupby('CB')['_count'].sum()
        summary['Intrapriming_prop_in_cell'] = safe_prop(intr.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        # Add per-category Intra-priming across all structural categories
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['perc_A_downstream_TTS'] >= 60)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{tag}_intrapriming_prop"] = safe_prop(numer, denom).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['perc_A_downstream_TTS'] >= 60)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_intrapriming_prop"] = safe_prop(numer, denom).fillna(0)

        tss_sup = (cls_valid['diff_to_gene_TSS'].abs() <= 50)
        sup_cnt = cls_valid[tss_sup].groupby('CB')['_count'].sum()
        summary['TSSAnnotationSupport_prop'] = safe_prop(sup_cnt.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        # Add per-category TSS support across all structural categories
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['diff_to_gene_TSS'].abs() <= 50)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{tag}_TSSAnnotationSupport"] = safe_prop(numer, denom).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['diff_to_gene_TSS'].abs() <= 50)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_TSSAnnotationSupport"] = safe_prop(numer, denom).fillna(0)

        summary['Annotated_genes_prop_in_cell'] = safe_prop(summary['Annotated_genes'], summary['Genes_in_cell']).fillna(0)

        anno_models = cls_valid[(cls_valid['exons'] > 1) & (~cls_valid['associated_transcript'].fillna('').str.startswith('novel'))] \
            .groupby('CB')['jxn_string'].nunique()
        
        if args.mode != 'isoforms':
            summary['Annotated_juction_strings_prop_in_cell'] = safe_prop(anno_models.reindex(summary.index, fill_value=0), summary['UJCs_in_cell']).fillna(0)
        else:
            summary['Annotated_juction_strings_prop_in_cell'] = 0

        canon_over = cls_valid[cls_valid['all_canonical'] == 'canonical'].groupby('CB')['_count'].sum()
        summary['Canonical_prop_in_cell'] = safe_prop(canon_over.reindex(summary.index, fill_value=0), summary['total_reads_no_monoexon']).fillna(0)
        # Add per-category Canonical across all structural categories (exclude mono-exon)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB')['_count'].sum()
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'canonical')].groupby('CB')['_count'].sum()
            summary[f"{tag}_canon_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB')['_count'].sum()
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'canonical')].groupby('CB')['_count'].sum()
            summary[f"{abbr}_canon_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)

        if not args.skipORF:
            nmd = cls_valid[cls_valid['predicted_NMD'].astype(str) == 'TRUE'].groupby('CB')['_count'].sum()
            summary['NMD_prop_in_cell'] = safe_prop(nmd.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM','ISM','NIC','NNC'] else tag.lower()
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['predicted_NMD'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum()
                summary[f"{coding_tag}_NMD_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                sub = cls_valid[cls_valid['structural_category'] == cat]
                denom = sub.groupby('CB')['_count'].sum()
                cod = sub[sub['coding'] == 'coding'].groupby('CB')['_count'].sum()
                ncod = sub[sub['coding'] == 'non_coding'].groupby('CB')['_count'].sum()
                coding_tag = tag if tag in ['FSM','ISM','NIC','NNC'] else tag.lower()
                summary[f"{coding_tag}_coding_prop"] = safe_prop(cod.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
                summary[f"{coding_tag}_non_coding_prop"] = safe_prop(ncod.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['NMD_prop_in_cell'] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM','ISM','NIC','NNC'] else tag.lower()
                summary[f"{coding_tag}_coding_prop"] = 0
                summary[f"{coding_tag}_non_coding_prop"] = 100

        if getattr(args, 'CAGE_peak', None):
            cage = cls_valid[cls_valid['within_CAGE_peak'].astype(str) == 'TRUE'].groupby('CB')['_count'].sum()
            summary['CAGE_peak_support_prop'] = safe_prop(cage.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            # Add per-category CAGE support across all structural categories
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['within_CAGE_peak'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                summary[f"{tag}_CAGE_peak_support_prop"] = safe_prop(numer, denom).fillna(0)
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['within_CAGE_peak'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum()
                summary[f"{abbr}_CAGE_peak_support_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['CAGE_peak_support_prop'] = 0
            # Initialize per-category CAGE support to 0 when unavailable
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_CAGE_peak_support_prop"] = 0

        if getattr(args, 'polyA_motif_list', None):
            pa = cls_valid[cls_valid['polyA_motif_found'].astype(str) == 'TRUE'].groupby('CB')['_count'].sum()
            summary['PolyA_motif_support_prop'] = safe_prop(pa.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            # Add per-category PolyA support across all structural categories
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['polyA_motif_found'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                summary[f"{tag}_PolyA_motif_support_prop"] = safe_prop(numer, denom).fillna(0)
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['polyA_motif_found'].astype(str) == 'TRUE')].groupby('CB')['_count'].sum()
                summary[f"{abbr}_PolyA_motif_support_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['PolyA_motif_support_prop'] = 0
            # Initialize per-category PolyA support to 0 when unavailable
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_PolyA_motif_support_prop"] = 0

        # Short reads coverage support (if available)
        if 'min_cov' in cls_valid.columns:
             # Ensure min_cov is valid numeric
            cls_valid['min_cov'] = pd.to_numeric(cls_valid['min_cov'], errors='coerce').fillna(0)
             # Define support using configurable threshold
            sr_sup = cls_valid[cls_valid['min_cov'] >= args.min_cov].groupby('CB')['_count'].sum()
            summary['srjunctions_support_prop'] = safe_prop(sr_sup.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)

            # Add per-category SR support
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['min_cov'] >= args.min_cov)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                summary[f"{tag}_srjunctions_support_prop"] = safe_prop(numer, denom).fillna(0)
            
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['min_cov'] >= args.min_cov)].groupby('CB')['_count'].sum()
                summary[f"{abbr}_srjunctions_support_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
             summary['srjunctions_support_prop'] = 0
             for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_srjunctions_support_prop"] = 0

        # TSS Ratio Validation Support (configurable threshold)
        # Check if 'ratio_TSS' column exists (and has valid data)
        tss_threshold = args.ratio_TSS_threshold
        if 'ratio_TSS' in cls_valid.columns:
             # Ensure numeric, fill NA with 0 (unvalidated)
            cls_valid['ratio_TSS'] = pd.to_numeric(cls_valid['ratio_TSS'], errors='coerce').fillna(0)
            
            # Identify validated transcripts
            tss_valid = cls_valid[cls_valid['ratio_TSS'] >= tss_threshold].groupby('CB')['_count'].sum()
            summary['TSS_ratio_validated_prop'] = safe_prop(tss_valid.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)

            # Add per-category TSS Validation
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['ratio_TSS'] >= tss_threshold)].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
                summary[f"{tag}_TSS_ratio_validated_prop"] = safe_prop(numer, denom).fillna(0)
            
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB')['_count'].sum()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['ratio_TSS'] >= tss_threshold)].groupby('CB')['_count'].sum()
                summary[f"{abbr}_TSS_ratio_validated_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
             summary['TSS_ratio_validated_prop'] = 0
             for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_TSS_ratio_validated_prop"] = 0

        summary = summary.reset_index()
        if args.mode == 'isoforms':
            summary = summary.rename(columns={'total_reads': 'Transcripts_in_cell', 'total_reads_no_monoexon': 'total_transcripts_no_monoexon', 'MT_reads_count': 'MT_transcripts_count'})
        else:
            summary = summary.rename(columns={'total_reads': 'Reads_in_cell', 'total_UMI': 'UMIs_in_cell'})

        for c in summary.columns[1:]:
            summary[c] = pd.to_numeric(summary[c], errors='coerce').fillna(0)

        try:
            summary.to_csv(out_summary, sep='\t', index=False, compression='gzip')
            print(f"**** Cell summary written: {out_summary}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_summary}: {e}", file=sys.stderr)


