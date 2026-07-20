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

    # ---- RT-switching junction column naming (shared by both mode paths) ----
    # The four splice-junction types, and the per-cell columns that feed the
    # report's two RT-switching violins. Totals already exist as
    # `{Known,Novel}_{canonical,non_canonical}_junctions`; here we add the
    # RTS-flagged numerators (all-junction plot) plus the distinct-intron totals
    # and RTS numerators (unique-junction plot). Column presence is the report's
    # capability signal: emitted only when the junctions file carries the needed
    # inputs (RTS_junction for the all-junction plot; RTS_junction + genomic
    # coords for the unique-junction plot), so the report skips a plot exactly
    # when its columns are absent -- mirroring the old runtime guards.
    junc_type_order = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
    junc_total_name = {
        'known_canonical': 'Known_canonical_junctions',
        'known_non_canonical': 'Known_non_canonical_junctions',
        'novel_canonical': 'Novel_canonical_junctions',
        'novel_non_canonical': 'Novel_non_canonical_junctions',
    }
    rts_junc_name = {jc: 'RTS_' + junc_total_name[jc] for jc in junc_type_order}
    unique_junc_name = {jc: 'unique_' + junc_total_name[jc] for jc in junc_type_order}
    rts_unique_junc_name = {jc: 'RTS_unique_' + junc_total_name[jc] for jc in junc_type_order}

    def _rts_true(series):
        return series.astype(str).str.lower().isin(['true', 't', '1', 'yes'])

    warnings.simplefilter('ignore', PerformanceWarning)

    def _isoforms_summary(cls, junc):
        """Isoforms-mode per-cell metrics without exploding the classification frame.

        The comma-separated CB/FL columns ARE a sparse (transcript x cell) matrix.
        Instead of materialising one row per (transcript, cell) -- which duplicates
        every attribute column ~110M times and OOMs at >1TB -- we parse CB/FL once
        into three flat arrays (row=transcript idx, col=cell idx, data=FL weight) and
        express every per-cell metric as a masked column-sum over those arrays. All
        transcript-level conditions are evaluated on the small T-row frame, so peak
        memory is a few GB (the nnz arrays) rather than the exploded frame. Output is
        byte-identical to the previous explode()-based implementation.
        """
        cls = cls.reset_index(drop=True)
        T = len(cls)

        # ---- parse CB/FL lists into COO arrays (row=transcript, col=cell, data=FL) ----
        cb_ser = cls['CB'].fillna('') if 'CB' in cls.columns else pd.Series([''] * T)
        if 'FL' in cls.columns:
            fl_ser = cls['FL'].fillna('1').astype(str)
            _tmp = pd.DataFrame({'CB': cb_ser.str.split(','), 'FL': fl_ser.str.split(',')})
            _tmp = _tmp.explode(['CB', 'FL'])
            row = _tmp.index.to_numpy()
            cb_flat = _tmp['CB'].to_numpy()
            data = pd.to_numeric(pd.Series(_tmp['FL'].to_numpy()), errors='coerce').to_numpy()
            data = np.where(np.isnan(data), 1.0, data)
            del _tmp
        else:
            _split = cb_ser.str.split(',')
            _lens = _split.str.len().fillna(0).astype(int).to_numpy()
            row = np.repeat(np.arange(T), _lens)
            cb_flat = np.concatenate(_split.to_numpy()) if _lens.sum() > 0 else np.array([], dtype=object)
            data = np.ones(len(cb_flat))

        keep = pd.notna(cb_flat) & (cb_flat != '')
        row = row[keep]
        cb_flat = cb_flat[keep]
        data = data[keep].astype(np.float64)
        # sorted unique cells + inverse mapping == searchsorted, matches groupby('CB') order
        cells, col = np.unique(cb_flat.astype(str), return_inverse=True)
        col = col.astype(np.int64)
        Cn = len(cells)
        cells_index = pd.Index(cells, name='CB')
        del cb_flat, keep

        # ---- transcript-level attribute arrays / masks (length T) ----
        def _num(colname):
            if colname in cls.columns:
                return pd.to_numeric(cls[colname], errors='coerce').to_numpy()
            return np.full(T, np.nan)

        def _str(colname):
            if colname in cls.columns:
                return cls[colname].astype(str).to_numpy()
            return np.full(T, 'nan', dtype=object)

        exons = _num('exons')
        length = _num('length')
        ref_length = _num('ref_length')
        percA = _num('perc_A_downstream_TTS')
        difftss = _num('diff_to_gene_TSS')
        mincov = np.nan_to_num(_num('min_cov'), nan=0.0)
        ratiotss = np.nan_to_num(_num('ratio_TSS'), nan=0.0)
        sc = _str('structural_category')
        subcat = _str('subcategory')
        allcanon = _str('all_canonical')
        chrom = _str('chrom')
        rts = _str('RTS_stage')
        nmd = _str('predicted_NMD')
        cage = _str('within_CAGE_peak')
        polya = _str('polyA_motif_found')
        coding = _str('coding')
        gene_s = cls['associated_gene'].fillna('').astype(str) if 'associated_gene' in cls.columns else pd.Series([''] * T)
        novel_gene = gene_s.str.startswith('novel').to_numpy()
        anno_mask = ~novel_gene
        assoc_tx = cls['associated_transcript'].fillna('').astype(str) if 'associated_transcript' in cls.columns else pd.Series([''] * T)
        gene_code = pd.factorize(cls['associated_gene'])[0] if 'associated_gene' in cls.columns else np.full(T, -1)

        m_all = np.ones(T, dtype=bool)
        m_ne1 = exons != 1
        m_gt1 = exons > 1
        m_eq1 = exons == 1

        def m_cat(cat):
            return sc == cat

        # ---- sparse per-cell reducers ----
        def msum(mask=None):
            if mask is None:
                s = np.bincount(col, weights=data, minlength=Cn)
            else:
                sel = mask[row]
                s = np.bincount(col[sel], weights=data[sel], minlength=Cn)
            return pd.Series(s, index=cells_index)

        def wsum(weight_T):
            w = weight_T[row].astype(np.float64) * data
            s = np.bincount(col, weights=w, minlength=Cn)
            return pd.Series(s, index=cells_index)

        def distinct(codes_T, mask=None):
            sel = (codes_T[row] >= 0)
            if mask is not None:
                sel = sel & mask[row]
            if not sel.any():
                return pd.Series(np.zeros(Cn), index=cells_index)
            g = codes_T[row][sel].astype(np.int64)
            c = col[sel]
            key = g * Cn + c
            uk = np.unique(key)
            cc = (uk % Cn).astype(np.int64)
            return pd.Series(np.bincount(cc, minlength=Cn).astype(float), index=cells_index)

        def weighted_median():
            wint = data.astype(np.int64)
            sel = np.isfinite(length[row]) & (wint > 0)
            if not sel.any():
                return pd.Series(np.zeros(Cn), index=cells_index)
            c = col[sel]
            L = length[row][sel].astype(np.float64)
            w = wint[sel]
            order = np.lexsort((L, c))
            c = c[order]; L = L[order]; w = w[order]
            Gcum = np.cumsum(w)
            tot = np.bincount(c, weights=w, minlength=Cn).astype(np.int64)
            Wbefore = np.cumsum(tot) - tot
            present = tot > 0
            p1 = (tot - 1) // 2
            p2 = tot // 2
            i1 = np.clip(np.searchsorted(Gcum, Wbefore + p1 + 1, side='left'), 0, len(L) - 1)
            i2 = np.clip(np.searchsorted(Gcum, Wbefore + p2 + 1, side='left'), 0, len(L) - 1)
            med = (L[i1] + L[i2]) / 2.0
            out = np.zeros(Cn)
            out[present] = med[present]
            return pd.Series(out, index=cells_index)

        if Cn == 0:
            return pd.DataFrame(index=pd.Index([], name='CB'))

        # ---- assemble summary (same column names/order as explode implementation) ----
        summary = pd.DataFrame({
            'total_reads': msum(None),
            'total_reads_no_monoexon': msum(m_ne1),
            'Median_length_per_cell': weighted_median(),
        }).fillna(0)

        for cat in structural_categories:
            summary[final_count_name(cat)] = msum(m_cat(cat))
        for cat in structural_categories:
            tag = final_count_name(cat)
            summary[f"{tag}_prop"] = safe_prop(summary[tag], summary['total_reads'])

        summary['Genes_in_cell'] = distinct(gene_code)

        mito_chroms = ['MT', 'chrM', 'chrMT', 'mt', 'Mt']
        mt_mask = pd.Series(chrom).isin(mito_chroms).to_numpy()
        summary['MT_reads_count'] = msum(mt_mask)
        summary['MT_perc'] = safe_prop(summary['MT_reads_count'], summary['total_reads'])

        summary['Annotated_genes'] = distinct(gene_code, anno_mask)
        summary['Novel_genes'] = distinct(gene_code, novel_gene)
        # FL (read) counts from annotated vs novel genes -- numerators the report
        # used to re-derive by re-exploding CB/FL (Annotated/Novel reads % plot).
        summary['Annotated_genes_reads'] = msum(anno_mask)
        summary['Novel_genes_reads'] = msum(novel_gene)

        # ---- junctions: per-isoform type counts weighted by per-cell FL ----
        junc_types = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
        junc_rename = {
            'known_canonical': 'Known_canonical_junctions',
            'known_non_canonical': 'Known_non_canonical_junctions',
            'novel_canonical': 'Novel_canonical_junctions',
            'novel_non_canonical': 'Novel_non_canonical_junctions',
        }
        have_junc = (junc is not None and not junc.empty and 'junction_category' in junc.columns
                     and 'canonical' in junc.columns and 'isoform' in junc.columns and 'isoform' in cls.columns)
        if have_junc:
            jt = junc[['isoform', 'junction_category', 'canonical']].copy()
            jt['jtype'] = jt['junction_category'].astype(str) + '_' + jt['canonical'].astype(str)
            jt_wide = jt.groupby(['isoform', 'jtype']).size().unstack(fill_value=0)
            per_tx = jt_wide.reindex(cls['isoform'].to_numpy())
            for jc in junc_types:
                w_T = per_tx[jc].fillna(0).to_numpy() if jc in per_tx.columns else np.zeros(T)
                summary[junc_rename[jc]] = wsum(w_T)
            summary['total_junctions'] = summary[[junc_rename[j] for j in junc_types]].sum(axis=1)
            for jc in junc_types:
                summary[f"{junc_rename[jc]}_prop"] = safe_prop(summary[junc_rename[jc]], summary['total_junctions'])
            # Per-(cell, structural category) junction-type counts: the same FL-weighted
            # totals restricted to transcripts of each category. Feeds the report's
            # "Junctions by Structural Category" violins. Stored as counts; the report
            # derives the per-category percentages (NA when a cell has no junctions of
            # that category).
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                cmask = m_cat(cat).astype(np.float64)
                for jc in junc_types:
                    base = per_tx[jc].fillna(0).to_numpy() if jc in per_tx.columns else np.zeros(T)
                    summary[f"{tag}_{jc}_junctions"] = wsum(base * cmask)

            # ---- RT-switching numerators (all-junction violin) ----
            # RTS-flagged per-transcript junction-type counts, FL-weighted like the
            # totals above. Emitted only if the junctions file has RTS_junction; the
            # denominators are the {Known,Novel}_*_junctions totals. The report
            # derives % = 100 * RTS / total per SJ type per cell.
            if 'RTS_junction' in junc.columns:
                jt_rts = junc.loc[_rts_true(junc['RTS_junction']),
                                  ['isoform', 'junction_category', 'canonical']].copy()
                jt_rts['jtype'] = jt_rts['junction_category'].astype(str) + '_' + jt_rts['canonical'].astype(str)
                per_tx_rts = (jt_rts.groupby(['isoform', 'jtype']).size().unstack(fill_value=0)
                              .reindex(cls['isoform'].to_numpy()))
                for jc in junc_types:
                    w_T = per_tx_rts[jc].fillna(0).to_numpy() if jc in per_tx_rts.columns else np.zeros(T)
                    summary[rts_junc_name[jc]] = wsum(w_T)

            # ---- unique (distinct-intron) junction counts per cell ----
            # A junction is "present" in a cell iff ANY transcript carrying it is
            # expressed there. Collapsing the SAME genomic intron across transcripts
            # makes this a set-union, not an FL-weighted sum, so it needs the
            # junction x transcript incidence times the transcript x cell presence:
            # P = Jt @ Xp, binarized. Distinct-per-cell counts (and their RTS subset)
            # are then column sums of P over each SJ type's rows. Kept fully sparse so
            # the (junction, cell) product is never densified. Emitted only when the
            # junctions file has RTS_junction + genomic coordinates.
            _chrom_col = 'chrom' if 'chrom' in junc.columns else ('chr' if 'chr' in junc.columns else None)
            have_coords = ('RTS_junction' in junc.columns and _chrom_col is not None
                           and {'genomic_start_coord', 'genomic_end_coord'}.issubset(junc.columns))
            if have_coords:
                from scipy import sparse
                iso_to_idx = pd.Series(np.arange(T), index=cls['isoform'].astype(str).to_numpy())
                jc_df = junc[['isoform', 'junction_category', 'canonical', _chrom_col,
                              'genomic_start_coord', 'genomic_end_coord', 'RTS_junction']].copy()
                if 'strand' in junc.columns:
                    jc_df['strand'] = junc['strand'].astype(str)
                else:
                    jc_df['strand'] = ''
                t_idx = iso_to_idx.reindex(jc_df['isoform'].astype(str).to_numpy()).to_numpy()
                gs = pd.to_numeric(jc_df['genomic_start_coord'], errors='coerce').to_numpy()
                ge = pd.to_numeric(jc_df['genomic_end_coord'], errors='coerce').to_numpy()
                chrom_s = jc_df[_chrom_col].astype(str).to_numpy()
                keep = (~np.isnan(t_idx)) & (~np.isnan(gs)) & (~np.isnan(ge)) & (chrom_s != '') & (chrom_s != 'nan')
                if keep.any():
                    t_idx = t_idx[keep].astype(np.int64)
                    jkey = (chrom_s[keep] + ':' + gs[keep].astype(np.int64).astype(str) + ':'
                            + ge[keep].astype(np.int64).astype(str) + ':' + jc_df['strand'].to_numpy()[keep])
                    jkey_code, _uniq = pd.factorize(jkey)
                    nJ = len(_uniq)
                    jtype_arr = (jc_df['junction_category'].astype(str) + '_'
                                 + jc_df['canonical'].astype(str)).to_numpy()[keep]
                    rts_arr = _rts_true(jc_df['RTS_junction']).to_numpy()[keep]
                    # per distinct junction: its SJ type (consistent across rows) and
                    # whether ANY supporting row is RTS.
                    jtype_J = np.empty(nJ, dtype=object)
                    jtype_J[jkey_code] = jtype_arr
                    rts_J = np.zeros(nJ, dtype=bool)
                    np.logical_or.at(rts_J, jkey_code, rts_arr)
                    # Jt: distinct junction x transcript incidence (binary)
                    Jt = sparse.coo_matrix((np.ones(len(jkey_code)), (jkey_code, t_idx)),
                                           shape=(nJ, T)).tocsr()
                    Jt.data[:] = 1.0
                    # Xp: transcript x cell presence (binary) from the parsed COO
                    Xp = sparse.coo_matrix((np.ones(len(row)), (row, col)), shape=(T, Cn)).tocsr()
                    Xp.data[:] = 1.0
                    P = Jt.dot(Xp)          # (nJ x Cn): supporting transcripts per (junction, cell)
                    P.data[:] = 1.0         # binarize -> present / absent
                    for jc in junc_types:
                        rmask = (jtype_J == jc)
                        if rmask.any():
                            u_tot = np.asarray(P[rmask, :].sum(axis=0)).ravel()
                            u_rts = np.asarray(P[rmask & rts_J, :].sum(axis=0)).ravel()
                        else:
                            u_tot = np.zeros(Cn)
                            u_rts = np.zeros(Cn)
                        summary[unique_junc_name[jc]] = pd.Series(u_tot, index=cells_index)
                        summary[rts_unique_junc_name[jc]] = pd.Series(u_rts, index=cells_index)
        else:
            for cn in list(junc_rename.values()) + ['total_junctions'] + [f"{v}_prop" for v in junc_rename.values()]:
                summary[cn] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                for jc in junc_types:
                    summary[f"{tag}_{jc}_junctions"] = 0

        # ---- subcategory props ----
        sublevels = {
            'full-splice_match': ['alternative_3end', 'alternative_3end5end', 'alternative_5end', 'reference_match', 'mono-exon'],
            'incomplete-splice_match': ['3prime_fragment', 'internal_fragment', '5prime_fragment', 'intron_retention', 'mono-exon'],
            'novel_in_catalog': ['combination_of_known_junctions', 'combination_of_known_splicesites', 'intron_retention', 'mono-exon_by_intron_retention', 'mono-exon'],
            'novel_not_in_catalog': ['at_least_one_novel_splicesite', 'intron_retention'],
            'genic': ['mono-exon', 'multi-exon'], 'antisense': ['mono-exon', 'multi-exon'],
            'fusion': ['intron_retention', 'multi-exon'], 'intergenic': ['mono-exon', 'multi-exon'],
            'genic_intron': ['mono-exon', 'multi-exon'],
        }
        for cat in structural_categories:
            denom = summary[final_count_name(cat)]
            for lv in sublevels[cat]:
                numer = msum(m_cat(cat) & (subcat == lv))
                summary[f"{cat_to_tag[cat]}_{lv.replace('-', '_')}_prop"] = safe_prop(numer, denom).fillna(0)

        # ---- gene bins (annotated genes only): abundance + unique isoforms ----
        # Two per-cell distributions over annotated genes, both from a single pass over
        # the (gene, cell) groups of the sparse COO -- no CB/FL explosion:
        #   * abundance = sum(FL) per gene -> half-open bins (0,1]..(9,10],(10,inf),
        #     named *_gene_reads_bin_* (renamed *_transcripts on output). Half-open so
        #     float/EM abundances never fall through a gap.
        #   * isoforms  = # distinct models per gene (group size) -> exact bins 1..10,>10,
        #     named anno_gene_isoforms_bin_* (always integer counts).
        # Novel-gene bins are intentionally not computed (no plot uses them).
        ab_labels = [f"anno_gene_reads_bin_{i}" for i in range(1, 11)] + ["anno_gene_reads_bin_10plus"]
        iso_labels = [f"anno_gene_isoforms_bin_{i}" for i in range(1, 11)] + ["anno_gene_isoforms_bin_10plus"]
        sel = anno_mask[row] & (gene_code[row] >= 0)
        if not sel.any():
            for lab in ab_labels + iso_labels:
                summary[lab] = 0.0
        else:
            g = gene_code[row][sel].astype(np.int64)
            c = col[sel]
            w = data[sel]
            key = g * Cn + c
            o = np.argsort(key, kind='stable')
            key = key[o]; w = w[o]
            first = np.ones(len(key), dtype=bool)
            first[1:] = key[1:] != key[:-1]
            idx_start = np.nonzero(first)[0]
            sums = np.add.reduceat(w, idx_start)                 # sum(FL) per (gene, cell)
            grp_sizes = np.diff(np.append(idx_start, len(key)))  # distinct models per (gene, cell)
            cc = (key[idx_start] % Cn).astype(np.int64)          # cell index per group

            def bincount_col(binmask, lab):
                summary[lab] = pd.Series(np.bincount(cc[binmask], minlength=Cn).astype(float), index=cells_index)

            # abundance: half-open (i-1, i] for i=1..10, then >10
            for i in range(1, 11):
                bincount_col((sums > i - 1) & (sums <= i), ab_labels[i - 1])
            bincount_col(sums > 10, ab_labels[10])
            # unique isoforms: exact model counts 1..10, then >10
            for i in range(1, 11):
                bincount_col(grp_sizes == i, iso_labels[i - 1])
            bincount_col(grp_sizes > 10, iso_labels[10])

        # ---- length bins ----
        def lenbin_masks(base_mask):
            return {
                'two_fifty': base_mask & (length <= 250),
                'five_hund': base_mask & (length > 250) & (length <= 500),
                'short': base_mask & (length > 500) & (length <= 1000),
                'mid': base_mask & (length > 1000) & (length <= 2000),
                'long': base_mask & (length > 2000),
            }

        lo = lenbin_masks(m_all)
        for dst, src in [('Total_250b_length_prop', 'two_fifty'), ('Total_500b_length_prop', 'five_hund'),
                         ('Total_short_length_prop', 'short'), ('Total_mid_length_prop', 'mid'), ('Total_long_length_prop', 'long')]:
            summary[dst] = safe_prop(msum(lo[src]), summary['total_reads']).fillna(0)
        lm = lenbin_masks(m_eq1)
        for dst, src in [('Total_250b_length_mono_prop', 'two_fifty'), ('Total_500b_length_mono_prop', 'five_hund'),
                         ('Total_short_length_mono_prop', 'short'), ('Total_mid_length_mono_prop', 'mid'), ('Total_long_length_mono_prop', 'long')]:
            summary[dst] = safe_prop(msum(lm[src]), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = summary[final_count_name(cat)]
            lc = lenbin_masks(m_cat(cat))
            for dst, src in [('250b_length_prop', 'two_fifty'), ('500b_length_prop', 'five_hund'),
                             ('short_length_prop', 'short'), ('mid_length_prop', 'mid'), ('long_length_prop', 'long')]:
                summary[f"{tag}_{dst}"] = safe_prop(msum(lc[src]), denom).fillna(0)
            lcm = lenbin_masks(m_cat(cat) & m_eq1)
            for dst, src in [('250b_length_mono_prop', 'two_fifty'), ('500b_length_mono_prop', 'five_hund'),
                             ('short_length_mono_prop', 'short'), ('mid_length_mono_prop', 'mid'), ('long_length_mono_prop', 'long')]:
                summary[f"{tag}_{dst}"] = safe_prop(msum(lcm[src]), denom).fillna(0)

        # ---- reference body coverage (FSM/ISM only) ----
        ref_cov_min = float(getattr(args, 'ref_cov_min_pct', 45.0))
        with np.errstate(divide='ignore', invalid='ignore'):
            ref_flag = (length / ref_length * 100.0) >= ref_cov_min
        for cat in ['full-splice_match', 'incomplete-splice_match']:
            tag = cat_to_tag[cat]
            denom = summary[final_count_name(cat)]
            summary[f"{tag}_ref_coverage_prop"] = safe_prop(msum(m_cat(cat) & ref_flag), denom).fillna(0)
        summary['ref_cov_min_pct'] = ref_cov_min

        # ---- RTS ----
        rts_true = rts == 'TRUE'
        summary['RTS_prop_in_cell'] = safe_prop(msum(rts_true), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_RTS_prop"] = safe_prop(msum(m_cat(cat) & rts_true), msum(m_cat(cat))).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_RTS_prop"] = safe_prop(msum(m_cat(cat) & rts_true), msum(m_cat(cat))).fillna(0)

        # ---- non-canonical ----
        nc_mask = allcanon == 'non_canonical'
        summary['Non_canonical_prop_in_cell'] = safe_prop(msum(nc_mask), summary['total_reads_no_monoexon']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_noncanon_prop"] = safe_prop(msum(m_cat(cat) & m_gt1 & nc_mask), msum(m_cat(cat) & m_gt1)).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_noncanon_prop"] = safe_prop(msum(m_cat(cat) & m_gt1 & nc_mask), msum(m_cat(cat) & m_gt1)).fillna(0)

        # ---- intra-priming ----
        intr_mask = percA >= 60
        summary['Intrapriming_prop_in_cell'] = safe_prop(msum(intr_mask), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_intrapriming_prop"] = safe_prop(msum(m_cat(cat) & intr_mask), msum(m_cat(cat))).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_intrapriming_prop"] = safe_prop(msum(m_cat(cat) & intr_mask), msum(m_cat(cat))).fillna(0)

        # ---- TSS annotation support ----
        tss_mask = np.abs(difftss) <= 50
        summary['TSSAnnotationSupport_prop'] = safe_prop(msum(tss_mask), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_TSSAnnotationSupport"] = safe_prop(msum(m_cat(cat) & tss_mask), msum(m_cat(cat))).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_TSSAnnotationSupport"] = safe_prop(msum(m_cat(cat) & tss_mask), msum(m_cat(cat))).fillna(0)

        summary['Annotated_genes_prop_in_cell'] = safe_prop(summary['Annotated_genes'], summary['Genes_in_cell']).fillna(0)
        summary['Annotated_juction_strings_prop_in_cell'] = 0

        # ---- canonical ----
        can_mask = allcanon == 'canonical'
        summary['Canonical_prop_in_cell'] = safe_prop(msum(can_mask), summary['total_reads_no_monoexon']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_canon_prop"] = safe_prop(msum(m_cat(cat) & m_gt1 & can_mask), msum(m_cat(cat) & m_gt1)).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_canon_prop"] = safe_prop(msum(m_cat(cat) & m_gt1 & can_mask), msum(m_cat(cat) & m_gt1)).fillna(0)

        # ---- NMD / coding ----
        if args.include_ORF:
            nmd_true = nmd == 'TRUE'
            summary['NMD_prop_in_cell'] = safe_prop(msum(nmd_true), summary['total_reads']).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM', 'ISM', 'NIC', 'NNC'] else tag.lower()
                summary[f"{coding_tag}_NMD_prop"] = safe_prop(msum(m_cat(cat) & nmd_true), msum(m_cat(cat))).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM', 'ISM', 'NIC', 'NNC'] else tag.lower()
                denom = msum(m_cat(cat))
                summary[f"{coding_tag}_coding_prop"] = safe_prop(msum(m_cat(cat) & (coding == 'coding')), denom).fillna(0)
                summary[f"{coding_tag}_non_coding_prop"] = safe_prop(msum(m_cat(cat) & (coding == 'non_coding')), denom).fillna(0)
        else:
            summary['NMD_prop_in_cell'] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM', 'ISM', 'NIC', 'NNC'] else tag.lower()
                summary[f"{coding_tag}_coding_prop"] = 0
                summary[f"{coding_tag}_non_coding_prop"] = 100

        # ---- CAGE ----
        if getattr(args, 'CAGE_peak', None):
            cage_true = cage == 'TRUE'
            summary['CAGE_peak_support_prop'] = safe_prop(msum(cage_true), summary['total_reads']).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_CAGE_peak_support_prop"] = safe_prop(msum(m_cat(cat) & cage_true), msum(m_cat(cat))).fillna(0)
            for abbr, cat in abbr_pairs:
                summary[f"{abbr}_CAGE_peak_support_prop"] = safe_prop(msum(m_cat(cat) & cage_true), msum(m_cat(cat))).fillna(0)
        else:
            summary['CAGE_peak_support_prop'] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_CAGE_peak_support_prop"] = 0

        # ---- PolyA ----
        if getattr(args, 'polyA_motif_list', None):
            pa_true = polya == 'TRUE'
            summary['PolyA_motif_support_prop'] = safe_prop(msum(pa_true), summary['total_reads']).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_PolyA_motif_support_prop"] = safe_prop(msum(m_cat(cat) & pa_true), msum(m_cat(cat))).fillna(0)
            for abbr, cat in abbr_pairs:
                summary[f"{abbr}_PolyA_motif_support_prop"] = safe_prop(msum(m_cat(cat) & pa_true), msum(m_cat(cat))).fillna(0)
        else:
            summary['PolyA_motif_support_prop'] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                summary[f"{tag}_PolyA_motif_support_prop"] = 0

        # ---- short-read junction support ----
        sr_mask = mincov >= args.min_cov
        summary['srjunctions_support_prop'] = safe_prop(msum(sr_mask), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_srjunctions_support_prop"] = safe_prop(msum(m_cat(cat) & sr_mask), msum(m_cat(cat))).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_srjunctions_support_prop"] = safe_prop(msum(m_cat(cat) & sr_mask), msum(m_cat(cat))).fillna(0)

        # ---- TSS ratio validation ----
        tss_v = ratiotss >= args.ratio_TSS_threshold
        summary['TSS_ratio_validated_prop'] = safe_prop(msum(tss_v), summary['total_reads']).fillna(0)
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            summary[f"{tag}_TSS_ratio_validated_prop"] = safe_prop(msum(m_cat(cat) & tss_v), msum(m_cat(cat))).fillna(0)
        for abbr, cat in abbr_pairs:
            summary[f"{abbr}_TSS_ratio_validated_prop"] = safe_prop(msum(m_cat(cat) & tss_v), msum(m_cat(cat))).fillna(0)

        return summary

    for _, r in df.iterrows():
        file_acc = r['file_acc']
        sampleID = r['sampleID']
        prefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{prefix}_classification.txt"
        junc_file = f"{prefix}_junctions.txt"

        out_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"

        # Read ONLY the columns the metrics actually use. The classification file has
        # ~80 columns; in isoforms mode the per-cell explode (below) turns each
        # transcript's comma-separated CB list into one row per (transcript, cell),
        # so for collapse (~8.9M transcripts) the frame becomes ~110M rows. Carrying
        # 80 object/str columns at that height needs >1TB and OOMs; restricting to the
        # ~23 needed columns is the single biggest memory saving and changes no output.
        _usecols = {'CB','isoform','associated_gene','structural_category','exons','length','ref_length',
                    'all_canonical','subcategory','chrom','UMI','jxn_string','associated_transcript',
                    'RTS_stage','predicted_NMD','within_CAGE_peak','polyA_motif_found','perc_A_downstream_TTS',
                    'diff_to_gene_TSS','coding','min_cov','ratio_TSS','FL'}
        try:
            cls = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                              usecols=lambda c: c in _usecols)
        except Exception:
            continue
        # In isoforms mode the per-cell CB column from junctions is NOT used: junctions
        # are aggregated per-isoform by (junction_category, canonical) and CB comes from
        # the classification table. That CB column is the bulk of the per-cell-enriched
        # junctions file (~70GB for Isosceles, dense EM output), so excluding it from the
        # read is the difference between OOM and a few GB. Reads mode still needs it.
        # RTS_junction + genomic coordinates feed the two RT-switching-by-SJ-type
        # violins. They are small next to the CB column that dominates the isoforms
        # junctions file (~70GB for dense EM output), so adding them keeps the read
        # cheap while letting the metrics be computed here instead of re-exploding
        # junctions x cells in the report.
        _junc_usecols = {'isoform','readID','read_id','ID','read_name','read',
                         'junction_category','canonical','RTS_junction',
                         'genomic_start_coord','genomic_end_coord','chrom','chr','strand'}
        if args.mode != 'isoforms':
            _junc_usecols.add('CB')
        try:
            junc = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False,
                               usecols=lambda c: c in _junc_usecols)
        except Exception:
            junc = pd.DataFrame()

        for col in ['CB','isoform','associated_gene','structural_category','exons','length','ref_length',
                    'all_canonical','subcategory','chrom','UMI','jxn_string','associated_transcript',
                    'RTS_stage','predicted_NMD','within_CAGE_peak','polyA_motif_found','perc_A_downstream_TTS','diff_to_gene_TSS','coding','min_cov','ratio_TSS']:
            if col not in cls.columns:
                cls[col] = np.nan
        for c in ['exons','length','ref_length','perc_A_downstream_TTS','diff_to_gene_TSS','min_cov','ratio_TSS']:
            cls[c] = pd.to_numeric(cls[c], errors='coerce').astype('float32')

        # Downcast low-cardinality string columns to 'category' BEFORE the explode so
        # the ~110M-row exploded frame stores 1-byte codes instead of full Python str
        # objects. Only columns that are either pure values (compared with ==/isin) or
        # tiny-cardinality groupby keys are converted: structural_category (9) and
        # subcategory (~25) are keys but so small that groupby expansion is negligible.
        # High-cardinality keys (CB, isoform, associated_gene) stay object to avoid the
        # categorical-groupby cartesian-product blow-up.
        for c in ['chrom','all_canonical','coding','RTS_stage','predicted_NMD',
                  'within_CAGE_peak','polyA_motif_found','structural_category','subcategory']:
            if c in cls.columns:
                cls[c] = cls[c].astype('category')

        if args.mode == 'isoforms':
            # Sparse per-cell path: avoids exploding the ~110M-row (transcript x cell)
            # frame. Produces the same columns as the reads-mode branch below.
            summary = _isoforms_summary(cls, junc)
            summary = summary.reset_index()
            # In isoforms mode the per-cell counter is transcript-level (sum of FL),
            # so the gene read-count columns follow the tool's convention and are
            # named *_transcripts / *_gene_transcripts_bin_* (cf. Transcripts_in_cell).
            gene_tx = {c: c.replace('_reads', '_transcripts') for c in summary.columns
                       if 'genes_reads' in c or 'gene_reads_bin' in c}
            summary = summary.rename(columns={'total_reads': 'Transcripts_in_cell', 'total_reads_no_monoexon': 'total_transcripts_no_monoexon', 'MT_reads_count': 'MT_transcripts_count', **gene_tx})
            for c in summary.columns[1:]:
                summary[c] = pd.to_numeric(summary[c], errors='coerce').fillna(0)
            try:
                summary.to_csv(out_summary, sep='\t', index=False, compression='gzip')
                print(f"**** Cell summary written: {out_summary}", file=sys.stdout)
            except Exception as e:
                print(f"[ERROR] Failed writing {out_summary}: {e}", file=sys.stderr)
            continue

        if False and args.mode == 'isoforms':
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

        # explode (isoforms mode) leaves DUPLICATE index labels. Any later .loc on a
        # non-unique index does a cartesian match (e.g. the weighted-median step's
        # index.repeat(_count) tried to allocate a 312 GiB indexer -> OOM). Reset to a
        # unique RangeIndex so all downstream positional repeats / lookups stay linear.
        cls = cls.reset_index(drop=True)
        cls_valid = cls[(cls['CB'].notna()) & (cls['CB'] != '')].copy()

        total_reads = cls_valid.groupby('CB')['_count'].sum().rename('total_reads')
        if args.mode == 'isoforms':
            total_umi = pd.Series(dtype=float) # No UMI
        else:
            total_umi = cls_valid.groupby('CB')['UMI'].nunique().rename('total_UMI')
        
        reads_no_mono = cls_valid[cls_valid['exons'] != 1].groupby('CB')['_count'].sum().rename('total_reads_no_monoexon')
        
        if 'length' in cls_valid.columns:
            length_df = cls_valid[['CB', 'length', '_count']].dropna(subset=['length']).copy()
            length_df['_count'] = pd.to_numeric(length_df['_count'], errors='coerce').fillna(1).astype(int)
            length_df = length_df[length_df['_count'] > 0]
            if not length_df.empty:
                weighted_length = length_df.loc[length_df.index.repeat(length_df['_count'])]
                median_length = weighted_length.groupby('CB')['length'].median().rename('Median_length_per_cell')
            else:
                median_length = pd.Series(dtype=float, name='Median_length_per_cell')
        else:
            median_length = pd.Series(dtype=float, name='Median_length_per_cell')

        if args.mode == 'isoforms':
             summary = pd.DataFrame(total_reads).join(reads_no_mono, how='outer').join(median_length, how='outer').fillna(0)
        else:
             summary = pd.DataFrame(total_reads).join(total_umi, how='outer').join(reads_no_mono, how='outer').join(median_length, how='outer').fillna(0)

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

        mito_chroms = ['MT', 'chrM', 'chrMT', 'mt', 'Mt']
        mt = cls_valid[cls_valid['chrom'].isin(mito_chroms)].groupby('CB')['_count'].sum().rename('MT_reads_count')
        summary = summary.join(mt, how='left').fillna({'MT_reads_count': 0})
        summary['MT_perc'] = safe_prop(summary['MT_reads_count'], summary['total_reads'])

        anno = (~cls_valid['associated_gene'].fillna('').str.startswith('novel'))
        summary['Annotated_genes'] = cls_valid[anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        summary['Novel_genes'] = cls_valid[~anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        # Read counts from annotated vs novel genes (report reads these directly).
        summary['Annotated_genes_reads'] = cls_valid[anno].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)
        summary['Novel_genes_reads'] = cls_valid[~anno].groupby('CB')['_count'].sum().reindex(summary.index, fill_value=0)

        if not junc.empty:
            junc_types = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
            junc_rename = {
                'known_canonical': 'Known_canonical_junctions',
                'known_non_canonical': 'Known_non_canonical_junctions',
                'novel_canonical': 'Novel_canonical_junctions',
                'novel_non_canonical': 'Novel_non_canonical_junctions'
            }

            if args.mode == 'isoforms':
                # Each junction gets weighted by the FL of its parent isoform in each
                # cell. Doing that by merging the raw junction table to the per-(isoform,
                # cell) table explodes to junctions x cells (tens of billions of rows for
                # collapse -> OOM). Instead, collapse junctions to per-isoform type COUNTS
                # first, then many-to-one join those onto the per-cell rows and multiply
                # by FL -- mathematically identical, no cartesian product.
                iso_col = next((c for c in ['isoform', 'readID', 'read_id', 'ID', 'read_name', 'read']
                                if c in junc.columns and c in cls_valid.columns), None)
                if iso_col is not None and 'junction_category' in junc.columns and 'canonical' in junc.columns:
                    jt = junc[[iso_col, 'junction_category', 'canonical']].copy()
                    jt['junction_type'] = jt['junction_category'].astype(str) + '_' + jt['canonical'].astype(str)
                    jt_per_iso = jt.groupby([iso_col, 'junction_type']).size().unstack(fill_value=0)
                    jtype_cols = list(jt_per_iso.columns)
                    merged = cls_valid[[iso_col, 'CB', '_count']].merge(
                        jt_per_iso, left_on=iso_col, right_index=True, how='inner')
                    cnt = merged['_count'].astype('float32')
                    for jc in jtype_cols:
                        merged[jc] = merged[jc].astype('float32') * cnt
                    counts = merged.groupby('CB')[jtype_cols].sum()
                else:
                    counts = pd.DataFrame(index=summary.index)
            else:
                # Reads mode: each junction row has a single CB; count rows.
                if 'CB' not in junc.columns or (junc['CB'].fillna('') == '').all():
                    iso_to_cb = cls_valid[['isoform', 'CB']].dropna().drop_duplicates() if 'isoform' in cls_valid.columns else pd.DataFrame()
                    if not iso_to_cb.empty and 'isoform' in junc.columns:
                        junc = pd.merge(junc, iso_to_cb, on='isoform', how='left')
                jv = junc[(junc['CB'].notna()) & (junc['CB'] != '')].copy()
                if not jv.empty:
                    jv['junction_type'] = jv['junction_category'].astype(str) + '_' + jv['canonical'].astype(str)
                    counts = jv.groupby(['CB', 'junction_type']).size().unstack(fill_value=0)
                else:
                    counts = pd.DataFrame(index=summary.index)

            if not counts.empty:
                for tp in junc_types:
                    if tp not in counts.columns:
                        counts[tp] = 0
                counts['total_junctions'] = counts[junc_types].sum(axis=1)
                counts = counts.rename(columns=junc_rename)
                for src, dst in [(v, f"{v}_prop") for v in junc_rename.values()]:
                    counts[dst] = safe_prop(counts[src].reindex(counts.index, fill_value=0), counts['total_junctions'])
                summary = summary.join(counts, how='left').fillna(0)
            else:
                for col in list(junc_rename.values()) + [f"{v}_prop" for v in junc_rename.values()] + ['total_junctions']:
                    summary[col] = 0
        else:
            summary[['Known_canonical_junctions', 'Known_non_canonical_junctions',
                      'Novel_canonical_junctions', 'Novel_non_canonical_junctions', 'total_junctions',
                      'Known_canonical_junctions_prop', 'Known_non_canonical_junctions_prop',
                      'Novel_canonical_junctions_prop', 'Novel_non_canonical_junctions_prop']] = 0

        # ---- per-(cell, structural category) junction-type counts ----
        # Each junction inherits its transcript's structural_category (from the
        # classification table) and CB. Feeds the report's "Junctions by Structural
        # Category" violins. Stored as counts; the report derives the per-category
        # percentages (NA when a cell has no junctions of that category).
        junc_cat_types = ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            for jc in junc_cat_types:
                summary[f"{tag}_{jc}_junctions"] = 0
        if not junc.empty and {'junction_category', 'canonical'}.issubset(junc.columns):
            iso_col_r = next((c for c in ['isoform', 'readID', 'read_id', 'ID', 'read_name', 'read']
                              if c in junc.columns and c in cls_valid.columns), None)
            jcat = junc.copy()
            if iso_col_r is not None:
                if 'CB' not in jcat.columns or (jcat['CB'].fillna('') == '').all():
                    jcat = jcat.merge(cls_valid[[iso_col_r, 'CB']].dropna().drop_duplicates(iso_col_r),
                                      on=iso_col_r, how='left')
                jcat = jcat.merge(cls_valid[[iso_col_r, 'structural_category']].drop_duplicates(iso_col_r),
                                  on=iso_col_r, how='left')
            if 'CB' in jcat.columns and 'structural_category' in jcat.columns:
                jcat = jcat[(jcat['CB'].notna()) & (jcat['CB'] != '') & (jcat['structural_category'].notna())].copy()
                if not jcat.empty:
                    jcat['junction_type'] = jcat['junction_category'].astype(str) + '_' + jcat['canonical'].astype(str)
                    g = (jcat.groupby(['CB', 'structural_category', 'junction_type'], observed=True)
                             .size().rename('n').reset_index())
                    for cat in structural_categories:
                        tag = cat_to_tag[cat]
                        sub = g[g['structural_category'] == cat]
                        if sub.empty:
                            continue
                        wide = sub.pivot_table(index='CB', columns='junction_type', values='n', fill_value=0)
                        for jc in junc_cat_types:
                            if jc in wide.columns:
                                summary[f"{tag}_{jc}_junctions"] = wide[jc].reindex(summary.index, fill_value=0)

        # ---- RT-switching all-junction / unique-junction counts per cell ----
        # Reads mode: each junction row already carries its read's CB (merged from
        # the classification table when the junctions file lacks CB). The
        # all-junction numerators count RTS rows per (CB, SJ type); the unique-junction
        # metric collapses identical genomic introns per cell first (a junction is
        # RTS if ANY supporting row is). Columns emitted only when RTS_junction (and,
        # for the unique plot, genomic coordinates) are present, so the report skips a
        # plot exactly when its columns are absent.
        if (not junc.empty and {'junction_category', 'canonical', 'RTS_junction'}.issubset(junc.columns)):
            jrt = junc.copy()
            iso_col_j = next((c for c in ['isoform', 'readID', 'read_id', 'ID', 'read_name', 'read']
                              if c in jrt.columns and c in cls_valid.columns), None)
            if ('CB' not in jrt.columns or (jrt['CB'].fillna('') == '').all()) and iso_col_j is not None:
                jrt = jrt.merge(cls_valid[[iso_col_j, 'CB']].dropna().drop_duplicates(iso_col_j),
                                on=iso_col_j, how='left')
            if 'CB' in jrt.columns:
                jrt = jrt[(jrt['CB'].notna()) & (jrt['CB'] != '')].copy()
                if not jrt.empty:
                    jrt['jtype'] = jrt['junction_category'].astype(str) + '_' + jrt['canonical'].astype(str)
                    jrt['rts'] = _rts_true(jrt['RTS_junction'])
                    # all-junction RTS numerators (denominators = *_junctions totals)
                    rts_all = (jrt[jrt['rts']].groupby(['CB', 'jtype']).size()
                               .unstack(fill_value=0))
                    for jc in junc_type_order:
                        summary[rts_junc_name[jc]] = (rts_all[jc].reindex(summary.index, fill_value=0)
                                                      if jc in rts_all.columns else 0)
                    # unique-junction counts: distinct genomic introns per (CB, SJ type)
                    _chrom_col = 'chrom' if 'chrom' in jrt.columns else ('chr' if 'chr' in jrt.columns else None)
                    if _chrom_col is not None and {'genomic_start_coord', 'genomic_end_coord'}.issubset(jrt.columns):
                        strand_s = jrt['strand'].astype(str) if 'strand' in jrt.columns else ''
                        jrt['jkey'] = (jrt[_chrom_col].astype(str) + ':'
                                       + jrt['genomic_start_coord'].astype(str) + ':'
                                       + jrt['genomic_end_coord'].astype(str) + ':' + strand_s)
                        valid = jrt['genomic_start_coord'].notna() & jrt['genomic_end_coord'].notna()
                        ju = jrt[valid]
                        if not ju.empty:
                            dd = (ju.groupby(['CB', 'jkey'])
                                    .agg(jtype=('jtype', 'first'), rts=('rts', 'any')).reset_index())
                            uniq_tot = dd.groupby(['CB', 'jtype']).size().unstack(fill_value=0)
                            uniq_rts = dd[dd['rts']].groupby(['CB', 'jtype']).size().unstack(fill_value=0)
                            for jc in junc_type_order:
                                summary[unique_junc_name[jc]] = (uniq_tot[jc].reindex(summary.index, fill_value=0)
                                                                 if jc in uniq_tot.columns else 0)
                                summary[rts_unique_junc_name[jc]] = (uniq_rts[jc].reindex(summary.index, fill_value=0)
                                                                     if jc in uniq_rts.columns else 0)

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

        # Gene read-count bins (annotated genes only) as raw gene COUNTS. Bins are the
        # 11 half-open intervals (0,1]..(9,10],(10,inf); the report labels them exactly
        # (1..10,>10) since reads-mode counts are always whole. Novel-gene bins are
        # intentionally not computed (no plot uses them).
        anno_counts = cls_valid[anno].groupby(['CB','associated_gene'])['_count'].sum().rename('read_count').reset_index()
        ab_labels = [f"anno_gene_reads_bin_{i}" for i in range(1, 11)] + ["anno_gene_reads_bin_10plus"]
        def _upper_half_open(lo, hi):
            return lambda s: ((s > lo) & (s <= hi)).sum()
        if not anno_counts.empty:
            agg_kwargs = {f"anno_gene_reads_bin_{i}": ('read_count', _upper_half_open(i - 1, i)) for i in range(1, 11)}
            agg_kwargs["anno_gene_reads_bin_10plus"] = ('read_count', lambda s: (s > 10).sum())
            abins = anno_counts.groupby('CB').agg(**agg_kwargs)
            summary = summary.join(abins.reindex(summary.index, fill_value=0))
        else:
            for lab in ab_labels:
                summary[lab] = 0

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
        # Only FSM and ISM have a meaningful associated reference transcript and ref_length;
        # other categories (NIC, NNC, Genic, etc.) should not have ref_coverage reported.
        ref_cov_categories = ['full-splice_match', 'incomplete-splice_match']
        for cat in ref_cov_categories:
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

        if args.include_ORF:
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
            gene_tx = {c: c.replace('_reads', '_transcripts') for c in summary.columns
                       if 'genes_reads' in c or 'gene_reads_bin' in c}
            summary = summary.rename(columns={'total_reads': 'Transcripts_in_cell', 'total_reads_no_monoexon': 'total_transcripts_no_monoexon', 'MT_reads_count': 'MT_transcripts_count', **gene_tx})
        else:
            summary = summary.rename(columns={'total_reads': 'Reads_in_cell', 'total_UMI': 'UMIs_in_cell'})

        for c in summary.columns[1:]:
            summary[c] = pd.to_numeric(summary[c], errors='coerce').fillna(0)

        try:
            summary.to_csv(out_summary, sep='\t', index=False, compression='gzip')
            print(f"**** Cell summary written: {out_summary}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_summary}: {e}", file=sys.stderr)


