"""
sc_export.py — Export SQANTI-sc outputs as AnnData .h5ad files.

Produces a single .h5ad per sample containing:
  .X             — gene × cell count matrix (sparse, both modes)
  .layers["isoform"] — isoform × cell count matrix (isoforms mode only)
  .obs           — per-cell QC metrics from _SQANTI_cell_summary.txt.gz
  .obs["cluster"] — cluster labels from umap_results.csv (if available)
  .obsm["X_umap"] — UMAP coordinates (if available)
  .var           — gene/feature metadata
"""

import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Classification preparation
# ---------------------------------------------------------------------------

def _prepare_classification(class_file, mode):
    """Read and prepare the classification file for matrix construction.

    In isoforms mode the comma-separated CB and FL columns are exploded
    so that each row represents one (isoform, cell) pair with its
    count stored in ``_count``.
    """
    cls = pd.read_csv(class_file, sep="\t", dtype=str, low_memory=False)

    required = ["CB", "associated_gene"]
    for col in required:
        if col not in cls.columns:
            print(
                f"[ERROR] Column '{col}' not found in {class_file}.",
                file=sys.stderr,
            )
            return None

    if mode == "isoforms":
        if "isoform" not in cls.columns:
            print(
                "[ERROR] Column 'isoform' not found in classification.",
                file=sys.stderr,
            )
            return None

        cls["CB"] = cls["CB"].fillna("")
        if "FL" in cls.columns:
            cls["FL"] = cls["FL"].fillna("1")
            cls["CB"] = cls["CB"].str.split(",")
            cls["FL"] = cls["FL"].astype(str).str.split(",")
            cls = cls.explode(["CB", "FL"])
            cls["_count"] = pd.to_numeric(cls["FL"], errors="coerce").fillna(1)
        else:
            cls = cls.assign(CB=cls["CB"].str.split(",")).explode("CB")
            cls["_count"] = 1
    else:
        cls["_count"] = 1

    # Filter to valid barcodes
    cls = cls[(cls["CB"].notna()) & (cls["CB"] != "") & (cls["CB"] != "NA")]

    if cls.empty:
        print("[WARNING] No valid cell barcodes found.", file=sys.stdout)
        return None

    return cls


# ---------------------------------------------------------------------------
# Count matrix builders
# ---------------------------------------------------------------------------

def _build_gene_matrix(cls_df, mode):
    """Build a cell × gene count matrix.

    Returns
    -------
    sp.csr_matrix, np.ndarray, np.ndarray
        Sparse matrix (cells × genes), barcode array, gene array.
    """
    if mode == "isoforms":
        counts = (
            cls_df.groupby(["CB", "associated_gene"])["_count"]
            .sum()
            .reset_index()
        )
    else:
        counts = (
            cls_df.groupby(["CB", "associated_gene"])
            .size()
            .reset_index(name="_count")
        )

    barcodes = sorted(counts["CB"].unique())
    genes = sorted(counts["associated_gene"].unique())
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    gene_idx = {g: i for i, g in enumerate(genes)}

    row = counts["CB"].map(bc_idx).values
    col = counts["associated_gene"].map(gene_idx).values
    data = counts["_count"].values.astype(np.float32)

    mat = sp.csr_matrix(
        (data, (row, col)), shape=(len(barcodes), len(genes))
    )
    return mat, np.array(barcodes), np.array(genes)


def _build_isoform_matrix(cls_df, barcodes):
    """Build a cell × isoform count matrix (isoforms mode only).

    Uses the same barcode ordering as the gene matrix so the two
    matrices are aligned on the obs axis.

    Returns
    -------
    sp.csr_matrix, np.ndarray, np.ndarray
        Sparse matrix (cells × isoforms), isoform ID array,
        associated gene array (parallel to isoform IDs).
    """
    counts = (
        cls_df.groupby(["CB", "isoform"])["_count"]
        .sum()
        .reset_index()
    )

    # Only keep barcodes present in the gene matrix
    counts = counts[counts["CB"].isin(barcodes)]

    # Build isoform → gene lookup (majority vote)
    iso_gene = (
        cls_df.groupby("isoform")["associated_gene"]
        .agg(lambda s: s.value_counts().idxmax())
    )

    isoforms = sorted(counts["isoform"].unique())
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    iso_idx = {g: i for i, g in enumerate(isoforms)}

    row = counts["CB"].map(bc_idx).values
    col = counts["isoform"].map(iso_idx).values
    data = counts["_count"].values.astype(np.float32)

    mat = sp.csr_matrix(
        (data, (row, col)), shape=(len(barcodes), len(isoforms))
    )
    iso_arr = np.array(isoforms)
    gene_arr = np.array([iso_gene.get(i, "") for i in isoforms])
    return mat, iso_arr, gene_arr


# ---------------------------------------------------------------------------
# AnnData builder
# ---------------------------------------------------------------------------

def _build_anndata(gene_mat, barcodes, genes, iso_mat=None,
                   isoforms=None, iso_genes=None,
                   cell_summary_file=None, clustering_file=None):
    """Create an AnnData object with gene matrix in .X and optional layers.

    Parameters
    ----------
    gene_mat : sp.csr_matrix
        Cell × gene count matrix.
    barcodes : np.ndarray
        Cell barcode strings.
    genes : np.ndarray
        Gene name strings.
    iso_mat : sp.csr_matrix or None
        Cell × isoform count matrix (isoforms mode only).
    isoforms : np.ndarray or None
        Isoform ID strings.
    iso_genes : np.ndarray or None
        Associated gene per isoform.
    cell_summary_file : str or None
        Path to _SQANTI_cell_summary.txt.gz.
    clustering_file : str or None
        Path to umap_results.csv.
    """
    var_df = pd.DataFrame({"gene_name": genes}, index=genes)
    var_df.index.name = "gene_id"

    adata = ad.AnnData(X=gene_mat, var=var_df)
    adata.obs_names = list(barcodes)

    # Isoform layer
    if iso_mat is not None and isoforms is not None:
        iso_var = pd.DataFrame(
            {"isoform_id": isoforms, "associated_gene": iso_genes},
            index=isoforms,
        )
        # Store as a separate AnnData in .uns for now, since layers
        # require matching .var dimensions. We use .obsm instead.
        # Actually, standard practice: store isoform matrix in .obsm
        # is not ideal either. Best approach: store as a varm/layer
        # with its own feature axis, which AnnData doesn't natively
        # support across different var axes.
        #
        # Practical solution: store isoform counts as a dense/sparse
        # matrix in .obsm with feature names in .uns.
        adata.obsm["isoform_counts"] = iso_mat
        adata.uns["isoform_features"] = iso_var.to_dict()

    # Embed cell summary as .obs
    if cell_summary_file and os.path.isfile(cell_summary_file):
        summary = pd.read_csv(cell_summary_file, sep="\t",
                              compression="gzip")
        bc_col = summary.columns[0]
        summary = summary.set_index(bc_col)
        common = summary.index.intersection(adata.obs_names)
        if len(common) > 0:
            for col in summary.columns:
                adata.obs[col] = summary[col].reindex(adata.obs_names)

    # Embed clustering results
    if clustering_file and os.path.isfile(clustering_file):
        clust = pd.read_csv(clustering_file)
        if "Barcode" in clust.columns:
            clust = clust.set_index("Barcode")
            common = clust.index.intersection(adata.obs_names)
            if len(common) > 0:
                if "Cluster" in clust.columns:
                    adata.obs["cluster"] = (
                        clust["Cluster"]
                        .reindex(adata.obs_names)
                        .astype(str)
                    )
                if "UMAP_1" in clust.columns and "UMAP_2" in clust.columns:
                    umap_df = clust[["UMAP_1", "UMAP_2"]].reindex(
                        adata.obs_names
                    )
                    adata.obsm["X_umap"] = umap_df.values.astype(np.float32)

    return adata


# ---------------------------------------------------------------------------
# Main entry point (called from pipeline)
# ---------------------------------------------------------------------------

def export_h5ad(args, df):
    """Export .h5ad files for each sample in the design table.

    Called from qc_pipeline.py when --export_h5ad is set.
    """
    for _, row in df.iterrows():
        file_acc = row["file_acc"]
        sampleID = row["sampleID"]
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"

        if not os.path.isfile(class_file):
            print(
                f"[WARNING] Classification file not found for {file_acc}. "
                f"Skipping .h5ad export.",
                file=sys.stdout,
            )
            continue

        print(f"**** Exporting .h5ad for {sampleID} ({file_acc})...",
              file=sys.stdout)

        cls = _prepare_classification(class_file, args.mode)
        if cls is None:
            continue

        # Gene-level matrix (.X)
        gene_mat, barcodes, genes = _build_gene_matrix(cls, args.mode)

        # Isoform-level matrix (isoforms mode only → .obsm)
        iso_mat = iso_arr = iso_genes = None
        if args.mode == "isoforms":
            iso_mat, iso_arr, iso_genes = _build_isoform_matrix(
                cls, barcodes
            )

        # Locate optional metadata files
        cell_summary = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
        clustering_dir = os.path.join(
            os.path.dirname(outputPathPrefix), "clustering"
        )
        clustering_file = os.path.join(clustering_dir, "umap_results.csv")

        adata = _build_anndata(
            gene_mat, barcodes, genes,
            iso_mat=iso_mat, isoforms=iso_arr, iso_genes=iso_genes,
            cell_summary_file=cell_summary,
            clustering_file=clustering_file,
        )

        h5ad_path = f"{outputPathPrefix}_sqanti_sc.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"**** AnnData written: {h5ad_path}", file=sys.stdout)
