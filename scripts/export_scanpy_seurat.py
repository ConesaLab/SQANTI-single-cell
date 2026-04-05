#!/usr/bin/env python3
"""
export_scanpy_seurat.py — Standalone script to export SQANTI-sc outputs
as AnnData .h5ad files.

Produces a single .h5ad containing:
  .X             — gene × cell count matrix (sparse, both modes)
  .obsm["isoform_counts"] — isoform × cell count matrix (isoforms mode only)
  .uns["isoform_features"] — feature metadata for the isoform matrix (isoforms mode only)
  .obs           — per-cell QC metrics from _SQANTI_cell_summary.txt.gz
  .obs["cluster"] — cluster labels from umap_results.csv (if available)
  .obsm["X_umap"] — UMAP coordinates (if available)
  .var           — gene/feature metadata

Usage:
    python export_scanpy_seurat.py --mode {reads,isoforms} --classification <class_file.tsv> \\
                                   [--cell_summary <summary_file.txt.gz>] \\
                                   [--clustering <umap_results.csv>] \\
                                   [-o ./results] [-p prefix]
"""

import argparse
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
    """Read and prepare the classification file for matrix construction."""
    cls = pd.read_csv(class_file, sep="\t", dtype=str, low_memory=False)

    required = ["CB", "associated_gene"]
    for col in required:
        if col not in cls.columns:
            print(f"[ERROR] Column '{col}' not found in {class_file}.", file=sys.stderr)
            sys.exit(1)

    if mode == "isoforms":
        if "isoform" not in cls.columns:
            print("[ERROR] Column 'isoform' not found in classification.", file=sys.stderr)
            sys.exit(1)

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
        print("[ERROR] No valid cell barcodes found.", file=sys.stderr)
        sys.exit(1)

    return cls


# ---------------------------------------------------------------------------
# Count matrix builders
# ---------------------------------------------------------------------------

def _build_gene_matrix(cls_df, mode):
    """Build a cell × gene count matrix."""
    if mode == "isoforms":
        counts = cls_df.groupby(["CB", "associated_gene"])["_count"].sum().reset_index()
    else:
        counts = cls_df.groupby(["CB", "associated_gene"]).size().reset_index(name="_count")

    barcodes = sorted(counts["CB"].unique())
    genes = sorted(counts["associated_gene"].unique())
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    gene_idx = {g: i for i, g in enumerate(genes)}

    row = counts["CB"].map(bc_idx).values
    col = counts["associated_gene"].map(gene_idx).values
    data = counts["_count"].values.astype(np.float32)

    mat = sp.csr_matrix((data, (row, col)), shape=(len(barcodes), len(genes)))
    return mat, np.array(barcodes), np.array(genes)


def _build_isoform_matrix(cls_df, barcodes):
    """Build a cell × isoform count matrix (isoforms mode only)."""
    counts = cls_df.groupby(["CB", "isoform"])["_count"].sum().reset_index()
    counts = counts[counts["CB"].isin(barcodes)]

    iso_gene = cls_df.groupby("isoform")["associated_gene"].agg(lambda s: s.value_counts().idxmax())

    isoforms = sorted(counts["isoform"].unique())
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    iso_idx = {g: i for i, g in enumerate(isoforms)}

    row = counts["CB"].map(bc_idx).values
    col = counts["isoform"].map(iso_idx).values
    data = counts["_count"].values.astype(np.float32)

    mat = sp.csr_matrix((data, (row, col)), shape=(len(barcodes), len(isoforms)))
    iso_arr = np.array(isoforms)
    gene_arr = np.array([iso_gene.get(i, "") for i in isoforms])
    return mat, iso_arr, gene_arr


# ---------------------------------------------------------------------------
# AnnData builder
# ---------------------------------------------------------------------------

def _build_anndata(gene_mat, barcodes, genes, iso_mat=None, isoforms=None, iso_genes=None,
                   cell_summary_file=None, clustering_file=None):
    """Create an AnnData object."""
    var_df = pd.DataFrame({"gene_name": genes}, index=genes)
    var_df.index.name = "gene_id"

    adata = ad.AnnData(X=gene_mat, var=var_df)
    adata.obs_names = list(barcodes)

    # Isoform data in .obsm
    if iso_mat is not None and isoforms is not None:
        iso_var = pd.DataFrame({"isoform_id": isoforms, "associated_gene": iso_genes}, index=isoforms)
        adata.obsm["isoform_counts"] = iso_mat
        adata.uns["isoform_features"] = iso_var.to_dict()

    # Embed cell summary
    if cell_summary_file and os.path.isfile(cell_summary_file):
        summary = pd.read_csv(cell_summary_file, sep="\t", compression="gzip")
        bc_col = summary.columns[0]
        summary = summary.set_index(bc_col)
        summary_aligned = summary.reindex(adata.obs_names)
        if summary_aligned.notna().any().any():
            adata.obs = pd.concat([adata.obs, summary_aligned], axis=1)

    # Embed clustering
    if clustering_file and os.path.isfile(clustering_file):
        clust = pd.read_csv(clustering_file)
        if "Barcode" in clust.columns:
            clust = clust.set_index("Barcode")
            clust_aligned = clust.reindex(adata.obs_names)
            if clust_aligned.notna().any().any():
                clust_cols = []
                if "Cluster" in clust.columns:
                    clust_cols.append(clust_aligned[["Cluster"]].rename(columns={"Cluster": "cluster"}).astype(str))
                if "UMAP_1" in clust.columns and "UMAP_2" in clust.columns:
                    adata.obsm["X_umap"] = clust_aligned[["UMAP_1", "UMAP_2"]].values.astype(np.float32)
                if clust_cols:
                    adata.obs = pd.concat([adata.obs] + clust_cols, axis=1)

    return adata


# ---------------------------------------------------------------------------
# Main export logic
# ---------------------------------------------------------------------------

def export(args):
    """Run the standalone export."""
    cls = _prepare_classification(args.classification, args.mode)

    print(f"**** Exporting .h5ad ({args.mode} mode)...", file=sys.stdout)
    os.makedirs(args.output_dir, exist_ok=True)

    gene_mat, barcodes, genes = _build_gene_matrix(cls, args.mode)

    iso_mat = iso_arr = iso_genes = None
    if args.mode == "isoforms":
        iso_mat, iso_arr, iso_genes = _build_isoform_matrix(cls, barcodes)

    adata = _build_anndata(
        gene_mat, barcodes, genes,
        iso_mat=iso_mat, isoforms=iso_arr, iso_genes=iso_genes,
        cell_summary_file=args.cell_summary,
        clustering_file=args.clustering,
    )

    h5ad_path = os.path.join(args.output_dir, f"{args.prefix}.h5ad")
    adata.write_h5ad(h5ad_path)
    print(f"**** AnnData written: {h5ad_path}", file=sys.stdout)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Export SQANTI-sc outputs to a Scanpy / Seurat compatible .h5ad file.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    req = p.add_argument_group("Required arguments")
    req.add_argument(
        "--mode", required=True, choices=["reads", "isoforms"],
        help="Pipeline mode used to generate the classification file.",
    )
    req.add_argument(
        "--classification", required=True,
        help="SQANTI-sc classification file (TSV).",
    )

    opt = p.add_argument_group("Optional inputs")
    opt.add_argument(
        "--cell_summary",
        help="SQANTI-sc cell summary file (_SQANTI_cell_summary.txt.gz).\n"
             "If provided, QC metrics are embedded in .obs.",
    )
    opt.add_argument(
        "--clustering",
        help="SQANTI-sc clustering results (umap_results.csv).\n"
             "If provided, UMAP coords and clusters are embedded in .obsm/.obs.",
    )

    out = p.add_argument_group("Output options")
    out.add_argument(
        "-o", "--output_dir", default=".",
        help="Output directory. Default: current directory.",
    )
    out.add_argument(
        "-p", "--prefix", default="sqanti_sc_export",
        help="Prefix for output .h5ad file. Default: sqanti_sc_export.",
    )

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not os.path.isfile(args.classification):
        print(f"[ERROR] Classification file not found: {args.classification}", file=sys.stderr)
        sys.exit(1)

    export(args)


if __name__ == "__main__":
    main()
