import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc

def prepare_anndata(args, row):
    """
    Load classification file and create an AnnData object.
    Aggregates counts by cell barcode (CB) and associated_gene.
    """
    file_acc = row['file_acc']
    sampleID = row['sampleID']
    outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
    class_file = f"{outputPathPrefix}_classification.txt"

    if not os.path.isfile(class_file):
        print(f"[ERROR] Classification file not found: {class_file}", file=sys.stderr)
        return None

    print(f"**** Loading data for clustering from {class_file}...", file=sys.stdout)
    
    try:
        # Read necessary columns
        cols = ['associated_gene', 'CB']
        if args.mode == 'isoforms':
            cols.append('FL')
        
        df = pd.read_csv(class_file, sep='\t', usecols=lambda c: c in cols, dtype=str)
        
        # Filter valid CBs
        df = df[(df['CB'].notna()) & (df['CB'] != '') & (df['CB'] != 'NA')]
        
        if df.empty:
            print(f"[WARNING] No valid cell barcodes found in {class_file}.", file=sys.stdout)
            return None

        # Expand data if needed (isoforms mode)
        if args.mode == 'isoforms':
            # CB and FL are comma-separated
            df['CB'] = df['CB'].str.split(',')
            if 'FL' in df.columns:
                df['FL'] = df['FL'].fillna('1').astype(str).str.split(',')
                df = df.explode(['CB', 'FL'])
                df['count'] = pd.to_numeric(df['FL'], errors='coerce').fillna(1)
            else:
                df = df.explode('CB')
                df['count'] = 1
        else:
            # Reads mode: each row is a read
            df['count'] = 1

        # Aggregate by Gene and Cell
        # Filter out novel genes if desired? User didn't specify, but usually we want all genes.
        # However, associated_gene might be "novelGene_X".
        
        # Create count matrix: Cell x Gene
        count_matrix = df.groupby(['CB', 'associated_gene'])['count'].sum().unstack(fill_value=0)
        
        # Create AnnData
        adata = sc.AnnData(count_matrix)
        adata.var_names_make_unique()
        
        return adata

    except Exception as e:
        print(f"[ERROR] Failed to prepare AnnData for {file_acc}: {e}", file=sys.stderr)
        return None

def run_clustering_analysis(args, row):
    """
    Run Scanpy clustering pipeline and save results.
    """
    file_acc = row['file_acc']
    sampleID = row['sampleID']
    outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
    clustering_dir = os.path.join(os.path.dirname(outputPathPrefix), "clustering")
    
    if not os.path.exists(clustering_dir):
        os.makedirs(clustering_dir)

    adata = prepare_anndata(args, row)
    
    if adata is None:
        return

    print(f"**** Running clustering analysis for {file_acc}...", file=sys.stdout)
    
    try:
        # NO filtering of cells/genes as requested
        
        # Normalization
        if args.normalization == 'log1p':
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        elif args.normalization == 'sqrt':
            adata.X = np.sqrt(adata.X)
        elif args.normalization == 'pearson':
            sc.experimental.pp.normalize_pearson_residuals(adata)

        # Highly Variable Genes
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes)
        
        # Scale
        sc.pp.scale(adata, max_value=10)
        
        # PCA
        sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_pc)
        
        # Neighbors
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pc)
        
        # UMAP
        sc.tl.umap(adata)
        
        # Clustering
        if args.clustering_method == 'leiden':
            sc.tl.leiden(adata, resolution=args.resolution)
            cluster_col = 'leiden'
        elif args.clustering_method == 'louvain':
            sc.tl.louvain(adata, resolution=args.resolution)
            cluster_col = 'louvain'
        elif args.clustering_method == 'kmeans':
            from sklearn.cluster import KMeans
            # K-means on PCA components
            kmeans = KMeans(n_clusters=args.n_clusters, random_state=0).fit(adata.obsm['X_pca'])
            adata.obs['kmeans'] = kmeans.labels_.astype(str)
            cluster_col = 'kmeans'
        
        # Save results
        # We need to export UMAP coordinates and Cluster IDs for R plotting
        
        umap_coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'], index=adata.obs_names)
        clusters = adata.obs[[cluster_col]].rename(columns={cluster_col: 'Cluster'})
        
        results = pd.concat([clusters, umap_coords], axis=1)
        results.index.name = 'Barcode'
        results = results.reset_index()
        
        out_file = os.path.join(clustering_dir, "umap_results.csv")
        results.to_csv(out_file, index=False)
        print(f"**** Clustering results saved to {out_file}", file=sys.stdout)
        
    except Exception as e:
        print(f"[ERROR] Clustering analysis failed for {file_acc}: {e}", file=sys.stderr)
