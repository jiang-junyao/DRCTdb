import matplotlib.pyplot as plt
import muon
import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import glob
import os
import re
import sys 




scvi.settings.seed = 0
sc.set_figure_params(dpi_save=300,figsize=(4, 4), frameon=False)
#os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:1024"


##load data
def usage():
    print('Usage: python script.py [input_h5ad]')

def run_scbasset(h5ad,outname,genome = 'hg38'):
    adata = sc.read_h5ad(h5ad)
    # compute the threshold: 5% of the cells
    min_cells = int(adata.shape[0] * 0.05)
    # in-place filtering of regions
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if "name" not in adata.var:
        adata.var["name"]  = adata.var_names
    split_interval = adata.var["name"].str.split("-", expand=True)
    adata.var["chr"] = split_interval[0]
    adata.var["start"] = split_interval[1].astype(int)
    adata.var["end"] = split_interval[2].astype(int)
    # Filter out non-chromosomal regions
    mask = adata.var["chr"].str.startswith("chr")
    adata = adata[:, mask].copy()
    if genome == 'hg38':
        scvi.data.add_dna_sequence(
            adata,
            genome_name="GRCh38",
            genome_dir="data",
            chr_var_key="chr",
            start_var_key="start",
            end_var_key="end",
            install_genome = False
        )
        print('Adding GRCh38 Genome')
        adata.varm["dna_sequence"]
    elif genome == 'hg19':
        scvi.data.add_dna_sequence(
            adata,
            genome_name="GRCh37",
            genome_dir="data",
            chr_var_key="chr",
            start_var_key="start",
            end_var_key="end",
            install_genome = False
        )
        print('Adding GRCh37 Genome')
        adata.varm["dna_sequence"]
    else:
        raise ValueError('Incorrect genome version')  
    ##start train
    bdata = adata.transpose()
    bdata.layers["binary"] = (bdata.X.copy() > 0).astype(float)
    scvi.external.SCBASSET.setup_anndata(bdata, layer="binary", dna_code_key="dna_code")
    bas = scvi.external.SCBASSET(bdata)
    bas.train(precision=16,batch_size=256)
    latent = bas.get_latent_representation()
    adata.obsm["X_scbasset"] = latent

    print('Successful learned latent representation form scbassset')
    # compute the k-nearest-neighbor graph that is used in both clustering and umap algorithms
    # sc.pp.neighbors(adata, use_rep="X_scbasset")
    # compute the umap
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, key_added="leiden_scbasset")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  

    ## output differential open chromtion regions
    sc.tl.rank_genes_groups(adata, 'cell_type', method='wilcoxon')
    DERs = sc.get.rank_genes_groups_df(adata, group = None, log2fc_min  = 1)    
    DERs.to_csv(f'{outname}_DERs.txt',sep = '\t')
    print('Write DERs')
    ##visualization
    # with plt.rc_context():  
    #     sc.pl.umap(adata, color="leiden_scbasset",title = '',show = False)
    #     plt.savefig(f"{outname}_umap_plot.svg", dpi = 300,bbox_inches="tight")


    motif_file = glob.glob("./data/motifs/shuffled_peaks_motifs/*.fasta")
    tf_name = [os.path.splitext(os.path.basename(i))[0] for i in motif_file]
    ##Calculate tf_activity
    def get_tf_activity(tf):
        tf_activ = bas.get_tf_activity(
                tf=tf,
                motif_dir="data/motifs",
            )
        return(tf_activ)
    tf_activity = []

    for tf in tqdm(tf_name):
        result = get_tf_activity(tf)
        tf_activity.append(result)

    tf_activity_df = pd.DataFrame({tf_name[i]: tf_activity[i] for i in range(len(tf_activity))}, index = adata.obs_names)
    tf_activity_df.to_csv(f'{outname}_tf_activity.txt',sep = '\t')
    print('Write TF activity')
    del adata.varm
    adata.write(f'{outname}_scbasset.h5ad', compression='gzip')
    

# def test(input,outname):
#     adata = sc.read_h5ad(input)
#     with plt.rc_context():  
#         sc.pl.umap(adata, color="leiden_scbasset",title = '',show = False)
#         plt.savefig(f"{outname}_umap_plot.svg", dpi = 300,bbox_inches="tight")  
#     print(f"{outname}_umap_plot.svg")
#     print('aa')

if __name__ == '__main__':
    try:      
        output_prefix =  os.path.splitext(sys.argv[1])[0]
        out_dir = re.search('\w+',output_prefix).group()

        if not os.path.lexists(out_dir):
            os.makedirs(out_dir)         
        run_scbasset(sys.argv[1],output_prefix,sys.argv[2])

    except:
        usage()

