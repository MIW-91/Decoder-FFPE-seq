import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import commot as ct
import seaborn as sns
import sys

sample = sys.argv[1]
databasename = "CellChat"
species = "human"
os.chdir('commot_filepath')


def assign_celltype_by_deconv(deconv_mat, top_n=2, avoid_types=['Epithelial']):
    celltype_assignments = []
    for _, row in deconv_mat.iterrows():
        sorted_types = row.sort_values(ascending=False).index.tolist()

        selected = sorted_types[:top_n]
        if any(typ in avoid_types for typ in selected) and len(sorted_types) > top_n:
            selected = sorted_types[:top_n+1]
        
        celltype_assignments.append(selected)
    return celltype_assignments

def create_pseudo_sc(adata, spot_celltypes):

    rep_indices = []
    celltype_labels = []
    for spot_idx, ct_list in enumerate(spot_celltypes):
        rep_indices.extend([spot_idx] * len(ct_list))
        celltype_labels.extend(ct_list)

    pseudo_adata = adata[rep_indices].copy()
    pseudo_adata.obs['celltype'] = celltype_labels

    return pseudo_adata

def replicate_communication_matrices(orig_adata, pseudo_adata, rep_indices):
    """
    expand the SRT to the pseudo single-cell data
    
    parameter:
        orig_adata : the original AnnData
        pseudo_adata : the pseudo single-cell AnnData
        rep_indices : the repetition of the original spot index 
    """
    for key in orig_adata.obsp.keys():
        if key.startswith('commot-'):
            orig_mat = orig_adata.obsp[key]
            pseudo_mat = orig_mat[rep_indices, :][:, rep_indices]
            pseudo_adata.obsp[key] = pseudo_mat

    if 'commot' in orig_adata.uns:
        pseudo_adata.uns['commot'] = orig_adata.uns['commot'].copy()

file_path='commot_filepath'
adata=sc.read_h5ad(file_path+sample+'_commot.h5ad')
adata.obs['region'] = np.where(adata.obs['region4'].isin(['T-M', 'M-P']), 'M', adata.obs['region4'])
adata = adata[adata.obs['region'] == 'M'].copy()

decon=pd.read_csv('deconvolution_path'+sample+'/'+sample+'_CARD.Proportion.csv',index_col=0)
newid=sample+"_"+decon.index
decon.index = newid
cell_ids = decon.index.tolist()

adata_use=adata[adata.obs_names.isin(cell_ids),:].copy()
decon=decon.loc[adata_use.obs.index]
decon.head(5)
spot_celltypes = assign_celltype_by_deconv(decon, top_n=2, avoid_types=['NormalEpith'])

pseudo_adata = create_pseudo_sc(adata_use,spot_celltypes=spot_celltypes)

rep_indices = []
for i, ct_list in enumerate(spot_celltypes):
    rep_indices.extend([i] * len(ct_list))

replicate_communication_matrices(adata_use, pseudo_adata, rep_indices)

def transform_to_tuples_per_group(df, group_col, col1, col2):
    grouped = df.groupby(group_col)
    result = {name: list(zip(group[col1], group[col2])) for name, group in grouped}
    return result

df = pseudo_adata.uns['commot-CellChat-info']['df_ligrec']
result = transform_to_tuples_per_group(df, 'pathway', 'ligand', 'receptor')
print(result)

a = 0 
for pathway,lrpairtplist in result.items():
    a += 1
    print("%d End of run command cluster_communication(), pathway: %s" % (a, pathway))
    for lrpairtp in lrpairtplist:
        print(pathway,lrpairtp)
        ct.tl.cluster_communication(pseudo_adata, database_name=databasename, lr_pair=lrpairtp, n_permutations=100,clustering='celltype')



lrpaircp=[]
for pathway,lrpairtplist in result.items():
    for lrpairtp in lrpairtplist:
        lrpairtpstr = '-'.join(str(x) for x in lrpairtp)
        lrcmatrix = pseudo_adata.uns['commot_cluster-celltype-%s-%s' % (databasename, lrpairtpstr)]['communication_matrix']
        lrcpvalue = pseudo_adata.uns['commot_cluster-celltype-%s-%s' % (databasename, lrpairtpstr)]['communication_pvalue']
        lrindices = [(row, col, lrcmatrix.loc[row, col], lrcpvalue.loc[row, col], lrpairtpstr, pathway) for row in lrcpvalue.index for col in lrcpvalue.columns if lrcpvalue.loc[row, col] < 0.05]
        lrpaircp.append(lrindices)

flattened_list = [item for sublist in lrpaircp for item in sublist]
df_lrpair = pd.DataFrame(flattened_list, columns=['celltype_source', 'celltype_target', 'Communication_score', 'p-values', 'LR-pair', 'pathway'])
df_lrpair.to_csv("./" + sample + "_lrpair_commot.csv", index=False)

print("the summary of %s is ok!" % sample)

adata_use.write("".join([sample,"_commot.h5ad"]))

print("the cluster_communication() analysis of %s is ok!" % sample)

