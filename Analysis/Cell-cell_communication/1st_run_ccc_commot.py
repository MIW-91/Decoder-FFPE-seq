
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import sys
import gc
import ot
import pickle
import anndata
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib as mpl
from matplotlib import rcParams

sample = sys.argv[1]


databasename = "CellChat"
species = "human"
os.chdir('commot_filepath')

file_path="my_filepath"+sample
# file_path=[path for path in paths if os.path.exists(path)]

adata=sc.read_visium(file_path)
adata.var_names_make_unique()
newid=sample+"_"+adata.obs_names
adata.obs_names = newid

###niche orig.ident
niche=pd.read_csv('my_filepath'+sample+'.spatialdomain.csv',index_col=0)
cell_ids = niche.index.tolist()
adata_use=adata[adata.obs_names.isin(cell_ids),:].copy()
niche=niche.loc[adata_use.obs.index]
adata_use.obs= adata_use.obs.join(niche, how='left')

sc.pp.normalize_total(adata_use, inplace=True)
sc.pp.log1p(adata_use)

dia=8 #spot diameter
conversion=dia/adata_use.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
print(conversion)

df_cellchat = ct.pp.ligand_receptor_database(species=species, signaling_type='Secreted Signaling', database=databasename)
df_cellchat = pd.concat([df_cellchat,ct.pp.ligand_receptor_database(species=species, signaling_type='Cell-Cell Contact', database=databasename)])
df_cellchat = pd.concat([df_cellchat,ct.pp.ligand_receptor_database(species=species, signaling_type='ECM-Receptor', database=databasename)])
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_use, min_cell_pct=0.01)
print(df_cellchat_filtered.shape)


# Spatial communication inference and Construct CCC networks
ct.tl.spatial_communication(adata_use,
    database_name=databasename, df_ligrec=df_cellchat_filtered, dis_thr=250/conversion, heteromeric=True, pathway_sum=True)

adata_use.write("".join([sample,"_commot.h5ad"]))
print("the summary of %s is ok!" % sample)
