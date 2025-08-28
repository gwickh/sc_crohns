import scanpy as sc
import os

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/"

gca_in = os.path.join(PATH, "gca_ref/Full_obj_raw_counts_nosoupx_v2.h5ad")
gca_out1 = os.path.join(PATH, "gca_ref/obj_healthy_adult_pediatric.h5ad")
gca_out2 = os.path.join(PATH, "gca_ref/obj_healthy_adult_pediatric_TIL.h5ad")

diagnoses_to_extract = ['Healthy adult','Pediatric healthy']

adata = sc.read_h5ad(gca_in)
adata = adata[adata.obs['Diagnosis'].isin(diagnoses_to_extract)]
adata.write_h5ad(gca_out1)

# regions_to_extract = ['DUO', 'DCL', 'TIL', 'ACL', 'APD', 'CAE', 'TCL', 'ILE2', 'SCL', 'ILE1', 'ILE', 'JEJ', 'REC', 'MLN'] 
regions_to_extract = ['TIL']

adata_filt = adata[adata.obs['Region code'].isin(regions_to_extract)]
adata_filt.write_h5ad(gca_out2)