import cellxgene_census, numpy as np, pandas as pd, scipy.sparse as sp
with cellxgene_census.open_soma(census_version="stable") as census:
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", measurement_name="RNA",
        var_value_filter="feature_name == 'IFIH1'",
        obs_value_filter="tissue_general == 'brain' and is_primary_data == True and disease in ['normal','multiple sclerosis']",
        obs_column_names=["cell_type","disease","donor_id","dataset_id","assay","suspension_type","raw_sum","tissue"],
        X_name="raw")
print("shape",ad.shape); print(ad.var)
x=np.asarray(sp.csr_matrix(ad.X).todense()).ravel()
o=ad.obs.copy(); o['ifih1']=x; o['cp10k']=1e4*o.ifih1/o.raw_sum
o.to_parquet('census_ifih1.parquet')
print("\n--- disease counts ---"); print(o.disease.value_counts().to_string())
print("\n--- datasets with MS ---"); print(o[o.disease=='multiple sclerosis'].groupby(['dataset_id','donor_id'],observed=True).size().to_string())
