import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks", palette="bright")

from scipy.stats import wilcoxon,ttest_ind,mannwhitneyu

gene_list = {
    # "ACBD4": "ENSG00000181513",
    # "ARMC3": "ENSG00000165309",
    # "BAIAP3": "ENSG00000007516",
    "BMP4": "ENSG00000125378",
    "BMP7": "ENSG00000101144",
    # "CALY":"ENSG00000130643",
    "CDKN2A":"ENSG00000147889",
    "CDKN2B":"ENSG00000147883",
    # "DNAH6":"ENSG00000115423",
    # "EFHB":"ENSG00000163576",
    "FLRT3": "ENSG00000125848",
    "GRIA1":"ENSG00000155511",
    "GRID1": "ENSG00000182771",
    "IRX2": "ENSG00000170561",
    "KCTD11":"ENSG00000213859",
    # "LRRC4C": "ENSG00000148948",
    "MARVELD1": "ENSG00000155254",
    # "MARVELD3": "ENSG00000140832",
    "MFSD3":"ENSG00000167700",
    # "MKX": "ENSG00000150051",
    # "NAP1L2":"ENSG00000186462",
    # "NKAIN4": "ENSG00000101198",
    "NTN1":"ENSG00000065320",
    # "NWD1":"ENSG00000188039",
    # "PRRG1": "ENSG00000130962",
    "PTPRU":"ENSG00000060656",
    "RALGPS2": "ENSG00000116191",
    # "RGL3":"ENSG00000205517",
    "RNF150": "ENSG00000170153",
    # "SLC20A2": "ENSG00000168575",
    "SLC2A12":"ENSG00000146411",
    # "SPAG6": "ENSG00000077327",
    "SPOCK2":"ENSG00000107742",
    "SULF2": "ENSG00000196562",
    # "SULT1C4":"ENSG00000198075",
    "TMEM121":"ENSG00000184986",
    # "TNXB":"ENSG00000168477",
    # "VWA3A":"ENSG00000175267",
    # "WNT3": "ENSG00000108379",
    "ZMAT3":"ENSG00000172667",
    "ZNRF3":"ENSG00000183579"
}

gene_ls = list(gene_list.values())

lgg_expr_df = pd.read_csv("../TCGA/data/lgg/htseq_fpkm-uq.tsv",index_col=0,header=0,sep="\t")
index_ls = [i[:15] for i in lgg_expr_df.index]
lgg_expr_df.index = index_ls
lgg_expr_genes = lgg_expr_df.loc[gene_ls,:]

gbm_expr_df = pd.read_csv("../TCGA/data/gbm/htseq_fpkm-uq.tsv",index_col=0,header=0,sep="\t")
index_ls = [i[:15] for i in gbm_expr_df.index]
gbm_expr_df.index = index_ls
gbm_expr_genes = gbm_expr_df.loc[gene_ls,:]

lgg_survival_df = pd.read_csv("../TCGA/data/lgg/survival.tsv",index_col=0,header=0,sep="\t")
gbm_survival_df = pd.read_csv("../TCGA/data/gbm/survival.tsv",index_col=0,header=0,sep="\t")
survival_df = lgg_survival_df.append(gbm_survival_df)
samples = survival_df.index.tolist()

combined_df = pd.merge(lgg_expr_genes,gbm_expr_genes,how="left",left_index=True,right_index=True)
combined_df = combined_df.loc[:,combined_df.columns.isin(samples)]

combined_df_zscore = combined_df.apply(lambda x: (x-x.mean())/x.std(),axis=1)
combined_df_zscore_avg = combined_df.mean(0)
combined_df_zscore_avg = combined_df_zscore_avg.to_frame(name="score")
score_median = combined_df_zscore_avg["score"].median()
combined_df_zscore_avg.loc[:,"Group"] = "High"
combined_df_zscore_avg.loc[combined_df_zscore_avg["score"]<=score_median,"Group"] = "Low"

survival_data = pd.merge(combined_df_zscore_avg,survival_df,how="left",right_index=True,left_index=True)
survival_data.to_csv("survival_input21.tsv",index=True,index_label="sample",sep="\t")
