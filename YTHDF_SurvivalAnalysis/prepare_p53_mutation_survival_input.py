import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks", palette="bright")

from scipy.stats import wilcoxon,ttest_ind,mannwhitneyu

gene_list = {"YTHDF2":"ENSG00000198492", 'TP53':"ENSG00000141510",'CDKN1A':"ENSG00000124762",'MDM2':"ENSG00000135679",
             'PUMA':"ENSG00000105327",'GAPDH':"ENSG00000111640","ACTB":"ENSG00000075624","TUBB":"ENSG00000196230",
             "LMNB1":"ENSG00000113368","RPL17":"ENSG00000265681","VCL":"ENSG00000035403",'KDM3A':"ENSG00000115548"}

gene_ensemble = ["ENSG00000198492", "ENSG00000141510","ENSG00000124762","ENSG00000135679",
                 "ENSG00000105327","ENSG00000111640","ENSG00000075624","ENSG00000196230",
                 "ENSG00000113368","ENSG00000265681","ENSG00000035403", "ENSG00000115548"]

gene_ensemble = ["ENSG00000198492"]

cancer_ls = ["lgg","gbm"]

cancer_mutation = {}
for cancer in cancer_ls:
    print(cancer)
    if cancer == "laml":
        mut_file = "../TCGA/data/"+cancer+"/mutation_wustl.tsv"
    else:
        mut_file = "../TCGA/data/"+cancer+"/mc3.txt.tsv"
    mutation_map = {}
    all_mut_sample = []
    mis_mut_sample = []
    GR_mut_sample = []

    with open(mut_file,"r") as f:
        next(f)
        for line in f:
            line_ls = line.strip().split("\t")
            sample = line_ls[0]
            gene = line_ls[6]
            mut_type = line_ls[7]
            aac = line_ls[8]

            if gene == "TP53":
                if sample not in all_mut_sample:
                    all_mut_sample.append(sample)

                if "missense" in mut_type.lower() and (sample not in mis_mut_sample):
                    mis_mut_sample.append(sample)

                if ("245" in aac) or ("175" in aac) or ("248" in aac) or ("249" in aac) or ("273" in aac) or ("282" in aac):
                    if sample not in GR_mut_sample:
                        GR_mut_sample.append(sample)

                if aac in mutation_map:
                    mutation_map[aac].append(sample)
                else:
                    mutation_map[aac] = [sample]
    # print(len(all_mut_sample),len(mis_mut_sample),len(GR_mut_sample))
    cancer_mutation[cancer]=[all_mut_sample,mis_mut_sample, GR_mut_sample]

for cancer in cancer_ls:
    expr_file = "../TCGA/data/"+cancer+"/htseq_fpkm-uq.tsv"
    RNAseq = pd.read_csv(expr_file, header=0, index_col=0, sep="\t")

    new_sample_name = {x: x[:-1] for x in RNAseq.columns}
    new_index_name = {x: x[:15] for x in RNAseq.index}
    RNAseq.rename(columns=new_sample_name, index=new_index_name, inplace=True)

    samples_tumor = RNAseq.loc[gene_ensemble,~RNAseq.columns.str.endswith("-11")]

    samples_tumor_ls = list(samples_tumor.columns)

    all_mut_samples = cancer_mutation[cancer][0]
    mis_mut_samples = cancer_mutation[cancer][1]
    GR_mut_samples = cancer_mutation[cancer][2]

    all_mut_sample_expr_ls = [x for x in all_mut_samples if x in samples_tumor_ls]
    mis_mut_sample_expr_ls = [x for x in mis_mut_samples if x in samples_tumor_ls]
    GR_mut_sample_expr_ls = [x for x in GR_mut_samples if x in samples_tumor_ls]

    wt_sample_expr = samples_tumor.drop(all_mut_sample_expr_ls, axis=1)
    wt_samples = list(wt_sample_expr.columns)

    ####################################
    used_mut_expr = all_mut_sample_expr_ls

    sub_wt_df = samples_tumor.loc[gene_ensemble,wt_samples]
    sub_mut_df = samples_tumor.loc[gene_ensemble, used_mut_expr]

    print("wt vs mut:", sub_wt_df.shape[1],sub_mut_df.shape[1])

    sub_wt_df.sort_values(by=gene_ensemble,axis=1,inplace=True)
    sub_mut_df.sort_values(by=gene_ensemble, axis=1, inplace=True)

    percent = 0.3
    sub_wt_high_df = sub_wt_df.iloc[:,-int(sub_wt_df.shape[1] * percent):]
    sub_wt_low_df = sub_wt_df.iloc[:, :int(sub_wt_df.shape[1] * percent)]
    print("WT high: WT low:", sub_wt_high_df.shape[1],sub_wt_low_df.shape[1])

    sub_mut_high_df = sub_mut_df.iloc[:, -int(sub_mut_df.shape[1] * percent):]
    sub_mut_low_df = sub_mut_df.iloc[:, :int(sub_mut_df.shape[1] * percent)]
    print("Mut high: Mut low:", sub_mut_high_df.shape[1], sub_mut_low_df.shape[1])

    wt_high_sample_ls = sub_wt_high_df.columns.tolist()
    wt_low_sample_ls = sub_wt_low_df.columns.tolist()

    mut_high_sample_ls = sub_mut_high_df.columns.tolist()
    mut_low_sample_ls = sub_mut_low_df.columns.tolist()

    fp = open("results/"+cancer+"_survival_input.tsv","w")
    with open("../TCGA/data/"+cancer+"/survival.tsv","r") as f:
        fp.write(f.readline())
        for line in f:
            line_ls = line.strip().split("\t")
            sample = line_ls[0][:-1]
            i = 0
            if sample in wt_high_sample_ls:
                line_ls += ["WT","High"]
                i+=1
            if sample in wt_low_sample_ls:
                line_ls += ["WT","Low"]
                i+=2
            if sample in mut_high_sample_ls:
                line_ls += ["Mut","High"]
                i+=4
            if sample in mut_low_sample_ls:
                line_ls += ["Mut","low"]
                i+=8
            if i not in [0,1,2,4,8]:
                print(i,sample)

            if i != 0:
                fp.write("\t".join(line_ls) + "\n")
    fp.close()





