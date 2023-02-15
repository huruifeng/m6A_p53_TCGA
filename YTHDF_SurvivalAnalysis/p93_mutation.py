import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks", palette="bright")

from scipy.stats import wilcoxon,ttest_ind,mannwhitneyu

# gene_list = {"YTHDF2":"ENSG00000198492", 'TP53':"ENSG00000141510",'CDKN1A':"ENSG00000124762",'MDM2':"ENSG00000135679",
#              'PUMA':"ENSG00000105327",'GAPDH':"ENSG00000111640","ACTB":"ENSG00000075624","TUBB":"ENSG00000196230",
#              "LMNB1":"ENSG00000113368","RPL17":"ENSG00000265681","VCL":"ENSG00000035403",'KDM3A':"ENSG00000115548"}
#
# gene_ensemble = ["ENSG00000198492", "ENSG00000141510","ENSG00000124762","ENSG00000135679",
#                  "ENSG00000105327","ENSG00000111640","ENSG00000075624","ENSG00000196230",
#                  "ENSG00000113368","ENSG00000265681","ENSG00000035403", "ENSG00000115548"]


#gene_list = {"YTHDF2":"ENSG00000198492", 'CDKN1A':"ENSG00000124762",'MDM2':"ENSG00000135679",'PUMA':"ENSG00000105327","CDKN2B":"ENSG00000147883"}
#gene_ensemble = ["ENSG00000198492","ENSG00000124762","ENSG00000135679","ENSG00000105327","ENSG00000147883"]

# gene_list = {"SPOCK2":"ENSG00000107742","NAP1L2":"ENSG00000186462","MARVELD3":"ENSG00000140832"}
# gene_ensemble = ["ENSG00000107742","ENSG00000186462","ENSG00000140832"]

gene_list = {"SOX2":"ENSG00000181449"}
gene_ensemble = ["ENSG00000181449"]

dir = "../TCGA/data"
cancer_ls = []
for dir_i in os.listdir(dir):
    cancer_ls.append(dir_i)

print(len(cancer_ls))

cancer_mutation = {}
for cancer in cancer_ls:
    if cancer.lower() not in ["lgg","gbm","brca","read"]:
        continue
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

#
# total_wt_df = pd.DataFrame(columns=gene_ensemble)
# total_mut_df = pd.DataFrame(columns=gene_ensemble)

# for cancer in cancer_mutation:

# cancer="lgg"
for cancer in ["lgg","gbm","brca","read"]:
    print("==============================================")
    total_wt_df = pd.DataFrame(columns=gene_ensemble)
    total_mut_df = pd.DataFrame(columns=gene_ensemble)

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

    all_mut_sample_expr = [x for x in all_mut_samples if x in samples_tumor_ls]
    mis_mut_sample_expr = [x for x in mis_mut_samples if x in samples_tumor_ls]
    GR_mut_sample_expr = [x for x in GR_mut_samples if x in samples_tumor_ls]

    sample_wt_expr = samples_tumor.drop(all_mut_sample_expr, axis=1)
    wt_samples_expr = sample_wt_expr.columns

    if float(len(all_mut_sample_expr)) / len(samples_tumor_ls) < 0.05:
        pass
    else:
        print(cancer,float(len(all_mut_sample_expr)) / len(samples_tumor_ls) )

    ####################################
    used_mut_expr = GR_mut_sample_expr

    sub_wt_df = samples_tumor.loc[:,wt_samples_expr]
    sub_mut_df = samples_tumor.loc[:, used_mut_expr]

    print("wt vs hotspots:", sub_wt_df.shape[1],sub_mut_df.shape[1])

    # total_wt_df = total_wt_df.append(sub_wt_df.T)
    # total_mut_df = total_mut_df.append(sub_mut_df.T)

    #################
    total_wt_df = sub_wt_df.T
    total_mut_df = sub_mut_df.T

    gene_p = {}
    d = {'gene': [], 'case_type': [], 'p': [], 'expr': []}
    df = pd.DataFrame(data=d)
    fout = open(cancer+"_hotspots_p_table.txt", "w")
    for gene_i in gene_list:
        # if gene_i != "YTHDF2":
        #     continue
        print(gene_i)
        gene_x_i = gene_list[gene_i]


        mut_expr = total_mut_df.loc[:,gene_x_i].values
        wt_expr = total_wt_df.loc[:, gene_x_i].values

        print("wt vs hotspots:", len(wt_expr),len(mut_expr))

        total_n = len(wt_expr) + len(mut_expr)

        case_type = ['TP53_WT'] * len(wt_expr) + ["TP53_hotspots"] * len(mut_expr)
        z_statistic, p_value = mannwhitneyu(wt_expr, mut_expr)
        p = [p_value] * total_n
        if p_value < 0.01:
            gene = [gene_i + "**"] * total_n
        elif p_value < 0.05:
            gene = [gene_i + "*"] * total_n
        else:
            gene = [gene_i] * total_n

        expression = list(np.append(wt_expr, mut_expr))
        d = {'gene': gene, 'case_type': case_type, 'p': p, 'expr': expression}
        df = pd.DataFrame(data=d)
        # df = df.append(sub_df, ignore_index=True)
        fout.write(gene_i + "\tWT_n="+str(len(wt_expr)) + "\tHotspots_n="+str(len(mut_expr))+"\t" + "Mann–Whitney U test p="+str(p_value) + "\n")


        # Draw a nested boxplot to show bills by day and time
        f, ax = plt.subplots(figsize=(7, 7))
        sns.boxplot(x="gene", y="expr",
                    hue="case_type", palette=["#A3CCDE", "#F58C30"],
                    data=df, )
        # sns.despine(trim=True) # remove frame

        ax.set(ylabel="Expression (log(FPKM+1))")
        if p_value < 0.01:
            ax.set(xlabel = gene_i+'**')#
        elif p_value < 0.05:
            ax.set(xlabel = gene_i + '*')  #
        else:
            ax.set(xlabel = gene_i)  #
        ax.set_xticks([])
        ax.set(title=gene_i + ' expression in TCGA - ' + cancer.upper())
        ax.legend(title="", bbox_to_anchor=(1.5, 0.5), frameon=False)

        # ax.set_xlabel('',fontsize = 20) #xlabel
        # ax.set_ylabel('', fontsize = 20)#ylabel
        # if p_value < 0.01:
        #     ax.set_title(gene_i+'**', fontsize = 20)#
        # elif p_value < 0.05:
        #     ax.set_title(gene_i + '*', fontsize=20)  #
        # else:
        #     ax.set_title(gene_i, fontsize=20)  #
        # ax.tick_params(axis='y', labelsize=16)
        # ax.set_xticks([])
        # ax.get_legend().remove()
        ##########


        # for i in range(0,2*len(gene_list),2):
        #     # Select which box you want to change
        #     mybox = ax.artists[i]
        #     # Change the appearance of that box
        #     mybox.set_facecolor('#A3CCDE')

        #     mybox = ax.artists[i+1]
        #     # Change the appearance of that box
        #     mybox.set_facecolor("#F58C30")
        plt.tight_layout()
        f.savefig(cancer + "_hotspots_"+gene_i+".pdf")
        plt.close()


    # # Draw a nested boxplot to show bills by day and time
    # f, ax = plt.subplots(figsize=(4.5, 6))
    # sns.boxplot(x="gene", y="expr",
    #             hue="case_type", palette=["#A3CCDE", "#F58C30"],
    #             data=df, )
    #
    # ax.set_xlabel('',fontsize = 20) #xlabel
    # ax.set_ylabel('', fontsize = 20)#ylabel
    # if p_value < 0.01:
    #     ax.set_title(gene_i+'**', fontsize = 20)#
    # elif p_value < 0.05:
    #     ax.set_title(gene_i + '*', fontsize=20)  #
    # else:
    #     ax.set_title(gene_i, fontsize=20)  #
    # ax.tick_params(axis='y', labelsize=16)
    # ax.set_xticks([])
    # # ax.get_legend().remove()
    # f.savefig(cancer + "_hotspots_"+gene_i+"_legend.pdf")
    # plt.close()



    fout.close()

