import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

sns.set(style="ticks", palette="bright")

from scipy.stats import wilcoxon,ttest_ind,mannwhitneyu

# gene_list = {"SPOCK2":"ENSG00000107742","NAP1L2":"ENSG00000186462","MARVELD3":"ENSG00000140832"}
# gene_ensemble = ["ENSG00000107742","ENSG00000186462","ENSG00000140832"]

gene_list = {"CCNG2":"CCNG2", "KMT5B":"SUV420H1", "DICER1":"DICER1","SWSAP1":"C19orf39"}
gene_ensemble = ["ENSG00000138764","ENSG00000110066","ENSG00000100697","ENSG00000173928"]
gene_symbol = ["CCNG2", "SUV420H1", "DICER1","C19orf39"]

dir = "F:/TCGA_Xena2/data"
cancer_ls = []
for dir_i in os.listdir(dir):
    cancer_ls.append(dir_i)

print(len(cancer_ls))

cancer_mutation = {}
for cancer in cancer_ls:
    print("=================================")
    print(cancer)
    if cancer == "laml":
        mut_file = dir+"/"+cancer+"/mutation_wustl.tsv"
    else:
        mut_file = dir+"/"+cancer+"/mc3.txt.tsv"

    all_mut_sample = []

    p53silent_sample = []

    p53r280t_sample = []
    p53del_sample = []
    p53othermut_sample = []
    p53wt_sample = []

    with open(mut_file,"r") as f:
        next(f)
        for line in f:
            line_ls = line.strip().split("\t")
            sample = line_ls[0]
            gene = line_ls[6]
            mut_type = line_ls[7]
            aac = line_ls[8]

            if gene == "TP53":
                if "silent" in mut_type.lower():
                    p53silent_sample.append(sample)
                    continue

                if sample not in all_mut_sample:
                    all_mut_sample.append(sample)

                if "r280t" in aac.lower():
                    p53r280t_sample.append(sample)
                elif "shift_ins" in mut_type.lower() or "shift_del" in mut_type.lower() or "large" in mut_type.lower() or "nonsense" in mut_type.lower():
                    p53del_sample.append(sample)
                else:
                    p53othermut_sample.append(sample)

    cancer_mutation[cancer]=[all_mut_sample,p53silent_sample, p53r280t_sample, p53del_sample, p53othermut_sample]

category = "WT"

fp = open("results/Mutation_count_TCGA.txt","w")
fp.write("cancer\tall_n\tnormal_n\tcancer_n\tp53wt_n\tp53mut_n\tp53silent_n\tp53r280t_n\tp53del_n\tp53othermut\n")

pan_df =pd.DataFrame(data= {'gene': [], 'case_type': [], 'p': [], 'expr': [],"cancer":[]})
with PdfPages('results/p53_R280Tvs'+category+'.pdf') as pdf:
    fout = open("results/p53_R280Tvs" + category + "_table.txt", "w")
    for cancer in cancer_ls:
        print("==============================================")
        expr_file = dir+"/"+cancer+"/HiSeqV2.tsv"
        RNAseq = pd.read_csv(expr_file, header=0, index_col=0, sep="\t")

        if cancer in ["laml"]:
            samples_tumor = RNAseq.loc[gene_symbol, RNAseq.columns.str.endswith("-03")]
            samples_normal = RNAseq.loc[gene_symbol, ~RNAseq.columns.str.endswith("-03")]
        elif cancer in ["skcm"]:
            samples_tumor = RNAseq.loc[gene_symbol, RNAseq.columns.str.endswith("-01") | RNAseq.columns.str.endswith("-06")]
            samples_normal = RNAseq.loc[
                gene_symbol, ~ (RNAseq.columns.str.endswith("-01") | RNAseq.columns.str.endswith("-06"))]
        else:
            samples_tumor = RNAseq.loc[gene_symbol, RNAseq.columns.str.endswith("-01")]
            samples_normal = RNAseq.loc[gene_symbol, ~RNAseq.columns.str.endswith("-01")]

        samples_tumor_ls = list(samples_tumor.columns)
        samples_normal_ls = list(samples_normal.columns)

        all_mut_sample = cancer_mutation[cancer][0]
        p53silent_sample = cancer_mutation[cancer][1]
        p53r280t_sample = cancer_mutation[cancer][2]
        p53del_sample = cancer_mutation[cancer][3]
        p53othermut_sample = cancer_mutation[cancer][4]

        all_mut_expr = [x for x in all_mut_sample if x in samples_tumor_ls]
        p53silent_expr = [x for x in p53silent_sample if x in samples_tumor_ls]
        p53r280t_expr = [x for x in p53r280t_sample if x in samples_tumor_ls]
        p53del_expr = [x for x in p53del_sample if x in samples_tumor_ls]
        p53othermut_expr = [x for x in p53othermut_sample if x in samples_tumor_ls]

        p53wt_df = samples_tumor.drop(all_mut_expr+p53silent_expr, axis=1)
        p53wt_expr = list(p53wt_df.columns)

        fp.write("\t".join([cancer,str(RNAseq.shape[1]), str(samples_normal.shape[1]),str(samples_tumor.shape[1]),
                           str(len(p53wt_expr)),str(len(all_mut_expr)),str(len(p53silent_expr)),
                           str(len(p53r280t_expr)),str(len(p53del_expr)),str(len(p53othermut_expr))])+"\n")

        if len(p53r280t_expr)==0: continue
        ####################################
        list1 = p53r280t_expr

        if category == "WT":
            list2 = p53wt_expr
        elif category == "NULL":
            list2 = p53del_expr
        else:
            list2 = p53othermut_expr

        sub_df1 = samples_tumor.loc[:, list1]
        sub_df2 = samples_tumor.loc[:, list2]

        print("p53_R280T vs p53_"+category+":", sub_df1.shape[1],sub_df2.shape[1])


        #################
        total_df1 = sub_df1.T
        total_df2 = sub_df2.T

        gene_p = {}
        d = {'gene': [], 'case_type': [], 'p': [], 'expr': [],"cancer":[]}
        df = pd.DataFrame(data=d)
        for gene_i in gene_list:
            # if gene_i != "YTHDF2":
            #     continue
            gene_x_i = gene_list[gene_i]

            print(gene_i)

            expr1 = total_df1.loc[:,gene_x_i].values
            expr2 = total_df2.loc[:, gene_x_i].values

            print("R280T vs "+category+":", len(expr1),len(expr2))

            total_n = len(expr1) + len(expr2)

            cance_n = [cancer] * total_n

            case_type = ['TP53_R280T'] * len(expr1) + ["TP53_"+category] * len(expr2)
            z_statistic, p_value = ttest_ind(expr1, expr2)
            p = [p_value] * total_n
            if p_value < 0.01:
                gene = [gene_i + "**"] * total_n
            elif p_value < 0.05:
                gene = [gene_i + "*"] * total_n
            else:
                gene = [gene_i] * total_n

            expression = list(np.append(expr1, expr2))
            d = {'gene': gene, 'case_type': case_type, 'p': p, 'expr': expression,"cancer":cance_n}
            sub_data = pd.DataFrame(data=d)
            df = df.append(sub_data, ignore_index=True)

            pan_d = {'gene': gene_i, 'case_type': case_type, 'p': p, 'expr': expression, "cancer": cance_n}
            pan_sub_data = pd.DataFrame(data=pan_d)
            pan_df = pan_df.append(pan_sub_data,ignore_index=True)

            fout.write(gene_i + "\tR280T_n="+str(len(expr1)) + "\t"+category+"_n="+str(len(expr2))+"\t" + "p(t-test)="+str(p_value) +  "\t"+cancer+"\n")

        f, ax = plt.subplots(figsize=(12, 7))
        sns.boxplot(x="gene", y="expr",
                    hue="case_type", palette=["#A3CCDE", "#F58C30"],
                    data=df, )
        # sns.despine(trim=True) # remove frame

        ax.set(ylabel="Expression (log2(norm_count+1))")
        ax.set(title='Expression in TCGA - ' + cancer.upper())
        ax.legend(title="", bbox_to_anchor=(1.2, 0.5), frameon=False)


        plt.tight_layout()
        pdf.savefig(f)
        plt.close()

    fout.close()
fp.close()

## pan cancer
fout = open("results/Pan_R280Tvs"+category+"_table.txt", "w")
for gene_i in gene_list:
    expr1 = pan_df.loc[(pan_df["gene"]==gene_i) & (pan_df["case_type"]=="TP53_R280T"),"expr"]
    expr2 = pan_df.loc[(pan_df["gene"] == gene_i) & (pan_df["case_type"] == "TP53_"+category), "expr"]
    z_statistic, p_value = ttest_ind(expr1, expr2)

    fout.write( gene_i + "\tR280T_n=" + str(len(expr1)) + "\t"+category+"_n=" + str(len(expr2)) + "\t" + "p(t-test)=" + str(
            p_value) + "\n")

    pan_df.loc[(pan_df["gene"] == gene_i), "p"] = p_value
    if p_value < 0.01:
        pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i + "**"
    elif p_value < 0.05:
        pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i + "*"
    else:
        pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i
fout.close()

f, ax = plt.subplots(figsize=(12, 7))
sns.boxplot(x="gene", y="expr",
            hue="case_type", palette=["#A3CCDE", "#F58C30"],
            data=pan_df, )
# sns.despine(trim=True) # remove frame

ax.set(ylabel="Expression (log2(norm_count+1))")
ax.set(title='Expression in Pan-TCGA')
ax.legend(title="", bbox_to_anchor=(1.2, 0.5), frameon=False)

plt.tight_layout()
f.savefig("results/Pan_R280Tvs"+category+".pdf")
plt.close()
