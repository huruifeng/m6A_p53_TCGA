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

gene_list = {"KDM3A":"ENSG00000115548"}
gene_ensemble = ["ENSG00000115548"]
gene_symbol = ["KDM3A"]

dir = "F:/TCGA_Xena2/data"
cancer_ls = []
for dir_i in os.listdir(dir):
    cancer_ls.append(dir_i)

print(len(cancer_ls))

category = "NULL"

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

    p53R175H_sample = []
    p53G245D_sample = []
    p53R282W_sample = []
    p53R273H_sample = []
    p53R248W_sample = []

    p53null_sample = []

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

            # if "silent" in mut_type.lower(): continue

            if gene == "TP53":
                if "silent" in mut_type.lower():
                    p53silent_sample.append(sample)
                    continue

                if sample not in all_mut_sample:
                    all_mut_sample.append(sample)

                if "r175h" in aac.lower():
                    p53R175H_sample.append(sample)
                elif "g245d" in aac.lower():
                    p53G245D_sample.append(sample)
                elif "r282w" in aac.lower():
                    p53R282W_sample.append(sample)
                elif "r273h" in aac.lower():
                    p53R273H_sample.append(sample)
                elif "r248w" in aac.lower():
                    p53R248W_sample.append(sample)
                elif "shift_ins" in mut_type.lower() or "shift_del" in mut_type.lower() or "large" in mut_type.lower() or "nonsense" in mut_type.lower():
                    p53null_sample.append(sample)
                else:
                    p53othermut_sample.append(sample)


    cancer_mutation[cancer]=[all_mut_sample, p53silent_sample,
                             p53R175H_sample, p53G245D_sample,p53R282W_sample,p53R273H_sample, p53R248W_sample,
                             p53null_sample, p53othermut_sample]


#####################
for mut_i in ["R175H","G245D","R282W","R273H","R248W"]:

    fp = open("results/Mutation_count_TCGA.txt", "w")
    fp.write("cancer\tall_n\tnormal_n\tcancer_n\tp53wt_n\tp53mut_n\tp53silient_n\t"
             "r175h\tg245d\tr282w\tr273h\tr248w\tp53_null(Ins/Del/Nonsense)\tp53Other_mut\n")

    pan_df = pd.DataFrame(data={'gene': [], 'case_type': [], 'p': [], 'expr': [], "cancer": []})
    with PdfPages('results/p53_'+mut_i+'vs'+category+'.pdf') as pdf:
        fout = open("results/p53_"+mut_i+"_vs"+category+"_table.txt", "w")
        for cancer in cancer_ls:
            print("==============================================")
            print(cancer)
            expr_file = dir+"/"+cancer+"/HiSeqV2.tsv"
            RNAseq = pd.read_csv(expr_file, header=0, index_col=0, sep="\t")

            if cancer in ["laml"]:
                samples_tumor = RNAseq.loc[gene_symbol, RNAseq.columns.str.endswith("-03")]
                samples_normal = RNAseq.loc[gene_symbol, ~RNAseq.columns.str.endswith("-03")]
            elif cancer in ["skcm"]:
                samples_tumor = RNAseq.loc[gene_symbol,  RNAseq.columns.str.endswith("-01") | RNAseq.columns.str.endswith("-06")]
                samples_normal = RNAseq.loc[gene_symbol, ~ (RNAseq.columns.str.endswith("-01") | RNAseq.columns.str.endswith("-06"))]
            else:
                samples_tumor = RNAseq.loc[gene_symbol, RNAseq.columns.str.endswith("-01")]
                samples_normal = RNAseq.loc[gene_symbol, ~RNAseq.columns.str.endswith("-01")]

            samples_tumor_ls = list(samples_tumor.columns)
            samples_normal_ls = list(samples_normal.columns)

            all_mut_sample = cancer_mutation[cancer][0]
            p53silent_sample = cancer_mutation[cancer][1]
            p53R175H_sample = cancer_mutation[cancer][2]
            p53G245D_sample = cancer_mutation[cancer][3]
            p53R282W_sample = cancer_mutation[cancer][4]
            p53R273H_sample = cancer_mutation[cancer][5]
            p53R248W_sample = cancer_mutation[cancer][6]

            p53null_sample = cancer_mutation[cancer][7]

            p53othermut_sample = cancer_mutation[cancer][8]


            all_mut_expr = [x for x in all_mut_sample if x in samples_tumor_ls]

            p53silent_expr = [x for x in p53silent_sample if x in samples_tumor_ls]

            p53R175H_expr = [x for x in p53R175H_sample if x in samples_tumor_ls]
            p53G245D_expr = [x for x in p53G245D_sample if x in samples_tumor_ls]
            p53R282W_expr = [x for x in p53R282W_sample if x in samples_tumor_ls]
            p53R273H_expr = [x for x in p53R273H_sample if x in samples_tumor_ls]
            p53R248W_expr = [x for x in p53R248W_sample if x in samples_tumor_ls]
            p53null_expr = [x for x in p53null_sample if x in samples_tumor_ls]

            p53othermut_expr = [x for x in p53othermut_sample if x in samples_tumor_ls]

            p53wt_df = samples_tumor.drop(all_mut_expr+p53silent_expr, axis=1)
            p53wt_expr = list(p53wt_df.columns)

            fp.write("\t".join([cancer,str(RNAseq.shape[1]), str(samples_normal.shape[1]),str(samples_tumor.shape[1]),
                               str(len(p53wt_expr)),str(len(all_mut_expr)),str(len(p53silent_expr)),
                               str(len(p53R175H_expr)),str(len(p53G245D_expr)),str(len(p53R282W_expr)),str(len(p53R273H_expr)),str(len(p53R248W_expr)),
                                str(len(p53null_sample)), str(len(p53othermut_expr))])+"\n")

            ####################################
            if mut_i == "R175H":
                if len(p53R175H_expr)==0: continue
                else:  list1 = p53R175H_expr
            if mut_i == "G245D":
                if len(p53G245D_expr)==0: continue
                else:  list1 = p53G245D_expr
            if mut_i == "R282W":
                if len(p53R282W_expr)==0: continue
                else:  list1 = p53R282W_expr
            if mut_i == "R273H":
                if len(p53R273H_expr)==0: continue
                else:  list1 = p53R273H_expr
            if mut_i == "R248W":
                if len(p53R248W_expr)==0: continue
                else:  list1 = p53R248W_expr

            if category == "WT":
                list2 = p53wt_expr
            else:
                list2 = p53null_expr

            sub_df1 = samples_tumor.loc[:, list1]
            sub_df2 = samples_tumor.loc[:, list2]

            print("p53_"+mut_i+" vs p53_"+category+":", sub_df1.shape[1],sub_df2.shape[1])


            #################
            total_df1 = sub_df1.copy()
            total_df2 = sub_df2.copy()

            gene_p = {}
            d = {'gene': [], 'case_type': [], 'p': [], 'expr': [],"cancer":[]}
            df = pd.DataFrame(data=d)

            for gene_i in gene_list:
                # if gene_i != "YTHDF2":
                #     continue
                print(gene_i)

                expr1 = total_df1.loc[gene_i,:].values
                expr2 = total_df2.loc[gene_i, :].values

                print(mut_i+" vs "+category+":", len(expr1),len(expr2))

                total_n = len(expr1) + len(expr2)

                cancer_n = [cancer] * total_n

                case_type = ['TP53_'+mut_i] * len(expr1) + ["TP53_"+category] * len(expr2)
                z_statistic, p_value = ttest_ind(expr1, expr2)
                p = [p_value] * total_n
                if p_value < 0.01:
                    gene = [gene_i + "**"] * total_n
                elif p_value < 0.05:
                    gene = [gene_i + "*"] * total_n
                else:
                    gene = [gene_i] * total_n

                expression = list(np.append(expr1, expr2))
                d = {'gene': gene, 'case_type': case_type, 'p': p, 'expr': expression,"cancer":cancer_n}
                sub_data = pd.DataFrame(data=d)
                df = df.append(sub_data, ignore_index=True)

                pan_d = {'gene': gene_i, 'case_type': case_type, 'p': p, 'expr': expression, "cancer": cancer_n}
                pan_sub_data = pd.DataFrame(data=pan_d)
                pan_df = pan_df.append(pan_sub_data,ignore_index=True)

                fout.write(gene_i + "\t"+mut_i+"_n="+str(len(expr1)) + "\t"+category+"_n="+str(len(expr2))+"\t" + "p(t-test)="+str(p_value) +  "\t"+cancer+"\n")

            f, ax = plt.subplots(figsize=(7, 7))
            sns.boxplot(x="gene", y="expr",
                        hue="case_type", palette=["#A3CCDE", "#F58C30"],
                        data=df)
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
    fout = open("results/PanTCGA_p53_"+mut_i+"vs"+category+"_table.txt", "w")
    for gene_i in gene_list:
        expr1 = pan_df.loc[(pan_df["gene"]==gene_i) & (pan_df["case_type"]=="TP53_"+mut_i),"expr"]
        expr2 = pan_df.loc[(pan_df["gene"] == gene_i) & (pan_df["case_type"] == "TP53_"+category), "expr"]
        z_statistic, p_value = ttest_ind(expr1, expr2)

        fout.write( gene_i + "\t"+mut_i+"_n=" + str(len(expr1)) + "\t"+category+"_n=" + str(len(expr2)) + "\t" + "p(t-test)=" + str(
                p_value) + "\n")

        pan_df.loc[(pan_df["gene"] == gene_i), "p"] = p_value
        if p_value < 0.01:
            pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i + "**"
        elif p_value < 0.05:
            pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i + "*"
        else:
            pan_df.loc[(pan_df["gene"] == gene_i) , "gene"] = gene_i
    fout.close()

    f, ax = plt.subplots(figsize=(7, 7))
    sns.boxplot(x="gene", y="expr",
                hue="case_type", palette=["#A3CCDE", "#F58C30"],
                data=pan_df, )
    # sns.despine(trim=True) # remove frame

    ax.set(ylabel="Expression (log2(norm_count+1))")
    ax.set(title='Expression in Pan-TCGA')
    ax.legend(title="", bbox_to_anchor=(1.2, 0.5), frameon=False)

    plt.tight_layout()
    f.savefig("results/PanTCGA_p53_"+mut_i+"vs"+category+".pdf")
    plt.close()
