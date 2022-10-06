import numpy as np
import scipy.stats
from scipy.stats import poisson
import pandas as pd


def ml_TPT(gene_count, total_transcripts):
    # total_transcripts = 100626
    bins = 10001
    expression_list = np.linspace(0, 1, num=bins)
    unnorm = poisson.pmf(k=gene_count, mu=expression_list * total_transcripts)
    normalized_dist = unnorm / unnorm.sum()
    _ = pd.DataFrame({"TPT": normalized_dist, 'probability': expression_list})
    TPT = round(sum(_['TPT'] * _['probability']) * 10000,3)
    return TPT


df = pd.read_csv(f'astrovirus_counts_nobystanders.csv', index_col=0)

gene_list = df.columns
counter = 0
# TPT_df = pd.DataFrame(columns=gene_list)
TPT_df = pd.DataFrame(columns=df.columns)
for index, row in df.iterrows():
    counter += 1
    print(counter, index)
    total_transcripts = row.sum()
    tmp = row[gene_list]
    # zero_TPT = round(ml_TPT(0, total_transcripts),2)
    TPTs = [round(ml_TPT(x, total_transcripts),2) if x > 0 else 0 for x in tmp]
    # TPTs.insert(0, int(total_transcripts))
#
    # df.loc[index][1:] = TPTs
    TPT_df.loc[index] = TPTs
    # print(TPTs)
    # if counter > 10:
    #     break
TPT_df.sort_values(by='Astrovirus', inplace=True)
#
# #
TPT_df.insert(loc = 0,
          column = 'log2_Astrovirus',
          value = np.log2(TPT_df['Astrovirus']+1).round(decimals=1))
# #
TPT_df.to_csv(f'Expected_TPT_matrix_Astrovirus_v2.csv')