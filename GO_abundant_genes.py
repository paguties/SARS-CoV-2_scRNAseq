import pandas as pd
import numpy as np
import scipy.stats as stats
from collections import Counter


GO_dictionary = {}
with open('GO_biological_process_Hsapiens.tsv') as file:
    next(file)
    for line in file:
        line = line.rstrip()
        gene = line.split('\t')[0]
        GO_terms = line.split('\t')[1:]
        GO_dictionary[gene] = GO_terms
        # print(line)



cells_list = ['hBECs','colon','ileum']


GO_counts = {}
for cells in cells_list:
    # cells = 'hBECs'
    # df = pd.read_csv(f'DEGS_{cells}.csv')
    df = pd.read_csv(f'gene_abundance_mock_{cells}.csv')
    # selection = (df['comparison'] == 'phaseI') & (df['log2FC'] > 0)
    df.sort_values(by='corrected_mean', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df = df.iloc[:20, :]
    # selected_ISGs = df.columns[df.columns.isin(ISG_list)]
    # selected_ISGs  = list(df[df['gene'].isin(ISG_list)]['gene'])


    for index, row in df.iterrows():

        gene = row['gene']
        # log2FC = row['log2FC']
        # regulation = ""

        try:
            for item in GO_dictionary[gene]:
                item = item.strip()
                if item not in GO_counts:
                    GO_counts[item] = [cells]
                else:
                    GO_counts[item].append(cells)
                # print(cells, item)
        except:
            pass

print('GO','hBECs','colon','ileum', sep="\t")
for k,v in GO_counts.items():
    hBECs_counts = v.count('hBECs')
    colon_counts = v.count('colon')
    ileum_counts = v.count('ileum')
    print(k,hBECs_counts,colon_counts, ileum_counts,  sep="\t")