import csv
import pandas as pd

df = pd.read_csv('nGRE_parse_output/ADX_F_DexvsADX_F_Veh/downregulated_gene_output.csv')

for column in range(10):
    test = df['mutations'][column]
    print(test)

print(df)
