import pandas as pd
import numpy as np

### filter out rows with more than 3 nan values

df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx',index_col=0)
df_sasa_clean = pd.DataFrame(columns=df_sasa.columns)

for row in df_sasa.iterrows():
    values = row[1].to_numpy()

    nan_number = np.count_nonzero(np.isnan(values))
    if nan_number >=3:
        continue
    else:
        df_sasa_clean.at[row[0],:] = values

### fill nan with nearby values

df_sasa_clean = df_sasa_clean.T.bfill().ffill().T # fill with forward and backward values
print (df_sasa_clean)
df_sasa_clean.to_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique_clean_fill.xlsx')