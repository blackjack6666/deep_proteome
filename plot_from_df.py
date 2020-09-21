import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_excel("D:/data/deep_proteome/9_20_identified_cov.xlsx")
df['protein_num'] = range(len(df))
df_split = np.split(df,[200],axis=0)[0]



fig, ax = plt.subplots(1,1, figsize=(50,25))
for each_cond in df_split.columns[1:-1]:
    print (each_cond)
    ax.plot(df_split['protein_num'],df_split[each_cond], '-', c=np.random.rand(3,), label=each_cond)

plt.show()