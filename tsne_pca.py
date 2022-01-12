import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
"""
df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/18_2_accumulated_cov_transpose.csv')

print (len(df['category'].unique()))
cate_number_dict = {category:num for category, num in zip(df['category'].unique(),range(6))}

sub_df = df.iloc[:,1:-1]
pca = PCA(n_components=3)
pca_result = pca.fit_transform(sub_df)
df['pca-one'] = pca_result[:,0]
df['pca-two'] = pca_result[:,1]
df['pca-three'] = pca_result[:,2]
df['color_number'] = [cate_number_dict[cat] for cat in df['category']]
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

plt.figure(figsize=(10,5))
sns.scatterplot(
    x="pca-one", y="pca-two",
    hue="category",
    palette=sns.color_palette("hls", 6),
    data=df,
    legend="full",
    alpha=0.3
)

# ax = plt.figure(figsize=(10,5)).gca(projection='3d')
# ax.scatter(
#     xs=df["pca-one"],
#     ys=df["pca-two"],
#     zs=df["pca-three"],
#     c=df["color_number"],
#     cmap='coolwarm'
# )
# ax.set_xlabel('pca-one')
# ax.set_ylabel('pca-two')
# ax.set_zlabel('pca-three')
plt.show()
"""
from os import rename
from glob import glob
folder_path = 'D:/data/Naba_deep_matrisome/proteomeXchange_submission/*/'
folders_peptide = glob(folder_path+'peptide.tsv')
folders_psms = glob(folder_path+'psm.tsv')
for each in folders_psms:
    sample = each.split('\\')[-2]
    rename(each,each.replace('psm.tsv',sample+'_'+'psm.tsv'))
