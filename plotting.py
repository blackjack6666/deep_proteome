import numpy as np


def label_stats(x_pos,y_pos, y_offset, p_val, ax, correction=False,number_compare=1):
    """
    label the statistics as stars
    :param x_pos: [x1, x2]
    :param y_pos: y
    :param p_val: test stats
    :return:
    """

    if correction and number_compare != 1:
        print ('Correct p value for multiple comparison, number of comparion: %i' % number_compare)
        p_val = np.min(np.array([p_val*number_compare,1]))
    print ("adjusted p val: %f" % p_val)

    if p_val <= 1.00e-04:
        text = '****'
    elif 1.00e-04 < p_val <= 1.00e-03:
        text = '***'
    elif 1.00e-03 < p_val <= 1.00e-02:
        text = '**'
    elif 1.00e-02 < p_val <= 5.00e-02:
        text = '*'
    else:
        text = 'ns'
    x1, x2 = x_pos
    line_x, line_y = [x1, x1, x2, x2], [y_pos, y_pos+y_offset, y_pos +y_offset, y_pos]
    ax.plot(line_x, line_y,color='0.2', lw=1.5)
    ann = ax.annotate(
        text, xy=(np.mean([x1, x2]), y_pos+y_offset),
        xytext=(0, 1), textcoords='offset points',
        xycoords='data', ha='center', va='bottom')
    return ax


def volcano_plot(df,lfc,pv,
                 lfc_thr=(1, 2), pv_thr=(0.05, 0.01),
                 color=("#00239CFF", "grey", "#E10600FF")):
    """

    :param df: a pandas dataframe
    :param lfc: column name for logfold change
    :param pv: column name for p value
    :return:
    """
    from bioinfokit import analys, visuz
    visuz.GeneExpression.volcano(df=df,
                                 lfc=lfc,
                                 pv=pv,
                                 lfc_thr=lfc_thr,
                                 pv_thr=pv_thr,
                                 color=color,
                                 sign_line=True,
                                 show=True)


def extract_clustered_table(res, data):
    """
    input
    =====
    res:     <sns.matrix.ClusterGrid>  the clustermap object
    data:    <pd.DataFrame>            input table

    output
    ======
    returns: <pd.DataFrame>            reordered input table
    """

    # if sns.clustermap is run with row_cluster=False:
    if res.dendrogram_row is None:
        print("Apparently, rows were not clustered.")
        return -1

    if res.dendrogram_col is not None:
        # reordering index and columns
        new_cols = data.columns[res.dendrogram_col.reordered_ind]
        new_ind = data.index[res.dendrogram_row.reordered_ind]

        return data.loc[new_ind, new_cols]

    else:
        # reordering the index
        new_ind = data.index[res.dendrogram_row.reordered_ind]

        return data.loc[new_ind, :]


def text_cloud(text_file,stop_words=None,out_put_png=None):
    from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator
    import matplotlib.pyplot as plt
    with open(text_file,'r',encoding="utf8") as f_o:
        f_string = f_o.read()
    wordcloud = WordCloud(stopwords=stop_words, background_color="white").generate(f_string).to_array()
    print (wordcloud)
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis("off")
    plt.show()

    if out_put_png:
        wordcloud.to_file(out_put_png)


if __name__=='__main__':
    from nltk.corpus import stopwords
    import nltk
    # nltk.download('stopwords')

    test_file = 'F:/matrisomedb2.0/test.txt'
    stopWords = set(stopwords.words('english'))

    text_cloud(test_file,stop_words=stopWords)
