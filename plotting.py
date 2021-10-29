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