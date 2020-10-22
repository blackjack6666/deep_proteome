"""
different classification model implementation
"""

import pandas as pd
import numpy as np
import pickle as ppp
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn import tree


def df_dummy_getter(polymer_label_dict):
    """
    get a pandas dataframe based on polymer label dict from different experiment results, values in df are binary
    :param polymer_label_dict:
    :return:
    """
    two_d_list = []
    for polymer in polymer_label_dict:
        one_d_list = []
        for aa in polymer:
            one_d_list.append(aa)
        one_d_list.append(polymer_label_dict[polymer])
        two_d_list.append(one_d_list)

    columns = [i for i in range(31)]
    columns.append('label')
    df = pd.DataFrame(two_d_list, columns=columns)
    df_dummy = pd.get_dummies(df)
    return df_dummy


def matrix_target_getter(df_dummy):
    target = df_dummy['label']
    matrix = df_dummy.drop('label', axis=1)
    return matrix, target


def train_test_data_split(matrix,target):
    """

    :param df_dummy:
    :return:
    """

    X_train, X_test, target_train, target_test = train_test_split(matrix, target)
    print ("X_train: %s, X_test: %s" % (X_train.shape,X_test.shape))
    return X_train, X_test, target_train, target_test


def dummy_clf(X_train, target_train, strategy='stratified'):
    """
    a dummy classifer used as baseline to compare with other classifiers
    :param X_train:
    :param target_train:
    :param strategy:
    :return:
    """
    from sklearn.dummy import DummyClassifier
    model = DummyClassifier(strategy=strategy)
    model.fit(X_train,target_train)
    return model


def decision_tree_classifer(X_train, target_train):
    treeclf = tree.DecisionTreeClassifier(criterion='entropy', min_samples_split=5,random_state=0)
    treeclf = treeclf.fit(X_train, target_train)
    return treeclf


def random_forest_classifer(X_train, target_train):
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier(random_state=0)
    clf.fit(X_train,target_train)
    return clf


def svm_classifer(X_train,target_train):
    from sklearn import svm
    svm_clf = svm.SVC(kernel='linear',probability=True)
    svm_clf.fit(X_train,target_train)
    return svm_clf


def model_accuracy_getter(trained_model, X_test, target_test):
    return trained_model.score(X_test, target_test)


def cross_validate(classifer, matrix, target, cv=10):
    cv_scores = cross_val_score(classifer, matrix, target, cv=cv)
    return cv_scores.mean(), cv_scores.std()*2


def plot_confusion_mtx(trained_clf, X_test, y_test):
    from sklearn.metrics import plot_confusion_matrix
    import matplotlib.pyplot as plt
    ax = plot_confusion_matrix(trained_clf,X_test,y_test)
    plt.show()
    return ax.confusion_matrix


def classifi_report(trained_clf, X_test, target_test):
    from sklearn.metrics import classification_report
    predicted = trained_clf.predict(X_test)
    return classification_report(target_test,predicted)


def precision_recall_curv(trained_clf,X_test,y_test):
    from sklearn.metrics import precision_recall_curve,auc
    import matplotlib.pyplot as plt

    yhat = trained_clf.predict_proba(X_test)
    # retrieve just the probabilities for the positive class
    pos_probs = yhat[:, 1]
    # calculate model precision-recall curve
    precision, recall, _ = precision_recall_curve(y_test, pos_probs)
    print ("AUC: ", auc(recall, precision))
    # plot the model precision-recall curve
    plt.plot(recall, precision, marker='.')
    # axis labels
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.show()

if __name__=='__main__':
    from collections import Counter

    t_37C_240min_dict = ppp.load(open('tryp_37C_1440min_cleavage_label.p','rb'))
    # print (Counter([t_37C_240min_dict[each] for each in t_37C_240min_dict]))
    df_dummy = df_dummy_getter(t_37C_240min_dict)
    matrix, target = matrix_target_getter(df_dummy)

    X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
    svm_clf = svm_classifer(X_train,target_train)
    # score = cross_validate(svm_clf,matrix,target)
    print (plot_confusion_mtx(svm_clf,X_test,target_test))
    print(classifi_report(svm_clf,X_test, target_test))
    precision_recall_curv(svm_clf,X_test,target_test)
# two_d_list = []
# for polymer in t_37C_240min:
#     one_d_list = []
#     for aa in polymer:
#         one_d_list.append(aa)
#     one_d_list.append(t_37C_240min[polymer])
#     two_d_list.append(one_d_list)
#
# pd.set_option('display.max_columns', None)
#
# columns = [i for i in range(31)]
# columns.append('label')
# df = pd.DataFrame(two_d_list, columns=columns)
#
# df_dummy = pd.get_dummies(df)
#
# target = df_dummy['label']
# matrix = df_dummy.drop('label',axis=1)
# X_train, X_test, target_train, target_test = train_test_split(matrix,target)
# print (X_train.shape, X_test.shape)
#
# treeclf = tree.DecisionTreeClassifier(criterion='entropy', min_samples_split=5)
# treeclf = treeclf.fit(X_train,target_train)
# treepreds_test = treeclf.predict(X_test)
# print (treeclf.score(X_test,target_test))
#
# cv_scores = cross_val_score(treeclf,matrix,target,cv=10)
# print (cv_scores)
