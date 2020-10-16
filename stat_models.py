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


def decision_tree_classifer(X_train, target_train):
    treeclf = tree.DecisionTreeClassifier(criterion='entropy', min_samples_split=5)
    treeclf = treeclf.fit(X_train, target_train)
    return treeclf


def random_forest_classifer(X_train, target_train):
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier()
    clf.fit(X_train,target_train)
    return clf


def svm_classifer(X_train,target_train):
    from sklearn import svm
    svm_clf = svm.SVC()
    svm_clf.fit(X_train,target_train)
    return svm_clf


def model_accuracy_getter(trained_model, X_test, target_test):
    return trained_model.score(X_test, target_test)


def cross_validate(classifer, matrix, target, cv=10):
    cv_scores = cross_val_score(classifer, matrix, target, cv=cv)
    return cv_scores.mean(), cv_scores.std()*2


if __name__=='__main__':

    t_37C_240min_dict = ppp.load(open('tryps_37C_240min_cleavage_label.p','rb'))
    df_dummy = df_dummy_getter(t_37C_240min_dict)
    matrix, target = matrix_target_getter(df_dummy)
    X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
    svm_clf = svm_classifer(X_train,target_train)
    score = cross_validate(svm_clf,matrix,target)
    print (score)

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
