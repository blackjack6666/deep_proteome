"""
different classification model implementation
"""

import pandas as pd
import numpy as np
import pickle as ppp
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn import tree
from sklearn.preprocessing import OneHotEncoder

def dump_data(polymer_label_dict):
    """
    dump 31mer and value into matrix before training, each
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
    matrix = np.array(two_d_list)
    matrix, target = [' '.join(each) for each in matrix[:,:-1]], matrix[:,-1].astype(np.int)

    return matrix,target


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


def matrix_target(polymer_label_dict):
    matrix = []
    for each in polymer_label_dict:
        one_d = []
        for each_aa in each:
            one_d.append(each_aa)
        one_d.append(polymer_label_dict[each])
        matrix.append(one_d)

    matrix = np.array(matrix)
    print (matrix.shape)
    matrix, target = matrix[:, :-1], matrix[:, -1].astype(np.int)
    return matrix, target


def ohe(matrix):
    """
    one hot encoder for 31-mers
    :param polymer_label_dict:
    :return:
    """

    # Create the encoder.
    encoder = OneHotEncoder(handle_unknown="ignore", sparse=False)
    encoder.fit(matrix)
    matrix = encoder.transform(matrix)
    return encoder, matrix

def matrix_target_getter(df_dummy):
    # from df_dummy
    target = df_dummy['label']
    matrix = df_dummy.drop('label', axis=1)
    return matrix, target


def train_test_data_split(matrix,target):
    """

    :param df_dummy:
    :return:
    """

    X_train, X_test, target_train, target_test = train_test_split(matrix, target,random_state=0)
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
    # print (yhat)
    # retrieve just the probabilities for the positive class
    pos_probs = yhat[:, 1]
    # calculate model precision-recall curve
    precision, recall, _ = precision_recall_curve(y_test, pos_probs)
    print ("AUC precision-recall: ", auc(recall, precision))
    # plot the model precision-recall curve
    plt.plot(recall, precision, marker='.')
    # axis labels
    # plt.ylim(0,1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.show()
    

def roc_curve(trained_clf,X_test,y_test):
    from sklearn.metrics import roc_curve,auc
    import matplotlib.pyplot as plt
    yhat = trained_clf.predict_proba(X_test)
    # print (yhat)
    # retrieve just the probabilities for the positive class
    pos_probs = yhat[:, 1]
    fpr, tpr, _ = roc_curve(y_test,pos_probs)
    roc_auc = auc(fpr,tpr)
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    plt.show()


if __name__=='__main__':
    from collections import Counter
    import time
    from parameters import custom_ohe

    t_37C_240min_dict = ppp.load(open('D:/data/non_specific_search/ecoli_non_specific_search_poly_dict.p','rb'))

    # test_dataset_dict = ppp.load(open('mouse_B_FT_31mer_dict.p','rb'))
    # predict_matrix = ppp.load(open('P62918_matrix_2d_array.p', 'rb'))

    print (Counter([t_37C_240min_dict[each] for each in t_37C_240min_dict]))
    # print (Counter([v for v in test_dataset_dict.values()]))
    pd.set_option('display.max_columns', 1000)
    df_dummy = df_dummy_getter(t_37C_240min_dict)
    # print (df_dummy.head())

    matrix,target = matrix_target(t_37C_240min_dict)
    matrix = custom_ohe(matrix)
    # encoder,matrix = ohe(matrix)
    # print (matrix.shape)
    # # test set from different dataset
    # test_maxtrix, test_target = matrix_target(test_dataset_dict)
    # print (test_maxtrix,test_target)
    # test_maxtrix = cust.transform(test_maxtrix)
    # print(test_maxtrix.shape)
    # predict_matrix = encoder.transform(predict_matrix)
    # matrix, target = matrix_target_getter(df_dummy)

    X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
    time_start = time.time()
    svm_clf = random_forest_classifer(X_train,target_train)
    ppp.dump(svm_clf, open('randomf_ecoli_non_spec_search.p','wb'))
    print('model trained time:',time.time() - time_start)
    # score = cross_validate(svm_clf,matrix,target)
    # print (score)
    print (plot_confusion_mtx(svm_clf,X_test,target_test))
    print(classifi_report(svm_clf,X_test, target_test))
    precision_recall_curv(svm_clf,X_test,target_test)
    roc_curve(svm_clf,X_test,target_test)
    # print (svm_clf.predict(predict_matrix))
    print (target_test[0])

