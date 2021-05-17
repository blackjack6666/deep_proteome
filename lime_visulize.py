import lime.lime_tabular
import pickle as ppp
from stat_models import matrix_target, train_test_data_split
from parameters import custom_ohe, array_to_seq
from sklearn.tree import export_graphviz, plot_tree
import matplotlib.pyplot as plt

t_37C_240min_dict = ppp.load(open('D:/data/deep_proteome/non_specfic_search/tryps_4h_polymer.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = custom_ohe(matrix)

test_matrix = [['Y','R','K','H','V','A','E','R','A','A','E','G','I','A','P','K','P','L','D','A','N','Q','M','A','A','L','V','E','L','L','K']]
test_matrix_ohe = custom_ohe(test_matrix)
print (array_to_seq(test_matrix_ohe[0]))
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
svm_clf = ppp.load(open('D:/data/deep_proteome/non_specfic_search/random_forest_tryps_4h.p','rb'))
print(array_to_seq(X_test[15]))
print (target_test[15])
print (svm_clf.predict_proba([X_test[15]]))
print (svm_clf.predict([X_test[15]]))

estimator = svm_clf.estimators_[80]
# Export as dot file
# export_graphviz(estimator, out_file='tree.dot',
#
#                 class_names = ['miss cleaved', 'cleaved'],
#                 rounded = True, proportion = False,
#                 precision = 2, filled = True)

# plt.figure(figsize=(40,20))
# plot_tree(estimator,
#           max_depth=3,
#                 rounded = True, proportion = False,
#                  filled = True)
# plt.show()
explainer = lime.lime_tabular.LimeTabularExplainer(X_train)
exp = explainer.explain_instance(X_test[15], svm_clf.predict_proba, num_features=20)
print (exp.as_list())
