import lime.lime_tabular
import pickle as ppp
from stat_models import matrix_target, train_test_data_split
from parameters import custom_ohe, array_to_seq


t_37C_240min_dict = ppp.load(open('D:/data/non_specific_search/ecoli_non_specific_search_poly_dict.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = custom_ohe(matrix)

test_matrix = [['Y','R','K','H','V','A','E','R','A','A','E','G','I','A','P','K','P','L','D','A','N','Q','M','A','A','L','V','E','L','L','K']]
test_matrix_ohe = custom_ohe(test_matrix)
print (array_to_seq(test_matrix_ohe[0]))
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
svm_clf = ppp.load(open('randomf_ecoli_non_spec_search.p','rb'))
print(array_to_seq(X_test[11]))
print (target_test[11])
print (svm_clf.predict_proba([test_matrix_ohe[0]]))
print (svm_clf.predict([test_matrix_ohe[0]]))
explainer = lime.lime_tabular.LimeTabularExplainer(X_train)
exp = explainer.explain_instance(test_matrix_ohe[0], svm_clf.predict_proba, num_features=6)
print (exp.as_list())