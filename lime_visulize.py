import lime.lime_tabular
import pickle as ppp
from stat_models import matrix_target, train_test_data_split
from parameters import custom_ohe, array_to_seq


t_37C_240min_dict = ppp.load(open('tryp_24h_label_dict_11_8.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = custom_ohe(matrix)
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
svm_clf = ppp.load(open('randomf_tryp37c_1440_11_8.p','rb'))
print(array_to_seq(X_test[11]))
print (target_test[11])
print (svm_clf.predict_proba([X_test[11]]))
explainer = lime.lime_tabular.LimeTabularExplainer(X_train)
exp = explainer.explain_instance(X_test[11], svm_clf.predict_proba, num_features=6)
print (exp.as_list())