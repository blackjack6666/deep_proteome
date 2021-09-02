"""
SHAP implementation
complete code is in lab jupyter hub xshao/dL_model_interpret.ipynb
"""
import pandas as pd
import shap
from tensorflow import keras
import pickle as ppp
from stat_models import df_dummy_getter, matrix_target_getter, train_test_data_split,dump_data, matrix_target,ohe
from parameters import custom_ohe
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from tcn import TCN
def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)

tf.compat.v1.disable_v2_behavior()

models = ['tryp_4h_gluc_ON','trypsin_4h_15mer_tcn','trypsin_CT_15mer_tcn','ct_4h_15mer_tcn','tryp_30_thermo_30','thermolysine_PRIDE']
datasets = ['tryp_Gluc_ON_15mer.p','tryps_4h_15mer.p','trypsin_CT_15mer.p', 'ct_37_15mer.p', 'tryp_30_thermo_30.p', 'thermolysin_PRIDE.p']

"""
model_weight_dict = {}

for each_model_path, data_path in zip(models,datasets):
    model = keras.models.load_model('D:/data/deep_proteome/deep_learning_models/'+each_model_path,custom_objects={'TCN': TCN})
    data = ppp.load(open('D:/data/deep_proteome/non_specfic_search/'+data_path,'rb'))
    matrix, target = matrix_target(data)
    matrix = custom_ohe(matrix,polymer_len=15)
    X_train, X_test, target_train, target_test = train_test_data_split(matrix, target)

# model = keras.models.load_model('D:/data/deep_proteome/bi_directional_lstm_trypsin4h')
# t_37C_240min_dict = ppp.load(open('D:/data/deep_proteome/non_specfic_search/tryps_4h_polymer.p','rb'))
# matrix,target = matrix_target(t_37C_240min_dict)
# matrix = custom_ohe(matrix)
# X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
# # X_train, X_test= X_train.reshape(X_train.shape[0],31,22), \
# #                                 X_test.reshape(X_test.shape[0],31,22)
# print (model.predict(X_train[:300]))
# print (target_train[:300])

    #initialize js methods for visualization
    shap.initjs()
    background = X_train[np.random.choice(X_train.shape[0], 1000, replace=False)]
    # print (background.shape)
    # create an instance of the DeepSHAP which is called DeepExplainer
    explainer_shap = shap.DeepExplainer(model=model,
                                     data=background)


    # Fit the explainer on a subset of the data (you can try all but then gets slower)
    shap_values = explainer_shap.shap_values(X=X_test[:500,:], ranked_outputs=True)
    model_weight_dict[each_model_path] = np.sum(shap_values[0][0],axis=0)
    print (each_model_path)

ppp.dump(model_weight_dict,open('D:/data/deep_proteome/deep_learning_models/models_weight.p','wb'))
"""
model_weight_dict = ppp.load(open('D:/data/deep_proteome/deep_learning_models/models_weight.p','rb'))
print (model_weight_dict['tryp_4h_gluc_ON'][154:176])
combined_array = []
# for i, j in zip(model_weight_dict['trypsin_4h_15mer_tcn'], model_weight_dict['ct_4h_15mer_tcn']):
#     if i>0 and j>0 and i>=j:
#         combined_array.append(i)
#     elif i>0 and j>0 and i<j:
#         combined_array.append(j)
#     elif i<0 and j<0 and i<=j:
#         combined_array.append(i)
#     elif i<0 and j<0 and i>j:
#         combined_array.append(j)
#     else:
#         combined_array.append((i+i)/2)
#
# print(cosine_sim_calculating(model_weight_dict['trypsin_CT_15mer_tcn'],combined_array ))
# print (shap_values)

"""
# now let's inspect some individual explanations inferred by DeepSHAP

shap.force_plot(explainer_shap.expected_value,
                shap_values[0][0][3],matplotlib=True
                )
plt.show()
# shap.summary_plot(shap_values[3], X_test[:10], plot_type="bar")
"""