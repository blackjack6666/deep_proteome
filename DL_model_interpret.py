"""
-----
protease prediction, get feature weights matrix (SHAP values), get max and min in each column (feature) and use
max/min array as input for cosine sim comparison
-----
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

models = ['gluc_4h_tryps_ON','tryp_4h_gluc_ON','trypsin_4h_15mer_tcn','trypsin_CT_15mer_tcn','ct_4h_15mer_tcn','tryp_30_thermo_30','thermolysine_PRIDE','ct_24h_15mer_tcn']
datasets = ['gluc_4h_tryps_ON.p','tryp_Gluc_ON_15mer.p','tryps_4h_15mer.p','trypsin_CT_15mer.p', 'ct_37_15mer.p', 'tryp_30_thermo_30.p', 'thermolysin_PRIDE.p','ct_24h_15mer.p']

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
    background = X_train[np.random.choice(X_train.shape[0], 500, replace=False)]
    # print (background.shape)
    # create an instance of the DeepSHAP which is called DeepExplainer
    explainer_shap = shap.DeepExplainer(model=model,
                                     data=background)


    # Fit the explainer on a subset of the data (you can try all but then gets slower)
    shap_values = explainer_shap.shap_values(X=X_test[np.random.choice(X_test.shape[0], 500, replace=False),:], ranked_outputs=True)
    max_array,min_array = np.amax(shap_values[0][0],axis=0), np.amin(shap_values[0][0],axis=0)
    model_weight_dict[each_model_path]=(max_array,min_array)
    print (each_model_path)

# ppp.dump(model_weight_dict,open('D:/data/deep_proteome/deep_learning_models/models_weight_shap_max_min_dict.p','wb'))
"""
model_weight_dict = ppp.load(open('D:/data/deep_proteome/deep_learning_models/models_weight_shap_max_min_dict.p','rb'))
trypsin_ct_combined = np.mean([np.concatenate((model_weight_dict['ct_4h_15mer_tcn'][0],model_weight_dict['ct_4h_15mer_tcn'][1]),axis=0),
                               np.concatenate((model_weight_dict['trypsin_CT_15mer_tcn'][0],model_weight_dict['trypsin_CT_15mer_tcn'][1]),axis=0)],axis=0)

trypsin_thermolysin_combined = np.mean([np.concatenate((model_weight_dict['trypsin_4h_15mer_tcn'][0],model_weight_dict['trypsin_4h_15mer_tcn'][1]),axis=0),
                               np.concatenate((model_weight_dict['thermolysine_PRIDE'][0],model_weight_dict['thermolysine_PRIDE'][1]),axis=0)],axis=0)

fig,axs = plt.subplots(1,1)
axs.plot(trypsin_thermolysin_combined,c='black')
plt.title('trypsin_4h_thermolysin_PRIDE_combined_npmean_+-_impact')
plt.show()
# print (cosine_sim_calculating(np.concatenate((model_weight_dict['tryp_30_thermo_30'][0],model_weight_dict['tryp_30_thermo_30'][1]),axis=0),
#                               trypsin_thermolysin_combined))

"""
# now let's inspect some individual explanations inferred by DeepSHAP

shap.force_plot(explainer_shap.expected_value,
                shap_values[0][0][3],matplotlib=True
                )
plt.show()
# shap.summary_plot(shap_values[3], X_test[:10], plot_type="bar")
"""