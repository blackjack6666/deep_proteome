"""
SHAP implementation
complete code is in lab jupyter hub xshao/dL_model_interpret.ipynb
"""

import shap
from tensorflow import keras
import pickle as ppp
from stat_models import df_dummy_getter, matrix_target_getter, train_test_data_split,dump_data, matrix_target,ohe
from parameters import custom_ohe
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

tf.compat.v1.disable_v2_behavior()

model = keras.models.load_model('D:/data/deep_proteome/bi_directional_lstm_trypsin4h')

t_37C_240min_dict = ppp.load(open('D:/data/deep_proteome/non_specfic_search/tryps_4h_polymer.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = custom_ohe(matrix)
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
# X_train, X_test= X_train.reshape(X_train.shape[0],31,22), \
#                                 X_test.reshape(X_test.shape[0],31,22)
print (model.predict(X_train[:300]))
print (target_train[:300])

#initialize js methods for visualization
shap.initjs()
background = X_train[np.random.choice(X_train.shape[0], 100, replace=False)]
print (background.shape)
# create an instance of the DeepSHAP which is called DeepExplainer
explainer_shap = shap.DeepExplainer(model=model,
                                 data=background)

# Fit the explainer on a subset of the data (you can try all but then gets slower)
shap_values = explainer_shap.shap_values(X=X_train[:300], ranked_outputs=True)
print (shap_values)
# now let's inspect some individual explanations inferred by DeepSHAP
shap.force_plot(explainer_shap.expected_value,
                shap_values[0][0][3],matplotlib=True
                )
plt.show()
# shap.summary_plot(shap_values[3], X_test[:10], plot_type="bar")
