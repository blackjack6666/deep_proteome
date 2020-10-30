import plotly.graph_objects as go
# categories = ['Precision','Recall','F1-score',
#               'Accuracy', 'AUC']
#
# fig = go.Figure()
#
# fig.add_trace(go.Scatterpolar(
#       r=[0.85,0.85,0.85,0.84,0.92],
#       theta=categories,
#       fill='none',
#       name='SVM'
# ))
#
# fig.add_trace(go.Scatterpolar(
#       r=[0.87,0.88,0.87,0.86,0.94],
#       theta=categories,
#       fill='none',
#       name='Random forest'
# ))
#
# fig.add_trace(go.Scatterpolar(
#       r=[0.53,0.55,0.54,0.51,0.66],
#       theta=categories,
#       fill='none',
#       name='Dummy classifer'
# ))
#
# fig.add_trace(go.Scatterpolar(
#       r=[0.87,0.83,0.85,0.84,0.90],
#       theta=categories,
#       fill='none',
#       name='CNN'
# ))
#
# fig.add_trace(go.Scatterpolar(
#       r=[0.83,0.87,0.85,0.83,0.90],
#       theta=categories,
#       fill='none',
#       name='LSTM'
# ))
#
#
#
# fig.update_layout(
#   polar=dict(
#     radialaxis=dict(
#       visible=True,
#       range=[0, 1]
#     )),
#   showlegend=True,
#
# )

import matplotlib.pyplot as plt
import pandas as pd
from math import pi

# Set data
df = pd.DataFrame({
    'group': ['SVM', 'Random forest', 'Dummy classifier', 'CNN', 'LSTM', 'Sequential'],
    'Precision': [0.80, 0.9, 0.58, 0.80, 0.73, 0.67],
    'Recall': [0.81, 0.82, 0.53, 0.78, 0.65, 0.72],
    'F1-score': [0.80, 0.86, 0.55, 0.79, 0.69, 0.70],
    'Accuracy': [0.79, 0.85, 0.50, 0.78, 0.69, 0.67],
    'AUC': [0.88, 0.93, 0.66, 0.84, 0.78, 0.72]
})

# ------- PART 1: Create background

# number of variable
categories = list(df)[1:]
N = len(categories)

# What will be the angle of each axis in the plot? (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# Initialise the spider plot
ax = plt.subplot(111, polar=True)

# If you want the first axis to be on top:
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

# Draw one axe per variable + add labels labels yet
plt.xticks(angles[:-1], categories)

# Draw ylabels
ax.set_rlabel_position(0)
plt.yticks([0.4,0.6,0.8,1.0], ["0.4", "0.6", "0.8", "1.0"], color="grey", size=7)
plt.ylim(0, 1)

# ------- PART 2: Add plots

# Plot each individual = each line of the data
# I don't do a loop, because plotting more than 3 groups makes the chart unreadable

color = ['r','b','y','p','g','a']
label_list = ['SVM', 'Random forest', 'Dummy classifier', 'CNN', 'LSTM', 'Sequential']
for i,c,label in zip(range(6),color,label_list):
    values = df.loc[i].drop('group').values.flatten().tolist()
    print (values)
    values += values[:1]
    ax.plot(angles, values, linewidth=1, linestyle='solid', label=label)
    # ax.fill(angles, values, c, alpha=0.1)

# Ind1
# values = df.loc[0].drop('group').values.flatten().tolist()
# print (values)
# values += values[:1]
# ax.plot(angles, values, linewidth=1, linestyle='solid', label="group A")
# ax.fill(angles, values, 'b', alpha=0.1)
#
# # Ind2
# values = df.loc[1].drop('group').values.flatten().tolist()
# values += values[:1]
# ax.plot(angles, values, linewidth=1, linestyle='solid', label="group B")
# ax.fill(angles, values, 'r', alpha=0.1)



# Add legend
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
plt.show()