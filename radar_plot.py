import plotly.graph_objects as go
categories = ['Precision','Recall','F1-score',
              'Accuracy', 'AUC']

fig = go.Figure()

fig.add_trace(go.Scatterpolar(
      r=[0.85,0.85,0.85,0.84,0.92],
      theta=categories,
      fill='none',
      name='SVM'
))

fig.add_trace(go.Scatterpolar(
      r=[0.87,0.88,0.87,0.86,0.94],
      theta=categories,
      fill='none',
      name='Random forest'
))

fig.add_trace(go.Scatterpolar(
      r=[0.53,0.55,0.54,0.51,0.66],
      theta=categories,
      fill='none',
      name='Dummy classifer'
))

fig.add_trace(go.Scatterpolar(
      r=[0.87,0.83,0.85,0.84,0.90],
      theta=categories,
      fill='none',
      name='CNN'
))

fig.add_trace(go.Scatterpolar(
      r=[0.83,0.87,0.85,0.83,0.90],
      theta=categories,
      fill='none',
      name='LSTM'
))



fig.update_layout(
  polar=dict(
    radialaxis=dict(
      visible=True,
      range=[0, 1]
    )),
  showlegend=True,

)

fig.show()