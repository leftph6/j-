import plotly.graph_objects as go
import json

# Data from the provided JSON
data = {"eigenvalues": [101.9956, 99.9934, 98.9864], "labels": ["λ₁", "λ₂", "λ₃"]}

eigenvalues = data["eigenvalues"]
labels = data["labels"]

# Create the bar chart
fig = go.Figure()

# Add bar trace with blue color
fig.add_trace(go.Bar(
    x=labels,
    y=eigenvalues,
    marker_color='#1FB8CD',  # Using the first brand color (blue/cyan)
    text=[f'{val:.4f}' for val in eigenvalues],  # Show exact values on bars
    textposition='outside',  # Position text above bars
    textfont=dict(size=12),
    name='Eigenvalues'
))

# Update layout
fig.update_layout(
    title='Top 3 Eigenvalues Comparison',
    xaxis_title='Eigenvalue No.',
    yaxis_title='Value',
    showlegend=False,  # Hide legend since we only have one series
)

# Update traces to prevent clipping
fig.update_traces(cliponaxis=False)

# Update y-axis to give some space above the bars for text
fig.update_yaxes(range=[95, max(eigenvalues) * 1.05])

# Save as both PNG and SVG
fig.write_image('eigenvalue_chart.png')
fig.write_image('eigenvalue_chart.svg', format='svg')

print("Chart created successfully!")
print(f"Eigenvalues: {eigenvalues}")
print(f"Labels: {labels}")