"""Export static & interactive plots"""

import os
import json
import plotly
import plotly.io as pio
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt

def get_plotly_template(name: str = "plotly_dark"):
    """
    Returns a Plotly template object by name.
    
    Args:
        name (str): Name of the Plotly template (e.g., 'plotly_dark', 'ggplot2').
    
    Returns:
        dict: The template object.
    
    Raises:
        ValueError: If the template name is not found.
    """
    if name not in pio.templates:
        raise ValueError(f"Template '{name}' not found. Available: {list(pio.templates.keys())}")
    return pio.templates[name]

def export_plotly_figure(fig: go.Figure, filename: str, template: str = "plotly_dark"):
    """
    Exports a Plotly figure to an HTML file with the specified template.
    
    Args:
        fig (go.Figure): The Plotly figure to export.
        filename (str): The name of the output HTML file.
        template (str): The name of the Plotly template to apply.
    
    Raises:
        ValueError: If the template name is not found.
    """
    try:
        fig.update_layout(template=get_plotly_template(template))
        fig.write_html(filename)
        print(f"Plotly figure exported successfully to {filename}")
    except Exception as e:
        print(f"Error exporting Plotly figure: {e}")

def export_matplotlib_figure(fig: plt.Figure, filename: str):
    """
    Exports a Matplotlib figure to a PNG file.
    
    Args:
        fig (plt.Figure): The Matplotlib figure to export.
        filename (str): The name of the output PNG file.
    """
    try:
        fig.savefig(filename, bbox_inches='tight')
        print(f"Matplotlib figure exported successfully to {filename}")
    except Exception as e:
        print(f"Error exporting Matplotlib figure: {e}")

def export_plotly_figure_json(fig: go.Figure, filename: str):
    """
    Exports a Plotly figure to a JSON file.
    
    Args:
        fig (go.Figure): The Plotly figure to export.
        filename (str): The name of the output JSON file.
    """
    try:
        fig_json = fig.to_json()
        with open(filename, 'w') as f:
            json.dump(json.loads(fig_json), f, indent=4)
        print(f"Plotly figure exported successfully to {filename}")
    except Exception as e:
        print(f"Error exporting Plotly figure to JSON: {e}")

def export_plotly_figure_image(fig: go.Figure, filename: str, format: str = 'png'):
    """
    Exports a Plotly figure to an image file (PNG, JPEG, etc.).
    
    Args:
        fig (go.Figure): The Plotly figure to export.
        filename (str): The name of the output image file.
        format (str): The image format (e.g., 'png', 'jpeg').
    
    Raises:
        ValueError: If the format is not supported.
    """
    supported_formats = ['png', 'jpeg', 'jpg', 'svg', 'pdf']
    if format not in supported_formats:
        raise ValueError(f"Format '{format}' not supported. Supported formats: {supported_formats}")
    
    try:
        fig.write_image(filename, format=format)
        print(f"Plotly figure exported successfully to {filename}")
    except Exception as e:
        print(f"Error exporting Plotly figure to image: {e}")




