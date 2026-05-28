"""Export visualization plots in multiple formats."""

from typing import Optional, Literal
import os
import json
import plotly.io as pio
import plotly.graph_objects as go
import matplotlib.pyplot as plt

def get_plotly_template(name: str = "plotly_dark") -> dict:
    """Get a Plotly template by name.
    
    Parameters
    ----------
    name : str, optional
        Name of the template (e.g., 'plotly_dark', 'ggplot2').
        
    Returns
    -------
    dict
        The template object.
        
    Raises
    ------
    ValueError
        If template name is not found.
    """
    if name not in pio.templates:
        available = list(pio.templates.keys())
        raise ValueError(f"Template '{name}' not found. Available: {available}")
    return pio.templates[name]


def export_plotly_figure(
    fig: go.Figure,
    filename: str,
    template: str = "plotly_dark"
) -> None:
    """Export a Plotly figure to HTML.
    
    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to export.
    filename : str
        Output HTML filename.
    template : str, optional
        Plotly template name (default: 'plotly_dark').
    """
    try:
        os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
        fig.update_layout(template=get_plotly_template(template))
        fig.write_html(filename)
        print(f"Plotly figure exported: {filename}")
    except Exception as e:
        print(f"Error exporting Plotly figure: {e}")


def export_matplotlib_figure(
    fig: plt.Figure,
    filename: str,
    dpi: int = 300
) -> None:
    """Export a Matplotlib figure to PNG.
    
    Parameters
    ----------
    fig : plt.Figure
        The Matplotlib figure to export.
    filename : str
        Output PNG filename.
    dpi : int, optional
        Resolution in DPI (default: 300).
    """
    try:
        os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Matplotlib figure exported: {filename}")
    except Exception as e:
        print(f"Error exporting Matplotlib figure: {e}")


def export_plotly_figure_json(
    fig: go.Figure,
    filename: str
) -> None:
    """Export a Plotly figure to JSON.
    
    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to export.
    filename : str
        Output JSON filename.
    """
    try:
        os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
        fig_json = fig.to_json()
        with open(filename, 'w') as f:
            json.dump(json.loads(fig_json), f, indent=2)
        print(f"Plotly figure exported to JSON: {filename}")
    except Exception as e:
        print(f"Error exporting to JSON: {e}")


def export_plotly_figure_image(
    fig: go.Figure,
    filename: str,
    format: Literal['png', 'jpeg', 'jpg', 'svg', 'pdf'] = 'png'
) -> None:
    """Export a Plotly figure to image format.
    
    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to export.
    filename : str
        Output image filename.
    format : str, optional
        Image format: 'png', 'jpeg', 'jpg', 'svg', or 'pdf' (default: 'png').
        
    Raises
    ------
    ValueError
        If format is not supported.
    """
    supported = ['png', 'jpeg', 'jpg', 'svg', 'pdf']
    if format not in supported:
        raise ValueError(f"Format '{format}' not supported. Use: {supported}")
    
    try:
        os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
        fig.write_image(filename, format=format)
        print(f"Plotly figure exported: {filename}")
    except Exception as e:
        print(f"Error exporting image: {e}")




