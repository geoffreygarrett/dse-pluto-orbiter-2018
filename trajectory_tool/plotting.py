import plotly
import plotly.graph_objs as go
import numpy as np


def create_soi(rsoi):
    return go.Surface(
                                        x=rsoi * np.outer(np.cos(np.linspace(0, 2 * np.pi, 100)),
                                                               np.sin(np.linspace(0, np.pi, 100))),
                                        y=rsoi * np.outer(np.sin(np.linspace(0, 2 * np.pi, 100)),
                                                               np.sin(np.linspace(0, np.pi, 100))),
                                        z=rsoi * np.outer(np.ones(100), np.cos(np.linspace(0, np.pi, 100))),
                                        colorscale=[[0, 'rgb(180, 40, 40)'],
                                                    [1, 'rgb(180, 40, 40)']],
                                        showscale=False,
                                        opacity=0.2
                    )