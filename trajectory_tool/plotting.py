import plotly
import plotly.graph_objs as go
import numpy as np
from poliastro.util import time_range
import astropy.units as u
from poliastro.plotting import OrbitPlotter, OrbitPlotter3D
from poliastro.twobody import Orbit
from copy import copy
from astropy import time

import matplotlib.pyplot as plt


def create_soi(rsoi):
    x = rsoi * np.outer(np.cos(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100)))
    y = rsoi * np.outer(np.sin(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100)))
    z = rsoi * np.outer(np.ones(100), np.cos(np.linspace(0, np.pi, 100)))

    proj_z = lambda x, y, z: z  # projection in the z-direction
    colorsurfz = proj_z(x, y, z)

    z_offset = 0 * np.ones(z.shape)  #
    x_offset = np.min(x) * np.ones(z.shape)
    y_offset = np.min(y) * np.ones(z.shape)

    sphere = go.Surface(
        x=x,
        y=y,
        z=z,
        colorscale=[[0, 'rgb(180, 40, 40)'],
                    [1, 'rgb(180, 40, 40)']],
        showscale=False,
        opacity=0.2,
        # showline=True
                    )

    tracez = go.Surface(z=z_offset,
                        x=x,
                        y=y,
                        colorscale=[[0, 'rgb(180, 40, 40)'],
                                    [1, 'rgb(180, 40, 40)']],
                        showlegend=False,
                        showscale=False,
                        # surfacecolor=colorsurfx,
                        text='testing',
                        hoverinfo='text',
                        opacity=0.4
                        )

    return [sphere, tracez]


def create_point(position):
    return [go.Scatter3d(
        x=[position[0].value],
        y=[position[1].value],
        z=[position[2].value],
        mode='marker'
    )]


def create_body(radius):
    return go.Surface(
        x=radius * np.outer(np.cos(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100))),
        y=radius * np.outer(np.sin(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100))),
        z=radius * np.outer(np.ones(100), np.cos(np.linspace(0, np.pi, 100))),
        colorscale=[[0, 'rgb(0, 0, 0)'],
                    [1, 'rgb(50, 50, 50)']],
        showscale=False,
        hoverinfo='Sphere of Influence (SOI)'
        # showline=True
                    )


def polytime_2_datetime(_time):
    temp = copy(_time)
    temp.format = 'datetime'
    return temp.value


def plot_pga_3D(_itinerary_data_indexed, rsoi, a_i, a_f, e_i, e_f, lan, aop, inclination, epoch_rp, r_entry, r_exit):
    op = OrbitPlotter3D()
    data = []
    body = _itinerary_data_indexed['b']
    op.set_attractor(body)
    ss_i = Orbit.from_classical(attractor=body, a=a_i, ecc=e_i * u.one, inc=inclination * u.rad,
                                raan=lan * u.rad, argp=aop * u.rad, nu=0 * u.rad, epoch=epoch_rp)

    ss_f = Orbit.from_classical(attractor=body, a=a_f, ecc=e_f * u.one, inc=inclination * u.rad,
                                raan=lan * u.rad, argp=aop * u.rad, nu=0 * u.rad, epoch=epoch_rp)

    epoch_entry_dt = polytime_2_datetime(_itinerary_data_indexed['d']['i'])
    epoch_exit_dt = polytime_2_datetime(_itinerary_data_indexed['d']['f'])
    epoch_rp_dt = polytime_2_datetime(_itinerary_data_indexed['d']['rp'])

    # Time frame arrays for entry and exit.
    tv_ent = time_range(epoch_entry_dt, periods=100, spacing=None, end=epoch_rp_dt)
    tv_ent_f = time_range(epoch_rp_dt, periods=100, spacing=None,
                          end=epoch_rp_dt + (epoch_rp_dt-epoch_entry_dt))
    tv_ext = time_range(epoch_rp_dt, periods=100, spacing=None, end=epoch_exit_dt)

    op.plot_trajectory(ss_i.sample(tv_ent)[-1], label='Entry hyperbole')
    op.plot_trajectory(ss_f.sample(tv_ext)[-1], label='Exit hyperbole')

    ss_i_after = ss_i.sample(tv_ent_f)[-1].get_xyz()

    ss_i_after_trace = go.Scatter3d(
        x = ss_i_after[0].value,
        y = ss_i_after[1].value,
        z = ss_i_after[2].value,
        mode='lines',
        line=dict(
            color='#1f77b4',
            width=7,
            dash="dash"
        ),
        text='Entry hyperbole extended',
        projection=dict(
            x=dict(
                show=True,
            ))
        )

    # trace1 = go.Surface(z=z,
    #                     x=x,
    #                     y=y,
    #                     colorscale=colorscale,
    #                     # text=textz,
    #                     hoverinfo='Equatorial Plane',
    #                     )

    # Low opacity visualisation of the SOI.
    # data.append()
    data.append(ss_i_after_trace)
    data.append(create_body(10*body.R.to(u.km).value))
    data += create_soi(rsoi)
    data += create_point(r_entry)
    data += create_point(r_exit)
    data = data + op._data

    layout = go.Layout(title="<b>{}</b> Flyby<br>Periapsis Epoch: {}".format(body, epoch_rp), width=800, height=800)
    fig = go.Figure(data=data, layout=layout)
    plotly.plotly.iplot(fig)

    # op.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)
    # op.savefig("{}_pga.png".format(str(body).split(' ')[0]), title="{} powered gravity assist".format(body))
