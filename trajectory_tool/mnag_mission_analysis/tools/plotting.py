import plotly
import plotly.graph_objs as go
import numpy as np
from poliastro.util import time_range
import astropy.units as u
from poliastro.plotting import OrbitPlotter, OrbitPlotter3D
from poliastro.twobody import Orbit
# from trajectory_tool.mnag_mission_analysis.planetary_flyby import PlanetaryFlyby
from copy import copy
from astropy import time


def polytime_2_datetime(_time):
    temp = copy(_time)
    temp.format = 'datetime'
    return temp.value


def trace_point(position):
    return [go.Scatter3d(
        x=[position[0]],
        y=[position[1]],
        z=[position[2]],
        mode='marker'
    )]


def trace_body(radius):
    return [go.Surface(
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
                    )]


def trace_soi(r_soi_magnitude, equatorial_plane=True):
    data = []

    x = r_soi_magnitude * np.outer(np.cos(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100)))
    y = r_soi_magnitude * np.outer(np.sin(np.linspace(0, 2 * np.pi, 100)),
                               np.sin(np.linspace(0, np.pi, 100)))
    z = r_soi_magnitude * np.outer(np.ones(100), np.cos(np.linspace(0, np.pi, 100)))

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
                    )
    data.append(sphere)

    if equatorial_plane:
        tracez = go.Surface(z=z_offset,
                            x=x,
                            y=y,
                            colorscale=[[0, 'rgb(180, 40, 40)'],
                                        [1, 'rgb(180, 40, 40)']],
                            showlegend=False,
                            showscale=False,
                            text='testing',
                            hoverinfo='text',
                            opacity=0.4
                            )
        data.append(tracez)
    return data


def plot_planetary_flyby(planetary_flyby):

    op = OrbitPlotter3D()
    data = []
    body = planetary_flyby.planetary_node.body
    op.set_attractor(body)
    ss_i = Orbit.from_classical(attractor=body,
                                a=planetary_flyby._sma_i,
                                ecc=np.linalg.norm(planetary_flyby._ecc_i) * u.one,
                                inc=planetary_flyby._inc * u.rad,
                                raan=planetary_flyby._raan * u.rad,
                                argp=planetary_flyby._aop * u.rad,
                                nu=0 * u.rad,
                                epoch=time.Time(planetary_flyby._epoch_periapsis))

    ss_f = Orbit.from_classical(attractor=body,
                                a=planetary_flyby._sma_f,
                                ecc=np.linalg.norm(planetary_flyby._ecc_f) * u.one,
                                inc=planetary_flyby._inc * u.rad,
                                raan=planetary_flyby._raan * u.rad,
                                argp=planetary_flyby._aop * u.rad,
                                nu=0 * u.rad,
                                epoch=time.Time(planetary_flyby._epoch_periapsis))

    epoch_entry_dt = time.Time(planetary_flyby.planetary_node.epoch_entry)
    epoch_exit_dt = time.Time(planetary_flyby.planetary_node.epoch_exit)
    epoch_rp_dt = time.Time(planetary_flyby.planetary_node.epoch_periapsis)

    # Time frame arrays for entry and exit.
    tv_ent = time_range(epoch_entry_dt, periods=100, spacing=None, end=epoch_rp_dt)
    tv_ent_f = time_range(epoch_rp_dt, periods=100, spacing=None,
                          end=epoch_rp_dt + (epoch_rp_dt - epoch_entry_dt))
    tv_ext = time_range(epoch_rp_dt, periods=100, spacing=None, end=epoch_exit_dt)
    op.plot_trajectory(ss_i.sample(tv_ent)[-1], label='Entry hyperbole')
    op.plot_trajectory(ss_f.sample(tv_ext)[-1], label='Exit hyperbole')

    ss_i_after = ss_i.sample(tv_ent_f)[-1].get_xyz()

    ss_i_after_trace = go.Scatter3d(
        x=ss_i_after[0].value,
        y=ss_i_after[1].value,
        z=ss_i_after[2].value,
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

    data.append(ss_i_after_trace)
    # data += trace_body(10 * body.R.to(u.km).value)
    # data += trace_soi(planetary_flyby.planetary_node.soi_periapsis_magnitude)
    data += trace_point(planetary_flyby.planetary_node.soi_entry_position_body_ecliptic)
    data += trace_point(planetary_flyby.planetary_node.soi_exit_position_body_ecliptic)

    print(data)

    layout = go.Layout(title="<b>{}</b> Flyby<br>Periapsis Epoch: {}".format(body, planetary_flyby.planetary_node.epoch_periapsis), width=800, height=800)
    fig = go.Figure(data=data, layout=layout)
    plotly.plotly.iplot(fig)

