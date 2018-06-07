"""
Created by Alejandro Daniel Noel
"""
import copy

import numpy as np
from ale_solar_system.solar_system import SolarSystem

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def angle_between(vector1, vector2):
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def true_anomaly(mu, a, e, t):
    m = np.sqrt(mu / (a ** 3)) * t
    E = copy.deepcopy(m)
    for _ in range(20):  # Do 20 fixed-point iterations (convergence guaranteed for e ∈ [0, 1) )
        E = m + e * np.sin(E)
    return E


def rotation_mat_between_vectors(vector1, vector2):
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    v_n = np.cross(v1_u, v2_u)
    c = np.dot(v1_u, v2_u)
    sk_v = np.array([[0, -v_n[2], v_n[1]],
                     [v_n[2], 0, -v_n[0]],
                     [-v_n[1], v_n[0], 0]])
    return np.eye(3) + sk_v + np.linalg.matrix_power(sk_v, 2) / (1 + c)


def orbit_coords_2d(mu, e, rp, t):
    a = rp / (1 - e)
    E = true_anomaly(mu, a, e, t)

    # p and q are the coordinates in the plane of the orbit with p+ pointing towards periapsis
    p = a * (np.cos(E) - e)
    q = a * np.sin(E) * np.sqrt(1 - e * e)
    return p, q


def get_eclipse_time(mu, lan, i, e, rp, epoch):
    """
    Returns duration of earth occultation during one orbit
    :param lan: longitude of ascending node of the orbiter
    :param i: inclination of the orbiter
    :param e: eccentricity of the orbiter
    :param rp: radius at periapsis
    :param epoch: epoch at which to do the analysis (assumed constant throughout 1 orbit)
    :return: tuple (orbital time, occultation time)
    """
    ss = SolarSystem()

    vec1es = ss.coordinates_of('earth', epoch)  # Earth position wrt Sun at epoch
    vec2es = ss.coordinates_of('earth', epoch + 30)  # Earth position wrt Sun 30 days after
    n_es = np.cross(vec1es, vec2es)  # Normal vector to the Earth-Sun orbital plane (direction according to right-hand rule)

    vec1cp = ss.coordinates_of('charon', epoch, ref_body='pluto')  # Charon position wrt Pluto at epoch
    vec2cp = ss.coordinates_of('charon', epoch + 1, ref_body='pluto')  # Charon position wrt Pluto 1 day after
    n_pc = np.cross(vec1cp, vec2cp)  # Normal vector to the Charon-Pluto orbital plane (direction according to right-hand rule)

    vec_ps = ss.coordinates_of('pluto', epoch)  # Pluto position wrt Sun at epoch
    vec_pe = ss.coordinates_of('pluto', epoch, ref_body='earth')  # Pluto position wrt Earth at epoch

    a_pc_es = angle_between(n_pc, n_es)  # note: order matters!
    a_ps_pc = angle_between(vec_ps, n_pc)  # note: order matters!

    a = rp / (1 - e)
    orbital_period = 2 * np.pi * np.sqrt(a ** 3 / mu)

    px = []
    py = []
    pz = []
    map_ocultation = []
    for t in np.linspace(0, orbital_period, 1000):
        x, y = orbit_coords_2d(mu, e, rp, t)

        c_o, s_o = np.cos(lan), np.sin(lan)
        c_i, s_i = np.cos(i), np.sin(i)
        rot_mat = np.array([[c_o, -s_o * c_i, s_o * s_i],
                            [s_o, c_o * c_i, -c_o * s_i],
                            [0, s_i, c_i]])

        vec_3d_p = np.matmul(rot_mat, (x, y, 0))    # apply inclination and longitude of ascending node

        # Rotate the pluto orbit ref system so that the x-axis points towards the sun
        v_sc = np.matmul(rotation_mat_between_vectors((0.0, 0.0, 1.0), n_pc), vec_3d_p)
        rad_dist_ocul = np.linalg.norm(np.cross(v_sc, vec_pe)) / np.linalg.norm(vec_pe)
        px.append(v_sc[0])
        py.append(v_sc[1])
        pz.append(v_sc[2])

    cx = []
    cy = []
    cz = []
    for t in np.linspace(0.0, 7.0, 1000):
        vcs = ss.coordinates_of('charon', epoch + t, ref_body='pluto') * 1000
        cx.append(vcs[0])
        cy.append(vcs[1])
        cz.append(vcs[2])

    minx = min(px + cx)
    miny = min(py + cy)
    minz = min(pz + cz)
    maxx = max(px + cx)
    maxy = max(py + cy)
    maxz = max(pz + cz)
    midx = (minx + maxx) * 0.5
    midy = (miny + maxy) * 0.5
    midz = (minz + maxz) * 0.5
    mid_range = max(maxx - minx, maxy - miny, maxz - minz) * 0.5

    vec_pe = unit_vector(vec_pe) * 10000000.0

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(px, py, pz, label='orbiter orbit')
    ax.plot(cx, cy, cz, label='Charon orbit')
    ax.plot((0, vec_pe[0]), (0, vec_pe[1]), (0, vec_pe[2]), label="Earth-Pluto vector")
    ax.scatter(*[0, 0, 0], label="Pluto", color='k')

    ax.set_xlim((midx - mid_range, midx + mid_range))
    ax.set_ylim((midy - mid_range, midy + mid_range))
    ax.set_zlim((midz - mid_range, midz + mid_range))
    ax.legend()
    plt.show()


if __name__ == '__main__':
    print(get_eclipse_time(mu=8.7e11,
                           lan=np.radians(90.0),
                           i=np.radians(30.0),
                           e=0.2,
                           rp=1588000.0,
                           epoch=2458275.5))
