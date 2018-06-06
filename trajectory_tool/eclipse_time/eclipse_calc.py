"""
Created by Alejandro Daniel Noel
"""
import numpy as np
from ale_solar_system.solar_system import SolarSystem


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def angle_between(vector1, vector2):
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def get_eclipse_time(ap, i, e, rp, epoch):
    ss = SolarSystem()

    vec1es = ss.coordinates_of('earth', epoch)  # Earth position wrt Sun at epoch
    vec2es = ss.coordinates_of('earth', epoch + 30)  # Earth position wrt Sun 30 days after
    n_es = np.cross(vec1es, vec2es)  # Normal vector to the Earth-Sun orbital plane (direction according to right-hand rule)

    vec1cp = ss.coordinates_of('charon', epoch, ref_body='pluto')  # Charon position wrt Pluto at epoch
    vec2cp = ss.coordinates_of('charon', epoch + 1, ref_body='pluto')  # Charon position wrt Pluto 1 day after
    n_pc = np.cross(vec1cp, vec2cp)  # Normal vector to the Charon-Pluto orbital plane (direction according to right-hand rule)

    vec_ps = ss.coordinates_of('pluto', epoch)

    a_pc_es = angle_between(n_pc, n_es)     # note: order matters!
    a_ps_pc = angle_between(vec_ps, n_pc)   # note: order matters!

    



if __name__ == '__main__':
    get_eclipse_time(0, 0, 0, 0, 2458275.5)
