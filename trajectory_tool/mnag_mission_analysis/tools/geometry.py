import math
import numpy as np
from copy import deepcopy


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def angle_between_cc(self, v1, v2):
    # TODO: Change the rotation check method from z to the normal vector.
    """ Returns the angle clockwise in z from v1 to v2."""
    a_i = self.angle_between(v1, v2)
    # rotation check
    v2_n = np.matmul(self.rotation_z(0.01), deepcopy(v2))
    a_new = self.angle_between(v1, v2_n)
    if a_new > a_i:  # Clockwise measurement
        angle = a_i
    else:
        angle = 2 * np.pi - a_i  # Adjust for clockwise measurement
    return angle
