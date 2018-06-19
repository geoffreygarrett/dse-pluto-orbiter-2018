# GENERAL & PROJECT-MADE
import datetime
import math
from trajectory_tool.helper import *
from trajectory_tool.plotting import *
from copy import deepcopy, copy
import matplotlib
from scipy import optimize

# ASTROPHYSICS & ORBITAL MECHANICS
import astropy.units as u
from astropy import time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from poliastro import iod
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
from poliastro.plotting import OrbitPlotter, OrbitPlotter3D
from poliastro.twobody import Orbit
from poliastro.util import time_range
import matplotlib.pyplot as plt


e0 = datetime.datetime.today()

solar_system_ephemeris.set("jpl")


def t_wait(r1, r2):
    # print(np.pi/np.sqrt(Sun.k) * ((r1*u.km + r2*u.km)/2) ** (3/2))
    return (np.pi/np.sqrt(Sun.k) * ((r1*u.km + r2*u.km)/2) ** (3/2)).to(u.s)

def delta_v(r1,r2):
    return np.sqrt(Sun.k.to(u.km**3/u.s**2).value/r1) * (np.sqrt((2*r2/(r1+r2)))-1) + np.sqrt(Sun.k.to(u.km**3/u.s**2).value/r2) * (1-np.sqrt(2*r1/(r1+r2)))

t_wait_list = []
dv_list = []
date_list = [e0]

r1 = np.linalg.norm(
    Orbit.from_body_ephem(Earth, time.Time(e0, scale='tdb')).r.to(u.km))
r2 = np.linalg.norm(
    Orbit.from_body_ephem(Pluto, time.Time(e0, scale='tdb')).r.to(u.km))

count=0
while True:

    t_w = 366 * 3600 * 24
    r1 = np.linalg.norm(Orbit.from_body_ephem(Earth, time.Time(date_list[-1], scale='tdb') + time.TimeDelta(t_w/(3600*24)*u.day)).r.to(u.km))
    r2 = np.linalg.norm(Orbit.from_body_ephem(Pluto, time.Time(date_list[-1], scale='tdb') +  time.TimeDelta(t_w/(3600*24)*u.day)).r.to(u.km))
    dv = delta_v(r1,r2)

    tof = t_wait(r1, r2).value
    date = date_list[-1] + datetime.timedelta(days=t_w/(3600*24))
    date_list.append(date)
    dv_list.append(dv)
    t_wait_list.append(tof/(3600)/24/365)
    if count>= 100:
        break
    count+=1

    # t_wait_list.append(t_wait(r1, r2).value)
    # dv_list.append(delta_v(r1,r2))
plt.grid(b=False, which='major', color='0.7', linestyle='--')
font = {'family' : 'normal',
        'size'   : 20,
        'serif': 'Times New Roman'}

matplotlib.rc('font', **font)

t_wait_ = np.array(t_wait_list)
# plt.style.use("seaborn")
plt.subplot(121)
plt.grid(b=False, which='major', color='0.7', linestyle='--')

plt.xlabel('Launch date [year]')
plt.ylabel('Time of flight [year]')
plt.plot(date_list[1:], t_wait_)
plt.subplot(122)
# plt.tight_layout()
plt.xlabel('Launch date [year]')
plt.grid(b=False, which='major', color='0.7', linestyle='--')

plt.ylabel('$\Delta{V}$ [km/s]')
plt.plot(date_list[1:], dv_list)
plt.rcParams["font.family"] = "serif"
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.tight_layout()
plt.show()


