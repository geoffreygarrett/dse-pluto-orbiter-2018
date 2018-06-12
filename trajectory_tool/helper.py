



import plotly
from collections import namedtuple
KEY = 'kHLfFnsPiyxyAfWXgLN6'
USER = 'Jones1311'
plotly.tools.set_credentials_file(username=USER, api_key=KEY)
lambert_parameters = namedtuple('lambert_parameters', 'r0 r1 v0 v1 tof attractor epoch0 epoch1 ss0 ss1 sst')
gravass_parameters = namedtuple('gravass_parameters', 'a_i a_f e_i e_f v_inf_i v_inf_f v_planet_i v_planet_f r_p_dv '
                                                      'v_p_i v_p_f t_p_i t_p_f')


factor = {
    'mercury':1.082,

    'venus': 1.047,

    'earth': 1.048,

    'mars': 1.076,

    'jupiter': 1.6,

    'saturn': 1.342,

    'uranus': 4.19,

    'neptune': 1.181,

    'pluto': 9.34,

}

body_d_domain = {
          'mercury':{
                            'upper':  1.12E5,
                            'lower': factor['mercury'] * 2539.7
  },
      'venus':{
                            'upper': 6.16E5,
                            'lower': factor['venus'] * 2440.
  },
      'earth':{
                            'upper': 9.25E5,
                            'lower': factor['earth'] * 2440.
  },
      'mars':{
                            'upper': 5.76E5,
                            'lower': factor['mars'] * 2440.
  },
      'jupiter':{
                            'upper': 4.82E7,
                            'lower': factor['jupiter'] * 2440.
  },
      'saturn':{
                            'upper': 5.48E7,
                            'lower': factor['saturn'] * 2440.
  },
      'uranus':{
                            'upper': 5.18E7,
                            'lower': factor['uranus'] * 2440.
  },
      'neptune':{
                            'upper': 8.66E7,
                            'lower': factor['neptune'] * 2440.
  },
      'pluto':{
                            'upper': 3.16E6,
                            'lower': factor['pluto'] * 2440.
  }
}