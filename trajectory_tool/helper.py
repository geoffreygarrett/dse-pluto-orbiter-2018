import plotly
from collections import namedtuple

KEY = 'kHLfFnsPiyxyAfWXgLN6'
USER = 'Jones1311'
plotly.tools.set_credentials_file(username=USER, api_key=KEY)
lambert_parameters = namedtuple('lambert_parameters', 'r0 r1 v0 v1 tof attractor epoch0 epoch1 ss0 ss1 sst')
gravass_parameters = namedtuple('gravass_parameters',
                                'a_i_mag a_f_mag e_i_mag e_f_mag v_inf_i_vec v_inf_f_vec v_planet_i_vec v_planet_f_vec r_p_dv_mag '
                                'v_p_i_vec v_p_f_vec t_p_i t_p_f e_i_vec e_f_vec aop lan inc type r_entry r_exit epoch_entry epoch_exit body_ga epoch_rp r_p')


factor = {
    'mercury': 1.082,

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
    'mercury': {
        'upper': 1.12E5,
        'lower': factor['mercury'] * 2439.7
    },
    'venus': {
        'upper': 6.16E5,
        'lower': factor['venus'] * 6051.8
    },
    'earth': {
        'upper': 9.25E5,
        'lower': factor['earth'] * 6371.
    },
    'mars': {
        'upper': 5.76E5,
        'lower': factor['mars'] * 3390.
    },
    'jupiter': {
        'upper': 4.82E7,
        'lower': factor['jupiter'] * 69911.
    },
    'saturn': {
        'upper': 5.48E7,
        'lower': factor['saturn'] * 58232.
    },
    'uranus': {
        'upper': 5.18E7,
        'lower': factor['uranus'] * 25362.
    },
    'neptune': {
        'upper': 8.66E7,
        'lower': factor['neptune'] * 24622.
    },
    'pluto': {
        'upper': 3.16E6,
        'lower': factor['pluto'] * 1188.
    }
}
