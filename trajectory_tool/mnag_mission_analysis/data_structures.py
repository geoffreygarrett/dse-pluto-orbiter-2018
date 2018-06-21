from collections import namedtuple


flyby_basic = namedtuple('flyby_basic', 'v_i, v_inf_i, v_planet_i, ecc_i, sma_i, v_p_i, r_p, r_p_dv,'
                                        'v_f, v_inf_f, v_planet_f, ecc_f, sma_f, v_p_f, a_req')

flyby_guess = namedtuple('flyby_guess', 'v_i, v_inf_i, v_planet_i, ecc_i, sma_i, v_p_i, r_p, r_p_dv, error, a_req,'
                                        'v_f, v_inf_f, v_planet_f, ecc_f, sma_f, v_p_f, argp, inc, raan, t_i, t_f, r_entry, r_exit')

flyby_refined = namedtuple('flyby_refined', 'v_i, v_inf_i, v_planet_i, ecc_i, sma_i, v_p_i, r_p, r_p_dv, a_req'
                                            'v_f, v_inf_f, v_planet_f, ecc_f, sma_f, v_p_f, argp, inc, raan, t_i, t_f,'
                                            'error_v, error_p, r_entry, r_exit')
