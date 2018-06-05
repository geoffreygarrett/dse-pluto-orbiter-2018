def gravity_assist2(self, body0, body1, body2, epoch0, epoch1, epoch2, plot=False):
    ss0 = Orbit.from_body_ephem(body0, epoch0)
    ss1 = Orbit.from_body_ephem(body1, epoch1)
    ss2 = Orbit.from_body_ephem(body2, epoch2)

    tof1 = epoch1 - epoch0
    tof2 = epoch2 - epoch1

    (v0_1, v1_1), = iod.lambert(Sun.k, ss0.r, ss1.r, tof1)
    (v0_2, v1_2), = iod.lambert(Sun.k, ss1.r, ss2.r, tof2)

    b0_str = body0.__str__().split(' ')[0].lower()
    b1_str = body1.__str__().split(' ')[0].lower()
    b2_str = body2.__str__().split(' ')[0].lower()

    d_limits = [i for i in [body_d_domain[b1_str]['lower']] + [body_d_domain[b1_str]['upper']]]

    tv_leg1 = time_range(start=epoch0, end=epoch1, periods=100)
    tv_leg2 = time_range(start=epoch1, end=epoch2, periods=100)

    # v_out, delta = grav_ass(v_s, v_p, curr_lambert_soln.v0, 10e2, mu)

    v_out1, _ = grav_ass(v1_1, ss1.state.v, v0_2, d_limits[0] + (d_limits[1] - d_limits[0]) * 0,
                         body1.k.to(u.km ** 3 / u.s ** 2))
    v_out2, _ = grav_ass(v1_1, ss1.state.v, v0_2, d_limits[1], body1.k.to(u.km ** 3 / u.s ** 2))

    # def py_ang(v1, v2):
    #     """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    #     cosang = np.dot(v1, v2)
    #     sinang = np.linalg.norm(np.cross(v1, v2))
    #     return np.arctan2(sinang, cosang)
    #
    # def f(x):
    #     v_out,_ = grav_ass(v1_1, ss1.state.v, v0_2, x, body1.k.to(u.km**3/u.s**2))
    #     print(v_out)
    #     print(v0_2)
    #     print(py_ang(v_out, v0_2))
    #     return py_ang(v_out, v0_2)
    #
    # res = minimize_scalar(f, method='bounded', bounds=d_limits)
    # print(d_limits)
    # print(res)

    x = np.linspace(start=d_limits[0], stop=d_limits[1], num=300)
    y = [grav_ass(v1_1, ss1.state.v, v0_2, x_i, body1.k.to(u.km ** 3 / u.s ** 2)) for x_i in x]
    y, _ = zip(*y)
    # print(y)
    #
    # print(y)
    # print(v1_1)
    # print(v0_2)
    y_new1 = [v0_2 - y_i for y_i in list(y)]
    y_new = [np.linalg.norm(v0_2 - y_i) for y_i in list(y)]
    # # y = np.linalg.norm( y-v0_2 )
    #
    plt.plot(x, y_new1)
    plt.show()
    plt.clf()

    print(v_out1)
    print(v0_2)

    ss = OrbitPlotter()
    ss.set_attractor(Sun)

    # def_optimal_v = res.x*(u.km/u.s)

    ss.plot(ss0)
    ss.plot(ss1)
    ss.plot(ss2)

    trajec_in = Orbit.from_vectors(Sun, ss0.r, v0_1, epoch0)
    trajec_out_needed = Orbit.from_vectors(Sun, ss1.r, v0_2, epoch1)

    trajec_out_low_b = Orbit.from_vectors(Sun, ss1.r, v_out1, epoch1)
    trajec_out_high_b = Orbit.from_vectors(Sun, ss1.r, v_out2, epoch1)

    # print(def_optimal_v)
    # trajec_opt = Orbit.from_vectors(Sun, ss1.r, def_optimal_v, epoch1)

    ss.plot_trajectory(trajec_in.sample(tv_leg1)[-1], label='Leg 1', color='blue')
    ss.plot_trajectory(trajec_out_needed.sample(tv_leg2)[-1], label='Leg 2', color='blue')

    ss.plot_trajectory(trajec_out_low_b.sample(tv_leg2)[-1], label='Leg 2', color='green')

    ss.plot_trajectory(trajec_out_high_b.sample(tv_leg2)[-1], label='Leg 2', color='red')

    # ss.plot_trajectory(trajec_opt.sample(tv_leg2)[-1], label = 'Leg 2 optimal', color='red')

    plt.show()

    # for i in range(len(_raw_itinerary['durations'])):
    #     frame.plot_trajectory(
    #         _itinerary_plot_data[i + 1]['tp'].sample(_itinerary_data[i + 1]['tv'])[-1],
    #         label="Leg {}".format(i + 1),
    #         color=color_trans)

    # _itinerary_data[i + 1]['tv'] = time_range(start=_itinerary_data[i]['d'],
    #                                           end=_itinerary_data[i + 1]['d'],
    #                                           periods=self.N)