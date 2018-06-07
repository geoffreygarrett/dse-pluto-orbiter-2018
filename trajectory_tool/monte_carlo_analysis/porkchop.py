
import numpy as np
from scipy.interpolate import Rbf
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import json
from trajectory_tool.monte_carlo_analysis.obj_def import Trajectory
from trajectory_tool.ale_solar_system.time_utils import datetime_to_jd, jd_to_datetime


def plot_porkchop(X, Y, Z, title):
    x_grid, y_grid = np.meshgrid(np.linspace(min(X), max(X), 700), np.arange(min(Y), max(Y), 700))
    z_grid = np.zeros_like(x_grid)
    f = Rbf(X, Y, Z, epsilon=1)
    for i in range(x_grid.shape[0]):
        for j in range(x_grid.shape[1]):
            z_grid[i, j] = f(x_grid[i, j], y_grid[i, j])

    mask = z_grid < 0.0
    z_grid = np.ma.array(z_grid, mask=mask)

    corner_mask = False
    levels = [z_grid.min(), z_grid.min() + 1, z_grid.min() + 2, z_grid.min() + 4, z_grid.min() + 8, z_grid.min() + 16, z_grid.min() + 25]
    cs = plt.contourf(np.vectorize(jd_to_datetime)(x_grid), y_grid / 365.0, z_grid, levels, corner_mask=corner_mask)
    cs2 = plt.contour(cs, colors='k')

    plt.xlabel("Launch date [year]")
    plt.ylabel("Trajectory duration [years]")

    cbar = plt.colorbar(cs)
    cbar.ax.set_ylabel('delta V')
    cbar.add_lines(cs2)

    # Plot grid.
    plt.grid(c='k', ls='-', alpha=0.3)
    plt.suptitle(title)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":

    base_path = os.path.dirname(os.path.realpath(__file__))
    itineraries = [folder for folder in os.listdir(base_path) if os.path.isdir(folder) and 'earth' in folder]
    for itinerary in itineraries:
        results_dir = os.path.join(base_path, itinerary, 'case_data')
        if os.path.exists(results_dir):
            results_files = [os.path.join(results_dir, file) for file in os.listdir(results_dir) if 'cases' in str(file)]
            trajectories = []
            for results_file in results_files:
                # with open(results_file, 'r') as f:
                #     file_string = f.read()
                # if '}\n=======\n[\n' in file_string:
                #     the_index = file_string.index('}\n=======\n[\n')
                #     print(file_string[the_index - 3: the_index + 40])
                #     file_string = file_string.replace('}\n=======\n[\n', '},\n')
                #     print(file_string[the_index - 3: the_index + 40])
                # if '<<<<<<< HEAD\n[' in file_string:
                #     file_string = file_string.replace('<<<<<<< HEAD\n[', '[')
                # if '>>>>>>> master\n' in file_string:
                #     file_string = file_string.replace('>>>>>>> master\n', '')
                # with open(results_file, 'w') as f:
                #     f.write(file_string)
                trajectories += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(results_file, 'r'))]
            launch_dates = []
            tof = []
            delta_v = []
            for case in trajectories:
                launch_dates.append(datetime_to_jd(case.planet_dates[0]))
                tof.append((case.planet_dates[-1] - case.planet_dates[0]).days)
                delta_v.append(case.delta_v)
                if not isinstance(delta_v[-1], float):
                    print(delta_v[-1])
                if not isinstance(tof[-1], int):
                    print(tof[-1])
                if not isinstance(launch_dates[-1], float):
                    print(launch_dates[-1])
            plot_porkchop(launch_dates, tof, delta_v, itinerary)

    exit()

    folder = 'earth-mars-earth-jupiter-pluto'

    results_dir = os.path.join(base_path, folder, 'case_data')
    boosted_dir = os.path.join(base_path, folder, 'boosted_data')
    results_files = [os.path.join(results_dir, file) for file in os.listdir(results_dir) if 'cases' in str(file)]
    boosted_files = [os.path.join(boosted_dir, file) for file in os.listdir(boosted_dir) if 'cases' in str(file)]
    itinerary = folder.split('-')

    cases = []
    for results_file in results_files:
        cases += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(results_file, 'r'))]

    # boosted_cases = []
    # for boosted_file in boosted_files:
    #     boosted_cases += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(boosted_file, 'r'))]

    cases_porkchop = cases

    launch_dates = []
    tof = []
    delta_v = []
    for case in cases_porkchop:
        launch_dates.append(datetime_to_jd(case.planet_dates[0]))
        tof.append((case.planet_dates[-1] - case.planet_dates[0]).days)
        delta_v.append(case.delta_v)
        # launch_dates.append(case.planet_dates[0].year + case.planet_dates[0].month / 12.0 + case.planet_dates[0].day / 365.0)
        # tof.append((case.planet_dates[-1] - case.planet_dates[0]).days / 365.0)
        # delta_v.append(case.delta_v)

    plot_porkchop(launch_dates, tof, delta_v)
    # x = np.random.rand(100)
    # y = np.random.rand(100)
    # z = np.sin(6.0 * x) * np.cos(6.0 * y)
    # plot_porkchop(x, y, z)

