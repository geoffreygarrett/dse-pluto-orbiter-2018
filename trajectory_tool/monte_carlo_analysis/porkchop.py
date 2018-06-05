
import numpy as np
from scipy.interpolate import interp2d, Rbf
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


x = np.random.rand(100)
y = np.random.rand(100)
z = np.sin(6.0*x) * np.cos(6.0*y)
f = Rbf(x, y, z, epsilon=1)
# f = interp2d(x, y, z, kind='cubic')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x, y, z)
# plt.show()
# exit()


x_grid, y_grid = np.meshgrid(np.arange(0, 1, 0.005), np.arange(0, 1, 0.005))
z_grid = np.zeros_like(x_grid)
for i in range(x_grid.shape[0]):
    for j in range(x_grid.shape[1]):
        z_grid[i, j] = np.clip(f(x_grid[i, j], y_grid[i, j]), -1.0, 1.0)

mask = z_grid < 0.0
z_grid = np.ma.array(z_grid, mask=mask)

corner_mask = True
cs = plt.contourf(x_grid, y_grid, z_grid, corner_mask=corner_mask)
cs2 = plt.contour(cs, colors='k')
plt.title('corner_mask = {0}'.format(corner_mask))

cbar = plt.colorbar(cs)
cbar.ax.set_ylabel('delta V')
# Add the contour line levels to the colorbar
cbar.add_lines(cs2)

# Plot grid.
plt.grid(c='k', ls='-', alpha=0.3)
plt.show()




