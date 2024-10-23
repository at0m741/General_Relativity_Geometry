
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("metric_results_3d.txt")

r = data[:, 0]
theta = data[:, 1]
g_rr = data[:, 2]

plt.contourf(r.reshape((256, 356)), theta.reshape((256, 356)), g_rr.reshape((256, 356)), cmap='inferno')
plt.colorbar(label='g_rr')
plt.xlabel('r')
plt.ylabel('theta')
plt.show()
