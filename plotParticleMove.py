#!/usr/bin/env /usr/local/anaconda3/bin/python3
#
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

# Open the netCDF file
ncfile = Dataset('particle-move.nc', 'r')

# Extract the data:
efieldx_data = ncfile.variables['Density rho'][:]

positions = ncfile.variables['Positions'][:]

xPositions = positions[:, 0]
yPositions = positions[:, 1]

# Close the netCDF file
ncfile.close()

plt.imshow(np.transpose(efieldx_data), cmap='seismic')

# Add a colorbar
plt.colorbar()
plt.xlabel('x axis')
plt.ylabel('y axis')

# Save the figure:
#plt.savefig('double-rho-erroneous.png', dpi=300)

plt.show()

# Positions scatterplot:
plt.scatter(xPositions[::10], yPositions[::10], s=4)
plt.xlabel('x position [m]')
plt.ylabel('y position [m]')
plt.grid()
# Show the cplot
plt.show()


