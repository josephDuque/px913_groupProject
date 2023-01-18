#!/usr/bin/env /usr/local/anaconda3/bin/python3
#
#Import modules
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

plt.imshow(efieldx_data[::-1, :], cmap='seismic')

# Add a colorbar
plt.colorbar()

# Save the figure:
plt.savefig('double-rho-erroneous.png', dpi=300)

plt.show()

# Positions scatterplot:
#plt.scatter(xPositions, yPositions)
# Show the cplot
#plt.show()


