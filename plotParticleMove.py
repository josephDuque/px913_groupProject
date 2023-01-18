#!/usr/bin/env /usr/local/anaconda3/bin/python3
#
#Import modules
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

# Open the netCDF file
ncfile = Dataset('particle-move.nc', 'r')

# Extract the data from the "Electric field" variable
efieldx_data = ncfile.variables['Electric field x-component'][:]

# Close the netCDF file
ncfile.close()

plt.imshow(efieldx_data[::-1, :], cmap='cool')

# Add a colorbar
plt.colorbar()

# Show the plot
plt.show()
