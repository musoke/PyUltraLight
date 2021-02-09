#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import math
import PyUltraLight
import pyfftw
import os
import sys
import multiprocessing
import numpy
import numba
import h5py
from IPython.core.display import clear_output, display

try:
    resol = int(sys.argv[1])
except:
    print("Resolution not supplied")
    exit()

# # Set Axion Mass

# In[2]:


axion_mass = 1e-22 *1.783e-36 #kg


# # Set Simulation Parameters

# In[3]:


# Set number of threads to target
if 'OMP_NUM_THREADS' in os.environ:
    # If slurm is used it sets this to the allocated number of cores
    num_threads = int(os.environ['OMP_NUM_THREADS'])
else:
    num_threads = multiprocessing.cpu_count()


# Set units for soliton parameters
s_mass_unit = ''     #Accepted units: 'kg', 'solar_masses', 'M_solar_masses', and '' for dimensionless units
s_position_unit = '' #Accepted units: 'm', 'km', 'pc', 'kpc', 'Mpc', 'ly', and '' for dimensionless units
s_velocity_unit = '' #Accepted units: 'm/s', 'km/s', 'km/h', and '' for dimensionless units

# Set box size and resolution
length = 10 # 1 code unit is ~38 kpc x (1e-22/m_a)^0.5
length_units = ''  # Accepted units: 'm', 'km', 'pc', 'kpc', 'Mpc', 'ly', and '' for dimensionless units.
# resol= 128 # It is recommended to check the upper bound on soliton mass for a given box size and resolution
duration = 0.1 # 1 code unit is ~70 Gyr (independent of axion mass assumption)
duration_units = ''  # Accepted units: 's', 'yr', 'kyr', 'Myr', and '' for dimensionless units
start_time = 0.0 # Should be given in the same units as duration. 
central_mass = 0. # Give this parameter in the same units as the soliton mass unit. i.e. units must match with s_mass_unit

#Data to save
save_rho = False # Saves density data for entire 3D simulation grid
save_psi = False # Saves full complex field data for entire 3D simulation grid
save_plane = False # Saves density data for plane z = 0
save_energies = False # Saves integrated gravitational, kinetic and total energies as lists
save_line = False # Saves density data for line y = 0, z = 0. Useful for examining intereference patterns. 

#Formats to save
hdf5 = False
npz = False
npy = True

step_factor = 1. # Change this to a larger number if velocities are sufficiently low that constraint on timestep can be relaxed. 
save_number = 10    # Choose number of 'frames' to save. Note that, depending on resolution, this could require significant disk space.
save_path = 'TestOutput'  # Set output directory

save_options = [save_rho,save_psi,save_plane,save_energies,save_line]


# # Set Initial Conditions:

# In[6]:


m = 8 #1 code unit is ~2.3e6 M_sol (1e-22/m_a)^1.5
r = 2 #1 code unit is ~38 kpc x (1e-22/m_a)^0.5
#v = np.sqrt(central_mass/r)

#Soliton parameters are mass, position, velocity and phase (radians)
soliton1 = [m, [r,0,0], [-20,0,0], 0]
soliton2 = [9, [-r,0,0], [20,0,0], 0]

solitons = [soliton1,soliton2]  


# # Run:

# In[8]:


time_per_step = PyUltraLight.evolve(central_mass, num_threads, length, length_units, resol, duration, duration_units, step_factor, save_number, save_options, save_path, npz, npy, hdf5, s_mass_unit, s_position_unit, s_velocity_unit, solitons, start_time)


print("{},{},{}".format(num_threads, resol, time_per_step))
