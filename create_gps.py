"""
 create_gps.py  -  description
 ---------------------------------------------------------------------------------
 create GPS source file
 ---------------------------------------------------------------------------------
 copyright            : (C) 2020 Valentina Fioretti
 email                : valentina.fioretti@inaf.it
 ----------------------------------------------
 Usage:
 python create_gps.py N_in N_runs target_type energy theta plist theta_bin angle_type 
 example:
 python create_gps.py 100000000 10 2 160 0.82 1 0.1 2
 ---------------------------------------------------------------------------------
 Parameters:
 - N_in = number of input particles
 - N_runs = number of simulation runs
 - target_type = target configuration [0 = eRosita, 1 = SPO]
 - energy = input energy at the target [keV]
 - theta = incident polar angle (from the slab plane)
 - plist = physics list on the target (0 = SPL, 1 = SS)
 - theta_bin = histogram bin
 - angle_type = angular distribution (0 = Beam, 1 = Iso, 2 = Real)
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2017/11/25: creation date
"""

from astropy.io import fits
import numpy as np
import math
import sys, os
import matplotlib.pyplot as plt


theta_deg = 0.0
phi_deg = 0.
theta = theta_deg*(np.pi/180.)
phi = phi_deg*(np.pi/180.)

# source height
h_s = 150.  # cm

# Global Geometry:
N_tray = 13  
N_layer = 2
N_strip = 3072
pitch = 0.121   # mm
Tray_side = 371.712  # mm

# Tracker geometry [mm]
Si_t = 0.410
Glue_t = 0.0
K_t = 0.05
CF_t = 0.5
Conv_t = 0.245

plane_distance = 18.7  # mm
dist_tray = 2.   # mm

Al_t = plane_distance - dist_tray - (Si_t + Glue_t + K_t + Conv_t) - (K_t + Glue_t + Si_t) - (CF_t*2.)
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '%                AGILE V2.0                    %'
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '% - Number of trays:', N_tray 
print '% - Number of strips:', N_strip
print '% - Pitch [mm]:', pitch
print '% - Tray side [mm]:', Tray_side 
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '% Tracker thicknesses:'
print '% - Silicon thickness [mm]:', Si_t
print '% - Glue thickness [mm]:', Glue_t
print '% - Kapton thickness [mm]:', K_t
print '% - Carbon fiber thickness [mm]:', CF_t
print '% - Converter (W) thickness [mm]:', Conv_t
print '% ----------------------------------------------'
print '% - Plane distance [mm]:', plane_distance
print '% - Trays distance [mm]:', dist_tray
print '% ----------------------------------------------'
print '% - Computed Al honeycomb thickness [mm]:', Al_t

Lower_module_t = Si_t + Glue_t + K_t + Conv_t
Lower_module_t_NoConv = Si_t + Glue_t + K_t

z_start = Lower_module_t_NoConv

Central_module_t = (CF_t*2.) + Al_t
Upper_module_t = K_t + Glue_t + Si_t

print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '% Tracker heights:'
print '% - Lower module height [mm]:', Lower_module_t
print '% - Central module height [mm]:', Central_module_t
print '% - Upper module height [mm]:', Upper_module_t
print '% - Tray height [mm]:', Lower_module_t + Central_module_t + Upper_module_t


TRK_t = Lower_module_t_NoConv + Central_module_t + Upper_module_t

for k in xrange(N_tray):
	if k > 2: 
		TRK_t = TRK_t + Lower_module_t + Central_module_t + Upper_module_t + dist_tray 
	else: 
		TRK_t = TRK_t + Lower_module_t_NoConv + Central_module_t + Upper_module_t + dist_tray


TRK_t = TRK_t - Lower_module_t_NoConv - Upper_module_t
z_end = TRK_t + z_start

print '% - Tracker height [mm]:', TRK_t
print '% - Tracker Z start [mm]:', z_start
print '% - Tracker Z end [mm]:', z_end

print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '% GPS Set-up for the point source position:'
print '% - theta [deg.]:', theta_deg
print '% - phi [deg.]:', phi_deg
print '% - source height [cm]:', h_s
print '% ----------------------------------------------'

# tracker height
h_t = z_end/10. # cm

# source height respect to tracker
h_r = h_s - h_t

# source distance from (0,0)
radius = h_r*np.tan(theta)
x_s = ((np.cos(phi))*radius)
y_s = ((np.sin(phi))*radius)


P_x = -np.sin(theta)*np.cos(phi)
P_y = -np.sin(theta)*np.sin(phi)
P_z = -np.cos(theta)

print '% Source position:'
print '% - X [cm]:', x_s
print '% - Y [cm]:', y_s
print '% - Z [cm]:', h_s
print '% - Source direction:'
print '% - P_x:', P_x
print '% - P_y:', P_y
print '% - P_z:', P_z
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

# Plane center distance from (0,0)
radius_plane = h_r*np.tan(theta)

# Plane center
c_x_plane = ((np.cos(phi))*radius_plane)
c_y_plane = ((np.sin(phi))*radius_plane)
halfx_plane = Tray_side/20. # cm
halfy_plane = Tray_side/20. # cm

# Plane Momenta

P_x = -np.sin(theta)*np.cos(phi)
P_y = -np.sin(theta)*np.sin(phi)
P_z = -np.cos(theta)

print '% Plane (square) center position:'
print '% - X [cm]:', c_x_plane
print '% - Y [cm]:', c_y_plane
print '% - Z [cm]:', h_s
print '% Plane (square) side position:'
print '% - Half X [cm]:', halfx_plane
print '% - Half Y [cm]:', halfy_plane
print '% Plane Source direction:'
print '% - P_x:', P_x
print '% - P_y:', P_y
print '% - P_z:', P_z
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
