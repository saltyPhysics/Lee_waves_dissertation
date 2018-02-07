#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:18:59 2018

@author: mdevana
"""

# Testing interpolation methods

# radial basis functions
from scipy.interpolate import Rbf
x, y, z, d = np.random.rand(4, 50)

lonpad = 1
latpad = 1
tpad = 1

n2_in = self.N2[lon_id-lonpad:lon_id+lonpad+1,
                lat_id-latpad:lat_id+latpad+1,
                :, time_id-tpad:time_id+tpad+1]

lonvec = gem.lon[lon_id-lonpad:lon_id+lonpad+1]

# map the grid onto uniform grids
lonmesh, latmesh, dmesh, tmesh = np.meshgrid(lonvec,
                                             latvec,
                                             gem.N2grid[:],
                                             tvec, sparse=False)

mask = ~np.isnan(n21)
 # radial basis function interpolator instance
rbfi = Rbf(lonmesh[mask],
           latmesh[mask],
           dmesh[mask],
           tmesh[mask], n21[mask]) 

test_coords = (-49.2, -53.7, 1500, 734555.2)
rbfi(-49.2, -53.7, 1500, 734555.2)

# test using rbf on whole satGEM subset (for speed and memory usage)
latmesh, lonmesh, dmesh, tmesh = np.meshgrid(gem.lat,
                                             gem.lon,
                                             gem.N2grid[:],
                                             gem.timevec[:])

mask = ~np.isnan(gem.N2[:])

rbfi = Rbf(lonmesh[mask],
           latmesh[mask],
           dmesh[mask],
           tmesh[mask],
           gem.N2[mask]) 

