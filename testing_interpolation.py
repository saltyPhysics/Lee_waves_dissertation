#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:18:59 2018

@author: mdevana
"""

# Testing interpolation methods

# radial basis functions
from scipy.interpolate import Rbf
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
x, y, z, d = np.random.rand(4, 50)

lonpad = .1
latpad = .1

lonid1 = np.argmin(np.abs((lon-lonpad) - gem.bathy['lon'][:]))
lonid2 = np.argmin(np.abs((lon+lonpad) - gem.bathy['lon'][:]))

latid1 = np.argmin(np.abs((lat - latpad) - gem.bathy['lat'][:]))
latid2 = np.argmin(np.abs((lat + latpad) - gem.bathy['lat'][:]))

lonvec = gem.bathy['lon'][lonid1:lonid2+1]
latvec = gem.bathy['lat'][latid1:latid2+1]
bathy_subgrid = gem.bathy['elevation'][latid1:latid2+1,
                          lonid1:lonid2+1]

lonmeshb, latmeshb = np.meshgrid(lonvec, latvec)

test = CGx()
F = interp2d(lonvec, latvec, bathy_subgrid, kind='cubic')
tpad = 1

n2_in = self.N2[lon_id-lonpad:lon_id+lonpad+1,
                lat_id-latpad:lat_id+latpad+1,
                :, time_id-tpad:time_id+tpad+1]

lonvec = gem.lon[lon_id-lonpad:lon_id+lonpad+1]

# map the grid onto uniform grids
latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             lonvec,
                                             gem.N2grid[:],
                                             tvec, sparse=False)

mask = ~np.isnan(n2temp)
 # radial basis function interpolator instance
F1 = LinearNDInterpolator(test, n2temp[mask]) 

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


distance = (1e-3 * np.sqrt(results['x']**2 + results['y']**2))
f = gsw.f(results['lat'])**2
N2_grid = results['N2']

latlims = np.array([np.nanmin(results['lat']) - buffer,
                    np.nanmax(results['lat']) + buffer])
latlims = [np.argmin(np.abs(lat_in - gem.bathy['lat'][:])) 
                    for lat_in in latlims]
latlims = np.arange(latlims[0], latlims[1])

lonlims = np.array([np.nanmin(results['lon']) - buffer,
                    np.nanmax(results['lon']) + buffer])
lonlims = [np.argmin(np.abs(lon_in - gem.bathy['lon'][:])) 
                    for lon_in in lonlims]
lonlims = np.arange(lonlims[0], lonlims[1])

bathy_rev = gem.bathy['elevation'][latlims, lonlims]
lat_b = gem.bathy['lat'][latlims]
lon_b = gem.bathy['lon'][lonlims]

clevels = np.linspace(np.nanmin(bathy_rev), np.nanmax(bathy_rev), cls)

fig = plt.figure(figsize=(16, 10))

# Map Plot
plt.subplot(221)
plt.contour(lon_b, lat_b, bathy_rev, colors='k', levels=clevels)
plt.pcolormesh(lon_b, lat_b, bathy_rev, shading='gouraud')
plt.plot(results['lon'], results['lat'], c='r')
plt.scatter(results['lon'][0], results['lat'][0],
            marker='*', c='#00ff32', s=ms)
plt.scatter(results['lon'][-1], results['lat'][-1], 
                        marker='*', c='r', s=ms)


plt.figure()
plt.contour(lon_b, lat_b, bathy_rev, colors='k', levels=clevels)
plt.contourf(gem.lon, gem.lat, gem.N2[:, :, 20, 18].T)
plt.ylim([-53.8488,-53.36])
plt.xlim([-49.44, -48.9742])

