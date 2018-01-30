"""
Author: Manish Devana
Date: December 31, 2017

Ray Tracing functions for internal waves

The goal is to use an object oriented approach to a ray tracing model.
CURRENT STATUS:
load in satGEM data and rewrite functions to remove all assumptions and run in a 4d field. 

- Figure out how to make k and l and m vary in all 4 dimensions (x,y,z, and t)
- need a solid conversion method for x and y distances to lat and long (long is the tricker one)




"""


import numpy as np
import scipy
import pandas as pd
import gsw
import oceans as oc
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as colors
import matplotlib.dates as mdates
import cmocean
import h5py
from datetime import datetime, timedelta
from netCDF4 import Dataset





def instructions():
    """
    Print instructions
    """
    text = '''
Ray Tracing Instructions:
------------------------
1. Generate a "wave" object : rt.wave(inputs)
    - enter wave frequency, horizontal and vertical wavenumber components, and initial depth
    - view properties to check if things were entered correctly
    - when loading and calculating N2, take chunks out at a time otherwise it will crash. (too big of a file)

satGEM Details
--------------
\t This ray tracing model utilizes the 4D velocity and density field constructed by Dr. Andrew Meijers. Full details are available in Meijers 2013. 
    '''

    print(text)

# Wave ray tracing equations

def Kvec(k, l, m):
    """
    Returns magnitude of wavenumber vector
    """

    return k**2 + l**2 + m**2

def dispersion(f, N2, k, l, m):
    """
    WKB Disperision Relation as a function of and K(k, l, m)
    """

    W = np.sqrt((f**2 * m**2 + N2*(k**2 + l**2))\
                / ((k**2 +l**2 + m**2)))

    return W


def CGz(Omega, k, l, m, f, N2, w=0):
    """
    Vertical Group Speed (includes vertical flow but is 0 by default)
    """
    K2 = k**2 + l**2 + m**2

    return (-1*(k**2 + l**2) * m * (N2 - f**2)) / (K2**2 * Omega)


def CGx(N2, Omega, k, l, m, u, f):
    """
    Horizontal group speed in x-direction in a flow
    """
#    K2 = k**2 + l**2 + m**2
    
    cgx = ((k * m**2 * (N2 - f**2))/((k**2 + l**2 + m**2)**2 * Omega)) + u
    
    return cgx


def CGy(N2, Omega, k, l, m, v, f):
    """
    Horizontal group speed in y-direction in a flow
    """
    K2 = k**2 + l**2 + m**2
    
    cgy = (l * m**2 * (N2 - f**2))/(K2**2 * Omega) + v

    return cgy


def EoZ(N2, w0, f, ):
    """
    Wave ray energy when variations can only occur in the vertical (i.e. N2 and
    flow only vary with depth not horizontally) - Olbers 1981
    """
    Ez = np.squeeze((w0**2 * (N2 - f**2))
                    / ((w0**2 - f**2)**(3 / 2) * (N2 - w0**2)**(1 / 2)))
    return Ez

def refraction(N, k, l, m, dN, di, Omega):
    """
    Refraction index of internal wave
    """
    K = k**2 + l**2 + m**2
    
    return ((N*(k**2 + l**2)) / (K * Omega)) * (dN/di)
    
    
def dk(dU, dV,dx, k, l , m, dN, N, Omega):
    """
    Change of wavenumber k in time
    """
    ri = refraction(N, k, l, m, dN, dx, Omega)

    dk = -1 * (ri + k * (dU/dx) + l * (dV/dx))

    return dk


def dl(dU, dV, dy, k, l, m, dN, N, Omega):
    """
    Change of wavenumber k in time
    """
    ri = refraction(N, k, l, m, dN, dy, Omega)

    dl = -1 * (ri + k * (dU / dy) + l * (dV / dy))

    return dl


def dm(dU, dV, dz, k, l, m, dN, N, Omega):
    """
    Discretized Change of wavenumber k in time
    """
    ri = refraction(N, k, l, m, dN, dz, Omega)

    dm = -1 * (ri + k * (dU / dz) + l * (dV / dz))

    return dm


def dOmega(rx, ry, rz, k, l, dU, dV):
    """
    Change in intrinsic frequency / dispersion relation
    """

    dW = (rx + ry + rx) + k * dU + l * dV

    return dW





def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def inverse_hav(x, y, lon1, lat1):
    """
    Uses the inverse haversine function to convert x and y distance to a new lat and long coordinate. (see ray tracing docs for full formula)

    Parameters
    ----------
    x: x distance traveled (east-west)
    y: y distance traveled (north-south)
    lon1: starting longitude (Degrees)
    lat1: starting latitude (Degrees)

    Returns
    -------
    lon2: final longitude (Degrees)
    lat2: final latitude (Degrees)

    
    """

    r = 6371e3 # radius of the earth
    d = np.sqrt(x**2 + y**2) # total distance traveled

    lat2 =  lat1 + (y/111.11e3) # convert y distance to a new latitude point

    # Convert to radians for use in trig functions
    latrev1 = np.deg2rad(lat1)
    latrev2 = np.deg2rad(lat2)

    # inverse haversine formula
    shift = 0.5 * np.rad2deg(np.arccos(1 - 2 * ((np.sin(d / (2 * r))**2 
                - np.sin((latrev2 - latrev1)/2)**2) /
                (np.cos(latrev1) * np.cos(latrev2)))))
        
    if x < 0:
        lon2 = lon1 - shift
    else:
        lon2 = lon1 + shift
    
    


    return lon2, lat2 # in degrees




class Wave(object):
    """
    Creates a wave which has varying functionality including:
    - time forward modelling
    - time reverse modelling
    - variable velocity and density field inputs
    - plotting and saving rays
    - HOPEFULLY: built in gui
    """

    # Add functionality for a default Buoyancy Frequncy and Velocity Profile

    def __init__(self, k=10*1000, l=10*1000, t0=datetime(2012, 11, 2, 3, 0, 0),
                 m=500, w0=-1.3e-4, z0=500, lat=-55, lon=-55):

        # Convert wavelengths into wavenumbers
        # Save initial values becuase running the model will change
        # the wave features.
        self.k = np.array([k], dtype='float')
        self.l = np.array([l], dtype='float')
        self.m = np.array([m], dtype='float')
        self.w0 = np.array([w0], dtype='float')
        self.kh = np.array([np.sqrt(self.k**2 + self.l**2)])
        self.z0 = np.array([z0], dtype='float')
        self.lat0 = np.array([lat], dtype='float')
        self.lon0 = np.array([lon], dtype='float')
        self.t0 = t0

        # These are empty for now- they get filled in during model runs. 
        self.x_all = []
        self.y_all = []
        self.z_all = []
        self.m_all = []
        self.w0_all = []
        self.E_all = []
        self.Ac_all = []


    def help(self):
        """
        Print instructions on how to use wave class
        """
        text = '''
Instructions for using ray tracing model.
\nGenerate a wave with chosen properties or use the defualt Parameters
'''
        print(text)


    def model_error_message(self, x, y, z, m, idx, idx2):
        error_message = '''
        current variables:
-----------------
x = {}
y = {}
z = {}
N2 = {}
U = {}
V = {}
m = {}
'''.format(x, y, z,
 self.N2[idx2], self.U[idx], self.V[idx], m)

        return error_message

    def properties(self):
        """
        Print wave properties
        """
        txt = '''Wave Properties:
---------------
    k: {}
    l: {}
    m: {}
    kh: {}
    Frequency: {}
'''.format(np.array2string(self.k), self.l, self.m, self.kh, self.w0)

        print(txt)
  




class satGEM_field(object):
    """
    load in the satGEM data as an object (this might be wierd though becuase the h5py module loads in each file as an object so not sure...)
    
    The objects built in functions can then be used to easily access the data set without ever having to load the whole thing in.
    
    Also Contains bathymetry data
    
    """

    def __init__(self):
        # Load satGEM data as h5py file objects
        gamma_file = h5py.File('DIMES_GAMMA_09_12_upd.mat')
        vel_file = h5py.File('DIMES_vel_09_12_upd.mat')
        ts_file = h5py.File('DIMES_TS_09_12_upd.mat')
        gamma = gamma_file['satGEM_gamma']

        self.u = vel_file['satGEM_east']
        self.v = vel_file['satGEM_north']
        self.temp = ts_file['satGEM_temp']
        self.sal = ts_file['satGEM_sal']

        # Data grids
        time = np.squeeze(np.array(gamma_file['time']))
        # convert from matlab to python date time.
        self.time = np.array([oc.matlab2datetime(timeIn) for timeIn in time])

        self.depth_grid = gamma_file['depthlvl']
        self.lon = gamma_file['lons']
        self.lat = gamma_file['lats']

        # The u and v grids are one point off each so I need
        # to figure out how to handle this
        self.centerlat = vel_file['centerlat']
        self.centerlon = vel_file['centerlon']
        
        ###################################
        # Bathymetry file
        self.bathy = Dataset('bathy.nc')
       
        

    def locate(self, lon, lat, depth, time):
        """
        Locate point/points within the satgem data set

        Parameters
        ----------
        lon: longitude of point
        lat: latitude of point
        depth: depth of point
        time: of point

        Returns
        -------
        lon_id: index along longitude axis
        lat_id: index along latitude axis
        depth_id: index along latitude axis
        time_id: index along time axis

        These are for the velocity grids
        centerlon_id: index along centerlon axis
        centerlat_id: index along centerlat axis


        """

        # Add warning for out of time and space boundaries. 

        lon_id = np.argmin(np.abs(self.lon[:] - lon))
        lat_id = np.argmin(np.abs(self.lat[:] - lat))
        depth_id = np.argmin(np.abs(self.depth_grid[:] - depth))
        time_id = np.argmin(np.abs(self.time[:] - time))

        centerlon_id = np.argmin(np.abs(self.centerlon[:] - lon))
        centerlat_id = np.argmin(np.abs(self.centerlat[:] - lat))

        return lon_id, lat_id, depth_id, time_id, centerlon_id, centerlat_id
            

    def subgrid(self, lon_c, lat_c, z_c, time_c, k0, l0, m0,
                x_pads=2, y_pads=2, z_pads=2, t_pads=1):
        """
        Generate a sub grid around a chosen point and return the indices of that grid
        """

        x_pad = (2 * ((x_pads * np.pi) / k0))  # pad with 2x wavelength on that axis
        
        y_pad = (2 * ((y_pads * np.pi) / l0))

        lon_pad1, lat_pad1 = inverse_hav(x_pad, y_pad, lon_c, lat_c)
        lon_pad2, lat_pad2 = inverse_hav(-x_pad, -y_pad, lon_c, lat_c)
        
        
        # Not sure if there is going to be problems near surface?
        z_pad1 = z_c + (2 * ((z_pads * np.pi) / m0))
        z_pad2 = z_c - (2 * ((z_pads * np.pi) / m0))

        # time padding (1 pad = 1 week - time grid of satGEM is weekly)
        tpad = time_c + dt.timedelta(days=7*t_pads) # if backwards in time use negative pad
        
        lon_id1, lat_id1,\
            depth_id1,time_id1, centerlon_id1,\
            centerlat_id1 = self.locate(lon_pad1, lat_pad1, z_pad1, time_c)
            
        lon_id2, lat_id2, depth_id2, time_id2, centerlon_id2, centerlat_id2 = self.locate(lon_pad2, lat_pad2, z_pad2, tpad)
            
        return np.array([lon_id1, lon_id2]), np.array([lat_id1, lat_id2]),\
                np.array([depth_id1, depth_id2]), np.array([time_id1, time_id2])








def run_tracing(wave, satGEM, time_direction='reverse', 
                duration=24, tstep=10, status=True):
    """
    Runs ray tracing using the wave 
    objects and gem field objects with option for
    forward and backwards time finite differenceing steps. 
    """

    if not isinstance(wave, Wave):
        raise ValueError('Wave input must be a Wave object')

    if not isinstance(satGEM, satGEM_field):
        raise ValueError('satGEM input must be a satGEM field object')
    

    # get initial values from wave object
    k = wave.k[:]
    l = wave.l[:]
    m = wave.m[:]
    Omega = wave.w0
    lat = wave.lat0
    lon = wave.lon0
    z = wave.z0[:]
    x = float(0)
    y = float(0)

    x_all = []
    y_all = []
    z_all = []
    k_all = []
    l_all = []
    m_all = []
    om_all = []
    lat_all = []
    lon_all = []
    cgx = []
    cgy = []
    cgz = []
    bathy = []
    N2_all = []
    N2_grid = []


    # add start values (theres probably a better way to do this)
    x_all.append(x)
    y_all.append(y)
    z_all.append(z)
    k_all.append(k)
    l_all.append(l)
    m_all.append(m)
    om_all.append(Omega)
    lat_all.append(lat)
    lon_all.append(lon)

    start_time = wave.t0 


    # Time arrays and start time in satgem field. 
    if time_direction == 'reverse':
        end_time = start_time - timedelta(hours=duration)
        tstep = -tstep
    elif time_direction == 'forward':
        end_time = start_time + timedelta(hours=duration)
        
    else:
        raise ValueError("Invalid time direction \
                            - accepts 'forward' or 'reverse'")

        
    time = np.arange(start_time, end_time, timedelta(seconds=tstep)).astype(datetime) # create time vector (seconds)

    # start time, depth, lat, and lon index in satGEM field
    
    # lon_idx, lat_idx, z_idx,\
    #         t_idx, clon_idx, clat_idx = satGEM.locate(lon, lat, z, time[0])
    # Run loop for tracing

    for i, t in enumerate(time[:-1]):
        
        f = gsw.f(lat)
        
        # list with [lon, lat, depth, time, centerlon, centerlat] indices
        lon_idx, lat_idx, z_idx,\
            t_idx, clon_idx, clat_idx = satGEM.locate(lon, lat, z, t)
        
        # satGEM values
        # Vertical N2 profile at current location 
        N2_vert, p_grid = gsw.Nsquared(satGEM.sal[lon_idx, lat_idx, :, t_idx],
                              satGEM.temp[lon_idx, lat_idx, :, t_idx],
                              satGEM.depth_grid[:,:].flatten(),lat[:],
                              axis=-1)
        N2_all.append(N2_vert.flatten())
        N2_grid.append(p_grid.flatten())
        # REMOVE THIS WHEN DONE!!!
        # plt.plot(p_grid, N2_vert)

        idx_n = np.argmin(np.abs(p_grid.flatten() - z))
        N2 = N2_vert.flatten()[idx_n]

        u = satGEM.u[lon_idx, clat_idx, z_idx, t_idx]
        v = satGEM.v[clon_idx, lat_idx, z_idx, t_idx]

        # Check 1 (these have to be done before calculations)
        if not np.isfinite(N2):
            print('N2 error')
            x_all.append(x)
            y_all.append(y)
            z_all.append(z)
            k_all.append(k)
            l_all.append(l)
            m_all.append(m)
            om_all.append(Omega)
            lat_all.append(lat)
            lon_all.append(lon)
            cgx.append(dx / tstep)
            cgy.append(dy / tstep)
            cgz.append(dz / tstep)
            bathy.append(bottom)
            break

        
        if not np.isfinite(u):
            print('u error')
            x_all.append(x)
            y_all.append(y)
            z_all.append(z)
            k_all.append(k)
            l_all.append(l)
            m_all.append(m)
            om_all.append(Omega)
            lat_all.append(lat)
            lon_all.append(lon)
            cgx.append(dx / tstep)
            cgy.append(dy / tstep)
            cgz.append(dz / tstep)
            bathy.append(bottom)
            break
        
        if not np.isfinite(v):
            print('v error')
            x_all.append(x)
            y_all.append(y)
            z_all.append(z)
            k_all.append(k)
            l_all.append(l)
            m_all.append(m)
            om_all.append(Omega)
            lat_all.append(lat)
            lon_all.append(lon)
            cgx.append(dx / tstep)
            cgy.append(dy / tstep)
            cgz.append(dz / tstep)
            bathy.append(bottom)
            break
        
        
        # X step
        dx = tstep * CGx(N2, Omega, k, l, m, u, f)
        x  = x + dx # use this form instead of x+= because it modifies old values
#        print(np.sqrt(u**2 + v**2))

#        print('x step: {}'.format(dx))
        
        # Y step
        dy = tstep * CGy(N2, Omega, k, l, m, v, f)
        y = y + dy
        
        # Z step
        dz = tstep * CGz(Omega, k, l, m, f, N2)
        z = z + dz
        
        
        # New position
        lon2, lat2 = inverse_hav(dx, dy, lon, lat)
        
        lon_idx2, lat_idx2, z_idx2,\
            t_idx2, clon_idx2, clat_idx2 = satGEM.locate(lon2, lat2, z, time[i+1])
        
        # change in satGEM properties (U, V, N)
        N2_vert2, p_grid2 = gsw.Nsquared(satGEM.sal[lon_idx2, lat_idx2,
                                :, t_idx2],
                                satGEM.temp[lon_idx2, lat_idx2, :, t_idx2],
                                satGEM.depth_grid[:,:].flatten(),
                                axis=-1)

        idx_n2 = np.argmin(np.abs(p_grid2 - z))
        N2_2 = np.abs(N2_vert2[idx_n2])
        dN = np.sqrt(N2_2) - np.sqrt(np.abs(N2))


        u2 = satGEM.u[lon_idx2, clat_idx2, z_idx2, t_idx2]
        v2 = satGEM.v[clon_idx2, lat_idx2, z_idx2, t_idx2]
        
        # Changes in U
        du = u2 - u
        
        # V changes
        dv = v2 - v

        # k step
        k = k + dk(du, dv, dx, k, l, m, dN, np.sqrt(N2_2), Omega)

        # l step
        l = l + dl(du, dv, dy, k, l, m, dN, np.sqrt(N2_2), Omega)
    
        # m step
        m = m + dm(du, dv, dz, k, l, m, dN, np.sqrt(N2_2), Omega)
        
        # omega step 

        # Refraction of internal wave through changing stratification
        rx = refraction(np.sqrt(N2_2), k, l, m, dN, dx, Omega)
        ry = refraction(np.sqrt(N2_2), k, l, m, dN, dy, Omega)
        rz = refraction(np.sqrt(N2_2), k, l, m, dN, dz, Omega)

        Omega = Omega + dOmega(rx, ry, rz, k, l, du, dv)

        
        # boundary checks
        
        lon = lon2 
        lat = lat2
        
        idx1 = np.argmin(np.abs(lon - satGEM.bathy['lon'][:]))
        idx2 = np.argmin(np.abs(lat - satGEM.bathy['lat'][:]))
        
        bottom = -1*satGEM.bathy['elevation'][idx2, idx1]
        
        
        # store data 
        x_all.append(x)
        y_all.append(y)
        z_all.append(z)
        k_all.append(k)
        l_all.append(l)
        m_all.append(m)
        om_all.append(Omega)
        lat_all.append(lat)
        lon_all.append(lon)
        cgx.append(dx/tstep)
        cgy.append(dy/tstep)
        cgz.append(dz/tstep)
        bathy.append(bottom)

        # Check Parameters before next step
        if z > bottom:
            print('Wave hit seafloor')
            break

        if z < 0 :
            print('Wave hit surface')
            break

        # if np.abs(Omega) < np.abs(f):
        #     print('Wave Frequency below inertial Frequency')
        #     break
        
        if status:
            print('\r{} % done'.format(100*(i/len(time))))

    # store all results in dictionary (keeps things concise when using)
    elapsed_time = np.vstack([(timeIn - time[0]).total_seconds() 
                    for timeIn in time[:i + 2]])
    results = {
            'x': np.vstack(x_all),
            'y': np.vstack(y_all),
            'z': np.vstack(z_all),
            'k': np.vstack(k_all),
            'l': np.vstack(l_all),
            'm': np.vstack(m_all),
            'omega': np.vstack(om_all),
            'lon': np.vstack(lon_all),
            'lat': np.vstack(lat_all),
            'time': time[:i+1],
            'elapsed_time': elapsed_time

    }
    
    
    if bathy:
        results['bathy'] = np.vstack(bathy)
        results['cgx'] = np.vstack(cgx)
        results['cgy'] = np.vstack(cgy)
        results['cgz'] = np.vstack(cgz)


    N2_all = np.vstack(N2_all)
    N2_grid = np.vstack(N2_grid)
    
    results['N2'] = N2_all
    results['N2_grid'] = N2_grid
    
    return results



def plot_tracing(results, gem, plot_satgem=True):
    """
    standard plotting function for viewing ray tracing results
    """
    
    if not isinstance(gem, satGEM_field):
        raise ValueError('satGEM input must be a satGEM field object')
        
    

    
    fig1 = plt.figure(figsize = (15,8))
    plt.subplot(421)
    plt.plot(results['x'], -results['z'])
    plt.title('x vs z')
    
    plt.subplot(422)
    plt.plot(results['time'], -results['z'])
    plt.title('time vs z')
    
    plt.subplot(423)
    plt.plot(results['time'], np.sqrt(results['x']**2 + results['y']**2))
    plt.title('time vs dist')
    
    plt.subplot(425)
    buffer = 0.2
    latlims = np.array([np.nanmin(results['lat'])-buffer, np.nanmax(results['lat'])+buffer])
    latlims = [np.argmin(np.abs(lat_in - gem.bathy['lat'][:])) for lat_in in latlims]
    latlims = np.arange(latlims[0], latlims[1])
    
    lonlims = np.array([np.nanmin(results['lon'])-buffer, np.nanmax(results['lon'])+buffer])
    lonlims = [np.argmin(np.abs(lon_in - gem.bathy['lon'][:])) for lon_in in lonlims]
    lonlims = np.arange(lonlims[0], lonlims[1])

    bathy_rev = gem.bathy['elevation'][latlims, lonlims]
    lat_b = gem.bathy['lat'][latlims]
    lon_b = gem.bathy['lon'][lonlims]
    plt.plot(results['lon'], results['lat'])
    
    plt.pcolormesh(lon_b, lat_b, bathy_rev)
    plt.title('lon vs lat')
 
    plt.subplot(426)
    plt.plot(results['elapsed_time'], np.log10(np.abs(results['k'])))
    plt.title('time vs log10(k)')
    
    plt.subplot(427)
    plt.plot(results['elapsed_time'], np.log10(np.abs(results['m'])))
    plt.title('time vs log10(m)')
    
    plt.subplot(428)
    plt.plot(np.log10(np.abs(results['k'])), np.log10(np.abs(results['m'])))
    plt.title('log k vs  log m')
    
    
    plt.tight_layout()
   


    plt.figure(figsize=(8,8))
    
    plt.pcolormesh(lon_b, lat_b, bathy_rev)
    plt.contour(lon_b, lat_b, bathy_rev, colors='k')
    plt.plot(results['lon'], results['lat'], c='r')
    plt.title('lon vs lat')
    
    
    plt.figure(figsize=(8,8))
    
    dist = np.sqrt(results['x']**2 + results['y']**2)
    plt.plot(dist, results['z'], c='r')
    plt.plot(dist[:-1], results['bathy'], c='k')
    plt.title('transect ')
 

def dashboard(results, gem, lc='#ff3f02',
                lw=1.5, ms=20, buffer=0.2, cls=20, plot_satgem=True):
    """
    Quick function for plotting a ray path over local bathymetry with
    start and ends marked

    Parameters
    ----------


    Returns
    -------


    """

    distance = 1e-3 * np.sqrt(results['x']**2 + results['y']**2)
    f = gsw.f(results['lat'])
    N2_grid = results['N2']
    p_grid = results['N2_grid']
    # if plot_satgem:
    #     # Get N2 field from satgem data for plotting in background
    #     N2_grid = []
    #     p_grid = []

    #     for i in range(len(results['time'])):
    #         lon_idx, lat_idx, z_idx,\
    #             t_idx, clon_idx, clat_idx = gem.locate(results['lon'][i],
    #                                                     results['lat'][i],
    #                                                     results['z'][i],
    #                                                     results['time'][i],
    #             )

    #         N2, p_mid = gsw.Nsquared(gem.sal[lon_idx, lat_idx, :, t_idx],
    #                         gem.temp[lon_idx, lat_idx, :, t_idx],
    #                         gem.depth_grid[:, 0],
    #                         axis=-1)

    #         N2_grid.append(N2)
    #         p_grid.append(p_mid)
            
    #     N2_grid = np.vstack(N2_grid)
    #     p_grid = np.vstack(p_grid)
        



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


    # Generate transec bathymetry (interpolate to smooth)
    

    bathy1, idx1 = np.unique(results['bathy'], return_index=True)
    bathy1 = np.vstack([results['bathy'][index] for index in sorted(idx1)])
    x1 = np.vstack([distance[index] for index in sorted(idx1)])
    bathy_rev = np.interp(distance.flatten(),  x1.flatten(), bathy1.flatten())

    plt.subplot(222)
    if len(distance) == N2_grid.shape[0]:
        x2 = np.tile(distance, (1, N2_grid.shape[1]))
    else:    
        x2 = np.tile(distance[:-1], (1, N2_grid.shape[1]))
    plt.contourf(x2, p_grid, np.log10(np.abs(N2_grid)),
                            cmap=cmocean.cm.tempo)
    c1 = plt.colorbar()
    c1.set_label('Log10 (N2)')
    plt.plot(distance, results['z'], c=lc, linewidth=lw)
    plt.fill_between(distance.flatten(), bathy_rev.flatten(), 6000,
                     facecolor='#606060')

    plt.scatter(distance[0], results['z'][0],
                                     marker='*', c='#00ff32', s=ms)
    plt.scatter(distance[-1], results['z'][-1],
                marker='*', c='r', s=ms)

    plt.xlabel('Distance Traveled (km)')
    plt.ylabel('Depth (m)')
    plt.title('Depth vs. Distance')
    plt.gca().invert_yaxis()

    plt.subplot(223)
    plt.plot(results['elapsed_time']/3600, np.log10(np.abs(results['m'])),
                                c=lc, linewidth=lw)
    plt.gca().format_xdata = mdates.DateFormatter('%h')
    plt.xlabel('Time (Hours)')
    plt.ylabel('log (m)')
    plt.title(' Vertical Wavenumber vs. Time')

    plt.subplot(224)
    plt.plot(results['elapsed_time']/3600, np.log10(np.abs(results['omega'])),
                    c=lc, linewidth=lw, label=r'$\omega$')
    plt.plot(results['elapsed_time'] / 3600, np.log10(np.abs(f)),
                    c='k', linewidth=lw, label='f')
    plt.legend()
    plt.xlabel('Time (Hours)')
    plt.ylabel(r'log ($\omega$)')
    plt.title(' Frequency vs. Time')

    title = """
        Ray Tracing Results: Runtime = {} hours 

        Initial Parameters (Observed)
        k0 = {}, l0 = {}, m0 = {}, w0 = {}  
        """.format(
            np.abs(results['elapsed_time'][-1] / 3600),
            results['k'][0], results['l'][0],
            results['m'][0], results['omega'][0]
            ).strip('[]')

    plt.suptitle(title, size=16)
    plt.tight_layout()
    fig.subplots_adjust(top=0.8)

                        

    return fig

def testing():
    """
    Random variables used for testing ray tracing
    
    """
    
    k =  0.000403364
    l = -0.000173485
    m = -0.01256
    w0 = -0.00014848
    
    towyo_date =  datetime(2012,	2,	28,	21,	33,	44)
    
    gem = satGEM_field()
    wave = Wave(k=k, l=l, m=m, w0=w0, t0=towyo_date)
    
    results = run_tracing(wave, gem, duration=10*24, tstep=30, time_direction='reverse')
    
    plot_tracing(results, gem)
    
