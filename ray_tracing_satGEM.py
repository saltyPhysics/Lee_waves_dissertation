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

    Parameters
    ----------
    Omega: intrinsic wave frequency
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    f: coriolis parameter
    N2: Buoyancy Frequency Squared
    w: vertical flow speed (default=0 because it is not available in satGEM)

    Returns
    -------
    cgz: Vertical Group Speed (m/s)

    """
    K2 = k**2 + l**2 + m**2

    cgz = (-1 * (k**2 + l**2) * m * (N2 - f**2)) / (K2**2 * Omega)

    return cgz



def CGx(N2, Omega, k, l, m, u, f):
    """
    Horizontal group speed in x-direction

    Parameters
    ----------
    N2: Buoyancy Frequency Squared
    Omega: intrinsic wave frequency
    k: zonal wavenumber
    l: meridional wavenumber 
    m: vertical wavenumber
    u: zonal flow speed
    f: coriolis parameter

    Returns
    -------
    cgx: Horizontal (x) Group Speed (m/s)


    """
#    K2 = k**2 + l**2 + m**2
    
    cgx = ((k * m**2 * (N2 - f**2))/((k**2 + l**2 + m**2)**2 * Omega)) + u
    
    return cgx


def CGy(N2, Omega, k, l, m, v, f):
    """
    Horizontal group speed in y-direction in a flow

    Parameters
    ----------
    N2: Buoyancy Frequency Squared
    Omega: intrinsic wave frequency
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    v: meridional flow speed
    f: coriolis parameter

    Returns
    -------
    cgy: Horizontal (y) Group Speed (m/s)

    """

    K2 = k**2 + l**2 + m**2
    
    cgy = (l * m**2 * (N2 - f**2))/(K2**2 * Omega) + v

    return cgy


def EoZ(N2, w0, f, ):
    """
    Wave ray energy when variations can only occur in the vertical (i.e. N2 and
    flow only vary with depth not horizontally) - Olbers 1981
    INCOMPLETE-----------
    Parameters
    ----------



    Returns
    -------

    """
    Ez = np.squeeze((w0**2 * (N2 - f**2))
                    / ((w0**2 - f**2)**(3 / 2) * (N2 - w0**2)**(1 / 2)))
    return Ez

def refraction(N, k, l, m, dN, di, Omega):
    """
    Refraction index of internal wave through stratification 

    Parameters
    ----------
    N: Buoyancy Frequency
    
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    u: zonal flow speed
    f: coriolis parameter
    Omega: intrinsic wave frequency

    Returns
    -------
    ri: refraction index for a single direction (x, y, or z)

    """

    K = k**2 + l**2 + m**2
    ri = ((N * (k**2 + l**2)) / (K * Omega)) * (dN / di)

    return ri
    
    
def dk(dU, dV,dx, k, l , m, dN, N, Omega):
    """
    Change of zonal wavenumber k in time

    Parameters
    ----------
    dU: change in U (zonal velocity)
    dV: change in V (meridional velocity)
    dx: x change (meters)
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    dN: change in buoyancy frequency
    N: Buoyancy frequency
    Omega: Intrinsic frequency

    Returns
    -------
    dk: Change in zonal wavenumber k

    """
    ri = refraction(N, k, l, m, dN, dx, Omega)

    dk = -1 * (ri + k * (dU/dx) + l * (dV/dx))

    return dk


def dl(dU, dV, dy, k, l, m, dN, N, Omega):
    """
    Change of meridional wavenumber l in time

    Parameters
    ----------
    dU: change in U (zonal velocity)
    dV: change in V (meridional velocity)
    dy: y change (meters)
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    dN: change in buoyancy frequency
    N: Buoyancy frequency
    Omega: Intrinsic frequency

    Returns
    -------
    dl: Change in meridional wavenumber l

    """
    ri = refraction(N, k, l, m, dN, dy, Omega)

    dl = -1 * (ri + k * (dU / dy) + l * (dV / dy))

    return dl


def dm(dU, dV, dz, k, l, m, dN, N, Omega):
    """
    Change of vertical wavenumber l in time

    Parameters
    ----------
    dU: change in U (zonal velocity)
    dV: change in V (meridional velocity)
    dz: z change (meters)
    k: zonal wavenumber
    l: meridional wavenumber
    m: vertical wavenumber
    dN: change in buoyancy frequency
    N: Buoyancy frequency
    Omega: Intrinsic frequency

    Returns
    -------
    dm: Change in vertical wavenumber m

    """
    ri = refraction(N, k, l, m, dN, dz, Omega)

    dm = -1 * (ri + k * (dU / dz) + l * (dV / dz))

    return dm


def dOmega(rx, ry, rz, k, l, dU, dV):
    """
    Change in intrinsic frequency / dispersion relation

    Paramaters
    ----------
    rx: wave refraction in x direction
    ry: wave refraction in y direction
    rz: wave refraction in z direction
    k: zonal wavenumber
    l: meridional wavenumber
    dU: change in U (zonal velocity)
    dV: change in V (meridional velocity)

    Returns
    -------
    dW: Change in frequency

    """

    dW = (rx + ry + rx) + k * dU + l * dV

    return dW



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
    Stores features of a wave as an object which can be input into 
    the ray tracing function
 
    """

    # Add functionality for a default Buoyancy Frequncy and Velocity Profile

    def __init__(self, k=(2*np.pi)/1000, l=(2*np.pi)/1000, 
                        t0=datetime(2012, 11, 2, 3, 0, 0),
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



    def help(self):
        """
        Print instructions on how to use wave class
        """
        text = '''
Instructions for using ray tracing model.
\nGenerate a wave with chosen properties or use the default Parameters
'''
        print(text)


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
    load in the satGEM data as an object (this might be wierd though because the h5py module loads in each file as an object so not sure...)
    
    The objects built in functions can then be used to easily access the data set without ever having to load the whole thing in.
    
    Also Contains bathymetry data from GEBCO which is also stored as an object
    rather than loaded in all at once. 
    
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
        # access bathy data like a dictionary with keys: elevation, lat, lon
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

        These are for the velocity grids which are offset
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
        Generate a sub grid around a chosen point and interpolate to get better resolution in space and time than satGEM
        INCOMPLETE --- need to add the interpolation
        - right now it returns the indicies of the subgrid
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




def run_tracing(wave, satGEM, time_direction='forward', 
                duration=24, tstep=10, status=False):
    """
    Runs ray tracing using the wave 
    objects and gem field objects with option for
    forward and backwards time finite differenceing steps. 

    Parameters
    ----------
    wave: Wave object with all initial wave
        features(i.e. wavenumbers, frequency, starting position)
    satGEM: satGEM object for accessing satGEM and bathymetry data
    time_direction: Direction in time to run ray tracing
                    ('reverse' or 'forward')
    duration: duration to run ray tracing (HOURS)
    tstep: time step (SECONDS)
    status: gives status updates (keep False unless you want a million messages)

    Returns
    -------
    wave.results: stores all results into a dictionary attached to 
        the initial wave object as 1D arrays
    
    results keys:
        x: x distance (meters) 
        y: y distance (meters)
        z: z distance (meters)
        k: zonal wavenumber
        l: meridional wavenumber
        m: vertical wavenumber
        omega: wave frequency
        lat: wave's latitude
        lon: wave's longitude
        cgx: horizontal (x direction) group speed
        cgy: horizontal (y direction) group speed
        cgz: vertical group speed
        bathy: bathymetry below wave
        N2: Buoyancy frequency along wave path
        N2_grid: pressure grid corresponding to N2
        time: time vector as datetime type
        elapsed_time: time vector as seconds elapsed 


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
        
        
        # Finite differencing 

        # X step
        dx = tstep * CGx(N2, Omega, k, l, m, u, f)
        x  = x + dx # use this form instead of x+= because it modifies old values
    
        
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

        
        # Update position
        lon = lon2 
        lat = lat2
        
        # find nearest location in bathymetry grid
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

        if np.abs(Omega) < np.abs(f*1):
            print('Wave Frequency below inertial Frequency')
            break
        
        # Keep as false it prints a million messages
        if status:
            print('\r{} % done'.format(100*(i/len(time))))

    # After ray tracing loop

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
    
    distance = 1e-3 * np.sqrt(results['x']**2 + results['y']**2)
    results['distance'] = distance

    # Had to add this condition for when the run quits out on first step
    if bathy:
        results['bathy'] = np.vstack(bathy)
        results['cgx'] = np.vstack(cgx)
        results['cgy'] = np.vstack(cgy)
        results['cgz'] = np.vstack(cgz)


    N2_all = np.vstack(N2_all)
    N2_grid = np.vstack(N2_grid)
    
    results['N2'] = N2_all
    results['N2_grid'] = N2_grid
    
    summary = """
Ray Tracing Summary
-------------------
Duration: {} hours
Time Step: {} seconds
Distance Traveled: {} km
Time Direction: {} 
""".format(
    duration,
    tstep,
    distance[-1],
    time_direction
).strip('[]')

    results['summary'] = summary

    # Store results onto wave object (just to see if this works for now)
    wave.results = results





def plot_results(results, gem, lc='#ff3f02',
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
        Time Step: {} Seconds
        Initial Parameters (Observed)
        k0 = {}, l0 = {}, m0 = {}, w0 = {}  
        """.format(
            np.abs(results['elapsed_time'][-1] / 3600),
            np.nanmean(np.abs(np.diff(results['elapsed_time'].flatten()))),
            results['k'][0], results['l'][0],
            results['m'][0], results['omega'][0]
            ).strip('[]')

    plt.suptitle(title, size=16, fontweight='bold')
    plt.tight_layout()
    fig.subplots_adjust(top=0.8)

                        

    return fig

def testing():
    """
    variables used for testing ray tracing (taken from transect observations)
    
    """
    
    towyo_date = datetime(2012, 2, 28, 21, 33, 44)
    gem = rt.satGEM_field()
    k0 = 0.000379449
    l0 = -0.000395896
    m0 = -0.0062492
    w0 = -0.00014730
    wave = Wave(k=k0, l=l0, m=m0, w0=w0, z0=1500, t0=towyo_date,
                lat=lat[:, 5], lon=lon[:, 5]
                )

    run_tracing(wave, gem, duration=24, tstep=30, time_direction='reverse')
    
  
    
