"""
Author: Manish Devana

Ray tracing using satGEM fields

CURRENT STATUS:
load in satGEM data and rewrite functions to remove all assumptions and run in a 4d field. 

- How to choose best interpolation method and how to choose best time step
- add more robust time stepping 




"""


import numpy as np
import gsw
import oceans as oc
from oceans import LinearNDInterpolatorExt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import h5py
from datetime import datetime, timedelta
from netCDF4 import Dataset
from scipy.interpolate import Rbf, LinearNDInterpolator, interp2d, NearestNDInterpolator
from mpl_toolkits.mplot3d import Axes3D





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

def datetime2ordinal(time):
    """
    convert datetime to ordinal number with decimals 
    
    
    """
    
    ordinal_time = time.toordinal() + time.hour/24 + time.second/(24*60*60)
    
    return ordinal_time


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


def action():
    """
    Adding Wave Action
    """


def refraction(N, k, l, m, dNdi, Omega):
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
    ri = ((N * (k**2 + l**2)) / (K * Omega)) * (dNdi)

    return ri
    
    
def del_wavenum(dudi, dvdi, k, l , m, N, dndi, Omega):
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
    ri = refraction(N, k, l, m, dndi, Omega)

    dk = -1 * (ri + k * (dudi) + l * (dvdi))

    return dk

def dk(dudx, dvdx, k, l , m, dndx, Omega):
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


def dOmega(rx, ry, rz, k, l, dudt, dvdt):
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

    dW = (rx + ry + rx) + k * dudt + l * dudt

    return dW


def plane_wave(k, l, m, omega, f, t, p0, N2):
    """
    Plane wave fits of velocity and density perturbations from polarization relations (Cusack 2017)
    """

    # Amplitudes (from polarization relations)
    u0 = p0*((k*omega + l*fj)/(omega**2 + f**2))
    v0 = p0*((l*omega + k*fj)/(omega**2 + f**2))
    w0 = p0 * ((-m*omega) / (N2 - omega**2))
    # b0 = 
    # p0 = 

    return None


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


def generate_rbf():
    """
    function for generating radial basis function within satgem fields
    using scipy's existing rbf function
    """
    
    

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
        self.z0 = np.asscalar(np.array([z0], dtype='float'))
        self.lat0 = np.asscalar(np.array([lat], dtype='float'))
        self.lon0 = np.asscalar(np.array([lon], dtype='float'))
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
  


class satGEM_mini(object):
    """
    --------------------TRIMMED VERSION OF SATGEM---------------------

    load in the satGEM data as an object (this might be wierd though because the h5py module loads in each file as an object so not sure...)
    
    The objects built in functions can then be used to easily access the data set without ever having to load the whole thing in.
    
    Also Contains bathymetry data from GEBCO which is also stored as an object
    rather than loaded in all at once. 
    
    """

    def __init__(self):
        # Load satGEM data as h5py file objects
        file = h5py.File('satGEM_update_md.mat')

        self.u = file['urev']
        self.v = file['vrev']
        self.N2 = file['n2']
        self.N2grid = file['n2grid']

        # Data grids
        time = np.squeeze(np.array(file['trev']))
        # convert from matlab to python date time.
        self.time = np.array([oc.matlab2datetime(timeIn) for timeIn in time])
        self.timevec = np.array([t.toordinal() for t in self.time])

        self.depth_grid = file['depthlvl']
        self.lon = file['lons2'][:].flatten()
        self.lat = file['lats2'][:].flatten()

        # The u and v grids are one point off each so I need
        # to figure out how to handle this
        # u uses centered lat
        self.centerlat = file['clats2'][:].flatten()
        # v uses centered lon
        self.centerlon = file['clons2'][:].flatten()

        # gradients
        self.dudx = file['dudx']
        self.dvdx = file['dvdx']
        self.dndx = file['dndx']
        
        self.dudy = file['dudy']
        self.dvdy = file['dvdy']
        self.dndy = file['dndy']
        
        self.dudz = file['dudz']
        self.dvdz = file['dvdz']
        self.dndz = file['dndz']

        self.dudt = file['dudt']
        self.dvdt = file['dvdt']

        # bathymetry data
        self.bathy = Dataset('bathy.nc')
        
        self.ngrids = [self.N2, self.dndx, self.dndy, self.dndz]
        self.ugrids = [self.u, self.dudx, self.dudy, self.dudz, self.dudt]
        self.vgrids = [self.v, self.dvdx, self.dvdy, self.dvdz, self.dvdt]

    def _locate_variables(self, lon, lat, depth, time,
                         tres=1, xres=10, yres=10, zres=5):
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
        
        # get indices
        lon_id = np.argmin(np.abs(self.lon[:] - lon))
        lat_id = np.argmin(np.abs(self.lat[:] - lat))
        depth_id = np.argmin(np.abs(self.depth_grid[:] - depth))
        time_id = np.argmin(np.abs(self.time[:] - time))

        clon_id = np.argmin(np.abs(self.centerlon[:] - lon))
        clat_id = np.argmin(np.abs(self.centerlat[:] - lat))
        
        # interpolate to position
        
        # y interps
        lonpad = 1
        latpad = 1
        tpad = 1
        
        n2_in = self.N2[lon_id-lonpad:lon_id+lonpad+1,
                        lat_id-latpad:lat_id+latpad+1,
                        :, time_id-tpad:time_id+tpad+1]
        
        lonvec = self.lon[lon_id-lonpad:lon_id+lonpad+1]
        latvec = self.lat[lat_id-latpad:lat_id+latpad+1]
        tvec = self.timevec[time_id-tpad:time_id+tpad+1]

        return lon_id, lat_id, depth_id, time_id, clon_id, clat_id


    def interpfield_rbf(self, lon, lat, depth, time,
                    lonpad=0.5, latpad=0.5, tpad=1,
                    function='linear'):
        """
        Generate 4D radial basis function for interpolating
        ray tracing points
        
        Status: The ray tracing works welll until reaching the bottom of the 
        satgem data where vertical group speed abruptly goes to zero.
        This happens whether using linearND interpolator or the radial basis function
        The best solution may be to fill the satgem datasets with a mean bottom 
        value based on the last available value for each cast. 
        
        
    
        """
        
        tpad  = timedelta(weeks=tpad)
            
        # get indices for center of interpolation field
        lon_id1, lat_id1, depth_id1, time_id1, clon_id1, clat_id1 = \
                    self._locate_variables(lon-lonpad,
                                           lat-latpad,
                                           depth,
                                           time-tpad)
        
        lon_id2, lat_id2, depth_id2, time_id2, clon_id2, clat_id2 = \
                    self._locate_variables(lon+lonpad,
                                           lat+latpad,
                                           depth,
                                           time+tpad)
        
        # make subgrid
        n2_sub = self.N2[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
        
        lonvec = self.lon[lon_id1:lon_id2+1]
        latvec = self.lat[lat_id1:lat_id2+1]
        tvec = self.timevec[time_id1:time_id2+1]
        
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             lonvec,
                                             self.N2grid[:],
                                             tvec, sparse=False)
#        mask = ~np.isnan(n2_sub)
#        # radial basis function interpolator instance
#        self.fn2 = Rbf(lonmesh[mask],
#                   latmesh[mask],
#                   dmesh[mask],
#                   tmesh[mask], 
#                   n2_sub[mask],
#                   function='linear')
        # order = [regular, dx, dy, dz, dt]
        names = ['Fn2', 'Fdndx', 'Fdndy', 'Fdndz']
        n2temp = self.N2[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
        
        for i, grid in enumerate(self.ngrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)
            
            F = Rbf(lonmesh[mask],
                   latmesh[mask],
                   dmesh[mask],
                   tmesh[mask], 
                   subgrid[mask],
                   function=function)
            
            setattr(self,names[i], F)
        
        # U functions -  u uses centered latitude grid
        latvec = self.centerlat[clat_id1:clat_id2+1]
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             lonvec,
                                             self.depth_grid[:],
                                             tvec, sparse=False)
        
        names = ['Fu', 'Fdudx', 'Fdudy', 'Fdudz','Fdudt']
        for i, grid in enumerate(self.ugrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)       
            F = Rbf(lonmesh[mask],
                   latmesh[mask],
                   dmesh[mask],
                   tmesh[mask], 
                   subgrid[mask],
                   function=function)
            
            setattr(self,names[i], F)
            
        # V functions -  u uses centered latitude grid
        clonvec = self.centerlon[clon_id1:clon_id2+1]
        
        # remake lat vector
        latvec = self.lat[lat_id1:lat_id2+1]
        
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             clonvec,
                                             self.depth_grid[:],
                                             tvec, sparse=False)
        # create rbf functions for each paramater F(lon, lat, depth, time)
        names = ['Fv', 'Fdvdx', 'Fdvdy', 'Fdvdz','Fdvdt']
        for i, grid in enumerate(self.vgrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)       
            F = Rbf(lonmesh[mask],
                   latmesh[mask],
                   dmesh[mask],
                   tmesh[mask], 
                   subgrid[mask],
                   function=function)
            
            setattr(self,names[i], F)
        
        # return functions' boundaries 
        latlims = [np.nanmin(latvec), np.nanmax(latvec)]
        lonlims = [np.nanmin(lonvec), np.nanmax(lonvec)]
        tlims = [np.nanmin(tvec), np.nanmax(tvec)]
        
        # bathy interpolation function
        lonid1 = np.argmin(np.abs((lon-lonpad) - self.bathy['lon'][:]))
        lonid2 = np.argmin(np.abs((lon+lonpad) - self.bathy['lon'][:]))
        
        latid1 = np.argmin(np.abs((lat - latpad) - self.bathy['lat'][:]))
        latid2 = np.argmin(np.abs((lat + latpad) - self.bathy['lat'][:]))
        
        lonvec = self.bathy['lon'][lonid1:lonid2+1]
        latvec = self.bathy['lat'][latid1:latid2+1]
        bathy_subgrid = self.bathy['elevation'][latid1:latid2+1,
                                  lonid1:lonid2+1]
        
        
        
        F = interp2d(lonvec, latvec, bathy_subgrid, kind='cubic')
            
        setattr(self,'bathyF', F)
        
        
        return lonlims, latlims, tlims
    
    def interpfield_lnd(self, lon, lat, depth, time,
                    lonpad=0.5, latpad=0.5, tpad=1):
        """
        Generate 4D radial basis function for interpolating
        ray tracing points
        
        """
        
        tpad  = timedelta(weeks=tpad)
            
        # get indices for center of interpolation field
        lon_id1, lat_id1, depth_id1, time_id1, clon_id1, clat_id1 = \
                    self._locate_variables(lon-lonpad,
                                           lat-latpad,
                                           depth,
                                           time-tpad)
        
        lon_id2, lat_id2, depth_id2, time_id2, clon_id2, clat_id2 = \
                    self._locate_variables(lon+lonpad,
                                           lat+latpad,
                                           depth,
                                           time+tpad)
        
        # make subgrid
        
        lonvec = self.lon[lon_id1:lon_id2+1]
        latvec = self.lat[lat_id1:lat_id2+1]
        tvec = self.timevec[time_id1:time_id2+1]
        
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             lonvec,
                                             self.N2grid[:],
                                             tvec, sparse=False)


        # order = [regular, dx, dy, dz, dt]
        names = ['Fn2', 'Fdndx', 'Fdndy', 'Fdndz']

        
        for i, grid in enumerate(self.ngrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)
            points = np.vstack([lonmesh[mask],
                                latmesh[mask],
                                dmesh[mask],
                                tmesh[mask]]).T
            
            F = LinearNDInterpolatorExt(points,
                   subgrid[mask])
            
            setattr(self,names[i], F)
        
        # U functions -  u uses centered latitude grid
        latvec = self.centerlat[clat_id1:clat_id2+1]
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             lonvec,
                                             self.depth_grid[:],
                                             tvec, sparse=False)
        
        names = ['Fu', 'Fdudx', 'Fdudy', 'Fdudz','Fdudt']
        for i, grid in enumerate(self.ugrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)
            points = np.vstack([lonmesh[mask],
                                latmesh[mask],
                                dmesh[mask],
                                tmesh[mask]]).T
            
            F = LinearNDInterpolatorExt(points,
                   subgrid[mask])
            
            setattr(self,names[i], F)
            
        # V functions -  u uses centered latitude grid
        clonvec = self.centerlon[clon_id1:clon_id2+1]
        
        # remake lat vector
        latvec = self.lat[lat_id1:lat_id2+1]
        
        latmesh, lonmesh, dmesh, tmesh = np.meshgrid(latvec,
                                             clonvec,
                                             self.depth_grid[:],
                                             tvec, sparse=False)
        # create rbf functions for each paramater F(lon, lat, depth, time)
        names = ['Fv', 'Fdvdx', 'Fdvdy', 'Fdvdz','Fdvdt']
        for i, grid in enumerate(self.vgrids):
            subgrid = grid[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
            
            mask = ~np.isnan(subgrid)
            points = np.vstack([lonmesh[mask],
                                latmesh[mask],
                                dmesh[mask],
                                tmesh[mask]]).T
            
            F = LinearNDInterpolatorExt(points,
                   subgrid[mask])
            
            setattr(self,names[i], F)
        
        # return functions' boundaries 
        latlims = [np.nanmin(latvec), np.nanmax(latvec)]
        lonlims = [np.nanmin(lonvec), np.nanmax(lonvec)]
        tlims = [np.nanmin(tvec), np.nanmax(tvec)]
        
        subgrid = self.N2[lon_id1:lon_id2+1,
                        lat_id1:lat_id2+1,
                        :, time_id1:time_id2+1]
        
        # create interpolation function for bathymetry
        lonid1 = np.argmin(np.abs((lon-lonpad) - self.bathy['lon'][:]))
        lonid2 = np.argmin(np.abs((lon+lonpad) - self.bathy['lon'][:]))
        
        latid1 = np.argmin(np.abs((lat - latpad) - self.bathy['lat'][:]))
        latid2 = np.argmin(np.abs((lat + latpad) - self.bathy['lat'][:]))
        
        lonvec = self.bathy['lon'][lonid1:lonid2+1]
        latvec = self.bathy['lat'][latid1:latid2+1]
        bathy_subgrid = self.bathy['elevation'][latid1:latid2+1,
                                  lonid1:lonid2+1]
        
        
        
        F = interp2d(lonvec, latvec, bathy_subgrid, kind='cubic')
            
        setattr(self,'bathyF', F)
        
        
        
        
        return lonlims, latlims, tlims
    
def ray_tracing_interp(wave, gem, time_direction='forward', 
                duration=24, tstep=10,
                latpad=1, lonpad=1, tpad=1,
                extend_bathy = 10000,
                interp_mode='lnd',
                interp_function='linear'):
    """
    Ray tracing using the satGEM density and velocity fields with
    local interpolation to avoid strange steps in ray property evolutions
    
    Parameters
    ----------
    wave: Wave object
    gem: satGEM (and bathymetry) data loaded as satGEM_mini object
    time_direction: time direction of 
                    ray tracing ('forward' (default) or 'reverse'))
    duration: duration of run (HOURS) 
    tstep: time step (SECONDS)
    latpad: Latitude padding for subgrid (bigger padding = more satGEM data 
            available but slower runtimes and possible memory issues)
    lonpad: Longitude padding for subgrid (bigger padding = more satGEM data 
            available but slower runtimes and possible memory issues)
    tpad: Time padding for subgrid (bigger padding = more satGEM data 
            available but slower runtimes and possible memory issues)
    extend_bathy: distance(METERS) to extend bathymetry data past end of run
    interp_mode: Method for interpolating and extrapolating through satGEM 
                default = 'lnd', 'rbf'
                Ray tracing is very sensitive to the methods with LND (linear interpolation) being the most consistent (other versions) need 
                further working
    
    Returns
    -------
    results: Dictionary of all results 
    
    
    """
    
    # argument checks
    if not isinstance(wave, Wave):
        raise ValueError('Wave input must be a Wave object')

    if not isinstance(gem, satGEM_mini):
        raise ValueError('satGEM input must be a satGEM field object')
        
        
        
        
    # get initial values from wave object
    k = wave.k[:]
    l = wave.l[:]
    m = wave.m[:]
    Omega = wave.w0
    lat = wave.lat0
    lon = wave.lon0
    z = wave.z0
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
    u_all = []
    v_all = []
    N2_all = []
    rx_all = []
    ry_all = []
    rz_all = []
    dudz_all = []
    dvdz_all = []
    
#    N2_all = []
#    N2_grid = []


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

        
    time = np.arange(start_time, end_time,
                     timedelta(seconds=tstep)).astype(datetime) # create time vector (seconds)
    N2subgrid=None
    if interp_mode == 'rbf':
        lonlims, latlims, tlims = gem.interpfield_rbf(lon, lat,
                                              z, time[0],
                                              function=interp_function)
    elif interp_mode == 'lnd':
        lonlims, latlims, tlims = \
                        gem.interpfield_lnd(lon, lat,
                                            z, time[0],
                                              )

    for i, t in enumerate(time[:-1]):
        
        f = gsw.f(lat)
        
        # list with [lon, lat, depth, time, centerlon, centerlat] indices
#        lon_idx, lat_idx, z_idx,\
#            t_idx, clon_idx, clat_idx = satGEM.locate(lon, lat, z, t)
        
        tnum = datetime2ordinal(t)
        
        # if depth is greater than 3000 and wave hasnt hit bottomr
        # use 3000 (max range of parameter functions)
        # this is assuming that below 3000m values are constant 
        latid = np.argmin(np.abs(lat -  gem.lat))
        lonid = np.argmin(np.abs(lon - gem.lon))
        timeid = np.argmin(np.abs(gem.time[:] - t))
        N2closest = gem.N2[lonid, latid, :, timeid]
        gembottom = gem.depth_grid[np.max(np.where(np.isfinite(N2closest)))]
        if z > gembottom:
            z1 = gembottom
        else:
            z1 = z
        u = gem.Fu(lon, lat, z1, tnum)
        dudx = gem.Fdudx(lon, lat, z1, tnum)

        dudy = gem.Fdudy(lon, lat, z1, tnum)
        dudz = gem.Fdudz(lon, lat, z1, tnum)
        dudt = gem.Fdudt(lon, lat, z1, tnum)
        if np.isnan(dudt):
            print('dudt ERROR')
            break 
        
        v = gem.Fv(lon, lat, z1, tnum)
        dvdx = gem.Fdvdx(lon, lat, z1, tnum)
        dvdy = gem.Fdvdy(lon, lat, z1, tnum)
        dvdz = gem.Fdvdz(lon, lat, z1, tnum)
        dvdt = gem.Fdvdt(lon, lat, z1, tnum)
        if np.isnan(dvdt):
            print('dvdt ERROR')
            print(tnum)
            print(lon)
            print(lat)
            print(z1)
            break 
        
        N2 = np.abs(gem.Fn2(lon, lat, z1, tnum))
        dndx = gem.Fdndx(lon, lat, z1, tnum)
        dndy = gem.Fdndy(lon, lat, z1, tnum)
        dndz = gem.Fdndz(lon, lat, z1, tnum)
        
        

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
        if np.isnan(CGx(N2, Omega, k, l, m, u, f)):
            print('dx Error')
            break
        x  = x + dx # use this form instead of x+= because it modifies old values
    
        
        # Y step
        dy = tstep * CGy(N2, Omega, k, l, m, v, f)
        if np.isnan(dy):
            print('dy Error')
            break
        y = y + dy
        
        # Z step
        dz = tstep * CGz(Omega, k, l, m, f, N2)
        z = z + dz
        z = np.asscalar(z)
        
        # k step
        k = k + del_wavenum(dudx, dvdx, k, l, m, np.sqrt(N2), dndx, Omega)

        # l step
        l = l + del_wavenum(dudy, dvdy, k, l, m, np.sqrt(N2), dndy, Omega)
    
        # m step
        m = m + del_wavenum(dudz, dvdz, k, l, m, np.sqrt(N2), dndz, Omega)
        


        # Refraction of internal wave through changing stratification
        rx = refraction(np.sqrt(N2), k, l, m, dndx, Omega)
        if np.isnan(rx):
            print['RX ERROR']
            break 
        ry = refraction(np.sqrt(N2), k, l, m, dndy, Omega)
        if np.isnan(ry):
            print['Ry ERROR']
            break 
        rz = refraction(np.sqrt(N2), k, l, m, dndz, Omega)
        if np.isnan(rz):
            print['Rz ERROR']
            break 
        

        Omega = Omega + dOmega(rx, ry, rz, k, l, dudt, dvdt)

        
        # Update position
        lon2, lat2 = inverse_hav(dx, dy, lon, lat)
        lon = np.asscalar(lon2) 
        lat = np.asscalar(lat2)
        
        # find nearest location in bathymetry grid
        bottom  = -gem.bathyF(lon, lat)
#        idx1 = np.argmin(np.abs(lon - gem.bathy['lon'][:]))
#        idx2 = np.argmin(np.abs(lat - gem.bathy['lat'][:]))
#        bottom = -1*gem.bathy['elevation'][idx2, idx1]
        
        
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
        u_all.append(u)
        v_all.append(v)
        N2_all.append(N2)
        rx_all.append(rx)
        ry_all.append(ry)
        rz_all.append(rz)
        dudz_all.append(dudz)
        dvdz_all.append(dvdz)

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
        
        # Boundary checks on rbf functions

        if lon <= lonlims[0] or lon >= lonlims[1]:
#            lonlims, latlims, tlims = gem.interpfield_lnd(lon, lat, z, t)
            print('WARNING: satGEM field functions must be regenerated')
            break

        
        if lat <= latlims[0] or lat >= latlims[1]:
#            lonlims, latlims, tlims = gem.interpfield_lnd(lon, lat, z, t)
            print('WARNING: satGEM field functions must be regenerated')
            break


        if tnum <= tlims[0] or lat >= tlims[1]:
#            lonlims, latlims, tlims = gem.interpfield_lnd(lon, lat, z, t)
            print('WARNING: satGEM field functions must be regenerated')
            break
        # print(t)


    # After ray tracing loop
    elapsed_time = np.vstack([(timeIn - time[0]).total_seconds() 
                    for timeIn in time[:i+2]])
    
    # Extend Bathymetry data so you can see whats just beyond end of ray trace
    ext = np.arange(0, extend_bathy, 10)
    x_extend = np.full_like(ext, np.nan)
    y_extend = np.full_like(ext, np.nan)
    b_extend = np.full_like(ext, np.nan)
    for i in range(len(ext)):
        x = x + dx
        y = y + dy
        lon2, lat2 = inverse_hav(dx, dy, lon, lat)
        lon = np.asscalar(lon2) 
        lat = np.asscalar(lat2)
        bottom  = -gem.bathyF(lon, lat)
        x_extend[i] = x
        y_extend[i] = y
        b_extend[i] = bottom
        
        

        
    
    # store all results in dictionary (keeps things concise when using)

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
            'elapsed_time': elapsed_time,
            'u': np.vstack(u_all)

    }
    
    distance = 1e-3 * np.sqrt(results['x']**2 + results['y']**2)
    results['distance'] = distance

    # Had to add this condition for when the run quits out on first step
    if bathy:
        results['bathy'] = np.vstack(bathy)
        results['bathy_ext'] = np.hstack((results['bathy'].flatten(), b_extend))
        results['x_ext'] = np.hstack((results['x'].flatten(), x_extend))
        results['y_ext'] = np.hstack((results['y'].flatten(), y_extend))
        results['cgx'] = np.vstack(cgx)
        results['cgy'] = np.vstack(cgy)
        results['cgz'] = np.vstack(cgz)
        results['N2'] = np.vstack(N2_all)
        results['rx'] = np.vstack(rx_all)
        results['ry'] = np.vstack(ry_all)
        results['rz'] = np.vstack(rz_all)
        results['dudz'] = np.vstack(dudz_all)
        results['dvdz'] = np.vstack(dvdz_all)


    if N2subgrid is not None:
        results['N2sub'] = N2subgrid
        results['lonvec'] = lonvec
        results['latvec'] = latvec


#    N2_all = np.vstack(N2_all)
#    N2_grid = np.vstack(N2_grid)
#    
#    results['N2'] = N2_all
#    results['N2_grid'] = N2_grid
    
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


















def testing():
    """
    variables used for testing ray tracing (taken from transect observations)
    
    """
    
    towyo_date = datetime(2012, 2, 28, 21, 33, 44)
    start_time = datetime(2012, 2, 28, 21, 33, 44)
    gem = satGEM_mini()
    lon = -49.24
    lat = -53.56
    k0 = 0.000379449
    omega = -0.000147305 	
    l0 = -0.000395896
    m0 = -0.0062492
    w0 = -0.00014730
    
    
    omega = -0.000147305
    k0 = 0.000379449
    l0 = -2.26563e-05
    m0 = -0.0062492
    z = 1536
    lat = -53.56228373
    lon = -49.23690843
    
    
    wave = Wave(k=k0, l=l0, m=m0, w0=w0, z0=1500, t0=towyo_date,
                lat=lat, lon=lon
                )

    
    ray_tracing_interp(wave, gem, duration=4*24,
                       tstep=5,
                       time_direction='reverse',
                       latpad=5, lonpad=5, tpad=2,
                       extend_bathy=10*1000,
                       interp_mode='lnd', interp_function='linear')
    
    results = wave.results
#    
    plot_results(results, gem)
    
    plt.figure()
    plt.plot(results['elapsed_time'][:-1]/3600, results['u']*results['k'][1:])
    
    plt.figure()
    plt.plot(results['elapsed_time'][:-1]/3600, results['N2'])
    
















def plot_results(results, gem, lc='#ff3f02',
                lw=1.5, ms=20, buffer=0.2, cls=30,
                 testing=False):
    """
    Quick function for plotting a ray path over local bathymetry with
    start and ends marked

    Parameters
    ----------


    Returns
    -------


    """

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


    # Generate transec bathymetry (interpolate to smooth)
    
    plt.subplot(222)
    distance1 = (1e-3 * np.sqrt(results['x_ext']**2 + results['y_ext']**2))
    bathy1, idx1 = np.unique(results['bathy_ext'], return_index=True)
    bathy1 = np.vstack([results['bathy_ext'][index] for index in sorted(idx1)])
    x1 = np.vstack([distance1[index] for index in sorted(idx1)])
    bathy_rev = np.interp(distance1.flatten(),  x1.flatten(), bathy1.flatten())
    plt.plot(distance, results['z'], c=lc, linewidth=lw)
    plt.fill_between(distance1.flatten(), bathy_rev.flatten(), 6000,
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
    plt.plot(results['elapsed_time']/3600, (2*np.pi)/ (np.abs(results['m'])),
                                c=lc, linewidth=lw)
    plt.gca().format_xdata = mdates.DateFormatter('%h')
    plt.xlabel('Time (Hours)')
    plt.ylabel('vertical Wavelength')
    plt.title(' Vertical Wavelength Evolution')

    plt.subplot(224)
    
    plt.plot(results['elapsed_time']/3600, np.log10(results['omega']**2),
                    c=lc, linewidth=lw, label=r'$\omega^2$')
    plt.plot(results['elapsed_time'] / 3600, np.log10(f),
                    c='k', linewidth=lw, label='f^2')

    plt.legend()
    plt.xlabel('Time (Hours)')
    plt.ylabel(r'$log (s^{-2})$')
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

    
