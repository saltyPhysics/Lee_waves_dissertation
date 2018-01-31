"""
Author: Manish Devana
Date: December 31, 2017

Ray Tracing functions for internal waves

The goal is to use an object oriented approach to a ray tracing model.
CURRENT STATUS:
load in satGEM data and rewrite functions to remove all assumptions and run in a 4d field. 

- Figure out how to make k and l and m vary in all 4 dimensions (x,y,z, and t)



"""


import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import gsw
import oceans as oc
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as colors
import cmocean


def instructions():
    """
    Print instructions
    """
    text = '''
Ray Tracing Instructions:
------------------------
1. Generate a "wave" object : rt.wave(inputs)
    - enter wave frequency, horizontal and vertical wavenumbner components, and initial depth
    - view properties to check if things were entered correctly

satGEM Details
--------------
\t This ray tracing model utlizes the 4D velocity and density field constructed by Dr. Andrew Meijers. Full details are available in Meijers 2013. 
    '''

    print(text)


def wavenumbers(lam_k, lam_l, lam_m):
    """
    Returns wavenumbers from wavelengths 
    """
    k = (2*np.pi)/lam_k/1000
    l = (2*np.pi)/lam_l/1000
    m = (2*np.pi)/lam_m

    return k, l, m


def cgz(w0, f, kh, m):
    """
    Vertical Group Speed
    """
    return np.squeeze((((w0**2 - f**2))/(w0*(kh**2 + m**2)))*m)


def cgx(N2, w0, k, kh, m):
    """
    Horizontal group speed in x-direction
    """
    return np.squeeze(((N2 - w0**2)/(w0*(kh**2 + m**2)))*k)


def cgy(N2, w0, l, kh, m):
    """
    Horizontal group speed in y-direction
    """
    return np.squeeze(((N2 - w0**2)/(w0*(kh**2 + m**2)))*l)


def EoZ(N2, w0, f, ):
    """
    Wave ray energy when variations can only occur in the vertical (i.e. N2 and
    flow only vary with depth not horizontally) - Olbers 1981
    """
    Ez = np.squeeze((w0**2 * (N2 - f**2))
                    / ((w0**2 - f**2)**(3 / 2) * (N2 - w0**2)**(1 / 2)))
    return Ez


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class wave(object):
    """
    Creates a wave which has varying functionality including:
    - time forward modelling
    - time reverse modelling
    - variable velocity and density field inputs
    - plotting and saving rays
    - HOPEFULLY: built in gui
    """

    # Add functionality for a default Buoyancy Frequncy and Velocity Profile

    def __init__(self, k=10*1000, l=10*1000, m=500, w0=8e-4, z0=500):

        # Convert wavelengths into wavenumbers
        # Save initial values becuase running the model will change
        # the wave features.
        self.k_init = np.array([k], dtype='float')
        self.l_init = np.array([l], dtype='float')
        self.m_init = np.array([m], dtype='float')
        self.w0_init = np.array([w0], dtype='float')
        self.kh_init = np.array([np.sqrt(self.k_init**2 + self.l_init**2)])
        self.x_init = np.array([0], dtype='float')
        self.y_init = np.array([0], dtype='float')
        self.z_init = np.array([z0], dtype='float')

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
'''.format(np.array2string(self.k_init), self.l_init, self.m_init, self.kh_init, self.w0_init)

        print(txt)
  
    def back_trace(self, duration=24, tstep=5,status=2,
                   seafloor=4000, print_run_report=False, updates=False):
        """
        3 dimensional ray tracing within the time evolving satGEM density and velocity fields


        Parameters
        ----------
        duration:  Run duration (in hours)
        tstep: time step of model in seconds
        status:
        seafloor: choose a seafloor depth in meters (maybe integrate bathymetry at some point using gebco data)
        
        
        Returns
        -------
        """

        # Set up model run
        tstep = float(tstep)
        duration = duration * 60 * 60  # Convert to seconds
        time = np.arange(0, duration, tstep)  # time grid for model run
        status = status*60*60
        f = gsw.f(self.lat)

        # Run model
        while steady_state:

            z = self.z_init[:]
            x = self.x_init[:]
            y = self.y_init[:]
            w0 = self.w0_init[:]
            m = self.m_init[:]
            k = self.k_init[:]
            l = self.l_init[:]
            kh = self.kh_init[:]



            bottom = 'NO'
            for t in range(len(time)):
                idx = np.nanargmin(np.abs(z - self.flow_grid))
                idx2 = np.nanargmin(np.abs(z - self.N2_grid))

                # X group step
                x = x - (cgx(self.N2[idx2], w0, k,\
                              kh, m) + self.U[idx])*tstep



                if not np.isfinite(x):
                    error = self.model_error_message(x, y, z, m, idx, idx2)
                    print(' X Error' + error)
                    break
                # Y group step
                y = y - (cgy(self.N2[idx2], w0, l,\
                              kh, m)+ self.V[idx])*tstep

                if not np.isfinite(y):
                    error = self.model_error_message(x, y, z, m, idx, idx2)
                    print(' Y Error' + error)
                    break

                # Z group step
                z = z + cgz(w0, f, kh, m)*tstep

                if not np.isfinite(z):
                    error = self.model_error_message(x, y, z, m, idx, idx2)
                    print(' X Error' + error)
                    break


                if z >= seafloor:
                    print('Wave ray hit seafloor')
                    bottom = 'YES'
                    break

                # vertical wavenumber change
                m = m - -1 * (k * self.dudz[idx] + l * self.dvdz[idx]) * tstep

                # Change in frequency (doppler shifting)
                w0 = w0 - (k * self.U[idx] + l * self.V[idx])

                # Wave Energy(z)
                Ez = EoZ(self.N2[idx2], w0, f)

                # Wave Action
                Ac = Ez / w0

                # Sometimes z turns to nan if the velocity or buoyancy profiles
                # are not complete, right now that breaks the model
                if not np.isfinite(z):
                    print('Error')
                    break

                # Store steps
                self.x_all.append(x)
                self.y_all.append(y)
                self.z_all.append(z)
                self.m_all.append(m)
                self.w0_all.append(w0)
                self.E_all.append(Ez)
                self.Ac_all.append(Ac)


                if time[t] % status == 0 and updates:
                    print("Progress: {} %".format(int(100
                                                      * time[t] / duration)))


            self.x_ray = np.vstack(self.x_all) / 1000  # Convert to Kilometers
            self.y_ray = np.vstack(self.y_all) / 1000  # Convert to Kilometers

            self.z_ray = np.vstack(self.z_all)
            self.m_ray = np.vstack(self.m_all)
            self.w0_ray = np.vstack(self.w0_all)
            self.E_ray = np.vstack(self.E_all)
            self.Ac_ray = np.vstack(self.Ac_all)
            self.bottom = bottom

            if updates:
                print('Run Complete!')
            Run_report = '''
    Ray Tracing Report:
    -------------------
    X-Distance: {} km
    Y-Distance: {} km
    Z-Distance: {} m
    Duration: {} seconds
    Time Step: {} seconds
    Hit Sea Floor : {}
    '''.format(self.x_ray[-1], self.y_ray[-1],
    self.z_ray[-1], duration, tstep, bottom)
            self.run_report = Run_report

            if print_run_report:
                print(Run_report)
            steady_state = False
            self.seafloor = seafloor





        
        
        



def testing():
    """
    This is where I load parameters to test functions
    """

    # These are parameters taken from the lee wave dissertation project
    l1 =  -0.000139543
    k1 = 0.000324448
    m1 = -0.00976433
    z0 = 1500
    w0 = -0.000132755
    wave1 = wave(k=k1, l=l1, m=m1, w0=w0, z0=z0)


    # run models
    duration = 48
    tstep = 10
    status = 6 # intervals to give run status
    wave1.backward3d(duration=duration, tstep=tstep,
                     status=status, print_run_report=True)


    plt.figure()
    plt.plot(wave1.x_ray, wave1.z_ray)
    plt.gca().invert_yaxis()
