"""
Author: Manish Devana
Date: December 31, 2017

Ray Tracing functions for internal waves

The goal is to use an object oriented approach to a ray tracing model.
"""


import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt




def cgz(w0, f, kh, m):
    """
    Vertical Group Speed
    """
    return ((w0**2 - f**2))/(w0*(kh**2 + m**2))*m

def cgx(N2, w0, k, kh, m):
    """
    Horizontal group speed in x-direction
    """
    return ((N2 - w0**2)/(w0*(kh**2 + m**2)))*k


def cgy(N2, w0, l, kh, m):
    """
    Horizontal group speed in y-direction
    """
    return ((N2 - w0**2)/(w0*(kh**2 + m**2)))*l


class wave(object):
    """
    Creates a wave which has varying functionality including:
    - time forward modelling
    - time reverse modelling
    - variable velocity and density field inputs
    - plotting and saving rays
    - HOPEFULLY: built in gui
    """
    lat = 55.0
    # Add functionality for a default Buoyancy Frequncy and Velocity Profile

    def __init__(self, k=10, l=10, m=500, w0=8e-4):

        # Convert wavelengths into wavenumbers
        self.k = (2*np.pi)/(k*1000)
        self.l = (2*np.pi)/(l*1000)
        self.m = (2*np.pi)/m
        self.w0 = w0
        self.kh = np.sqrt(self.k**2 + self.l**2)

        # Save initial values becuase running the model will change
        # the wave features.
        self.k_init = (2*np.pi)/(k*1000)
        self.l_init = (2*np.pi)/(l*1000)
        self.m_init = (2*np.pi)/m
        self.w0_init = w0
        self.kh_init = np.sqrt(self.k**2 + self.l**2)

    def add_all_fields(self, U, V, N2, rho):
        """
        Allows user to add buoyancy, velocity, and density fields
        """
        self.U = U
        self.V = V
        self.N2 = N2
        self.rho = rho


    def help(self):
        """
        Print instructions on how to use wave class
        """
        text = '''
Instructions for using ray tracing model.
\nGenerate a wave with chosen properties or use the defualt Parameters
'''
        print(text)


    def forward2d(steady_state=True):
        """
        Time forward 2d modelling with option for steady state or time evolving
        fields.
        """


    def back2d():
        """
        Time reverese 2d modelling
        """

    def back3d(self, x0=0, y0=0, z0=500, duration=24, tstep=5,
                    steady_state=True, status=2):
        """
        Time forward 2d modelliing
        """
        duration = duration*60*60 # Convert to seconds
        time = np.arange(0, duration, tstep) # time grid for model run
        status = status*60*60


        if steady_state:
            x_all = []
            y_all = []
            z_all = []
            m_all = []
            w0 = []
            x = x0
            y = y0
            z = z0
            w0 = self.w0_init
            m = self.m_init
            k = self.k_init
            l = self.l_init



            for t in range(len(time)):
                x = cgx(N2, w0, k, kh, m)*tstep
                y = cgy(N2, w0, l, kh, m)*tstep
                z = cgz(w0, f, kh, m)*tstep
                x_all.append(x)
                y_all.append(y)
                z_all.append(z)
                if time[t]%status == 0:
                    print("Progress: {} %".format(int(100*time[t]/duration)))

            self.x = np.vstack(x_all)
            self.y = np.vstack(y_all)
            self.z = np.vstack(z_all)
            print('Run Complete!')

    def save_run(self, fname=None):

        if not fname:
            fname = 'ray_trace.csv'
            np.savetxt(fname)
