"""
Created on January 11th 2018. 

@author: manishdevana
Analyzing and cleaning sat gem data - from Andrew Meijers (BAS)
 - For loading data use hdf python library... works SOOO MUCH BETTER

"""

import matplotlib.pyplot as plt
import numpy as np
import h5py



def load_data():
    """ 
    Load satGEM and transect data
    """
    
    # load transect data
    jc69 = h5py.File('JC69_Devana.mat')
    
  
    # Keep data in h5py file objects otherwise the whole thing goes to hell
    gamma_file = h5py.File('DIMES_GAMMA_09_12_upd.mat')
    gamma_variables = list(gamma_file.keys())
    vel_file =  h5py.File('DIMES_vel_09_12_upd.mat')
    vel_variables = list(vel_file.keys())
    TS_file =  h5py.File('DIMES_TS_09_12_upd.mat')
    TS_variables = list(TS_file.keys())
    
    
    
    # THESE ARE NOT NUMPY ARRAYS THEY ARE H5PY OBJECT BUT BEHAVE SIMILARLY
    # if there needs t be
    gamma = gamma_file['satGEM_gamma']
    u = vel_file['satGEM_east']
    v = vel_file['satGEM_north']
    temp = TS_file['satGEM_temp']
    sal = TS_file['satGEM_sal']
    
    # Data grids
    time = gamma_file['time']
    depth_grid = gamma_file['depthlvl']
    lon = gamma_file['lons']
    lat = gamma_file['lats']
    u_lat = vel_file['centerlat']
    v_lon = vel_file['centerlon']
    
    # data is structured in the shape [380, 401, 60, 180]
    # which equals [lon, lat, depth_grid, time]
    
    






