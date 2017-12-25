#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 21:02:14 2017

@author: manishdevana
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import gsw
from scipy import signal, stats
import ctd_functions as ctdf
import mat4py as m4p
from itertools import chain

def load_data():
    fname = 'JC69_Comb.mat'
    
    all_data = sio.loadmat(fname)
    ctd = all_data['ctd']
    ctd_dict = {}
    ctd_dict['p'] = ctd['p'][0][0]
    ctd_dict['s'] = ctd['s'][0][0]
    ctd_dict['t'] = ctd['t'][0][0]
    ctd_dict['pt'] = ctd['pt'][0][0]
    ctd_dict['lat'] = ctd['lat'][0][0]
    ctd_dict['lon'] = ctd['lon'][0][0]
    ctd_dict['ctd'] = ctd['ctd'][0][0]
    
    # When extracting ladcp data I found a much more efficient way by importing a 
    # all the field names from matlab. go back and change it for ctd data later
    ladcp = all_data['ladcp']
    ladcp_names = m4p.loadmat('ladcp_names.mat')
    ladcp_names = ladcp_names['names']
    ladcp_names = list(chain.from_iterable(ladcp_names))
    ladcp_dict = {name:ladcp[name][0][0] for name in ladcp_names}
    
    
    # Previously determined in matlab and verified with alberto. the range is 40-63 
    # since python range doesnt include the last value. Using the loop below, The 
    # data for the desired stations are extracted
    
    station_idx = list(range(41,62)) 
    
    for key, value in ctd_dict.items():
        unsorted = ctd_dict[key]
        ctd_dict[key] = unsorted[:,station_idx]
        
    ladcp_grids = {}
    for key, value in ladcp_dict.items():
        unsorted = ladcp_dict[key]
        if unsorted.shape == (100,1000):
            ladcp_grids[key] = unsorted.T[:,station_idx]
        elif unsorted.shape == (1, 100):
            ladcp_grids[key] = unsorted.T[station_idx]
        elif unsorted.shape == (1,1000):
            ladcp_grids[key] = unsorted.T
    
    
    
    bathy = -1*(np.loadtxt('bathy_transect.txt'))
    
    return ladcp_grids, ctd_dict, bathy

            

""" Use GSW- python package for density (and other hydrography) calculations. The 
link can be found here for the docs ---> https://teos-10.github.io/GSW-Python/density.html """



""" Binning CTD data for spectral calculations
        This loop bins data according to the bin_size in half-overlapping segments
        descending. In order to deal with the length of the casts not matching 
        and integer number of bins, the remainder of the cast length and bin size
        subtracted out and trimmed off the top. Each binned set is detrended
        before being stored. ***I'm not sure if this is the right
        way to deal with this*** 
        
        For CTD - bin size = 256 (half overlapping bins)
        --> this means 22 bins (256 points in each bin) in each cast with 
        total 23 casts.
        
        Try to rewrite the binning as a function 
        
        """ 


