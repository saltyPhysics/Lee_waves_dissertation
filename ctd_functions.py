#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 13:40:35 2017

@author: manishdevana
function set for ctd microstructure processing. 
"""

import numpy as np
from scipy import signal, stats


def data_binner(data_dict, p_grid, bin_size=256):
    """ 
    This function bins data from top down in specified bin size and returns the
    binned pressure grid associated.
    """
    binned_p = []
    binned_data = {}
    for data_key in data_dict.keys():
        data_in = data_dict[data_key]
        binned_casts = []
        for cast in data_in.T:
            count2 = 0
            if count2 <= 1:
                # Creates associated pressure grid
                cutoff = len(data_in) % bin_size
                usable_range = list(range(cutoff, len(data_in)+1, int(.5*bin_size)))
                for i in usable_range[:-2]:
                    binned_p.append(p_grid[i:i+bin_size])

            # Bins data and detrends using linear regression and p-grid
            for i in usable_range[:-2]:
                binned = cast[i:i+bin_size]
                nn_idx = ~np.isnan(binned)
                if np.sum(nn_idx) <= .5*bin_size:
                    continue
                m, b, r_val, p_val, std_err = stats.linregress(binned_p[count2][nn_idx], binned[nn_idx])
                detrend_cast = binned - (m*binned_p[count2] + b)
            binned_casts.append(detrend_cast)
            count2 += 1
        binned_data[data_key] = binned_casts
#    binned_data['p'] = binned_p
    return binned_data
