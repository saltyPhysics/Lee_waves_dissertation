"""
Created on January 11th 2018. 

@author: manishdevana
Analyzing and cleaning sat gem data - from Andrew Meijers (BAS)


"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


mat = sio.loadmat('DIMES_vel_09_12_upd.mat', squeeze_me=True)

