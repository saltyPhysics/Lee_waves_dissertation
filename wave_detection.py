"""

Created January 12th  2018

Attempt at improving wave detection diagnostics from transect of hydrographic sections.

PLAN:
Instead of segmenting casts and choosing integration limits and seeing if things match up, this will treat the transect as a 2D image. A "patch" will be determined by vertical and horizontal distances. Within the path, the PSD at a specfic frequency (inverse wavelength) will fill each cell of the path and try to isolate patches which have a high level of similarity in psd at a given wavelength.


"""

import matplotlib.pyplot as plt 
import numpy as np
import scipy
import gsw
import cmocean