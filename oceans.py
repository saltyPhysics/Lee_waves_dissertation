#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:11:30 2017

General Oceanographic tools


NOTES: geostrophy is wrong.. need to go back and try with geopotentials?


@author: manishdevana
"""
import numpy as np
import pandas as pd
import scipy.interpolate as interp
import scipy
import matplotlib.pyplot as plt
import gsw
import cmocean
from netCDF4 import Dataset
from matplotlib.collections import LineCollection

default_params = {
    'bin_size': 200,
    'overlap': 100,
    'm_0': 1. / 300.,  # 1/lamda limits
    'm_c': 1. / 10.,  # 1/lamda limits
    'order': 1,
    'nfft': 256,
    'plot_data': 'on',
    'transect_fig_size': (6, 4),
    'reference pressure': 0,
    'plot_check_adiabatic': False,
    'plot_spectrums': False,
        }

def loadCTD(ctd):
    """
    Loads standard ctd data from a dictionary of ctd values
    """

    S = ctd['s']
    T = ctd['t']
    p = ctd['p']
    lat = ctd['lat']
    lon = ctd['lon']


    return S, T, p, lat, lon

def loadLADCP(ladcp, full_set=False):
    """
    Loads standard ctd data from a dictionary of ctd values
    """

    U = ladcp['u']
    V = ladcp['v']
    p = ladcp['p_grid']
    uup = ladcp['uup']
    vup = ladcp['vup']
    udo = ladcp['udo']
    vdo = ladcp['vdo']





    if full_set:
        return U, V, p, uup, vup, udo, vdo
    else:
        return U, V, p

def depthTrim(data, z, maxDepth):
    """
    Trims array to a max depth with corresponding depth or pressure grid
    """
    idx = np.squeeze(z <= maxDepth)
    dataTrimmed = data[idx,:]

    return dataTrimmed, z[idx]

def rhoFromCTD(S, T, p, lon, lat):
    """
    Rho from measured salinity and temperature values
    """

    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_t(SA, T, p)
    rho = gsw.density.rho(SA, CT, p)

    return rho

def sigma0FromCTD(S, T, p, lon, lat):
    """
    Rho from measured salinity and temperature values
    """

    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_t(SA, T, p)
    rho = gsw.density.sigma0(SA, CT)

    return rho



def gswN2(S, T, p, lat, lon, axis=0):


    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_t(SA, T, p)

    N2, dump = gsw.stability.Nsquared(SA, CT, p, lat, axis=axis)
    N2poly = []
#    for cast in N2.T:
#        fitrev = vert_polyFit(cast, p[:,1], 100, deg=1)
#        N2poly.append(fitrev)
#
#    N2poly = np.vstack(N2poly).T
    bottom_row = np.full((1,21), np.nan)
    N2 = np.vstack((N2, bottom_row))


    return N2


def bfrq(T, S, z, lat, lon):
    """
    Calculate Bouyancy frequency from practical salinity and temperature.
    Converts to potential temperature and absolute salinity to calculate density
    """


    SA = gsw.conversions.SA_from_SP(S, z, lon, lat)
    g = gsw.grav(lat, z)
    pdens = gsw.rho_t_exact(SA, T, z)
    dpdens = np.gradient(gsw.rho_t_exact(SA, T, z), axis=0)
    dz = np.gradient(z, axis=0)
    N2 = (g/pdens)*(dpdens/dz)

    return N2

def sliding_bins(data, p, bin_size, step=100):
    """
    Bin data in sliding vertical bins with a specified overlap
    """

    dz = np.nanmean(np.gradient(p))

    bin_step = int(np.ceil(bin_size/dz))
    step = int(np.ceil(step/dz))

    # start point of each bin
    starts = np.arange(0, data.shape[0], step)

    idx = []
    for start in starts:
        if start+bin_step <= data.shape[0]:
            idx.append(np.arange(start, start+bin_step, dtype=int))

    return np.vstack(idx)

def binData(data, p, bin_size):
    """
    This function takes an array of hydrographic casts (each column representing
    an indivdual cast) and bins the data into half overlapping bins, returning
    only the indices i.e. no modifications are made to the original data
    """

    # Pressure steps
    dp = np.nanmean(np.gradient(p))

    steps = int(np.ceil(bin_size/dp))

    starts = np.arange(0, data.shape[0], steps/2)

    idx = []
    for i in starts:
        if i+steps <= data.shape[0]:
            idx.append(np.arange(i,i+steps, dtype=int))

    idx = np.vstack(idx)

    return idx





def depthAvgFlow(U, V, dz, depthRange=500):
    """
    Depth averaged flow for flow from bottom up to depthrange.
    """
    steps = int(np.ceil(depthRange/dz))
    Umeans = np.empty(0)
    Vmeans = np.empty(0)

    for Ui, Vi in zip(U.T, V.T):
        mask = np.where(np.isfinite(Ui))[0]
        mask = mask[-1:np.max(mask)-steps:-1]
        Umeans = np.append(Umeans, np.nanmean(Ui[mask]))
        Vmeans = np.append(Vmeans, np.nanmean(Vi[mask]))


    return Umeans, Vmeans


def bathyLoadNc(fname, lat, lon, add_buffer=True, buffer=.1):
    """
    Loads bathymetric data from Netcdf file and extracts box around transect
    """
    file = Dataset(fname)

    bathyDict = {key:file.variables[key][:] for key in file.variables.keys()}
    file.close()
    lat2 = bathyDict['lat']
    lon2 = bathyDict['lon']
    if add_buffer:
        # go back and actually add the buffer in
        latidx = np.where(np.logical_and(lat2 < np.nanmax(lat)+buffer\
                                , lat2 > np.nanmin(lat)-buffer))[0]

        lonidx = np.where(np.logical_and(lon2 < np.nanmax(lon)+buffer\
                                , lon2 > np.nanmin(lon)-buffer))[0]
    else:
        latidx = np.where(np.logical_and(lat2 < np.nanmax(lat)\
                                , lat2 > np.nanmin(lat)))[0]

        lonidx = np.where(np.logical_and(lon2 < np.nanmax(lon)\
                                , lon2 > np.nanmin(lon)))[0]

    bathyrev = bathyDict['elevation'][latidx,:]
    bathyrev = bathyrev[:,lonidx]
    longrid, latgrid = np.meshgrid(lon2[lonidx], lat2[latidx])



    return bathyrev, longrid, latgrid


def matchGrids(data, oldGrid, newGrid, interp_method='pchip'):
    """
    General function for interpolating one grid of hydographic data onto a
    desired grid

    Parameters
    ----------
    data : the data to be gridded onto a new grid
    griddata : the datagrid to interpolate onto

    Returns
    -------
    gridded_data : the data gridded onto chosen grid

    """

    newGrid = np.squeeze(newGrid)
    oldGrid = np.squeeze(oldGrid)

    # find which axis to interpolate along (and which counts the stations)
    data_dims = np.array(data.shape)
    axis_choose = int(np.where(data_dims != len(oldGrid))[0])

    gridded_data = np.empty((len(newGrid), data.shape[axis_choose]))
    # Interpolate vertically through each cast
    for i in range(data.shape[axis_choose]):
        f = interp.PchipInterpolator(oldGrid, data[:,i])
        gridded_data[:,i] = f(newGrid)


    return gridded_data

def verticalBoxFilter1(data, z, box=100):
    """
    Function for filtering a vertical cast of data using a box average method
    with fixed box size

    """
    dz = np.nanmean(np.gradient(np.squeeze(z)))
    window = int(np.ceil(box/dz))
    filtered = []
    for i in range(len(data)):
        if i < len(data) - window:
            filtered.append(np.nanmean(data[i:i+window]))
        else:
            filtered.append(np.nanmean(data[i:-1]))

    filtered = np.hstack(filtered)

    return filtered


def vert_polyFit(data, z, bin0, step=1, deg=2):
    """
    Trying to use the vertical polynomial fit to clean up the data
    not reallly sure about what im doing though
    """

    data = np.squeeze(data)
    z = np.squeeze(z)

    dz = np.nanmean(np.gradient(np.squeeze(z)))
    bin1 = int(np.ceil(bin0/dz))

    fits = []
    zFits = []
    bins = []
    for i in range(len(z)):

        if 2*i+bin1 < len(z):
            bins.append(np.arange(i,2*i+bin1+1))
            mask = np.isfinite(data[i:2*i+bin1])
            dataIn = data[i:2*i+bin1]
            zIn = z[i:2*i+bin1]
            dataIn = dataIn[mask]
            if dataIn.size == 0:
                fits.append(np.nan)
                zFits.append(np.nan)
            else:
                zIn = zIn[mask]
                zFits.append(np.nanmean(z[i:2*i+bin1]))
                p = scipy.polyfit(zIn, dataIn, deg)
                fits.append(np.nanmean(scipy.polyval(p, z[i:2*i+bin1])))

    fits = np.hstack(fits)
    zFits = np.hstack(zFits)
    P = scipy.interpolate.interp1d(zFits, fits, fill_value='extrapolate')
    fitrev = P(z)


    return fitrev


def vert_polyFit2(data, z, bin0, step=1, deg=2):
    """
    Trying to use the vertical polynomial fit to clean up the data
    not reallly sure about what im doing though
    """

    data = np.squeeze(data)
    z = np.squeeze(z)

    dz = np.nanmean(np.gradient(np.squeeze(z)))
    bin1 = int(np.ceil(bin0/dz))

    fits = []
    zFits = []
    bins = []
    for i in range(len(z)):

        if 2*i+bin1 < len(z):
            bins.append(np.arange(i,2*i+bin1+1))
            mask = np.isfinite(data[i:2*i+bin1])
            dataIn = data[i:2*i+bin1]
            zIn = z[i:2*i+bin1]
            dataIn = dataIn[mask]
            if dataIn.size == 0:
                fits.append(np.nan)
                zFits.append(np.nan)
            else:
                zIn = zIn[mask]
                zFits.append(np.nanmean(z[i:2*i+bin1]))
                p = scipy.polyfit(zIn, dataIn, deg)
                fits.append(np.nanmean(scipy.polyval(p, z[i:2*i+bin1][mask])))

    fits = np.hstack(fits)
    zFits = np.hstack(zFits)
    mask2 = np.isfinite(fits)
    P = scipy.interpolate.interp1d(zFits[mask2], fits[mask2], fill_value='extrapolate')

    fitrev = P(z)


    return fitrev


def velocityFiltered(U, V, p_ladcp, deg=1):

    Upoly = []

    for cast in U.T:
        fitrev = vert_polyFit(cast, p_ladcp, 100, deg=deg)
        Upoly.append(fitrev)

    Upoly = np.vstack(Upoly).T
    U = U - Upoly

    Vpoly = []

    for cast in V.T:
        fitrev = vert_polyFit(cast, p_ladcp, 100, deg=deg)
        Vpoly.append(fitrev)

    Vpoly = np.vstack(Vpoly).T
    V = V - Vpoly

    return U, V


def SpectrumGenerator_vertical(data, dz, length_sig, binned=False,\
                               detrend=True,params=default_params):
    """
    Function for getting the spectrum and associated wavenumber and wavelength
    axes
    """

    nfft = int(2**np.ceil(np.log2(length_sig)))

    # Sampling frequency
    fs = 1./dz
    Nc = .5*fs

    mx = Nc*(np.arange(0, nfft))/(nfft/2)


    spectrum = []
    #Assumes binned data
    if binned:
        for station in data:
            spectrum.append(scipy.fftpack.fft(station, n=nfft, axis=1))

    else:
        #Assume gridded profiles with each column being a new cast
        for i in range(data.shape[1]):
            mask = np.isfinite(data[:,i])
            if detrend:
                data[mask,i] = scipy.signal.detrend(data[mask,i])
                data[mask,i] = data[mask,i]-np.nanmean(data[mask,i])
            spectrum.append(scipy.fftpack.fft(data[mask,i], n=nfft))
        spectrum = np.vstack(spectrum).T



    # Convert 1/lambda to k (wavenumber)
    kx = 2*np.pi*mx

    return spectrum, mx, kx


def verticalBandPass(data, z, m_min, m_max, return_sig=True):
    """
    filter vertical profiles using band pass filters to look for wavelike wiggles
    """

    # invert wavelengths since fft returns on frequency grid (1/lamda)
    m1 = 1/m_min
    m2 = 1/m_max


    # get spectra of each vertical cast
    dz = np.nanmean(np.gradient(np.squeeze(z)))
    spectra, mx, kx = SpectrumGenerator_vertical(data, dz, data.shape[0])

    # Normalize Power
    power = np.abs(scipy.fftpack.fftshift(spectra, axes=0))
    power = power/len(mx)


    # Filter on shifted spectrum
    midpoint = int(len(mx)/2)
    pos_half = mx[1:midpoint+1]
    neg_half = np.flipud(-pos_half)
    mxShift = np.hstack((neg_half, pos_half))

    mask1 = np.logical_and(np.abs(mxShift)>=m2, np.abs(mxShift)<=m1)
    bandpass1 = np.logical_not(mask1)

    filtShift = scipy.fftpack.fftshift(spectra)
    filtShift[bandpass1,:] = 0
    powerFilt = np.abs(filtShift)
    powerFilt = 2*powerFilt/len(mx)



    # create band bass filters using min and max lamdas
    mask = np.logical_and(mx>=m2, mx<=m1)
    bandpass = np.logical_not(mask)

    # Apply filter be turning all non desired values to zero

    filtered = spectra[:]
    filtered[bandpass,:] = 0


    # shift mx grid
    midpoint = int(len(mx)/2)
    pos_half = mx[1:midpoint+1]
    neg_half = np.flipud(-pos_half)
    mxShift = np.hstack((neg_half, pos_half))


    # retur wavnumber and wavelength grids along with the spectra and filter
    return mx, kx, spectra, bandpass, filtShift, power, mxShift, powerFilt



def velocityInLine(U, V, lat, lon):
    """
    Returns velocity in line with transect
    """

    # Angle between velocity vector and transect
#    theta = np.arctan((lat[-1]-lat[0])/(lon[-1]-lon[0])) - \
#                    np.arctan(U/V)

    theta = np.full_like(U, np.nan)
    for i in range(len(lat)):
        if i != 0:
            theta[:,i] = np.arctan((lat[i]-lat[i-1])/(lon[i]-lon[i-1])) - \
                                np.arctan(U[:,i]/V[:,i])
        else:
            theta[:,i] = np.arctan((lat[i]-lat[i-1])/(lon[i]-lon[i-1])) - \
                                np.arctan(U[:,i]/V[:,i])

    # Velocity Magnitude
    Umag = np.sqrt(U**2 + V**2)


    # Velocity in direction of transect
    X = Umag*np.cos(theta)

    # U and V breakdown of X in direction of transect line
    phi = np.full_like(U, np.nan)
    for i in range(len(lat)):
        phi[:,i]  = np.arctan((lat[i]-lat[i-1])/(lon[i]-lon[i-1]))

    Ux = np.cos(phi)*Umag
    Vx = np.sin(phi)*Umag

    return X, Ux, Vx


def speedAtZ(U, V, z, depth, bin_width=100):
    """
    Extracts mean velocity at depth with average of bin width
    """

    z_upper = np.nanargmin(np.abs(z - (depth+.5*bin_width)))
    z_lower = np.nanargmin(np.abs(z - (depth-.5*bin_width)))
    zmean = np.nanmean(z[z_lower:z_upper])

    Urev = U[z_lower:z_upper,:]
    Urev = np.nanmean(Urev, axis=0)
    Vrev = V[z_lower:z_upper,:]
    Vrev = np.nanmean(Vrev, axis=0)

    return Urev, Vrev


def adiabatic_level(S, T, p, lat, pref=0,
                    pressure_range=400,
                    order=1, axis=0,
                    eta_calc=True):
    """
    Adiabatic Leveling from Bray and Fofonoff (1981) - or at least my best
    attempt at this procedure.

    -------------
    Data should be loaded in bins - might try to add a way to check that later

    -------------
    CTD and ladcp data are on two different grids so this function uses the finer
    CTD data to generate the regressions which are then used on the pressure grid

    PARAMETERS
    ----------
    S: Salinity
    T: Temperature
    p: Pressure grid
    lat: latitude grid
    pref: reference pressure
    pressure_range: window of pressure for adiabtic Leveling
    order: polynomial order for fitting
    axis: axis to perform calculations on
    eta_calc: If true, calculates eta via (rho - rho_ref)/(drho_ref/dz) in same
            windoes as adiabatic leveling calculations


    Returns
    -------
    N2: Bouyancy Frequency Squared
    N2_ref: Buoyancy Frequency Squared reference field from adiabtic Leveling
    p_mid: midpoint pressure grid from N2 calculations
    strain: strain calculated using N2 (vertical gradient of vertical isopycnal displacement)


    """

    # Pressure window - See Bray and Fofonoff 1981 for details
    plev = pressure_range

    # Calculate Buoyancy Frequency using GSW toolbox
    N2, p_mid = gsw.stability.Nsquared(S, T, p, lat, axis=axis)

    # Keeps as rank 2 array making it compatible with casts and arrays
    if len(N2.shape) < 2:
        N2 = np.expand_dims(N2, 1)
        p_mid = np.expand_dims(p_mid, 1)
        S = np.expand_dims(S, 1)
        T = np.expand_dims(T, 1)
        p = np.expand_dims(p, 1)

    # Calculate reference N2 Field
    alpha = np.full((p_mid.shape[0], p_mid.shape[1]), np.nan)
    g = np.full((p_mid.shape[0], p_mid.shape[1]), np.nan)
    rho_bar = np.full((p_mid.shape[0], p_mid.shape[1]), np.nan)
    drho_dz = np.full((p_mid.shape[0], p_mid.shape[1]), np.nan)


#    N2_ref = np.full_like(N2, np.nan)
    for i in range(p_mid.shape[axis]):
        # Bottom and top for window - the max and min of the 2-element array
        # allows windows to run to ends of profile (I think)

        for k in range(p_mid.shape[1]):
            p_min = np.max([p_mid[i,k] - 0.5*plev, p[0,k]])
            p_max = np.min([p_mid[i,k] + 0.5*plev, p[-1,k]])
            data_in = np.logical_and(p >= p_min, p<= p_max)

            if np.sum(data_in) > 0:
                p_bar = np.nanmean(p[data_in[:,k],k])
                t_bar = np.nanmean(T[data_in[:,k],k])
                s_bar = np.nanmean(S[data_in[:,k],k])
                rho_bar[i,k] = gsw.density.rho(s_bar, t_bar, p_bar)

                # Potential Temperature referenced to pbar
                theta = gsw.conversions.pt_from_t(S[data_in[:,k],k],
                                                  T[data_in[:,k],k],
                                                  p[data_in[:,k],k],
                                                  p_bar)

                sv = 1. / gsw.pot_rho_t_exact(S[data_in[:,k],k],
                                              theta,
                                              p[data_in[:,k],k],
                                              p_bar)

                # Regress pressure onto de-meaned specific volume and store coefficients
                Poly = np.polyfit(p[data_in[:, k], k], sv - np.nanmean(sv), order)

                # Do something with coefficients that I dont understand
                alpha[i, k] = Poly[order - 1]

                # calculate reference N2 reference field
                g[i, k] = gsw.grav(lat.T[k], p_bar)

    # Calculate N2 grid (10^-4 to convert from db to pascals)
    N2_ref = -1e-4 * rho_bar**2 * g**2 * alpha

    # strain calcuations
    strain = (N2 - N2_ref) / N2_ref

    # Append row of nans onto the top or bottom because computing N2 chops a row
    p_mid = np.vstack((np.zeros((1, p_mid.shape[1])), p_mid))
    rho_bar = np.vstack((np.full((1, p_mid.shape[1]), np.nan), rho_bar))

    return N2_ref, N2, strain, p_mid, rho_bar


def isopycnal_displacements(rho, rho_ref, p, lat, diff_window=400, axis=0):
    """
    Compute vertical isopycnal displacements
    (eta) as (rho - rho_ref)/(drho_ref/dz) with drho_ref/dz computed over a
    vertical window (diff_window)

    Parameters
    ----------
    rho: density (nd-array)
    rho_ref: reference density (see adiabatic leveling for computing ref rho)
    p: pressure array
    diff_window: vertical window for computing drho_ref/dz
    axis: axis of vertical direction

    Returns
    -------
    eta: vertical isopycnal displacements


    """

    win = diff_window
    flip = False

    # makes the vertical be along dimension 0 if not already
    if axis == 1:
        rho = rho.T
        rho_ref = rho_ref.T
        p = p.T
        flip = True

    elif axis > 1:
        raise ValueError('Only handles 2-D Arrays')

    z = -1 * gsw.z_from_p(p[:, 0], lat[:, 0])
    dz = np.nanmean(np.diff(z))
    step = int(np.floor(.5 * win / dz))
    eta = np.full_like(rho, np.nan)
    for i in range(rho.shape[0]):

        # If in the TOP half of the profile the window needs to be adjusted
        if i - step < 0:
            lower = 0
            upper = int(2 * step)

        # If in the BOTTOM half of the profile the window needs to be adjusted
        elif i + step > (rho.shape[0] - 1):
            lower = int(rho.shape[0] - 2*step)
            upper = -1

        else:
            upper = i + step
            lower = i - step

        drefdz = (rho_ref[upper] - rho_ref[lower])/win

        eta[i, :] = (rho[i, :] - rho_ref[i]) / drefdz

    return eta


def thorpe_scales(S, T, p, lat, lon, axis=-1):
    """
    Thorpe scales simplified for estimating dissipation rates

    Parameters
    ----------
    S : Practical Salinity
    T : Practical Temperature
    P : Pressure Measurements -- NOT A NORMALIZED PRESSURE GRID
    lat : Latitude
    lon : Longitude

    Returns
    -------
    thorpe: Thorpe Scale estimated from vertical displacement RMS

    """
    rho = rhoFromCTD(S, T, p, lat, lon)
    order = np.argsort(rho, axis=axis)

    displacements = []
    for i, p_in in enumerate(p):
        displacements.append(p[order[i]] - p)

    displacements = np.vstack(displacements)
    thorpe =  np.sqrt(np.mean(displacements**2))


    return thorpe


def geoflow(S, T, p, lat, lon, plot=False):
    """
    Function for calculating the geostrophic flow for an array of ctd data
    using geopotential heights


    Parameters:
        S : Practical Salinity
        T : Practical Temperature
        P : Pressure Measurements -- NOT A NORMALIZED PRESSURE GRID
        lat : Latitude
        lon : Longitude
        plot : switch to run test plots of velocity ---> Default=False

    Returns:
        U : Geostrophic flow in the U Component (Perpendicular to transect)
        V : Geostrophic flow in the V Component (Perpendicular to transect)

    """


    # More precise gravity calculations and density calculations
    g = gsw.grav(lat, p)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_t(SA, T, p)
    rho = gsw.rho(SA, CT, p)

    # Specific Volume
    sv = gsw.density.specvol_anom_standard(SA, CT, p)
    f = gsw.geostrophy.f(lat)

    # calculate distance traveled in meters
    dist = gsw.distance(lon, lat)
    dist = np.cumsum(dist)/1000
    dist = np.append(0,dist)
    dist = dist*1000

    dist = np.tile(dist, (np.shape(rho)[0], 1))
    f = np.tile(f, (np.shape(rho)[0], 1))

    # Geopotential Anomaly Calculation
    phiSum = np.nancumsum(sv, axis=0)

    # loop over grids to make calculations for vertical Shear
    v = np.full((sv.shape[0], sv.shape[1]-1), np.nan)
    for i in range(sv.shape[1]-1):
        for m in range(sv.shape[0]-1):
            dphiB = np.trapz([sv[m+1,i+1]*p[m+1,i+1], sv[m,i+1]*p[m,i+1]])
            dphiA = np.trapz([sv[m+1,i]*p[m+1,i],  sv[m,i]*p[m,i]])
            slope = (dphiB - dphiA)/dist[m,i+1]
            v[m, i] = (1/f[m,i])*slope

    # Integrate vertically to calculate absolute velocity (assume no flow at
    # and integrate upwards)
    geoMag = np.flipud(np.nancumsum(np.flipud(v), axis=0))

    # Calculate U and V components of geostrophic flow
    dlat = np.diff(np.squeeze(lat))
    dlon = np.diff(np.squeeze(lon))
    theta = np.arctan(dlat/dlon)

    # U and V components
    U = np.vstack([-geoMag[:,i]*np.sin(theta[i]) for i in range(len(theta))]).T
    V = np.vstack([geoMag[:,i]*np.cos(theta[i]) for i in range(len(theta))]).T

    # revised distance with dropped first column (distance = 0) to match size
    # of Velocity array
    distRev = dist[:,1:]

    return U, V, geoMag, distRev


def VectorOverheadPlot(U, V, lat, lon, z, depth,\
                       bathy_file, average_window=100,\
                       x_scale=.3, y_scale=.85):
    """
    Function for plotting the flow vectors at a target depth averaged in a
    window (default 100 meters)
    """
    X = x_scale
    Y = y_scale
    # Load bathymetry data from gebco netCDF4 file
    bathy, longrid, latgrid = bathyLoadNc(bathy_file,\
                                          lat,\
                                          lon,\
                                          add_buffer=True)
    #Take Mean flow at chose depth
    U, V = speedAtZ(U, V, z, depth)

    fig = plt.figure()
    plt.pcolormesh(longrid, latgrid,\
                   bathy, cmap=cmocean.cm.deep_r,\
                   shading='gouraud')
    plt.colorbar(label='Depth (meters)')
    plt.plot(lon, lat)
    q1 = plt.quiver(lon, lat,\
                    U, V,\
                    scale=2,\
                    color='red',\
                    width=.004)
    plt.quiverkey(q1, X, Y, .5, "1 m/s", labelpos='W')
    plt.title("Flow at " + str(depth) + "m")

    return fig

def transect_flow_vector(U, V, lat, lon,\
                       bathy_file, title, arrow_scale=1.5, x_scale=.3, y_scale=.85):
    """
    Function for plotting the flow vectors at a target depth averaged in a
    window (default 100 meters)
    """
    X = x_scale
    Y = y_scale
    # Load bathymetry data from gebco netCDF4 file
    bathy, longrid, latgrid = bathyLoadNc(bathy_file,\
                                          lat,\
                                          lon,\
                                          add_buffer=True)

    fig = plt.figure()
    plt.pcolormesh(longrid, latgrid,\
                   bathy, cmap=cmocean.cm.deep_r,\
                   shading='gouraud')
    plt.colorbar(label='Depth (meters)')
#    plt.scatter(np.squeeze(lon), np.squeeze(lat))
    for i, dump in enumerate(np.squeeze(lat)):
        plt.annotate(str(i+1), (np.squeeze(lon)[i],np.squeeze(lat)[i]))
    q1 = plt.quiver(np.squeeze(lon), np.squeeze(lat),\
                    U, V,\
                    scale=arrow_scale,\
                    color='red',\
                    width=.004)
    plt.quiverkey(q1, X, Y, .5, "1 m/s", labelpos='W')
    plt.title(title)

    return fig

def VectorOverheadPlot_compare(U, V, Ui, Vi, lat, lon, z, depth,\
                       bathy_file, average_window=100,\
                       x_scale=.3, y_scale=.85, scale=1.5):
    """
    Function for plotting the flow vectors at a target depth averaged in a
    window (default 100 meters)
    """
    X = x_scale
    Y = y_scale
    # Load bathymetry data from gebco netCDF4 file
    bathy, longrid, latgrid = bathyLoadNc(bathy_file,\
                                          lat,\
                                          lon,\
                                          add_buffer=True)
    lat = np.squeeze(lat)
    lon = np.squeeze(lon)
    #Take Mean flow at chose depth
    U, V = speedAtZ(U, V, z, depth)
    Ui, Vi = speedAtZ(Ui, Vi, z, depth)

    fig, ax = plt.subplots(1,2,figsize=(10,8))
    c1 = ax[0].pcolormesh(longrid, latgrid,\
                   bathy, cmap=cmocean.cm.deep_r,\
                   shading='gouraud')
    ax[0].plot(lon, lat, color='orange', linewidth=1)
    for i, dump in enumerate(lat):
        ax[0].annotate(str(i+1), (np.squeeze(lon)[i],np.squeeze(lat)[i]))
    q1 = ax[0].quiver(lon, lat,\
                    U, V,\
                    scale=scale,\
                    color='red',\
                    width=.004)
    plt.quiverkey(q1, X, Y, .5, "1 m/s", labelpos='W')
    ax[0].set_title("Flow at " + str(depth) + "m - residual")

    c1 = ax[1].pcolormesh(longrid, latgrid,\
                   bathy, cmap=cmocean.cm.deep_r,\
                   shading='gouraud')
    plt.colorbar(c1, label='Depth (meters)')
    ax[1].plot(lon, lat, color='orange', linewidth=1)
    for i, dump in enumerate(lat):
        ax[1].annotate(str(i+1), (np.squeeze(lon)[i],np.squeeze(lat)[i]))
    q1 = ax[1].quiver(lon, lat,\
                    U, V,\
                    scale=scale,\
                    color='red',\
                    width=.004)
    plt.quiverkey(q1, X, Y, .5, "1 m/s", labelpos='W')
    ax[1].set_title("Flow at " + str(depth) + "m - original")

    return fig

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def display(array, cmap=None, index=None, caption=' ', precision=5):
    """
    Function for quick displaying tables in a jupyter notebook
    """

    disp = pd.DataFrame(array, index=index, columns=np.arange(1, array.shape[1] + 1))
    disp.style.background_gradient(cmap=cmap, axis=1)\
        .set_properties(**{'max-width': '300px', 'font-size': '12pt'})\
        .set_caption('buts')\
        .set_precision(precision)

    return disp

def magnify():
    """
    Magnifier for interactive tables in Pandas for jupyter notebook
    """

    return [dict(selector="th",
                 props=[("font-size", "12pt")]),
            dict(selector="td",
                 props=[('padding', "0em 0em")]),
            dict(selector="th:hover",
                 props=[("font-size", "18pt")]),
            dict(selector="tr:hover td:hover",
                 props=[('max-width', '200px'),
                        ('font-size', '18pt')])
]
