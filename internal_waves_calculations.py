"""
Created on December 24th 15:44:35 2017

@author: manishdevana
This toolbox calculates internal wave properties



"""


import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import data_load
import gsw
import oceans as oc

def reset_test():
    """
    This loads the data for testing the functions
    (its also all the data for the project)
    """


    ladcp, ctd, bathy = data_load.load_data()
    rho_neutral =  np.genfromtxt('neutral_rho.csv', delimiter=',')
    strain = np.genfromtxt('strain.csv', delimiter=',')
    N2 = np.genfromtxt('ref_N2.csv', delimiter=',')
    wl_max=1000
    wl_min=300
    ctd_bin_size=1024
    ladcp_bin_size=1024

def PowerDens(data, dz, wlmax, wlmin, axis=0, grid=False, detrend='constant'):
    """
    Using periodogram function to estimate power spectral density

    PARAMETERS:
    data: data array for input (assumes column wise data (axis=0))
    dz: vertical sampling rate
    wlmax: maximimum vertical wavelength integration limit
    wlmin: minimum vertical wavelength integration limit

    RETURN:
    variance: The integrated psd of the profiles


    """

    mgrid, Psd = sig.periodogram(data, fs=1/dz, axis=axis,\
                                 nfft=len(data), detrend=detrend,\
                                 window='hanning', scaling='density')

    # Integration limits set by minimum and maximum vertical wavelength
    int_limit = np.logical_and(mgrid <= (1)/wlmin, mgrid >= (1)/wlmax)


    # Integrate between set limits
    variance = np.trapz(Psd[int_limit], x=mgrid[int_limit])

    if grid:
        return variance, mgrid, Psd
    else:
        return variance



def PE_isopycnal(N2, z, rho, strain, wl_min, wl_max,\
                     bin_idx, window=400, detrend=False):
    """
    Calculate internal wave potential energy based on isopycnal displacements
    and using neutral densities. (optional to do this) The function should work
    the same with normal density and an accompanying reference density profile.

    update:
    Right now this uses strain x dz with strain calculated from alex's code on
    the bray and fofonoff leveling method. it seems to work better than when I
    do it with my own density calculations (gets crazy numbers)

    """
    rho = rho - np.nanmean(rho)
    ref_rho = []
    counter = 0
    for cast in rho_neutral.T:
        rho_i = oc.vert_polyFit2(cast, z[:,0], 100, deg=1)
        ref_rho.append(rho_i)
        counter += 1

    ref_rho = (np.vstack(ref_rho)).T

#    rho = rho - np.nanmean(rho)
#    ref_rho = np.nanmean(rho, axis=1)


    eta = np.full_like(rho_neutral, np.nan)
    for i, cast in enumerate(rho_neutral.T):
        for k, elem in enumerate(cast):
            mask = np.abs(rho[k,i] - ref_rho[:,i])
            if np.nansum(mask) == 0:
                eta[k,i] = np.nan
            else:
#                print(np.nanargmin(mask))
                eta[k,i] = z[np.nanargmin(np.abs(rho_neutral[k,i] - ref_rho[:,i])),i] - z[k,i]
    drhodz = []
#    step = np.floor(ref_rho.shape[0]/window)
#    for i in zip(rho, ref_rho):
#
#        drhodz.append((np.nanmax(i)-np.nanmin(i))/\
#                      (z[:,0][np.argmax(i)] - z[:,0][np.argmin(i)]))
#    drhodz= np.squeeze(np.vstack(drhodz))
#    eta = (rho-ref_rho.T)
#    eta = np.vstack([etaIn/drhodzIn for etaIn, drhodzIn in zip(eta.T, drhodz)]).T

    eta = (rho-ref_rho)
    eta = eta**2

    # Assumes that strain is the gradient version of isopycnal displacements
    dz = np.nanmean(diff(z, axis=0))
    eta_psd = []

    # Use periodogram and integrate between target wavelengths
    eta1 = np.full((bin_idx.shape[0],eta.shape[1]), np.nan)
    for k, cast in enumerate(eta.T):
        for i, binIn in enumerate(bin_idx):
            eta1[i,k], f_grid, psd = PowerDens(cast[binIn], dz, wl_max,\
                                                wl_min, grid=True)
            eta_psd.append(psd)

    eta_psd = np.vstack(eta_psd)


    # Calculate mean Buoyancy Frequency for each bin
    N2mean = []
    for binIn in bin_idx:
        N2mean.append(np.nanmean(N2[binIn,:], axis=0))

    N2mean = np.vstack(N2mean)

    PE = 0.5*N2mean*eta1

    return PE, f_grid, eta_psd, N2mean


def PE_strain(N2, z, strain, wl_min, wl_max, bin_idx, detrend='constant'):
    """
    Calculate internal wave potential energy based on isopycnal displacements
    and using neutral densities. (optional to do this) The function should work
    the same with normal density and an accompanying reference density profile.

    update:
    Right now this uses strain x dz with strain calculated from alex's code on
    the bray and fofonoff leveling method. it seems to work better than when I
    do it with my own density calculations (gets crazy numbers)

    """
    # Assumes that strain is the gradient version of isopycnal displacements
    dz = np.nanmean(diff(z, axis=0))
    m = (2*np.pi)/(np.nanmean([wl_max, wl_min]))
    eta = (strain/m)
    eta_psd = []


    # Use periodogram and integrate between target wavelengths
    eta1 = np.full((bin_idx.shape[0],eta.shape[1]), np.nan)
    for k, cast in enumerate(eta.T):
        for i, binIn in enumerate(bin_idx):
            eta1[i,k], f_grid, psd = PowerDens(cast[binIn], dz, wl_max,\
                                                wl_min, grid=True,detrend=detrend)
            eta_psd.append(psd)

    eta_psd = np.vstack(eta_psd)


    # Calculate mean Buoyancy Frequency for each bin
    N2mean = []
    for binIn in bin_idx:
        N2mean.append(np.nanmean(N2[binIn,:], axis=0))

    N2mean = np.vstack(N2mean)


    PE = 0.5*eta1*N2mean

    return PE, f_grid, eta_psd, N2mean



def KE_UV(U, V, z, bin_idx, wl_min, wl_max, detrend=False):
    """

    Internal wave kinetic energy

    """

    Upoly = []
    for cast in U.T:
        fitrev = oc.vert_polyFit(cast, z[:,0], 100, deg=1)
        Upoly.append(fitrev)

    Upoly = np.vstack(Upoly).T
    U = U - Upoly

    Vpoly = []

    for cast in V.T:
        fitrev = oc.vert_polyFit(cast, z[:,0], 100, deg=1)
        Vpoly.append(fitrev)

    Vpoly = np.vstack(Vpoly).T
    V = V - Vpoly
    dz = np.nanmean(np.gradient(z, axis=0))
    KE_psd = []

    # Use periodogram and integrate between target wavelengths
    Pu = np.full((bin_idx.shape[0],U.shape[1]), np.nan)
    Pv = np.full((bin_idx.shape[0],U.shape[1]), np.nan)
    for k, (Ui, Vi) in enumerate(zip(U.T, V.T)):
        for i, binIn in enumerate(bin_idx):
            Pu[i,k], f_grid, psd = PowerDens(Ui[binIn], dz, wl_max,\
                                      wl_min, grid=True, detrend=detrend)
            Pv[i,k], f_grid, psd = PowerDens(Ui[binIn], dz, wl_max,\
                                      wl_min, grid=True, detrend=detrend)
            KE_psd.append(psd)

    KE_psd = np.vstack(KE_psd)


    # New Version
    KE = 0.5*(Pu + Pv)
    clean  = KE_psd < 1e-8
    KE_psd[clean] = np.nan

    return KE, f_grid, KE_psd

def wave_components_with_strain(ctd, ladcp, strain,\
                         ctd_bin_size=1024, ladcp_bin_size=1024,\
                         wl_min=300, wl_max=1000, plots=False):

    """
    Calculating Internal Wave Energy

    Internal wave energy calcuations following methods in waterman et al 2012.

    """

    # Load Hydrographic Data
    U, V, p_ladcp = oc.loadLADCP(ladcp)
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)
    SA = gsw.SA_from_SP(S, p_ctd, lon, lat)
    CT = gsw.CT_from_t(SA, T, p_ctd)
    N2, dump = gsw.stability.Nsquared(SA, CT, p_ctd, lat)


    maxDepth = 4000
    idx_ladcp = p_ladcp[:,-1] <= maxDepth
    idx_ctd = p_ctd[:,-1] <= maxDepth

    strain = strain[idx_ctd, :]
    S = S[idx_ctd,:]
    p_ctd = p_ctd[idx_ctd,:]
    U = U[idx_ladcp, :]
    V = V[idx_ladcp, :]
    p_ladcp = p_ladcp[idx_ladcp,:]
    # Bin CTD data
    ctd_bins = oc.binData(S, p_ctd[:,0], ctd_bin_size)
    # Bin Ladcp Data
    ladcp_bins = oc.binData(U, p_ladcp[:,0], ladcp_bin_size)

    # Depth and lat/long grids
    depths = np.vstack([np.nanmean(p_ctd[binIn]) for binIn in ctd_bins])
    dist = gsw.distance(lon, lat)
    dist = np.cumsum(dist)/1000
    dist = np.append(0,dist)


    # Calculate Potential Energy
    z = -1*gsw.z_from_p(p_ctd, lat)
    PE, PE_grid, eta_psd, N2mean = PE_strain(N2, z, strain, wl_min, wl_max, ctd_bins)

    # Calculate Kinetic Energy
    z = -1*gsw.z_from_p(p_ladcp, lat)
    KE, KE_grid, KE_psd = KE_UV(U, V, z, ladcp_bins, wl_min, wl_max, detrend=False)

    # Total Kinetic Energy
    Etotal = (KE + PE) * 1025 # Multiply by density to get Joules

    # wave components
    f = np.nanmean(gsw.f(lat))

    # version 2 omega calculation
    omega = f*np.sqrt((KE+PE)/(KE-PE))

    # version 2 omega calculation
    omega2 = np.abs((f**2)*((KE+PE)/(KE-PE)))
    rw = KE/PE
    w0 = ((f**2)*((rw+1)/(rw-1)))
    m = (2*np.pi)/np.mean((wl_min, wl_max))

    # version 1 kh calculation
    khi = m*np.sqrt(((f**2 - omega**2)/(omega**2 - N2mean)))

    # version 2 kh calculation
    kh = (m/np.sqrt(N2mean))*(np.sqrt(omega2 - f**2))
    mask = kh == 0
    kh[mask]= np.nan
    lambdaH = 1e-3*(2*np.pi)/kh

    if plots:
        plt.figure()
        plt.subplot(211)
        plt.loglog(KE_grid, KE_psd.T)
        plt.subplot(212)
        plt.loglog(PE_grid, .5*np.nanmean(N2)*eta_psd.T)

    if save_data:

        file2save = pd.DataFrame(lambdaH)
        file2save.index = np.squeeze(depths)
        file2save.to_excel('lambdaH_dec24.xlsx')
        file2save = pd.DataFrame(Etotal)
        file2save.index = np.squeeze(depths)
        file2save.to_excel('E_total.xlsx')



    return PE, KE, omega, m, kh, lambdaH, Etotal

def doppler_shifts(kh, ladcp, avg=1000, bin_size = 512):
    """
    Doppler shift the internal frequency to test for lee waves
    using the depth averaged floww
    """
    U, V, p_ladcp = oc.loadLADCP(ladcp)
    maxDepth = 4000
    idx_ladcp = p_ladcp[:,-1] <= maxDepth
    dz = int(np.nanmean(np.gradient(p_ladcp, axis=0)))
    window = int(np.ceil(avg/dz))
    Ubar = []
    for u, v in zip(U.T, V.T):
        mask = np.isfinite(u)
        u = u[mask]
        v = v[mask]
        u = np.nanmean(u[-window:])
        v = np.nanmean(v[-window:])
        Ubar.append(np.sqrt(u**2 + v**2))

    Ubar = np.vstack(Ubar)
    dshift = []
    for cast, ubar in zip(kh.T, Ubar):
        dshift.append(cast*ubar)

    dshift = np.vstack(dshift).T

    return dshift

def momentumFluxes(kh, m, N2,ladcp, z_ctd, bin_size=512, h0=1000):
    """
    Calculating vertical and horizontal momentum fluxes of internal waves within
    safe range
    """

    dshift = doppler_shifts(kh, ladcp)

    # Test whether K*U is between N and f
    f = (np.nanmean(gsw.f(lat)))

    dshiftTest = np.full_like(kh, np.nan)

    for i, dump in enumerate(kh.T):
        N2mean = np.nanmean(N2[:,i])
        dshiftTest[:,i] = np.logical_and(dshift[:,i]**2 >= f**2, dshift[:,i]**2<= N2mean)

    dshiftTest2 = dshiftTest == 0
    kh[dshiftTest2] = np.nan

    maxDepth = 4000
    idx_ladcp = z[:,-1] <= maxDepth
    idx_ctd = z_ctd[:,-1] <= maxDepth
    z = z[idx_ladcp,:]
    U = U[idx_ladcp,:]
    V = V[idx_ladcp,:]
    rho = rho_neutral[idx_ctd,:]
    bins = oc.binData(U, z[:,0], bin_size)
    Umean = np.vstack([np.nanmean(U[binIn,:], axis=0) for binIn in bins])
    Vmean = np.vstack([np.nanmean(V[binIn,:], axis=0) for binIn in bins])
    Umag = np.sqrt(Umean**2 + Vmean**2)

    bins = oc.binData(rho, z_ctd[:,0], bin_size)
    rho0 = np.vstack([np.nanmean(rho[binIn,:], axis=0) for binIn in bins])



    tau = .5*kh*rho0*Umag*np.sqrt((N2/Umag**2)-kh**2)
