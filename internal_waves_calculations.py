"""
Created on December 24th 15:44:35 2017

@author: manishdevana
This toolbox calculates internal wave properties and energetics
"""


import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import data_load
import gsw
import oceans as oc
from scipy import interpolate
import cmocean
default_params = {
        'nfft': 2048,
        'plots': True,
        'rho0': 1025
        }

def reset_test():
    """
    This loads the data for testing the functions test
    (its also all the data for the project)
    """

    ladcp, ctd = data_load.load_data()
#    rho_neutral =  np.genfromtxt('neutral_rho.csv', delimiter=',')
    strain = np.genfromtxt('strain.csv', delimiter=',')
#    N2 = np.genfromtxt('ref_N2.csv', delimiter=',')
    wl_max = 1500
    wl_min = 500
    ctd_bin_size = 1024
    ladcp_bin_size = 1024
    nfft = 2048
    U, V, p_ladcp = oc.loadLADCP(ladcp)
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)

    return ladcp, ctd, strain,wl_max, wl_min, ctd_bin_size, ladcp_bin_size, nfft, S, T, p_ctd, U, V, p_ladcp, lat, lon

def PowerDens(data, dz, wlmax, wlmin, axis=0, grid=False,
              nfft=None, window='hanning', detrend='constant'):
    """
    Using periodogram function to estimate power spectral density

    PARAMETERS
    ----------
    data: data array for input (assumes column wise data (axis=0))
    dz: vertical sampling rate
    wlmax: maximimum vertical wavelength integration limit
    wlmin: minimum vertical wavelength integration limit
    axis: axis to perform PSD calculations on
    grid: if True, returns frequency grid
    nfft: default nfft is length(data) unless specified otherwise
    window: window function to reduce loss on FFT, default='hanning', refer to
                scipy docs for other options
    detrend: detrending options , defaults to 'constant', other option='linear'

    RETURN
    ------
    variance: The integrated psd of the profiles
    mgrid: frequency grid as 1/vertical wavelength
    Psd: full power spectral density spectrum for input data
    peaks: Frequency of maximum Power spectral density

    """

    if not nfft:
        nfft = len(data)
    mgrid, Psd = sig.periodogram(data, fs=1/dz, axis=axis,
                                 nfft=nfft, detrend=detrend,
                                 window=window, scaling='density')

    # Integration limits set by minimum and maximum vertical wavelength
    int_limit = np.logical_and(mgrid <= (1)/wlmin, mgrid >= (1)/wlmax)

    # Integrate between set limits using trapezoid rule
    variance = np.trapz(Psd[int_limit], x=mgrid[int_limit])

    # find wavelength of max energy density
    idx = np.argmax(Psd)
    peaks = mgrid[idx]

    if grid:
        return variance, mgrid, Psd, peaks
    else:
        return variance


def PE_isopycnal(N2, z, b, wl_min, wl_max,
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
            eta1[i,k], f_grid, psd = PowerDens(cast[binIn], dz, wl_max,
                                                wl_min, grid=True)
            eta_psd.append(psd)

    eta_psd = np.vstack(eta_psd)


    # Calculate mean Buoyancy Frequency for each bin
    N2mean = []
    for binIn in bin_idx:
        N2mean.append(np.nanmean(N2[binIn,:], axis=0))

    N2mean = np.vstack(N2mean)

    PE = 1027*0.5*N2mean*eta1

    return PE, f_grid, eta_psd, N2mean


def potential_energy(S, T, p, lat,
                     min_w=400,
                     max_w=1500,
                     nfft=1024,
                     axis=0,
                     window='hanning',
                     detrend='constant',
                     with_strain=True):
    """
    Calculate internal wave potential energy following waterman et al 2012 and
    Meyer et al. 2016. isopyncal displacement is calculated using a reference
    densiity field constructed by adiabtic leveling following Bray and Fofonoff
    1981. Gives the option to calculate eta directly or calculate strain and
    integrated between depth intervals (not sure if this is valid but seems to
    work much better).

    PARAMETERS
    ----------
    S: Salinity
    T: Temperature
    p: pressure
    lat: latitude
    min_w: lower vertical wavelength limit for integration
    max_w: upper vertical wavelength limit for Integration
    nfft: number of points in fft calculation (default=length of data)
    axis: axis to perform calculations default=0
    window: window function to reduce variance loss in FFT (default='hanning')
    detrend: detrend method for fft (default='constant')
    with_strain: option to calculate eta from strain instead of from density
                    surfaces. (default=True)

    RETURNS
    -------
    PE: Potential energy in joules/m^3
    PE_psd: periodogram used in calculations
    peaks: wavelenths of max power in periodograms
    f_grid: frequency grid (inverse vertical wavelengths
    """

    # Use adiabatic leveling to generate reference fields
    N2_ref, N2, strain, p_mid = oc.adiabatic_level(S, T, p, lat)


def PE(N2, z, eta,
              wl_min,
              wl_max,
              bin_idx,
              nfft=2048,
              detrend='linear'):
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
    dz = np.nanmean(np.diff(z, axis=0))
    # Use periodogram and integrate between target wavelengths
    eta1 = np.full((bin_idx.shape[0], eta.shape[1]), np.nan)

    peaks = []
    eta_psd = []

    for k, cast in enumerate(eta.T):
        for i, binIn in enumerate(bin_idx):
            data_in = cast[binIn]
            mask = np.isfinite(data_in)
            good = data_in[mask]
            eta1[i, k], f_grid, psd, peaks_i = PowerDens(good, dz, wl_max,
                                            wl_min, grid=True, nfft=nfft, detrend=detrend)
            eta_psd.append(psd)
            peaks.append(peaks_i)

    eta_psd = np.vstack(eta_psd)


    # Calculate mean Buoyancy Frequency for each bin using a mean vertical
    # buoyancy profile for the entire grid
    N2ref = np.nanmean(N2, axis=1)
    N2mean = []
    for binIn in bin_idx:
        N2mean.append(np.nanmean(N2ref[binIn], axis=0))



    N2mean = np.vstack(N2mean)
    N2mean2 = np.tile(N2mean, [1,N2.shape[1]])

    PE = 0.5 * eta1 * N2mean

    return PE, f_grid, eta_psd, N2mean2, np.vstack(peaks)

def PE_strain(N2, z, strain, wl_min,
              wl_max, bin_idx, nfft=2048,
              detrend='constant'):
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
    dz = np.nanmean(np.diff(z, axis=0))
    m = (2*np.pi)/(np.nanmean([wl_max, wl_min]))
    eta = strain
#    eta = (strain/m)
    eta_psd = []

    # Use periodogram and integrate between target wavelengths
    eta1 = np.full((bin_idx.shape[0], eta.shape[1]), np.nan)
    peaks = []
    for k, cast in enumerate(eta.T):
        for i, binIn in enumerate(bin_idx):
            eta1[i, k], f_grid, psd, peaks_i = PowerDens(cast[binIn], dz, wl_max,
                                            wl_min, grid=True, nfft=nfft, detrend=detrend)
            eta_psd.append(psd)
            peaks.append(peaks_i)

    eta_psd = np.vstack(eta_psd)


    # Calculate mean Buoyancy Frequency for each bin using a mean vertical
    # buoyancy profile for the entire grid
    N2ref = np.nanmean(N2, axis=1)
    N2mean = []
    for binIn in bin_idx:
        N2mean.append(np.nanmean(N2ref[binIn], axis=0))



    N2mean = np.vstack(N2mean)
    N2mean2 = np.tile(N2mean, [1,N2.shape[1]])

    PE = 0.5*1027*eta1*N2mean

    return PE, f_grid, eta_psd, N2mean2, np.vstack(peaks)


def KE_UV(U, V, z, bin_idx, wl_min, wl_max, lc=400, nfft=2048, detrend='constant'):
    """

    Internal wave kinetic energy

    """


    # Clean Up velocity data (u' = U - u_bar)
    Upoly = []
    for cast in U.T:
        fitrev = oc.vert_polyFit(cast, z[:, 0], 100, deg=1)
        Upoly.append(fitrev)

    Upoly = np.vstack(Upoly).T
    U = U - Upoly


    dz = 8  # This is the vertical spacing between measurements in metres.
    lc = lc  # This is the cut off vertical scale in metres, the filter will remove variability smaller than this.
    mc = 1./lc  # Cut off wavenumber.
    normal_cutoff = mc*dz*2.  # Nyquist frequency is half 1/dz.
    a1, a2 = sig.butter(4, normal_cutoff, btype='lowpass')  # This specifies you use a lowpass butterworth filter of order 4, you can use something else if you want
    for i in range(U.shape[1]):
        mask = ~np.isnan(U[:,i])
        U[mask,i] = sig.filtfilt(a1, a2, U[mask,i])

    Vpoly = []

    for cast in V.T:
        fitrev = oc.vert_polyFit(cast, z[:, 0], 100, deg=1)
        Vpoly.append(fitrev)

    Vpoly = np.vstack(Vpoly).T
    V = V - Vpoly


    for i in range(V.shape[1]):
        mask = ~np.isnan(U[:,i])
        V[mask,i] = sig.filtfilt(a1, a2, V[mask,i])
    dz = np.nanmean(np.gradient(z, axis=0))
    KE_psd = []
    peaks = []
    # Use periodogram and integrate between target wavelengths
    Pu = np.full((bin_idx.shape[0],U.shape[1]), np.nan)
    Pv = np.full((bin_idx.shape[0],U.shape[1]), np.nan)
    for k, (Ui, Vi) in enumerate(zip(U.T, V.T)):
        for i, binIn in enumerate(bin_idx):
            Pu[i, k], f_grid, psd, u_peaks = PowerDens(Ui[binIn], dz, wl_max,
                                      wl_min, grid=True, nfft=nfft, detrend=detrend)
            Pv[i, k], f_grid, psd1, v_peaks = PowerDens(Vi[binIn], dz, wl_max,
                                      wl_min, grid=True, nfft=nfft, detrend=detrend)
            KE_psd.append(.5 * (psd + psd1))
            peaks.append([u_peaks, v_peaks])

    KE_psd = np.vstack(KE_psd)



    # New Version
    KE = 0.5*(Pu + Pv)
    clean  = KE_psd < 1e-8
    KE_psd[clean] = np.nan

    return KE, f_grid, KE_psd, U, V, np.vstack(peaks)

def momentumFlux(kh, m, N2mean, f):
    """
    Calculating internal wave momentum fluxes (N22f 2)a2
    """
    a = kh/m

    cgz = (a**2)*(N2mean - f)/(m*((1+a**2)**1.5)*(f**2 + a*N2mean)**.5)

    Ev = Etotal*cgz

    plt.figure()
    plt.contourf(Ev)
    plt.contour(m, color='k')
    plt.gca().invert_yaxis()






def wave_components_with_strain(ctd, ladcp, strain,
                                rho0=default_params['rho0'],
                                ctd_bin_size=1024, ladcp_bin_size=1024,
                                wl_min=300, wl_max=1000,
                                nfft=default_params['nfft'],
                                plots=default_params['plots'], save_data=False):

    """
    Calculating Internal Wave Energy

    Internal wave energy calcuations following methods in waterman et al 2012.

    """

    # Load Hydrographic Data
    g = 9.8
    U, V, p_ladcp = oc.loadLADCP(ladcp)
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)
    SA = gsw.SA_from_SP(S, p_ctd, lon, lat)
    CT = gsw.CT_from_t(SA, T, p_ctd)
    N2, dump = gsw.stability.Nsquared(SA, CT, p_ctd, lat)

    maxDepth = 4000
    idx_ladcp = p_ladcp[:, -1] <= maxDepth
    idx_ctd = p_ctd[:, -1] <= maxDepth

    strain = strain[idx_ctd, :]
    S = S[idx_ctd,:]
    T = T[idx_ctd,:]
    p_ctd = p_ctd[idx_ctd, :]
    U = U[idx_ladcp, :]
    V = V[idx_ladcp, :]
    p_ladcp = p_ladcp[idx_ladcp, :]
    rho = oc.rhoFromCTD(S, T, p_ctd, lon, lat)
    # Bin CTD data
    ctd_bins = oc.binData(S, p_ctd[:, 0], ctd_bin_size)
    # Bin Ladcp Data
    ladcp_bins = oc.binData(U, p_ladcp[:, 0], ladcp_bin_size)

    # Depth and lat/long grids
    depths = np.vstack([np.nanmean(p_ctd[binIn]) for binIn in ctd_bins])
    dist = gsw.distance(lon, lat)
    dist = np.cumsum(dist)/1000
    dist = np.append(0,dist)


    # Calculate Potential Energy
    z = -1*gsw.z_from_p(p_ctd, lat)
    PE, PE_grid, eta_psd, N2mean, pe_peaks = PE_strain(N2, z, strain,
                                             wl_min, wl_max, ctd_bins, nfft=2048)

    # Calculate Kinetic Energy
    z = -1*gsw.z_from_p(p_ladcp, lat)
    KE, KE_grid, KE_psd, Uprime, Vprime, ke_peaks = KE_UV(U, V, z, ladcp_bins,
                                wl_min, wl_max, lc=wl_min-50,
                                nfft=2048, detrend='constant')

    # Total Kinetic Energy
    Etotal = 1027*(KE + PE) # Multiply by density to get Joules

    # wave components
    f = np.nanmean(gsw.f(lat))

    # version 2 omega calculation
    omega = f*np.sqrt((KE+PE)/(KE-PE))

    # version 2 omega calculation
    omega2 = np.abs((f**2)*((KE+PE)/(KE-PE)))
    rw = KE/PE
    w0 = ((f**2)*((rw+1)/(rw-1)))
#    m = (2*np.pi)/np.mean((wl_min, wl_max))
    m = np.nanmean(ke_peaks, axis=1)
    m = ke_peaks[:,0]
    m = m.reshape(omega.shape)
    m = (2*np.pi)*m

    # version 1 kh calculation
    khi = m*np.sqrt(((f**2 - omega**2)/(omega**2 - N2mean)))

    # version 2 kh calculation
    kh = (m/np.sqrt(N2mean))*(np.sqrt(omega2 - f**2))
    mask = khi == 0
    khi[mask]= np.nan
    lambdaH = 1e-3*(2*np.pi)/khi

    # Get coherence of u'b' and v'b' and use to estimate horizontal wavenumber
    # components. This uses the previously binned data but now regrids velocity
    # onto the density grid so there are the same number of grid points
    b = (-g*rho)/rho0
    b_poly = []
    z = -1*gsw.z_from_p(p_ctd, lat)
    fs = 1/np.nanmean(np.diff(z, axis=0))
    for cast in b.T:
        fitrev = oc.vert_polyFit(cast, z[:, 0], 100, deg=1)
        b_poly.append(fitrev)

    b_poly = np.vstack(b_poly).T
    b_prime = b - b_poly

    dz = 1/fs  # This is the vertical spacing between measurements in metres.
    lc = wl_min-50  # This is the cut off vertical scale in metres, the filter will remove variability smaller than this.
    mc = 1./lc  # Cut off wavenumber.
    normal_cutoff = mc*dz*2.  # Nyquist frequency is half 1/dz.
    a1, a2 = sig.butter(4, normal_cutoff, btype='lowpass')  # This specifies you use a lowpass butterworth filter of order 4, you can use something else if you want
    for i in range(b_prime.shape[1]):
        mask = ~np.isnan(b_prime[:,i])
        b_prime[mask,i] = sig.filtfilt(a1, a2, b_prime[mask,i])

    ub = []
    vb = []

    for i in range(ctd_bins.shape[0]):

        Uf = interpolate.interp1d(p_ladcp[ladcp_bins[i,:]].squeeze(),
                                        Uprime[ladcp_bins[i, :], :],
                                        axis=0, fill_value='extrapolate')

        Vf = interpolate.interp1d(p_ladcp[ladcp_bins[i,:]].squeeze(),
                                        Vprime[ladcp_bins[i, :], :],
                                        axis=0, fill_value='extrapolate')
        new_z = p_ctd[ctd_bins[i,:],0]
        u_f, ub_i = sig.coherence(b_prime[ctd_bins[i,:],:],
                                   Uf(new_z), nfft=nfft, fs=fs, axis=0)
        v_f, vb_i = sig.coherence(b_prime[ctd_bins[i,:],:],
                                   Vf(new_z), nfft=nfft, fs=fs, axis=0)

        ub.append(ub_i)
        vb.append(vb_i)

    ub = np.hstack(ub).T
    vb = np.hstack(vb).T


    # Random plots (only run if youre feeling brave)
    m_plot = np.array([(2*np.pi)/wl_max,
                       (2*np.pi)/wl_max, (2*np.pi)/wl_min,
                       (2*np.pi)/wl_min])

    if plots:
        plt.figure(figsize=[12,6])
        plt.subplot(121)
        plt.loglog(KE_grid, KE_psd.T, linewidth=.6, c='b', alpha=.1)
        plt.loglog(KE_grid, np.nanmean(KE_psd, axis=0).T, lw=1.5, c='k')
        ylims = plt.gca().get_ylim()
        ylim1 = np.array([ylims[0], ylims[1]])
        plt.plot(m_plot[2:], ylim1, lw=1,
                 c='k', alpha=.5,
                 linestyle='dotted')
        plt.plot(m_plot[:2], ylim1, lw=1,
                 c='k', alpha=.5,
                 linestyle='dotted')
        plt.ylim(ylims)
        plt.ylabel('Kinetic Energy Density')
        plt.xlabel('Vertical Wavenumber')
        plt.gca().grid(True, which="both", color='k', linestyle='dotted', linewidth=.2)
        plt.subplot(122)
        plt.loglog(PE_grid, .5*np.nanmean(N2)*eta_psd.T,
                   lw=.6, c='b', alpha=.1)
        plt.loglog(KE_grid, .5*np.nanmean(N2)*np.nanmean(eta_psd, axis=0).T,
                   lw=1.5, c='k')
        plt.plot(m_plot[2:], ylim1, lw=1,
                 c='k', alpha=.5,
                 linestyle='dotted')
        plt.plot(m_plot[:2], ylim1, lw=1,
                 c='k', alpha=.5,
                 linestyle='dotted')
        plt.ylim(ylims)
        plt.gca().grid(True, which="both", color='k', linestyle='dotted', linewidth=.2)
        plt.ylabel('Potential Energy Density')
        plt.xlabel('Vertical Wavenumber')

        plt.figure()
        Kemax = np.nanmax(KE_psd, axis=1)
        kespots = np.nanargmax(KE_psd, axis=1)
        ax = plt.gca()
        ax.scatter(KE_grid[kespots],Kemax , c='blue', alpha=0.3, edgecolors='none')
        ax.set_yscale('log')
        ax.set_xscale('log')

        plt.figure(figsize=[12,6])
        plt.subplot(121)
        plt.semilogx(u_f, ub.T, linewidth=.5, alpha=.5)
        plt.gca().grid(True, which="both", color='k', linestyle='dotted', linewidth=.2)
        plt.subplot(122)
        plt.semilogx(v_f, vb.T, linewidth=.5)
        plt.gca().grid(True, which="both", color='k', linestyle='dotted', linewidth=.2)
#        plt.xlim([10**(-2.5), 10**(-2)])

        plt.figure()
        ub_max = np.nanmax(ub, axis=1)
        kespots = np.argmax(ub, axis=1)
        ax = plt.gca()
        ax.scatter(u_f[kespots],ub_max , c='blue', alpha=0.3, edgecolors='none')
        ax.set_xscale('log')
        ax.set_xlim([1e-3, 1e-5])

        Kemax = np.nanmax(.5*np.nanmean(N2)*eta_psd.T, axis=1)
        kespots = np.nanargmax(.5*np.nanmean(N2)*eta_psd.T, axis=1)
        ax = plt.gca()
        ax.scatter(PE_grid[kespots],Kemax , c='red', alpha=0.3, edgecolors='none')
        ax.set_yscale('log')
        ax.set_xscale('log')


        # Peaks lots
        plt.figure()
        mask = np.isfinite(Etotal)
        Etotal[~mask]= 0
        distrev = np.tile(dist, [kh.shape[0],1])
        depthrev = np.tile(depths, [1, kh.shape[1]])
        plt.pcolormesh(distrev, depthrev, Etotal, shading='gouraud')
        plt.gca().invert_yaxis()

        plt.figure()
        plt.pcolormesh(dist, p_ladcp.squeeze(),
                       Uprime, cmap=cmocean.cm.balance,
                       shading='flat')
        levels = np.arange(np.nanmin(Etotal), np.nanmax(Etotal)+.5,.05)
        plt.contour(distrev, depthrev, Etotal)
        plt.gca().invert_yaxis()

    if save_data:

        file2save = pd.DataFrame(lambdaH)
        file2save.index = np.squeeze(depths)
        file2save.to_excel('lambdaH_dec24.xlsx')
        file2save = pd.DataFrame(Etotal)
        file2save.index = np.squeeze(depths)
        file2save.to_excel('E_total.xlsx')

    return PE, KE, omega, m, kh, lambdaH,\
            Etotal, khi, Uprime, Vprime, b_prime,\
            ctd_bins, ladcp_bins, KE_grid, PE_grid,\
            ke_peaks, pe_peaks, dist, depths, KE_psd,\
            eta_psd, N2, N2mean

def horizontal_azimuth(Uprime, Vprime, dz, wl_min, wl_max, axis=0, nfft=1024):
    """
    Attempt to decompose horizontal wave vector
    Following methods used in Polzin 2007 (internal waves in eddies or something like that)
    """

    # U'* x b'
    Uspec = np.fft.fft(Uprime, n=nfft, axis=axis)
    Vspec = np.fft.fft(Vprime, n=nfft, axis=axis)
    fs = 1./dz
    fmin = 0
    fi = 1/nfft
    fmax = .5*fs
    mx = np.linspace(fmin, fmax, num=nfft)

    int_limit = np.logical_and(mx <= (1)/wl_min, mx >= (1)/wl_max)
    Uspec = np.nanmean(Uspec[int_limit,:], axis=axis)
    Vspec = np.nanmean(Vspec[int_limit,:], axis=axis)

    theta = []
    for Uin, Vin in zip(Uspec.T, Vspec.T):
        u_conj = np.conj(Uin)
        v_prime = Vin
        u_prime = Uin
        v_conj = np.conj(Vin)
        theta.append(np.arctan(2*np.real((u_conj*v_prime)/(u_conj*u_prime - v_conj*v_prime)))/2)

    theta = np.vstack(theta).T


    return theta

#
#def doppler_shifts(kh, ladcp, avg=1000, bin_size = 512):
#    """
#    Doppler shift the internal frequency to test for lee waves
#    using the depth averaged floww
#    """
#    U, V, p_ladcp = oc.loadLADCP(ladcp)
#    maxDepth = 4000
#    idx_ladcp = p_ladcp[:,-1] <= maxDepth
#    dz = int(np.nanmean(np.gradient(p_ladcp, axis=0)))
#    window = int(np.ceil(avg/dz))
#    Ubar = []
#    for u, v in zip(U.T, V.T):
#        mask = np.isfinite(u)
#        u = u[mask]
#        v = v[mask]
#        u = np.nanmean(u[-window:])
#        v = np.nanmean(v[-window:])
#        Ubar.append(np.sqrt(u**2 + v**2))
#
#    Ubar = np.vstack(Ubar)
#    dshift = []
#    for cast, ubar in zip(kh.T, Ubar):
#        dshift.append(cast*ubar)
#
#    dshift = np.vstack(dshift).T
#
#    return dshift
#
#
#def energyFlux(kh, m, KE, PE, N2, U, V):
#    """
#    Calculates the energy flux of an internal wave
#    """
#
#
#
#def momentumFluxes(kh, m, N2,ladcp, z_ctd, bin_size=512, h0=1000):
#    """
#    Calculating vertical and horizontal momentum fluxes of internal waves within
#    safe range
#    """
#
#    dshift = doppler_shifts(kh, ladcp)
#
#    # Test whether K*U is between N and f
#    f = (np.nanmean(gsw.f(lat)))
#
#    dshiftTest = np.full_like(kh, np.nan)
#
#    for i, dump in enumerate(kh.T):
#        N2mean = np.nanmean(N2[:,i])
#        dshiftTest[:,i] = np.logical_and(dshift[:,i]**2 >= f**2, dshift[:,i]**2<= N2mean)
#
#    dshiftTest2 = dshiftTest == 0
#    kh[dshiftTest2] = np.nan
#
#    maxDepth = 4000
#    idx_ladcp = z[:,-1] <= maxDepth
#    idx_ctd = z_ctd[:,-1] <= maxDepth
#    z = z[idx_ladcp,:]
#    U = U[idx_ladcp,:]
#    V = V[idx_ladcp,:]
#    rho = rho_neutral[idx_ctd,:]
#    bins = oc.binData(U, z[:,0], bin_size)
#    Umean = np.vstack([np.nanmean(U[binIn,:], axis=0) for binIn in bins])
#    Vmean = np.vstack([np.nanmean(V[binIn,:], axis=0) for binIn in bins])
#    Umag = np.sqrt(Umean**2 + Vmean**2)
#
#    bins = oc.binData(rho, z_ctd[:,0], bin_size)
#    rho0 = np.vstack([np.nanmean(rho[binIn,:], axis=0) for binIn in bins])
#
#
#
#    tau = .5*kh*rho0*Umag*np.sqrt((N2/Umag**2)-kh**2)
