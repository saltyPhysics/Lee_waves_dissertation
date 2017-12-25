#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 12:05:57 2017

@author: manishdevana

Further experiments with ray tracing after doing wave property calculations
(runs calculations in internal_wave_properties)

"""

import numpy as np
import matplotlib.pyplot as plt
import data_load
import gsw
import oceans as oc
import Internal_wave_properties_REV as iwp
import pandas as pd


#   Testdata load
def internal_wave_properties(save_only=True):
    ladcp, ctd, bathy = data_load.load_data()
    rho_neutral = np.genfromtxt('neutral_rho.csv', delimiter=',')
    strain = np.genfromtxt('strain.csv', delimiter=',')
    wl_max=350
    wl_min=100
    lambdaH, kh, omega, N2, dist, depths,\
        U2, V2, p_ladcp, Uspec, Vspec,\
        etaSpec, aspect, Ek, Ep, Etotal = iwp.frequencyEstimator(ctd, ladcp, bathy,\
                                            rho_neutral,strain, wl_max=500, full_set=True)

    if not save_only:
        return lambdaH, kh, omega, N2, dist, depths,\
                U2, V2, p_ladcp, Uspec, Vspec, etaSpec, aspect


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



def horizontal_wave_vector_decomposition(Uspec, Vspec):
    """
    Attempt to decompose horizontal wave vector
    Following methods used in Polzin 2007 (internal waves in eddies or something like that)
    """

    # U'* x b'
    theta = []
    for Uin, Vin in zip(Uspec, Vspec):
        u_conj = np.conj(Uin[:,1])
        v_prime = Vin[:,1]
        u_prime = Uin[:,1]
        v_conj = np.conj(Vin[:,1])
        theta.append(np.arctan(2*np.real((u_conj*v_prime)/(u_conj*u_prime - v_conj*v_prime)))/2)

    theta = np.vstack(theta).T

    k = -kh*np.cos(theta)
    k_wave = (2*np.pi/k)*1e-3
    l = kh*np.sin(theta)
    l_wave = (2*np.pi/l)*1e-3

    return k, l

def lee_wave_tests(kh, omega, N2, ctd, ladcp, dist, depths, error_factor=2, plots=False):
    """
    Testing whether or not the observations can be attributed to lee waves
    """
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)

#    k, l = horizontal_wave_vector_decomposition(Uspec, Vspec)


    dshift = doppler_shifts(kh, ladcp)

    # Test whether K*U is between N and f
    f = (np.nanmean(gsw.f(lat)))

    dshiftTest = np.full_like(kh, np.nan)
    test_final = np.full_like(kh, np.nan)

    for i, dump in enumerate(kh.T):
        N2mean = np.nanmean(N2[:,i])
        testA = np.abs(dshift[:,i]**2-f**2) <= f**2
        testB = np.abs(dshift[:,i]**2-f**2) <= N2mean
        test_final[:,i] = np.logical_and(testA, testB)
        dshiftTest[:,i] = np.logical_and(dshift[:,i]**2 >= f**2, dshift[:,i]**2<= N2mean)

    dshiftTest2 = dshiftTest == 0
    mask = np.logical_not(test_final)
    kh[mask] = np.nan
    omega[mask] = np.nan
    lambdaH[mask] = np.nan
    k[mask] = np.nan
    l[mask] = np.nan
#    kh[dshiftTest2] = np.nan
#    omega[dshiftTest2] = np.nan
#    lambdaH[dshiftTest2] = np.nan
#    k[dshiftTest2] = np.nan
#    l[dshiftTest2] = np.nan

    file2save = pd.DataFrame(kh)
    file2save.index = np.squeeze(depths)
    file2save.to_excel('Kh_masked.xlsx')
    file2save = pd.DataFrame(omega)
    file2save.index = np.squeeze(depths)
    file2save.to_excel('omega_masked.xlsx')
    file2save = pd.DataFrame(lambdaH)
    file2save.index = np.squeeze(depths)
    file2save.to_excel('lambda_masked.xlsx')
    np.savetxt('kh_masked2.csv', kh)
    np.savetxt('k_masked2.csv', k)
    np.savetxt('l_masked2.csv', l)
    np.savetxt('omega_masked.csv', omega)



    # Test phase of velocity and isopycnal perturbations
    stns = np.arange(1, 1+S.shape[1])
    stns = np.expand_dims(stns, 1)
    stns = np.repeat(stns.T, kh.shape[0], axis=0)

    np.savetxt('dshift_mask.csv',dshiftTest)

    if plots:
        fig = plt.figure()
        plt.contourf(dist, np.squeeze(p_ladcp), U, cmap='seismic')
        plt.colorbar()
        plt.pcolormesh(dist, np.squeeze(depths), dshiftTest, cmap='binary', alpha=.2)
        plt.fill_between(dist, bathy, 4000, color = '#B4B4B4')
        plt.ylim(0, 4000)
        plt.gca().invert_yaxis()
        plt.title("u' with bins with f < Kh*U < N")
        plt.xlabel('Distance Along Transect (km)')
        plt.ylabel('Pressure (dB)')
        for i in range(U.shape[1]):
            plt.annotate(str(i+1), (dist[i], 3600))






def energy_flux(Uspec, Vspec, Ek, Ep, ctd, depths, m=300):
    """
    Use to calculate energy flux for a given internal wave following the equations
    shown in Cusack et al. 2017

    energy flux in the direction of the wave (direction of Cg)
    This function will return the total energy flux and the x,y, and z components
    of energy flux


    """
    m = 2*np.pi/m
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)
    rho = oc.rhoFromCTD(S, T, p_ctd, lon, lat)
    maxDepth = 4000
    idx_ctd = p_ctd[:,-1] <= maxDepth
    p_ctd = p_ctd[idx_ctd,:]
    rho = rho[idx_ctd,:]

    ctd_bins = oc.binData(rho, p_ctd[:,0], bin_size=512)
    rho2 = np.vstack([np.nanmean(rho[binIn]) for binIn in ctd_bins])
    # Energy Density = Kinetic Energy + Potential Energy
    E = (Ek + Ep)*np.nanmean(rho2)
    k, l = horizontal_wave_vector_decomposition(Uspec, Vspec)





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

    np.savetxt('wave_momentum.csv', tau)
















def perturbation_comparisons(ctd, ladcp, bin_size=512, plots=False):
    """
    Comparison of velocity and density perturbations to see if they are in or
    out of phase
    """
    # Load ctd data
    S, T, p_ctd, lat, lon = oc.loadCTD(ctd)

    # Calculate sigma0 and filter it
    rho = oc.rhoFromCTD(S, T, p_ctd, lon, lat)

    rho_poly = []
    for cast in rho.T:
        fitrev = oc.vert_polyFit2(cast, p_ctd[:,0], 100, deg=2)
        rho_poly.append(fitrev)

    rho_poly = np.vstack(rho_poly).T
    rho2 = rho - rho_poly

    # load ladcp data and filter it
    U, V, p_ladcp, uup, vup, udo, vdo = oc.loadLADCP(ladcp, full_set=True)

    uup, vup = oc.velocityFiltered(uup, vup, p_ladcp, deg=1)
    udo, vdo = oc.velocityFiltered(udo, vdo, p_ladcp, deg=1)

    U2, V2 = oc.velocityFiltered(U, V, p_ladcp, deg=1)

    # Apply 100m wide box filter to clean up perturbations
    for i, dump in enumerate(rho2.T):
        rho2[:,i] = oc.verticalBoxFilter1(rho2[:,i], p_ctd[:,0], box=100)
        U2[:,i] = oc.verticalBoxFilter1(U2[:,i], p_ladcp[:,0], box=100)
        V2[:,i] = oc.verticalBoxFilter1(V2[:,i], p_ladcp[:,0], box=100)
        uup[:,i] = oc.verticalBoxFilter1(uup[:,i], p_ladcp[:,0], box=100)
        vup[:,i] = oc.verticalBoxFilter1(vup[:,i], p_ladcp[:,0], box=100)
        udo[:,i] = oc.verticalBoxFilter1(udo[:,i], p_ladcp[:,0], box=100)
        vdo[:,i] = oc.verticalBoxFilter1(vdo[:,i], p_ladcp[:,0], box=100)



    if plots:

        count = 0
        for rho0, u1, u2 in zip(rho2.T, uup.T, udo.T):
            count += 1
            fig = plt.figure(figsize=(3.88, 7.62))
            ax1 = fig.add_subplot(111)
            ax2 = ax1.twiny()
            l1 = ax1.plot(rho0, p_ctd[:,0],'r', label='\rho')
            l2 = ax2.plot(u1, p_ladcp,'b', label="up'")
            l3 = ax2.plot(u2, p_ladcp,'k', label="down'")
            ax2.set_xlabel("u'")
            ax1.set_xlabel("rho'")
            ax1.set_ylim((0, 4000))
            plt.gca().invert_yaxis()
            plt.title('station ' + str(count), y=1.1)
            plt.savefig('figures/perturbation_comparisons/station ' + str(count) + '.png',\
                        bbox_inches='tight', dpi=400)
            plt.close()
