#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 14:17:31 2017

This uses my ray tracing model to run different experiments on ray tracing.

@author: manishdevana
"""

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import data_load
import gsw
import oceans as oc
from scipy import interpolate
import cmocean
import internal_waves_calculations as iwc
import ray_tracing as rt
import os
import plotly.plotly as py
import plotly.graph_objs as go

from mpl_toolkits.mplot3d import Axes3D
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
#
#
#
#ladcp, ctd, strain, wl_max, wl_min,\
# ctd_bin_size, ladcp_bin_size, nfft = iwc.reset_test()

def load():
    ladcp, ctd, strain, wl_max, wl_min,\
        ctd_bin_size, ladcp_bin_size, nfft = iwc.reset_test()

    PE, KE, omega, m, kh, lambdaH, Etotal,\
     khi, Uprime, Vprime, b_prime, ctd_bins,\
     ladcp_bins, KE_grid, PE_grid, ke_peaks,\
     pe_peaks, dist, depths, KE_psd, eta_psd, N2, N2mean = iwc.wave_components_with_strain(ctd,\
     ladcp, eta, wl_min=wl_min, wl_max=wl_max, plots=False)


def plots():


    fig = plt.figure(figsize=[10,10])
    plt.subplot(121)
    plt.imshow(Uprime, cmap=cmocean.cm.amp, aspect=.03)
    plt.colorbar()
    plt.title("U'")
    plt.subplot(122)
    plt.imshow(Vprime, cmap=cmocean.cm.amp, aspect=.03)
    plt.colorbar()
    plt.title("V'")


    u2  = Uprime - np.nanmean(Uprime)
    test2 = np.fft.fft2(u2, norm='ortho')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(dist, p_ladcp, np.abs(np.fft.fftshift(test2)))


    plt.figure()
    plt.hist(np.log10(kh.flatten()), bins=10,
             normed=False, range=(np.nanmin(np.log10(kh)),
                                  np.nanmax(np.log10(kh))))

    plt.figure()
    plt.hist(lambdaH.flatten(), bins=20,
             normed=False, range=(np.nanmin(lambdaH),
                                  np.nanmax(lambdaH)))




def make_waves(test=False):
    k = []
    l = []

    theta = []
    dz = 8

    for i in ladcp_bins:
        theta.append(iwc.horizontal_azimuth(Uprime[i,:], Vprime[i,:], dz,\
                                            wl_min=wl_min,
                                            wl_max=wl_max,
                                            nfft=1024))
    theta = np.vstack(theta)

    k = kh*np.cos(theta)
    l = kh*np.sin(theta)

    k_wave = (2*np.pi)/k/1000
    l_wave = (2*np.pi)/l/1000

    waves = []
    depthsrev = np.tile(depths, [21,1])
    m = (-2*np.pi)/np.mean([wl_min, wl_max])

    for i in range(len(depthsrev)):
        waves.append(rt.wave(k=k.flatten(order='F')[i],
                             l=l.flatten(order='F')[i],
                             m=m, w0=omega.flatten(order='F')[i],
                             z0=depthsrev[i]))

    if test:
        i=7
        depthsrev = np.tile(depths, [21,1])
        m = (-2*np.pi)/np.mean([wl_min, wl_max])
        waves = []
        waves.append(rt.wave(k=k.flatten(order='F')[i],
                             l=l.flatten(order='F')[i],
                             m=m, w0=omega.flatten(order='F')[i],
                             z0=depthsrev[i]))
        duration = 48
        status = 12
        runcheck = []
        for i, wave in enumerate(waves):
            if np.isfinite(wave.k_init) and np.isfinite(wave.w0_init):
                wave.back3d(duration=duration, status=status,
                            print_run_report=True)
                runcheck.append('OK')
            else:
                runcheck.append('Bad Bin')
            print('Completed {}/{} bins'.format(i+1, len(waves)))

        fig1 = plt.figure()
    #    _, Ugrid = np.meshgrid(waves[7].x_ray, waves[7].U)
    #    plt.contourf(waves[7].x_ray.flatten(), waves[7].flow_grid, Ugrid)
        for i, flag in enumerate(runcheck):
            if flag == 'OK':
                waves[i].x_m_plot(fig=fig1, contour=True)

        plt.gca().invert_yaxis()




    duration = 24
    status = 12
    tstep = 10
    runcheck = []
    for i, wave in enumerate(waves):
        if np.isfinite(wave.k_init) and np.isfinite(wave.w0_init):
            wave.back3d(duration=duration,
                        tstep=tstep, status=status,
                        print_run_report=True)
            runcheck.append('OK')
        else:
            runcheck.append('Bad Bin')
        print('Completed {}/{} bins'.format(i+1, len(waves)))



    fig1 = plt.figure()
    _, Ugrid = np.meshgrid(waves[7].x_ray, waves[7].U)
    plt.contourf(waves[7].x_ray.flatten(), waves[7].flow_grid, Ugrid)
    cbar_flag = True
    for i, flag in enumerate(runcheck):
        if flag == 'OK':
            if cbar_flag:
                waves[i].plot_flow_x(fig=fig1)
                cbar_flag = None
            else:
                waves[i].plot_flow_x(fig=fig1)

    plt.gca().invert_yaxis()


    trace = go.Scatter(x=waves[7].x_ray.flatten(),
                       y=waves[7].z_ray.flatten(),
                       name = """
Wavenumbers:
   <br> k: {} <br>
    <br> l: {} <br>
<br> Frequency : {} <br>
<br> starting depth: {}  <br>""".format(waves[7].k_init, waves[7].k_init,
    waves[7].w0_init, waves[7].z_init))

    data = [trace]

    py.plot(data)



#
#for i in ladcp_bins:
#    k_i, l_i = iwc.horizontal_wave_vector_decomposition(Uprime[i,:],\
#                                                        Vprime[i,:],
#                                                        axis=0, nfft=1024)
#    k.append(k_i)
#    l.append(v_i)
