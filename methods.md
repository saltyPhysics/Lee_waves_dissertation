---
title: "Methods"
date: "May 2017"
author: "Manish Devana"
---
## Spectral analysis

Wave properties are inferred here from internal wave energy densities calculated from velocity and buoyancy perturbations using the parameterizations derived above. Using the dispersion relations of internal waves from internal wave theory. Estimation of intrinsic frequency and wavenumber components can help separate lee waves from other types of internal waves. Velocity and density measurements are separated into mean and perturbation components to isolate the wave driven effects.
$$ U = \bar u + u' $$
$$ \rho = \bar \rho + \rho'$$

Kinetic and potential energy densities are estimated using the power spectral density of the perturbations integrated between a range of vertical wavelengths, following the work of several previous internal wave investigations [@Waterman2013; @Meyer2016]. The Fourier transforms are performed on half overlapping segments of 1024 meters. Integration limits are determined qualitatively through examining the velocity perturbation profiles and identifying coherent structures. Fig (velocity profiles) shows that along the transect, these coherent structures have markedly different vertical sizes. Therefore, calculations are repeated several times with varying integration limits, Fourier transform sizes, and segment lengths. Hanning windows are applied to reduce variance loss. Kinetic and Potential energies are estimated as follows.
$$ E_{total} = \frac{1}{2} \bigg[ \rho \big(\langle u'^2 \rangle + \langle v'^2 \rangle ) + \rho N^2 \langle \eta^2 \rangle )\bigg]  $$

### Kinetic Energy.
$u$ and $v$ mean components are estimated using a sliding 2^nd^ order polynomial and subtracted from the measured flow. The polynomial fit starts with a 100 meter vertical segment which increases by 100 meters every 8 meters, or one step on the normalized depth grid. The power spectral density of these profiles are calculated for each bin and integrated between a chosen wavelength band. Kinetic energy per bin is calculated as
$$ KE =\frac{1}{2}\rho \int_{m_{c}}^{m_{0}}\langle u'^2 \rangle + \langle v'^2 \rangle dm $$

where $m_0$ and $m_c$ are the vertical wavenumber integration limits.
The ranges of vertical wavelengths utilized were 200-500 and 500-1000 meters.

### Potential Energy
Potential energy is estimated from isopycnal displacement $\eta$. Strain, $\xi$, is the depth-wise gradient of isopycnal displacement which is calculated following Bray and Fofonoff's adiabatic leveling method .


## **References**
