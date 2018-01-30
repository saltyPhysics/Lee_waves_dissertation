# Dissertation Outline

## Methods
- Describe data - i.e. towyo transect ~110km long over complex shag rock topography
- Get Anomalies - (u’, v’, rho’) 
  - 2nd order polynomial sliding polynomial fit subtracted from each profile (subtracting mean flow). 
  - bandpass filter applied to filter out signals with anything shorter than 300m vertical wavelengths
  - Data binned into half overlapping vertical segments each 1024 meters long.
  - Power spectral density (welch method with hanning windowing) was used to estimate the energy (integrated between target vertical wavelengths 
    - sensitivity to integration limits tested by varying limits (little effect found (quantify this))


- Kinetic and Potential energies estimated as 
