## Ray tracing with satGEM build notes:

        pandoc --variable geometry:margin=2cm -s ray_tracing.md -o ray_tracing.pdf
### Loading satGEM:
* Since the satGEM files are large (>10 gb), use the h5py python library which loads the files as objects which can be indexed and searched without actually loading the whole file. 
* all of the satGEM files are loaded into a single object defined in the ray tracing module. This makes accessing the data easier because it is all centralized into that one object without ever having to load the whole data set. 
* the velocity data grids are offset by one cell so the **U** grid corresponds to the normal lat and **centerlon**
* the **V** grid corresponds to normal lon and **centerlat**


### Interpolating and Extrapolating through satGEM and bathymetry
All the gradients are precalculated (dndx, dudx, dvdz, etc...) and saved as an additional dataset along side satGEM files (with corresponding modified grids from gradient calculations).

Interpolating through satGEM data is very sensitive to the method. The 3 main approaches have been to use a linear radial basis function, linear ND interpolation, and nearest neighbor interpolation. 

* Nearest Neighbor interpolation is by far the fastest but also least accurate
* radial basis function works but fails at extrapolating 
* linear ND gives more consistent results than radial basis function but waaaay slower

**FINAL METHOD**: Use linear ND on a subset of the dataset (adjustable) and in regions below satGEM depth, switches to nearest neighbor. This is probably ok since below 3000 meters, conditions are less likely to vary significantly 

### accessing satGEM while running ray tracing:
* the satGEM object has a locate function which uses a distance formul along each axis to get the index along each axis. 
    * **only the indices are returned** 
* the indices are tested against previous steps indices, and if they haven't changed there is no need to reload the data. 


-------------------------------------------------------------------------------------------------------------------------------------------

## Ray Tracing Equations:

- **These are adapted from the buhler and olbers ray tracing papers to work in 4 dimensions with no assumptions about flow or the density field in time and space**

- Go back and look at the original derivations from the **Gill** book. How the equations used in ray tracing papers starts to make sense when you look at how aspect ratio is substituted into the dispersion relation
    - When doing this, aspect can be assummed small and negligible compared to m because of the order of magntitude difference between them (validated in my observations) .
    - What happens if I use the full equations without doing this?
    - **what is the fast marching method and eikonal equation .... Seems super important to solving ray tracing**
    - When loading grid in for N2, calculations instead of loading the whole box, load in a box (4d) around the rays current location which has enough buffer to make sure the ray doesnt move out of it. Also, add a check to say if the ray is getting close to the edge of the box, re center the box. **(see diagram)**

\pagebreak

## Ray Equations


### Group Speeds

$$ C_gx = \frac{k m^2 (N^2 - f^2)}{\omega(k^2 + l^2 + m^2)^2} + U(x,y,z,t) $$

$$ C_gy = \frac{l m^2 (N^2 - f^2)}{\omega(k^2 + l^2 + m^2)^2} + V(x,y,z,t) $$

$$ C_gz = \frac{-m (k^2 + l^2) (N^2 - f^2)}{\omega(k^2 + l^2 + m^2)^2} $$

### Wave Property Changes
$$ \frac{\partial k}{\partial t}  = -\bigg[\frac{N (k^2 +l^2)}{\omega(k^2 +l^2 + m^2)} \frac{\partial N}{\partial x }\bigg] - k\frac{\partial U}{\partial x } - l\frac{\partial V}{\partial x } $$

$$ \frac{\partial l}{\partial t} = -\bigg[\frac{N (k^2 +l^2)}{\omega(k^2 +l^2 + m^2)} \frac{\partial N}{\partial y }\bigg]  - k\frac{\partial U}{\partial y } - l\frac{\partial V}{\partial y }$$

$$ \frac{\partial m}{\partial t} = -\bigg[\frac{N (k^2 +l^2)}{\omega(k^2 +l^2 + m^2)} \frac{\partial N}{\partial z } \bigg] - k\frac{\partial U}{\partial z } - l\frac{\partial V}{\partial z }$$

$$ \frac{\partial \omega}{\partial t} = \bigg[\frac{N (k^2 +l^2)}{\omega(k^2 +l^2 + m^2)}\bigg] \bigg[\frac{\partial N}{\partial x } + \frac{\partial N}{\partial y }+  \frac{\partial N}{\partial z } \bigg] + k\frac{\partial U}{\partial t } + l\frac{\partial V}{\partial t }$$

### Plane Wave Fitting
For all relevant parameters ($u',v',w',\rho', p'$), a plane wave solution is assumed in the form:

$$ X' = X_0cos(kx +ly +mz - \omega t) $$


### Wave Action and Momentum Fluxes 
Wave action is the ratio of energy density to frequency which remains constant along the ray path. Consequently, changes in local frequency must result in directly proportional changes in energy density. The changes in energy density represent momentum fluxes between a ray _tube_ and the mean flow.


Waves transport momentum via reynolds stresses in the form:

$$ \rho\overline{u'w'} $$ 
**which represents vertical flux of horizontal momentum  in x direction**

The stress on mean flow (work done on mean flow) is :

$$ -\rho \overline{u'w'} \frac{\partial  U}{\partial z} $$  


$$ A = \frac{E}{\omega} = \frac{\frac{1}{2}\rho(u'^2 + v'^2 + w'^2) + \frac{1}{2}\rho(N^2 \overline{b'^2})}{\omega_{r}}  $$

$$ \frac{\partial A}{\partial t} =  $$ 

$$ \frac{dE}{dt} = -E\bigg[\frac{\partial  C_gx}{\partial x} +\frac{\partial  C_gy}{\partial y} +\frac{\partial  C_gz}{\partial z} \bigg] - \mathbf{\frac{E}{\omega}\frac{\partial \omega}{\partial t}} $$

$$ \frac{dE}{dt} = -E\bigg[\frac{\partial  C_gx}{\partial x} +\frac{\partial  C_gy}{\partial y} +\frac{\partial  C_gz}{\partial z}\bigg] - \mathbf{A\frac{\partial \omega}{\partial t}} $$

The second term describes the momentum fluxes between the mean flow and wave -  it works out to be momentum flux per time (N/m^2) /s so can I integrate this over time to get the total momentum flux along the ray path?

modify equation to solve for reynolds stress term

$$ A\frac{\partial \omega}{\partial t} = \frac{\partial F}{\partial t} $$

Where $F$ is the Momentum flux with mean flow 

$$ \frac{\partial F}{\partial t}=  -E\bigg[\frac{\partial  C_gx}{\partial x} +\frac{\partial  C_gy}{\partial y} +\frac{\partial  C_gz}{\partial z}\bigg] - \frac{dE}{dt}   $$ 

$$  \frac{\partial F}{\partial t}  = \frac{kg}{ms^3} = \frac{N/m^2}{s} $$ 

$$ \int_{0}^{t_f} \frac{\partial F}{\partial t} dt = F  $$



Where $F$ is the total  flux of  momentum between mean flow and the (from $t=0$ to $t_f$ where wave either propogation ends)
wave. 

Through ray tracing, the total time of wave propagation, and changes to 
#### Get change in position and new field values: _Inverse Haversine Forumala for spherical distance_

$$ \Delta longitude  = (\Delta x)(**conversion factor**)  $$

$$ \Delta latitude  = (\Delta y)(111.11 km)  $$eff


#### Changes in wave parameters:



$$ \Delta k  =   \Delta t \bigg [  r_x  + k\frac{\partial U}{\partial x} + l\frac{\partial V}{\partial x} \bigg ] $$

$$ \Delta l  =   \Delta t \bigg [  r_y  + k\frac{\partial U}{\partial y} + l\frac{\partial V}{\partial y} \bigg ] $$

$$ \Delta m  =   \Delta t \bigg [  r_z  + k\frac{\partial U}{\partial z} + l\frac{\partial V}{\partial z} \bigg ] $$


where $r_{x,y,z}$ is the refraction index:

$$ r = \frac{N(k^2 +l^2)}{m^2 \Omega} \frac{\partial N}{\partial (x,y,z)} $$

#### Wave dispersion relation change

Changes to intrinsic frequency/dispersion relation along the wave ray is equal to the sum of the changes to the wave number components :

$$ \Delta \Omega = \Delta t \bigg [ (r_x + r_y + r_z) + k\Delta U + l\Delta V \bigg ]   $$


## Wave Action and Momentum fluxes
Wave action is conserved along the ray path but frequency and energy density can vary (proportionally). In regions where frequency increases (such as when wave moves into stronger flows), energy density increases proportionally , removing energy/momentum from the mean flow. When a wave moves into weaker flow and frequency decreases, energy is fed into the mean flow from the wave. This momentum flux acts as a drag on the mean flow (and is ultimately a reynolds stress between mean flow and a wave). In order to quantify the momentum flux over the entire lifetime of a lee wave (generation to breaking or suppression), the wave action flux is multiplied by the **WHATTTTTTT** to get the momentum flux along the group path.

One way to do it is directly calculate the u' and w' components via a plane wave fit but then what is the point of doing wave action calcuations? Additionally, regardless of whether using wave action or direction reynolds stress calculations, there needs to be integration limits - **WHICH ARE???** 
