## Ray tracing with satGEM build notes:

        pandoc --variable geometry:margin=2cm -s ray_tracing_w_satGEM_notes.md -o ray_tracing.pdf
### Loading satGEM:
* Since the satGEM files are large (>10 gb), use the h5py python library which loads the files as objects which can be indexed and searched without actually loading the whole file. 
* all of the satGEM files are loaded into a single object defined in the ray tracing module. This makes accessing the data easier because it is all centralized into that one object without ever having to load the whole data set. 
* the velocity data grids are offset by one cell so the **U** grid corresponds to the normal lat and **centerlon**
* the **V** grid corresponds to normal lon and **centerlat**

 
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

**SETUP**

#### x, y, and z steps:

$$ \Delta x  = \Delta t \bigg [ \bigg (\frac{N^2 k}{m^2 \Omega} \bigg )  + U(x,y,z,t) \bigg ] $$

$$ \Delta y  = \Delta t \bigg [ \bigg ( \frac{N^2 l}{m^2 \Omega} \bigg ) + V(x,y,z,t) \bigg ] $$

$$ \Delta z  = \Delta t \bigg [ \frac{N^2(k^2 + l^2)}{m^3 \Omega}   \bigg ] $$

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