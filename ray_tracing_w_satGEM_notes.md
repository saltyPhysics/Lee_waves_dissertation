## Ray tracing with satGEM build notes:

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



