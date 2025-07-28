
// File: wrap__deposit3D_8h.xml

%feature("docstring") deposit3D_InContainer "

Computes a new deposit under gravity in a container.  

i_shape = 0 : box  

*   a point (x, y, z) is in the box iff x is in [-lx/2, lx/2], y is in [-ly/2,
    ly/2] and z is in [0, lz] i_shape = 1 : cylinder  
*   a point (x, y, z) is in the cylinder iff x^2 + y^2 is in [0, R^2] and z is
    in [0, lz] i_shape = 2 : sphere  
*   a point (x, y, z) is in the sphere iff x^2 + y^2 + z^2 is in [0, R^2]  

python call: radii, coor = deposit3D_InContaier(in_radii, shape, p1, p2, p3[,
dradii, dcoor, seed, with_log])  

Parameters
----------
* `in_radii` :  
    (double array): given radii list (i.e. granulometry)  
* `shape(integer)` :  
    of container (0->box, 1->cylinder, 2->sphere)  
* `p1` :  
    (double): box-> lx, cylinder->R, sphere->R  
* `p2` :  
    (double): box-> ly, cylinder->lz, sphere->ignored  
* `p3` :  
    (double): box-> lz, cylinder->ignored, sphere->ignored  
* `dradii` :  
    (double array) (optional) : a list of already deposited radii  
* `dcoor` :  
    (double array) (optional) : a list of already deposited coor (must be of
    size [nb_dradii,3])  
* `seed` :  
    (integer array) (optional) : an input seed to control randomness  
* `with_log(integer)` :  
    de/activate log message  

Returns
-------
radii (double array): list of deposited radii coor (double array): coordinates
of deposited radii (shape [nb_radii,3]) PYDOC  
";

