
// File: wrap__deposit2D_8h.xml

%feature("docstring") deposit2D_Potential "

Computes a new deposit under potential with or without big particles.  

python call: coor = deposit2D_Potential(in_radii, lx, potential[, dradii,
dcoor])  

Parameters
----------
* `in_radii` :  
    (double array): given radii list (i.e. granulometry)  
* `lx` :  
    (double): width of the box in which to deposit  
* `potential` :  
    (integer): for deposit (1->gravity, 2->wall, 3->big_particles)  
* `dradii` :  
    (double array) (optional) : a list of already deposited radii  
* `dcoor` :  
    (double array) (optional) : a list of already deposited coor (must be of
    size [nb_dradii,3])  

Returns
-------
coor (double array): coordinates of deposited radii (shape [nb_radii,2]) PYDOC  
";

