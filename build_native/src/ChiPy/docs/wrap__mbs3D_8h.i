
// File: wrap__mbs3D_8h.xml

%feature("docstring") MBS3D_setNb "

Set the number of MBS.  

python usage : MBS3D_setNb(nb)  

Parameters
----------
* `nb` :  
    (integer) : set the number of MBS  
";

%feature("docstring") MBS3D_getNb "

Get the number of MBS.  

python usage : nb = MBS3D_getNb()  

Returns
-------
nb (integer) : the number of MBS  
";

%feature("docstring") MBS3D_setNbNodes "

Set the number of nodes of a MBS.  

python usage : MBS3D_setNbNodes(ibdyty, nb)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  
* `nb` :  
    (integer) : the number of nodes of the MBS  
";

%feature("docstring") MBS3D_setNbTactors "

Set the number contactors of a MBS.  

python usage : MBS3D_setNbTactors(ibdyty, nb)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  
* `nb` :  
    (integer) : the number of contactor of the MBS  
";

%feature("docstring") MBS3D_getPtrCoor "

Get a pointer on the coor of a MBS.  

usage : coor = MBS3D_GetPtrCoor(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered MBS  

Returns
-------
coor (double 2D-array) : reference on the coordinates of the nodes  
";

%feature("docstring") MBS3D_getPtrCoorTT "

Set the array of coordinates of nodes of a MBS.  

python usage : coor = MBS3D_getPtrCoorTT(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  

Returns
-------
coor (double array) : coordinates of nodes of a MBS (in contact configuration)  
";

%feature("docstring") MBS3D_getPtrLocalFrame "

Get a pointer on the coor of a MBS.  

usage : frame = MBS3D_GetPtrLocalFrame(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered MBS  

Returns
-------
frame (double 2D-array) : local frame  
";

%feature("docstring") MBS3D_getPtrLocalFrameTT "

Set the array of coordinates of nodes of a MBS.  

python usage : frameTT = MBS3D_GetPtrLocalFrameTT(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  

Returns
-------
frameTT (double array) : local frame (in contact configuration)  
";

%feature("docstring") MBS3D_addContactor "

Add a new contactor to a MBS.  

Available contactor types are :  

*   PLANx: inputs are:
    -   rdata must hold [axe_x, axe_y, axe_z]  
*   POLYR: inputs are:
    -   rdata must hold the coordinates of the vertices [x_1, y_1, z_1, ... x_n,
        y_n, z_n]  
    -   idata must hold the connecivity of each triangle defining the surface  

python usage : MBS3D_addContactor(ibdyty, inodty, itacty, tacttype, color,
rdata, idata=None)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the MBS  
* `inodty` :  
    (integer) : rank of the node of the MBS the contactor is tied to  
* `itacty` :  
    (integer) : rank of the contactor of MBS  
* `tactype` :  
    (string [5]) : the type of contactor  
* `color` :  
    (string [5]) : the color of the contactor  
* `rdata` :  
    (double array) : the new value of the vector  
* `idata` :  
    (integer array) : the new value of the vector  
";

%feature("docstring") MBS3D_initialize "

Initialize MBS module once loading is done.  

python usage : MBS3D_initialize()  
";

%feature("docstring") MBS3D_finalize "

Finalize MBS module.  

python usage : MBS3D_finalize()  
";

%feature("docstring") MBS3D_IncrementStep "

compute the current velocity and displacement  

python usage : MBS3D_IncrementStep()  
";

%feature("docstring") MBS3D_ComputeFreeVelocity "

compute free velocity  

python usage : MBS3D_ComputeFreeVelocity()  
";

%feature("docstring") MBS3D_ComputeDof "

update current position and velocity  

python usage : MBS3D_ComputeDof()  
";

%feature("docstring") MBS3D_UpdateDof "

save d.o.f. of the end of the time step to d.o.f. of the begining of the next
one  

python usage : MBS3D_UpdateDof()  
";

