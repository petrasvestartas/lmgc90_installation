
// File: wrap__mbs2D_8h.xml

%feature("docstring") MBS2D_setNb "

Set the number of MBS.  

python usage : MBS2D_setNb(nb)  

Parameters
----------
* `nb` :  
    (integer) : set the number of MBS  
";

%feature("docstring") MBS2D_getNb "

Get the number of MBS.  

python usage : nb = MBS2D_getNb()  

Returns
-------
nb (integer) : the number of MBS  
";

%feature("docstring") MBS2D_setNbNodes "

Set the number of nodes of a MBS.  

python usage : MBS2D_setNbNodes(ibdyty, nb)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  
* `nb` :  
    (integer) : the number of nodes of the MBS  
";

%feature("docstring") MBS2D_setNbTactors "

Set the number contactors of a MBS.  

python usage : MBS2D_setNbTactors(ibdyty, nb)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  
* `nb` :  
    (integer) : the number of contactor of the MBS  
";

%feature("docstring") MBS2D_getPtrCoor "

Get a pointer on the coor of a MBS.  

python usage : coor = MBS2D_getPtrCoor(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered MBS  

Returns
-------
coor (double 2D-array) : reference on the coordinates of the nodes  
";

%feature("docstring") MBS2D_getPtrCoorTT "

Set the array of coordinates of nodes of a MBS.  

python usage : coorTT = MBS2D_GetPtrCoorTT(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : id of the MBS  

Returns
-------
coorTT (double array) : coordinates of nodes of a MBS  
";

%feature("docstring") MBS2D_addContactor "

Add a new contactor to a MBS.  

Available contactor types are :  

*   JONCx: inputs are:
    -   rdata must hold [axe_x, axe_y]  
*   POLYR: inputs are:
    -   rdata must hold the coordinates of the vertices [x_1, y_1, ... x_n, y_n]  
    -   idata must hold the number of vertices  

python usage : MBS2D_addContactor(ibdyty, inodty, itacty, tacttype, color,
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

%feature("docstring") MBS2D_initialize "

Initialize MBS module once loading is done.  

python usage : MBS2D_initialize()  
";

%feature("docstring") MBS2D_finalize "

Finalize MBS module.  

python usage : MBS2D_finalize()  
";

%feature("docstring") MBS2D_IncrementStep "

compute the current velocity and displacement  

python usage : MBS2D_IncrementStep()  
";

%feature("docstring") MBS2D_ComputeFreeVelocity "

compute free velocity  

python usage : MBS2D_ComputeFreeVelocity()  
";

%feature("docstring") MBS2D_ComputeDof "

update current position and velocity  

python usage : MBS2D_ComputeDof()  
";

%feature("docstring") MBS2D_UpdateDof "

save d.o.f. of the end of the time step to d.o.f. of the begining of the next
one  

python usage : MBS2D_UpdateDof()  
";

