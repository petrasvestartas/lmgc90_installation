
// File: wrap__DISKx_8h.xml

%feature("docstring") DISKx_LoadTactors "

load DISKx from RBDY2 file and initialize existing_entites  

python usage : DISKx_LoadTactors()  
";

%feature("docstring") DISKx_GetNbDISKx "

Get the number of DISKx in the container.  

python usage : nb_diskx = DISKx_GetNbDISKx()  

Returns
-------
nb_DISKx (integer) : the number of DISKx in container  
";

%feature("docstring") DISKx_GetDISKx2BDYTY "

Get a copy of map DISKx2bdyty.  

usage : polyr2bdyty = DISKx_GetDISKx2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") DISKx_GetPtrDISKx2BDYTY "

return a pointer onto the map diskx2rbdy2  

python usage : diskx2bdyty = DISKx_GetPtrDISKx2BDYTY()  

Returns
-------
diskx2bdyty (integer array) : reference on map between diskx rank and body rank  
";

%feature("docstring") DISKx_IsVisible "

return if a body visible  

python usage : visible = DISKx_IsVisible(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of DISKx  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") DISKx_GetContactorRadius "

Get the radius of a given DISKx.  

python usage : radius = DISKx_GetContactorRadius(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of a DISKx (in the list of all the DISKx)  

Returns
-------
radius (double) : the radius of the DISKx of rank itact  
";

%feature("docstring") DISKx_GetMeanRadius "

Get the mean radius of DISKx in the container.  

python usage : radius = DISKx_GetMeanRadius()  

Returns
-------
radius (double) : the mean radius of DISKx in the container  
";

%feature("docstring") DISKx_GetMaxRadius "

Get the max radius of DISKx in the container.  

python usage : radius = DISKx_GetMaxRadius()  

Returns
-------
radius (double) : the max radius of DISKx in the contactor  
";

%feature("docstring") DISKx_GetMinRadius "

Get the min radius of DISKx in the container.  

python usage : radius = DISKx_GetMinRadius()  

Returns
-------
radius (double) : the min radius of DISKx in the container  
";

%feature("docstring") DISKx_GetContactorColor "

Get the color of a given DISKx.  

python usage : color = DISKx_GetContactorColor(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of a DISKx  

Returns
-------
color (string) : the color of the DISKx itact  
";

%feature("docstring") DISKx_GetRadius "

get radius of a DISKx  

python usage : radius = DISKx_GetRadius(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of DISKx  

Returns
-------
radius (double) : the radius of DISKx of body ibdyty  
";

%feature("docstring") DISKx_GetContactorCoor "

get coordinates of the center of a given DISKx  

python usage : vector = DISKx_GetContactorCoor(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") DISKx_InitOutlines "

Get a reference on the outlines of all DISKx.  

python usage : outlines = DISKx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_DISKx  
";

%feature("docstring") DISKx_InitScalarFields "

Get a reference on the scalar fields of all DISKx.  

python usage : scalarfields = DISKx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_DISKx array  
";

%feature("docstring") DISKx_UpdatePostdata "

Update values of outlines_DISKx and scalarfields_DISKx pointers.  

python usage : DISKx_UpdatePostdata()  
";

%feature("docstring") DISKx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = DISKx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
DISKx  
";

%feature("docstring") DISKx_GetNbScalarFields "

Get the number of scalar fields of a DISKx.  

python usage : nb_scalarfields = DISKx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a DISKx  
";

%feature("docstring") DISKx_CleanMemory "

Free all memory allocated within DISKx module.  

python usage : DISKx_CleanMemory()  
";

%feature("docstring") DISKx_SetXdilation "

set increase of radius of a DISKx due to expansion  

python usage : DISKx_SetXdilation(itacty,x)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  
* `x` :  
    (float) : increase of radius  
";

%feature("docstring") DISKx_SetVdilation "

set increase rate of radius of a DISKx due to expansion  

python usage : DISKx_SetVdilation(itacty, v)  

Parameters
----------
* `itacty` :  
    (integer) : rank of contactor  
* `v` :  
    (float) : radius increase rate  
";

