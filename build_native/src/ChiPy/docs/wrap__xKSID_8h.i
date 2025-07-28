
// File: wrap__xKSID_8h.xml

%feature("docstring") xKSID_LoadTactors "

load xKSID from RBDY2 and initialize existing_entites  

python usage : xKSID_LoadTactors()  
";

%feature("docstring") xKSID_GetNbxKSID "

Get the number of xKSID in the container.  

python usage : nb_diskx = xKSID_GetNbxKSID()  

Returns
-------
nb_xKSID (integer) : the number of xKSID in container  
";

%feature("docstring") xKSID_GetPtrxKSID2BDYTY "

return a pointer onto the map xksid2rbdy2  

python usage : xksid2rbdy2 = xKSID_GetPtrxKSID2BDYTY()  

Returns
-------
xksid2rbdy2 (integer array) : reference on map between xksid rank and body/tact
rank  
";

%feature("docstring") xKSID_IsVisible "

return if a body visible  

usage : visible = xKSID_IsVisible(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of xKSID  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") xKSID_GetContactorRadius "

Get the radius of a given xKSID.  

python usage : radius = xKSID_GetContactorRadius(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of a xKSID (in the list of all the xKSID)  

Returns
-------
radius (double) : the radius of the xKSID of rank itact  
";

%feature("docstring") xKSID_GetContactorCoor "

get coordinates of the center of a given xKSID  

usage : vector = xKSID_GetContactorCoor(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") xKSID_InitOutlines "

Get a reference on the outlines of all xKSID.  

usage : outlines = xKSID_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_xKSID  
";

%feature("docstring") xKSID_InitScalarFields "

Get a reference on the scalar fields of all xKSID.  

usage : scalarfields = xKSID_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_xKSID array  
";

%feature("docstring") xKSID_UpdatePostdata "

Update values of outlines_xKSID and scalarfields_xKSID pointers.  

usage : xKSID_UpdatePostdata  
";

%feature("docstring") xKSID_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = xKSID_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
xKSID  
";

%feature("docstring") xKSID_GetNbScalarFields "

Get the number of scalar fields of a xKSID.  

python usage : nb_scalarfields = xKSID_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a xKSID  
";

%feature("docstring") xKSID_CleanMemory "

Free all memory allocated within xKSID module.  

python usage : xKSID_CleanMemory()  
";

%feature("docstring") xKSID_SetXdilation "

set increase of radius of a xKSID due to expansion  

python usage : xKSID_SetXdilation(itacty,x)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  
* `x` :  
    (float) : increase of radius  
";

%feature("docstring") xKSID_SetVdilation "

set increase rate of radius of a xKSID due to expansion  

python usage : xKSID_SetVdilation(itacty, v)  

Parameters
----------
* `itacty` :  
    (integer) : rank of contactor  
* `v` :  
    (float) : radius increase rate  
";

