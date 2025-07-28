
// File: wrap__CYLND_8h.xml

%feature("docstring") CYLND_LoadTactors "

load CYLND from RBDY3 and initialize existing_entites  

python usage : CYLND_LoadTactors()  
";

%feature("docstring") CYLND_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = CYLND_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") CYLND_GetNbCYLND "

Get the number of CYLND.  

python usage : nb_CYLND = CYLND_GetNbCYLND()  

Returns
-------
nb_CYLND (integer) : the number of CYLND  
";

%feature("docstring") CYLND_GetShape "

Get the shape of a CYLND.  

usage : shape = CYLND_GetShape(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of CYLND  

Returns
-------
shape (double array) : axis length of the CYLND  
";

%feature("docstring") CYLND_GetPtrCYLND2BDYTY "

return a pointer onto the map cylnd2bdyty  

python usage : cylnd2bdyty = CYLND_GetPtrCYLND2BDYTY()  

Returns
-------
cylnd2bdyty (integer array) : reference on map between cylnd rank and body rank  
";

%feature("docstring") CYLND_InitOutlines "

Get a reference on the outlines of all CYLND.  

python usage : outlines = CYLND_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_CYLND  
";

%feature("docstring") CYLND_InitScalarFields "

Get a reference on the scalar fields of all CYLND.  

python usage : scalarfields = CYLND_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_CYLND array  
";

%feature("docstring") CYLND_UpdatePostdata "

Update values of outlines_CYLND and scalarfields_CYLND pointers.  

python usage : CYLND_UpdatePostdata  
";

%feature("docstring") CYLND_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = CYLND_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
CYLND  
";

%feature("docstring") CYLND_GetNbScalarFields "

Get the number of scalar fields of a CYLND.  

python usage : nb_scalarfields = CYLND_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a CYLND  
";

%feature("docstring") CYLND_GetPtrAllConnectivities "

Get a reference on the connectivities of all CYLND.  

python usage : connec = CYLND_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") CYLND_CleanMemory "

Free all memory allocated within CYLND module.  

python usage : CYLND_CleanMemory()  
";

