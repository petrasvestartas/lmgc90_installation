
// File: wrap__PLANx_8h.xml

%feature("docstring") PLANx_LoadTactors "

load PLANx from RBDY3 and initialize existing_entites  

python usage : PLANx_LoadTactors()  
";

%feature("docstring") PLANx_GetNbPLANx "

Get the number of PLANx.  

python usage : nb_PLANx = PLANx_GetNbPLANx()  

Returns
-------
nb_PLANx (integer) : the number of PLANx  
";

%feature("docstring") PLANx_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = PLANx_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PLANx_GetPtrPLANx2BDYTY "

return a pointer onto the map planx2bdyty  

python usage : planx2bdyty = PLANx_GetPtrPLANx2BDYTY()  

Returns
-------
planx2bdyty (integer array) : reference on map between planx rank and body rank  
";

%feature("docstring") PLANx_InitOutlines "

Get a reference on the outlines of all PLANx.  

usage : outlines = PLANx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PLANx  
";

%feature("docstring") PLANx_InitScalarFields "

Get a reference on the scalar fields of all PLANx.  

usage : scalarfields = PLANx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PLANx array  
";

%feature("docstring") PLANx_UpdatePostdata "

Update values of outlines_PLANx and scalarfields_PLANx pointers.  

usage : PLANx_UpdatePostdata  
";

%feature("docstring") PLANx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PLANx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PLANx  
";

%feature("docstring") PLANx_GetNbScalarFields "

Get the number of scalar fields of a PLANx.  

python usage : nb_scalarfields = PLANx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PLANx  
";

%feature("docstring") PLANx_GetPtrAllConnectivities "

Get a reference on the connectivities of all PLANx.  

usage : connec = PLANx_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") PLANx_CleanMemory "

Free all memory allocated within PLANx module.  

python usage : PLANx_CleanMemory()  
";

