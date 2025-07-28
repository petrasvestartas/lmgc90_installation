
// File: wrap__PT3Dx_8h.xml

%feature("docstring") PT3Dx_LoadTactors "

load PT3Dx from RBDY3 and initialize existing_entites  

python usage : PT3Dx_LoadTactors()  
";

%feature("docstring") PT3Dx_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = PT3Dx_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PT3Dx_GetNbPT3Dx "

Get the number of PT3Dx.  

python usage : nb_PT3Dx = PT3Dx_GetNbPT3Dx()  

Returns
-------
nb_PT3Dx (integer) : the number of PT3Dx  
";

%feature("docstring") PT3Dx_SetDisplayRadius "

set the size of the glyph representing the PT3Dx  

python usage : PT3Dx_SetDisplayRadius(radius)  

Parameters
----------
* `radius` :  
    (double): radius of the PT3Dx contactors  
";

%feature("docstring") PT3Dx_GetPtrPT3Dx2BDYTY "

return a pointer onto the map pt3dx2bdyty  

python usage : pt3dx2bdyty = PT3Dx_GetPtrPT3Dx2BDYTY()  

Returns
-------
pt3dx2bdyty (integer array) : reference on map between pt3dx rank and body rank  
";

%feature("docstring") PT3Dx_InitOutlines "

Get a reference on the outlines of all PT3Dx.  

usage : outlines = PT3Dx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PT3Dx  
";

%feature("docstring") PT3Dx_InitScalarFields "

Get a reference on the scalar fields of all PT3Dx.  

usage : scalarfields = PT3Dx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PT3Dx array  
";

%feature("docstring") PT3Dx_UpdatePostdata "

Update values of outlines_PT3Dx and scalarfields_PT3Dx pointers.  

usage : PT3Dx_UpdatePostdata  
";

%feature("docstring") PT3Dx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PT3Dx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PT3Dx  
";

%feature("docstring") PT3Dx_GetNbScalarFields "

Get the number of scalar fields of a PT3Dx.  

python usage : nb_scalarfields = PT3Dx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PT3Dx  
";

%feature("docstring") PT3Dx_GetPtrAllConnectivities "

Get a reference on the connectivities of all PT3Dx.  

usage : connec = PT3Dx_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") PT3Dx_CleanMemory "

Free all memory allocated within PT3Dx module.  

python usage : PT3Dx_CleanMemory()  
";

