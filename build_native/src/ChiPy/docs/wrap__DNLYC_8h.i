
// File: wrap__DNLYC_8h.xml

%feature("docstring") DNLYC_LoadTactors "

load DNLYC from RBDY3 and initialize existing_entites  

python usage : DNLYC_LoadTactors()  
";

%feature("docstring") DNLYC_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = DNLYC_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") DNLYC_GetNbDNLYC "

Get the number of DNLYC.  

python usage : nb_DNLYC = DNLYC_GetNbDNLYC()  

Returns
-------
nb_DNLYC (integer) : the number of DNLYC  
";

%feature("docstring") DNLYC_GetPtrDNLYC2BDYTY "

return a pointer onto the map dnlyc2bdyty  

python usage : dnlyc2bdyty = DNLYC_GetPtrDNLYC2BDYTY()  

Returns
-------
dnlyc2bdyty (integer array) : reference on map between dnlyc rank and body rank  
";

%feature("docstring") DNLYC_InitOutlines "

Get a reference on the outlines of all DNLYC.  

usage : outlines = DNLYC_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_DNLYC  
";

%feature("docstring") DNLYC_InitScalarFields "

Get a reference on the scalar fields of all DNLYC.  

usage : scalarfields = DNLYC_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_DNLYC array  
";

%feature("docstring") DNLYC_UpdatePostdata "

Update values of outlines_DNLYC and scalarfields_DNLYC pointers.  

usage : DNLYC_UpdatePostdata  
";

%feature("docstring") DNLYC_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = DNLYC_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
DNLYC  
";

%feature("docstring") DNLYC_GetNbScalarFields "

Get the number of scalar fields of a DNLYC.  

python usage : nb_scalarfields = DNLYC_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a DNLYC  
";

%feature("docstring") DNLYC_GetPtrAllConnectivities "

Get a reference on the connectivities of all DNLYC.  

usage : connec = DNLYC_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") DNLYC_CleanMemory "

Free all memory allocated within DNLYC module.  

python usage : DNLYC_CleanMemory()  
";

