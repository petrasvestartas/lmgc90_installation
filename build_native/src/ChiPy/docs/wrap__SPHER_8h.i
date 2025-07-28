
// File: wrap__SPHER_8h.xml

%feature("docstring") SPHER_LoadTactors "

load SPHER from RBDY3 and initialize existing_entites  

python usage : SPHER_LoadTactors()  
";

%feature("docstring") SPHER_SetRadiusCorrection "

set a radius correction  

python usage : SPHER_SetRadiusCorrection(corr)  

Parameters
----------
* `corr` :  
    (real) :  
";

%feature("docstring") SPHER_GetNbSPHER "

Get the number of SPHER.  

python usage : nb_SPHER = SPHER_GetNbSPHER()  

Returns
-------
nb_SPHER (integer) : the number of SPHER  
";

%feature("docstring") SPHER_GetSPHER2BDYTY "

Get a copy of map SPHER2bdyty.  

usage : polyr2bdyty = SPHER_GetSPHER2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") SPHER_GetPtrSPHER2BDYTY "

return a pointer onto the map spher2bdyty  

python usage : spher2bdyty = SPHER_GetPtrSPHER2BDYTY()  

Returns
-------
spher2bdyty (integer array) : reference on map between spher rank and body rank  
";

%feature("docstring") SPHER_GetContactorRadius "

Get the radius of a SPHER contactor.  

python usage : radius = SPHER_GetContactorRadius(itact)  

Parameters
----------
* `itact` :  
    (integer) : id of a SPHER  

Returns
-------
radius (double) : the radius of the SPHER number itact  
";

%feature("docstring") SPHER_GetContactorCoor "

get coordinates of the center of a given SPHER  

usage : vector = SPHER_GetContactorCoor(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") SPHER_GetContactorCoorb "

get coordinates at the begin of the time step of the center of a given SPHER  

usage : vector = SPHER_GetContactorCoorb(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") SPHER_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = SPHER_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") SPHER_InitOutlines "

Get a reference on the outlines of all SPHER.  

usage : outlines = SPHER_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_SPHER  
";

%feature("docstring") SPHER_InitScalarFields "

Get a reference on the scalar fields of all SPHER.  

usage : scalarfields = SPHER_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_SPHER array  
";

%feature("docstring") SPHER_UpdatePostdata "

Update values of outlines_SPHER and scalarfields_SPHER pointers.  

usage : SPHER_UpdatePostdata  
";

%feature("docstring") SPHER_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = SPHER_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
SPHER  
";

%feature("docstring") SPHER_GetNbScalarFields "

Get the number of scalar fields of a SPHER.  

python usage : nb_scalarfields = SPHER_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a SPHER  
";

%feature("docstring") SPHER_GetPtrAllConnectivities "

Get a reference on the connectivities of all SPHER.  

usage : connec = SPHER_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") SPHER_CleanMemory "

Free all memory allocated within SPHER module.  

python usage : SPHER_CleanMemory()  
";

