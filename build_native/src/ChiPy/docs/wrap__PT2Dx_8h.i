
// File: wrap__PT2Dx_8h.xml

%feature("docstring") PT2Dx_LoadTactors "

load PT2Dx from RBDY2 and initialize existing_entites  

python usage : PT2Dx_LoadTactors()  
";

%feature("docstring") PT2Dx_GetNbPT2Dx "

Get the number of PT2Dx in the container.  

python usage : nb_pt2d = PT2Dx_GetNbPT2Dx()  

Returns
-------
nb_pt2d (integer) : the number of PT2Dx in container  
";

%feature("docstring") PT2Dx_SetDisplayRadius "

Set a radius to display a pt2dx.  

python usage : PT2Dx_SetDisplayRadius(radius)  

Parameters
----------
* `radius` :  
    (double) : value of the radius which should be used for display  
";

%feature("docstring") PT2Dx_GetPtrPT2Dx2BDYTY "

return a pointer onto the map pt2dx2rbdy2  

python usage : ptd2x2rbdy2 = PT2Dx_GetPtrPT2Dx2BDYTY()  

Returns
-------
pt2dx2rbdy2 (integer array) : reference on map between pt2dx rank and body/tact
rank  
";

%feature("docstring") PT2Dx_IsVisible "

return if a body visible  

usage : visible = PT2Dx_IsVisible(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of PT2Dx  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PT2Dx_InitOutlines "

Get a reference on the outlines of all PT2Dx.  

usage : outlines = PT2Dx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PT2Dx  
";

%feature("docstring") PT2Dx_InitScalarFields "

Get a reference on the scalar fields of all PT2Dx.  

usage : scalarfields = PT2Dx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PT2Dx array  
";

%feature("docstring") PT2Dx_UpdatePostdata "

Update values of outlines_PT2Dx and scalarfields_PT2Dx pointers.  

usage : PT2Dx_UpdatePostdata  
";

%feature("docstring") PT2Dx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PT2Dx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PT2Dx  
";

%feature("docstring") PT2Dx_GetNbScalarFields "

Get the number of scalar fields of a PT2Dx.  

python usage : nb_scalarfields = PT2Dx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PT2Dx  
";

%feature("docstring") PT2Dx_CleanMemory "

Free all memory allocated within PT2Dx module.  

python usage : PT2Dx_CleanMemory()  
";

