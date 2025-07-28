
// File: wrap__JONCx_8h.xml

%feature("docstring") JONCx_LoadTactors "

load JONCx from RBDY2 and initialize existing_entites  

python usage : JONCx_LoadTactors()  
";

%feature("docstring") JONCx_GetNbJONCx "

Get the number of JONCx in container.  

python usage : nb_joncx = JONCx_GetNbJONCx()  

Returns
-------
nb_joncx (integer) : the number of JONCx in container  
";

%feature("docstring") JONCx_GetBodyId "

Get the body rank of a given JONCx.  

python usage : ibdyty = JONCx_GetBodyId(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : JONCx rank  

Returns
-------
ibdyty (integer) : body rank  
";

%feature("docstring") JONCx_GetShape "

Get the shape of a JONCx.  

usage : shape = JONCx_GetShape(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of JONCx  

Returns
-------
shape (double array) : axis length of the JONCx  
";

%feature("docstring") JONCx_GetCoor "

Get the coor of a JONCx.  

usage : coor = JONCx_GetCoor(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of JONCx  

Returns
-------
coor (double array) : coordinates of the JONCx  
";

%feature("docstring") JONCx_GetPtrJONCx2BDYTY "

return a pointer onto the map joncx2rbdy2  

python usage : joncx2rbdy2 = JONCx_GetPtrJONCx2BDYTY()  

Returns
-------
joncx2rbdy2 (integer array) : reference on map between joncx rank and body/tact
rank  
";

%feature("docstring") JONCx_IsVisible "

return if a body visible  

usage : visible = JONCx_IsVisible(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of JONCx  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") JONCx_InitOutlines "

Get a reference on the outlines of all JONCx.  

usage : outlines = JONCx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_JONCx  
";

%feature("docstring") JONCx_InitScalarFields "

Get a reference on the scalar fields of all JONCx.  

usage : scalarfields = JONCx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_JONCx array  
";

%feature("docstring") JONCx_UpdatePostdata "

Update values of outlines_JONCx and scalarfields_JONCx pointers.  

usage : JONCx_UpdatePostdata()  
";

%feature("docstring") JONCx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = JONCx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
JONCx  
";

%feature("docstring") JONCx_GetNbScalarFields "

Get the number of scalar fields of a JONCx.  

python usage : nb_scalarfields = JONCx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a JONCx  
";

%feature("docstring") JONCx_CleanMemory "

Free all memory allocated within JONCx module.  

python usage : JONCx_CleanMemory()  
";

