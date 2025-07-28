
// File: wrap__POLYG_8h.xml

%feature("docstring") POLYG_LoadTactors "

load POLYG from RBDY2 and initialize existing_entites  

python usage : POLYG_LoadTactors()  
";

%feature("docstring") POLYG_GetMinRadius "

give min radius used during detection  

python usage : POLYG_GetMinRadius()  
";

%feature("docstring") POLYG_GetMaxRadius "

give max radius used during detection  

python usage : POLYG_GetMaxRadius()  
";

%feature("docstring") POLYG_GetNbPOLYG "

Get the number of POLYG in container.  

python usage : nb_polyg = POLYG_GetNbPOLYG()  

Returns
-------
nb_polyg (integer) : the number of POLYG in container  
";

%feature("docstring") POLYG_GetPOLYG2BDYTY "

Get a copy of map POLYG2bdyty.  

usage : polyr2bdyty = POLYG_GetPOLYG2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") POLYG_GetPtrPOLYG2BDYTY "

return a pointer onto the map polyg2rbdy2  

python usage : polyg2rbdy2 = POLYG_GetPtrPOLYG2BDYTY()  

Returns
-------
polyg2rbdy2 (integer array) : reference on map between polyg rank and
body/tactor rank  
";

%feature("docstring") POLYG_IsVisible "

return if a body visible  

usage : visible = POLYG_IsVisible(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of POLYG  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") POLYG_GetContactorRadius "

Get the radius of a given POLYG.  

python usage : radius = POLYG_GetContactorRadius(itact)  

Parameters
----------
* `itact` :  
    (integer) : rank of a POLYG  

Returns
-------
radius (double) : the radius of the POLYG of rank itact  
";

%feature("docstring") POLYG_GetNbVertices "

Get the number of vertices of the first POLYG of a body.  

python usage : nb_vertices = POLYG_GetNbVertices(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of a body  

Returns
-------
nb_vertices (integer) : the number of vertices of the first POLYG of the body  
";

%feature("docstring") POLYG_GetVertices "

Get the coordinates of the vertices of the first POLYG of a body.  

usage : vertices = POLYG_GetVertices(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `vertices` :  
    (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYG_GetNbVertex "

Get the number of vertices of a POLYG.  

usage : nb_vertex = POLYG_GetNpVertex(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the POLYG contactor  

Returns
-------
nb_vertex (int) : the number of vertices of the POLYG  
";

%feature("docstring") POLYG_GetVertex "

Get the outline of a POLYG.  

usage : vertex = POLYG_GetVertex(itacty, length)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `length` :  
    (integer) : 2 * number of vertices  
* `vertex` :  
    (double array) : the coordinates of the vertices  
";

%feature("docstring") POLYG_GetBodyId "

Get the id of the body which the tactor belongs.  

python usage : id = POLYG_GetBodyId(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of a POLYG contactor  

Returns
-------
id (integer) : the id of the body  
";

%feature("docstring") POLYG_InitOutlines "

Get a reference on the outlines of all POLYG.  

usage : outlines = POLYG_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_POLYG  
";

%feature("docstring") POLYG_InitScalarFields "

Get a reference on the scalar fields of all POLYG.  

usage : scalarfields = POLYG_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_POLYG array  
";

%feature("docstring") POLYG_UpdatePostdata "

Update values of outlines_POLYG and scalarfields_POLYG pointers.  

usage : POLYG_UpdatePostdata  
";

%feature("docstring") POLYG_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = POLYG_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
POLYG  
";

%feature("docstring") POLYG_GetNbScalarFields "

Get the number of scalar fields of a POLYG.  

python usage : nb_scalarfields = POLYG_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a POLYG  
";

%feature("docstring") POLYG_SetXdilation "
";

%feature("docstring") POLYG_SetVdilation "
";

%feature("docstring") POLYG_CleanMemory "

Free all memory allocated within POLYG module.  

python usage : POLYG_CleanMemory()  
";

