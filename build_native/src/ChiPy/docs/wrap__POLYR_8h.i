
// File: wrap__POLYR_8h.xml

%feature("docstring") POLYR_LoadTactors "

load POLYR from RBDY3 or MAILx and initialize existing_entites  

python usage : POLYR_LoadTactors()  
";

%feature("docstring") POLYR_GetContactorColor "

Get the color of a given POLYR.  

python usage : color = POLYR_GetContactorColor(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of POLYR  

Returns
-------
color (string) : the color of the POLYR itact  
";

%feature("docstring") POLYR_SaveVertex "

write position of vertex in a file  

python usage : POLYR_SaveVertex()  
";

%feature("docstring") POLYR_ModifyRadius "

apply an amplification/reduction size factor  

python usage : POLYR_ModifyRadius(ratio)  

Parameters
----------
* `ratio` :  
    (real) : ratio factor  
* `ratio` :  
    (double) : ratio factor  
";

%feature("docstring") POLYR_SetThresholdBigPolyr "

define the threshold between a plain and a big polyr. big polyr are such that
radius > threshold*mean_radius. default threshold = 4. Must be defined before
the load of tactors.  

python usage : POLYR_SetThresholdBigPolyr(ratio)  

Parameters
----------
* `ratio` :  
    (real) : ratio factor  
* `ratio` :  
    (double) : ratio factor  
";

%feature("docstring") POLYR_SetBigPolyr "

impose explicitly that an object is big. Must be set after the load of tactors.  

python usage : POLYR_SetBigPolyr(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of the polyr  
* `itacty` :  
    (int) : rank of the polyr  
";

%feature("docstring") POLYR_SetNbBigPolyr "

impose explicitly the number of big POLYR. Must be set after the load of
tactors.  

python usage : POLYR_SetNbBigPolyr(nb)  

Parameters
----------
* `nb` :  
    (integer) : number of polyr  
* `number` :  
    (int) : number of polyr  
";

%feature("docstring") POLYR_SkipTopoBigPolyr "

skip the topological decomposition of a big POLYR. its surface is considered as
a soup of triangle. usefull with complicated surface using Cundall CP detection  

python usage : POLYR_SkipTopoBigPolyr()  
";

%feature("docstring") POLYR_SkipAutomaticReorientation "

disable automatic reorientation (which works only with convex POLYR).  

python usage : POLYR_SkipAutomaticReorientation()  

  
 Disable the automatic reorientation of normals performed by lmgc90.  
 This is necessary when using non-convex objects.  
";

%feature("docstring") POLYR_SkipHEBuild "

disable Half-Edge structure generation (HE is necessary for non convex contact
detection)  

python usage : POLYR_SkipHEBuild()  

  
 Disable the Half-Edge structure generation performed by lmgc90.  
 This is necessary when testing the import of strange object.  
";

%feature("docstring") POLYR_TopologyAngle "

set the maximum angle (between 0 and 180 degree) threshold between 2 elements to
declare them as belonging to the same topological face  

python usage : POLYR_TopologyAngle(angle)  
";

%feature("docstring") POLYR_FlatnessAngle "

set the maximum angle (between 0 and 180 degree) variation between elements of a
topological face to declare it as flat  

python usage : POLYR_FlatnessAngle(angle)  
";

%feature("docstring") POLYR_GetWireframe "

Get wireframe of a POLYR.  

python usage : coor,connectivity = POLYR_GetWireframe(itacty, angle)  

Parameters
----------
* `itacty` :  
    (integer) : rank of the POLYR  
* `angle` :  
    (double) : threshold angle to skip some nodes on boundary of faces of the
    POLYR  

Returns
-------
coor (double array) : reference on the coor vector seen as a numpy array of size
[nb_point,3] connectivity (integer array) : reference on the connectivity vector
seen as a numpy array  
";

%feature("docstring") POLYR_GetVertex "

Get the outline of a POLYR in almost current configuration.  

If the POLYR is a real POLYR the current position of the center of the POLYR is
used but the local frame for the orientation is the on in detection
configuration.  

If the POLYR is in fact a POLYD the current position of nodes of the mesh are
used.  

usage : vertex = POLYR_GetVertex(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetPtrVertexTT "

Get a pointer on the outline of a POLYR in detection configuration.  

usage : vertex = POLYR_GetPtrVertexTT(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetPtrNormalTT "

Get a pointer on the outline of a POLYR in detection configuration - be carefull
to move polyr.  

usage : normal = POLYR_GetPtrNormalTT(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered POLYR  

Returns
-------
normal (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_MoveToConfigurationTT "

move the polyr in the configuration TT ; mandatory to get the wireframe in
deformed configuration  

python usage : POLYR_MoveToConfigurationTT()  
";

%feature("docstring") POLYR_GetPOLYR2BDYTY "

Get a copy of map POLYR2bdyty.  

usage : polyr2bdyty = POLYR_GetPOLYR2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") POLYR_GetPtrPOLYR2BDYTY "

Get a pointer on map POLYR2bdyty.  

usage : polyr2bdyty = POLYR_GetPtrPOLYR2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : a pointer in the polyr2bdyty map  
";

%feature("docstring") POLYR_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = POLYR_IsVisible(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") POLYR_GetNbPOLYR "

Get the number of POLYR.  

python usage : nb_POLYR = POLYR_GetNbPOLYR()  

Returns
-------
nb_POLYR (integer) : the number of POLYR  
";

%feature("docstring") POLYR_InitOutlines "

Get a reference on the outlines of all POLYR.  

usage : outlines = POLYR_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_POLYR  
";

%feature("docstring") POLYR_InitScalarFields "

Get a reference on the scalar fields of all POLYR.  

usage : scalarfields = POLYR_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_POLYR array  
";

%feature("docstring") POLYR_UpdatePostdata "

Update values of outlines_POLYR and scalarfields_POLYR pointers.  

usage : POLYR_UpdatePostdata()  
";

%feature("docstring") POLYR_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = POLYR_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
POLYR  
";

%feature("docstring") POLYR_GetNbScalarFields "

Get the number of scalar fields of a POLYR.  

python usage : nb_scalarfields = POLYR_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a POLYR  
";

%feature("docstring") POLYR_GetPtrAllConnectivities "

Get a reference on the connectivities of all POLYR.  

usage : connec = POLYR_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") POLYR_GetPtrConnectivity "

Get a reference on the connectivity of one POLYR.  

usage : connec = POLYR_GetPtrConnectivity(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : POLYR number  

Returns
-------
connec (integer 2D-array) : reference on connectivity  
";

%feature("docstring") POLYR_GetPtrVertexRef "

Get the position of the vertices of a POLYR in its inertia frame.  

usage : vertex = POLYR_GetPtrVertexRef(itacty)  

Parameters
----------
* `itacty` :  
    (integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetTopoData "

Get for each face of all POLYR : contactor id, topo id, face id and face status.  

usage : topo_data = POLYR_GetTopoData()  

Returns
-------
topt_data (int 2D-array) : topology data of all faces of all POLYR  
";

%feature("docstring") POLYR_CleanMemory "

Free all memory allocated within POLYR module.  

python usage : POLYR_CleanMemory()  
";

