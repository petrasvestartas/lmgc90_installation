
// File: wrap__MAILx_8h.xml

%feature("docstring") MAILx_ReadBodies "

read MAILx from DATBOX/BODIES.DAT  

Input string is of form vX.Y where X is major version number and Y is minor one.  
 If not specified, last available version is used.  

python usage : MAILx_ReadBodies(version)  

param[in] version (string) : file format version to use  
";

%feature("docstring") MAILx_WriteBodies "

write MAILx to OUTBOX/BODIES.OUT  

Input string is of form vX.Y where X is major version number and Y is minor one.  
 If not specified, last available version is used.  

python usage : MAILx_WriteBodies(version)  

param[in] version (string) : file format version to use  
";

%feature("docstring") MAILx_WriteLastGPV "

write OUTBOX/GPV.LAST  

python usage : MAILx_WriteLastGPV()  
";

%feature("docstring") MAILx_WriteOutGPV "

write OUTBOX/GPV.OUT.x  

python usage : MAILx_WriteOutGPV()  
";

%feature("docstring") MAILx_DisplayOutGPV "

display GPV values  

python usage : MAILx_DisplayOutGPV()  
";

%feature("docstring") MAILx_AddDof2InBodies "

set cooref = cooref + X  

python usage : MAILx_AddDof2InBodies()  
";

%feature("docstring") MAILx_GetNbMAILx "

Get the number of MAILx.  

python usage : nb_MAILx = GetNbMAILx()  

Returns
-------
nb_MAILx (integer) : number of MAILx  
";

%feature("docstring") MAILx_GetNbCell "

Get the number of Cells of a given MAILx.  

python usage : nb_MAILx = GetNbCell(IdBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  

Returns
-------
nb_cell (integer) : number of cell  
";

%feature("docstring") MAILx_SetCoorRef "

set reference coordinates on a given body  

python usage : MAILx_SetCoorRef(IdBody, f, length)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f` :  
    (double array) : value of the vitesse  
* `length` :  
    (integer) : length of vector  
";

%feature("docstring") MAILx_GetCoordNodty "

Get one coordinate of a node of a body.  

python usage : x = MAILx_GetCoordNodty(int ibdty,int inodty,int icomp)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `inodty` :  
    (integer) : the node  
* `icomp` :  
    (integer) : the component  

Returns
-------
x (double) : coordinate of node  
";

%feature("docstring") MAILx_GetCoordsNodty "

Get the coordinates of a node of a body.  

python usage : x = MAILx_GetCoordsNodty(int ibdty, int inodty, int length)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `inodty` :  
    (integer) : the node  
* `length` :  
    (integer) : the number of component  

Returns
-------
x (double array) : the desired vector  
";

%feature("docstring") MAILx_GetNbNodes "

Get the number of nodes of a given MAILx.  

python usage : nb_nodes = MAILx_GetNbNodes(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : body id  

Returns
-------
nb_nodes (integer) : number of nodes of the body  
";

%feature("docstring") MAILx_InitNodalFields "

Set the number of nodal_fields for a given body.  

python usage : MAILx_InitNodalFields(ibdyty,nb_nodal_fields)  

Parameters
----------
* `ibdyty` :  
    (integer) : body id  
* `nb_nodal_fields` :  
    (integer) : number of nodal fields required  
";

%feature("docstring") MAILx_InitNodalField "

Set name and size of a nodal_field of a given body.  

python usage : MAILx_InitNodalField(ibdyty,name,rank,sz)  

Parameters
----------
* `ibdyty` :  
    (integer) : body id  
* `name` :  
    (string) : field name  
* `rank` :  
    (integer) : field rank  
* `sz` :  
    (integer) : size of the field  
";

%feature("docstring") MAILx_SetNodalField "

Set a nodal_field of a given body.  

python usage : MAILx_SetNodalField(ibdyty,rank,field)  

Parameters
----------
* `ibdyty` :  
    (integer) : body id  
* `rank` :  
    (integer) : field rank  
* `field` :  
    (double vector) : field  
";

%feature("docstring") MAILx_CleanMemory "

Free all memory allocated within MAILx module.  

python usage : MAILx_CleanMemory()  
";

