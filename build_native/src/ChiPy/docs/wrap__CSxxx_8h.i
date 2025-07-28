
// File: wrap__CSxxx_8h.xml

%feature("docstring") CSxxx_LoadTactors "

Load CSxxx from MAILx and Initialize existing_entities.  

python usage : CSxxx_LoadTactors()  
";

%feature("docstring") CSxxx_PushPreconNodes "

set CSxxx supporting nodes as precon  

python usage : CSxxx_PushPreconNodes()  
";

%feature("docstring") CSxxx_FlipOrientation "

Flip normal of all CSxxx of a given MAILx body.  

python usage : CSxxx_FlipOrientation(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of desired body  
";

%feature("docstring") CSxxx_FlipOrientationOnePatch "

Flip normal of CSxxx belonging to given patch of a given MAILx body.  

python usage : CSxxx_FlipOrientationOnePatch(ibdyty,icspxx)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of desired body  
* `icspxx` :  
    (integer) : rank of desired patch  
";

%feature("docstring") CSxxx_SetShrink "

shrink position of nodes in CSxxx contactors  

python usage : CSxxx_SetShrink(shrink)  

Parameters
----------
* `shrink` :  
    (real) : shrink value  
";

%feature("docstring") CSxxx_SetQuadrature "

Set the contact quadrature rule of a CSxxx face. OBSOLETE FUNCTION !!!! To
remove in the future.  

python usage : CSxxx_SetQuadrature(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : degree on CSxxx contactor  
";

%feature("docstring") CSxxx_AddReac "

Apply an external reaction on a CSxxx.  

python usage : CSxxx_AddReac(datatype, iCSxxx, reac)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `iCSxxx` :  
    (integer) : id of the CSpxx  
* `reac` :  
    (double array) : the value to add  
";

%feature("docstring") CSpxx_ApplySurfaceLoad "
";

%feature("docstring") CSpxx_ApplyPressure "

Apply an external pressure on a CSpxx.  

python usage : CSpxx_ApplyPressure(ivalue,rvalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the CSpxx  
* `rvalue` :  
    (real) : pressure  
";

%feature("docstring") CSxxx_GetNbCSxxx "

Get the number of CSxxx.  

usage : nb_CSxxx = CSxxx_GetNbCSxxx()  

Parameters
----------
* `nb_CSxxx` :  
    (integer) : number of CSxxx in container  
";

%feature("docstring") CSpxx_GetAllConnec "

return connectivity of all CS in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = CSxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of CSxxx elements  
";

%feature("docstring") CSpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all CSxxx  

python usage : idata, rdata = CSxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") CSxxx_CleanMemory "

Free all memory allocated within CSxxx module.  

python usage : CSxxx_CleanMemory()  
";

