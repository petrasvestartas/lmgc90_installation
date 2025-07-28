
// File: wrap__CLxxx_8h.xml

%feature("docstring") CLxxx_LoadTactors "

load CLxxx from MAILx and Initialize existing_entities  

python usage : CLxxx_LoadTactors()  
";

%feature("docstring") CLxxx_SetNbNodesByCLxxx "

Set the number of CL nodes by edges. It helps to compute the length associated
to a contact node. Default is 2.  

python usage : CLxxx_SetNbNodesByCLxxx(nb_nodes)  

Parameters
----------
* `nb_nodes` :  
    (integer) : number of CLxxx contactors by edges  
";

%feature("docstring") CLxxx_PushPreconNodes "

set CLxxx supporting nodes as precon  

python usage : CLxxx_PushPreconNodes()  
";

%feature("docstring") CLxxx_GetNbCLxxx "

Get the number of CLxxx.  

usage : nb_CLxxx = CLxxx_GetNbCLxxx()  

Parameters
----------
* `nb_CLxxx` :  
    (integer) : number of CLxxx in container  
";

%feature("docstring") CLpxx_GetAllConnec "

return connectivity of all CL in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = CLxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of CLxxx elements  
";

%feature("docstring") CLpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all CLxxx  

python usage : idata, rdata = CLxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") CLxxx_CleanMemory "

Free all memory allocated within CLxxx module.  

python usage : CLxxx_CleanMemory()  
";

