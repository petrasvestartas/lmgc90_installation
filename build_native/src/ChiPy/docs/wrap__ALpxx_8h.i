
// File: wrap__ALpxx_8h.xml

%feature("docstring") ALpxx_LoadTactors "

load ALpxx from MAILx and initialize existing_entities  

python usage : ALpxx_LoadTactors()  
";

%feature("docstring") ALpxx_PushPreconNodes "

set ALpxx supporting nodes as precon  

python usage : ALpxx_PushPreconNodes()  
";

%feature("docstring") ALpxx_GetAllConnec "

return connectivity of all AL in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = ALxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of ALxxx elements  
";

%feature("docstring") ALpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all ALxxx  

python usage : idata, rdata = ALxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") ALpxx_CleanMemory "

Free all memory allocated within ALpxx module.  

python usage : ALpxx_CleanMemory()  
";

