
// File: wrap__ASpxx_8h.xml

%feature("docstring") ASpxx_LoadTactors "

Load ASpxx from MAILx and Initialize existing_entities.  

python usage : ASpxx_LoadTactors()  
";

%feature("docstring") ASpxx_PushPreconNodes "

set ASpxx supporting nodes as precon  

python usage : ASpxx_PushPreconNodes()  
";

%feature("docstring") ASpxx_GetAllConnec "

return connectivity of all AS in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = ASxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of ASxxx elements  
";

%feature("docstring") ASpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all ASxxx  

python usage : idata, rdata = ASxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") ASpxx_CleanMemory "

Free all memory allocated within ASpxx module.  

python usage : ASpxx_CleanMemory()  
";

%feature("docstring") ASpxx_ExplodePatch "

Explode ASpxx patch in singleton.  

python usage : ASpxx_ExplodePatch()  
";

