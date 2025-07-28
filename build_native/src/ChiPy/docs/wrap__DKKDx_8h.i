
// File: wrap__DKKDx_8h.xml

%feature("docstring") DKKDx_SelectProxTactors "

contact detection between DISKx and xKSID tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKKDx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKKDx_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKKDx_SmoothForceComputation()  
";

%feature("docstring") DKKDx_WriteLastVlocRloc "

write last local values of all DKKDx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKKDx_WriteLastVlocRloc()  
";

%feature("docstring") DKKDx_WriteOutVlocRloc "

write local values of all DKKDx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKKDx_WriteOutVlocRloc()  
";

%feature("docstring") DKKDx_DisplayOutVlocRloc "

display local values of all DKKDx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKKDx_DisplayOutVlocRloc()  
";

%feature("docstring") DKKDx_DisplayProxTactors "

display contacts  

python usage : DKKDx_DisplayProxTactors()  
";

%feature("docstring") DKKDx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKKDx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKKDx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKKDx_SetSurfaceSectors(nbsect)  

Parameters
----------
* `nbsect` :  
    (integer) : number of sectors  
";

%feature("docstring") DKKDx_CleanMemory "

Free all memory allocated within DKKDx module.  

python usage : DKKDx_CleanMemory()  
";

