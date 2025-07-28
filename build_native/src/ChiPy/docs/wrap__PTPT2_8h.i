
// File: wrap__PTPT2_8h.xml

%feature("docstring") PTPT2_SelectProxTactors "

contact detection between PT2Dx tactors  

python usage : PTPT2_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PTPT2_WriteLastVlocRloc "

write last local values of all PTPT2 contacts  

python usage : PTPT2_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_WriteOutVlocRloc "

write local values of all PTPT2 contacts  

python usage : PTPT2_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_DisplayOutVlocRloc "

display local values of all PTPT2 contacts  

python usage : PTPT2_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_DisplayProxTactors "

display contacts  

python usage : PTPT2_DisplayProxTactors()  
";

%feature("docstring") PTPT2_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PTPT2_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PTPT2_LoadNetwork "

read a PTPT2 network from a file  

python usage : PTPT2_LoadNetwork()  
";

%feature("docstring") PTPT2_SetTolerance "

set the maximum violation for a point to point link  

python usage : PTPT2_SetTolerance(tol)  
";

%feature("docstring") PTPT2_SetExplicitLocalFrame "

local frame is computed only once at the first step  

python usage : PTPT2_SetExplicitLocalFrame()  
";

%feature("docstring") PTPT2_LoadParams "

read a PTPT2 surface and l0 from a file  

python usage : PTPT2_LoadParams()  
";

%feature("docstring") PTPT2_UseCurrentNonuc0 "

Use GetCoor or value given from file insted of computing nonuc0 from reference
coordinates.  

python usage : PTPT2_UseCurrentNonuc0(to_use) param[in] to_use (integer) : 1 to
activate, 0 to deactivate feature  
";

%feature("docstring") PTPT2_CleanMemory "

Free all memory allocated within PTPT2 module.  

python usage : PTPT2_CleanMemory()  
";

