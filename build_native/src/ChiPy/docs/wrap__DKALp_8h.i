
// File: wrap__DKALp_8h.xml

%feature("docstring") DKALp_SelectProxTactors "

contact detection between DISKx and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKALp_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKALp_WriteLastVlocRloc "

write last local values of all DKALp contacts  

The values written are relative velocity, forces and local frame  

python usage : DKALp_WriteLastVlocRloc()  
";

%feature("docstring") DKALp_WriteOutVlocRloc "

write local values of all DKALp contacts  

The values written are relative velocity, forces and local frame  

python usage : DKALp_WriteOutVlocRloc()  
";

%feature("docstring") DKALp_DisplayOutVlocRloc "

display local values of all DKALp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKALp_DisplayOutVlocRloc()  
";

%feature("docstring") DKALp_DisplayProxTactors "

display contacts  

python usage : DKALp_DisplayProxTactors()  
";

%feature("docstring") DKALp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

python usage : DKALp_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKALp_CleanMemory "

Free all memory allocated within DKALp module.  

python usage : DKALp_CleanMemory()  
";

