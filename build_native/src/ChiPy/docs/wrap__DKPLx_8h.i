
// File: wrap__DKPLx_8h.xml

%feature("docstring") DKPLx_SelectProxTactors "

contact detection between DISKx and POLYG tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKPLx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKPLx_WriteLastVlocRloc "

write last local values of all DKPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKPLx_WriteLastVlocRloc()  
";

%feature("docstring") DKPLx_WriteOutVlocRloc "

write local values of all DKPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKPLx_WriteOutVlocRloc()  
";

%feature("docstring") DKPLx_DisplayOutVlocRloc "

display local values of all DKPLx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKPLx_DisplayOutVlocRloc()  
";

%feature("docstring") DKPLx_DisplayProxTactors "

display contacts  

python usage : DKPLx_DisplayProxTactors()  
";

%feature("docstring") DKPLx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKPLx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKPLx_SetPeriodicCondition(period)  

Parameters
----------
* `period` :  
    (double) : value of the period  
";

%feature("docstring") DKPLx_CleanMemory "

Free all memory allocated within DKPLx module.  

python usage : DKPLx_CleanMemory()  
";

