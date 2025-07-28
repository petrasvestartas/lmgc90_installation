
// File: wrap__DKDKL_8h.xml

%feature("docstring") DKDKL_SelectProxTactors "

contact detection between DISKx and DISKL tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKDKL_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKDKL_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKDKL_SmoothForceComputation  
";

%feature("docstring") DKDKL_WriteLastVlocRloc "

write last local values of all DKDKL contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKL_WriteLastVlocRloc()  

The values written are relative velocity, forces and local frame  
";

%feature("docstring") DKDKL_WriteOutVlocRloc "

write local values of all DKDKL contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKL_WriteOutVlocRloc()  
";

%feature("docstring") DKDKL_DisplayOutVlocRloc "

display local values of all DKDKL contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKDKL_DisplayOutVlocRloc()  
";

%feature("docstring") DKDKL_DisplayProxTactors "

display contacts  

python usage : DKDKL_DisplayProxTactors()  
";

%feature("docstring") DKDKL_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKDKL_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKDKL_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKDKL_SetPeriodicCondition(period)  

Parameters
----------
* `period` :  
    (double) : value of the period  
";

%feature("docstring") DKDKL_CleanMemory "

Free all memory allocated within DKDKL module.  

python usage : DKDKL_CleanMemory()  
";

