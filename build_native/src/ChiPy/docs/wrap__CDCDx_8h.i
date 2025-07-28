
// File: wrap__CDCDx_8h.xml

%feature("docstring") CDCDx_SelectProxTactors "

contact detection between CYLND and CYLND tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list  

python usage : CDCDx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CDCDx_SmoothForceComputation "

computes smooth contact forces (if any)  

python usage : CDCDx_SmoothForceComputation()  
";

%feature("docstring") CDCDx_WriteLastVlocRloc "

write last local values of all CDCDx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDCDx_WriteLastVlocRloc()  
";

%feature("docstring") CDCDx_WriteOutVlocRloc "

write local values of all CDCDx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDCDx_WriteOutVlocRloc()  
";

%feature("docstring") CDCDx_DisplayOutVlocRloc "

display local values of all CDCDx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CDCDx_DisplayOutVlocRloc()  
";

%feature("docstring") CDCDx_DisplayProxTactors "

display detected contacts  

python usage : CDCDx_DisplayProxTactors()  
";

%feature("docstring") CDCDx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read,
    -   num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : CDCDx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CDCDx_SetXPeriodicCondition "

initialise data for simulation using periodic condition along X  

python usage : CDCDx_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (double) : period on x axis  
";

%feature("docstring") CDCDx_SetYPeriodicCondition "

initialise data for simulation using periodic condition along Y  

python usage : CDCDx_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (double) : period on y axis  
";

%feature("docstring") CDCDx_SetNumberInterByContact "

define the number of interaction by contact (experimental)  

python usage : CDCDx_SetNumberInterByContact(nb_interactions)  

Parameters
----------
* `nb_interactions` :  
    (integer) : number of interactions per contact  
";

%feature("docstring") CDCDx_SetContactRadius "

define the contact radius (experimental)  

python usage : CDCDx_SetContactRadius(radius)  

Parameters
----------
* `radius` :  
    (double) : contact radius  
";

%feature("docstring") CDCDx_CleanMemory "

Free all memory allocated within CDCDx module.  

python usage : CDCDx_CleanMemory()  
";

