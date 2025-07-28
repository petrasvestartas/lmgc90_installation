
// File: wrap__PTPT3_8h.xml

%feature("docstring") PTPT3_SelectProxTactors "

contact detection between PTxxx and PTxxx tactors  

python usage : PTPT3_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PTPT3_SmoothForceComputation "

computes smooth forces (if any)  

python usage : PTPT3_SmoothForceComputation()  
";

%feature("docstring") PTPT3_WriteLastVlocRloc "

write last local values of all PTPT3 contacts  

python usage : PTPT3_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_WriteOutVlocRloc "

write local values of all PTPT3 contacts  

python usage : PTPT3_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_DisplayOutVlocRloc "

display local values of all PTPT3 contacts  

python usage : PTPT3_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_DisplayProxTactors "

display contacts  

python usage : PTPT3_DisplayProxTactors()  
";

%feature("docstring") PTPT3_LoadNetwork "

read a PTPT3 network from a file  

python usage : PTPT3_LoadNetwork()  
";

%feature("docstring") PTPT3_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PTPT3_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PTPT3_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PTPT3_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (real) : period on x axis  
";

%feature("docstring") PTPT3_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PTPT3_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (real) : period on y axis  
";

%feature("docstring") PTPT3_SetExplicitLocalFrame "

local frame is computed only once at the first step  

python usage : PTPT3_SetExplicitLocalFrame()  
";

%feature("docstring") PTPT3_LoadParams "

read a PTPT3 surface and l0 from a file  

python usage : PTPT3_LoadParams()  
";

%feature("docstring") PTPT3_UseCurrentNonuc0 "

Use GetCoor or value given from file insted of computing nonuc0 from reference
coordinates.  

python usage : PTPT3_UseCurrentNonuc0(to_use) param[in] to_use (integer) : 1 to
activate, 0 to deactivate feature  
";

%feature("docstring") PTPT3_CleanMemory "

Free all memory allocated within PTPT3 module.  

python usage : PTPT3_CleanMemory()  
";

