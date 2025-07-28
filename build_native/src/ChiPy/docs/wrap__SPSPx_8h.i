
// File: wrap__SPSPx_8h.xml

%feature("docstring") SPSPx_SelectProxTactors "

contact detection between SPxxx and SPxxx tactors  

python usage : SPSPx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPSPx_SmoothForceComputation "

recup values of local contact forces of the last time step  

python usage : SPSPx_SmoothForceComputation()  
";

%feature("docstring") SPSPx_WriteLastVlocRloc "

write last local values of all SPSPx contacts  

python usage : SPSPx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_WriteOutVlocRloc "

write local values of all SPSPx contacts  

python usage : SPSPx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_DisplayOutVlocRloc "

display local values of all SPSPx contacts  

python usage : SPSPx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_DisplayProxTactors "

display contacts  

python usage : SPSPx_DisplayProxTactors()  
";

%feature("docstring") SPSPx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPSPx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") SPSPx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPSPx_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (real) : period on x axis  
";

%feature("docstring") SPSPx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPSPx_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (real) : period on y axis  
";

%feature("docstring") SPSPx_SetNumberInterByContact "

define the number of interaction by contact (experimental)  

python usage : SPSPx_SetNumberInterByContact(nb_interactions)  

Parameters
----------
* `nb_interactions` :  
    (integer) : number of interactions per contact  
";

%feature("docstring") SPSPx_SetContactRadius "

define the contact radius (experimental)  

python usage : SPSPx_SetContactRadius(radius)  

Parameters
----------
* `radius` :  
    (real) : contact radius  
";

%feature("docstring") SPSPx_FdSelectProxTactors "

contact detection between SPHER and SPHER tactors  

python usage : SPSPx_FdSelectProxTactors()  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPSPx_CleanMemory "

Free all memory allocated within SPSPx module.  

python usage : SPSPx_CleanMemory()  
";

