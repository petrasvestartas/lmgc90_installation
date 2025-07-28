
// File: wrap__SPPRx_8h.xml

%feature("docstring") SPPRx_SelectProxTactors "

contact detection between SPxxx and PLxxx tactors  

python usage : SPPRx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPPRx_SmoothForceComputation "

compute smooth contact law (in any)  

python usage : SPPRx_SmoothForceComputation()  
";

%feature("docstring") SPPRx_WriteLastVlocRloc "

write last local values of all SPPRx contacts  

python usage : SPPRx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_WriteOutVlocRloc "

write local values of all SPPRx contacts  

python usage : SPPRx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_DisplayOutVlocRloc "

display local values of all SPPRx contacts  

python usage : SPPRx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_DisplayProxTactors "

display contacts  

python usage : SPPRx_DisplayProxTactors()  
";

%feature("docstring") SPPRx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") SPPRx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPPRx_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (real) : period on x axis  
";

%feature("docstring") SPPRx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPPRx_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (real) : period on y axis  
";

%feature("docstring") SPPRx_CleanMemory "

Free all memory allocated within SPPRx module.  

python usage : SPPRx_CleanMemory()  
";

