
// File: wrap__SPCDx_8h.xml

%feature("docstring") SPCDx_SelectProxTactors "

contact detection between SPxxx and CDxxx tactors  

python usage : SPCDx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPCDx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : SPCDx_SmoothForceComputation()  
";

%feature("docstring") SPCDx_WriteLastVlocRloc "

write last local values of all SPCDx contacts  

python usage : SPCDx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_WriteOutVlocRloc "

write local values of all SPCDx contacts  

python usage : SPCDx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_DisplayOutVlocRloc "

display local values of all SPCDx contacts  

python usage : SPCDx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_DisplayProxTactors "

display contacts  

python usage : SPCDx_DisplayProxTactors()  
";

%feature("docstring") SPCDx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPCDx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") SPCDx_CleanMemory "

Free all memory allocated within SPCDx module.  

python usage : SPCDx_CleanMemory()  
";

