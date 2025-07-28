
// File: wrap__SPDCx_8h.xml

%feature("docstring") SPDCx_SelectProxTactors "

contact detection between SPxxx and DCxxx tactors  

python usage : SPDCx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPDCx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : SPDCx_SmoothForceComputation()  
";

%feature("docstring") SPDCx_WriteLastVlocRloc "

write last local values of all SPDCx contacts  

python usage : SPDCx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_WriteOutVlocRloc "

write local values of all SPDCx contacts  

python usage : SPDCx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_DisplayOutVlocRloc "

display local values of all SPDCx contacts  

python usage : SPDCx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_DisplayProxTactors "

display contacts  

python usage : SPDCx_DisplayProxTactors()  
";

%feature("docstring") SPDCx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPDCx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") SPDCx_CleanMemory "

Free all memory allocated within SPDCx module.  

python usage : SPDCx_CleanMemory()  
";

