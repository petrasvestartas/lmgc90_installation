
// File: wrap__CDPLx_8h.xml

%feature("docstring") CDPLx_SelectProxTactors "

contact detection between CYLND and PLANx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list  

python usage : CDPLx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CDPLx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : CDPLx_SmoothForceComputation()  
";

%feature("docstring") CDPLx_WriteLastVlocRloc "

write last local values of all CDPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDPLx_WriteLastVlocRloc()  
";

%feature("docstring") CDPLx_WriteOutVlocRloc "

write local values of all CDPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDPLx_WriteOutVlocRloc()  
";

%feature("docstring") CDPLx_DisplayOutVlocRloc "

display local values of all CDPLx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CDPLx_DisplayOutVlocRloc()  
";

%feature("docstring") CDPLx_DisplayProxTactors "

display contacts  

python usage : CDPLx_DisplayProxTactors()  
";

%feature("docstring") CDPLx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : CDPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CDPLx_CleanMemory "

Free all memory allocated within CDPLx module.  

python usage : CDPLx_CleanMemory()  
";

