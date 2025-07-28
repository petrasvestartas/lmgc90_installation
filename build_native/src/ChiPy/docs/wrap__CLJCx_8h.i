
// File: wrap__CLJCx_8h.xml

%feature("docstring") CLJCx_SelectProxTactors "

contact detection between CLxxx and JCxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : CLJCx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CLJCx_WriteLastVlocRloc "

write last local values of all CLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : CLJCx_WriteLastVlocRloc()  
";

%feature("docstring") CLJCx_WriteOutVlocRloc "

write local values of all CLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : CLJCx_WriteOutVlocRloc()  
";

%feature("docstring") CLJCx_DisplayOutVlocRloc "

display local values of all CLJCx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CLJCx_DisplayOutVlocRloc()  
";

%feature("docstring") CLJCx_DisplayProxTactors "

display contacts  

python usage : CLJCx_DisplayProxTactors()  
";

%feature("docstring") CLJCx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CLJCx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CLJCx_CleanMemory "

Free all memory allocated within CLJCx module.  

python usage : CLJCx_CleanMemory()  
";

