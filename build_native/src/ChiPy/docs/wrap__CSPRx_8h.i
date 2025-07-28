
// File: wrap__CSPRx_8h.xml

%feature("docstring") CSPRx_SelectProxTactors "

contact detection between CSxxx and PRxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : CSPRx_SelectProxTactors(int reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CSPRx_WriteLastVlocRloc "

write last local values of all CSPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : CSPRx_WriteLastVlocRloc()  

The values written are relative velocity, forces and local frame  
";

%feature("docstring") CSPRx_WriteOutVlocRloc "

write local values of all CSPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : CSPRx_WriteOutVlocRloc()  
";

%feature("docstring") CSPRx_DisplayOutVlocRloc "

display local values of all CSPRx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CSPRx_DisplayOutVlocRloc()  
";

%feature("docstring") CSPRx_DisplayProxTactors "

display contacts  

python usage : CSPRx_DisplayProxTactors()  
";

%feature("docstring") CSPRx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CSPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CSPRx_Trim "

trim contact (only node face contact)  

python usage : CSPRx_Trim()  
";

%feature("docstring") CSPRx_GetInfo "

return contact info for the icdan CSPRx contact  

python usage : a = CSPRx_GetInfo(icdan)  

Parameters
----------
* `icdan` :  
    (integer) : contact identifiant  

Returns
-------
a (array integer) : info array  
";

%feature("docstring") CSPRx_Smoothing "

smooth contact reaction  

python usage : CSPRx_Smmothing()  
";

%feature("docstring") CSPRx_AddReac "

add contact force to body Reac  

python usage : CSPRx_AddReac()  
";

%feature("docstring") CSPRx_CleanMemory "

Free all memory allocated within CSPRx module.  

python usage : CSPRx_CleanMemory()  
";

