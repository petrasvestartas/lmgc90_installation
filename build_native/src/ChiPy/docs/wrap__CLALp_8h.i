
// File: wrap__CLALp_8h.xml

%feature("docstring") CLALp_SelectProxTactors "

contact detection between CLxxx and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
If reset not equal to 0, the initialization flag is reset and detection skipped  

python usage : CLALp_SelectProxTactors(reset=0, use_external=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
* `use_external` :  
    (integer) : if not 0, external detection is used  
";

%feature("docstring") CLALp_UpdateWear "

python usage : CLALp_UpdateWear()  
";

%feature("docstring") CLALp_WriteLastVlocRloc "

write last local values of all CLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : CLALp_WriteLastVlocRloc()  
";

%feature("docstring") CLALp_WriteOutVlocRloc "

write local values of all CLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : CLALp_WriteOutVlocRloc()  
";

%feature("docstring") CLALp_DisplayOutVlocRloc "

display local values of all CLALp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CLALp_DisplayOutVlocRloc()  
";

%feature("docstring") CLALp_DisplayProxTactors "

display contacts  

python usage : CLALp_DisplayProxTactors()  
";

%feature("docstring") CLALp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

python usage : CLALp_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CLALp_SetNonSymmetricDetection "

this function allows non symmetric detection i.e. only one interaction is kept
when two bodies with candidate and antagonist contactors see each other  

python usage : CLALp_SetNonSymmetricDetection()  
";

%feature("docstring") CLALp_Trim "

trim contact (only contact within a line - not with extremities)  

python usage : CLALp_Trim()  
";

%feature("docstring") CLALp_CleanMemory "

Free all memory allocated within CLALp module.  

python usage : CLALp_CleanMemory()  
";

