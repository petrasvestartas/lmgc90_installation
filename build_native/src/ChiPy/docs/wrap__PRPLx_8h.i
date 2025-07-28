
// File: wrap__PRPLx_8h.xml

%feature("docstring") PRPLx_SelectProxTactors "

contact detection between PRxxx and PLxxx tactors  

python usage : PRPLx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PRPLx_WriteLastVlocRloc "

write last local values of all PRPLx contacts  

python usage : PRPLx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_WriteOutVlocRloc "

write local values of all PRPLx contacts  

python usage : PRPLx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_DisplayOutVlocRloc "

display local values of all PRPLx contacts  

python usage : PRPLx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_DisplayProxTactors "

display contacts  

python usage : PRPLx_DisplayProxTactors()  
";

%feature("docstring") PRPLx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PRPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PRPLx_CleanMemory "

Free all memory allocated within PRPLx module.  

python usage : PRPLx_CleanMemory()  
";

