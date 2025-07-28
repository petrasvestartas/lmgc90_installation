
// File: wrap__PRASp_8h.xml

%feature("docstring") PRASp_SelectProxTactors "

contact detection between PRxxx and ASpxx tactors  

python usage : PRASp_SelectProxTactors(int reset=0) param[in] reset (integer) :
if not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PRASp_WriteLastVlocRloc "

write last local values of all PRASp contacts  

python usage : PRASp_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRASp_WriteOutVlocRloc "

write local values of all PRASp contacts  

python usage : PRASp_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRASp_DisplayOutVlocRloc "

display local values of all PRASp contacts  

python usage : PRASp_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PRASp_DisplayProxTactors "

display contacts  

python usage : PRASp_DisplayProxTactors()  
";

%feature("docstring") PRASp_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PRASp_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PRASp_CleanMemory "

Free all memory allocated within PRASp module.  

python usage : PRASp_CleanMemory()  
";

