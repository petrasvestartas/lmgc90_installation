
// File: wrap__PLALp_8h.xml

%feature("docstring") PLALp_SelectProxTactors "

contact detection between POLYG and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLALp_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") PLALp_WriteLastVlocRloc "

write last local values of all PLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : PLALp_WriteLastVlocRloc()  
";

%feature("docstring") PLALp_WriteOutVlocRloc "

write local values of all PLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : PLALp_WriteOutVlocRloc()  
";

%feature("docstring") PLALp_DisplayOutVlocRloc "

display local values of all PLALp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : PLALp_DisplayOutVlocRloc()  
";

%feature("docstring") PLALp_DisplayProxTactors "

display contacts  

python usage : PLALp_DisplayProxTactors()  
";

%feature("docstring") PLALp_ReadIniVlocRloc "

Read VlocRloc file.  

-If num <= 0 : DATBOX/VlocRloc.INI file is read -Else : OUTBOX/VlocRloc.OUT.num
is read, num being  

*   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : PLALp_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PLALp_CleanMemory "

Free all memory allocated within PLALp module.  

python usage : PLALp_CleanMemory()  
";

