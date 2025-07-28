
// File: wrap__PLJCx_8h.xml

%feature("docstring") PLJCx_SelectProxTactors "

contact detection between POLYG and JONCx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLJCx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") PLJCx_WriteLastVlocRloc "

write last local values of all PLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : PLJCx_WriteLastVlocRloc()  
";

%feature("docstring") PLJCx_WriteOutVlocRloc "

write local values of all PLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : PLJCx_WriteOutVlocRloc()  
";

%feature("docstring") PLJCx_DisplayOutVlocRloc "

display local values of all PLJCx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : PLJCx_DisplayOutVlocRloc()  
";

%feature("docstring") PLJCx_DisplayProxTactors "

display contacts  

python usage : PLJCx_DisplayProxTactors()  
";

%feature("docstring") PLJCx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : PLJCx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PLJCx_ComputeStress "

compute the PLJC contribution to the equivalent stress tensor  

python usage : PLJCx_ComputeStress()  
";

%feature("docstring") PLJCx_CleanMemory "

Free all memory allocated within PLJCx module.  

python usage : PLJCx_CleanMemory()  
";

%feature("docstring") PLJCx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : PLJCx_SetFrictionModel(cflag)  

Parameters
----------
* `cflag` :  
    (char) : model to use ('min', 'max' or 'ave')  
";

