
// File: wrap__DKJCx_8h.xml

%feature("docstring") DKJCx_SelectProxTactors "

contact detection between DISKx and JONCx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKJCx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKJCx_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKJCx_SmoothForceComputation()  
";

%feature("docstring") DKJCx_WriteLastVlocRloc "

write last local values of all DKJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKJCx_WriteLastVlocRloc()  
";

%feature("docstring") DKJCx_WriteOutVlocRloc "

write local values of all DKJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKJCx_WriteOutVlocRloc()  
";

%feature("docstring") DKJCx_DisplayOutVlocRloc "

display local values of all DKJCx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKJCx_DisplayOutVlocRloc()  
";

%feature("docstring") DKJCx_DisplayProxTactors "

display contacts  

python usage : DKJCx_DisplayProxTactors()  
";

%feature("docstring") DKJCx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKJCx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKJCx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKJCx_SetSurfaceSectors(nbsect)  

Parameters
----------
* `nbsect` :  
    (integer) : number of sectors  
";

%feature("docstring") DKJCx_ComputeStress "

compute the DKJC contribution to the equivalent stress tensor  

python usage : DKJCx_ComputeStress()  
";

%feature("docstring") DKJCx_CleanMemory "

Free all memory allocated within DKJCx module.  

python usage : DKJCx_CleanMemory()  
";

%feature("docstring") DKJCx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : DKJCx_SetFrictionModel(cflag)  

Parameters
----------
* `cflag` :  
    (char) : model to use ('min', 'max' or 'ave')  
";

