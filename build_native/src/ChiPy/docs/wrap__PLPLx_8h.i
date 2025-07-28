
// File: wrap__PLPLx_8h.xml

%feature("docstring") PLPLx_SelectProxTactors "

contact detection between POLYG tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLPLx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") PLPLx_WriteLastVlocRloc "

write last local values of all PLPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : PLPLx_WriteLastVlocRloc()  
";

%feature("docstring") PLPLx_WriteOutVlocRloc "

write local values of all PLPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : PLPLx_WriteOutVlocRloc()  
";

%feature("docstring") PLPLx_DisplayOutVlocRloc "

display local values of all PLPLx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : PLPLx_DisplayOutVlocRloc()  
";

%feature("docstring") PLPLx_DisplayProxTactors "

display contacts  

python usage : PLPLx_DisplayProxTactors()  
";

%feature("docstring") PLPLx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : PLPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PLPLx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : PLPLx_SetPeriodicCondition(period)  

Parameters
----------
* `period` :  
    (double) : value of the period  
";

%feature("docstring") PLPLx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : PLPLx_SetFrictionModel(cflag)  

Parameters
----------
* `cflag` :  
    (char) : model to use ('min', 'max' or 'ave')  
";

%feature("docstring") PLPLx_SetBigPolygTolerance "

python usage : PLPLx_SetBigPolygTolerance(tol)  

Parameters
----------
* `period` :  
    (double) : value of the tolerance  
";

%feature("docstring") PLPLx_ComputeStress "

compute stress  

python usage : PLPLx_ComputeStress()  
";

%feature("docstring") PLPLx_ComputeBetai "

compute equivalent damage parameter  

python usage : PLPLx_ComputeBetai()  
";

%feature("docstring") PLPLx_ComputeCZMEnergy "

compute and decompose local contact energy with CZM law  

python usage : PLPLx_ComputeCZMEnergy()  
";

%feature("docstring") PLPLx_CleanMemory "

Free all memory allocated within PLPLx module.  

python usage : PLPLx_CleanMemory()  
";

%feature("docstring") PLPLx_GetCZMEnergy "

Get the CZM energy of a given contact.  

python usage energy = PLPLx_GetCZMEnergy(icdan)  

Parameters
----------
* `icdan(int)` :  
    : index of the PLPLx contact  

Returns
-------
energy(double[4]) : energy value  
";

%feature("docstring") PLPLx_UseNcDetection "

chooses contact detection methode between non-convex shapes  

python usage : PLPLx_UseNcDetection()  
";

%feature("docstring") PLPLx_ShrinkPolygFaces "

Shrink the face of the polygon for the detection.  

python usage : PLPLx_ShrinkPolygFaces(shrink)  

Parameters
----------
* `shrink` :  
    (real) :  
";

