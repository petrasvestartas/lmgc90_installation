
// File: wrap__DKDKx_8h.xml

%feature("docstring") DKDKx_SelectProxTactors "

contact detection between CLxxx and JCxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKDKx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKDKx_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKDKx_SmoothForceComputation()  
";

%feature("docstring") DKDKx_UseVaVDetection "

allow to increase the number of contact for a pair cd/an.  

python usage : DKDKx_UseVaVDetection(nb)  

Parameters
----------
* `nb` :  
    (integer) : number of contact points for a couple (cd,an)  
";

%feature("docstring") DKDKx_WriteLastVlocRloc "

write last local values of all DKDKx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKx_WriteLastVlocRloc()  
";

%feature("docstring") DKDKx_WriteOutVlocRloc "

write local values of all DKDKx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKx_WriteOutVlocRloc()  
";

%feature("docstring") DKDKx_DisplayOutVlocRloc "

display local values of all DKDKx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKDKx_DisplayOutVlocRloc()  
";

%feature("docstring") DKDKx_DisplayProxTactors "

display contacts  

python usage : DKDKx_DisplayProxTactors()  
";

%feature("docstring") DKDKx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKDKx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") DKDKx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKDKx_SetPeriodicCondition(period)  

Parameters
----------
* `period` :  
    (double) : value of the period  
";

%feature("docstring") DKDKx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : DKDKx_SetFrictionModel(cflag)  

Parameters
----------
* `cflag` :  
    (char) : model to use ('min', 'max' or 'ave')  
";

%feature("docstring") DKDKx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKDKx_SetSurfaceSectors(nbsect)  

Parameters
----------
* `nbsect` :  
    (integer) : number of sectors  
";

%feature("docstring") DKDKx_UpdateSurfaceEnergySector "

update surface energy sector  

python usage : DKDKx_UpdateSurfaceEnergySector()  
";

%feature("docstring") DKDKx_ComputeStress "

update surface energy sector  

python usage : DKDKx_ComputeStress()  
";

%feature("docstring") DKDKx_ComputeBetai "

compute equivalent damage parameter  

python usage : DKDKx_ComputeBetai()  
";

%feature("docstring") DKDKx_ComputeCZMEnergy "

compute and decompose local contact energy with CZM law  

python usage : DKDKx_ComputeCZMEnergy()  
";

%feature("docstring") DKDKx_CleanMemory "

Free all memory allocated within DKDKx module.  

python usage : DKDKx_CleanMemory()  
";

%feature("docstring") DKDKx_GetCZMEnergy "

Get the CZM energy of a given contact.  

python usage : energy = DKDKx_GetCZMEnergy(icdan)  

Parameters
----------
* `icdan(int)` :  
    : index of the DKDKx contact  

Returns
-------
energy(double[4]) : energy value  
";

