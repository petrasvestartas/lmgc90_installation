
// File: wrap__poroMAILx_8h.xml

%feature("docstring") poroMAILx_LoadModels "

load from MAILx and models  

python usage : poroMAILx_LoadModels()  
";

%feature("docstring") poroMAILx_LoadBehaviours "

load from bulk_behav  

python usage : pordMAILx_LoadBehaviours()  
";

%feature("docstring") poroMAILx_PushProperties "

declares to models couple (model,behav)  

python usage : poroMAILx_PushProperties()  
";

%feature("docstring") poroMAILx_ReadDrivenDof "

Read DRV_DOF.DAT.  

python usage : poroMAILx_ReadDrivenDof()  
";

%feature("docstring") poroMAILx_WriteDrivenDof "

Write DRV_DOF.OUT.  

python usage : poroMAILx_WriteDrivenDof()  
";

%feature("docstring") poroMAILx_ReadIniDof "

Read DOF.INI.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

python usage : poroMAILx_ReadIniDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : which DOF file to read  
";

%feature("docstring") poroMAILx_ReadIniMecaDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniMecaDof last call  

python usage : poroMAILx_ReadIniMecaDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : which DOF file to read  
";

%feature("docstring") poroMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : poroMAILx_ReadIniGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : which GPV file to read  
";

%feature("docstring") poroMAILx_ReadIniMecaGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read Else : OUTBOX/GPV.OUT.num is read, num
being the parameter used in TimeEvolution_ReadIniMecaGPV last call  

python usage : poroMAILx_ReadIniMecaGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : which GPV file to read  
";

%feature("docstring") poroMAILx_WriteLastDof "

Write ascii DOF.LAST file.  

python usage : poroMAILx_WriteLastDof()  
";

%feature("docstring") poroMAILx_ComputeMass "

compute elementary mass and inertia of bodies  

python usage : poroMAILx_ComputeMass()  
";

%feature("docstring") poroMAILx_ComputeFext "

compute elementary external forces  

python usage : poroMAILx_ComputeFext()  
";

%feature("docstring") poroMAILx_ComputeBulk "

compute elementary stiffness  

python usage : poroMAILx_ComputeBulk()  
";

%feature("docstring") poroMAILx_ComputeDamping "

compute elemenatry damping  

python usage : poroMAILx_ComputeDamping()  
";

%feature("docstring") poroMAILx_AssembKT "

assembles matrice  

python usage : poroMAILx_AssembKT()  
";

%feature("docstring") poroMAILx_AssembRHS "

assembles RHS  

python usage : poroMAILx_AssembRHS()  
";

%feature("docstring") poroMAILx_ComputeFreeVelocity "

computes free motion (without contact contribution)  

python usage : poroMAILx_ComputeFreeVelocity()  
";

%feature("docstring") poroMAILx_ComputeDof "

computes motion (free + contact)  

python usage : poroMAILx_ComputeDof()  
";

%feature("docstring") poroMAILx_DisplayOutDof "

Display body degrees of freedom.  

python usage : poroMAILx_DisplayOutDof()  
";

%feature("docstring") poroMAILx_UpdateDof "

update begin dof with current dof  

python usage : poroMAILx_UpdateDof()  
";

%feature("docstring") poroMAILx_UpdateBulk "

update begin elementary fields with current elementary fields  

python usage : poroMAILx_UpdateBulk(void)  
";

%feature("docstring") poroMAILx_ComputeGrad "

apply elementary fields gradient  

python usage : poroMAILx_ComputeGrad(void)  
";

%feature("docstring") poroMAILx_GetBodyVector "

Get a copy of a vector of a given body.  

*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"VbALE\": fluid velocity at beginning of time step  
*   \"V_ALE\": fluid velocity at beginning of time step  
*   \"Vaux_\": working array for velocity  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Fext_\": external forces  
*   \"Fint_\": internal forces  
*   \"Pbeg_\": pressure at beginning of time step  
*   \"P____\": pressure in computed configuration  
*   \"Qext_\": external fluxes  
*   \"Qint_\": internal luxces  
*   \"NodId\":  

Possible values for datatype field are \"X____\", \"Xbeg_\", \"V____\",
\"Vbeg_\", \"Vaux_\", \"Reac_\", \"Vfree\", \"Fext_\", \"Fint_\"  

Python usage : vector = poroMAILx_GetBodyVector(datatype, ibdyty)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to get  
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired vector  
";

%feature("docstring") poroMAILx_GetNbNodes "

Get the number of nodes of a poroMAILx.  

python usage : nb_nodes = poroMAILx_GetNbNodes(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the poroMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a poroMAILx  
";

%feature("docstring") poroMAILx_GetNbElements "

Get the number of elements of a poroMAILx.  

python usage : nb_elements = poroMAILx_GetNbElements(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the poroMAILx  

Returns
-------
nb_nodes (integer) : number of elements of a poroMAILx  
";

%feature("docstring") poroMAILx_IncrementStep "

correction of the configuration parameter using the theta-method  

python usage : poroMAILx_IncrementStep()  
";

%feature("docstring") poroMAILx_WithoutRenumbering "

skip renumbering of the unknowns using a rcc method  

python usage : poroMAILx_WithoutRenumbering()  
";

%feature("docstring") poroMAILx_BandStorage "

use band matrix  

python usage : poroMAILx_BandStorage()  
";

%feature("docstring") poroMAILx_SparseStorage "

use sparse matrix  

python usage : poroMAILx_SparseStorage()  
";

%feature("docstring") poroMAILx_ExplodedStorage "

use element by element matrix  

python usage : poroMAILx_ExplodedStorage()  
";

%feature("docstring") poroMAILx_DiagonalStorage "

use diagonal matrix  

python usage : poroMAILx_DiagonalStorage()  
";

%feature("docstring") poroMAILx_SkylineStorage "

use skyline matrix  

python usage : poroMAILx_SkylineStorage()  
";

%feature("docstring") poroMAILx_FullStorage "

use full matrix  

python usage : poroMAILx_FullStorage()  
";

%feature("docstring") poroMAILx_SymmetricShape "

assume matrix is symmetrical  

python usage : poroMAILx_SymmetricShape()  
";

%feature("docstring") poroMAILx_UnspecifiedShape "

does not assume any thing on matrix shape  

python usage : poroMAILx_UnspecifiedShape()  
";

%feature("docstring") poroMAILx_SetMecaScalarFieldByNode "

Update an external field on a given body.  

python usage : poroMAILx_SetMecaScalarFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
  
 You need to set this field in your models.dat  
";

%feature("docstring") poroMAILx_SetTherScalarFieldByNode "

Update an external field on a given body.  

python usage : poroMAILx_SetTherieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
  
 You need to set this field in your models.dat  
";

%feature("docstring") poroMAILx_SetMecaScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetMecaScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") poroMAILx_SetTherScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetTherScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") poroMAILx_GetMecaScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetMecaScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the concern body  
* `iblmty` :  
    (integer) : id of the concern element  
* `name` :  
    (string) : name of the desired scalar field  

Returns
-------
f_rank (integer) : rank of the corresponding scalar field  
";

%feature("docstring") poroMAILx_GetMecaVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetMecaVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the concern body  
* `iblmty` :  
    (integer) : id of the concern element  
* `name` :  
    (string) : name of the desired vector field  

Returns
-------
f_rank (integer) : rank of the corresponding vector field  
";

%feature("docstring") poroMAILx_GetTherScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetTherScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the concern body  
* `iblmty` :  
    (integer) : id of the concern element  
* `name` :  
    (string) : name of the desired scalar field  

Returns
-------
f_rank (integer) : rank of the corresponding scalar field  
";

%feature("docstring") poroMAILx_GetTherVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetTherVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the concern body  
* `iblmty` :  
    (integer) : id of the concern element  
* `name` :  
    (string) : name of the desired vector field  

Returns
-------
f_rank (integer) : rank of the corresponding vector field  
";

%feature("docstring") poroMAILx_SetMecaVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetMecaVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetVectorFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetTherVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetTherVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") poroMAILx_LoadALE "

Apply an ALE Formulation in Fluid zone.  

python usage : poroMAILx_LoadALE(IdBody)  
";

%feature("docstring") poroMAILx_PutBodyVector "

Set a vector of a given body.  

*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"VbALE\": fluid velocity at beginning of time step  
*   \"V_ALE\": fluid velocity at beginning of time step  
*   \"Raux_\": working array for reaction  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Fext_\": external forces  
*   \"Fint_\": internal forces  
*   \"Qext_\": external fluxes  
*   \"Qint_\": internal luxces  
*   \"Pbeg_\": pressure at beginning of time step  
*   \"P____\": pressure in computed configuration  

python usage : poroMAILx_PutBodyVector(datatype, ibdyty, matrix)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `ibdyty` :  
    (integer) : rank of body  
* `matrix` :  
    (double array) : the new values  
";

%feature("docstring") poroMAILx_ComputeResidueNorm "

computes the norm of the residue  

python usage : norm = poroMAILx_ComputeResidueNorm()  

Returns
-------
norm (double) : Residue Norm  
";

%feature("docstring") poroMAILx_GetStress "

Get a copy of a stress of a given body.  

Python usage : stress = poroMAILx_GetStress(ibdyty,required_field=0)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `required_field` :  
    (integer) : required additional field  

Returns
-------
matrix_out (double 2D-array) : the desired stress  
";

%feature("docstring") poroMAILx_GetStrain "

Get a copy of a strain of a given body.  

Python usage : strain = poroMAILx_GetStrain(ibdyty, required_field=0)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `required_field` :  
    (integer) : required additional field  

Returns
-------
strain (double 2D-array) : the desired strain  
";

%feature("docstring") poroMAILx_ComputeContactDetectionConfiguration "

compute the contact detection configuration  

python usage : poroMAILx_ComputeContactDetectionConfiguration()  
";

%feature("docstring") poroMAILx_SetPreconAllBodies "

ask for precomputation of the W matrix on support node dofs of contactors for
all bodies. Assumes bulk behaviour is linear.  

python usage : poroMAILx_SetPreconAllBodies()  
";

%feature("docstring") poroMAILx_ComputePreconW "

compute the precon W on precon bodies  

python usage : poroMAILx_ComputePreconW()  
";

%feature("docstring") poroMAILx_GetNbPoroMAILx "

Get the number of poroMAILx.  

python usage : nb_poroMAILx = poroMAILx_GetNbPoroMAILx()  

Returns
-------
nb_poroMAILx (integer) : number of poroMAILx  
";

%feature("docstring") poroMAILx_GetCoor "

return node coordinates of idBody  

python usage : array = poroMAILx_GetCoor(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") poroMAILx_GetAll "

return poro mechanical data computed for idBody  

python usage : array = poroMAILx_GetAll(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : poro mechanical data  
";

%feature("docstring") poroMAILx_GetGrad "

Get a copy of a grad P of a given body.  

Python usage : grad = poroMAILx_GetGrad(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
grad (double 2D-array) : the desired grad  
";

%feature("docstring") poroMAILx_GetFlux "

Get a copy of a Darcy Flux of a given body.  

Python usage : flux = poroMAILx_GetFlux(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
flux (double 2D-array) : the desired flux  
";

%feature("docstring") poroMAILx_GetInternal "

return internal mechanical data computed for idBody  

python usage : array = poroMAILx_GetInternal(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical internal data  
";

%feature("docstring") poroMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = poroMAILx_GetConnectivity(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") poroMAILx_SetVlocyDrivenDof "

Apply Drv Dof on a given body.  

python usage : poroMAILx_SetVlocyDrivenDof(IdBody, f_dof, f_node, f_value)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_dof` :  
    (integer) : dof of the concern node  
* `f_node` :  
    (integer) : node  
* `f_value` :  
    (double) : value of the drvdof  
";

%feature("docstring") poroMAILx_AddFieldLoad "

Add elementary load through a nodal external field on a given body.  

python usage : poroMAILx_AddFieldLoad(IdBody, Ideriv, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") poroMAILx_WriteOutDof "

Write ascii DOF.OUT file. Can be activate only each N step.  

python usage : poroMAILx_WriteOutDof()  
";

%feature("docstring") poroMAILx_PostModels "

load from MAILx and models for post  

python usage : poroMAILx_PostModels()  
";

%feature("docstring") poroMAILx_CleanMemory "

Free all memory allocated within poroMAILx module.  

python usage : poroMAILx_CleanMemory()  
";

%feature("docstring") poroMAILx_CheckProperties "

check if model and material are matching ; set material parameter if external
model  

python usage : poroMAILx_CheckProperties()  
";

%feature("docstring") poroMAILx_GetNbGpByElem "

Get the list of finite elements for porox models and the associated number of
Gauss Points for MECA and THER physics.  

Here memory is allocated within lmgc90 so that the pointer can be freely
modified by third parties without nasty effect on lmgc90 functioning.  

python usage : names, meca_nb, ther_nb = poroMAILx_GetNbGpByElem()  

Returns
-------
names (string list) : list of the finite elements meca_nb (integer list): list
of the number of Gauss Points for MECA ther_nb (integer list): list of the
number of Gauss Points for THER  
";

