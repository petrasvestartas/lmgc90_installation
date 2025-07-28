
// File: wrap__mecaMAILx_8h.xml

%feature("docstring") mecaMAILx_WithoutRenumbering "

skip renumbering of the unknowns using a rcc method  

python usage : mecaMAILx_WithoutRenumbering()  
";

%feature("docstring") mecaMAILx_BandStorage "

use band matrix  

python usage : mecaMAILx_BandStorage()  
";

%feature("docstring") mecaMAILx_SparseStorage "

use sparse matrix  

python usage : mecaMAILx_SparseStorage()  
";

%feature("docstring") mecaMAILx_ExplodedStorage "

use element by element matrix  

python usage : mecaMAILx_ExplodedStorage()  
";

%feature("docstring") mecaMAILx_DiagonalStorage "

use diagonal matrix  

python usage : mecaMAILx_DiagonalStorage()  
";

%feature("docstring") mecaMAILx_SkylineStorage "

use skyline matrix  

python usage : mecaMAILx_SkylineStorage()  
";

%feature("docstring") mecaMAILx_FullStorage "

use full matrix  

python usage : mecaMAILx_FullStorage()  
";

%feature("docstring") mecaMAILx_SymmetricShape "

assume matrix is symmetrical  

python usage : mecaMAILx_SymmetricShape()  
";

%feature("docstring") mecaMAILx_UnspecifiedShape "

does not assume any thing on matrix shape  

python usage : mecaMAILx_UnspecifiedShape()  
";

%feature("docstring") mecaMAILx_GetNbMecaMAILx "

Get the number of mecaMAILx.  

python usage : nb_mecaMAILx = mecaMAILx_GetNbMecaMAILx()  

Returns
-------
nb_mecaMAILx (integer) : number of mecaMAILx  
";

%feature("docstring") mecaMAILx_GetNbNodes "

Get the number of nodes of a mecaMAILx.  

python usage : nb_nodes = mecaMAILx_GetNbNodes(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the mecaMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a mecaMAILx  
";

%feature("docstring") mecaMAILx_GetNbElements "

Get the number of elements of a mecaMAILx.  

python usage : nb_elements = mecaMAILx_GetNbElements(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the mecaMAILx  

Returns
-------
nb_nodes (integer) : number of elements of a mecaMAILx  
";

%feature("docstring") mecaMAILx_GetNbGp "

Get the number of Gauss points of an element of a mecaMAILx.  

python usage : nb_gp = mecaMAILx_GetNbElements(ibdyty, iblmty)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the mecaMAILx  
* `iblmty` :  
    (integer) : id of the element  

Returns
-------
nb_gp (integer) : number of Gauss point of an element of a mecaMAILx  
";

%feature("docstring") mecaMAILx_SetPreconBody "

ask for precomputation of the W matrix on support node dofs of contactors for
one body. Assumes bulk behaviour is linear.  

python usage : mecaMAILx_SetPreconBody(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of body to set precon  
";

%feature("docstring") mecaMAILx_SetPreconAllBodies "

ask for precomputation of the W matrix on support node dofs of contactors for
all bodies. Assumes bulk behaviour is linear.  

python usage : mecaMAILx_SetPreconAllBodies()  
";

%feature("docstring") mecaMAILx_ComputePreconW "

compute the precon W on precon bodies  

python usage : mecaMAILx_ComputePreconW()  
";

%feature("docstring") mecaMAILx_InitPreconW "

initialize an empty precon W  

python usage : mecaMAILx_InitPreconW()  
";

%feature("docstring") mecaMAILx_PutPreconW "

push a column of precon W  

python usage : mecaMAILx_PutPreconW(ivalue1, ivalue2, ivalue3, vect)  

Parameters
----------
* `ivalue1` :  
    (integer) : body number  
* `ivalue2` :  
    (integer) : node number  
* `ivalue3` :  
    (integer) : dof number  
* `vect` :  
    (double) : column  
";

%feature("docstring") mecaMAILx_GetNodesPrecon "

Get the list of preconditionned nodes of a mecaMAILx body.  

Here memory is allocated within lmgc90 so that the pointer can be freely
modified by third parties without nasty effect on lmgc90 functioning.  

python usage : precon_list = mecaMAILx_GetNodesPrecon(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the desired mecaMAILx  

Returns
-------
precon_list (integer list) : list of the preconditionned nodes  
";

%feature("docstring") mecaMAILx_SetCoroAllBodies "

ask for corotationnal computation of the W matrix. Assumes bulk behaviour is
linear.  

python usage : mecaMAILx_SetCoroAllBodies()  
";

%feature("docstring") mecaMAILx_SetCoroBody "

ask for corotationnal computation of the W matrix of a given body. Assumes bulk
behaviour is linear.  

python usage : mecaMAILx_SetCoroBody(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of body to set coro  
";

%feature("docstring") mecaMAILx_SetTolCoro "

set the admssible tolerance on rigid body velocity computed by deformable model  

python usage : mecaMAILx_SetTolCoro(tol)  

Parameters
----------
* `tol` :  
    (double) : tolerance  
";

%feature("docstring") mecaMAILx_SetRigidAllBodies "

ask for rigid computation of the W matrix. Assumes bulk behaviour is linear.  

python usage : mecaMAILx_SetRigidAllBodies()  
";

%feature("docstring") mecaMAILx_SetRigidBody "

ask for rigid computation of the W matrix of a given body. Assumes bulk
behaviour is linear.  

python usage : mecaMAILx_SetRigidBody(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of body to compute as a rigid  
";

%feature("docstring") mecaMAILx_SkipDeformableComputationAllBodies "

avoid deformable part computation of a deformable body declared as rigid  

python usage : mecaMAILx_SkipDeformableComputationAllBodies()  
";

%feature("docstring") mecaMAILx_SkipDeformableComputationBody "

avoid deformable part computation of a given deformable body declared as rigid  

python usage : mecaMAILx_SkipDeformableComputationBody(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of body to compute without deformation  
";

%feature("docstring") mecaMAILx_BuildRigidBodies "

computes internal matrices for rigid description  

python usage : mecaMAILx_BuildRigidBodies()  
";

%feature("docstring") mecaMAILx_IsRigid "

return 1 if a given body is rigid/coro, 0 otherwize  

python usage : rigid = mecaMAILx_IsRigid(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body we want visibility  

Returns
-------
rigid (integer) : 1 if body is visible, 0 otherwize  
";

%feature("docstring") mecaMAILx_GetRigidFrame "

return an inertia frame matrix  

Possible values for datatype field are \"RFbeg\", \"RF___\", \"RFTT_ (stands for
Rigid Frame)  

python usage : mat = mecaMAILx_GetRigidFrame(datatype, ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body  

Returns
-------
vec (float matrix) : frame matrix (beg, current or TT)  
";

%feature("docstring") mecaMAILx_GetRigidCoorTT "

return TT center of inertia coordinates  

python usage : vec = mecaMAILx_GetRigidCoorTT(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body  

Returns
-------
vec (float vector) : TT center of inertia coordinates  
";

%feature("docstring") mecaMAILx_GetRigidCooref "

return ref center of inertia coordinates  

python usage : vec = mecaMAILx_GetRigidCooref(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body  

Returns
-------
vec (float vector) : ref center of inertia coordinates  
";

%feature("docstring") mecaMAILx_SetRVDrivenDofs "

declares rigid velocity dof as driven  

python usage : mecaMAILx_SetRVDrivenDofs(idbody,vector_in)  

Parameters
----------
* `idbody` :  
    (integer) : id of the body  
* `vector` :  
    (integer) : list of driven dofs  
";

%feature("docstring") mecaMAILx_SetRVDrivenDofValue "

set the value of rigid velocity dof value  

python usage : mecaMAILx_SetRVDrivenDofValue(idbody,iddof,rv)  

Parameters
----------
* `idbody` :  
    (integer) : id of the body  
* `iddof` :  
    (integer) : id of dof  
* `rv` :  
    (float) : value  
";

%feature("docstring") mecaMAILx_PutBodyRVector "

Set a vector of a coro or rigid mecaMAILx body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Raux_\": working array for reaction force  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  

uses copy, and in case fo Fext, the operation is not just setting but adding  

python usage : mecaMAILx_PutBodyRVector(datatype, ibdyty, vector)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to set  
* `ibdyty` :  
    (integer) : rank of the RBDY3  
* `vector` :  
    (double array) : the new value of the vector  
";

%feature("docstring") mecaMAILx_GetBodyRVector "

Get a copy of a vector of a mecaMAILx body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coorb\": coordinates at beginning of time step  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"XTT__\": cumulated displacements over time in detection configuration  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Raux_\": working array for reaction force  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  
*   \"Fint_\": internal forces  

python usage : vector = mecaMAILx_GetBodyRVector(datatype, ibdyty)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to get  
* `ibdyty` :  
    (integer) : rank of the RBDY3  

Returns
-------
vector (double array) : output vector  
";

%feature("docstring") mecaMAILx_PutBodyVector "

Set a vector of a given body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coor_\": coordinates in computed configuration  
*   \"Coorb\": coordinates at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Raux_\": working array for reaction force  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  
*   \"Fint_\": internal forces  

python usage : mecaMAILx_PutBodyVector(datatype, ibdyty, matrix)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `ibdyty` :  
    (integer) : rank of body  
* `matrix` :  
    (double array) : the new value  
";

%feature("docstring") mecaMAILx_GetBodyVector "

Get a copy of a vector of a given body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coor_\": coordinates in computed configuration  
*   \"Coorb\": coordinates at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Vaux_\": working array for velocity  
*   \"Vfree\": velocity free of contacts  
*   \"Reac_\": contact reaction force  
*   \"Raux_\": working array for reaction force  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  
*   \"Fint_\": internal forces  

Python usage : vector = mecaMAILx_GetBodyVector(datatype, ibdyty)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to get  
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired data  
";

%feature("docstring") mecaMAILx_GetMaterials "

Get a copy of a the elements' material vector of a given body.  

Python usage : materials = mecaMAILx_GetMaterials(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
vector (double 1D-array) : the material index of elements  
";

%feature("docstring") mecaMAILx_GetStress "

Get a copy of the smoothed nodal stress (Cauchy) of a given body: 2D
Sxx,Syy,Sxy,Szz,Svm | 3D Sxx,Sxy,Syy,Sxz,Syz,Szz,Svm.  

Python usage : stress = mecaMAILx_GetStress(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
stress (double 2D-array) : nodal stress of the desired body  
";

%feature("docstring") mecaMAILx_GetStrain "

Get a copy of the smoothed nodal strain (Almansi) of a given body: 2D
Exx,Eyy,Exy,Ezz,J | 3D Exx,Exy,Eyy,Exz,Eyz,Ezz,J.  

Python usage : strain = mecaMAILx_GetStrain(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
strain (double 2D-array) : nodal strain of the desired body  
";

%feature("docstring") mecaMAILx_GetInternalVariables "

Get a copy of the smoothed nodal internal variables (2D:10 ; 3D:57)  

Python usage : strain = mecaMAILx_GetInternalVariables(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
strain (double 2D-array) : nodal internal variables of the desired body  
";

%feature("docstring") mecaMAILx_GetElementStress "

Get a copy of the mean stress (Cauchy) of a given body: 2D Sxx,Syy,Sxy,Szz,Svm |
3D Sxx,Sxy,Syy,Sxz,Syz,Szz,Svm.  

Python usage : stress = mecaMAILx_GetElementStress(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
stress (double 2D-array) : nodal stress of the desired body  
";

%feature("docstring") mecaMAILx_PushProperties "

gives to model the couple of model,behavior used at gauss point  

python usage : mecaMAILx_PushProperties()  
";

%feature("docstring") mecaMAILx_UseNewPPSet "

each gauss point will have its own property set (necessary in multi physics)  

python usage : mecaMAILx_UseNewPPSet()  
";

%feature("docstring") mecaMAILx_ComputeFreeVelocity "

computes free velocity of a list of bodies  

python usage : mecaMAILx_ComputeFreeVelocity(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute free velocity if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_AssembKT "

assemble pseudo mass matrix and apply drvdof of a list of bodies  

python usage : mecaMAILx_AssembKT(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to assemble pseudo mass matrix and apply
    drvdof if omitted works on all objects  
";

%feature("docstring") mecaMAILx_OnlyAssembKT "

assemble pseudo mass matrix of a list of bodies  

python usage : mecaMAILx_OnlyAssembKT(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to assemble pseudo mass matrix if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ApplyDrvDofKT "

apply drvdof pseudo mass matrix  

python usage : mecaMAILx_ApplyDrvDofKT(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to apply drvdof on pseudo mass matrix if
    omitted works on all objects  
";

%feature("docstring") mecaMAILx_AssembRHS "

assembles right hand side of a list of bodies  

python usage : mecaMAILx_AssembRHS(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to assemble right hand side if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeResidueNorm "

computes the norm of the residue of a list of bodies  

python usage : norm = mecaMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute the norm of the residue if
    omitted works on all objects  

Returns
-------
norm (double) : Residue Norm  
";

%feature("docstring") mecaMAILx_ComputeBulk "

computes elementary stiffness and viscosity matrices and internal forces of a
list of bodies  

python usage : mecaMAILx_ComputeBulk(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute stiffness and viscosity
    matrices and internal forces if omitted works on all objects  
";

%feature("docstring") mecaMAILx_ComputeField "

computes elementary fields of a list of bodies  

python usage : mecaMAILx_ComputeField(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeFint "

computes elementary internal forces of a list of bodies  

python usage : mecaMAILx_ComputeFint(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute internal forces if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_UpdateBulk "

update begin elementary fields with current elementary fields of a list of
bodies  

python usage : mecaMAILx_UpdateBulk(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_UpdateDof "

update begin d.o.f. with current d.o.f. of a list of bodies  

python usage : mecaMAILx_UpdateDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to update current d.o.f if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_ComputeDof "

computes the current d.o.f knowing all the forces (free + contact) of a list of
bodies  

python usage : mecaMAILx_ComputeDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute current d.o.f if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_IncrementStep "

initializes the current d.o.f and some driven d.o.f values  

python usage : mecaMAILx_IncrementStep()  
";

%feature("docstring") mecaMAILx_ComputeFext "

compute elementary external forces of a list of bodies  

python usage : mecaMAILx_ComputeFext(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute external forces if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeMass "

compute elementary mass and inertia of a list of bodies  

python usage : mecaMAILx_ComputeMass(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute mass and inertia if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_FatalDamping "

set to 0 current velocities of a list of bodies  

This keyword must be between the ComputeDof and UpdateDof ones.  

python usage : mecaMAILx_FatalDamping(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to reset current velocity if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_CheckEquilibriumState "

Check if the bodies riches an equilibrium state (velocities almost equal to 0)  

python usage : iconv = mecaMAILx_CheckEquilibriumState()  

Returns
-------
iconv (boolean) : True if in equilibrium state  
";

%feature("docstring") mecaMAILx_SetEquilibriumNorm "

set the norm for CheckEquilibriumState  

Type of check test:  

*   Qvlcy : quadratic norm of velocy  
*   Maxm : maximum norm of velocy  

python usage : mecaMAILx_SetEquilibriumNorm(checktype, tol)  

Parameters
----------
* `checktype` :  
    (char[5]) : type of check test  
* `tol` :  
    (double) : tolerance  
";

%feature("docstring") mecaMAILx_ReadDrivenDof "

Read DRV_DOF.DAT.  

python usage : mecaMAILx_ReadDrivenDof()  
";

%feature("docstring") mecaMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : mecaMAILx_ReadIniGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : which GPV file to read  
";

%feature("docstring") mecaMAILx_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

python usage : mecaMAILx_ReadIniDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : which DOF file to read  
";

%feature("docstring") mecaMAILx_LoadBehaviours "

load behaviours from bulk_behav  

python usage : mecaMAILx_LoadBehaviours()  
";

%feature("docstring") mecaMAILx_LoadModels "

load models from models  

python usage : mecaMAILx_LoadModels()  
";

%feature("docstring") mecaMAILx_WriteDrivenDof "

Write DRV_DOF.OUT.  

python usage : mecaMAILx_WriteDrivenDof()  
";

%feature("docstring") mecaMAILx_WriteLastDof "

Write ascii DOF.LAST file.  

python usage : mecaMAILx_WriteLastDof()  
";

%feature("docstring") mecaMAILx_WriteOutDof "

Write ascii DOF.OUT file. Can be activate only each N step.  

python usage : mecaMAILx_WriteOutDof()  
";

%feature("docstring") mecaMAILx_DisplayOutDof "

Display body degrees of freedom.  

python usage : mecaMAILx_DisplayOutDof()  
";

%feature("docstring") mecaMAILx_DisplayBulkElement "

Display fields of a bulk element.  

Parameters
----------
* `IdBody` :  
    (int) : id of the concern body  
* `IdElem` :  
    (int) : id of the concern element python usage :
    mecaMAILx_DisplayBulkElement(IdBody,IdElem)  
";

%feature("docstring") mecaMAILx_WriteLastRnod "

Write ascii Rnod.LAST file of a list of bodies.  

python usage : mecaMAILx_WriteLastRnod(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write in Rnod.LAST if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_WriteOutRnod "

Write ascii Rnod.OUT file of a list of bodies. Can be activat only each N step.  

python usage : mecaMAILx_WriteOutRnod(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write in Rnod.OUT if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_DisplayOutRnod "

Display body forces of a list of bodies.  

python usage : mecaMAILx_DisplayOutRnod(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to display body forces if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_WriteLastNodalForces "

Write ascii Rnod.LAST file of a list of bodies.  

This function is almost like WriteLastRnod, but write also internal and inertial
forces.  

python usage : mecaMAILx_WriteLastNodalForces(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write in Rnod.LAST if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_WriteOutNodalForces "

Write ascii Rnod.OUT file of a list of bodies. Can be activat only each N step.  

This function is almost like WriteOutRnod, but write also internal and inertial
forces.  

python usage : mecaMAILx_WriteOutNodalForces(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write in Rnod.OUT if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_DisplayOutNodalForces "

Display computed nodal forces of a list of bodies.  

python usage : mecaMAILx_DisplayOutNodalForces(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to display body forces if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_GetScalarFieldRank "

Get the rank of scalar field of an element of a body from its name.  

python usage : f_rank = mecaMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the concern body  
* `iblmty` :  
    (integer) : id of the concern element  
* `name` :  
    (string) : name of the desired field  

Returns
-------
f_rank (integer) : rank of the corresponding field  
";

%feature("docstring") mecaMAILx_SetScalarFieldByNode "

Update elementary scalar field through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

python usage : mecaMAILx_SetScalarFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
  
 You need to declare this field in your MODELS.DAT  
";

%feature("docstring") mecaMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

python usage : mecaMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
  
 You need to declare this field in your MODELS.DAT  
";

%feature("docstring") mecaMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = mecaMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

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

%feature("docstring") mecaMAILx_SetVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : mecaMAILx_SetVectorFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") mecaMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : mecaMAILx_SetVectorFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") mecaMAILx_Terminate "

Stop job properly.  

python usage : mecaMAILx_Terminate()  
";

%feature("docstring") mecaMAILx_ComputeOrthoFrame "

Use user routine to compute the ortho frame of a list of bodies.  

This method uses a routine define by the user in user.f90  

python usage : mecaMAILx_ComputeOrthoFrame(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute ortho frame with user routine
    if omitted works on all objects  
";

%feature("docstring") mecaMAILx_ComputeUserField "

Use user routine to compute a field at gp.  

python usage : mecaMAILx_ComputeUserField(ifield, i_list)  

Parameters
----------
* `ifield` :  
    (integer) : id of the field to compute  
* `i_list` :  
    (list of integer) : list of bodies to compute user fields on if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_SetVisible "

set visible a given mecaMAILx  

python usage : mecaMAILx_SetVisible(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the mecaMAILx  
";

%feature("docstring") mecaMAILx_SetInvisible "

rended a given mecaMAILx invisible  

python usage : mecaMAILx_SetInvisible(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the mecaMAILx  
";

%feature("docstring") mecaMAILx_IsVisible "

return if a given body visible  

python usage : visible = mecaMAILx_IsVisible(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") mecaMAILx_ComputeRayleighDamping "

compute the Rayleigh damping: C=alpha*M+beta*K of a list of bodies  

python usage : mecaMAILx_ComputeRayleighDamping(alpha,beta,i_list)  

Parameters
----------
* `alpha` :  
    (real) : damping value  
* `beta` :  
    (real) : damping value  
* `i_list` :  
    (list of integer) : list of bodies to compute Rayleigh damping if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeRayleighDampingDiscreteElement "

set damping for discrete FE element of a list of bodies  

python usage : mecaMAILx_ComputeRayleighDampingDiscreteElement(damp, i_list)  

Parameters
----------
* `ref_size` :  
    (real) : damping value  
* `i_list` :  
    (list of integer) : list of bodies to compute damping for discrete FE
    element if omitted works on all objects  
";

%feature("docstring") mecaMAILx_GetNodeCoorTT "

return TT node coordinates  

python usage : vec = mecaMAILx_GetNodeCoorTT(ibdyty,inodty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body  
* `inodty(integer)` :  
    : id of the node  

Returns
-------
vec (float vector) : TT node coordinates  
";

%feature("docstring") mecaMAILx_GetNodeCooref "

return ref node coordinates  

python usage : vec = mecaMAILx_GetNodeCoorref(ibdyty,inodty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body  
* `inodty(integer)` :  
    : id of the node  

Returns
-------
vec (float vector) : ref node coordinates  
";

%feature("docstring") mecaMAILx_GetBodyMatrix "

Get a copy of a matrix of a given body.  

Possible values for datatype field are \"mass_\", \"stiff\", \"damp_\"  

Python usage : matrix = mecaMAILx_GetBodyMatrix(datatype, ibdyty)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the matrix to get  
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
matrix (double array) : the desired matrix  
";

%feature("docstring") mecaMAILx_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = mecaMAILx_getDrvVlocy(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the mecaMAILx  
* `drvdof_indices` :  
    (integer array) : indices list of driven dof  
* `drvdof_values` :  
    (real array) : values of the driven dof  
";

%feature("docstring") mecaMAILx_computeDrvVlocy "

Compute the value of the driven velocity of a body a current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : mecaMAILx_computeDrvVlocy(ibdyty, values)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the mecaMAILx  
* `values` :  
    (double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") mecaMAILx_SetVlocyDrivenDof "

Apply Drv Dof on a given body.  

python usage : mecaMAILx_SetVlocyDrivenDof(IdBody, f_dof, f_node, f_value)  

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

%feature("docstring") mecaMAILx_ComputeContactDetectionConfiguration "

compute the contact detection configuration of a list of bodies  

python usage : mecaMAILx_ComputeContactDetectionConfiguration(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute contact detection
    configuration if omitted works on all objects  
";

%feature("docstring") mecaMAILx_NullifyReac "

set to 0 the reac of the IdBody mecaMAILx  

python usage : mecaMAILx_NullifyReac(datatype, IdBody)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `IdBody` :  
    (integer) : id of the concerned body  
";

%feature("docstring") mecaMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = mecaMAILx_GetAll(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") mecaMAILx_GetCooref "

return node coordinates of idBody  

python usage : array = mecaMAILx_GetCooref(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") mecaMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = mecaMAILx_GetConnectivity(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") mecaMAILx_GetElementsVolume "

return volume of elements  

python usage : volumes = mecaMAILx_GetElementsVolume(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
volumes[nb_ele] (double) : volume  
";

%feature("docstring") mecaMAILx_GetGpCoor "

return Gauss points coordinates of idBody  

python usage : array = mecaMAILx_GetGpCoor(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates of all Gauss points  
";

%feature("docstring") mecaMAILx_GetGpStrain "

return strain values stored at a gp  

python usage : strain = mecaMAILx_GetGpStrain(idBody,idEle,idGp)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  
* `IdGp` :  
    (integer) : id of the concerned gauss point  

Returns
-------
strain[size] (double) : value of strain  
";

%feature("docstring") mecaMAILx_GetGpStress "

return stress values stored at a gp  

python usage : stress = mecaMAILx_GetGpStress(idBody,idEle,idGp)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  
* `IdGp` :  
    (integer) : id of the concerned gauss point  

Returns
-------
stress[size] (double) : value of stress  
";

%feature("docstring") mecaMAILx_GetGpInternals "

return internal values stored at a gp  

python usage : internals = mecaMAILx_GetGpInternals(idBody,idEle,idGp)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  
* `IdGp` :  
    (integer) : id of the concerned gauss point  

Returns
-------
internals[nb_internals] (double) : value of internals  
";

%feature("docstring") mecaMAILx_GetGpPrincipalField "

return principal field (strain or stress) at a gp  

python usage : field = mecaMAILx_GetGpPrincipalField(idBody,idEle,idGp,idField)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  
* `IdGp` :  
    (integer) : id of the concerned gauss point  
* `IdField(integer)` :  
    : id of the field (1: strain, 2: stress)  

Returns
-------
field (double array): tensor field with principal values  
";

%feature("docstring") mecaMAILx_GetElementsInternal "

return a value over elements of an internal stored at gp  

python usage : internals = mecaMAILx_GetElementsInternal(idBody,id,f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `Id` :  
    (integer) : id of the internal  
* `f` :  
    (integer) : flag 1: mean, 2: sum, 3:max, 4: min  

Returns
-------
internals[nb_ele] (double) : value of internal  
";

%feature("docstring") mecaMAILx_GetElementsInternalIntegral "

return integral over elements of an internal stored at gp  

python usage : internals = mecaMAILx_GetElementsInternalIntegral(idBody,id)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `Id` :  
    (integer) : id of the internal  

Returns
-------
internals[nb_ele] (double) : value of internal  
";

%feature("docstring") mecaMAILx_GetElementsCenter "

return center of elements  

python usage : centers = mecaMAILx_GetElementsCenter(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
centers[3*nb_ele] (double) : center  
";

%feature("docstring") mecaMAILx_GetElementsJacobian "

return jacobian of elements  

python usage : jacobians = mecaMAILx_GetElementsJacobian(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
jacobians[nb_ele] (double) : jacobian  
";

%feature("docstring") mecaMAILx_ComputeElementsEnergy "

return energy of elements  

python usage : mecaMAILx_ComputeElementsEnergy()  
";

%feature("docstring") mecaMAILx_GetPtrElementsEnergy "

return energy of elements  

python usage : energies = mecaMAILx_GetElementsEnergy(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
energies (double array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_GetElementsNeighbor "

return elements in the tol-neighbor of an element of idBody  

python usage : neighbors = mecaMAILx_GetElementsNeighbor(idBody,tol)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `tol` :  
    (double) : tolerance  

Returns
-------
array (double 2D-array) : neighbor[nb_ele,max_neighbors]  
";

%feature("docstring") mecaMAILx_GetPtrElementsVisibility "

Get a pointer on the elements visibility vector.  

python usage : eviz = mecaMAILx_GetPtrElementsVisibility(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the mecaMAILx  

Returns
-------
eviz (int array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_AddNodalFieldDivergence "

Add the divergence of a diagonal field to external forces.  

python usage : mecaMAILx_AddNodalFieldDivergence(ibdyty, ifield)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `ifield` :  
    (integer) : rank of field  
";

%feature("docstring") mecaMAILx_CleanMemory "

Free all memory allocated within mecaMAILx module.  

python usage : mecaMAILx_CleanMemory()  
";

%feature("docstring") mecaMAILx_ComputeInfoPrincipalStressField "

Get info on the principal stress field: min,mean,max.  

Python usage : info = mecaMAILx_ComputeInfoPrincipalStressField(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
info (double array) : the desired info  
";

%feature("docstring") mecaMAILx_ComputePDFPressure "

Get pdf on the pressure.  

Python usage : pdf = mecaMAILx_ComputePDFPressure()  

Returns
-------
pdf (double array) : the desired info  
";

%feature("docstring") mecaMAILx_GetDeformationEnergy "

Get the deformation energy of a given displacement field.  

python usage : energy = mecaMAILx_GetDeformationEnergy(id,displacement)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `displacement` :  
    (double vector) : displacement field  

Returns
-------
energy (double) : deformation energy  
";

%feature("docstring") mecaMAILx_GetKineticEnergy "

Get the kinetic energy of a given velocity field.  

python usage : energy = mecaMAILx_GetKineticEnergy(id,velocity)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `velocity` :  
    (double vector) : velocity field  

Returns
-------
energy (double) : kinetic energy  
";

%feature("docstring") mecaMAILx_GetNeighborElementsToElement "

return neighbor elements to element idEle of body idBody  

python usage : vector = mecaMAILx_GetNeighborElementsToElement(idBody,idEle)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  

Returns
-------
vector (integer) : list of elements  
";

%feature("docstring") mecaMAILx_GetNeighborElementsToNode "

return neighbor elements to node idNode of body idBody  

python usage : vector = mecaMAILx_GetNeighborElementsToNode(idBody,idNode)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdNode` :  
    (integer) : id of the concerned node  

Returns
-------
vector (integer) : list of elements  
";

%feature("docstring") mecaMAILx_GetBoundaryElements "

return boundary elements  

python usage : vector = mecaMAILx_GetBoundaryElements(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : for each element =0 no boundary, otherwise gives the number
of free edge/face  
";

%feature("docstring") mecaMAILx_LoadWPreconBody "

load the precomputed W matrix on support node dofs of contactors for one body.
Assumes bulk behaviour is linear.  

python usage : mecaMAILx_LoadWPreconBody(ivalue)  

Parameters
----------
* `ivalue` :  
    (integer) : id of body to set precon  
";

%feature("docstring") mecaMAILx_GetPtrPreconW "

Get a pointer on the preconW Matrix of a given body.  

python usage : pcW = mecaMAILx_GetPtrPreconW(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the mecaMAILx  

Returns
-------
pcW (double array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_GetInternalVariable "

Get a copy of the internal variable of a given body.  

Python usage : internal = mecaMAILx_GetInternalVariable(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
internal (double array) : internal variable of desired body  
";

%feature("docstring") mecaMAILx_GetNbInternal "

Get the number of internal variable of a given body.  

python usage : nb_internal = mecaMAILx_GetNbInternal(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : rank of the body  

Returns
-------
nb_internal (integer) : number of internal variable of a body  
";

%feature("docstring") mecaMAILx_GetPtrBodyVector "

return pointer on body vector cvalue1_c of body IdBody  

Reac and Raux are impulsions (and not forces)  

python usage : vector_ptr = mecaMAILx_GetPtrBodyVector( cvalue1_c, IdBody )  

Parameters
----------
* `cvalue1_c` :  
    (string of size 5) : name of the body vector  
* `IdBody` :  
    (integer) : id of the body  

Returns
-------
vector_ptr (double array) : reference on the desired body vector  
";

%feature("docstring") mecaMAILx_GetDofStatus "

Get the status of nodes: 0 free, 1 x, 10 y.  

Python usage : vector = mecaMAILx_GetDofStatus(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired data  
";

%feature("docstring") mecaMAILx_PrepGlobalSolver "

computes free velocity of a list of bodies  

python usage : mecaMAILx_PrepGlobalSolver(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute free velocity if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_PostGlobalSolver "

computes the current d.o.f knowing all the forces (free + contact) of a list of
bodies  

python usage : mecaMAILx_PostGlobalSolver(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute current d.o.f if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_AddBodyForceToFext "

Add a body force (M*gamma) to Fext for a given body.  

python usage : mecaMAILx_AddBodyForceToFext(ibdyty, matrix)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `matrix` :  
    (double array) : the new value  
";

%feature("docstring") mecaMAILx_CheckProperties "

check if model and material are matching ; set material parameter if external
model  

python usage : mecaMAILx_CheckProperties()  
";

%feature("docstring") mecaMAILx_GetNbGpByElem "

Get the list of finite elements for MECAx models and the associated number of
Gauss Points.  

Here memory is allocated within lmgc90 so that the pointer can be freely
modified by third parties without nasty effect on lmgc90 functioning.  

python usage : names, nb_gps = mecaMAILx_GetNbGpByElem()  

Returns
-------
names (string list) : list of the finite elements nb_gps (integer list): list of
the number of Gauss Points  
";

%feature("docstring") mecaMAILx_MassScaling "

set mass scaling (default 1.d0)  

python usage : mecaMAILx_MassScaling(scale)  

Parameters
----------
* `scale` :  
    (double) : scaling  
";

%feature("docstring") mecaMAILx_GetGpAllJoint "

return GP value for joints  

python usage : vec = mecaMAILx_GetGpAllJoint()  

Returns
-------
vec (float matrix) : value at GP  
";

%feature("docstring") mecaMAILx_SetVisibleVlocyDrivenDof "

allows to (re)activate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : mecaMAILx_SetVisibleVlocyDrivenDof(ibdyty, inod, idof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the mecaMAILx  
* `inod(integer)` :  
    : index of the node to set visible  
* `idof(integer)` :  
    : index of the dof of the node to set visible  
";

%feature("docstring") mecaMAILx_SetInvisibleVlocyDrivenDof "

allows to deactivate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : mecaMAILx_SetInvisibleVlocyDrivenDof(ibdyty, inod, idof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the mecaMAILx  
* `inod(integer)` :  
    : index of the node to set invisible  
* `idof(integer)` :  
    : index of the dof of the node to set invisible  
";

%feature("docstring") mecaMAILx_UpdateVlocyDrivenDofStructures "

takes into account modifications on Vlocy driven dof status  

python usage : mecaMAILx_UpdateVlocyDrivenDofStructures(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the mecaMAILx  
";

