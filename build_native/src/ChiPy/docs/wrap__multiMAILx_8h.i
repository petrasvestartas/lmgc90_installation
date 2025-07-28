
// File: wrap__multiMAILx_8h.xml

%feature("docstring") multiMAILx_UsePicardScheme "

use Picard scheme (fixed point method)  

python usage : multiMAILx_UsePicardScheme()  
";

%feature("docstring") multiMAILx_UseNewtonScheme "

use Newton scheme  

python usage : multiMAILx_UseNewtonScheme()  
";

%feature("docstring") multiMAILx_WithoutRenumbering "

skip renumbering of the unknowns using a rcc method  

python usage : multiMAILx_WithoutRenumbering()  
";

%feature("docstring") multiMAILx_BandStorage "

use band matrix  

python usage : multiMAILx_BandStorage()  
";

%feature("docstring") multiMAILx_SparseStorage "

use sparse matrix  

python usage : multiMAILx_SparseStorage()  
";

%feature("docstring") multiMAILx_ExplodedStorage "

use element by element matrix  

python usage : multiMAILx_ExplodedStorage()  
";

%feature("docstring") multiMAILx_DiagonalStorage "

use diagonal matrix  

python usage : multiMAILx_DiagonalStorage()  
";

%feature("docstring") multiMAILx_SkylineStorage "

use skyline matrix  

python usage : multiMAILx_SkylineStorage()  
";

%feature("docstring") multiMAILx_FullStorage "

use full matrix  

python usage : multiMAILx_FullStorage()  
";

%feature("docstring") multiMAILx_SymmetricShape "

assume matrix is symmetrical  

python usage : multiMAILx_SymmetricShape()  
";

%feature("docstring") multiMAILx_UnspecifiedShape "

does not assume any thing on matrix shape  

python usage : multiMAILx_UnspecifiedShape()  
";

%feature("docstring") multiMAILx_GetNb "

Get the number of multiMAILx.  

python usage : nb_multiMAILx = multiMAILx_GetNb()  

Returns
-------
nb_multiMAILx (integer) : number of multiMAILx  
";

%feature("docstring") multiMAILx_GetNbNodes "

Get the number of nodes of a multiMAILx.  

python usage : nb_nodes = multiMAILx_GetNbNodes(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the multiMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a multiMAILx  
";

%feature("docstring") multiMAILx_GetNbElements "

Get the number of elements of a multiMAILx.  

python usage : nb_elements = multiMAILx_GetNbElements(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the multiMAILx  

Returns
-------
nb_nodes (integer) : number of elements of a multiMAILx  
";

%feature("docstring") multiMAILx_IsVisible "

return if a given body visible  

python usage : visible = multiMAILx_IsVisible(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") multiMAILx_GetBodyVector "

Get a copy of a vector of a given body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Pcbeg\": pressure of 1st fluid at beginning of time step  
*   \"Pc___\": pressure of 1st fluid in computed configuration  
*   \"Pnbeg\": pressure of 2nd fluid at beginning of time step  
*   \"Pn___\": pressure of 2nd fluid in computed configuration  
*   \"U_Fex\": external forces  
*   \"PcFex\": external pressure for 1st fluid  
*   \"PnFex\": external pressure for 2nd fluid  
*   \"U_Fin\": internal forces  
*   \"PcFin\": internal pressure for 1st fluid  
*   \"PnFin\": internal pressure for 2nd fluid  
*   \"U_Fdp\":  
*   \"PcFdp\":  
*   \"PnFdp\":  
*   \"U_Fdy\":  
*   \"PcFdy\":  
*   \"PnFdy\":  

Possible values for datatype field are \"X____\", \"Xbeg_\", \"V____\",
\"Vbeg_\", \"Pw___\", \"Pwbeg\", \"Pn___\", \"Pnbeg\"  

Python usage : vector = multiMAILx_GetBodyVector(datatype, ibdyty)  

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

%feature("docstring") multiMAILx_PutBodyVector "

Set a vector of a given body.  

Possible values for datatype field are:  

*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Pcbeg\": pressure of 1st fluid at beginning of time step  
*   \"Pc___\": pressure of 1st fluid in computed configuration  
*   \"Pnbeg\": pressure of 2nd fluid at beginning of time step  
*   \"Pn___\": pressure of 2nd fluid in computed configuration  

python usage : multiMAILx_PutBodyVector(datatype, ibdyty, matrix)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `ibdyty` :  
    (integer) : rank of body  
* `matrix` :  
    (double array) : the new values  
";

%feature("docstring") multiMAILx_ReadDrivenDof "

Read DRV_DOF.DAT.  

python usage : multiMAILx_ReadDrivenDof()  
";

%feature("docstring") multiMAILx_WriteDrivenDof "

Write DRV_DOF.OUT.  

python usage : multiMAILx_WriteDrivenDof()  
";

%feature("docstring") multiMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : multiMAILx_ReadIniGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : which GPV file to read  
";

%feature("docstring") multiMAILx_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

python usage : multiMAILx_ReadIniDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : which DOF file to read  
";

%feature("docstring") multiMAILx_WriteLastDof "

Write DOF.LAST file.  

python usage : multiMAILx_WriteLastDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write dof if omitted works on all
    objects  
";

%feature("docstring") multiMAILx_WriteOutDof "

Write DOF.OUT file.  

python usage : multiMAILx_WriteOutDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to write dof if omitted works on all
    objects  
";

%feature("docstring") multiMAILx_LoadBehaviours "

load behaviours from bulk_behav  

python usage : multiMAILx_LoadBehaviours()  
";

%feature("docstring") multiMAILx_LoadModels "

load models from models  

python usage : multiMAILx_LoadModels()  
";

%feature("docstring") multiMAILx_PushProperties "

gives to model the couple of model,behavior used at gauss point  

python usage : multiMAILx_PushProperties()  
";

%feature("docstring") multiMAILx_IncrementStep "

initializes the current d.o.f and some driven d.o.f values  

python usage : multiMAILx_IncrementStep()  
";

%feature("docstring") multiMAILx_ComputeMass "

compute elementary mass and inertia of a list of bodies  

python usage : multiMAILx_ComputeMass(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute mass and inertia if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_ComputeBulk "

computes elementary stiffness and viscosity matrices of a list of bodies  

python usage : multiMAILx_ComputeBulk(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute stiffness and viscosity
    matrices if omitted works on all objects  
";

%feature("docstring") multiMAILx_ComputeFext "

compute elementary external forces of a list of bodies  

python usage : multiMAILx_ComputeFext(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute external forces if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_AssembKT "

assemble pseudo mass matrix and apply drvdof of a list of bodies  

python usage : multiMAILx_AssembKT(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to assemble pseudo mass matrix and apply
    drvdof if omitted works on all objects  
";

%feature("docstring") multiMAILx_AssembRHS "

assembles right hand side of a list of bodies  

python usage : multiMAILx_AssembRHS(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to assemble right hand side if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_ComputeResidueNorm "

computes the norm of the residue of a list of bodies  

python usage : norm = multiMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute the norm of the residue if
    omitted works on all objects  

Returns
-------
norm (double) : Residue Norm  
";

%feature("docstring") multiMAILx_ComputeFreeState "

computes free (of interactions) state of a list of bodies  

python usage : multiMAILx_ComputeFreeState(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute free state if omitted works on
    all objects  
";

%feature("docstring") multiMAILx_ComputeDof "

computes the current d.o.f knowing all the forces/fluxses (free + contact) of a
list of bodies  

python usage : multiMAILx_ComputeDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute current d.o.f if omitted works
    on all objects  
";

%feature("docstring") multiMAILx_ComputeField "

computes elementary fields of a list of bodies  

python usage : multiMAILx_ComputeField(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_UpdateBulk "

update begin elementary fields with current elementary fields of a list of
bodies  

python usage : multiMAILx_UpdateBulk(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_UpdateDof "

update begin d.o.f. with current d.o.f. of a list of bodies  

python usage : multiMAILx_UpdateDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to update current d.o.f if omitted works
    on all objects  
";

%feature("docstring") multiMAILx_GetScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = multiMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

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

%feature("docstring") multiMAILx_SetScalarFieldByNode "

Update elementary fields through a nodal external field on a given body.  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetScalarFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") multiMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") multiMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = multiMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

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

%feature("docstring") multiMAILx_SetVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") multiMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") multiMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = multiMAILx_GetConnectivity(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") multiMAILx_GetCoor "

return node coordinates of idBody  

python usage : array = multiMAILx_GetCoor(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") multiMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = multiMAILx_GetAll(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") multiMAILx_GetElementsVolume "

return volume of elements  

python usage : volumes = multiMAILx_GetElementsVolume(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
volumes[nb_ele] (double) : volume  
";

%feature("docstring") multiMAILx_GetElementsNeighbor "

return elements in the tol-neighbor of an element of idBody  

python usage : neighbors =
multiMAILx_GetElementsNeighbor(idBody,tol,max_neighbors)  

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

%feature("docstring") multiMAILx_GetPtrElementsEnergy "

return pointer on energy of elements  

python usage : energies = multiMAILx_GetPtrElementsEnergy(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
energies[nb_ele] (double) : energy  
";

%feature("docstring") multiMAILx_ComputeElementsEnergy "

compute energy of elements  

python usage : multiMAILx_ComputeElementsEnergy(idBody)  
";

%feature("docstring") multiMAILx_GetPtrElementsJacobian "

return jacobian of elements  

python usage : jacobians = multiMAILx_GetPtrElementsJacobian(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
jacobians[nb_ele] (double) : jacobian  
";

%feature("docstring") multiMAILx_ComputeElementsJacobian "

compute jacobian of elements  

python usage : multiMAILx_ComputeElementsJacobian(idBody)  
";

%feature("docstring") multiMAILx_GetPtrElementsVisibility "

Get a pointer on the elements visibility vector.  

python usage : eviz = multiMAILx_GetPtrElementsVisibility(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the multiMAILx  

Returns
-------
eviz (int array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") multiMAILx_GetDeformationEnergy "

Get the deformation energy of a given displacement field.  

python usage : energy = multiMAILx_GetDeformationEnergy(id,displacement)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `displacement` :  
    (double matrix) : displacement field  

Returns
-------
energy (double) : deformation energy  
";

%feature("docstring") multiMAILx_GetPtrBoundaryElements "

return boundary elements  

python usage : vector = multiMAILx_GetPtrBoundaryElements(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : for each element =0 no boundary, otherwise gives the number
of free edge/face  
";

%feature("docstring") multiMAILx_CleanMemory "

Free all memory allocated within multiMAILx module.  

python usage : multiMAILx_CleanMemory()  
";

