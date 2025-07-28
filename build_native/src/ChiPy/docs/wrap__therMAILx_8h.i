
// File: wrap__therMAILx_8h.xml

%feature("docstring") therMAILx_GetNbTherMAILx "

Get the number of therMAILx.  

python usage : nb_therMAILx = therMAILx_GetNbTherMAILx()  

Returns
-------
nb_therMAILx (integer) : number of therMAILx  
";

%feature("docstring") therMAILx_GetNbNodes "

Get the number of nodes of a therMAILx.  

python usage : nb_nodes = therMAILx_GetNbNodes(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the therMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a therMAILx  
";

%feature("docstring") therMAILx_GetNbElements "

Get the number of nodes of a therMAILx.  

python usage : nb_nodes = therMAILx_GetNbElements(ibdyty)  

Parameters
----------
* `ivalue` :  
    (integer) : id of the therMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a therMAILx  
";

%feature("docstring") therMAILx_GetNbDofs "

Get the number of dofs for the therMAILX.  

python usage : nb_dofs = therMAILx_GetNbDofs(int ibdyty)  

Returns
-------
nb_dofs (integer) : number of dofs of the body for the model  
";

%feature("docstring") therMAILx_IncrementStep "

initializes current dof  

python usage : therMAILx_IncrementStep()  
";

%feature("docstring") therMAILx_ComputeConductivity "

computes the elementary conductivity matrices of a list of bodies  

If the input list is empty, the conductivities of all bodies will be computed  

python usage : therMAILx_ComputeConductivity(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeCapacity "

computes the elemetary capacity matrices  

python usage : therMAILx_ComputeCapacity(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeConvection "

compute elementary convection terms  

python usage : therMAILx_ComputeConvection(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeInternalFlux "

compute elementary internal flux  

python usage : therMAILx_ComputeInternalFlux(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeExternalFlux "

compute elementary external flux  

python usage : therMAILx_ComputeExternalFlux(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_AssembThermKT "

assembles elementary matrices  

python usage : therMAILx_AssembKT(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_AssembThermRHS "

assembles elementary vectors  

python usage : therMAILx_AssembRHS(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeThermDof "

computes current dof  

python usage : therMAILx_ComputeThermDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeThermFields "

computes elementary fields  

python usage : therMAILx_ComputeThermFields(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_UpdateThermDof "

update begin dof with current dof  

python usage : therMAILx_UpdateThermDof(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_UpdateThermBulk "

update begin elementary fields with current elementary fields  

python usage : therMAILx_UpdateThermBulk(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeResidueNorm "

compute the residue of the thermal equation  

python usage : norm = therMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  

Returns
-------
norm (double) : value of the norm  
";

%feature("docstring") therMAILx_ReadDrivenDof "

Read DRV_DOF.DAT.  

python usage : therMAILx_ReadDrivenDof()  
";

%feature("docstring") therMAILx_WriteDrivenDof "

Write DRV_DOF.OUT.  

python usage : therMAILx_WriteDrivenDof()  
";

%feature("docstring") therMAILx_LoadModels "

loads models frol models module  

python usage : therMAILx_LoadModels()  
";

%feature("docstring") therMAILx_LoadBehaviours "

loads bulk behaviors parameters from bulk_behav module  

python usage : therMAILx_LoadBehaviours()  
";

%feature("docstring") therMAILx_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read Else : OUTBOX/DOF.OUT.num is read, num
being the parameter used in TimeEvolution_ReadIniDof last call  

python usage : therMAILx_ReadIniDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : which DOF file to read  
";

%feature("docstring") therMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : therMAILx_ReadIniGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : which GPV file to read  
";

%feature("docstring") therMAILx_WriteLastDof "

Write ascii DOF.LAST file.  

python usage : therMAILx_WriteLastDof()  
";

%feature("docstring") therMAILx_WriteOutDof "

Write ascii DOF.OUT file. Can be activate only each N step.  

python usage : therMAILx_WriteOutDof()  
";

%feature("docstring") therMAILx_DisplayOutDof "

Display body degrees of freedom.  

python usage : therMAILx_DisplayOutDof()  
";

%feature("docstring") therMAILx_PutBodyVector "

Set a vector of a given body.  

Possible values for datatype field are:  

*   \"T____\": Temperature in computed configuration  
*   \"Tbeg_\": Temperature at beginning of time step  
*   \"Taux_\": Temperature in working array  
*   \"Fext_\": external flux  
*   \"Fint_\": internal flux  

Uses copy  

python usage : therMAILx_PutBodyVector(datatype, ibdyty, matrix)  

Parameters
----------
* `datatype` :  
    (string of size 5) : the vector to set  
* `ibdyty` :  
    (integer) : rank of body  
* `matrix` :  
    (double array) : the new values  
";

%feature("docstring") therMAILx_GetBodyVector "

Get a copy of a vector of a given body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"T____\": Temperature in computed configuration  
*   \"Tbeg_\": Temperature at beginning of time step  
*   \"Taux_\": Temperature in working array  
*   \"Fext_\": external flux  
*   \"Fint_\": internal flux  

Uses copy  

python usage : vector = therMAILx_GetBodyVector(datatype, ibdyty)  

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

%feature("docstring") therMAILx_GetScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = therMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

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

%feature("docstring") therMAILx_SetScalarFieldByNode "

Update an external field on a given body.  

You need to set this field in your models.dat  

python usage : therMAILx_SetScalarFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") therMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : therMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the field to set  
* `f` :  
    (double array) : value of the field  
";

%feature("docstring") therMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = therMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

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

%feature("docstring") therMAILx_SetVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : therMAILx_SetFieldByNode(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") therMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : therMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concern body  
* `f_rank` :  
    (integer) : rank of the vector field to set  
* `f` :  
    (double array) : value of the vector field  
";

%feature("docstring") therMAILx_AddSource "

Add a volumic source into a given body.  

python usage : therMAILx_AddSource(ibdyty, ifield)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `ifield` :  
    (integer) : rank of field  
";

%feature("docstring") therMAILx_AddNodalFieldDivergence "

Add the divergence of a field to external flux.  

python usage : therMAILx_AddNodalFieldDivergence(ibdyty, ifield)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `ifield` :  
    (integer) : rank of field  
";

%feature("docstring") therMAILx_PushProperties "

declares to module model the couples (model,behavior) used  

python usage : therMAILx_PushProperties()  
";

%feature("docstring") therMAILx_WithoutRenumbering "

skip renumbering of the unknowns using a rcc method  

python usage : therMAILx_WithoutRenumbering()  
";

%feature("docstring") therMAILx_BandStorage "

use band matrix  

python usage : therMAILx_BandStorage()  
";

%feature("docstring") therMAILx_SparseStorage "

use sparse matrix  

python usage : therMAILx_SparseStorage()  
";

%feature("docstring") therMAILx_ExplodedStorage "

use element by element matrix  

python usage : therMAILx_ExplodedStorage()  
";

%feature("docstring") therMAILx_DiagonalStorage "

use diagonal matrix  

python usage : therMAILx_DiagonalStorage()  
";

%feature("docstring") therMAILx_SkylineStorage "

use skyline matrix  

python usage : therMAILx_SkylineStorage()  
";

%feature("docstring") therMAILx_FullStorage "

use full matrix  

python usage : therMAILx_FullStorage()  
";

%feature("docstring") therMAILx_SymmetricShape "

assume matrix is symmetrical  

python usage : therMAILx_SymmetricShape()  
";

%feature("docstring") therMAILx_UnspecifiedShape "

does not assume any thing on matrix shape  

python usage : therMAILx_UnspecifiedShape()  
";

%feature("docstring") therMAILx_GetGrad "

Get a copy of a gradient of a given body.  

Python usage : grad_T = therMAILx_GetGrad(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
grad_T (double 2D-array) : the desired gradient  
";

%feature("docstring") therMAILx_GetFlux "

Get a copy of a gradient of a given body.  

Python usage : Flux_T = therMAILx_GetFlux(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  

Returns
-------
Flux_T (double array) : the desired flux  
";

%feature("docstring") therMAILx_InitializeElementaryFlux "

set elementary flux to 0  

python usage : therMAILx_InitializeElementaryFlux()  
";

%feature("docstring") therMAILx_GetCoor "

return node coordinates of idBody  

python usage : array = therMAILx_GetCoor(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") therMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = therMAILx_GetConnectivity(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") therMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = therMAILx_GetAll(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") therMAILx_GetGpCoor "

return Gauss points coordinates of idBody  

python usage : array = therMAILx_GetGpCoor(idBody)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates of all Gauss points  
";

%feature("docstring") therMAILx_GetGpField "

return field values stored at a gp  

python usage : field = therMAILx_GetGpField(idBody,idEle,idGp,idField)  

Parameters
----------
* `IdBody` :  
    (integer) : id of the concerned body  
* `IdEle` :  
    (integer) : id of the concerned element  
* `IdGp` :  
    (integer) : id of the concerned gauss point  
* `IdField` :  
    (integer) : id of the concerned field  

Returns
-------
field (double array) : field value  
";

%feature("docstring") therMAILx_TrialAssembThermKT "

[experimental] assembles elementary matrices  

python usage : therMAILx_AssembThermKT()  
";

%feature("docstring") therMAILx_TrialAssembThermRHS "

[experimental] assembles elementary vectors  

python usage : therMAILx_AssembThermRHS()  
";

%feature("docstring") therMAILx_CleanMemory "

Free all memory allocated within therMAILx module.  

python usage : therMAILx_CleanMemory()  
";

%feature("docstring") therMAILx_CheckProperties "

check if model and material are matching ; set material parameter if external
model  

python usage : therMAILx_CheckProperties()  
";

%feature("docstring") therMAILx_GetNbGpByElem "

Get the list of finite elements for therx models and the associated number of
Gauss Points.  

Here memory is allocated within lmgc90 so that the pointer can be freely
modified by third parties without nasty effect on lmgc90 functioning.  

python usage : names, nb_gps = therMAILx_GetNbGpByElem()  

Returns
-------
names (string list) : list of the finite elements nb_gps (integer list): list of
the number of Gauss Points  
";

%feature("docstring") therMAILx_GetNbGp "

Get the number of Gauss points of an element of a therMAILx.  

python usage : nb_gp = therMAILx_GetNbElements(ibdyty, iblmty)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of the therMAILx  
* `iblmty` :  
    (integer) : id of the element  

Returns
-------
nb_gp (integer) : number of Gauss point of an element of a therMAILx  
";

