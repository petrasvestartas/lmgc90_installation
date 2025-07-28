
// File: wrap__models_8h.xml

%feature("docstring") models_ReadModels "

read models from DATBOX/MODELS.DAT  

python usage : models_ReadModels()  
";

%feature("docstring") models_WriteModels "

write models to OUTBOX/MODELS.OUT  

python usage : models_WriteModels()  
";

%feature("docstring") models_InitModels "

initialize models  

python usage : models_InitModels()  
";

%feature("docstring") models_InitProperties "

initialize properties  

In face re-initialize properties (since it is done in InitModels). Necessary if
a Store has been done and it is wanted again to LoadModel  

python usage : models_InitProperties()  
";

%feature("docstring") models_StoreProperties "

create properties (couple of model and models)  

python usage : models_StoreProperties()  
";

%feature("docstring") models_CleanMemory "

Free all memory allocated within models module.  

python usage : models_CleanMemory()  
";


// File: wrap__mp__solver_8h.xml

%feature("docstring") mp_solver_ReadMpBehaviour "

python usage : mp_solver_ReadMpBehaviour()  
";

%feature("docstring") mp_solver_WriteMpBehaviour "

python usage : mp_solver_WriteMpBehaviour()  
";

%feature("docstring") mp_solver_ReadIniMpValues "

Read MP_VALUES file.  

If num <= 0 : DATBOX/MP_VALUES.INI file is read Else : OUTBOX/MP_VALUES.OUT.num
is read, num being the parameter used in TimeEvolution_ReadIniDof last call  

usage : mp_solver_ReadIniMpValues(num=0)  

Parameters
----------
num(integer) : which file to read  
";

%feature("docstring") mp_solver_WriteOutMpValues "

python usage : mp_solver_WriteOutMpValues()  
";

%feature("docstring") mp_solver_WriteLastMpValues "

python usage : mp_solver_WriteLastMpValues()  
";

%feature("docstring") mp_solver_SolveElectro1G "

python usage : mp_solver_SolveElectro1G()  
";

%feature("docstring") mp_solver_SolveNlElectro1G "

python usage : mp_solver_SolveNlElectro1G()  
";

%feature("docstring") mp_solver_SolveThermoProblem "

python usage : mp_solver_SolveThermoProblem()  
";

%feature("docstring") mp_solver_UpdateThermoProblem "

python usage : mp_solver_UpdateThermoProblem()  
";

%feature("docstring") mp_solver_RecupTemperature "

python usage : mp_solver_RecupTemperature()  
";

%feature("docstring") mp_solver_RecupPotential "

python usage : mp_solver_RecupPotential()  
";

%feature("docstring") mp_solver_UpdateConductivity "

python usage : mp_solver_UpdateConductivity()  
";

%feature("docstring") mp_solver_InitThermalConductivity "

python usage : mp_solver_InitThermalConductivity()  
";

%feature("docstring") mp_solver_GetBrancheValues "
";

%feature("docstring") mp_solver_PutHeatGenerationFactor "

python usage : value = mp_solver_PutHeatGenerationFactor(ivalue)  
";

%feature("docstring") mp_solver_PutHeatConductionContinueFactor "

python usage : value = mp_solver_PutHeatConductionContinueFactor(ivalue)  
";


// File: wrap__ALpxx_8h.xml

%feature("docstring") ALpxx_LoadTactors "

load ALpxx from MAILx and initialize existing_entities  

python usage : ALpxx_LoadTactors()  
";

%feature("docstring") ALpxx_PushPreconNodes "

set ALpxx supporting nodes as precon  

python usage : ALpxx_PushPreconNodes()  
";

%feature("docstring") ALpxx_GetAllConnec "

return connectivity of all AL in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = ALxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of ALxxx elements  
";

%feature("docstring") ALpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all ALxxx  

python usage : idata, rdata = ALxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") ALpxx_CleanMemory "

Free all memory allocated within ALpxx module.  

python usage : ALpxx_CleanMemory()  
";


// File: wrap__xKSID_8h.xml

%feature("docstring") xKSID_LoadTactors "

load xKSID from RBDY2 and initialize existing_entites  

python usage : xKSID_LoadTactors()  
";

%feature("docstring") xKSID_GetNbxKSID "

Get the number of xKSID in the container.  

python usage : nb_diskx = xKSID_GetNbxKSID()  

Returns
-------
nb_xKSID (integer) : the number of xKSID in container  
";

%feature("docstring") xKSID_GetPtrxKSID2BDYTY "

return a pointer onto the map xksid2rbdy2  

python usage : xksid2rbdy2 = xKSID_GetPtrxKSID2BDYTY()  

Returns
-------
xksid2rbdy2 (integer array) : reference on map between xksid rank and body/tact
rank  
";

%feature("docstring") xKSID_IsVisible "

return if a body visible  

usage : visible = xKSID_IsVisible(itact)  

Parameters
----------
itact(integer) : rank of xKSID  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") xKSID_GetContactorRadius "

Get the radius of a given xKSID.  

python usage : radius = xKSID_GetContactorRadius(itact)  

Parameters
----------
itact(integer) : rank of a xKSID (in the list of all the xKSID)  

Returns
-------
radius (double) : the radius of the xKSID of rank itact  
";

%feature("docstring") xKSID_GetContactorCoor "

get coordinates of the center of a given xKSID  

usage : vector = xKSID_GetContactorCoor(itacty)  

Parameters
----------
itacty(integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") xKSID_InitOutlines "

Get a reference on the outlines of all xKSID.  

usage : outlines = xKSID_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_xKSID  
";

%feature("docstring") xKSID_InitScalarFields "

Get a reference on the scalar fields of all xKSID.  

usage : scalarfields = xKSID_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_xKSID array  
";

%feature("docstring") xKSID_UpdatePostdata "

Update values of outlines_xKSID and scalarfields_xKSID pointers.  

usage : xKSID_UpdatePostdata  
";

%feature("docstring") xKSID_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = xKSID_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
xKSID  
";

%feature("docstring") xKSID_GetNbScalarFields "

Get the number of scalar fields of a xKSID.  

python usage : nb_scalarfields = xKSID_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a xKSID  
";

%feature("docstring") xKSID_CleanMemory "

Free all memory allocated within xKSID module.  

python usage : xKSID_CleanMemory()  
";

%feature("docstring") xKSID_SetXdilation "

set increase of radius of a xKSID due to expansion  

python usage : xKSID_SetXdilation(itacty,x)  

Parameters
----------
itacty(integer) : rank of considered contactor  
x(float) : increase of radius  
";

%feature("docstring") xKSID_SetVdilation "

set increase rate of radius of a xKSID due to expansion  

python usage : xKSID_SetVdilation(itacty, v)  

Parameters
----------
itacty(integer) : rank of contactor  
v(float) : radius increase rate  
";


// File: wrap__a__EF_8h.xml

%feature("docstring") a_EF_InterpolateField "
";

%feature("docstring") a_EF_ComputeCenter "

Compute the geometric center of an element.  

python usage : center = a_EF_ComputeCenter(coor)  

Parameters
----------
coor(double array) : coordinates the nodes of the element  

Returns
-------
center (double array) : computed center of the element  
";


// File: wrap__postpro__3D_8h.xml

%feature("docstring") postpro_3D_PostproDuringComputation "

Scan postprocessing function which should be call during the computation
process.  

python usage : postpro_3D_PostproDuringComputation()  
";

%feature("docstring") postpro_3D_FlushDuringComputation "

Flush all postpro files.  

python usage : postpro_3D_FlushDuringComputation()  
";

%feature("docstring") postpro_3D_ReadCommands "

Scan postprocessing functions which should be call during the computation
process.  

python usage : postpro_3D_ReadCommands()  
";

%feature("docstring") postpro_3D_PostproBeforeComputation "

Data initialization.  

python usage : postpro_3D_PostproBeforeComputation(restart=False) param[in]
restart (integer) : if the Postpro file must append to existing ones and
starting index of CONTACT_FORCE_DISTRIBUTION files  
";

%feature("docstring") postpro_3D_ClosePostproFiles "

Close all postpro files.  

python usage : postpro_3D_ClosePostproFiles()  
";

%feature("docstring") postpro_3D_GetKineticEnergy "

Compute Kinetic Energy for all bodies (rigids and defo)  

python usage : KE = postpro_3D_GetKineticEnergy()  
";

%feature("docstring") postpro_3D_GetRBDY3PrincStress "

Return the principal stresses on each RBDY3.  

python usage : pstress = postpro_3D_GetRBDY3PrincStress()  

Returns
-------
pstress (double 2D-array) : the interactions  
";

%feature("docstring") postpro_3D_CleanMemory "

Free all memory allocated within postpro_3D module.  

python usage : postpro_3D_CleanMemory()  
";


// File: wrap__mbs2D_8h.xml

%feature("docstring") MBS2D_setNb "

Set the number of MBS.  

python usage : MBS2D_setNb(nb)  

Parameters
----------
nb(integer) : set the number of MBS  
";

%feature("docstring") MBS2D_getNb "

Get the number of MBS.  

python usage : nb = MBS2D_getNb()  

Returns
-------
nb (integer) : the number of MBS  
";

%feature("docstring") MBS2D_setNbNodes "

Set the number of nodes of a MBS.  

python usage : MBS2D_setNbNodes(ibdyty, nb)  

Parameters
----------
ibdyty(integer): id of the MBS  
nb(integer) : the number of nodes of the MBS  
";

%feature("docstring") MBS2D_setNbTactors "

Set the number contactors of a MBS.  

python usage : MBS2D_setNbTactors(ibdyty, nb)  

Parameters
----------
ibdyty(integer): id of the MBS  
nb(integer) : the number of contactor of the MBS  
";

%feature("docstring") MBS2D_getPtrCoor "

Get a pointer on the coor of a MBS.  

python usage : coor = MBS2D_getPtrCoor(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered MBS  

Returns
-------
coor (double 2D-array) : reference on the coordinates of the nodes  
";

%feature("docstring") MBS2D_getPtrCoorTT "

Set the array of coordinates of nodes of a MBS.  

python usage : coorTT = MBS2D_GetPtrCoorTT(ibdyty)  

Parameters
----------
ibdyty(integer): id of the MBS  

Returns
-------
coorTT (double array) : coordinates of nodes of a MBS  
";

%feature("docstring") MBS2D_addContactor "

Add a new contactor to a MBS.  

Available contactor types are :  

*   JONCx: inputs are:
    -   rdata must hold [axe_x, axe_y]  
*   POLYR: inputs are:
    -   rdata must hold the coordinates of the vertices [x_1, y_1, ... x_n, y_n]  
    -   idata must hold the number of vertices  

python usage : MBS2D_addContactor(ibdyty, inodty, itacty, tacttype, color,
rdata, idata=None)  

Parameters
----------
ibdyty(integer) : rank of the MBS  
inodty(integer) : rank of the node of the MBS the contactor is tied to  
itacty(integer) : rank of the contactor of MBS  
tactype(string [5]) : the type of contactor  
color(string [5]) : the color of the contactor  
rdata(double array) : the new value of the vector  
idata(integer array) : the new value of the vector  
";

%feature("docstring") MBS2D_initialize "

Initialize MBS module once loading is done.  

python usage : MBS2D_initialize()  
";

%feature("docstring") MBS2D_finalize "

Finalize MBS module.  

python usage : MBS2D_finalize()  
";

%feature("docstring") MBS2D_IncrementStep "

compute the current velocity and displacement  

python usage : MBS2D_IncrementStep()  
";

%feature("docstring") MBS2D_ComputeFreeVelocity "

compute free velocity  

python usage : MBS2D_ComputeFreeVelocity()  
";

%feature("docstring") MBS2D_ComputeDof "

update current position and velocity  

python usage : MBS2D_ComputeDof()  
";

%feature("docstring") MBS2D_UpdateDof "

save d.o.f. of the end of the time step to d.o.f. of the begining of the next
one  

python usage : MBS2D_UpdateDof()  
";


// File: wrap__DDM__2D_8h.xml

%feature("docstring") DDM_2D_SetDDWorkingDirectory "
";

%feature("docstring") DDM_2D_Initialize "

Initialize 2D DDM module.  

ddm_type may be :  

*   1 for Feti (DDM without overlap)  
*   2 for Schwarz (DDM with overlap)  

python usage : DDM_2D_Initialize(nb_sdmx, nb_sdmy, ddm_type)  

Parameters
----------
nb_sdmx(integer) : number of domains on x-axis  
nb_sdmy(integer) : number of domains on y-axis  
ddm_type(integer) : type of DDM to use  
";

%feature("docstring") DDM_2D_Partitioning "
";

%feature("docstring") DDM_2D_AddToFext "

Add external forces due to DDM.  

python usage : DDM_2D_AddToFext()  

To add after RBDY2_ComputeFext  
";

%feature("docstring") DDM_2D_ExSolver "

Solve fully the local contact problem with DDM.  

python usage : DDM_2D_ExSolver(storage, checktype, tol, relax, nb_iter_check,
nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";

%feature("docstring") DDM_2D_ComputeDof "

Compute degrees of freedom.  

python usage : DDM_2D_ComputeDof()  
";

%feature("docstring") DDM_2D_Post "

Does postpro operation related to DDM.  

python usage : DDM_2D_Post()  
";

%feature("docstring") DDM_2D_SetParameters "

Set frequencies parameters of DDM.  

python usage : DDM_2D_SetParameters(f_ddm, f_out, f_last, f_postpro, f_display)  

Parameters
----------
f_ddm(integer) : frequency of ddm partionning  
f_out(integer) : frequency of output file writing  
f_last(integer) : frequency of last file writing  
f_postpro(integer) : frequency of postpro file writing  
f_display(integer) : frequency of display file writing  
";

%feature("docstring") DDM_2D_WriteLast "

Write some data at the end of computation.  

python usage : DDM_2D_WriteLast()  
";

%feature("docstring") DDM_2D_Finalize "

End of computation/module's life.  

python usage : DDM_2D_Finalize()  
";


// File: wrap__cut2D_8h.xml

%feature("docstring") cut2D_Cut "
";


// File: wrap__deposit3D_8h.xml

%feature("docstring") deposit3D_InContainer "

Computes a new deposit under gravity in a container.  

i_shape = 0 : box  

*   a point (x, y, z) is in the box iff x is in [-lx/2, lx/2], y is in [-ly/2,
    ly/2] and z is in [0, lz] i_shape = 1 : cylinder  
*   a point (x, y, z) is in the cylinder iff x^2 + y^2 is in [0, R^2] and z is
    in [0, lz] i_shape = 2 : sphere  
*   a point (x, y, z) is in the sphere iff x^2 + y^2 + z^2 is in [0, R^2]  

python call: radii, coor = deposit3D_InContaier(in_radii, shape, p1, p2, p3[,
dradii, dcoor, seed, with_log])  

Parameters
----------
in_radii(double array): given radii list (i.e. granulometry)  
shape(integer)of container (0->box, 1->cylinder, 2->sphere)  
p1(double): box-> lx, cylinder->R, sphere->R  
p2(double): box-> ly, cylinder->lz, sphere->ignored  
p3(double): box-> lz, cylinder->ignored, sphere->ignored  
dradii(double array) (optional) : a list of already deposited radii  
dcoor(double array) (optional) : a list of already deposited coor (must be of
    size [nb_dradii,3])  
seed(integer array) (optional) : an input seed to control randomness  
with_log(integer)de/activate log message  

Returns
-------
radii (double array): list of deposited radii coor (double array): coordinates
of deposited radii (shape [nb_radii,3]) PYDOC  
";


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
ivalue(integer) : id of the mecaMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a mecaMAILx  
";

%feature("docstring") mecaMAILx_GetNbElements "

Get the number of elements of a mecaMAILx.  

python usage : nb_elements = mecaMAILx_GetNbElements(ibdyty)  

Parameters
----------
ivalue(integer) : id of the mecaMAILx  

Returns
-------
nb_nodes (integer) : number of elements of a mecaMAILx  
";

%feature("docstring") mecaMAILx_GetNbGp "

Get the number of Gauss points of an element of a mecaMAILx.  

python usage : nb_gp = mecaMAILx_GetNbElements(ibdyty, iblmty)  

Parameters
----------
ibdyty(integer) : id of the mecaMAILx  
iblmty(integer) : id of the element  

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
ivalue(integer) : id of body to set precon  
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
ivalue1(integer) : body number  
ivalue2(integer) : node number  
ivalue3(integer) : dof number  
vect(double) : column  
";

%feature("docstring") mecaMAILx_GetNodesPrecon "

Get the list of preconditionned nodes of a mecaMAILx body.  

Here memory is allocated within lmgc90 so that the pointer can be freely
modified by third parties without nasty effect on lmgc90 functioning.  

python usage : precon_list = mecaMAILx_GetNodesPrecon(ibdyty)  

Parameters
----------
ibdyty(integer) : index of the desired mecaMAILx  

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
ivalue(integer) : id of body to set coro  
";

%feature("docstring") mecaMAILx_SetTolCoro "

set the admssible tolerance on rigid body velocity computed by deformable model  

python usage : mecaMAILx_SetTolCoro(tol)  

Parameters
----------
tol(double) : tolerance  
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
ivalue(integer) : id of body to compute as a rigid  
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
ivalue(integer) : id of body to compute without deformation  
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
idbdy(integer): id of the body we want visibility  

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
idbdy(integer): id of the body  

Returns
-------
vec (float matrix) : frame matrix (beg, current or TT)  
";

%feature("docstring") mecaMAILx_GetRigidCoorTT "

return TT center of inertia coordinates  

python usage : vec = mecaMAILx_GetRigidCoorTT(ibdyty)  

Parameters
----------
idbdy(integer): id of the body  

Returns
-------
vec (float vector) : TT center of inertia coordinates  
";

%feature("docstring") mecaMAILx_GetRigidCooref "

return ref center of inertia coordinates  

python usage : vec = mecaMAILx_GetRigidCooref(ibdyty)  

Parameters
----------
idbdy(integer): id of the body  

Returns
-------
vec (float vector) : ref center of inertia coordinates  
";

%feature("docstring") mecaMAILx_SetRVDrivenDofs "

declares rigid velocity dof as driven  

python usage : mecaMAILx_SetRVDrivenDofs(idbody,vector_in)  

Parameters
----------
idbody(integer) : id of the body  
vector(integer) : list of driven dofs  
";

%feature("docstring") mecaMAILx_SetRVDrivenDofValue "

set the value of rigid velocity dof value  

python usage : mecaMAILx_SetRVDrivenDofValue(idbody,iddof,rv)  

Parameters
----------
idbody(integer) : id of the body  
iddof(integer) : id of dof  
rv(float) : value  
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
datatype(string [5]) : the vector to set  
ibdyty(integer) : rank of the RBDY3  
vector(double array) : the new value of the vector  
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
datatype(string [5]) : the vector to get  
ibdyty(integer) : rank of the RBDY3  

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
datatype(string of size 5) : the vector to set  
ibdyty(integer) : rank of body  
matrix(double array) : the new value  
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
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired data  
";

%feature("docstring") mecaMAILx_GetMaterials "

Get a copy of a the elements' material vector of a given body.  

Python usage : materials = mecaMAILx_GetMaterials(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

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
ibdyty(integer) : rank of considered body  

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
ibdyty(integer) : rank of considered body  

Returns
-------
strain (double 2D-array) : nodal strain of the desired body  
";

%feature("docstring") mecaMAILx_GetInternalVariables "

Get a copy of the smoothed nodal internal variables (2D:10 ; 3D:57)  

Python usage : strain = mecaMAILx_GetInternalVariables(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

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
ibdyty(integer) : rank of considered body  

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
i_list(list of integer) : list of bodies to compute free velocity if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_AssembKT "

assemble pseudo mass matrix and apply drvdof of a list of bodies  

python usage : mecaMAILx_AssembKT(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to assemble pseudo mass matrix and apply
    drvdof if omitted works on all objects  
";

%feature("docstring") mecaMAILx_OnlyAssembKT "

assemble pseudo mass matrix of a list of bodies  

python usage : mecaMAILx_OnlyAssembKT(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to assemble pseudo mass matrix if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ApplyDrvDofKT "

apply drvdof pseudo mass matrix  

python usage : mecaMAILx_ApplyDrvDofKT(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to apply drvdof on pseudo mass matrix if
    omitted works on all objects  
";

%feature("docstring") mecaMAILx_AssembRHS "

assembles right hand side of a list of bodies  

python usage : mecaMAILx_AssembRHS(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to assemble right hand side if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeResidueNorm "

computes the norm of the residue of a list of bodies  

python usage : norm = mecaMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute the norm of the residue if
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
i_list(list of integer) : list of bodies to compute stiffness and viscosity
    matrices and internal forces if omitted works on all objects  
";

%feature("docstring") mecaMAILx_ComputeField "

computes elementary fields of a list of bodies  

python usage : mecaMAILx_ComputeField(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeFint "

computes elementary internal forces of a list of bodies  

python usage : mecaMAILx_ComputeFint(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute internal forces if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_UpdateBulk "

update begin elementary fields with current elementary fields of a list of
bodies  

python usage : mecaMAILx_UpdateBulk(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_UpdateDof "

update begin d.o.f. with current d.o.f. of a list of bodies  

python usage : mecaMAILx_UpdateDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to update current d.o.f if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_ComputeDof "

computes the current d.o.f knowing all the forces (free + contact) of a list of
bodies  

python usage : mecaMAILx_ComputeDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute current d.o.f if omitted works
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
i_list(list of integer) : list of bodies to compute external forces if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeMass "

compute elementary mass and inertia of a list of bodies  

python usage : mecaMAILx_ComputeMass(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute mass and inertia if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_FatalDamping "

set to 0 current velocities of a list of bodies  

This keyword must be between the ComputeDof and UpdateDof ones.  

python usage : mecaMAILx_FatalDamping(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to reset current velocity if omitted
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
checktype(char[5]) : type of check test  
tol(double) : tolerance  
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
num(integer) : which GPV file to read  
";

%feature("docstring") mecaMAILx_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

python usage : mecaMAILx_ReadIniDof(num=0)  

Parameters
----------
num(integer) : which DOF file to read  
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
IdBody(int) : id of the concern body  
IdElem(int) : id of the concern element python usage :
    mecaMAILx_DisplayBulkElement(IdBody,IdElem)  
";

%feature("docstring") mecaMAILx_WriteLastRnod "

Write ascii Rnod.LAST file of a list of bodies.  

python usage : mecaMAILx_WriteLastRnod(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write in Rnod.LAST if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_WriteOutRnod "

Write ascii Rnod.OUT file of a list of bodies. Can be activat only each N step.  

python usage : mecaMAILx_WriteOutRnod(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write in Rnod.OUT if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_DisplayOutRnod "

Display body forces of a list of bodies.  

python usage : mecaMAILx_DisplayOutRnod(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to display body forces if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_WriteLastNodalForces "

Write ascii Rnod.LAST file of a list of bodies.  

This function is almost like WriteLastRnod, but write also internal and inertial
forces.  

python usage : mecaMAILx_WriteLastNodalForces(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write in Rnod.LAST if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_WriteOutNodalForces "

Write ascii Rnod.OUT file of a list of bodies. Can be activat only each N step.  

This function is almost like WriteOutRnod, but write also internal and inertial
forces.  

python usage : mecaMAILx_WriteOutNodalForces(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write in Rnod.OUT if omitted works on
    all objects  
";

%feature("docstring") mecaMAILx_DisplayOutNodalForces "

Display computed nodal forces of a list of bodies.  

python usage : mecaMAILx_DisplayOutNodalForces(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to display body forces if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_GetScalarFieldRank "

Get the rank of scalar field of an element of a body from its name.  

python usage : f_rank = mecaMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
  
 You need to declare this field in your MODELS.DAT  
";

%feature("docstring") mecaMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

python usage : mecaMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
  
 You need to declare this field in your MODELS.DAT  
";

%feature("docstring") mecaMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = mecaMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired vector field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") mecaMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : mecaMAILx_SetVectorFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
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
i_list(list of integer) : list of bodies to compute ortho frame with user routine
    if omitted works on all objects  
";

%feature("docstring") mecaMAILx_ComputeUserField "

Use user routine to compute a field at gp.  

python usage : mecaMAILx_ComputeUserField(ifield, i_list)  

Parameters
----------
ifield(integer) : id of the field to compute  
i_list(list of integer) : list of bodies to compute user fields on if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_SetVisible "

set visible a given mecaMAILx  

python usage : mecaMAILx_SetVisible(ibdyty)  

Parameters
----------
ibdyty(integer): index of the mecaMAILx  
";

%feature("docstring") mecaMAILx_SetInvisible "

rended a given mecaMAILx invisible  

python usage : mecaMAILx_SetInvisible(ibdyty)  

Parameters
----------
ibdyty(integer): index of the mecaMAILx  
";

%feature("docstring") mecaMAILx_IsVisible "

return if a given body visible  

python usage : visible = mecaMAILx_IsVisible(ibdyty)  

Parameters
----------
idbdy(integer): id of the body we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") mecaMAILx_ComputeRayleighDamping "

compute the Rayleigh damping: C=alpha*M+beta*K of a list of bodies  

python usage : mecaMAILx_ComputeRayleighDamping(alpha,beta,i_list)  

Parameters
----------
alpha(real) : damping value  
beta(real) : damping value  
i_list(list of integer) : list of bodies to compute Rayleigh damping if omitted
    works on all objects  
";

%feature("docstring") mecaMAILx_ComputeRayleighDampingDiscreteElement "

set damping for discrete FE element of a list of bodies  

python usage : mecaMAILx_ComputeRayleighDampingDiscreteElement(damp, i_list)  

Parameters
----------
ref_size(real) : damping value  
i_list(list of integer) : list of bodies to compute damping for discrete FE
    element if omitted works on all objects  
";

%feature("docstring") mecaMAILx_GetNodeCoorTT "

return TT node coordinates  

python usage : vec = mecaMAILx_GetNodeCoorTT(ibdyty,inodty)  

Parameters
----------
idbdy(integer): id of the body  
inodty(integer): id of the node  

Returns
-------
vec (float vector) : TT node coordinates  
";

%feature("docstring") mecaMAILx_GetNodeCooref "

return ref node coordinates  

python usage : vec = mecaMAILx_GetNodeCoorref(ibdyty,inodty)  

Parameters
----------
idbdy(integer): id of the body  
inodty(integer): id of the node  

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
datatype(string of size 5) : the matrix to get  
ibdyty(integer) : rank of considered body  

Returns
-------
matrix (double array) : the desired matrix  
";

%feature("docstring") mecaMAILx_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = mecaMAILx_getDrvVlocy(ibdyty)  

Parameters
----------
ibdyty(integer) : index of the mecaMAILx  
drvdof_indices(integer array) : indices list of driven dof  
drvdof_values(real array) : values of the driven dof  
";

%feature("docstring") mecaMAILx_computeDrvVlocy "

Compute the value of the driven velocity of a body a current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : mecaMAILx_computeDrvVlocy(ibdyty, values)  

Parameters
----------
ibdyty(integer) : index of the mecaMAILx  
values(double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") mecaMAILx_SetVlocyDrivenDof "

Apply Drv Dof on a given body.  

python usage : mecaMAILx_SetVlocyDrivenDof(IdBody, f_dof, f_node, f_value)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_dof(integer) : dof of the concern node  
f_node(integer) : node  
f_value(double) : value of the drvdof  
";

%feature("docstring") mecaMAILx_ComputeContactDetectionConfiguration "

compute the contact detection configuration of a list of bodies  

python usage : mecaMAILx_ComputeContactDetectionConfiguration(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute contact detection
    configuration if omitted works on all objects  
";

%feature("docstring") mecaMAILx_NullifyReac "

set to 0 the reac of the IdBody mecaMAILx  

python usage : mecaMAILx_NullifyReac(datatype, IdBody)  

Parameters
----------
datatype(string of size 5) : the vector to set  
IdBody(integer) : id of the concerned body  
";

%feature("docstring") mecaMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = mecaMAILx_GetAll(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") mecaMAILx_GetCooref "

return node coordinates of idBody  

python usage : array = mecaMAILx_GetCooref(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") mecaMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = mecaMAILx_GetConnectivity(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") mecaMAILx_GetElementsVolume "

return volume of elements  

python usage : volumes = mecaMAILx_GetElementsVolume(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
volumes[nb_ele] (double) : volume  
";

%feature("docstring") mecaMAILx_GetGpCoor "

return Gauss points coordinates of idBody  

python usage : array = mecaMAILx_GetGpCoor(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates of all Gauss points  
";

%feature("docstring") mecaMAILx_GetGpStrain "

return strain values stored at a gp  

python usage : strain = mecaMAILx_GetGpStrain(idBody,idEle,idGp)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  
IdGp(integer) : id of the concerned gauss point  

Returns
-------
strain[size] (double) : value of strain  
";

%feature("docstring") mecaMAILx_GetGpStress "

return stress values stored at a gp  

python usage : stress = mecaMAILx_GetGpStress(idBody,idEle,idGp)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  
IdGp(integer) : id of the concerned gauss point  

Returns
-------
stress[size] (double) : value of stress  
";

%feature("docstring") mecaMAILx_GetGpInternals "

return internal values stored at a gp  

python usage : internals = mecaMAILx_GetGpInternals(idBody,idEle,idGp)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  
IdGp(integer) : id of the concerned gauss point  

Returns
-------
internals[nb_internals] (double) : value of internals  
";

%feature("docstring") mecaMAILx_GetGpPrincipalField "

return principal field (strain or stress) at a gp  

python usage : field = mecaMAILx_GetGpPrincipalField(idBody,idEle,idGp,idField)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  
IdGp(integer) : id of the concerned gauss point  
IdField(integer): id of the field (1: strain, 2: stress)  

Returns
-------
field (double array): tensor field with principal values  
";

%feature("docstring") mecaMAILx_GetElementsInternal "

return a value over elements of an internal stored at gp  

python usage : internals = mecaMAILx_GetElementsInternal(idBody,id,f)  

Parameters
----------
IdBody(integer) : id of the concerned body  
Id(integer) : id of the internal  
f(integer) : flag 1: mean, 2: sum, 3:max, 4: min  

Returns
-------
internals[nb_ele] (double) : value of internal  
";

%feature("docstring") mecaMAILx_GetElementsInternalIntegral "

return integral over elements of an internal stored at gp  

python usage : internals = mecaMAILx_GetElementsInternalIntegral(idBody,id)  

Parameters
----------
IdBody(integer) : id of the concerned body  
Id(integer) : id of the internal  

Returns
-------
internals[nb_ele] (double) : value of internal  
";

%feature("docstring") mecaMAILx_GetElementsCenter "

return center of elements  

python usage : centers = mecaMAILx_GetElementsCenter(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
centers[3*nb_ele] (double) : center  
";

%feature("docstring") mecaMAILx_GetElementsJacobian "

return jacobian of elements  

python usage : jacobians = mecaMAILx_GetElementsJacobian(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

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
IdBody(integer) : id of the concerned body  

Returns
-------
energies (double array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_GetElementsNeighbor "

return elements in the tol-neighbor of an element of idBody  

python usage : neighbors = mecaMAILx_GetElementsNeighbor(idBody,tol)  

Parameters
----------
IdBody(integer) : id of the concerned body  
tol(double) : tolerance  

Returns
-------
array (double 2D-array) : neighbor[nb_ele,max_neighbors]  
";

%feature("docstring") mecaMAILx_GetPtrElementsVisibility "

Get a pointer on the elements visibility vector.  

python usage : eviz = mecaMAILx_GetPtrElementsVisibility(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the mecaMAILx  

Returns
-------
eviz (int array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_AddNodalFieldDivergence "

Add the divergence of a diagonal field to external forces.  

python usage : mecaMAILx_AddNodalFieldDivergence(ibdyty, ifield)  

Parameters
----------
ibdyty(integer) : rank of body  
ifield(integer) : rank of field  
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
ibdyty(integer) : rank of considered body  

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
ibdyty(integer) : rank of considered body  
displacement(double vector) : displacement field  

Returns
-------
energy (double) : deformation energy  
";

%feature("docstring") mecaMAILx_GetKineticEnergy "

Get the kinetic energy of a given velocity field.  

python usage : energy = mecaMAILx_GetKineticEnergy(id,velocity)  

Parameters
----------
ibdyty(integer) : rank of considered body  
velocity(double vector) : velocity field  

Returns
-------
energy (double) : kinetic energy  
";

%feature("docstring") mecaMAILx_GetNeighborElementsToElement "

return neighbor elements to element idEle of body idBody  

python usage : vector = mecaMAILx_GetNeighborElementsToElement(idBody,idEle)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  

Returns
-------
vector (integer) : list of elements  
";

%feature("docstring") mecaMAILx_GetNeighborElementsToNode "

return neighbor elements to node idNode of body idBody  

python usage : vector = mecaMAILx_GetNeighborElementsToNode(idBody,idNode)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdNode(integer) : id of the concerned node  

Returns
-------
vector (integer) : list of elements  
";

%feature("docstring") mecaMAILx_GetBoundaryElements "

return boundary elements  

python usage : vector = mecaMAILx_GetBoundaryElements(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

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
ivalue(integer) : id of body to set precon  
";

%feature("docstring") mecaMAILx_GetPtrPreconW "

Get a pointer on the preconW Matrix of a given body.  

python usage : pcW = mecaMAILx_GetPtrPreconW(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the mecaMAILx  

Returns
-------
pcW (double array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") mecaMAILx_GetInternalVariable "

Get a copy of the internal variable of a given body.  

Python usage : internal = mecaMAILx_GetInternalVariable(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

Returns
-------
internal (double array) : internal variable of desired body  
";

%feature("docstring") mecaMAILx_GetNbInternal "

Get the number of internal variable of a given body.  

python usage : nb_internal = mecaMAILx_GetNbInternal(ibdyty)  

Parameters
----------
ivalue(integer) : rank of the body  

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
cvalue1_c(string of size 5) : name of the body vector  
IdBody(integer) : id of the body  

Returns
-------
vector_ptr (double array) : reference on the desired body vector  
";

%feature("docstring") mecaMAILx_GetDofStatus "

Get the status of nodes: 0 free, 1 x, 10 y.  

Python usage : vector = mecaMAILx_GetDofStatus(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired data  
";

%feature("docstring") mecaMAILx_PrepGlobalSolver "

computes free velocity of a list of bodies  

python usage : mecaMAILx_PrepGlobalSolver(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute free velocity if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_PostGlobalSolver "

computes the current d.o.f knowing all the forces (free + contact) of a list of
bodies  

python usage : mecaMAILx_PostGlobalSolver(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute current d.o.f if omitted works
    on all objects  
";

%feature("docstring") mecaMAILx_AddBodyForceToFext "

Add a body force (M*gamma) to Fext for a given body.  

python usage : mecaMAILx_AddBodyForceToFext(ibdyty, matrix)  

Parameters
----------
ibdyty(integer) : rank of body  
matrix(double array) : the new value  
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
scale(double) : scaling  
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
ibdyty(integer): index of the mecaMAILx  
inod(integer): index of the node to set visible  
idof(integer): index of the dof of the node to set visible  
";

%feature("docstring") mecaMAILx_SetInvisibleVlocyDrivenDof "

allows to deactivate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : mecaMAILx_SetInvisibleVlocyDrivenDof(ibdyty, inod, idof)  

Parameters
----------
ibdyty(integer): index of the mecaMAILx  
inod(integer): index of the node to set invisible  
idof(integer): index of the dof of the node to set invisible  
";

%feature("docstring") mecaMAILx_UpdateVlocyDrivenDofStructures "

takes into account modifications on Vlocy driven dof status  

python usage : mecaMAILx_UpdateVlocyDrivenDofStructures(ibdyty)  

Parameters
----------
ibdyty(integer): index of the mecaMAILx  
";


// File: wrap__PRPRx_8h.xml

%feature("docstring") PRPRx_SelectProxTactors "

contact detection between PRxxx and PRxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PRPRx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") PRPRx_UseCpCundallDetection "

chooses the Cundall iterative detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpCundallDetection(nb_iter, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
nb_iter(integer) : max number of iterations  
cd_shrink(real) : shrink parameter (length) in clipper for candidate  
an_shrink(real) : shrink parameter (length) in clipper for antagonist  
delta(real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseCpF2fExplicitDetection "

chooses the face 2 face combinatory detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpF2fExplicitDetection(tol, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
tol(real) : tolerance on normal orientations  
cd_shrink(real) : shrink parameter (length) in clipper for candidate  
an_shrink(real) : shrink parameter (length) in clipper for antagonist  
delta(real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseCpF2fDetection "

chooses a mix of the face 2 face and Cundall detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpF2fDetection(tol, iter, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
tol(real) : tolerance on normal orientations  
iter(integer) : max number of iterations  
cd_shrink(real) : shrink parameter (length) in clipper for candidate  
an_shrink(real) : shrink parameter (length) in clipper for antagonist  
delta(real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseNcDetection "

chooses contact detection methode between non-convex shapes  

python usage : PRPRx_UseNcDetection(gdist)  

Parameters
----------
gdist(real) : global distance  
";

%feature("docstring") PRPRx_UseNcF2fDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

python usage : PRPRx_UseNcF2fDetection(gdist,tol)  

Parameters
----------
gdist(real) : global distance  
tol(real) : tolerance on normal orientations  
";

%feature("docstring") PRPRx_UseNcF2fExplicitDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

python usage : PRPRx_UseNcF2fExplicitDetection(gdist,tol)  

Parameters
----------
gdist(real) : global distance  
tol(real) : tolerance on normal orientations  
";

%feature("docstring") PRPRx_UseTrianglesIntersectionDetection "

chooses contact detection finding intersection in a soup of triangles.  

The number of point provided is an internal parameter of the algorithm which
control the maximum number of intersection points stored when looking for the
triangles intersection before restricting it to only 4 of them. So it must be
strictly superior to 4.  

python usage : PRPRx_UseTrianglesIntersectionDetection(nb_max_pt=16)  

Parameters
----------
nb_max_pt(integer): maximum contact points to store/check during detection  
";

%feature("docstring") PRPRx_SetF2fMinimalSurfaceSize "

set the minimum contact surface size with f2f algo otherwize contact is not
computed  

python usage : PRPRx_SetF2fMinimalSurfaceSize(tol)  

Parameters
----------
tol(real) : minimum surface size  
";

%feature("docstring") PRPRx_UseExternalDetection "

chooses external contact detection (bindings)  

python usage : PRPRx_UseExternalDetection()  
";

%feature("docstring") PRPRx_WriteLastVlocRloc "

write last local values of all PRPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : PRPRx_WriteLastVlocRloc()  
";

%feature("docstring") PRPRx_WriteOutVlocRloc "

write local values of all PRPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : PRPRx_WriteOutVlocRloc()  
";

%feature("docstring") PRPRx_DisplayOutVlocRloc "

display local values of all PRPRx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : PRPRx_DisplayOutVlocRloc()  
";

%feature("docstring") PRPRx_DisplayProxTactors "

display contacts  

python usage : PRPRx_DisplayProxTactors()  
";

%feature("docstring") PRPRx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : PRPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PRPRx_ShrinkPolyrFaces "

Shrink the face of the candidate polyhedron for the detection.  

May be conflicting with the shrink parameters of the detections functions used
by clipper library. The difference is that clipper use a single length for all
sample, whereas this function use a scale factor to retract the vertices of the
candidate polyhedron inside the the surface.  

python usage : PRPRx_ShrinkPolyrFaces(shrink)  

Parameters
----------
shrink(real) : scale factor allowing to shrink candidate surface  

    1.  no shrink, 1. no surface  
";

%feature("docstring") PRPRx_LowSizeArrayPolyr "

abscons parameter to manage memory allocation  

python usage : PRPRx_LowSizeArrayPolyr(sfactor)  

Parameters
----------
sfactor(integer) :  
";

%feature("docstring") PRPRx_SaveProxTactorsToFile "

write selected contacts to file  

python usage : PRPRx_SaveProxTactorsToFile()  
";

%feature("docstring") PRPRx_LoadProxTactorsFromFile "

load selected contact from files  

python usage : PRPRx_LoadProxTactorsFromFile()  
";

%feature("docstring") PRPRx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PRPRx_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(real) : periode on x axis  
";

%feature("docstring") PRPRx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PRPRx_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(real) : period on y axis  
yperiod(double) : period on y axis  
";

%feature("docstring") PRPRx_VerboseF2F "

ask for verbose comment concerning contact detection between cd and an  

python usage : PRPRx_VerboseF2F(cd,an)  

Parameters
----------
cd(integer) : candidate  
an(integer) : antagoniste  
";

%feature("docstring") PRPRx_GetNbF2f "

Get the number of f2f structures stored This is the real size of the array, and
not the number of active f2f structure.  

python usage : nb_f2f = PRPRx_GetNbF2f()  

Returns
-------
nb_f2f (integer) : the size of the f2f array  
";

%feature("docstring") PRPRx_GetF2f2Inters "

Get the list of interactions for each face-to-face structure Array of integer
with number of f2f, then for each f2f, the number of interactions then the list
of interaction id.  

python usage : f2f_inters = PRPRx_GetF2f2Inters()  

Returns
-------
f2f_inters (integer array) : the integer array  
";

%feature("docstring") PRPRx_GetF2fOutlines "

Get the connectivity of all intersection polytopes of all face2face and the
corresponding coordinates.  

The connectivity containes first the number of f2f, then for each, the number of
polytope, then for each the number of vertices.  

The coordinates must be counted from this ordering...  

python usage : connec, points = PRPRx_GetF2fOutlines()  

Returns
-------

*   connec (integer array) : the connectivities  
*   points (double array) : the coordinates  
";

%feature("docstring") PRPRx_GetF2fAllIdata "

Get topological face id of cd/an for all F2f structure.  

python usage : idata = PRPRx_GetF2fAllIdata()  

Returns
-------

*   idata (integer array) : size [nb_f2fx2] with the face id  
";

%feature("docstring") PRPRx_GetF2fCentralKernel "

Give the central kernel coordinates, the equivalent normal stress and if the
center of pressure is inside.  

python usage : ck_coor, sn, is_in = PRPRx_GetF2fStress(i_f2f)  

Returns
-------  
";

%feature("docstring") PRPRx_GetF2fStress "

Give the polygons of the compressed and decompressed part and linear stress
repartition.  

In the case when the minimization algorithm failed, the decompression value is
set to -99. so that when writing the vtk files, the 'ids' numbering is kept
consistent.  

python usage : coorC, sizeC, coorD, sizeD, sigma, decomp =
PRPRx_GetF2fStress(i_f2f)  

Returns
-------  
";

%feature("docstring") PRPRx_SetCundallNeighbor "

set a neighbor distance around common plane to select projected nodes  

python usage : PRPRx_SetCundallNeighbor(neighbor)  

Parameters
----------
neighbor(real) : ratio of a reference size  
";

%feature("docstring") PRPRx_CpUseOldCcpm "

use the old method for computing contact point position  

python usage : PRPRx_CpUseOldCcpm()  
";

%feature("docstring") PRPRx_SetReactionTrackingLength "

function which makes possible to set the length of the hexaedra glyph
representing the visavis reaction  

python usage : PRPRx_SetReactionTrackingLength(length)  

Parameters
----------
length(real) : length the hexaedra glyph  
";

%feature("docstring") PRPRx_SetTolRecupRloc "

set the distance tolerance used in PRPRx_RecupRloc  

python usage : PRPRx_SetTolRecupRloc(tol)  

Parameters
----------
tol(double) : tolerance  
";

%feature("docstring") PRPRx_GetInteractionVector "

Get a copy of a vector of a PRPRx.  

possible values for datatype field are \"Coor_\", \"N____\"  

python usage : vector = PRPRx_GetInteractionVector(datatype, icdan)  

Parameters
----------
datatype(string [5]) : the vector to get  
icdan(integer) : rank of the PRPRx  

Returns
-------
vector (double array) : output vector  
";

%feature("docstring") PRPRx_SetInteractionInternal "

Set a value of the internal vector of a PRPRx.  

python usage : PRPRx_SetInteractionInternal(i, icdan, value)  

Parameters
----------
i(integer) : rank of internal  
icdan(integer) : rank of the PRPRx  
value(double) : value to set  
";

%feature("docstring") PRPRx_GetInteractionInternal "

Get a value from the internal vector of a PRPRx.  

python usage : value = PRPRx_GetInteractionInternal(i, icdan)  

Parameters
----------
i(integer) : rank of internal  
icdan(integer) : rank of the PRPRx  
value(double) : value to get  
";

%feature("docstring") PRPRx_GetInteractionInternalComment "

Get internal comment of a given interaction.  

python usage : comment=PRPRx_GetInteractionInternalComment(icdan)  

Parameters
----------
icdan(integer) : rank of the PRPRx  

Returns
-------
comment (char[100]) : the string to get  
";

%feature("docstring") PRPRx_WithNodalContact "

use cd contact points at nodes instead at faces with NcDetection  

python usage : PRPRx_WithNodalContact()  
";

%feature("docstring") PRPRx_SetInternalSurface "

Set the value of a surface type (point, line or surf) for wti detection.  

For surface, if the value is left to 0., then the surface of the triangle is
computed To select the type of surface : 1->point, 2->line, 3->surface  

python usage : PRPRx_SetInternalSurface(itype, value)  

Parameters
----------
itype(integer) : the type of surface to set  
value(double) : value to set  
";

%feature("docstring") PRPRx_UseStoDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

Face to face detection implemented by Stono which can mix between the standard
f2f detection and the non convex one. Furthermor the decompression parameter can
help with putting the contact points either near the  

python usage : PRPRx_UseFCDetection(explicite, decompression, tol, kappa)  

Parameters
----------
explicite(boolean) : use explicit detection  
decompression(double) : surface decompression (value in [-1., 1.])  
tol(real) : tolerance on normal orientations  
kappa(boolean) : compute kappas coefficient  
";

%feature("docstring") PRPRx_ForceF2fDetection "

force f2f detection method even for non-convex surfaces  

python usage : PRPRx_ForceF2fDetection()  
";

%feature("docstring") PRPRx_ForceNcDetection "

force nc detection method even for flat surfaces  

python usage : PRPRx_ForceNcDetection()  
";

%feature("docstring") PRPRx_CleanMemory "

Free all memory allocated within PRPRx module.  

python usage : PRPRx_CleanMemory()  
";


// File: wrap__deposit2D_8h.xml

%feature("docstring") deposit2D_Potential "

Computes a new deposit under potential with or without big particles.  

python call: coor = deposit2D_Potential(in_radii, lx, potential[, dradii,
dcoor])  

Parameters
----------
in_radii(double array): given radii list (i.e. granulometry)  
lx(double): width of the box in which to deposit  
potential(integer): for deposit (1->gravity, 2->wall, 3->big_particles)  
dradii(double array) (optional) : a list of already deposited radii  
dcoor(double array) (optional) : a list of already deposited coor (must be of
    size [nb_dradii,3])  

Returns
-------
coor (double array): coordinates of deposited radii (shape [nb_radii,2]) PYDOC  
";


// File: wrap__SPHER_8h.xml

%feature("docstring") SPHER_LoadTactors "

load SPHER from RBDY3 and initialize existing_entites  

python usage : SPHER_LoadTactors()  
";

%feature("docstring") SPHER_SetRadiusCorrection "

set a radius correction  

python usage : SPHER_SetRadiusCorrection(corr)  

Parameters
----------
corr(real) :  
";

%feature("docstring") SPHER_GetNbSPHER "

Get the number of SPHER.  

python usage : nb_SPHER = SPHER_GetNbSPHER()  

Returns
-------
nb_SPHER (integer) : the number of SPHER  
";

%feature("docstring") SPHER_GetSPHER2BDYTY "

Get a copy of map SPHER2bdyty.  

usage : polyr2bdyty = SPHER_GetSPHER2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") SPHER_GetPtrSPHER2BDYTY "

return a pointer onto the map spher2bdyty  

python usage : spher2bdyty = SPHER_GetPtrSPHER2BDYTY()  

Returns
-------
spher2bdyty (integer array) : reference on map between spher rank and body rank  
";

%feature("docstring") SPHER_GetContactorRadius "

Get the radius of a SPHER contactor.  

python usage : radius = SPHER_GetContactorRadius(itact)  

Parameters
----------
itact(integer) : id of a SPHER  

Returns
-------
radius (double) : the radius of the SPHER number itact  
";

%feature("docstring") SPHER_GetContactorCoor "

get coordinates of the center of a given SPHER  

usage : vector = SPHER_GetContactorCoor(itacty)  

Parameters
----------
itacty(integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") SPHER_GetContactorCoorb "

get coordinates at the begin of the time step of the center of a given SPHER  

usage : vector = SPHER_GetContactorCoorb(itacty)  

Parameters
----------
itacty(integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") SPHER_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = SPHER_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") SPHER_InitOutlines "

Get a reference on the outlines of all SPHER.  

usage : outlines = SPHER_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_SPHER  
";

%feature("docstring") SPHER_InitScalarFields "

Get a reference on the scalar fields of all SPHER.  

usage : scalarfields = SPHER_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_SPHER array  
";

%feature("docstring") SPHER_UpdatePostdata "

Update values of outlines_SPHER and scalarfields_SPHER pointers.  

usage : SPHER_UpdatePostdata  
";

%feature("docstring") SPHER_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = SPHER_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
SPHER  
";

%feature("docstring") SPHER_GetNbScalarFields "

Get the number of scalar fields of a SPHER.  

python usage : nb_scalarfields = SPHER_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a SPHER  
";

%feature("docstring") SPHER_GetPtrAllConnectivities "

Get a reference on the connectivities of all SPHER.  

usage : connec = SPHER_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") SPHER_CleanMemory "

Free all memory allocated within SPHER module.  

python usage : SPHER_CleanMemory()  
";


// File: wrap__CLJCx_8h.xml

%feature("docstring") CLJCx_SelectProxTactors "

contact detection between CLxxx and JCxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : CLJCx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CLJCx_WriteLastVlocRloc "

write last local values of all CLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : CLJCx_WriteLastVlocRloc()  
";

%feature("docstring") CLJCx_WriteOutVlocRloc "

write local values of all CLJCx contacts  

The values written are relative velocity, forces and local frame  

python usage : CLJCx_WriteOutVlocRloc()  
";

%feature("docstring") CLJCx_DisplayOutVlocRloc "

display local values of all CLJCx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CLJCx_DisplayOutVlocRloc()  
";

%feature("docstring") CLJCx_DisplayProxTactors "

display contacts  

python usage : CLJCx_DisplayProxTactors()  
";

%feature("docstring") CLJCx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CLJCx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CLJCx_CleanMemory "

Free all memory allocated within CLJCx module.  

python usage : CLJCx_CleanMemory()  
";


// File: wrap__PTPT3_8h.xml

%feature("docstring") PTPT3_SelectProxTactors "

contact detection between PTxxx and PTxxx tactors  

python usage : PTPT3_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PTPT3_SmoothForceComputation "

computes smooth forces (if any)  

python usage : PTPT3_SmoothForceComputation()  
";

%feature("docstring") PTPT3_WriteLastVlocRloc "

write last local values of all PTPT3 contacts  

python usage : PTPT3_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_WriteOutVlocRloc "

write local values of all PTPT3 contacts  

python usage : PTPT3_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_DisplayOutVlocRloc "

display local values of all PTPT3 contacts  

python usage : PTPT3_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PTPT3_DisplayProxTactors "

display contacts  

python usage : PTPT3_DisplayProxTactors()  
";

%feature("docstring") PTPT3_LoadNetwork "

read a PTPT3 network from a file  

python usage : PTPT3_LoadNetwork()  
";

%feature("docstring") PTPT3_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PTPT3_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PTPT3_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PTPT3_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(real) : period on x axis  
";

%feature("docstring") PTPT3_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PTPT3_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(real) : period on y axis  
";

%feature("docstring") PTPT3_SetExplicitLocalFrame "

local frame is computed only once at the first step  

python usage : PTPT3_SetExplicitLocalFrame()  
";

%feature("docstring") PTPT3_LoadParams "

read a PTPT3 surface and l0 from a file  

python usage : PTPT3_LoadParams()  
";

%feature("docstring") PTPT3_UseCurrentNonuc0 "

Use GetCoor or value given from file insted of computing nonuc0 from reference
coordinates.  

python usage : PTPT3_UseCurrentNonuc0(to_use) param[in] to_use (integer) : 1 to
activate, 0 to deactivate feature  
";

%feature("docstring") PTPT3_CleanMemory "

Free all memory allocated within PTPT3 module.  

python usage : PTPT3_CleanMemory()  
";


// File: wrap__PLJCx_8h.xml

%feature("docstring") PLJCx_SelectProxTactors "

contact detection between POLYG and JONCx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLJCx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
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
num(integer) : which VlocRloc file to read  
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
cflag(char) : model to use ('min', 'max' or 'ave')  
";


// File: wrap__SPPLx_8h.xml

%feature("docstring") SPPLx_SelectProxTactors "

contact detection between SPxxx and PLxxx tactors  

python usage : SPPLx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPPLx_SmoothForceComputation "

compute smooth contact law (in any)  

python usage : SPPLx_SmoothForceComputation()  
";

%feature("docstring") SPPLx_WriteLastVlocRloc "

write last local values of all SPPLx contacts  

python usage : SPPLx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_WriteOutVlocRloc "

write local values of all SPPLx contacts  

python usage : SPPLx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_DisplayOutVlocRloc "

display local values of all SPPLx contacts  

python usage : SPPLx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_DisplayProxTactors "

display contacts  

python usage : SPPLx_DisplayProxTactors()  
";

%feature("docstring") SPPLx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") SPPLx_CleanMemory "

Free all memory allocated within SPPLx module.  

python usage : SPPLx_CleanMemory()  
";


// File: wrap__global__thermal__solver_8h.xml

%feature("docstring") gts_Initialize "

Initialize global solver module.  

python usage : gts_Initialize()  
";

%feature("docstring") gts_AssembleSystem "

Assembling of the global system.  

python usage : gts_AssembleSystem()  
";

%feature("docstring") gts_PrepSystem "

Preparing the global system.  

python usage gts_PrepSystem()  
";

%feature("docstring") gts_AssembleLHS "

Assembling the lhs of the global system.  

python usage : gts_AssemblerLHS()  
";

%feature("docstring") gts_AssembleRHS "

Assembling the rhs of the global system.  

i  

python usage : gts_AssembleRHS()  
";

%feature("docstring") gts_Solve "

Solving of the global system.  

i  

python usage : gts_Solve()  
";

%feature("docstring") gts_Finalize "

Clean memory of global solver module.  

python usage : gts_Finalize()  
";


// File: wrap__mbs3D_8h.xml

%feature("docstring") MBS3D_setNb "

Set the number of MBS.  

python usage : MBS3D_setNb(nb)  

Parameters
----------
nb(integer) : set the number of MBS  
";

%feature("docstring") MBS3D_getNb "

Get the number of MBS.  

python usage : nb = MBS3D_getNb()  

Returns
-------
nb (integer) : the number of MBS  
";

%feature("docstring") MBS3D_setNbNodes "

Set the number of nodes of a MBS.  

python usage : MBS3D_setNbNodes(ibdyty, nb)  

Parameters
----------
ibdyty(integer): id of the MBS  
nb(integer) : the number of nodes of the MBS  
";

%feature("docstring") MBS3D_setNbTactors "

Set the number contactors of a MBS.  

python usage : MBS3D_setNbTactors(ibdyty, nb)  

Parameters
----------
ibdyty(integer): id of the MBS  
nb(integer) : the number of contactor of the MBS  
";

%feature("docstring") MBS3D_getPtrCoor "

Get a pointer on the coor of a MBS.  

usage : coor = MBS3D_GetPtrCoor(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered MBS  

Returns
-------
coor (double 2D-array) : reference on the coordinates of the nodes  
";

%feature("docstring") MBS3D_getPtrCoorTT "

Set the array of coordinates of nodes of a MBS.  

python usage : coor = MBS3D_getPtrCoorTT(ibdyty)  

Parameters
----------
ibdyty(integer): id of the MBS  

Returns
-------
coor (double array) : coordinates of nodes of a MBS (in contact configuration)  
";

%feature("docstring") MBS3D_getPtrLocalFrame "

Get a pointer on the coor of a MBS.  

usage : frame = MBS3D_GetPtrLocalFrame(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered MBS  

Returns
-------
frame (double 2D-array) : local frame  
";

%feature("docstring") MBS3D_getPtrLocalFrameTT "

Set the array of coordinates of nodes of a MBS.  

python usage : frameTT = MBS3D_GetPtrLocalFrameTT(ibdyty)  

Parameters
----------
ibdyty(integer): id of the MBS  

Returns
-------
frameTT (double array) : local frame (in contact configuration)  
";

%feature("docstring") MBS3D_addContactor "

Add a new contactor to a MBS.  

Available contactor types are :  

*   PLANx: inputs are:
    -   rdata must hold [axe_x, axe_y, axe_z]  
*   POLYR: inputs are:
    -   rdata must hold the coordinates of the vertices [x_1, y_1, z_1, ... x_n,
        y_n, z_n]  
    -   idata must hold the connecivity of each triangle defining the surface  

python usage : MBS3D_addContactor(ibdyty, inodty, itacty, tacttype, color,
rdata, idata=None)  

Parameters
----------
ibdyty(integer) : rank of the MBS  
inodty(integer) : rank of the node of the MBS the contactor is tied to  
itacty(integer) : rank of the contactor of MBS  
tactype(string [5]) : the type of contactor  
color(string [5]) : the color of the contactor  
rdata(double array) : the new value of the vector  
idata(integer array) : the new value of the vector  
";

%feature("docstring") MBS3D_initialize "

Initialize MBS module once loading is done.  

python usage : MBS3D_initialize()  
";

%feature("docstring") MBS3D_finalize "

Finalize MBS module.  

python usage : MBS3D_finalize()  
";

%feature("docstring") MBS3D_IncrementStep "

compute the current velocity and displacement  

python usage : MBS3D_IncrementStep()  
";

%feature("docstring") MBS3D_ComputeFreeVelocity "

compute free velocity  

python usage : MBS3D_ComputeFreeVelocity()  
";

%feature("docstring") MBS3D_ComputeDof "

update current position and velocity  

python usage : MBS3D_ComputeDof()  
";

%feature("docstring") MBS3D_UpdateDof "

save d.o.f. of the end of the time step to d.o.f. of the begining of the next
one  

python usage : MBS3D_UpdateDof()  
";


// File: wrap__DKJCx_8h.xml

%feature("docstring") DKJCx_SelectProxTactors "

contact detection between DISKx and JONCx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKJCx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
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
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKJCx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKJCx_SetSurfaceSectors(nbsect)  

Parameters
----------
nbsect(integer) : number of sectors  
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
cflag(char) : model to use ('min', 'max' or 'ave')  
";


// File: wrap__DDM__3D_8h.xml

%feature("docstring") DDM_3D_SetDDWorkingDirectory "
";

%feature("docstring") DDM_3D_Initialize "

Initialize 3D DDM solver.  

ddm_type may be :  

*   1 for Feti (DDM without overlap)  
*   2 for Schwarz (DDM with overlap)  

python usage : DDM_3D_Initialize(nb_sdmx, nb_sdmy, nb_sdmz, ddm_type)  

Parameters
----------
nb_sdmx(integer) : number of domains on x-axis  
nb_sdmy(integer) : number of domains on y-axis  
nb_sdmz(integer) : number of domains on z-axis  
ddm_type(integer) : type of DDM to use  
";

%feature("docstring") DDM_3D_Partitioning "
";

%feature("docstring") DDM_3D_ExperimentalPartitioning "
";

%feature("docstring") DDM_3D_AddToFext "

Add external forces due to DDM.  

python usage : DDM_3D_AddToFext()  

To add after RBDY2_ComputeFext  
";

%feature("docstring") DDM_3D_SelectProxTactors "

DDM way to compute contact.  

python usage : DDM_3D_SelectProxTactors()  

To add after RBDY3_ComputeFreeVelocity()  
";

%feature("docstring") DDM_3D_ExSolver "

Solve fully the local contact problem with DDM.  

python usage : DDM_3D_ExSolver(storage, checktype, tol, relax, nb_iter_check,
nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";

%feature("docstring") DDM_3D_ComputeDof "

Compute degrees of freedom.  

python usage : DDM_3D_ComputeDof()  
";

%feature("docstring") DDM_3D_Post "

Does postpro operation related to DDM.  

python usage : DDM_3D_Post()  
";

%feature("docstring") DDM_3D_SetParameters "

Set frequencies parameters of DDM.  

python usage : DDM_3D_SetParameters(f_ddm, f_out, f_last, f_postpro, f_display)  

Parameters
----------
f_ddm(integer) : frequency of ddm partionning  
f_out(integer) : frequency of output file writing  
f_last(integer) : frequency of last file writing  
f_postpro(integer) : frequency of postpro file writing  
f_display(integer) : frequency of display file writing  
";

%feature("docstring") DDM_3D_WriteLast "

Write some data at the end of computation.  

python usage : DDM_3D_WriteLast()  
";

%feature("docstring") DDM_3D_Finalize "

End of computation/module's life.  

python usage : DDM_3D_Finalize()  
";


// File: wrap__utilities_8h.xml

%feature("docstring") utilities_logMes "

ask to write a message  

If the message is too long (more than 256 characters) it will be truncated  

python usage : utilities_logMes(message)  

Parameters
----------
message(string) : log message to add  
length(integer) : length of the message string  
";

%feature("docstring") utilities_DisableLogMes "

disable printing of messages  

python usage : utilities_DisableLogMes()  
";

%feature("docstring") utilities_EnableLogMes "

enable priting of messages  

python usage : utilities_EnableLogMes()  
";

%feature("docstring") utilities_setIoUnitLimits "

set the interval of unit numbers lmgc90 can use to open file  
";

%feature("docstring") utilities_setStopMode "

Decide to stop or store a message in case of fatal error.  

python usage : utilities_setStopMode()  
";

%feature("docstring") utilities_resetFatal "

Clean fatal error state.  

This function is not intended to be used in python but by swig to throw an
excpetion  
";

%feature("docstring") utilities_checkFatal "
";

%feature("docstring") utilities_OpenFileStandardOutput "

Select the file for standard and errors outputs.  

If the filename is too long (more than 256 characters) it will be truncated  

python usage : utilities_OpenFileStandardOutput(filename)  

Parameters
----------
filename(string) : the name of file  
length(integer) : length the name of the file  
";

%feature("docstring") utilities_CloseFileStandardOutput "

Close the file for standard and errors outputs.  

python usage : utilities_CloseFileStandardOutput()  
";

%feature("docstring") utilities_InitRandomSeed "

Re-initialize the seed of the build-in random function.  

python usage : utilities_InitRandomSeed([seed])  

Parameters
----------
seed(integer array) : an optional desired input seed  
";

%feature("docstring") utilities_Finalize "

End of simulation operations.  

Only close all possibly opened units by the program.  

python usage : utilities_Finalize()  
";


// File: wrap__SPPRx_8h.xml

%feature("docstring") SPPRx_SelectProxTactors "

contact detection between SPxxx and PLxxx tactors  

python usage : SPPRx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPPRx_SmoothForceComputation "

compute smooth contact law (in any)  

python usage : SPPRx_SmoothForceComputation()  
";

%feature("docstring") SPPRx_WriteLastVlocRloc "

write last local values of all SPPRx contacts  

python usage : SPPRx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_WriteOutVlocRloc "

write local values of all SPPRx contacts  

python usage : SPPRx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_DisplayOutVlocRloc "

display local values of all SPPRx contacts  

python usage : SPPRx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPPRx_DisplayProxTactors "

display contacts  

python usage : SPPRx_DisplayProxTactors()  
";

%feature("docstring") SPPRx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") SPPRx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPPRx_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(real) : period on x axis  
";

%feature("docstring") SPPRx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPPRx_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(real) : period on y axis  
";

%feature("docstring") SPPRx_CleanMemory "

Free all memory allocated within SPPRx module.  

python usage : SPPRx_CleanMemory()  
";


// File: wrap__PRPLx_8h.xml

%feature("docstring") PRPLx_SelectProxTactors "

contact detection between PRxxx and PLxxx tactors  

python usage : PRPLx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PRPLx_WriteLastVlocRloc "

write last local values of all PRPLx contacts  

python usage : PRPLx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_WriteOutVlocRloc "

write local values of all PRPLx contacts  

python usage : PRPLx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_DisplayOutVlocRloc "

display local values of all PRPLx contacts  

python usage : PRPLx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PRPLx_DisplayProxTactors "

display contacts  

python usage : PRPLx_DisplayProxTactors()  
";

%feature("docstring") PRPLx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PRPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PRPLx_CleanMemory "

Free all memory allocated within PRPLx module.  

python usage : PRPLx_CleanMemory()  
";


// File: wrap__mp__solver__3D_8h.xml

%feature("docstring") mp_solver_3D_ReadMpBehaviour "

python usage : mp_solver_3D_ReadMpBehaviour()  
";

%feature("docstring") mp_solver_3D_WriteMpBehaviour "

python usage : mp_solver_3D_WriteMpBehaviour()  
";

%feature("docstring") mp_solver_3D_ReadIniMpValues "

Read MP_VALUES file.  

If num <= 0 : DATBOX/MP_VALUES.INI file is read Else : OUTBOX/MP_VALUES.OUT.num
is read, num being the parameter used in TimeEvolution_ReadIniDof last call  

usage : mp_solver_3D_ReadIniMpValues(num=0)  

Parameters
----------
num(integer) : which file to read  
";

%feature("docstring") mp_solver_3D_WriteOutMpValues "

python usage : mp_solver_3D_WriteOutMpValues()  
";

%feature("docstring") mp_solver_3D_WriteLastMpValues "

python usage : mp_solver_3D_WriteLastMpValues()  
";

%feature("docstring") mp_solver_3D_SolveElectro1G "

python usage : mp_solver_3D_SolveElectro1G()  
";

%feature("docstring") mp_solver_3D_SolveNlElectro1G "

python usage : mp_solver_3D_SolveNlElectro1G()  
";

%feature("docstring") mp_solver_3D_SolveThermoProblem "

python usage : mp_solver_3D_SolveThermoProblem()  
";

%feature("docstring") mp_solver_3D_UpdateThermoProblem "

python usage : mp_solver_3D_UpdateThermoProblem()  
";

%feature("docstring") mp_solver_3D_RecupTemperature "

python usage : mp_solver_3D_RecupTemperature()  
";

%feature("docstring") mp_solver_3D_RecupPotential "

python usage : mp_solver_3D_RecupPotential()  
";

%feature("docstring") mp_solver_3D_UpdateConductivity "

python usage : mp_solver_3D_UpdateConductivity()  
";


// File: wrap__RBDY3_8h.xml

%feature("docstring") RBDY3_IncrementStep "

compute the current velocity and displacement  

python usage : RBDY3_IncrementStep()  
";

%feature("docstring") RBDY3_SetVlocyDrivenDof "

Override the value of an existing velocity driven dof.  

usage : RBDY3_SetVlocyDrivenDof(ibdyty, idrvdof, value)  

Parameters
----------
ibdyty(integer) : rank of considered  
idrvdof(integer) : index of velocity driven dof to set  
value(real) : new value of the velocity driven dof  
";

%feature("docstring") RBDY3_FatalDamping "

Nullify body velocities (current and initial) of a list of bodies.  

python usage : RBDY3_FatalDamping(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to reset current velocity if omitted
    works on all objetcs  
";

%feature("docstring") RBDY3_ComputeFext "

compute external forces  

python usage : RBDY3_ComputeFext()  
";

%feature("docstring") RBDY3_ComputeBulk "

compute internal forces  

python usage : RBDY3_ComputeBulk()  
";

%feature("docstring") RBDY3_ComputeFreeVelocity "

compute free velocity  

python usage : RBDY3_ComputeFreeVelocity()  
";

%feature("docstring") RBDY3_ComputeDof "

update current position and velocity  

python usage : RBDY3_ComputeDof()  
";

%feature("docstring") RBDY3_UpdateDof "

save d.o.f. of the end of the time step to d.o.f. of the begining of the next
one  

python usage : RBDY3_UpdateDof()  
";

%feature("docstring") RBDY3_ComputeContactDetectionConfiguration "

compute the contact detection configuration  

python usage : RBDY3_ComputeContactDetectionConfiguration()  
";

%feature("docstring") RBDY3_WriteLastDof "

write ascii DOF.LAST file  

python usage : RBDY3_WriteLastDof()  
";

%feature("docstring") RBDY3_WriteOutDof "

write ascii DOF.OUT file. Can be activate only each N step  

If 0 for ifrom and ito, dofs of all bodies are written.  

python usage : RBDY3_WriteOutDof(ifrom=0, ito=0)  

Parameters
----------
ifrom(integer) : begining of bodys' index that will be written  
ito(integer) : end of bodys'index that will be written  
";

%feature("docstring") RBDY3_DisplayOutDof "

display bodies degrees of freedom  

python usage : RBDY3_DisplayOutDof()  
";

%feature("docstring") RBDY3_WriteLastRnod "

write ascii Rnod.LAST file  

python usage : RBDY3_WriteLastRnod()  
";

%feature("docstring") RBDY3_WriteOutRnod "

write ascii Rnod.OUT file. Can be activate only each N step.  

python usage : RBDY3_WriteOutRnod()  
";

%feature("docstring") RBDY3_DisplayOutRnod "

display body forces.  

python usage : RBDY3_DisplayOutRnod()  
";

%feature("docstring") RBDY3_WriteBodies "

write BODIES.OUT file  

python usage : RBDY3_WriteBodies()  
";

%feature("docstring") RBDY3_WriteDrivenDof "

write DRV_DOF.OUT file  

python usage : RBDY3_WriteDrivenDof()  
";

%feature("docstring") RBDY3_ReadBodies "

read BODIES.DAT file  

Initializes existing_entities variable in RBDY3  

Adds the number of found bodies to entity  

python usage : RBDY3_ReadBodies()  
";

%feature("docstring") RBDY3_ReadCompressedBodies "

read BODIES.DAT file without any comment  

Initializes existing_entities variable in RBDY3  

Adds the number of found bodies to entity  

python usage : RBDY3_ReadCompressedBodies()  
";

%feature("docstring") RBDY3_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

usage : RBDY3_ReadIniDof(num=0)  

Parameters
----------
num(integer) : which DOF file to read  
";

%feature("docstring") RBDY3_ReadDrivenDof "

read DRV_DOF.DAT file  

python usage : RBDY3_ReadDrivenDof()  
";

%feature("docstring") RBDY3_LoadBehaviours "

Load bulk behaviour id from bulk_behav module.  

python usage : RBDY3_LoadBehaviours()  
";

%feature("docstring") RBDY3_ComputeMass "

compute mass and inertia of bodies  

python usage : RBDY3_ComputeMass()  
";

%feature("docstring") RBDY3_NewRotationScheme "

active new rotation scheme FLAG  

python usage : RBDY3_NewRotationScheme()  
";

%feature("docstring") RBDY3_SetZminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetZminBoundary(Zmin)  

Parameters
----------
Zmin(real) : inferior boundary value  
";

%feature("docstring") RBDY3_SetZmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetZmaxBoundary(Zmax)  

Parameters
----------
Zmax(real) : superior boundary value  
";

%feature("docstring") RBDY3_SetYminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetYminBoundary(Ymin)  

Parameters
----------
Ymin(real) : left boundary value  
";

%feature("docstring") RBDY3_SetYmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetYmaxBoundary(Ymax)  

Parameters
----------
Ymax(real) : right boundary value  
";

%feature("docstring") RBDY3_SetXminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetXminBoundary(Xmin)  

Parameters
----------
Xmin(real) : inferior boundary value  
";

%feature("docstring") RBDY3_SetXmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetXmaxBoundary(Xmax)  

Parameters
----------
Xmax(real) : front boundary value  
";

%feature("docstring") RBDY3_SetXPeriodicCondition "

set the period on X axis  

python usage : RBDY3_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(real) : period on x axis  
";

%feature("docstring") RBDY3_SetYPeriodicCondition "

set the periode on Y axis  

python usage : RBDY3_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(real) : period on y axis  
";

%feature("docstring") RBDY3_AvoidBodyRotation "

kill rotation effect for RBDY3  

python usage : RBDY3_AvoidBodyRotation()  
";

%feature("docstring") RBDY3_SkipInvisible "

if a body is invisible, il will not be written in bodies.out and dof.out  

python usage : RBDY3_SkipInvisible()  
";

%feature("docstring") RBDY3_KeepIniDofOrder "

numbering information as they are read  

python usage : RBDY3_KeepIniDofOrder()  
";

%feature("docstring") RBDY3_SetVisible "

rended a given RBDY3 visible  

python usage : RBDY3_SetVisible(ibdyty)  

Parameters
----------
ibdyty(integer): index of the RBDY3  
";

%feature("docstring") RBDY3_SetInvisible "

rended a given RBDY3 invisible  

python usage : RBDY3_SetInvisible(ibdyty)  

Parameters
----------
ibdyty(integer): index of the RBDY3  
";

%feature("docstring") RBDY3_IsVisible "

return if a given body visible  

python usage : visible = RBDY3_IsVisible(ibdyty)  

Parameters
----------
idbdy(integer): id of the body we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") RBDY3_CompCoor "

Compute the position of bodies.  

python usage : RBDY3_CompCoor()  
";

%feature("docstring") RBDY3_GetBodyDensity "

Get the density of a given body.  

python usage : density = RBDY3_GetBodyDensity(ibdyty)  

Parameters
----------
ibdyty(integer): rank of the RBDY3  

Returns
-------
density(double) : density of the RBDY3  
";

%feature("docstring") RBDY3_GetBodyInertia "

Get the principal inertia of a given RBDY3.  

python usage : inertia = RBDY3_GetBodyInertia(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3  

Returns
-------
inertia (double array) : inertia vector of the desired RBDY3  
";

%feature("docstring") RBDY3_GetAllInertia "

Get the inertia of a all RBDY3 body.  

usage : inertia = RBDY3_GetAllInertia()  

Parameters
----------
inertia(double array): the inertia of all bodies  
";

%feature("docstring") RBDY3_CollectBodiesDotOUT "

python usage : RBDY3_CollectBodiesDotOUT()  
";

%feature("docstring") RBDY3_AppendToBodiesDotOUT "

python usage : RBDY3_AppendToBodiesDotOUT()  
";

%feature("docstring") RBDY3_RebuildBodiesDotDAT "

python usage : RBDY3_RebuildBodiesDotDAT()  
";

%feature("docstring") RBDY3_PutBodyVector "

Set a vector of a RBDY3 body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
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

uses copy, and in case of Fext, the operation is not just setting but adding  

python usage : RBDY3_PutBodyVector(datatype, ibdyty, vector)  

Parameters
----------
datatype(string [5]) : the vector to set  
ibdyty(integer) : rank of the RBDY3  
vector(double array) : the new value of the vector  
";

%feature("docstring") RBDY3_PutAllBodyVector "

Put an array of a vector of all RBDY3 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY3_PutBodyVector  

python usage : RBDY3_PutAllBodyVector(datatype, matrix)  

Parameters
----------
datatype(string [5]) : the vector to set  
matrix(double array) : input matrix  
";

%feature("docstring") RBDY3_GetBodyVector "

Get a copy of a vector of a RBDY3 body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coor_\": coordinates in computed configuration  
*   \"Coorb\": coordinates at beginning of time step  
*   \"Coorm\": coordinates in detection configuration  
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

python usage : vector = RBDY3_GetBodyVector(datatype, ibdyty)  

Parameters
----------
datatype(string [5]) : the vector to get  
ibdyty(integer) : rank of the RBDY3  

Returns
-------
vector (double array) : output vector  
";

%feature("docstring") RBDY3_GetAllBodyVector "

Get an array of a vector of all RBDY3 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY3_GetBodyVector  

python usage : matrix = RBDY3_GetBodyVector(datatype, ibdyty)  

Parameters
----------
datatype(string [5]) : the vector to get  

Returns
-------
matrix (double array) : output matrix  
";

%feature("docstring") RBDY3_GetPtrBodyVector "

Get a pointer on a vector of a RBDY3 body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"Vaux_\": working array for velocity  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  

python usage : vector_ptr = RBDY3_GetPtrBodyVector(datatype, ibdyty)  

Parameters
----------
datatype(string [5]) : the vector to set  
ibdyty(integer) : rank of the RBDY3  

Returns
-------
vector_ptr (double array) : reference on the desired vector seen as a numpy
array  
";

%feature("docstring") RBDY3_PutBodyMatrix "

Set a matrix of a RBDY3 body.  

Possible values for datatype field are:  

*   \"IFbeg\": inertia frame at beginning of time step  
*   \"IFTT_\": inertia frame in detection configuration  
*   \"IF___\": inertia frame in computed configuration  

Uses copy  

python usage : RBDY3_PutBodyMatrix(datatype, ibdyty, matrix)  

Parameters
----------
datatype(string [5]) : the vector to set  
ibdyty(integer) : rank of the RBDY3  
matrix(double array) : a matrix  
";

%feature("docstring") RBDY3_GetBodyMatrix "

Get a copy of a matrix of a RBDY3 body.  

Possible values for datatype field are:  

*   \"IFref\": inertia frame in reference configuration  
*   \"IFbeg\": inertia frame at beginning of time step  
*   \"IFTT_\": inertia frame in detection configuration  
*   \"IF___\": inertia frame in computed configuration  

Uses copy  

python usage : matrix = RBDY3_GetBodyMatrix(datatype, ibdyty)  

Parameters
----------
datatype(string [5]) : the vector to get  
ibdyty(integer) : rank of the RBDY3  

Returns
-------
matrix (double array) : output matrix  
";

%feature("docstring") RBDY3_GetAllRData "

Get a copy of a real data of all rbdy3.  

In this order : coor, frame, vlocy, spin, fext, reac  

python usage : rdata = RBDY3_GetAllRData()  

Returns
-------
rdata (double array) : output matrix  
";

%feature("docstring") RBDY3_GetNbRBDY3 "

get the number of RBDY3  

python usage : nb_RBDY3 = RBDY3_GetNbRBDY3()  

Returns
-------
nb_RBDY3 (integer) : number of RBDY3 in container  
";

%feature("docstring") RBDY3_GetMass "

Get the mass of a body.  

python usage : mass = RBDY3_GetMass(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3  

Returns
-------
mass (double) : mass of the RBDY3  
";

%feature("docstring") RBDY3_GetAllMass "

Get the mass of a all body (visible and invisible)  

python usage : masses = RBDY3_GetAllMass()  

Returns
-------
masses (double array) : masses of all RBDY3  
";

%feature("docstring") RBDY3_GetPtrMass "

Get a pointer onto the mass matrix of a body.  

Parameters
----------
ibdyty(int): index of the RBDY3  
mass(double**): mass matrix of the RBDY3  
";

%feature("docstring") RBDY3_GetVelocity "

Get the velocity of a body.  

Parameters
----------
ibdyty(int): index of the RBDY3  
velocity(double[6]): velocity of the RBDY3  
";

%feature("docstring") RBDY3_GetGlobInertia "

Get the global inertia.  

usage : inertia = RBDY3_GetGlobInertia(ibdyty)  

Parameters
----------
ibdyty(integer) : id of desired RBDY3  

Returns
-------
inertia (double 2D array) : the inertia matrix  
";

%feature("docstring") RBDY3_GetBehavior "

Get the type of the nickname of the behavior.  

usage name = RBDY3_GetBehavior(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3 in container  

Returns
-------
type (string) : nickname  
";

%feature("docstring") RBDY3_GetNbContactor "

get the number of contactor of RBDY3  

python usage : nb = RBDY3_GetNbContactor(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3 in container  

Returns
-------
nb (integer) : number of contactor attached to a RBDY3  
";

%feature("docstring") RBDY3_GetContactorType "

Get the type of the itacty contactor of a body ibdyty.  

usage type = RBDY3_GetContactorType(ibdyty,itacty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3 in container  
itacty(integer) : rank of the contactor in the RBDY3  

Returns
-------
type (string) : type of the contactor of the body  
";

%feature("docstring") RBDY3_SetContactorColor "

Set the color of a given contactor of a body.  

usage : RBDY3_SetContactorColor(ibdyty, itacty, color)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3  
itacty(integer) : rank of the contactor in the RBDY3  
color(string of size 5) : the color  
";

%feature("docstring") RBDY3_GetContactorColor "

Get the color of the itacty contactor of a body ibdyty.  

usage color = RBDY3_GetContactorColor(ibdyty,itacty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY3 in container  
itacty(integer) : rank of the contactor in the RBDY3  

Returns
-------
color (string) : color of the contactor of the body  
";

%feature("docstring") RBDY3_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = RBDY3_getDrvVlocy(ibdyty)  

Parameters
----------
ibdyty(integer) : index of the RBDY3  
drvdof_indices(integer array) : indices list of driven dof  
drvdof_values(real array) : values of the driven dof  
";

%feature("docstring") RBDY3_computeDrvVlocy "

Compute the value of the driven velocity of a body at current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : RBDY3_computeDrvVlocy(ibdyty, values)  

Parameters
----------
ibdyty(integer) : index of the RBDY3  
values(double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") RBDY3_WriteOutOneBody "

write a bdyty to BODIES.OUT with a given rank  

python usage : RBDY3_WriteOutOneBody(ibdyty, new_ibdyty)  

Parameters
----------
ibdyty(integer) : index of the RBDY3  
new_ibdyty(integer): new index of the RBDY3  
";

%feature("docstring") RBDY3_WriteOutDofOneBody "

write a bdyty dof to DOF.OUT with a given rank  

python usage : RBDY3_WriteOutDofOneBody(ibdyty, new_ibdyty)  

Parameters
----------
ibdyty(integer) : index of the RBDY3  
new_ibdyty(integer): new index of the RBDY3  
";

%feature("docstring") RBDY3_LoadThreadNetwork "

read thread structure for textile structure  

python usage : RBDY3_LoadThreadNetwork(void);  
";

%feature("docstring") RBDY3_SetInvisibleSmallObjects "

Set the objects to invisible if their average radius is less than radius.  

python usage : RBDY3_SetInvisibleSmallObjects(radius)  

Parameters
----------
radius(double) : radius threshold  
";

%feature("docstring") RBDY3_SetVisibleVlocyDrivenDof "

rended a given Velocy DOF visible  

python usage : RBDY3_SetVisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
ibdyty(integer): index of the RBDY3  
iccdof(integer): index of the DOF to set visible  
";

%feature("docstring") RBDY3_SetInvisibleVlocyDrivenDof "

rended a given Velocy DOF invisible  

python usage : RBDY3_SetInvisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
ibdyty(integer): index of the RBDY3  
iccdof(integer): index of the DOF to set invisible  
";

%feature("docstring") RBDY3_PartialDamping "

Limit body velocity to Vmax value.  

usage : RBDY3_PartialDamping(nb_steps, Vmax)  

Parameters
----------
nb_steps(integer) : periodicity @parma[in] Vmax (double) : Vmax  
";

%feature("docstring") RBDY3_GetVolume "

Get volume of a body.  

usage : volume = RBDY3_GetVolume(ibdyty)  

Parameters
----------
ibdyty(integer) : RBDY3 id  

Returns
-------
volume (double) : volume  
";

%feature("docstring") RBDY3_GetAllVolume "

Get the area of a all body (visible and invisible)  

python usage : area = RBDY3_GetAllVolume()  

Returns
-------
area (double array) : masses of all RBDY2  
";

%feature("docstring") RBDY3_RenumVisibleBodies "

give a new numerotation of visible bodies  

python usage : RBDY3_RenumVisibleBodies()  
";

%feature("docstring") RBDY3_GetBulkBehavNumber "

return the bulk number of a given RBDY3  

python usage : ibehav = RBDY3_GetBulkBehavNumber(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of a RBDY3  

Returns
-------
ibehav (integer) : the bulk behav number  
";

%feature("docstring") RBDY3_CleanMemory "

Free all memory allocated within RBDY3 module.  

python usage : RBDY3_CleanMemory()  
";

%feature("docstring") RBDY3_LoadMpBehaviours "

read extra physical behaviour in BULK_BEHAV.DAT file.  

Must be used for THERMO_RIGID ELECTRO_RIGID and THERMO_ELECTRO_RIGID behaviour  

python usage : RBDY3_LoadMpBehaviours(disper)  

Parameters
----------
disper(double): some dispersion coefficient  
";

%feature("docstring") RBDY3_IncrementWSvsT "

python usage : RBDY3_IncrementWSvsT()  
";

%feature("docstring") RBDY3_UpdateGAMMAvsT "

python usage : RBDY3_UpdateGAMMAvsT()  
";

%feature("docstring") RBDY3_GetThermalValue "

Get temperature of rigid particle.  

usage : T = RBDY3_GetThermalValu(ibdyty, itacty)  

Parameters
----------
ibdyty(integer) : rank of body  
itacty(integer) : rank of tacty \"  
";

%feature("docstring") RBDY3_SetEquilibriumNorm "

Initialization of data for the equilibrium state check.  

You must precise the type of check test :  

*   Qvlcy : quadratic norm velocy  
*   Mvlcy : maximum norm velocy  

usage : RBDY3_CheckEquilibrium(norm_type , tolerance)  

Parameters
----------
norm_type(string of size 5) : norm type use for the equilibrium check  
tolerance(double) : norm tolerance  
";

%feature("docstring") RBDY3_CheckEquilibriumState "

check if all the RBDY3 rich an equilibrium state (velocity is almost equal to
zero)  

usage : isBalanced = RBDY3_CheckEquilibriumState()  

Returns
-------
isBalanced (boolean) : True if in equilibrium state  
";

%feature("docstring") RBDY3_SetSourcePoint "

create an assembly by source point deposit  

python usage : RBDY3_SetSourcePoint(first_RBDY3, radius, Xshift, Yshift, Zshift)  

Parameters
----------
first_RBDY3(int): number of first invisible body  
radius: source point area radius  
Xshift: X translation of deposited object from reference coordinate  
Yshift: Y translation of deposited object from reference coordinate  
Zshift: Z translation of deposited object from reference coordinate  
";

%feature("docstring") RBDY3_SetSourcePointWithIni "

create an assembly by source point deposit  

python usage : RBDY3_SetSourcePointWithIni(first_RBDY3, radius, Xshift, Yshift,
Zshift)  

Parameters
----------
first_RBDY3(int): number of first invisible body  
radius: source point area radius  
Xshift: X coordinate of deposited object  
Yshift: Y coordinate of deposited object  
Zshift: Z coordinate of deposited object  
";

%feature("docstring") RBDY3_InitializeProgressiveActivation "

set the progression of altitude  

python usage : RBDY3_InitializeProgressiveActivation(zini, dz)  

Parameters
----------
zini(real) : initial altitude  
dz(real) : increment of altitude  
";

%feature("docstring") RBDY3_ApplyProgressiveActivation "

set occurence of activation  

python usage : RBDY3_ApplyProgressiveActivation(freq)  

Parameters
----------
freq(integer) : activation frequence of progression  
";

%feature("docstring") RBDY3_InitFreeBoundary "

python usage : RBDY3_InitFreeBoundary(xmin, xmax, ymin, ymax, radius)  

Parameters
----------
xmin(real) :  
xmax(real) :  
ymin(real) :  
ymax(real) :  
radius(real) :  
";

%feature("docstring") RBDY3_TriaxialLoading "

Triaxial load of a sample using a rigid box.  

python usage : TriaxialLoading(num_down, num_right, num_up, num_left, num_front,
num_rear, nb_loads, loads)  

Parameters
----------
num_down(integer) :  
num_right(integer) :  
num_up(integer) :  
num_left(integer) :  
num_front(integer) :  
num_rear(integer) :  
nb_loads(integer) : the number of walls you want to load with a pressure (1 to 6)  
loads(array) : loads(2,nb_loads): load(1,i) contains which wall is loaded
    (1==down, 2==right, 3==up, 4==left, 5==front, 6==rear) and load(2,i)
    contains the amplitude of the stress (a positive value means compression).  
";

%feature("docstring") RBDY3_GetDofStatus "

Get dof status.  

python usage : status = RBDY3_GetDofStatus(ibdyty)  

Parameters
----------
ibdyty(integer): rank of the RBDY3  

Returns
-------
status(integer) : dof status of the RBDY3  
";


// File: wrap__nlgs__3D_8h.xml

%feature("docstring") nlgs_3D_ExIter "

Executes nb_iter NLGS iterations.  

python usage : nlgs_3D_ExIter(nb_iter) param[in] nb_iter (integer) : number of
iterations to do  
";

%feature("docstring") nlgs_3D_ExIterJacobi "

Executes nb_iter NLJacobi iterations.  

python usage : nlgs_3D_ExIterJacobi(nb_iter) param[in] nb_iter (integer) :
number of iterations to do  
";

%feature("docstring") nlgs_3D_AfterIterCheck "

Control NLGS convergence.  

python usage : convergence = nlgs_3D_AfterIterCheck()  

Returns
-------
convergence (integer) :  
";

%feature("docstring") nlgs_3D_AfterIterCheckJacobi "

Control NLGS convergence.  

python usage : convergence = nlgs_3D_AfterIterCheckJacobi()  

Returns
-------
convergence (integer) :  
";

%feature("docstring") nlgs_3D_ScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_3D_ScrambleContactOrder()  
";

%feature("docstring") nlgs_3D_QuickScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_3D_QuickScrambleContactOrder()  
";

%feature("docstring") nlgs_3D_ReverseContactOrder "

reverse the numbering of the contact list  

python usage : nlgs_3D_ReverseContactOrder()  
";

%feature("docstring") nlgs_3D_DisplayAfterIterCheck "

display NLGS convergence results  

python usage : nlgs_3D_DisplayAfterIterCheck()  
";

%feature("docstring") nlgs_3D_ScaleRloc "

scale all local contact forces of a factor equal to 0.9 < f < 1.1  

python usage : nlgs_3D_ScaleRloc()  
";

%feature("docstring") nlgs_3D_ComputeRnod "

mapping from local contact forces to global ones  

python usage : nlgs_3D_ComputeRnod()  
";

%feature("docstring") nlgs_3D_ExPost "

run a jacobi iteration with the solution obtain with the NLGS algorithm  

python usage : nlgs_3D_ExPost()  
";

%feature("docstring") nlgs_3D_ExPostJacobi "

run a jacobi iteration with the solution obtain with the NLGS algorithm  

python usage : nlgs_3D_ExPostJacobi()  
";

%feature("docstring") nlgs_3D_SetCheckType "

define numerical convergence of the NLGS algorithm  

python usage : nlgs_SetCheckType(check_type, tolerance, relaxation)  

Parameters
----------
chekctype_c(char[5]) : type of convergence check  
tol(double) : norm tolerance  
relax(double) : relaxation factor  
  
 convergence check keywords:  
 Quad : quadratic norm (faulty contacts are redeemed by accurate contacts;
laxist norm)  
 Maxm : maximum norm (faulty contacts must comply; severe norm)  
 QM/16 : maximum of Quad and Maxm/16 norms (a compromise). For large dense
collections Quad ranges usually around 1/16 Maxm  
 where Quad,Maxm,QM/16 are keywords for the check test, and the following real
number is the tolerance value.  
";

%feature("docstring") nlgs_3D_ExPrep "

Prepare matrix storage.  

python usage : nlgs_ExPrep(storage)  

Parameters
----------
storage_c(char[30]): matrix storage  
  
 prepare the matrix and the RHS of the contact problem in regards of the
selected matrix storage:  

*   Exchange_Local_Global (the standard case) only the diagonal blocks are
    computed and stored.  
*   Stored_Delassus_Loops (faster but memory expensive) the complete Delassus
    matrix is computed.  
";

%feature("docstring") nlgs_3D_WriteNormCheck "

write norm to file  

python usage : nlgs_3D_WriteNormCheck()  
";

%feature("docstring") nlgs_3D_DiagonalResolution "

python usage : nlgs_3D_DiagonalResolution()  
";

%feature("docstring") nlgs_3D_SetWithQuickScramble "

Activate quick scramble in macro function ExSolver.  

python usage : nlgs_3D_SetWithQuickScramble()  
";

%feature("docstring") nlgs_3D_SetWithReverseContactOrder "

Activate reverse order in macro function ExSolver.  

python usage : nlgs_3D_SetWithReverseContactOrder()  
";

%feature("docstring") nlgs_3D_UseJacobiSolver "

Use a Jacobi solver instead of Gauss Seidel solver.  

usage : nlgs_3D_UseJacobiSolver(True) or nlgs_UseJacobiSolver(False)  
";

%feature("docstring") nlgs_3D_ExSolver "

Solve fully the local contact problem.  

python usage : nlgs_3D_ExSolver(storage, checktype, tol, relax, nb_iter_check,
nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";

%feature("docstring") nlgs_3D_UpdateTactBehav "

update internal parameters of contact laws for each contact  

python usage : nlgs_3D_UpdateTactBehav()  
";

%feature("docstring") nlgs_3D_IsInitialized "

In case of restart say that nlgs is initialized.  

python usage : nlgs_3D_IsInitialized(is_init=1)  
";

%feature("docstring") nlgs_3D_DisplayTacInfo "

Display information concerning one contact.  

python usage : nlgs_3D_DsplayTacInfo(itac) param[in] itac (integer) : contact
rank  
";

%feature("docstring") nlgs_3D_UseRegularization "

use some regularization heuristics on interaction laws  

python usage : nlgs_3D_UseRegularization(krn, krt)  

Parameters
----------
krn(double) : normal penality (default 1e14)  
krt(double) : tangential penality (default 1e14)  
";

%feature("docstring") nlgs_3D_CutOpenCZM "

If some czm contact have a gap greater than the given they are considered as
broken ; works only with EXPO_CZM or IQS_EXPO_CZM.  

python usage : nlgs_3D_CutOpenCZM(tol)  

Parameters
----------
tol(double) : threshold on positive distance (default 1e-6)  
";

%feature("docstring") nlgs_3D_ManageInterpenetratedCZM "

Apply a g0 strategy if gap is negative and if gap is positive (without using
nlgs_3D_CutOpenCZM) ; works only with EXPO_CZM or IQS_EXPO_CZM.  

python usage : nlgs_3D_ManageInterpenetratedCZM()  
";


// File: wrap__tact__behav_8h.xml

%feature("docstring") tact_behav_OpenBehavContainer "

open the container (access as a linked list) in order to add/remove objects  

python usage : tact_behav_OpenBehavContainer()  
";

%feature("docstring") tact_behav_CloseBehavContainer "

close the container (access as an array)  

python usage : tact_behav_TactBehavContainer()  
";

%feature("docstring") tact_behav_OpenSeeContainer "

open the container (access as a linked list) in order to add/remove objects  

python usage : tact_behav_OpenSeeContainer()  
";

%feature("docstring") tact_behav_CloseSeeContainer "

close the container (access as an array)  

python usage : tact_behav_CloseSeeContainer()  
";

%feature("docstring") tact_behav_FillContainersFromFile "

read DATBOX/TACT_BEHAV.DAT and fill the containers (see and tact)  

python usage : tact_behav_FillContainersFromFile()  
";

%feature("docstring") tact_behav_AddToSeeContainer "

add a see table to the container  

python usage :
tact_behav_AddToSeeContainer(cdbdy,cdtac,cdcol,behav,anbdy,antac,ancol,alert,global_alert)  
";

%feature("docstring") tact_behav_ReadBehaviours "

open + fill + close  

python usage : tact_behav_ReadBehaviours()  
";

%feature("docstring") tact_behav_CollectOutTactBehav "

old fashion read from OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_CollectOutTactBehav()  
";

%feature("docstring") tact_behav_WriteBehaviours "

write (replace) tact and see to OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_WriteBehaviours()  
";

%feature("docstring") tact_behav_AppendOutTactBehav "

write (append) tact and see to OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_AppendOutTactBehav()  
";

%feature("docstring") tact_behav_RebuildInTactBehav "

write (replace) tact and see to DATBOX/TACT_BEHAV.DAT  

python usage : tact_behav_RebuildInTactBehav()  
";

%feature("docstring") tact_behav_CleanOutTactBehav "

erase OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_CleanOutTactBehav()  
";

%feature("docstring") tact_behav_GetNbTactBehav "

get the number of tact laws  

python usage : nb_tact_behav = tact_behav_GetNbTactBehav()  

Parameters
----------
nb_tact_behav(integer) : number of contact behaviour in lmgc90  
";

%feature("docstring") tact_behav_GetTactBehav "

get information related to a given tact law  

python usage : [lawty, behav, param] = tact_behav_GetTactBehav(i_tb)  

Parameters
----------
i_tb(integer) : rank (in the contact laws list) of the desired tact_behav  
lawty(string) : type of the contact law  
behav(string) : name of the contact law  
param(real vector) : parameters of the law  
";

%feature("docstring") tact_behav_GetInternalComment "

Get internal variables comment of a given interaction law.  

python usage : comment = tact_behav_GetInternalComment(ilaw)  

Parameters
----------
ilaw(integer) : rank of the interaction law  

Returns
-------
comment (char[100]) : the string to get  
";

%feature("docstring") tact_behav_SetCZMwithInitialFriction "

define the way friction evolve with damage: =0. constant value, (1. - beta)**pow
otherwize  

python usage : tact_behav_SetCZMwithInitialFriction(pow)  

Parameters
----------
pow(real) : parameter of power law evlution for friction
    mu(beta)=mu_s*(1-beta)**pow  
";

%feature("docstring") tact_behav_initFrictionEvolution "

[experimental] read a friction time evolution map  

python usage : tact_behav_initFrictionEvolution()  
";

%feature("docstring") tact_behav_setRandomFriction "

Active variation of local friction.  

python usage : tact_behav_setRandomFriction(r8)  
";

%feature("docstring") tact_behav_GetTactBehavRankFromName "

get the rank (in the list of tact laws) of a tact behav law  

python usage : rank = tact_behav_GetTactBehavRankFromName(c5)  
";

%feature("docstring") tact_behav_GetParamRankFromName "

get the rank of a param for a given tact behav law  

python usage : rank = tact_behav_GetParamRankFromName(i_tact,c5)  
";

%feature("docstring") tact_behav_GetParam "

get the value of a parameter  

python usage : param = tact_behav_GetParam(i_tact,i_param)  

Parameters
----------
i_tact(integer) : rank of the interaction law  
i_param(integer) : rank of the parameter  
param(real ) : value of the parameter  
";

%feature("docstring") tact_behav_SetParam "

set the value ...  

python usage : tact_behav_SetParam(i_tact, i_param, param)  

Parameters
----------
i_tact(integer) : rank of the interaction law  
i_param(integer) : rank of the parameter  
param(real ) : value of the parameter  
";

%feature("docstring") tact_behav_GetLawInternalComment "
";

%feature("docstring") tact_behav_SetRNcap "

set a maximal compression value  

python usage : tact_behav_SetRNcap(param)  
";

%feature("docstring") tact_behav_SetDilatancyParameters "

set dilatancy parameters  

python usage : tact_behav_SetDilatancyParameters(fric,height)  
";

%feature("docstring") tact_behav_SetPressureParameters "

set pressure parameters  

Parameters
----------
ibehav(integer) : rank of the tact behav  
flag(integer) : kind of build-in pressure law (0 no pressure, 1: time dependent,
    2: linearly progressive since crack starts, 3: exponentially progressive
    since crack starts, 4: external)  
params(double array) : the new value of the params [p0,dp,tau,alpha]  

python usage : tact_behav_SetPressureParameters(ibehav,flag,params)  
";

%feature("docstring") tact_behav_CleanMemory "

Free all memory allocated within tact_behav module.  

python usage : tact_behav_CleanMemory()  
";


// File: wrap__cpg_8h.xml

%feature("docstring") cpg_ExIter "

Execute one CPG iteration over the contact loop.  

python usage cpg_ExIter()  
";

%feature("docstring") cpg_AfterIterCheck "

Control CPG convergence.  

python usage cpg_AfterIterCheck()  
";

%feature("docstring") cpg_ExPost "

Transfer local solution.  

python usage cpg_ExPost()  
";

%feature("docstring") cpg_ExPrep "

prepare the matrix and the RHS of the contact problem  

python usage cpg_ExPrep()  
";

%feature("docstring") cpg_ScaleRloc "

scale all local contact forces of a factor equal to 0.9 < f < 1.1  

python usage cpg_ScaleRloc()  
";

%feature("docstring") cpg_SetDiagonalPrecond "

active diagonal preconditioner  

python usage cpg_SetDiagonalPrecond()  
";

%feature("docstring") cpg_SetFrictionless "

active frictionless solver  

python usage cpg_SetFrictionless()  
";

%feature("docstring") cpg_SetNoConjugaison "

desactive conjugaison  

python usage cpg_SetNoConjugaison()  
";

%feature("docstring") cpg_SetCheckType "

define numerical convergence of the NLGS algorithm  

python usage cpg_SetCheckType(checktype, tol)  

Parameters
----------
chekctype(char[5]) : type of convergence check  
tol(double) : norm tolerance  
  
 convergence check keywords:  
 Quad : quadratic norm (faulty contacts are redeemed by accurate contacts;
laxist norm)  
 Maxm : maximum norm (faulty contacts must comply; severe norm)  
 QM/16 : maximum of Quad and Maxm/16 norms (a compromise). For large dense
collections Quad ranges usually around 1/16 Maxm  
 where Quad,Maxm,QM/16 are keywords for the check test, and the following real
number is the tolerance value.  
";

%feature("docstring") cpg_NormCheck "

Active one step norm evolution.  

python usage : cpg_norm_check()  
";

%feature("docstring") cpg_ExSolver "

Solve fully the local contact problem.  

python usage : cpg_ExSolver(checktype, tol, nb_iter_check, nb_block_iter)  

Parameters
----------
checktype(char[5]) c : convergentce test keyword  
tol(double) : tolerance value  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";


// File: wrap__SPCDx_8h.xml

%feature("docstring") SPCDx_SelectProxTactors "

contact detection between SPxxx and CDxxx tactors  

python usage : SPCDx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPCDx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : SPCDx_SmoothForceComputation()  
";

%feature("docstring") SPCDx_WriteLastVlocRloc "

write last local values of all SPCDx contacts  

python usage : SPCDx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_WriteOutVlocRloc "

write local values of all SPCDx contacts  

python usage : SPCDx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_DisplayOutVlocRloc "

display local values of all SPCDx contacts  

python usage : SPCDx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPCDx_DisplayProxTactors "

display contacts  

python usage : SPCDx_DisplayProxTactors()  
";

%feature("docstring") SPCDx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPCDx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") SPCDx_CleanMemory "

Free all memory allocated within SPCDx module.  

python usage : SPCDx_CleanMemory()  
";


// File: wrap__SiconosNumerics_8h.xml

%feature("docstring") SiconosNumerics_SetParameters "
";

%feature("docstring") SiconosNumerics_ExSolver "

Solve fully the local contact problem.  

python usage : SiconosNumerics_ExSolver()  
";

%feature("docstring") SiconosNumerics_IsInitialized "

In case of restart say that nlgs is initialized.  

python usage : SiconosNumerics_IsInitialized()  
";


// File: wrap__PLPLx_8h.xml

%feature("docstring") PLPLx_SelectProxTactors "

contact detection between POLYG tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLPLx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
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
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PLPLx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : PLPLx_SetPeriodicCondition(period)  

Parameters
----------
period(double) : value of the period  
";

%feature("docstring") PLPLx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : PLPLx_SetFrictionModel(cflag)  

Parameters
----------
cflag(char) : model to use ('min', 'max' or 'ave')  
";

%feature("docstring") PLPLx_SetBigPolygTolerance "

python usage : PLPLx_SetBigPolygTolerance(tol)  

Parameters
----------
period(double) : value of the tolerance  
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
icdan(int): index of the PLPLx contact  

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
shrink(real) :  
";


// File: wrap__DKPLx_8h.xml

%feature("docstring") DKPLx_SelectProxTactors "

contact detection between DISKx and POLYG tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKPLx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKPLx_WriteLastVlocRloc "

write last local values of all DKPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKPLx_WriteLastVlocRloc()  
";

%feature("docstring") DKPLx_WriteOutVlocRloc "

write local values of all DKPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKPLx_WriteOutVlocRloc()  
";

%feature("docstring") DKPLx_DisplayOutVlocRloc "

display local values of all DKPLx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKPLx_DisplayOutVlocRloc()  
";

%feature("docstring") DKPLx_DisplayProxTactors "

display contacts  

python usage : DKPLx_DisplayProxTactors()  
";

%feature("docstring") DKPLx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKPLx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKPLx_SetPeriodicCondition(period)  

Parameters
----------
period(double) : value of the period  
";

%feature("docstring") DKPLx_CleanMemory "

Free all memory allocated within DKPLx module.  

python usage : DKPLx_CleanMemory()  
";


// File: wrap__P2P2L_8h.xml

%feature("docstring") P2P2L_SelectProxTactors "

contact detection between PT2DL tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : P2P2L_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") P2P2L_WriteLastVlocRloc "

write last local values of all P2P2L contacts  

The values written are relative velocity, forces and local frame  

python usage : P2P2L_WriteLastVlocRloc()  
";

%feature("docstring") P2P2L_WriteOutVlocRloc "

write local values of all P2P2L contacts  

The values written are relative velocity, forces and local frame  

python usage : P2P2L_WriteOutVlocRloc()  
";

%feature("docstring") P2P2L_DisplayOutVlocRloc "

display local values of all P2P2L contacts  

The values displayed are relative velocity, forces and local frame  

python usage : P2P2L_DisplayOutVlocRloc()  
";

%feature("docstring") P2P2L_DisplayProxTactors "

display contacts  

python usage : P2P2L_DisplayProxTactors()  
";

%feature("docstring") P2P2L_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : P2P2L_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") P2P2L_CleanMemory "

Free all memory allocated within P2P2L module.  

python usage : P2P2L_CleanMemory()  
";


// File: wrap__DDM__ExternalFEM_8h.xml

%feature("docstring") DDM_ExternalFEM_SetDDWorkingDirectory "

Working directories for each subdomain.  

python usage : DDM_ExternalFEM_SetDDWorkingDirectory()  
";

%feature("docstring") DDM_ExternalFEM_ExSolver "

Solve fully the local contact problem in DDM.  

python usage : DDM_ExternalFEM_ExSolver(storage, checktype, tol, relax,
nb_iter_check, nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";

%feature("docstring") DDM_ExternalFEM_ExSolver_3D "

Solve fully the local contact problem.  

python usage : DDM_ExternalFEM_ExSolver_3D(storage, checktype, tol, relax,
nb_iter_check, nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";


// File: wrap__ASpxx_8h.xml

%feature("docstring") ASpxx_LoadTactors "

Load ASpxx from MAILx and Initialize existing_entities.  

python usage : ASpxx_LoadTactors()  
";

%feature("docstring") ASpxx_PushPreconNodes "

set ASpxx supporting nodes as precon  

python usage : ASpxx_PushPreconNodes()  
";

%feature("docstring") ASpxx_GetAllConnec "

return connectivity of all AS in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = ASxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of ASxxx elements  
";

%feature("docstring") ASpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all ASxxx  

python usage : idata, rdata = ASxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") ASpxx_CleanMemory "

Free all memory allocated within ASpxx module.  

python usage : ASpxx_CleanMemory()  
";

%feature("docstring") ASpxx_ExplodePatch "

Explode ASpxx patch in singleton.  

python usage : ASpxx_ExplodePatch()  
";


// File: wrap__surface__T3_8h.xml

%feature("docstring") surface_T3_compute_volume_inertia "

Computes the volume of an object described by a triangulated surface.  

**Warning**: 1) we assume size_coor is three times the number of nodes and
    size_connec is three times the number of elements python call: x_G, I,
    vol=surface_T3_compute_volume_inertia(coor, connec, 3, 9)  

Parameters
----------
coor_size(int): size of coor  
coor(double *): node coordinates  
connec_size(int): size of connec  
vol(double *): computed volume  
x_G(double *): mass center coordinates  
x_G_size(int): size of x_G  
I(double *): inertia matrix, stored a a vector  
I_size(int): size of I  
";

%feature("docstring") surface_T3_identify_entities "

Attributes an entity number to triangles, by computing connected components.  

**Warning**: 1) we assume size_connec is three times the number of elements and
    size_ele2entity is the number of elements python call:
    ele2entity=surface_T3_identify_entities(nbnode, max_adj_ele_2_node, connec,
    nbele)  

Parameters
----------
nbnode(int): the number of nodes  
max_adj_ele_2_node(int): the maximal number of adjacent elements per node  
connec(int *): connecivity of elements  
connec_size(int): size of connec  
ele2entity(double *): entity number for each element  
ele2entity_size(int): size of ele2entity  
";


// File: wrap__CSPRx_8h.xml

%feature("docstring") CSPRx_SelectProxTactors "

contact detection between CSxxx and PRxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : CSPRx_SelectProxTactors(int reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CSPRx_WriteLastVlocRloc "

write last local values of all CSPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : CSPRx_WriteLastVlocRloc()  

The values written are relative velocity, forces and local frame  
";

%feature("docstring") CSPRx_WriteOutVlocRloc "

write local values of all CSPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : CSPRx_WriteOutVlocRloc()  
";

%feature("docstring") CSPRx_DisplayOutVlocRloc "

display local values of all CSPRx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CSPRx_DisplayOutVlocRloc()  
";

%feature("docstring") CSPRx_DisplayProxTactors "

display contacts  

python usage : CSPRx_DisplayProxTactors()  
";

%feature("docstring") CSPRx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CSPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CSPRx_Trim "

trim contact (only node face contact)  

python usage : CSPRx_Trim()  
";

%feature("docstring") CSPRx_GetInfo "

return contact info for the icdan CSPRx contact  

python usage : a = CSPRx_GetInfo(icdan)  

Parameters
----------
icdan(integer) : contact identifiant  

Returns
-------
a (array integer) : info array  
";

%feature("docstring") CSPRx_Smoothing "

smooth contact reaction  

python usage : CSPRx_Smmothing()  
";

%feature("docstring") CSPRx_AddReac "

add contact force to body Reac  

python usage : CSPRx_AddReac()  
";

%feature("docstring") CSPRx_CleanMemory "

Free all memory allocated within CSPRx module.  

python usage : CSPRx_CleanMemory()  
";


// File: wrap__POLYR_8h.xml

%feature("docstring") POLYR_LoadTactors "

load POLYR from RBDY3 or MAILx and initialize existing_entites  

python usage : POLYR_LoadTactors()  
";

%feature("docstring") POLYR_GetContactorColor "

Get the color of a given POLYR.  

python usage : color = POLYR_GetContactorColor(itacty)  

Parameters
----------
itacty(integer) : rank of POLYR  

Returns
-------
color (string) : the color of the POLYR itact  
";

%feature("docstring") POLYR_SaveVertex "

write position of vertex in a file  

python usage : POLYR_SaveVertex()  
";

%feature("docstring") POLYR_ModifyRadius "

apply an amplification/reduction size factor  

python usage : POLYR_ModifyRadius(ratio)  

Parameters
----------
ratio(real) : ratio factor  
ratio(double) : ratio factor  
";

%feature("docstring") POLYR_SetThresholdBigPolyr "

define the threshold between a plain and a big polyr. big polyr are such that
radius > threshold*mean_radius. default threshold = 4. Must be defined before
the load of tactors.  

python usage : POLYR_SetThresholdBigPolyr(ratio)  

Parameters
----------
ratio(real) : ratio factor  
ratio(double) : ratio factor  
";

%feature("docstring") POLYR_SetBigPolyr "

impose explicitly that an object is big. Must be set after the load of tactors.  

python usage : POLYR_SetBigPolyr(itacty)  

Parameters
----------
itacty(integer) : rank of the polyr  
itacty(int) : rank of the polyr  
";

%feature("docstring") POLYR_SetNbBigPolyr "

impose explicitly the number of big POLYR. Must be set after the load of
tactors.  

python usage : POLYR_SetNbBigPolyr(nb)  

Parameters
----------
nb(integer) : number of polyr  
number(int) : number of polyr  
";

%feature("docstring") POLYR_SkipTopoBigPolyr "

skip the topological decomposition of a big POLYR. its surface is considered as
a soup of triangle. usefull with complicated surface using Cundall CP detection  

python usage : POLYR_SkipTopoBigPolyr()  
";

%feature("docstring") POLYR_SkipAutomaticReorientation "

disable automatic reorientation (which works only with convex POLYR).  

python usage : POLYR_SkipAutomaticReorientation()  

  
 Disable the automatic reorientation of normals performed by lmgc90.  
 This is necessary when using non-convex objects.  
";

%feature("docstring") POLYR_SkipHEBuild "

disable Half-Edge structure generation (HE is necessary for non convex contact
detection)  

python usage : POLYR_SkipHEBuild()  

  
 Disable the Half-Edge structure generation performed by lmgc90.  
 This is necessary when testing the import of strange object.  
";

%feature("docstring") POLYR_TopologyAngle "

set the maximum angle (between 0 and 180 degree) threshold between 2 elements to
declare them as belonging to the same topological face  

python usage : POLYR_TopologyAngle(angle)  
";

%feature("docstring") POLYR_FlatnessAngle "

set the maximum angle (between 0 and 180 degree) variation between elements of a
topological face to declare it as flat  

python usage : POLYR_FlatnessAngle(angle)  
";

%feature("docstring") POLYR_GetWireframe "

Get wireframe of a POLYR.  

python usage : coor,connectivity = POLYR_GetWireframe(itacty, angle)  

Parameters
----------
itacty(integer) : rank of the POLYR  
angle(double) : threshold angle to skip some nodes on boundary of faces of the
    POLYR  

Returns
-------
coor (double array) : reference on the coor vector seen as a numpy array of size
[nb_point,3] connectivity (integer array) : reference on the connectivity vector
seen as a numpy array  
";

%feature("docstring") POLYR_GetVertex "

Get the outline of a POLYR in almost current configuration.  

If the POLYR is a real POLYR the current position of the center of the POLYR is
used but the local frame for the orientation is the on in detection
configuration.  

If the POLYR is in fact a POLYD the current position of nodes of the mesh are
used.  

usage : vertex = POLYR_GetVertex(itacty)  

Parameters
----------
itacty(integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetPtrVertexTT "

Get a pointer on the outline of a POLYR in detection configuration.  

usage : vertex = POLYR_GetPtrVertexTT(itacty)  

Parameters
----------
itacty(integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetPtrNormalTT "

Get a pointer on the outline of a POLYR in detection configuration - be carefull
to move polyr.  

usage : normal = POLYR_GetPtrNormalTT(itacty)  

Parameters
----------
itacty(integer) : rank of considered POLYR  

Returns
-------
normal (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_MoveToConfigurationTT "

move the polyr in the configuration TT ; mandatory to get the wireframe in
deformed configuration  

python usage : POLYR_MoveToConfigurationTT()  
";

%feature("docstring") POLYR_GetPOLYR2BDYTY "

Get a copy of map POLYR2bdyty.  

usage : polyr2bdyty = POLYR_GetPOLYR2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") POLYR_GetPtrPOLYR2BDYTY "

Get a pointer on map POLYR2bdyty.  

usage : polyr2bdyty = POLYR_GetPtrPOLYR2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : a pointer in the polyr2bdyty map  
";

%feature("docstring") POLYR_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = POLYR_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") POLYR_GetNbPOLYR "

Get the number of POLYR.  

python usage : nb_POLYR = POLYR_GetNbPOLYR()  

Returns
-------
nb_POLYR (integer) : the number of POLYR  
";

%feature("docstring") POLYR_InitOutlines "

Get a reference on the outlines of all POLYR.  

usage : outlines = POLYR_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_POLYR  
";

%feature("docstring") POLYR_InitScalarFields "

Get a reference on the scalar fields of all POLYR.  

usage : scalarfields = POLYR_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_POLYR array  
";

%feature("docstring") POLYR_UpdatePostdata "

Update values of outlines_POLYR and scalarfields_POLYR pointers.  

usage : POLYR_UpdatePostdata()  
";

%feature("docstring") POLYR_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = POLYR_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
POLYR  
";

%feature("docstring") POLYR_GetNbScalarFields "

Get the number of scalar fields of a POLYR.  

python usage : nb_scalarfields = POLYR_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a POLYR  
";

%feature("docstring") POLYR_GetPtrAllConnectivities "

Get a reference on the connectivities of all POLYR.  

usage : connec = POLYR_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") POLYR_GetPtrConnectivity "

Get a reference on the connectivity of one POLYR.  

usage : connec = POLYR_GetPtrConnectivity(itacty)  

Parameters
----------
itacty(integer) : POLYR number  

Returns
-------
connec (integer 2D-array) : reference on connectivity  
";

%feature("docstring") POLYR_GetPtrVertexRef "

Get the position of the vertices of a POLYR in its inertia frame.  

usage : vertex = POLYR_GetPtrVertexRef(itacty)  

Parameters
----------
itacty(integer) : rank of considered POLYR  

Returns
-------
vertex (double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYR_GetTopoData "

Get for each face of all POLYR : contactor id, topo id, face id and face status.  

usage : topo_data = POLYR_GetTopoData()  

Returns
-------
topt_data (int 2D-array) : topology data of all faces of all POLYR  
";

%feature("docstring") POLYR_CleanMemory "

Free all memory allocated within POLYR module.  

python usage : POLYR_CleanMemory()  
";


// File: wrap__io__hdf5__hl_8h.xml

%feature("docstring") io_hdf5_initOutFile "

Init HDF5 file in which to write results.  

python usage : io_hdf5_initOutFile(filename)  

Parameters
----------
filename(string) : file in which to write  
";

%feature("docstring") io_hdf5_write "

write output data in HDF5 file (GPV, DOF and VlocRloc).  

python usage : io_hdf5_write()  
";

%feature("docstring") io_hdf5_write_last "

write output data in HDF5 file (GPV, DOF and VlocRloc) in file  

python usage : io_hdf5_write_last(filename)  

Parameters
----------
filename(string) : file in which to write  
";

%feature("docstring") io_hdf5_read "

read output data from HDF5 file (DOF and VlocRloc).  

python usage : io_hdf5_read(filename, step)  

Parameters
----------
filename(string) : file to read  
step(integer) : step number to read  
";

%feature("docstring") io_hdf5_cleanMemory "

cleanMemory of io_hdf5 module  

python usage : io_hdf5_cleanMemory()  
";

%feature("docstring") io_hdf5_fixVersion "

Will try to fix the file when reading it.  

Because the parameters changed within version 0, this flag is needed to fix the
file whend reading it.  

python usage : io_hdf5_fixVersion(version)  
";


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
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PRASp_CleanMemory "

Free all memory allocated within PRASp module.  

python usage : PRASp_CleanMemory()  
";


// File: wrap__cpg__3D_8h.xml

%feature("docstring") cpg_3D_ExIter "

Execute one CPG iteration over the contact loop.  

python usage cpg_3D_ExIter()  
";

%feature("docstring") cpg_3D_AfterIterCheck "

Control CPG convergence.  

python usage cpg_3D_AfterIterCheck()  
";

%feature("docstring") cpg_3D_ExPost "

Transfer local solution.  

python usage cpg_3D_ExPost()  
";

%feature("docstring") cpg_3D_ExPrep "

prepare the matrix and the RHS of the contact problem  

python usage cpg_3D_ExPrep()  
";

%feature("docstring") cpg_3D_ScaleRloc "

scale all local contact forces of a factor equal to 0.9 < f < 1.1  

python usage cpg_3D_ScaleRloc()  
";

%feature("docstring") cpg_3D_SetDiagonalPrecond "

active diagonal preconditioner  

python usage cpg_3D_SetDiagonalPrecond()  
";

%feature("docstring") cpg_3D_SetFrictionless "

active frictionless solver  

python usage cpg_3D_SetFrictionless()  
";

%feature("docstring") cpg_3D_BimodalContactOrder "

active bimodal list  

python usage : cpg_3D_BimodalContactOrder()  
";

%feature("docstring") cpg_3D_SetCheckType "

define numerical convergence of the NLGS algorithm  

python usage cpg_3D_SetCheckType(checktype, tol, idproj)  

Parameters
----------
chekctype(char[5]) : type of convergence check  
tol(double) : norm tolerance  
idproj(integer) :  
  
 convergence check keywords:  
 Quad : quadratic norm (faulty contacts are redeemed by accurate contacts;
laxist norm)  
 Maxm : maximum norm (faulty contacts must comply; severe norm)  
 where Quad,Maxm,QM/16 are keywords for the check test, and the following real
number is the tolerance value.  
 The identifiant projection parameter corrsponds to :  
 PYRAMIDAL APPROXIMATION (1)  
 Efficient but no more isotropic friction  
 NORMAL PROJECTION (2)  
 The basic projection but not really efficient  
 HYBRID CORRECTION (3)  
 Efficient for sphere but not really sense for other bodies.  
";

%feature("docstring") cpg_3D_NormCheck "

Active one step norm evolution.  

python usage : cpg_3D_norm_check()  
";

%feature("docstring") cpg_3D_ExSolver "

Solve fully the local contact problem.  

python usage : cpg_3D_ExSolver(checktype, tol, idpoj, nb_iter_check,
nb_block_iter)  

Parameters
----------
checktype(char[5]) : convergentce test keyword  
tol(double) : tolerance value  
idproj(integer) :  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";


// File: wrap__CDPLx_8h.xml

%feature("docstring") CDPLx_SelectProxTactors "

contact detection between CYLND and PLANx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list  

python usage : CDPLx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CDPLx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : CDPLx_SmoothForceComputation()  
";

%feature("docstring") CDPLx_WriteLastVlocRloc "

write last local values of all CDPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDPLx_WriteLastVlocRloc()  
";

%feature("docstring") CDPLx_WriteOutVlocRloc "

write local values of all CDPLx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDPLx_WriteOutVlocRloc()  
";

%feature("docstring") CDPLx_DisplayOutVlocRloc "

display local values of all CDPLx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CDPLx_DisplayOutVlocRloc()  
";

%feature("docstring") CDPLx_DisplayProxTactors "

display contacts  

python usage : CDPLx_DisplayProxTactors()  
";

%feature("docstring") CDPLx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : CDPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CDPLx_CleanMemory "

Free all memory allocated within CDPLx module.  

python usage : CDPLx_CleanMemory()  
";


// File: wrap__SPSPx_8h.xml

%feature("docstring") SPSPx_SelectProxTactors "

contact detection between SPxxx and SPxxx tactors  

python usage : SPSPx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPSPx_SmoothForceComputation "

recup values of local contact forces of the last time step  

python usage : SPSPx_SmoothForceComputation()  
";

%feature("docstring") SPSPx_WriteLastVlocRloc "

write last local values of all SPSPx contacts  

python usage : SPSPx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_WriteOutVlocRloc "

write local values of all SPSPx contacts  

python usage : SPSPx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_DisplayOutVlocRloc "

display local values of all SPSPx contacts  

python usage : SPSPx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPSPx_DisplayProxTactors "

display contacts  

python usage : SPSPx_DisplayProxTactors()  
";

%feature("docstring") SPSPx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPSPx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") SPSPx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPSPx_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(real) : period on x axis  
";

%feature("docstring") SPSPx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : SPSPx_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(real) : period on y axis  
";

%feature("docstring") SPSPx_SetNumberInterByContact "

define the number of interaction by contact (experimental)  

python usage : SPSPx_SetNumberInterByContact(nb_interactions)  

Parameters
----------
nb_interactions(integer) : number of interactions per contact  
";

%feature("docstring") SPSPx_SetContactRadius "

define the contact radius (experimental)  

python usage : SPSPx_SetContactRadius(radius)  

Parameters
----------
radius(real) : contact radius  
";

%feature("docstring") SPSPx_FdSelectProxTactors "

contact detection between SPHER and SPHER tactors  

python usage : SPSPx_FdSelectProxTactors()  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPSPx_CleanMemory "

Free all memory allocated within SPSPx module.  

python usage : SPSPx_CleanMemory()  
";


// File: wrap__meca__polygon_8h.xml

%feature("docstring") MecaPolyg_CentralKernel "

Compute the central kernel of an input surface (which may be composed of several
polygons)  

python usage : ck_pts = MecaPolyg_CentralKernel(points, sizes)  

Parameters
----------
points(double array) : coordinates of the points of the surface (2D)  
sizes(integer array) : number of vertices of each polygons of the surface  

Returns
-------
ck_pts (double array) : coordinates of the points of the central kernel  
";

%feature("docstring") MecaPolyg_StressField "
";


// File: wrap__CSxxx_8h.xml

%feature("docstring") CSxxx_LoadTactors "

Load CSxxx from MAILx and Initialize existing_entities.  

python usage : CSxxx_LoadTactors()  
";

%feature("docstring") CSxxx_PushPreconNodes "

set CSxxx supporting nodes as precon  

python usage : CSxxx_PushPreconNodes()  
";

%feature("docstring") CSxxx_FlipOrientation "

Flip normal of all CSxxx of a given MAILx body.  

python usage : CSxxx_FlipOrientation(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of desired body  
";

%feature("docstring") CSxxx_FlipOrientationOnePatch "

Flip normal of CSxxx belonging to given patch of a given MAILx body.  

python usage : CSxxx_FlipOrientationOnePatch(ibdyty,icspxx)  

Parameters
----------
ibdyty(integer) : rank of desired body  
icspxx(integer) : rank of desired patch  
";

%feature("docstring") CSxxx_SetShrink "

shrink position of nodes in CSxxx contactors  

python usage : CSxxx_SetShrink(shrink)  

Parameters
----------
shrink(real) : shrink value  
";

%feature("docstring") CSxxx_SetQuadrature "

Set the contact quadrature rule of a CSxxx face. OBSOLETE FUNCTION !!!! To
remove in the future.  

python usage : CSxxx_SetQuadrature(ivalue)  

Parameters
----------
ivalue(integer) : degree on CSxxx contactor  
";

%feature("docstring") CSxxx_AddReac "

Apply an external reaction on a CSxxx.  

python usage : CSxxx_AddReac(datatype, iCSxxx, reac)  

Parameters
----------
datatype(string of size 5) : the vector to set  
iCSxxx(integer) : id of the CSpxx  
reac(double array) : the value to add  
";

%feature("docstring") CSpxx_ApplySurfaceLoad "
";

%feature("docstring") CSpxx_ApplyPressure "

Apply an external pressure on a CSpxx.  

python usage : CSpxx_ApplyPressure(ivalue,rvalue)  

Parameters
----------
ivalue(integer) : id of the CSpxx  
rvalue(real) : pressure  
";

%feature("docstring") CSxxx_GetNbCSxxx "

Get the number of CSxxx.  

usage : nb_CSxxx = CSxxx_GetNbCSxxx()  

Parameters
----------
nb_CSxxx(integer) : number of CSxxx in container  
";

%feature("docstring") CSpxx_GetAllConnec "

return connectivity of all CS in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = CSxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of CSxxx elements  
";

%feature("docstring") CSpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all CSxxx  

python usage : idata, rdata = CSxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") CSxxx_CleanMemory "

Free all memory allocated within CSxxx module.  

python usage : CSxxx_CleanMemory()  
";


// File: wrap__SPDCx_8h.xml

%feature("docstring") SPDCx_SelectProxTactors "

contact detection between SPxxx and DCxxx tactors  

python usage : SPDCx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPDCx_SmoothForceComputation "

computes smooth forces (if any)  

python usage : SPDCx_SmoothForceComputation()  
";

%feature("docstring") SPDCx_WriteLastVlocRloc "

write last local values of all SPDCx contacts  

python usage : SPDCx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_WriteOutVlocRloc "

write local values of all SPDCx contacts  

python usage : SPDCx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_DisplayOutVlocRloc "

display local values of all SPDCx contacts  

python usage : SPDCx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPDCx_DisplayProxTactors "

display contacts  

python usage : SPDCx_DisplayProxTactors()  
";

%feature("docstring") SPDCx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPDCx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") SPDCx_CleanMemory "

Free all memory allocated within SPDCx module.  

python usage : SPDCx_CleanMemory()  
";


// File: wrap__DNLYC_8h.xml

%feature("docstring") DNLYC_LoadTactors "

load DNLYC from RBDY3 and initialize existing_entites  

python usage : DNLYC_LoadTactors()  
";

%feature("docstring") DNLYC_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = DNLYC_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") DNLYC_GetNbDNLYC "

Get the number of DNLYC.  

python usage : nb_DNLYC = DNLYC_GetNbDNLYC()  

Returns
-------
nb_DNLYC (integer) : the number of DNLYC  
";

%feature("docstring") DNLYC_GetPtrDNLYC2BDYTY "

return a pointer onto the map dnlyc2bdyty  

python usage : dnlyc2bdyty = DNLYC_GetPtrDNLYC2BDYTY()  

Returns
-------
dnlyc2bdyty (integer array) : reference on map between dnlyc rank and body rank  
";

%feature("docstring") DNLYC_InitOutlines "

Get a reference on the outlines of all DNLYC.  

usage : outlines = DNLYC_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_DNLYC  
";

%feature("docstring") DNLYC_InitScalarFields "

Get a reference on the scalar fields of all DNLYC.  

usage : scalarfields = DNLYC_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_DNLYC array  
";

%feature("docstring") DNLYC_UpdatePostdata "

Update values of outlines_DNLYC and scalarfields_DNLYC pointers.  

usage : DNLYC_UpdatePostdata  
";

%feature("docstring") DNLYC_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = DNLYC_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
DNLYC  
";

%feature("docstring") DNLYC_GetNbScalarFields "

Get the number of scalar fields of a DNLYC.  

python usage : nb_scalarfields = DNLYC_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a DNLYC  
";

%feature("docstring") DNLYC_GetPtrAllConnectivities "

Get a reference on the connectivities of all DNLYC.  

usage : connec = DNLYC_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") DNLYC_CleanMemory "

Free all memory allocated within DNLYC module.  

python usage : DNLYC_CleanMemory()  
";


// File: wrap__DKKDx_8h.xml

%feature("docstring") DKKDx_SelectProxTactors "

contact detection between DISKx and xKSID tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKKDx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKKDx_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKKDx_SmoothForceComputation()  
";

%feature("docstring") DKKDx_WriteLastVlocRloc "

write last local values of all DKKDx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKKDx_WriteLastVlocRloc()  
";

%feature("docstring") DKKDx_WriteOutVlocRloc "

write local values of all DKKDx contacts  

The values written are relative velocity, forces and local frame  

python usage : DKKDx_WriteOutVlocRloc()  
";

%feature("docstring") DKKDx_DisplayOutVlocRloc "

display local values of all DKKDx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKKDx_DisplayOutVlocRloc()  
";

%feature("docstring") DKKDx_DisplayProxTactors "

display contacts  

python usage : DKKDx_DisplayProxTactors()  
";

%feature("docstring") DKKDx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKKDx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKKDx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKKDx_SetSurfaceSectors(nbsect)  

Parameters
----------
nbsect(integer) : number of sectors  
";

%feature("docstring") DKKDx_CleanMemory "

Free all memory allocated within DKKDx module.  

python usage : DKKDx_CleanMemory()  
";


// File: wrap__ExternalModels_8h.xml

%feature("docstring") ExternalModels_InitModels "

Initialize the external models (if any)  

python usage : ExternalModels_InitModels()  
";

%feature("docstring") ExternalModels_StoreProperties "

Store external models (if any)  

python usage : ExternalModels_StoreProperties()  
";

%feature("docstring") ExternalModels_CleanMemory "

Free all memory allocated within ExternalModels module.  

python usage : ExternalModels_CleanMemory()  
";


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
ivalue(integer) : id of the multiMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a multiMAILx  
";

%feature("docstring") multiMAILx_GetNbElements "

Get the number of elements of a multiMAILx.  

python usage : nb_elements = multiMAILx_GetNbElements(ibdyty)  

Parameters
----------
ivalue(integer) : id of the multiMAILx  

Returns
-------
nb_nodes (integer) : number of elements of a multiMAILx  
";

%feature("docstring") multiMAILx_IsVisible "

return if a given body visible  

python usage : visible = multiMAILx_IsVisible(ibdyty)  

Parameters
----------
idbdy(integer): id of the body we want visibility  

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
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  

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
datatype(string of size 5) : the vector to set  
ibdyty(integer) : rank of body  
matrix(double array) : the new values  
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
num(integer) : which GPV file to read  
";

%feature("docstring") multiMAILx_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniDof last call  

python usage : multiMAILx_ReadIniDof(num=0)  

Parameters
----------
num(integer) : which DOF file to read  
";

%feature("docstring") multiMAILx_WriteLastDof "

Write DOF.LAST file.  

python usage : multiMAILx_WriteLastDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write dof if omitted works on all
    objects  
";

%feature("docstring") multiMAILx_WriteOutDof "

Write DOF.OUT file.  

python usage : multiMAILx_WriteOutDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to write dof if omitted works on all
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
i_list(list of integer) : list of bodies to compute mass and inertia if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_ComputeBulk "

computes elementary stiffness and viscosity matrices of a list of bodies  

python usage : multiMAILx_ComputeBulk(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute stiffness and viscosity
    matrices if omitted works on all objects  
";

%feature("docstring") multiMAILx_ComputeFext "

compute elementary external forces of a list of bodies  

python usage : multiMAILx_ComputeFext(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute external forces if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_AssembKT "

assemble pseudo mass matrix and apply drvdof of a list of bodies  

python usage : multiMAILx_AssembKT(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to assemble pseudo mass matrix and apply
    drvdof if omitted works on all objects  
";

%feature("docstring") multiMAILx_AssembRHS "

assembles right hand side of a list of bodies  

python usage : multiMAILx_AssembRHS(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to assemble right hand side if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_ComputeResidueNorm "

computes the norm of the residue of a list of bodies  

python usage : norm = multiMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute the norm of the residue if
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
i_list(list of integer) : list of bodies to compute free state if omitted works on
    all objects  
";

%feature("docstring") multiMAILx_ComputeDof "

computes the current d.o.f knowing all the forces/fluxses (free + contact) of a
list of bodies  

python usage : multiMAILx_ComputeDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute current d.o.f if omitted works
    on all objects  
";

%feature("docstring") multiMAILx_ComputeField "

computes elementary fields of a list of bodies  

python usage : multiMAILx_ComputeField(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_UpdateBulk "

update begin elementary fields with current elementary fields of a list of
bodies  

python usage : multiMAILx_UpdateBulk(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute elementary fields if omitted
    works on all objects  
";

%feature("docstring") multiMAILx_UpdateDof "

update begin d.o.f. with current d.o.f. of a list of bodies  

python usage : multiMAILx_UpdateDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to update current d.o.f if omitted works
    on all objects  
";

%feature("docstring") multiMAILx_GetScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = multiMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") multiMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") multiMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = multiMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired vector field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") multiMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : multiMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") multiMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = multiMAILx_GetConnectivity(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") multiMAILx_GetCoor "

return node coordinates of idBody  

python usage : array = multiMAILx_GetCoor(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") multiMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = multiMAILx_GetAll(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") multiMAILx_GetElementsVolume "

return volume of elements  

python usage : volumes = multiMAILx_GetElementsVolume(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

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
IdBody(integer) : id of the concerned body  
tol(double) : tolerance  

Returns
-------
array (double 2D-array) : neighbor[nb_ele,max_neighbors]  
";

%feature("docstring") multiMAILx_GetPtrElementsEnergy "

return pointer on energy of elements  

python usage : energies = multiMAILx_GetPtrElementsEnergy(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

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
IdBody(integer) : id of the concerned body  

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
ibdyty(integer) : rank of the multiMAILx  

Returns
-------
eviz (int array) : reference on the desired vector seen as a numpy array  
";

%feature("docstring") multiMAILx_GetDeformationEnergy "

Get the deformation energy of a given displacement field.  

python usage : energy = multiMAILx_GetDeformationEnergy(id,displacement)  

Parameters
----------
ibdyty(integer) : rank of considered body  
displacement(double matrix) : displacement field  

Returns
-------
energy (double) : deformation energy  
";

%feature("docstring") multiMAILx_GetPtrBoundaryElements "

return boundary elements  

python usage : vector = multiMAILx_GetPtrBoundaryElements(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
vector (integer) : for each element =0 no boundary, otherwise gives the number
of free edge/face  
";

%feature("docstring") multiMAILx_CleanMemory "

Free all memory allocated within multiMAILx module.  

python usage : multiMAILx_CleanMemory()  
";


// File: wrap__JONCx_8h.xml

%feature("docstring") JONCx_LoadTactors "

load JONCx from RBDY2 and initialize existing_entites  

python usage : JONCx_LoadTactors()  
";

%feature("docstring") JONCx_GetNbJONCx "

Get the number of JONCx in container.  

python usage : nb_joncx = JONCx_GetNbJONCx()  

Returns
-------
nb_joncx (integer) : the number of JONCx in container  
";

%feature("docstring") JONCx_GetBodyId "

Get the body rank of a given JONCx.  

python usage : ibdyty = JONCx_GetBodyId(itacty)  

Parameters
----------
itacty(integer) : JONCx rank  

Returns
-------
ibdyty (integer) : body rank  
";

%feature("docstring") JONCx_GetShape "

Get the shape of a JONCx.  

usage : shape = JONCx_GetShape(itacty)  

Parameters
----------
itacty(integer) : rank of JONCx  

Returns
-------
shape (double array) : axis length of the JONCx  
";

%feature("docstring") JONCx_GetCoor "

Get the coor of a JONCx.  

usage : coor = JONCx_GetCoor(itacty)  

Parameters
----------
itacty(integer) : rank of JONCx  

Returns
-------
coor (double array) : coordinates of the JONCx  
";

%feature("docstring") JONCx_GetPtrJONCx2BDYTY "

return a pointer onto the map joncx2rbdy2  

python usage : joncx2rbdy2 = JONCx_GetPtrJONCx2BDYTY()  

Returns
-------
joncx2rbdy2 (integer array) : reference on map between joncx rank and body/tact
rank  
";

%feature("docstring") JONCx_IsVisible "

return if a body visible  

usage : visible = JONCx_IsVisible(itact)  

Parameters
----------
itact(integer) : rank of JONCx  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") JONCx_InitOutlines "

Get a reference on the outlines of all JONCx.  

usage : outlines = JONCx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_JONCx  
";

%feature("docstring") JONCx_InitScalarFields "

Get a reference on the scalar fields of all JONCx.  

usage : scalarfields = JONCx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_JONCx array  
";

%feature("docstring") JONCx_UpdatePostdata "

Update values of outlines_JONCx and scalarfields_JONCx pointers.  

usage : JONCx_UpdatePostdata()  
";

%feature("docstring") JONCx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = JONCx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
JONCx  
";

%feature("docstring") JONCx_GetNbScalarFields "

Get the number of scalar fields of a JONCx.  

python usage : nb_scalarfields = JONCx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a JONCx  
";

%feature("docstring") JONCx_CleanMemory "

Free all memory allocated within JONCx module.  

python usage : JONCx_CleanMemory()  
";


// File: wrap__PTPT2_8h.xml

%feature("docstring") PTPT2_SelectProxTactors "

contact detection between PT2Dx tactors  

python usage : PTPT2_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") PTPT2_WriteLastVlocRloc "

write last local values of all PTPT2 contacts  

python usage : PTPT2_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_WriteOutVlocRloc "

write local values of all PTPT2 contacts  

python usage : PTPT2_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_DisplayOutVlocRloc "

display local values of all PTPT2 contacts  

python usage : PTPT2_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") PTPT2_DisplayProxTactors "

display contacts  

python usage : PTPT2_DisplayProxTactors()  
";

%feature("docstring") PTPT2_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : PTPT2_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PTPT2_LoadNetwork "

read a PTPT2 network from a file  

python usage : PTPT2_LoadNetwork()  
";

%feature("docstring") PTPT2_SetTolerance "

set the maximum violation for a point to point link  

python usage : PTPT2_SetTolerance(tol)  
";

%feature("docstring") PTPT2_SetExplicitLocalFrame "

local frame is computed only once at the first step  

python usage : PTPT2_SetExplicitLocalFrame()  
";

%feature("docstring") PTPT2_LoadParams "

read a PTPT2 surface and l0 from a file  

python usage : PTPT2_LoadParams()  
";

%feature("docstring") PTPT2_UseCurrentNonuc0 "

Use GetCoor or value given from file insted of computing nonuc0 from reference
coordinates.  

python usage : PTPT2_UseCurrentNonuc0(to_use) param[in] to_use (integer) : 1 to
activate, 0 to deactivate feature  
";

%feature("docstring") PTPT2_CleanMemory "

Free all memory allocated within PTPT2 module.  

python usage : PTPT2_CleanMemory()  
";


// File: wrap__overall_8h.xml

%feature("docstring") overall_Initialize "

Initialize LMGC90.  

python usage : overall_Initialize()  
";

%feature("docstring") overall_Finalize "

Finalize LMGC90.  

python usage : overall_Finalize()  
";

%feature("docstring") overall_InitEntityList "

Initialize entity list : must be done after LoadTactors.  

python usage : overall_InitEntityList()  
";

%feature("docstring") TimeEvolution_SetTimeStep "

Set value of the time step.  

python usage : TimeEvolution_SetTimeStep(dt)  

Parameters
----------
dt(double) : value of time step  
";

%feature("docstring") TimeEvolution_IncrementStep "

Increment curent time, time step and eventually initialize NR loop counter.  

python usage : TimeEvolution_IncrementStep()  
";

%feature("docstring") TimeEvolution_UpdateStep "

update the initial time to the current time  

python usage : TimeEvolution_UpdateStep()  
";

%feature("docstring") TimeEvolution_DisplayStep "

Display time evolution step informations.  

python usage : TimeEvolution_DisplayStep()  
";

%feature("docstring") TimeEvolution_SetInitialStep "

Set the rank of the first time step.  

python usage : TimeEvolution_SetInitialStep(first_step)  

Parameters
----------
first_step(integer) : rank of the first time step  
";

%feature("docstring") TimeEvolution_SetInitialTime "

Set initial time.  

python usage : TimeEvolution_SetInitialTime(t_init)  

Parameters
----------
t_init(double) : initial time  
";

%feature("docstring") TimeEvolution_GetTime "

get current time  

python usage : time = TimeEvolution_GetTime()  

Returns
-------
time (double) : current time  
";

%feature("docstring") TimeEvolution_GetTimeStep "

get current time step  

python usage : dt = TimeEvolution_GetTimeStep()  

Returns
-------
dt (double) : time step  
";

%feature("docstring") TimeEvolution_GetStep "

get current step number  

python usage : it = TimeEvolution_GetStep()  

Returns
-------
it (int) : current step number  
";

%feature("docstring") TimeEvolution_WriteLastDof "

python usage : TimeEvolution_WriteLastDof()  
";

%feature("docstring") TimeEvolution_WriteOutDof "

python usage : TimeEvolution_WriteOutDof(Nstep_writeDof)  

Parameters
----------
Nstep_writeDof(integer) : periodicity of DOF write  
";

%feature("docstring") TimeEvolution_DisplayOutDof "

python usage : TimeEvolution_DisplayOutDof()  
";

%feature("docstring") TimeEvolution_WriteLastRnod "

python usage : TimeEvolution_WriteLastRnod()  
";

%feature("docstring") TimeEvolution_WriteOutRnod "

python usage : TimeEvolution_WriteOutRnod(nstep)  

Parameters
----------
nstep(integer) : a freq of writing  
";

%feature("docstring") TimeEvolution_DisplayOutRnod "

python usage : TimeEvolution_DisplayOutRnod()  
";

%feature("docstring") TimeEvolution_WriteLastVlocRloc "

python usage : TimeEvolution_WriteLastVlocRloc()  
";

%feature("docstring") TimeEvolution_WriteOutVlocRloc "

python usage : TimeEvolution_WriteOutVlocRloc(nstep)  

Parameters
----------
nstep(integer) : a freq of writing  
";

%feature("docstring") TimeEvolution_DisplayOutVlocRloc "

python usage : TimeEvolution_DisplayOutVlocRloc()  
";

%feature("docstring") TimeEvolution_WriteLastGPV "

python usage : TimeEvolution_WriteLastGPV()  
";

%feature("docstring") TimeEvolution_WriteOutGPV "

python usage : TimeEvolution_WriteOutGPV(nstep)  

Parameters
----------
nstep(integer) : a freq of writing  
";

%feature("docstring") TimeEvolution_ReadIniDof "

Read header of a DOF file.  

python usage : TimeEvolution_ReadIniDof(num=0)  

Parameters
----------
num(integer) : num of file to read  
";

%feature("docstring") TimeEvolution_ReadIniVlocRloc "

Read header of a VlocRloc file.  

python usage : TimeEvolution_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : num of file to read  
";

%feature("docstring") TimeEvolution_ReadIniGPV "

Read header of a GPV file.  

python usage : TimeEvolution_ReadIniGPV(num=0)  

Parameters
----------
num(integer) : num of file to read  
";

%feature("docstring") NewtonRaphson_Initialize "

initialize Newton Raphson Loop  

python usage : NewtonRaphson_Initialize(tol)  

Parameters
----------
tol(double) : tolerance  
";

%feature("docstring") NewtonRaphson_CheckConvergence "

check if Newton Raphson loop converges  

python usage : iconv = NewtonRaphson_CheckConvergence(norm)  

Parameters
----------
norm(double) : value to check  

Returns
-------
iconv (integer) : convergence status  

*   iconv = 0 : converges  
*   iconv = 1 : unknown  
*   iconv = 2 : diverges  
";

%feature("docstring") NewtonRaphson_ComputeTimeStep "

manages time step evolution depending on newton raphson convergence  

python usage : itodo = NewtonRaphson_ComputeTimeStep()  

Returns
-------
itodo (integer) : what to do now  

*   itodo = 0 : just keep going (time step may have been modified)  
*   itodo = 1 : redo time step (time step has been decreased)  
*   itodo = 2 : it's hopeless just stop where you are  
";

%feature("docstring") NewtonRaphson_SetMinTimeStep "

Set value of the mininum possible time step.  

python usage : NewtonRaphson_SetMinTimeStep(dt)  

Parameters
----------
dt(double) : minimum value of time step  
  
 Needed only if adaptive time step feature is used  
";

%feature("docstring") NewtonRaphson_SetMaxTimeStep "

Set value of the maximum possible time step.  

python usage : NewtonRaphson_SetMaxTimeStep(dt)  

Parameters
----------
dt(double) : maximum value of time step  
  
 Needed only if adaptive time step feature is used  
";

%feature("docstring") NewtonRaphson_SetFinalTime "

Set final time.  

python usage : NewtonRaphson_SetFinalTime(t_final)  

Parameters
----------
t_final(double) : final time  
";

%feature("docstring") NewtonRaphson_SetMaxIter "

Max number of iterations - default is 50.  

python usage : NewtonRaphson_SetMaxIter(max_iter)  

Parameters
----------
max_iter(integer) :  
";

%feature("docstring") NewtonRaphson_SetGoodIter "

Set the max number of iterations for good convergence - default is 10.  

python usage : NewtonRaphson_SetGoodIter(good_iter)  

Parameters
----------
good_iter(integer) :  
";

%feature("docstring") NewtonRaphson_SetBadIter "

Set the max number of iterations for bad convergence - default is 30.  

python usage : NewtonRaphson_SetBadIter(bad_iter)  

Parameters
----------
bad_iter(integer) :  
";

%feature("docstring") NewtonRaphson_SetIncPatience "

Set the number of increments to adapt the time step when successive good
convergence (increase time step) or bad convergence (decrease time step) -
default is 3.  

python usage : NewtonRaphson_SetIncPatience(patience)  

Parameters
----------
patience(integer) :  
";

%feature("docstring") overall_SelectProxTactors "

Prepare contact detection.  

python usage : overall_SelectProxTactors(Nstep_rough_seek)  

Parameters
----------
Nstep_rough_seek(integer) : periodicity of rough detection  
";

%feature("docstring") overall_DisplayProxTactors "

python usage : overall_DisplayProxTactors()  
";

%feature("docstring") overall_DIME "

set space dimension and in 2D the modelling assumption  

python usage : overall_DIME(idim, imod)  

Parameters
----------
idim(integer) : dimension (2 or 3)  
imod(integer) : kind of model (2D only)  

*   imod = 1 => plane strain  
*   imod = 2 => plane stress  
*   imod = 3 => axisymmetric  
";

%feature("docstring") Integrator_InitTheta "

python usage : Integrator_InitTheta(theta)  

Parameters
----------
theta(double) : value of theta in integrator  
";

%feature("docstring") Integrator_InitQS "

python usage : Integrator_InitQS()  
";

%feature("docstring") Integrator_InitCrankNickolson "

python usage : Integrator_InitCrankNickolson(theta)  

Parameters
----------
theta(double) : value of theta in integrator  
";

%feature("docstring") Integrator_InitGear "

python usage : Integrator_InitGear()  
";

%feature("docstring") Integrator_InitVerlet "

python usage : Integrator_InitVerlet()  
";

%feature("docstring") Integrator_InitBeta2 "

python usage : Integrator_InitBeta2(value)  

Parameters
----------
value(double) : numeric diffusion ([0.5,1] and 0.5 is conservative)  
";

%feature("docstring") Integrator_SetContactDetectionConfiguration "

set the parameters necessary to define the contact detection configuration
(default: 1-theta, 0.)  

python usage : Integrator_SetContactDetectionConfiguration(alpha_b,alpha_e)  

Parameters
----------
alpha_b(double) : value of the V_begin weight  
alpha_e(double) : value of the V weight  
";

%feature("docstring") overall_RequireXxlComputation "

python usage : overall_RequireXxlComputation()  
";

%feature("docstring") overall_UpdatePostData "

python usage : overall_UpdatePostData()  
";

%feature("docstring") overall_InitPostData "

python usage : overall_InitPostData(ifirst, ilast)  

Parameters
----------
ifirst(integer) :  
ilast(integer) :  
";

%feature("docstring") overall_SetWorkingDirectory "

python usage : overall_SetWorkingDirectory(path)  

Parameters
----------
path(string) : set path to DATBOX directory  
";

%feature("docstring") overall_GetWorkingDirectory "

python usage : path = overall_GetWorkingDirectory()  

Returns
-------
path (string) : working directory  
";

%feature("docstring") overall_WriteDrivenDof "

python usage : overall_WriteDrivenDof()  
";

%feature("docstring") overall_WriteOutDisplayFile "

python usage : overall_WriteOutDisplayFile(freq_display)  

Parameters
----------
freq_display(integer) : periodicity of display write  
";

%feature("docstring") TimeEvolution_ReadIniMpValues "

Read header of a MP_VALUES file.  

python usage : TimeEvolution_ReadIniMpValues(num=0)  

Parameters
----------
num(integer) : num of file to read  
";

%feature("docstring") TimeEvolution_WriteOutMpValues "

python usage : TimeEvolution_WriteOutMpValues(nstep)  

Parameters
----------
nstep(integer) : a freq of writing  
";

%feature("docstring") TimeEvolution_WriteLastMpValues "

python usage : TimeEvolution_WriteLastMpValues()  
";

%feature("docstring") overall_WriteBodies "

python usage : overall_WriteBodies()  
";

%feature("docstring") overall_CleanOutBodies "

python usage : overall_CleanOutBodies()  
";

%feature("docstring") overall_RebuildInBodies "

python usage : overall_RebuildInBodies()  
";

%feature("docstring") overall_CleanWriteOutFlags "

python usage : overall_CleanWriteOutFlags()  
";

%feature("docstring") overall_UseExperimentalDev "

Activate some unstable devs.  

python usage : overall_UseExperimentalDev()  
";

%feature("docstring") overall_UseExternalFem "

Allow to use the externalFem library instead of lmgc90 Fem lib.  

python usage : overall_UseExternalFem()  
";

%feature("docstring") overall_GetMaxInternalTact "

get max of internal for tact  

python usage : nb = overall_GetMaxInternalTact()  

Returns
-------
nb (integer) : maximum number of internal for interactions  
";


// File: wrap__mesh2D_8h.xml

%feature("docstring") mesh2D_GetIndicesMeshQ4 "

this function gives the couple (i, j) of indices coresponding to a given node n  

**Warning**: python call: [i, j]=mesh2D_GetIndicesMeshQ4(n)  

Parameters
----------
n(int): the given node  
i(int *): index in the u direction  
j(int *): index in the v direction  
";

%feature("docstring") mesh2D_SizeMeshQ4 "

this function computes the sizes of vectors used to store a mesh made of Q4 in
the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4 6  *---*---*  
     | 1 | 2 |  *---*---*  
     1 3 5  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]  
*   nb_node_per_ele = [4, 4]  
*   conn = [1, 3, 4, 2, 3, 5, 6, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMeshQ4(nb_elem_x, nb_elem_y)  

    Parameters:  
    nb_elem_x(int): number of elements in the horizontal direction  
    nb_elem_y(int): number of elements in the vertical direction  
    size_coor(int *): size of coor  
    size_nb_node_per_ele(int *): size of nb_node_per_ele  
    size_conn(int *): size of conn  
";

%feature("docstring") mesh2D_SizeMesh2T3 "

this function computes the sizes of vectors used to store a mesh made of T3 ---
obtained by splitting a Q4 in two T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *----*  
     | 1 /|  
     | / |  
     | / |  
     |/ 2 |  *----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4]  
*   nb_node_per_ele = [3, 3]  
*   conn = [1, 3, 4, 2, 1, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh2T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int *): size of coor  
    size_nb_node_per_ele(int *): size of nb_node_per_ele  
    size_conn(int *): size of conn  
";

%feature("docstring") mesh2D_SizeMesh4T3 "

this function computes the sizes of vectors used to store a mesh made of T3 ---
obtained by splitting a Q4 in four T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *-----*  
     |\\ 4 /|  
     | \\ / |  
     |1 5 3|  
     | / \\ |  
     |/ 2 |  *-----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5]  
*   nb_node_per_ele = [3, 3, 3, 3]  
*   conn = [1, 5, 2, 1, 3, 5, 3, 4, 5, 2, 5, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh4T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int *): size of coor  
    size_nb_node_per_ele(int *): size of nb_node_per_ele  
    size_conn(int *): size of conn  
";

%feature("docstring") mesh2D_SizeMeshQ8 "

this function computes the sizes of vectors used to store a mesh made of Q8
following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     3 7 3  *---*---*  
     | |  
     | |  
     8 * 1 * 6  
     | |  
     | |  *---*---*  
     1 5 2  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8]  
*   nb_node_per_ele = [8]  
*   conn = [1, 2, 3, 4, 5, 6, 7, 8]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh4T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int *): size of coor  
    size_nb_node_per_ele(int *): size of nb_node_per_ele  
    size_conn(int *): size of conn  
";

%feature("docstring") mesh2D_MeshQ4 "

this function computes and returns a mesh made of Q4 in the following generic
format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4 6  *---*---*  
     | 1 | 2 |  *---*---*  
     1 3 5  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]  
*   nb_node_per_ele = [4, 4]  
*   conn = [1, 3, 4, 2, 3, 5, 6, 4]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_MeshQ4(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    x0(double): abscissa of the lower left corner of the rectangle  
    y0(double): ordinate of the lower left corner of the rectangle  
    lx(double): length of the mesh, following the axis (Ox)  
    ly(double): length of the mesh, following the axis (Oy)  
    nb_elem_x(int): number of elements in the horizontal direction  
    nb_elem_y(int): number of elements in the vertical direction  
    size_coor(int): size of coor  
    size_nb_node_per_ele(int): size of nb_node_per_ele  
    size_conn(int): size of conn  
    coor(double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    nb_node_per_ele(int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    conn(int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_Mesh2T3 "

this function computes an returns a mesh made of T3 --- obtained by splitting a
Q4 in two T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *----*  
     | 1 /|  
     | / |  
     | / |  
     |/ 2 |  *----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4]  
*   nb_node_per_ele = [3, 3]  
*   conn = [1, 3, 4, 2, 1, 4]  
*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_Mesh2T3(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    x0(double): abscissa of the lower left corner of the rectangle  
    y0(double): ordinate of the lower left corner of the rectangle  
    lx(double): length of the mesh, following the axis (Ox)  
    ly(double): length of the mesh, following the axis (Oy)  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int): size of coor  
    size_nb_node_per_ele(int): size of nb_node_per_ele  
    size_conn(int): size of conn  
    coor(double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    nb_node_per_ele(int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    conn(int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_Mesh4T3 "

this function computes and return a mesh made of T3 --- obtained by splitting a
Q4 in four T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *-----*  
     |\\ 4 /|  
     | \\ / |  
     |1 5 3|  
     | / \\ |  
     |/ 2 |  *-----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5]  
*   nb_node_per_ele = [3, 3, 3, 3]  
*   conn = [1, 5, 2, 1, 3, 5, 3, 4, 5, 2, 5, 4]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_Mesh4T3(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    x0(double): abscissa of the lower left corner of the rectangle  
    y0(double): ordinate of the lower left corner of the rectangle  
    lx(double): length of the mesh, following the axis (Ox)  
    ly(double): length of the mesh, following the axis (Oy)  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int): size of coor  
    size_nb_node_per_ele(int): size of nb_node_per_ele  
    size_conn(int): size of conn  
    coor(double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    nb_node_per_ele(int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    conn(int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_MeshQ8 "

this function computes and returns a mesh made of Q8 in the following generic
format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     3 7 3  *---*---*  
     | |  
     | |  
     8 * 1 * 6  
     | |  
     | |  *---*---*  
     1 5 2  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8]  
*   nb_node_per_ele = [8]  
*   conn = [1, 2, 3, 4, 5, 6, 7, 8]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_MeshQ8(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    x0(double): abscissa of the lower left corner of the rectangle  
    y0(double): ordinate of the lower left corner of the rectangle  
    lx(double): length of the mesh, following the axis (Ox)  
    ly(double): length of the mesh, following the axis (Oy)  
    nb_elem_x(int): number of elements Q4 in the horizontal direction  
    nb_elem_y(int): number of elements Q4 in the vertical direction  
    size_coor(int): size of coor  
    size_nb_node_per_ele(int): size of nb_node_per_ele  
    size_conn(int): size of conn  
    coor(double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    nb_node_per_ele(int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    conn(int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";


// File: wrap__PT2DL_8h.xml

%feature("docstring") PT2DL_LoadTactors "

Initialize existing_entities variable for PT2DL contactors.  

python usage : PT2DL_LoadTactors()  
";

%feature("docstring") PT2DL_PushPreconNodes "

python usage : PT2DL_PushPreconNodes()  
";

%feature("docstring") PT2DL_GetNbPT2DL "

Get the number of PT2DL.  

usage : nb_PT2DL = PT2DL_GetNbPT2DL()  

Parameters
----------
nb_PT2DL(integer) : number of PT2DL in container  
";

%feature("docstring") PT2DL_GetNbPT2TL "

Get the number of PT2TL of a body.  

usage : nb_PT2DL = PT2DL_GetNbPT2TL(ibdyty)  

Parameters
----------
nb_PT2TL(integer) : number of PT2TL in container  
";

%feature("docstring") PT2DL_ComputeConvectiveFlux "

python usage : PT2DL_ComputeConvectiveFlux()  
";

%feature("docstring") PT2DL_AssembThermKT "

python usage : PT2DL_AssembThermKT()  
";

%feature("docstring") PT2DL_AssembThermRHS "

python usage : PT2DL_AssembThermRHS()  
";

%feature("docstring") PT2DL_GetBody "

return corresponding body  

python usage : ibdy = PT2DL_GetBody(itacty)  
";

%feature("docstring") PT2DL_CleanMemory "

Free all memory allocated within PT2DL module.  

python usage : PT2DL_CleanMemory()  
";


// File: wrap__user_8h.xml

%feature("docstring") user_getWoodFrame "
";


// File: wrap__DISKL_8h.xml

%feature("docstring") DISKL_LoadTactors "

Load DISKL from MAILx and Initialize existing_entities.  

python usage : DISKL_LoadTactors()  
";

%feature("docstring") DISKL_PushPreconNodes "

declare the DISKL supporting nodes as precon  

python usage : DISKL_PushPreconNodes()  
";

%feature("docstring") DISKL_CleanMemory "

Free all memory allocated within DISKL module.  

python usage : DISKL_CleanMemory()  
";


// File: wrap__postpro_8h.xml

%feature("docstring") postpro_PostproDuringComputation "

Scan postprocessing function which should be call during the computation
process.  

python usage : postpro_PostproDuringComputation()  
";

%feature("docstring") postpro_ReadCommands "

Scan postprocessing function which should be call during the computation
process.  

python usage : postpro_ReadCommands()  
";

%feature("docstring") postpro_PostproBeforeComputation "

Data initialization and scan postprocessing function which should be called
before the computation process.  

python usage : postpro_PostproBeforeComputation(restart=0) param[in] restart
(integer) : if the Postpro file must append to existing ones and starting index
of CONTACT_FORCE_DISTRIBUTION files  
";

%feature("docstring") postpro_FlushDuringComputation "

Flush all postpro files.  

python usage : postpro_FlushDuringComputation()  
";

%feature("docstring") postpro_ClosePostproFiles "

Close all postpro files.  

python usage : postpro_ClosePostproFiles()  
";

%feature("docstring") postpro_SetCircularSelectionZone "

Initialize data for postreatment using a circular selection.  

python usage : postpro_SetCircularSelectionZone(rvalue1, rvalu2, rvalue3)  

Parameters
----------
rvalue1(double) : X coordinate  
rvalue2(double) : Y coordinate  
rvalue3(double) : radius selection  
";

%feature("docstring") postpro_MoveCircularSelectionZone "

Increment the position of the circular selection defined with
CIRCULAR_SELECTION.  

python usage : postpro_MoveCircularSelectionZone(rvalue1, rvalu2)  

Parameters
----------
rvalue1(double) : X translational velocity  
rvalue2(double) : Y translational velocity  
";

%feature("docstring") postpro_CleanMemory "

Free all memory allocated within postpro module.  

python usage : postpro_CleanMemory()  
";

%feature("docstring") postpro_2D_GetKineticEnergy "

Compute Kinetic Energy for all bodies (rigids and defo)  

python usage : KE = postpro_2D_GetKineticEnergy()  
";


// File: wrap__RBDY2_8h.xml

%feature("docstring") RBDY2_PutBodyInvMass "

Set inv mass diagonal matrix of a given body. Overwrites the computed values.  

usage : RBDY2_PutBodyInvMass(ibdyty, inv_mass)  

Parameters
----------
ibdyty(integer) : rank of RBDY2  
inv_mass(double array) : inv_mass of RBDY2 (size 3)  
";

%feature("docstring") RBDY2_PutBodyPreconW "

Put preconW of a given body.  

usage : RBDY2_PutBodyPreconW(ibdyty, idof, W)  

Parameters
----------
ibdyty(integer) : rank of RBDY2  
idof(integer) : corresponding dof to set  
W(double array) :  
";

%feature("docstring") RBDY2_PutBodyVector "

Set a vector of a given body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coorb\": coordinates at beginning of time step  
*   \"Coor_\": coordinates in computed configuration  
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

Uses copy, and in case fo Fext, the operation is not just setting but adding  

usage : RBDY2_PutBodyVector(datatype, ibdyty, vector)  

Parameters
----------
datatype(string of size 5) : the vector to set  
ibdyty(integer) : rank of body  
vector(double array) : the new value  
";

%feature("docstring") RBDY2_PutAllBodyVector "

Put an array of a vector of all RBDY2 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY2_PutBodyVector  

python usage : RBDY2_PutAllBodyVector(datatype, matrix)  

Parameters
----------
datatype(string [5]) : the vector to set  
matrix(double array) : input matrix  
";

%feature("docstring") RBDY2_GetBodyVector "

Get a copy of a vector of a given RBDY2 body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Coorb\": coordinates at beginning of time step  
*   \"Coorm\": coordinates in detection configuration  
*   \"Coor_\": coordinates in computed configuration  
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

usage : vector = RBDY2_GetBodyVector(datatype, ibdyty)  

Parameters
----------
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  
vector(double array) : the desired vector  
";

%feature("docstring") RBDY2_GetAllBodyVector "

Get an array of a vector of all RBDY2 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY2_GetBodyVector  

python usage : matrix = RBDY2_GetBodyVector(datatype, ibdyty)  

Parameters
----------
datatype(string [5]) : the vector to get  

Returns
-------
matrix (double array) : output matrix  
";

%feature("docstring") RBDY2_GetPtrBodyVector "

Get a pointer on a vector of a given RBDY2 body.  

Possible values for datatype field are:  

*   \"Coor0\": reference coordinates  
*   \"Xbeg_\": cumulated displacements over time at beginning of time step  
*   \"X____\": cumulated displacements over time in computed configuration  
*   \"Vbeg_\": velocity at beginning of time step  
*   \"V____\": velocity in computed configuration  
*   \"Vaux_\": working array for velocity  
*   \"Ireac\": contact impulse  
*   \"Iaux_\": working array for impulste  
*   \"Fext_\": external forces  

usage : vector_ptr = RBDY2_GetPtrVector(datatype, ibdyty)  

Parameters
----------
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  
vector_ptr(double array) : reference on the desired vector viewed as a numpy array  
";

%feature("docstring") RBDY2_GetBodyInertia "

Get the inertia of a given RBDY2 body.  

usage : inertia = RBDY2_GetBodyInertia(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  
inertia(double) : the inertia of desired body  
";

%feature("docstring") RBDY2_GetAllInertia "

Get the inertia of a all RBDY2 body.  

usage : inertia = RBDY2_GetAllInertia()  

Parameters
----------
inertia(double array): the inertia of all bodies  
";

%feature("docstring") RBDY2_IncrementStep "

increment values at the current time step (prediction)  

usage : RBDY2_IncrementStep()  
";

%feature("docstring") RBDY2_SetVlocyDrivenDof "

Override the value of an existing velocity driven dof ; use after IncrementStep.  

usage : RBDY2_SetVlocyDrivenDof(ibdyty, idrvdof, value)  

Parameters
----------
ibdyty(integer) : rank of considered  
idrvdof(integer) : index of velocity driven dof to set  
value(real) : new value of the velocity driven dof  
";

%feature("docstring") RBDY2_ComputeDof "

Compute current DOF of bodies in container.  

usage : RBDY2_ComputeDof()  
";

%feature("docstring") RBDY2_UpdateDof "

set current DOF as initial DOF of bodies in container  

usage : RBDY2_UpdateDof()  
";

%feature("docstring") RBDY2_ComputeFreeVelocity "

Compute free velocity of bodies in container.  

usage : RBDY2_ComputeFreeVelocity()  
";

%feature("docstring") RBDY2_ComputeFext "

Compute impulse of external forces of bodies in container.  

usage : RBDY2_ComputeFext()  
";

%feature("docstring") RBDY2_ComputeBulk "

Compute impulse of internal forces of bodies in container.  

usage : RBDY2_ComputeBulk()  
";

%feature("docstring") RBDY2_CheckEquilibriumState "

check if all the RBDY2 rich an equilibrium state (velocity is almost equal to
zero)  

usage : isBalanced = RBDY2_CheckEquilibriumState()  

Returns
-------
isBalanced (boolean) : True if in equilibrium state  
";

%feature("docstring") RBDY2_GhostToInvisible "

set bodies with ghost behaviour nickname as invisible  

usage : RBDY2_GhostToInvisible()  
";

%feature("docstring") RBDY2_FatalDamping "

Nullify body current and initial velocities of a list of bodies.  

This keyword must be between the ComputeDof and UpdateDof ones.  

usage : RBDY2_FatalDamping(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to reset current velocity if omitted
    works on all objetcs  
";

%feature("docstring") RBDY2_PartialDamping "

Limit body velocity to Vmax value.  

usage : RBDY2_PartialDamping(nb_steps, Vmax)  

Parameters
----------
nb_steps(integer) : periodicity @parma[in] Vmax (double) : Vmax  
";

%feature("docstring") RBDY2_WriteLastDof "

Write ascii DOF.LAST file.  

usage : RBDY2_WriteLastDof()  
";

%feature("docstring") RBDY2_WriteOutDof "

Write ascii DOF.OUT file. Can be activated only each N steps.  

If 0 for ivalue1 and ivalue2, dofs of all bodies are written.  

usage : RBDY2_WriteOutDof(ivalue1=0, ivalue2=0)  

Parameters
----------
ivalue1(integer) : first body  
ivalue2(integer) : last body  
";

%feature("docstring") RBDY2_DisplayOutDof "

Display bodies degrees of freedom.  

usage : RBDY2_DisplayOutDof()  
";

%feature("docstring") RBDY2_WriteLastRnod "

Write ascii Rnod.LAST file.  

usage : RBDY2_WriteLastRnod()  
";

%feature("docstring") RBDY2_WriteOutRnod "

Write ascii Rnod.OUT file. Can be activated only each N steps.  

usage : RBDY2_WriteOutRnod()  
";

%feature("docstring") RBDY2_DisplayOutRnod "

display body forces  

usage : RBDY2_DisplayOutRnod()  
";

%feature("docstring") RBDY2_WriteBodies "

Write BODIES.OUT file.  

Write DRV_DOF.OUT file.  

usage : RBDY2_WriteBodies()  
";

%feature("docstring") RBDY2_ClearedWriteBodies "

...  

usage : RBDY2_ClearedWriteBodies()  
";

%feature("docstring") RBDY2_WriteDrivenDof "
";

%feature("docstring") RBDY2_ReadBodies "

Read BODIES.DAT file.  

usage : RBDY2_ReadBodies()  

  
 Initialize existing_entities variable in RBDY2  
 Adds the number of found bodies to entity  
";

%feature("docstring") RBDY2_ReadIniDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read Else : OUTBOX/DOF.OUT.num is read, num
being the parameter used in TimeEvolution_ReadIniDof last call  

usage : RBDY2_ReadIniDof(num=0)  

Parameters
----------
num(integer) : which DOF file to read  
";

%feature("docstring") RBDY2_ReadDrivenDof "

Read DRV_DOF.DAT file.  

usage : RBDY2_ReadDrivenDof()  
";

%feature("docstring") RBDY2_LoadBehaviours "

Load bulk behaviour id from bulk_behav module.  

usage : RBDY2_LoadBehaviours()  
";

%feature("docstring") RBDY2_MP_LoadBehaviours "

Load extra physical behaviour read in BULK_BEHAV.DAT file.  

Must be used with THERMO_RIGID ELECTRO_RIGID THERMO_ELECTRO_RIGID behaviours  

usage : RBDY2_MP_LoadBehaviours(disper)  

Parameters
----------
disper(double) : dispersion variable  
";

%feature("docstring") RBDY2_UpdateWSvsT "

update surface energy with temperature  

Must be used with THERMO_RIGID ELECTRO_RIGID THERMO_ELECTRO_RIGID behaviours  

usage : RBDY2_UpdateWSvsT()  
";

%feature("docstring") RBDY2_UpdateWSvsTime "

update surface energy with time  

Must be used with THERMO_RIGID ELECTRO_RIGID THERMO_ELECTRO_RIGID behaviours  

usage : RBDY2_UpdateWSvsTime()  
";

%feature("docstring") RBDY2_ComputeMass "

Compute mass and inertia of bodies.  

usage : RBDY2_ComputeMass()  
";

%feature("docstring") RBDY2_SetPeriodicCondition "

define the space is X periodic [0,periode]  

The X variable reaches a value between 0 and periode  

usage : RBDY2_SetPeriodicCondition(periode)  

Parameters
----------
period(double) : periode  
";

%feature("docstring") RBDY2_ResizeBodies "

resize body radius by a factor  

usage : RBDY2_ResizeBodies(homo)  

Parameters
----------
homo(double) : resize factor  
";

%feature("docstring") RBDY2_NullifyDisplacements "

Set displacements equal to 0.  

usage : RBDY2_NullifyDisplacements()  
";

%feature("docstring") RBDY2_NullifyVelocities "

Set velocity to 0.  

usage : RBDY2_NullifyVelocities()  
";

%feature("docstring") RBDY2_SetSourcePoint "

Create an assembly by source point deposit.  

usage : RBDY2_SetSourcePoint(ibdyty, radius, x_coor, y_coor)  

Parameters
----------
ibdyty(integer): rank of first invisible body  
radius(double) : radius of source point area  
x_coor(double) : X translation from the set of grains  
y_coor(double) : Y translation from the set of grains  
";

%feature("docstring") RBDY2_CheckSourcePoint "

check if it possible to deposit a new particle  

usage : RBDY2_CheckSourcePoint()  
";

%feature("docstring") RBDY2_MembraneBiaxialLoading "

Biaxial load of a sample using pseudo membrane.  

usage : RBDY2_MembraneBiaxialLoading(down, up, thickness, stress)  

Parameters
----------
down(integer) : rank of the lower body  
up(integer) : rank of the upper body  
thickness(double) : thickness of the membrane  
stress(double) : pressure on the membrane  
";

%feature("docstring") RBDY2_BiaxialLoading "

Biaxial load of a sample using a rigid box.  

usage : RBDY2_BiaxialLoading(down,f_down,right,f_right,up,f_up,left,f_left)  

Parameters
----------
down(integer) : rank of the lower body  
f_down(double) : pressure on the lower body  
right(integer) : rank of the right body  
f_right(double) : pressure on the right body  
up(integer) : rank of the upper body  
f_up(double) : pressure on the upper body  
left(integer) : rank of the left body  
f_left(double) : pressure on the left body  
";

%feature("docstring") RBDY2_SetYminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetYminBoundary(inf_boundary)  

Parameters
----------
inf_boundary(double) : inferior boundary value  
";

%feature("docstring") RBDY2_SetYmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetYmaxBoundary(up_boundary)  

Parameters
----------
up_boundary(double) : superior boundary value  
";

%feature("docstring") RBDY2_SetXminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetXminBoundary(left_boundary)  

Parameters
----------
left_boundary(double) : left boundary value  
";

%feature("docstring") RBDY2_SetXmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetXmaxBoundary(right_boundary)  

Parameters
----------
right_boundary(double) : right boundary value  
";

%feature("docstring") RBDY2_SetEquilibriumNorm "

Initialization of data for the equilibrium state check.  

You must precise the type of check test :  

*   Qvlcy : quadratic norm velocy  
*   Mvlcy : maximum norm velocy  

usage : RBDY2_CheckEquilibrium(norm_type , tolerance)  

Parameters
----------
norm_type(string of size 5) : norm type use for the equilibrium check  
tolerance(double) : norm tolerance  
";

%feature("docstring") RBDY2_AddDof2InBodies "

Create a new BODIES.OUT file as combination of the last one and of the last
DOF.OUT file.  

usage : RBDY2_AddDof2InBodies()  
";

%feature("docstring") RBDY2_InitFreeBoundary "

usage : RBDY2_InitFreeBoundary(xmin, xmax, radius)  

Parameters
----------
xmin(double) :  
xmax(double) :  
radius(double) :  
";

%feature("docstring") RBDY2_UpdateThermalStrain "

usage : RBDY2_UpdateThermalStrain()  
";

%feature("docstring") RBDY2_GetNbRBDY2 "

Get the number of RBDY2.  

usage : nb_rbdy2 = RBDY2_GetNbRBDY2()  

Parameters
----------
nb_rbdy2(integer) : number of RBDY2 in container  
";

%feature("docstring") RBDY2_GetBodyArea "

Get the area (2D volume equivalent) of a given body.  

usage : area = GetBodyArea(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the body  
area(double) : area  
";

%feature("docstring") RBDY2_GetAllArea "

Get the area of a all body (visible and invisible)  

python usage : area = RBDY2_GetAllArea()  

Returns
-------
area (double array) : masses of all RBDY2  
";

%feature("docstring") RBDY2_ComputePartialEquilibriumState "

Compute norms used to check if a part of the sample is in a equilibrium state.  

Compute norms used to test if there is an equilibrium state between abs_min and
abs_max  

Usefull in case of silos to access norms used to decide if the arch research
must be activated  

usage : Qnorm, Mnorm = RBDY2_CheckPartialEquilibriumState(abs_min, abs_max)  

Parameters
----------
abs_min(double) : min abscisse of sub domaine tested  
abs_max(double) : max abscisse of sub domaine tested  
Qnorm(double) : quadratric norm of the velocities of the grains  
Mnorm(double) : quadratric norm of the velocities of the grains  
";

%feature("docstring") RBDY2_CheckPartialEquilibriumState "

Check if a part of the sample is in a equilibrium state.  

Test if there is an equilibrium state between abs_min and abs_max  

Usefull in case of silos to decide if the arch research must be activated  

usage : isPartiallyEquilibriumed = RBDY2_CheckPartialEquilibriumState(abs_min,
abs_max)  

Parameters
----------
abs_min(double) : min abscisse of sub domaine tested  
abs_max(double) : max abscisse of sub domaine tested  
isPartiallyEquilibriumed(boolean) : true if at partial equlibrium state, else false  
";

%feature("docstring") RBDY2_SetBodiesInvisible "

Set a list of body to invisible state.  

usage : RBDY2_SetBodiesInvisible(list_bdy)  

Parameters
----------
list_bdy(integer array) : list of rank of bodies of the container  
";

%feature("docstring") RBDY2_IsVisible "

return if a body visible  

usage : visible = RBDY2_IsVisible(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of body  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") RBDY2_GetBodyMass "

Get the mass of a body.  

usage : mass = RBDY2_GetBodyMass(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of desired body  
mass(double) : mass of body  
";

%feature("docstring") RBDY2_GetAllMass "

Get the mass of a all body (visible and invisible)  

python usage : masses = RBDY2_GetAllMass()  

Returns
-------
masses (double array) : masses of all RBDY2  
";

%feature("docstring") RBDY2_CompCoor "

Compute the position of bodies.  

usage : RBDY2_CompCoor()  
";

%feature("docstring") RBDY2_GetBodyDensity "

Get the density of a body.  

usage density = RBDY2_GetBodyDensity(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY2 in container  
density(double) : density of the RBDY2  
";

%feature("docstring") RBDY2_GetNbContactor "

get the number of contactor of RBDY2  

python usage : nb = RBDY2_GetNbContactor(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY2 in container  

Returns
-------
nb (integer) : number of contactor attached to a RBDY2  
";

%feature("docstring") RBDY2_GetContactorType "

Get the type of the first contactor of a body.  

usage type = RBDY2_GetContactorType(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY2 in container  

Returns
-------
type (string) : type of the first contactor of the body  
";

%feature("docstring") RBDY2_GetContactorColor "

Get the color of the itacty contactor of a body ibdyty.  

usage color = RBDY2_GetContactorColor(ibdyty,itacty)  

Parameters
----------
ibdyty(integer) : rank of the RBDY2 in container  
itacty(integer) : rank of the contactor in the RBDY2  

Returns
-------
color (string) : color of the contactor of the body  
";

%feature("docstring") RBDY2_SetContactorColor "

Set the color of a given contactor of a body.  

usage : RBDY2_SetContactorColor(ibdyty, itacty, color)  

Parameters
----------
ibdyty(integer) : rank of the RBDY2  
itacty(integer) : rank of the contactor in the RBDY2  
color(string of size 5) : the color  
";

%feature("docstring") RBDY2_GetPtrMass "

Get a pointer onto the mass matrix of a body.  

Parameters
----------
ibdyty(int): index of the RBDY2  
mass(double**): mass matrix of the RBDY2  
";

%feature("docstring") RBDY2_GetVelocity "

Get the velocity of a body.  

Parameters
----------
ibdyty(int): index of the RBDY2  
velocity(double[6]): velocity of the RBDY2  
";

%feature("docstring") RBDY2_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = RBDY2_getDrvVlocy(ibdyty)  

Parameters
----------
ibdyty(integer) : index of the RBDY2  
drvdof_indices(integer array) : indices list of driven dof  
drvdof_values(real array) : values of the driven dof  
";

%feature("docstring") RBDY2_computeDrvVlocy "

Compute the value of the driven velocity of a body a current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : RBDY2_computeDrvVlocy(ibdyty, values)  

Parameters
----------
ibdyty(integer) : index of the RBDY2  
values(double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") RBDY2_SetVisible "

Set a given body as visible.  

python usage : RBDY2_SetVisible(ibdyt)  

Parameters
----------
ibdyty(integer) : index of the RBDY2  
";

%feature("docstring") RBDY2_SetInvisible "

Set a given body as invisible.  

python usage : RBDY2_SetInvisible(ibdyty)  

Parameters
----------
ibdyty(integer) : index of the RBDY2  
";

%feature("docstring") RBDY2_SetVisibleVlocyDrivenDof "

allows to (re)activate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : RBDY2_SetVisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
ibdyty(integer): index of the RBDY2  
iccdof(integer): index of the DOF to set visible  
";

%feature("docstring") RBDY2_SetInvisibleVlocyDrivenDof "

allows to deactivate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : RBDY2_SetInvisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
ibdyty(integer): index of the RBDY2  
iccdof(integer): index of the DOF to set invisible  
";

%feature("docstring") RBDY2_GetBulkBehavID "

return the ID of a given bulk of a given body  

python usage : blmID = DISKx_GetBulkBehavID(ibdyty, iblmty)  

Parameters
----------
ibdyty(integer) : rank of a RBDY2  
iblmty(integer) : rank of a bulk of the giben RBDY2 (typically 1!)  

Returns
-------
blmID (string) : the bulk behav ID  
";

%feature("docstring") RBDY2_GetBulkBehavNumber "

return the bulk ID of a given RBDY2  

python usage : ibehav = RBDY2_GetBulkBehavNumber(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of a RBDY2  

Returns
-------
ibehav (integer) : the bulk behav number  
";

%feature("docstring") RBDY2_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : RBDY2_SetSurfaceSectors(nbsect)  

Parameters
----------
nbsect(integer) : number of sectors  
";

%feature("docstring") RBDY2_GetStress "

Get the mean stress field of a rigid object.  

python usage : matrix = RBDY2_GetStress(ibdyty)  

Parameters
----------
ibdyty(integer) : body to get stress of  
matrix(double array) : stress matrix  
";

%feature("docstring") RBDY2_ModifyBody "

Modify a body tactor.  

usage : RBDY2_ModifyBody(ibdyty, itacty, vector)  

Parameters
----------
ibdyty(integer) : rank of body  
itacty(integer) : rank of tacty  
vector(double array) : the new value  
";

%feature("docstring") RBDY2_SkipInvisible "

skip invisible objects when writing BODIES.OUT  

usage : RBDY2_SkipInvisible()  
";

%feature("docstring") RBDY2_InitializeStresses "

initialize stress for rigid bodies  

usage : RBDY2_InitializeStresses()  
";

%feature("docstring") RBDY2_InitializeWS "

initialize WS for rigid bodies with a value between wsmin and wsmax ponderate by
rvalue1  

python usage : RBDY2_InitializeWS(double rvalue1)  
";

%feature("docstring") RBDY2_CleanMemory "

Free all memory allocated within RBDY2 module.  

python usage : RBDY2_CleanMemory()  
";

%feature("docstring") RBDY2_GetThermalValue "

Get temperature of rigid particle.  

usage : T = RBDY2_GetThermalValu(ibdyty, itacty)  

Parameters
----------
ibdyty(integer) : rank of body  
itacty(integer) : rank of tacty  
";

%feature("docstring") RBDY2_GetElectricalPotential "

Get electrical potential of rigid particle.  

usage : ep = RBDY2_GetElectricalPotential(ibdyty, itacty)  

Parameters
----------
ibdyty(integer) : rank of body  

Returns
-------
ep (double) : electrical potential  
";

%feature("docstring") RBDY2_GetElectricalCurrent "

Get electrical potential of rigid particle.  

usage : ep = RBDY2_GetElectricalCurrent(ibdyty, itacty)  

Parameters
----------
ibdyty(integer) : rank of body  

Returns
-------
ep (double) : electrical current  
";

%feature("docstring") RBDY2_GetBetai "

Get equivalent damage related to CZM interaction for rigid particle.  

usage : betai = RBDY2_GetBetai(ibdyty, itacty)  

Parameters
----------
ibdyty(integer) : rank of body  
itacty(integer) : rank of tacty  

Returns
-------
betai (double) : equivalent damage  
";

%feature("docstring") RBDY2_GetPeriode "

Get the periode id (0, 1 or -1) for rigid particles.  

usage : iper = RBDY2_GetPeriode(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of body  

Returns
-------
iper (integer) : periode id  
";

%feature("docstring") RBDY2_GetAverageSurfaceEnergy "
";


// File: wrap__PT2Dx_8h.xml

%feature("docstring") PT2Dx_LoadTactors "

load PT2Dx from RBDY2 and initialize existing_entites  

python usage : PT2Dx_LoadTactors()  
";

%feature("docstring") PT2Dx_GetNbPT2Dx "

Get the number of PT2Dx in the container.  

python usage : nb_pt2d = PT2Dx_GetNbPT2Dx()  

Returns
-------
nb_pt2d (integer) : the number of PT2Dx in container  
";

%feature("docstring") PT2Dx_SetDisplayRadius "

Set a radius to display a pt2dx.  

python usage : PT2Dx_SetDisplayRadius(radius)  

Parameters
----------
radius(double) : value of the radius which should be used for display  
";

%feature("docstring") PT2Dx_GetPtrPT2Dx2BDYTY "

return a pointer onto the map pt2dx2rbdy2  

python usage : ptd2x2rbdy2 = PT2Dx_GetPtrPT2Dx2BDYTY()  

Returns
-------
pt2dx2rbdy2 (integer array) : reference on map between pt2dx rank and body/tact
rank  
";

%feature("docstring") PT2Dx_IsVisible "

return if a body visible  

usage : visible = PT2Dx_IsVisible(itact)  

Parameters
----------
itact(integer) : rank of PT2Dx  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PT2Dx_InitOutlines "

Get a reference on the outlines of all PT2Dx.  

usage : outlines = PT2Dx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PT2Dx  
";

%feature("docstring") PT2Dx_InitScalarFields "

Get a reference on the scalar fields of all PT2Dx.  

usage : scalarfields = PT2Dx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PT2Dx array  
";

%feature("docstring") PT2Dx_UpdatePostdata "

Update values of outlines_PT2Dx and scalarfields_PT2Dx pointers.  

usage : PT2Dx_UpdatePostdata  
";

%feature("docstring") PT2Dx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PT2Dx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PT2Dx  
";

%feature("docstring") PT2Dx_GetNbScalarFields "

Get the number of scalar fields of a PT2Dx.  

python usage : nb_scalarfields = PT2Dx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PT2Dx  
";

%feature("docstring") PT2Dx_CleanMemory "

Free all memory allocated within PT2Dx module.  

python usage : PT2Dx_CleanMemory()  
";


// File: wrap__MAILx_8h.xml

%feature("docstring") MAILx_ReadBodies "

read MAILx from DATBOX/BODIES.DAT  

Input string is of form vX.Y where X is major version number and Y is minor one.  
 If not specified, last available version is used.  

python usage : MAILx_ReadBodies(version)  

param[in] version (string) : file format version to use  
";

%feature("docstring") MAILx_WriteBodies "

write MAILx to OUTBOX/BODIES.OUT  

Input string is of form vX.Y where X is major version number and Y is minor one.  
 If not specified, last available version is used.  

python usage : MAILx_WriteBodies(version)  

param[in] version (string) : file format version to use  
";

%feature("docstring") MAILx_WriteLastGPV "

write OUTBOX/GPV.LAST  

python usage : MAILx_WriteLastGPV()  
";

%feature("docstring") MAILx_WriteOutGPV "

write OUTBOX/GPV.OUT.x  

python usage : MAILx_WriteOutGPV()  
";

%feature("docstring") MAILx_DisplayOutGPV "

display GPV values  

python usage : MAILx_DisplayOutGPV()  
";

%feature("docstring") MAILx_AddDof2InBodies "

set cooref = cooref + X  

python usage : MAILx_AddDof2InBodies()  
";

%feature("docstring") MAILx_GetNbMAILx "

Get the number of MAILx.  

python usage : nb_MAILx = GetNbMAILx()  

Returns
-------
nb_MAILx (integer) : number of MAILx  
";

%feature("docstring") MAILx_GetNbCell "

Get the number of Cells of a given MAILx.  

python usage : nb_MAILx = GetNbCell(IdBody)  

Parameters
----------
IdBody(integer) : id of the concern body  

Returns
-------
nb_cell (integer) : number of cell  
";

%feature("docstring") MAILx_SetCoorRef "

set reference coordinates on a given body  

python usage : MAILx_SetCoorRef(IdBody, f, length)  

Parameters
----------
IdBody(integer) : id of the concern body  
f(double array) : value of the vitesse  
length(integer) : length of vector  
";

%feature("docstring") MAILx_GetCoordNodty "

Get one coordinate of a node of a body.  

python usage : x = MAILx_GetCoordNodty(int ibdty,int inodty,int icomp)  

Parameters
----------
ibdyty(integer) : rank of considered body  
inodty(integer) : the node  
icomp(integer) : the component  

Returns
-------
x (double) : coordinate of node  
";

%feature("docstring") MAILx_GetCoordsNodty "

Get the coordinates of a node of a body.  

python usage : x = MAILx_GetCoordsNodty(int ibdty, int inodty, int length)  

Parameters
----------
ibdyty(integer) : rank of considered body  
inodty(integer) : the node  
length(integer) : the number of component  

Returns
-------
x (double array) : the desired vector  
";

%feature("docstring") MAILx_GetNbNodes "

Get the number of nodes of a given MAILx.  

python usage : nb_nodes = MAILx_GetNbNodes(ibdyty)  

Parameters
----------
ibdyty(integer) : body id  

Returns
-------
nb_nodes (integer) : number of nodes of the body  
";

%feature("docstring") MAILx_InitNodalFields "

Set the number of nodal_fields for a given body.  

python usage : MAILx_InitNodalFields(ibdyty,nb_nodal_fields)  

Parameters
----------
ibdyty(integer) : body id  
nb_nodal_fields(integer) : number of nodal fields required  
";

%feature("docstring") MAILx_InitNodalField "

Set name and size of a nodal_field of a given body.  

python usage : MAILx_InitNodalField(ibdyty,name,rank,sz)  

Parameters
----------
ibdyty(integer) : body id  
name(string) : field name  
rank(integer) : field rank  
sz(integer) : size of the field  
";

%feature("docstring") MAILx_SetNodalField "

Set a nodal_field of a given body.  

python usage : MAILx_SetNodalField(ibdyty,rank,field)  

Parameters
----------
ibdyty(integer) : body id  
rank(integer) : field rank  
field(double vector) : field  
";

%feature("docstring") MAILx_CleanMemory "

Free all memory allocated within MAILx module.  

python usage : MAILx_CleanMemory()  
";


// File: wrap__PT3Dx_8h.xml

%feature("docstring") PT3Dx_LoadTactors "

load PT3Dx from RBDY3 and initialize existing_entites  

python usage : PT3Dx_LoadTactors()  
";

%feature("docstring") PT3Dx_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = PT3Dx_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PT3Dx_GetNbPT3Dx "

Get the number of PT3Dx.  

python usage : nb_PT3Dx = PT3Dx_GetNbPT3Dx()  

Returns
-------
nb_PT3Dx (integer) : the number of PT3Dx  
";

%feature("docstring") PT3Dx_SetDisplayRadius "

set the size of the glyph representing the PT3Dx  

python usage : PT3Dx_SetDisplayRadius(radius)  

Parameters
----------
radius(double): radius of the PT3Dx contactors  
";

%feature("docstring") PT3Dx_GetPtrPT3Dx2BDYTY "

return a pointer onto the map pt3dx2bdyty  

python usage : pt3dx2bdyty = PT3Dx_GetPtrPT3Dx2BDYTY()  

Returns
-------
pt3dx2bdyty (integer array) : reference on map between pt3dx rank and body rank  
";

%feature("docstring") PT3Dx_InitOutlines "

Get a reference on the outlines of all PT3Dx.  

usage : outlines = PT3Dx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PT3Dx  
";

%feature("docstring") PT3Dx_InitScalarFields "

Get a reference on the scalar fields of all PT3Dx.  

usage : scalarfields = PT3Dx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PT3Dx array  
";

%feature("docstring") PT3Dx_UpdatePostdata "

Update values of outlines_PT3Dx and scalarfields_PT3Dx pointers.  

usage : PT3Dx_UpdatePostdata  
";

%feature("docstring") PT3Dx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PT3Dx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PT3Dx  
";

%feature("docstring") PT3Dx_GetNbScalarFields "

Get the number of scalar fields of a PT3Dx.  

python usage : nb_scalarfields = PT3Dx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PT3Dx  
";

%feature("docstring") PT3Dx_GetPtrAllConnectivities "

Get a reference on the connectivities of all PT3Dx.  

usage : connec = PT3Dx_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") PT3Dx_CleanMemory "

Free all memory allocated within PT3Dx module.  

python usage : PT3Dx_CleanMemory()  
";


// File: wrap__DISKx_8h.xml

%feature("docstring") DISKx_LoadTactors "

load DISKx from RBDY2 file and initialize existing_entites  

python usage : DISKx_LoadTactors()  
";

%feature("docstring") DISKx_GetNbDISKx "

Get the number of DISKx in the container.  

python usage : nb_diskx = DISKx_GetNbDISKx()  

Returns
-------
nb_DISKx (integer) : the number of DISKx in container  
";

%feature("docstring") DISKx_GetDISKx2BDYTY "

Get a copy of map DISKx2bdyty.  

usage : polyr2bdyty = DISKx_GetDISKx2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") DISKx_GetPtrDISKx2BDYTY "

return a pointer onto the map diskx2rbdy2  

python usage : diskx2bdyty = DISKx_GetPtrDISKx2BDYTY()  

Returns
-------
diskx2bdyty (integer array) : reference on map between diskx rank and body rank  
";

%feature("docstring") DISKx_IsVisible "

return if a body visible  

python usage : visible = DISKx_IsVisible(itact)  

Parameters
----------
itact(integer) : rank of DISKx  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") DISKx_GetContactorRadius "

Get the radius of a given DISKx.  

python usage : radius = DISKx_GetContactorRadius(itact)  

Parameters
----------
itact(integer) : rank of a DISKx (in the list of all the DISKx)  

Returns
-------
radius (double) : the radius of the DISKx of rank itact  
";

%feature("docstring") DISKx_GetMeanRadius "

Get the mean radius of DISKx in the container.  

python usage : radius = DISKx_GetMeanRadius()  

Returns
-------
radius (double) : the mean radius of DISKx in the container  
";

%feature("docstring") DISKx_GetMaxRadius "

Get the max radius of DISKx in the container.  

python usage : radius = DISKx_GetMaxRadius()  

Returns
-------
radius (double) : the max radius of DISKx in the contactor  
";

%feature("docstring") DISKx_GetMinRadius "

Get the min radius of DISKx in the container.  

python usage : radius = DISKx_GetMinRadius()  

Returns
-------
radius (double) : the min radius of DISKx in the container  
";

%feature("docstring") DISKx_GetContactorColor "

Get the color of a given DISKx.  

python usage : color = DISKx_GetContactorColor(itact)  

Parameters
----------
itact(integer) : rank of a DISKx  

Returns
-------
color (string) : the color of the DISKx itact  
";

%feature("docstring") DISKx_GetRadius "

get radius of a DISKx  

python usage : radius = DISKx_GetRadius(itacty)  

Parameters
----------
itacty(integer) : rank of DISKx  

Returns
-------
radius (double) : the radius of DISKx of body ibdyty  
";

%feature("docstring") DISKx_GetContactorCoor "

get coordinates of the center of a given DISKx  

python usage : vector = DISKx_GetContactorCoor(itacty)  

Parameters
----------
itacty(integer) : rank of considered contactor  

Returns
-------
vector (double array) : the desired vector  
";

%feature("docstring") DISKx_InitOutlines "

Get a reference on the outlines of all DISKx.  

python usage : outlines = DISKx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_DISKx  
";

%feature("docstring") DISKx_InitScalarFields "

Get a reference on the scalar fields of all DISKx.  

python usage : scalarfields = DISKx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_DISKx array  
";

%feature("docstring") DISKx_UpdatePostdata "

Update values of outlines_DISKx and scalarfields_DISKx pointers.  

python usage : DISKx_UpdatePostdata()  
";

%feature("docstring") DISKx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = DISKx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
DISKx  
";

%feature("docstring") DISKx_GetNbScalarFields "

Get the number of scalar fields of a DISKx.  

python usage : nb_scalarfields = DISKx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a DISKx  
";

%feature("docstring") DISKx_CleanMemory "

Free all memory allocated within DISKx module.  

python usage : DISKx_CleanMemory()  
";

%feature("docstring") DISKx_SetXdilation "

set increase of radius of a DISKx due to expansion  

python usage : DISKx_SetXdilation(itacty,x)  

Parameters
----------
itacty(integer) : rank of considered contactor  
x(float) : increase of radius  
";

%feature("docstring") DISKx_SetVdilation "

set increase rate of radius of a DISKx due to expansion  

python usage : DISKx_SetVdilation(itacty, v)  

Parameters
----------
itacty(integer) : rank of contactor  
v(float) : radius increase rate  
";


// File: wrap__bulk__behav_8h.xml

%feature("docstring") bulk_behav_ReadBehaviours "

read gravity and behaviors from DATBOX/BULK_BEHAV.DAT file  

python usage : bulk_behav_ReadBehaviours()  
";

%feature("docstring") bulk_behav_WriteBehaviours "

write gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_WriteBehaviours()  
";

%feature("docstring") bulk_behav_CollectOutBulkBehav "

read gravity and behaviors from OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_CollectOutBulkBehav()  
";

%feature("docstring") bulk_behav_CleanOutBulkBehav "

write (replacing) gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_CleanOutBulkBehav()  
";

%feature("docstring") bulk_behav_AppendOutBulkBehav "

write (appending) gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_AppendOutBulkBehav()  
";

%feature("docstring") bulk_behav_RebuildInBulkBehav "

write (replace) gravity and behaviors to DATBOX/BULK_BEHAV.DAT file  

python usage : bulk_behav_RebuildInBulkBehav()  
";

%feature("docstring") bulk_behav_GetGravity "

get the gravity acceleration used  

python usage : gravity = bulk_behav_GetGravity()  

Returns
-------
gravity (double array) : gravity vector  
";

%feature("docstring") bulk_behav_SetGravity "

set the gravity acceleration to be used  

python usage : bulk_behav_SetGravity(gravity)  

Parameters
----------
gravity(double array) : gravity vector (size 3)  
";

%feature("docstring") bulk_behav_SetConductivity "

set the conductivity parameter to be used  

python usage : bulk_behav_SetConductivity(cvalue ,ivalue, rvalue)  

Parameters
----------
cvalue(string of size 5) : nickname of bulk behaviour  
ivalue(integer) : type of parameter: 0 = constant, 1 = field  
rvalue(real) : conductivity value  
";

%feature("docstring") bulk_behav_SetCapacity "

set the Capacity parameter to be used  

python usage : bulk_behav_SetCapacity(cvalue ,ivalue, rvalue)  

Parameters
----------
cvalue(string of size 5) : nickname of bulk behaviour  
ivalue(integer) : type of parameter: 0 = constant, 1 = field  
rvalue(real) : Capacity value  
";

%feature("docstring") bulk_behav_SetBiot "

set the Biot parameter to be used  

python usage : bulk_behav_SetBiot(cvalue ,ivalue, rvalue)  

Parameters
----------
cvalue(string of size 5) : nickname of bulk behaviour  
ivalue(integer) : type of parameter: 0 = constant, 1 = field  
rvalue(real) : Biot value  
";

%feature("docstring") bulk_behav_SetExternalFlux "

set the External Flux parameter to be used  

python usage : bulk_behav_SetExternalFlux(cvalue ,ivalue, rvalue)  

Parameters
----------
cvalue(string of size 5) : nickname of bulk behaviour  
ivalue(integer) : type of parameter: 0 = constant, 1 = field  
rvalue(real) : External Flux value  
";

%feature("docstring") bulk_behav_SetDensity "

set the Density parameter to be used  

python usage : bulk_behav_SetDensity(cvalue , rvalue)  

Parameters
----------
cvalue(string of size 5) : nickname of bulk behaviour  
rvalue(real) : Density value  
";

%feature("docstring") bulk_behav_GetNbBulkBehav "

get the number of bulk laws  

python usage : nb_bulk_behav = bulk_behav_GetNbBulkBehav()  

Parameters
----------
nb_bulk_behav(integer) : number of bulk behaviour in lmgc90  
";

%feature("docstring") bulk_behav_GetBulkBehav "

get a given bulk law  

python usage : lawty, behav = bulk_behav_GetBulkBehav(i_bb)  

Parameters
----------
i_bb(integer) : index of the desired bulk_behav  
lawty(string) : type of the bulk law  
behav(string) : name of the bulk law  
param(real vector) : parameters of the law  
";

%feature("docstring") bulk_behav_CleanMemory "

Free all memory allocated within bulk_behav module.  

python usage : bulk_behav_CleanMemory()  
";


// File: wrap__DKDKL_8h.xml

%feature("docstring") DKDKL_SelectProxTactors "

contact detection between DISKx and DISKL tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKDKL_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKDKL_SmoothForceComputation "

explicit computation of contact forces  

python usage : DKDKL_SmoothForceComputation  
";

%feature("docstring") DKDKL_WriteLastVlocRloc "

write last local values of all DKDKL contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKL_WriteLastVlocRloc()  

The values written are relative velocity, forces and local frame  
";

%feature("docstring") DKDKL_WriteOutVlocRloc "

write local values of all DKDKL contacts  

The values written are relative velocity, forces and local frame  

python usage : DKDKL_WriteOutVlocRloc()  
";

%feature("docstring") DKDKL_DisplayOutVlocRloc "

display local values of all DKDKL contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKDKL_DisplayOutVlocRloc()  
";

%feature("docstring") DKDKL_DisplayProxTactors "

display contacts  

python usage : DKDKL_DisplayProxTactors()  
";

%feature("docstring") DKDKL_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : DKDKL_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKDKL_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKDKL_SetPeriodicCondition(period)  

Parameters
----------
period(double) : value of the period  
";

%feature("docstring") DKDKL_CleanMemory "

Free all memory allocated within DKDKL module.  

python usage : DKDKL_CleanMemory()  
";


// File: wrap__nlgs_8h.xml

%feature("docstring") nlgs_ExPrep "

Prepare matrix storage.  

python usage : nlgs_ExPrep(storage)  

Parameters
----------
sotrage(char[30]) : matrix storage  
  
 prepare the matrix and the RHS of the contact problem in regards of the
selected matrix storage:  

*   Exchange_Local_Global (the standard case) only the diagonal blocks are
    computed and stored.  
*   Stored_Delassus_Loops (faster but memory expensive) the complete Delassus
    matrix is computed.  
";

%feature("docstring") nlgs_ExIter "

Execute NLGS iterations over the contact loop.  

python usage : nlgs_ExIter(nb_iter) param[in] nb_iter (integer) : number of
iterations to do  
";

%feature("docstring") nlgs_ExPost "

Run a jacobi iteration with the solution obtained with the NLGS algorithm.  

python usage : nlgs_ExPost()  
";

%feature("docstring") nlgs_AfterIterCheck "

Control NLGS convergence.  

python usage : convergence = nlgs_AfterIterCheck()  

Returns
-------
convergence (integer) :  
";

%feature("docstring") nlgs_DisplayAfterIterCheck "

Display NLGS convergence results.  

python usage : nlgs_DisplayAfterIterCheck()  
";

%feature("docstring") nlgs_NormCheck "

Active one step norm evolution.  

python usage : nlgs_NormCheck()  
";

%feature("docstring") nlgs_UpdateTactBehav "

Update internal parameters of contact lawz for each contact.  

python usage : nlgs_UpdateTactBehav()  
";

%feature("docstring") nlgs_SetCheckType "

Define numerical convergence of the NLGS algorithm.  

python usage : nlgs_SetCheckType(check_type, tolerance, relaxation)  

Parameters
----------
check_type(char[5]) : type of convergence check  
tolerance(double) : norm tolerance  
relaxation(double) : relaxation factor  
  
 convergence check keywords:  
 Quad : quadratic norm (faulty contacts are redeemed by accurate contacts;
laxist norm)  
 Maxm : maximum norm (faulty contacts must comply; severe norm)  
 QM/16 : maximum of Quad and Maxm/16 norms (a compromise). For large dense
collections Quad ranges usually around 1/16 Maxm  
 where Quad,Maxm,QM/16 are keywords for the check test, and the following real
number is the tolerance value.  
";

%feature("docstring") nlgs_ScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_ScrambleContactOrder()  
";

%feature("docstring") nlgs_QuickScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_QuickScrambleContactOrder()  
";

%feature("docstring") nlgs_SetWithQuickScramble "

active quick scramble in macro function ExSolver  

python usage : nlgs_SetWithQuickScramble()  
";

%feature("docstring") nlgs_ReverseContactOrder "

Reverse the numbering of the contact list.  

python usage : nlgs_ReverseContactOrder()  
";

%feature("docstring") nlgs_BimodalContactOrder "

Renumbering of the contact list using the definition of weak and strong network
in granular assemblies.  

python usage : nlgs_BimodalContactOrder()  
";

%feature("docstring") nlgs_ScaleRloc "

Scale all local contact forces of a factor equal to * 0.9 < f < 1.1.  

python usage : nlgs_ScaleRloc()  
";

%feature("docstring") nlgs_ComputeRnod "

mapping from local contact forces to global ones  

python usage : nlgs_ComputeRnod()  
";

%feature("docstring") nlgs_DisplayRlocNSum "

Display the sum of normal contact forces.  

python usage : nlgs_DisplayRlocNSum()  
";

%feature("docstring") nlgs_ExSolver "

Solve fully the local contact problem.  

python usage : nlgs_ExSolver(storage, checktype, tol, relax, nb_iter_check,
nb_block_iter)  

Parameters
----------
storage(char[30]) : matrix storage (cf nlgs_ExPrep)  
checktype(char[5]) : convergentce test keyword  
tolerance(double) : tolerance value  
relaxation(double) : relaxation number  
nb_iter_check(integer) : number of iteration between convergence test  
nb_block_iter(integer) : number of block iterations  
";

%feature("docstring") nlgs_UpdateCohesiveBehav "

update internal parameters of contact laws for each contact  

python usage : nlgs_UpdateCohesiveBehav(void)  
";

%feature("docstring") nlgs_UpdateFrictionalBehav "

update internal parameters of contact laws for each contact  

python usage : nlgs_UpdateFrictionalBehav(void)  
";

%feature("docstring") nlgs_GetAllThis "

Get all interactions in \"this\" array.  

Each interaction has (in this order): coor, tuc, nuc, rlt, rln, vlt, vln  

usage : interactions = nlgs_GetAllThis()  

Returns
-------
interactions (double 2D-array) : the interactions  
";

%feature("docstring") nlgs_UseJacobiSolver "

Use a Jacobi solver instead of Gauss Seidel solver.  

usage : nlgs_UseJacobiSolver(True) or nlgs_UseJacobiSolver(False)  
";

%feature("docstring") nlgs_UseRegularization "

use some regularization heuristics on interaction laws  

python usage : nlgs_UseRegularization(krn, krt)  

Parameters
----------
krn(double) : normal penality (default 1e14)  
krt(double) : tangential penality (default 1e14)  
";

%feature("docstring") nlgs_SetTemporaryVariable "

set temporary variables used in nlgs ; ivalue2 == 3 gives access to post crack
pressure  

python usage : nlgs_SetTemporaryVariable(icdan,id,val)  

Parameters
----------
icdan(int) : interaction rank  
id(int) : value rank  
val(double) : value  
";

%feature("docstring") nlgs_GetTemporaryVariable "

get temporary variables used in nlgs ; ivalue2 == 3 gives access to post crack
pressure  

python usage : val = nlgs_GetTemporaryVariable(icdan,id)  

Parameters
----------
icdan(int) : interaction rank  
id(int) : value rank  
val(double) : value  
";

%feature("docstring") nlgs_IsInitialized "

In case of restart say that nlgs is initialized or reset it.  

python usage : nlgs_IsInitialized(is_init=1)  
";


// File: wrap__timer_8h.xml

%feature("docstring") timer_InitializeTimers "

Set all timers to 0.  

python usage : timer_InitializeTimers()  
";

%feature("docstring") timer_WriteOutTimers "

write the cumulated times of all the timers  

python usage : timer_WriteOutTimers()  
";

%feature("docstring") timer_GetNewTimer "

create a new timer  

python usage : id = timer_GetNewTimer(name)  

Parameters
----------
name(string) : name of new timer  

Returns
-------
id (integer) : id of the timer created  
";

%feature("docstring") timer_StartTimer "

start a given timer  

python usage : timer_StartTimer(timer_id)  

Parameters
----------
timer_id(integer) : id of the timer to start  
";

%feature("docstring") timer_StopTimer "

stop a given timer, and add the elapsed time since start to the time  

python usage : timer_StopTimer(timer_id)  

Parameters
----------
timer_id(integer) : id of the timer to stop  
";

%feature("docstring") timer_ClearAll "

clear all timers (internal, external and user)  

python usage : timer_ClearAll()  
";


// File: wrap__DKDKx_8h.xml

%feature("docstring") DKDKx_SelectProxTactors "

contact detection between CLxxx and JCxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKDKx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
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
nb(integer) : number of contact points for a couple (cd,an)  
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
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKDKx_SetPeriodicCondition "

initialize data for simulation using periodic condition  

python usage : DKDKx_SetPeriodicCondition(period)  

Parameters
----------
period(double) : value of the period  
";

%feature("docstring") DKDKx_SetFrictionModel "

initialize data for simulation using evolutive local friction  

python usage : DKDKx_SetFrictionModel(cflag)  

Parameters
----------
cflag(char) : model to use ('min', 'max' or 'ave')  
";

%feature("docstring") DKDKx_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : DKDKx_SetSurfaceSectors(nbsect)  

Parameters
----------
nbsect(integer) : number of sectors  
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
icdan(int): index of the DKDKx contact  

Returns
-------
energy(double[4]) : energy value  
";


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
num(integer) : which DOF file to read  
";

%feature("docstring") poroMAILx_ReadIniMecaDof "

Read DOF file.  

If num <= 0 : DATBOX/DOF.INI file is read  

Else : OUTBOX/DOF.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniMecaDof last call  

python usage : poroMAILx_ReadIniMecaDof(num=0)  

Parameters
----------
num(integer) : which DOF file to read  
";

%feature("docstring") poroMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : poroMAILx_ReadIniGPV(num=0)  

Parameters
----------
num(integer) : which GPV file to read  
";

%feature("docstring") poroMAILx_ReadIniMecaGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read Else : OUTBOX/GPV.OUT.num is read, num
being the parameter used in TimeEvolution_ReadIniMecaGPV last call  

python usage : poroMAILx_ReadIniMecaGPV(num=0)  

Parameters
----------
num(integer) : which GPV file to read  
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
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired vector  
";

%feature("docstring") poroMAILx_GetNbNodes "

Get the number of nodes of a poroMAILx.  

python usage : nb_nodes = poroMAILx_GetNbNodes(ibdyty)  

Parameters
----------
ivalue(integer) : id of the poroMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a poroMAILx  
";

%feature("docstring") poroMAILx_GetNbElements "

Get the number of elements of a poroMAILx.  

python usage : nb_elements = poroMAILx_GetNbElements(ibdyty)  

Parameters
----------
ivalue(integer) : id of the poroMAILx  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
  
 You need to set this field in your models.dat  
";

%feature("docstring") poroMAILx_SetTherScalarFieldByNode "

Update an external field on a given body.  

python usage : poroMAILx_SetTherieldByNode(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
  
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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") poroMAILx_SetTherScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetTherScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") poroMAILx_GetMecaScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetMecaScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired scalar field  

Returns
-------
f_rank (integer) : rank of the corresponding scalar field  
";

%feature("docstring") poroMAILx_GetMecaVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetMecaVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired vector field  

Returns
-------
f_rank (integer) : rank of the corresponding vector field  
";

%feature("docstring") poroMAILx_GetTherScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetTherScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired scalar field  

Returns
-------
f_rank (integer) : rank of the corresponding scalar field  
";

%feature("docstring") poroMAILx_GetTherVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = poroMAILx_GetTherVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired vector field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetMecaVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetVectorFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetTherVectorFieldByNode "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetFieldByNode(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") poroMAILx_SetTherVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : poroMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
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
datatype(string of size 5) : the vector to set  
ibdyty(integer) : rank of body  
matrix(double array) : the new values  
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
ibdyty(integer) : rank of considered body  
required_field(integer) : required additional field  

Returns
-------
matrix_out (double 2D-array) : the desired stress  
";

%feature("docstring") poroMAILx_GetStrain "

Get a copy of a strain of a given body.  

Python usage : strain = poroMAILx_GetStrain(ibdyty, required_field=0)  

Parameters
----------
ibdyty(integer) : rank of considered body  
required_field(integer) : required additional field  

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
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") poroMAILx_GetAll "

return poro mechanical data computed for idBody  

python usage : array = poroMAILx_GetAll(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : poro mechanical data  
";

%feature("docstring") poroMAILx_GetGrad "

Get a copy of a grad P of a given body.  

Python usage : grad = poroMAILx_GetGrad(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

Returns
-------
grad (double 2D-array) : the desired grad  
";

%feature("docstring") poroMAILx_GetFlux "

Get a copy of a Darcy Flux of a given body.  

Python usage : flux = poroMAILx_GetFlux(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

Returns
-------
flux (double 2D-array) : the desired flux  
";

%feature("docstring") poroMAILx_GetInternal "

return internal mechanical data computed for idBody  

python usage : array = poroMAILx_GetInternal(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical internal data  
";

%feature("docstring") poroMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = poroMAILx_GetConnectivity(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") poroMAILx_SetVlocyDrivenDof "

Apply Drv Dof on a given body.  

python usage : poroMAILx_SetVlocyDrivenDof(IdBody, f_dof, f_node, f_value)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_dof(integer) : dof of the concern node  
f_node(integer) : node  
f_value(double) : value of the drvdof  
";

%feature("docstring") poroMAILx_AddFieldLoad "

Add elementary load through a nodal external field on a given body.  

python usage : poroMAILx_AddFieldLoad(IdBody, Ideriv, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f(double array) : value of the field  
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


// File: wrap__inter__handler__2D_8h.xml

%feature("docstring") inter_handler_2D_tgetNb "

return the number of interactions of the selected type stored in this data
structure  

python usage : nb_inter = inter_handler_2D_tgetNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_2D_tgetTactLawNb "

return the contact law number of an interaction stored this data structure  

python usage : tact_law = inter_handler_2D_tgetTactLawNb(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
tact_law (integer) : contact law number  
";

%feature("docstring") inter_handler_2D_tgetIdBodies "

return the serial numbers of contacting objects for an interaction stored in
this data structure  

python usage : idBodies = inter_handler_2D_tgetIdBodies(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
idBodies (integer) : array with cd and an bodies serial number  
";

%feature("docstring") inter_handler_2D_tgetIData "

Get the integer data of an interaction stored in this data structure.  

idata vector holds cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

usage : idata = inter_handler_2D_tgetIData(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
idata (integer array) : the values array  
";

%feature("docstring") inter_handler_2D_tgetRData "

return the real data associated with an interactions  

Get an output array with, in this order, : coor, tuc, nuc, vlt, vln, rlt, rln,
gapTT  

python usage : rdata = inter_handler_2D_tgetRData(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
rdata (double array) : array with real data of the interaction  
";

%feature("docstring") inter_handler_2D_tsetInternal "

Set the internal of an interaction (either the array or a single value) stored
in this data structure.  

Uses copy. If internal array is provided, the whole array is set. Otherwise
index and value must be provided and a single value is set.  

usage : inter_handler_2D_tsetInternal(inter_id, icdan, internal) or
inter_handler_2D_tsetInternal(inter_id, icdan, index, value)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  
internal(double array) : the new values  
index(integer) : the index where to set single value  
value(double ) : the new value to put at index  
";

%feature("docstring") inter_handler_2D_tgetInternal "

Get the internal of an interaction stored in this data structure.  

usage : internal = inter_handler_2D_tgetInternal(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  
internal(double array) : the new values array  
";

%feature("docstring") inter_handler_2D_getNbRecup "

return the number of recuped interactions of the selected type  

python usage : nb_recup = inter_handler_2D_getNbRecup(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction recuped of selected type  
";

%feature("docstring") inter_handler_2D_getNb "

return the number of interactions of the selected type stored in verlet data
structure  

python usage : nb_inter = inter_handler_2D_getNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_2D_getAllTactLawNb "

return the tact law number of all interactions stored in verlet data structure  

python usage : vector = inter_handler_2D_getAllTactLawNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
vector (int 1D-array) : mechanical data  
";

%feature("docstring") inter_handler_2D_getAll "

return coorx,coory,tx,ty,nx,ny,rlt,rln,vlt,vln,gaptt for all interactions stored
in verlet data structure  

python usage : array = inter_handler_2D_getAll(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_2D_getAllInternal "

return internal variables for all interactions stored in verlet data structure  

python usage : array = inter_handler_2D_getAllInternal(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_2D_getAllIdata "

return all integer data of all 'verlet' interaction  

Which are in order cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

python usage : array = inter_handler_2D_getAllIdata(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (int 2D-array) : identification data  
";

%feature("docstring") inter_handler_2D_getVerletAdjsz "

return integer number of verlet interaction of a candidate  

python usage : iantac = inter_handler_2D_getVerletAdjsz(inter_id, icdtac)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdtac(integer) : candidate contactor id  

Returns
-------
iantac (integer) : number of verlet interactions on candidate  
";

%feature("docstring") inter_handler_2D_getVerletIantac "

return integer antagonist contact of a verlet interaction  

python usage : iantac = inter_handler_2D_getVerletIantac(inter_id, icdtac, iadj)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdtac(integer) : candidate contactor id  
iadj(integer) : id of adjacent of candidate  

Returns
-------
iantac (integer) : id of antagonist contactor corresponding to verlet
interaction  
";

%feature("docstring") inter_handler_2D_computeRnod "

Put back the reac value of bodies from (this) interactions.  
";

%feature("docstring") inter_handler_2D_stockRloc "

stock from this to verlet  

python usage : inter_handler_2D_stockRloc(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_2D_recupRloc "

recup from verlet to this  

python usage : inter_handler_2D_recupRloc(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_2D_recupRlocByPos "

recup from verlet to this using position as criteria  

Only available for CLALp and PLPLx inter_id  

python usage : inter_handler_2D_recupRloc(inter_id, rtol)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
rtol(real) : tolerance to decide if contact is recup  
";


// File: wrap__CSASp_8h.xml

%feature("docstring") CSASp_SelectProxTactors "

contact detection between CSxxx and ASpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

If reset not equal to 0, the initialization flag is reset and detection skipped  

python usage : CSASp_SelectProxTactors(reset=0,use_external=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
use_external(integer) : if not 0, external detection is used  
";

%feature("docstring") CSASp_WriteLastVlocRloc "

write last local values of all CSASp contacts  

The values written are relative velocity, forces and local frame  

python usage : CSASp_WriteLastVlocRloc()  
";

%feature("docstring") CSASp_WriteOutVlocRloc "

write local values of all CSASp contacts  

The values written are relative velocity, forces and local frame  

python usage : CSASp_WriteOutVlocRloc()  
";

%feature("docstring") CSASp_DisplayOutVlocRloc "

display local values of all CSASp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CSASp_DisplayOutVlocRloc()  
";

%feature("docstring") CSASp_DisplayProxTactors "

display contacts  

python usage : CSASp_DisplayProxTactors()  
";

%feature("docstring") CSASp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CSASp_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CSASp_SkipAutoContact "

avoid CSxxx/ASpxx contact detection when they belong to the same entity  

python usage : CSASp_SkipAutoContact()  
";

%feature("docstring") CSASp_SetNonSymmetricDetection "

this function allows non symetric detection i.e. only one interaction is kept
when two bodies with candidate and antagonist contactors see each other  

python usage : CSASp_SetNonSymmetricDetection()  
";

%feature("docstring") CSASp_Trim "

trim contact (only contact within surface - not with extremities)  

python usage : CSASp_Trim()  
";

%feature("docstring") CSASp_SetTrimAngle "

set the trim angle (only contact within surface - not with extremities)  

python usage : CSASp_SetTrimAngle(angle)  

Parameters
----------
angle(real) : angle in degree - default 87 deg  
";

%feature("docstring") CSASp_AddReac "

add contact force to body Reac  

python usage : CSASp_AddReac()  
";

%feature("docstring") CSASp_AssumeOldFiles "

to read file with the CSpxx rank instead of CSxxx one  

python usage : CSASp_AssumeOldFiles()  
";

%feature("docstring") CSASp_CleanMemory "

Free all memory allocated within CSASp module.  

python usage : CSASp_CleanMemory()  
";


// File: wrap__CDCDx_8h.xml

%feature("docstring") CDCDx_SelectProxTactors "

contact detection between CYLND and CYLND tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list  

python usage : CDCDx_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") CDCDx_SmoothForceComputation "

computes smooth contact forces (if any)  

python usage : CDCDx_SmoothForceComputation()  
";

%feature("docstring") CDCDx_WriteLastVlocRloc "

write last local values of all CDCDx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDCDx_WriteLastVlocRloc()  
";

%feature("docstring") CDCDx_WriteOutVlocRloc "

write local values of all CDCDx contacts  

The values written are relative velocity, forces and local frame  

python usage : CDCDx_WriteOutVlocRloc()  
";

%feature("docstring") CDCDx_DisplayOutVlocRloc "

display local values of all CDCDx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CDCDx_DisplayOutVlocRloc()  
";

%feature("docstring") CDCDx_DisplayProxTactors "

display detected contacts  

python usage : CDCDx_DisplayProxTactors()  
";

%feature("docstring") CDCDx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read,
    -   num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : CDCDx_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CDCDx_SetXPeriodicCondition "

initialise data for simulation using periodic condition along X  

python usage : CDCDx_SetXPeriodicCondition(xperiod)  

Parameters
----------
xperiod(double) : period on x axis  
";

%feature("docstring") CDCDx_SetYPeriodicCondition "

initialise data for simulation using periodic condition along Y  

python usage : CDCDx_SetYPeriodicCondition(yperiod)  

Parameters
----------
yperiod(double) : period on y axis  
";

%feature("docstring") CDCDx_SetNumberInterByContact "

define the number of interaction by contact (experimental)  

python usage : CDCDx_SetNumberInterByContact(nb_interactions)  

Parameters
----------
nb_interactions(integer) : number of interactions per contact  
";

%feature("docstring") CDCDx_SetContactRadius "

define the contact radius (experimental)  

python usage : CDCDx_SetContactRadius(radius)  

Parameters
----------
radius(double) : contact radius  
";

%feature("docstring") CDCDx_CleanMemory "

Free all memory allocated within CDCDx module.  

python usage : CDCDx_CleanMemory()  
";


// File: wrap__parameters_8h.xml

%feature("docstring") parameters_getPhysicTypeId "

Get the id a body type from its name.  

usage i_param = parameters_getPhysicTypeId(bodyName)  

Parameters
----------
bodyName(string): body type name  

Returns
-------
i_param (int) : body type parameter  
";

%feature("docstring") parameters_getPhysicTypeNames "

Get the list of body types.  

usage bodyNames = parameters_getPhysicTypeName()  

Returns
-------
bodyName (string array) : body type names  
";

%feature("docstring") parameters_getBodyModelId "

Get the id a body model from its name.  

usage i_param = parameters_getBodyModelId(bodyName)  

Parameters
----------
bodyName(string): body model name  

Returns
-------
i_param (int) : body model parameter  
";

%feature("docstring") parameters_getBodyModelNames "

Get the list of body types.  

usage bodyNames = parameters_getBodyModelName()  

Returns
-------
bodyName (string array) : body model names  
";

%feature("docstring") parameters_getContactorId "

Get the id a contactor from its name.  

usage i_param = parameters_getContactorId(bodyName)  

Parameters
----------
bodyName(string): contactor name  

Returns
-------
i_param (int) : contactor parameter  
";

%feature("docstring") parameters_getContactorNames "

Get the list of contactors.  

usage bodyNames = parameters_getContactorName()  

Returns
-------
bodyName (string array) : contactor names  
";

%feature("docstring") parameters_getInteractionId "

Get the id a interaction from its name.  

usage i_param = parameters_getInteractionId(bodyName)  

Parameters
----------
bodyName(string): interaction name  

Returns
-------
i_param (int) : interaction parameter  
";

%feature("docstring") parameters_getInteractionNames "

Get the list of interactions.  

usage bodyNames = parameters_getInteractionName()  

Returns
-------
bodyName (string array) : interaction names  
";

%feature("docstring") parameters_getMatrixStorageId "

Get the id a matrix storage from its name.  

usage i_param = parameters_getMatrixStorageId(bodyName)  

Parameters
----------
bodyName(string): matrix storage name  

Returns
-------
i_param (int) : matrix storage parameter  
";

%feature("docstring") parameters_getMatrixStorageNames "

Get the list of matrix storages.  

usage bodyNames = parameters_getMatrixStorageName()  

Returns
-------
bodyName (string array) : matrix storage names  
";

%feature("docstring") parameters_getMatrixShapeId "

Get the id a matrix shape from its name.  

usage i_param = parameters_getMatrixShapeId(bodyName)  

Parameters
----------
bodyName(string): matrix shape name  

Returns
-------
i_param (int) : matrix shape parameter  
";

%feature("docstring") parameters_getMatrixShapeNames "

Get the list of matrix shapes.  

usage bodyNames = parameters_getMatrixShapeName()  

Returns
-------
bodyName (string array) : matrix shape names  
";

%feature("docstring") parameters_getGeneralizedCoordinatesId "

Get the id a generalized coordinates from its name.  

usage i_param = parameters_getGeneralizedCoordinatesId(bodyName)  

Parameters
----------
bodyName(string): generalized coordinates name  

Returns
-------
i_param (int) : generalized coordinates parameter  
";

%feature("docstring") parameters_getGeneralizedCoordinatesNames "

Get the list of generalized coordinatess.  

usage bodyNames = parameters_getGeneralizedCoordinatesName()  

Returns
-------
bodyName (string array) : generalized coordinates names  
";

%feature("docstring") parameters_getSurfaceEnergyStatusId "

Get the id a surface energy status from its name.  

usage i_param = parameters_getSurfaceEnergyStatusId(bodyName)  

Parameters
----------
bodyName(string): surface energy status name  

Returns
-------
i_param (int) : surface energy status parameter  
";

%feature("docstring") parameters_getSurfaceEnergyStatusNames "

Get the list of surface energy statuss.  

usage bodyNames = parameters_getSurfaceEnergyStatusName()  

Returns
-------
bodyName (string array) : surface energy status names  
";

%feature("docstring") parameters_getInterLawId "

Get the id a inter law from its name.  

usage i_param = parameters_getInterLawId(bodyName)  

Parameters
----------
bodyName(string): inter law name  

Returns
-------
i_param (int) : inter law parameter  
";

%feature("docstring") parameters_getInterLawNames "

Get the list of inter laws.  

usage bodyNames = parameters_getInterLawName()  

Returns
-------
bodyName (string array) : inter law names  
";

%feature("docstring") parameters_getIntegratorId "

Get the id a integrator from its name.  

usage i_param = parameters_getIntegratorId(bodyName)  

Parameters
----------
bodyName(string): integrator name  

Returns
-------
i_param (int) : integrator parameter  
";

%feature("docstring") parameters_getIntegratorNames "

Get the list of integrators.  

usage bodyNames = parameters_getIntegratorName()  

Returns
-------
bodyName (string array) : integrator names  
";

%feature("docstring") parameters_getNodeId "

Get the id a node from its name.  

usage i_param = parameters_getNodeId(bodyName)  

Parameters
----------
bodyName(string): node name  

Returns
-------
i_param (int) : node parameter  
";

%feature("docstring") parameters_getNodeNames "

Get the list of nodes.  

usage bodyNames = parameters_getNodeName()  

Returns
-------
bodyName (string array) : node names  
";

%feature("docstring") parameters_getDimeModeId "

Get the id a dime mode from its name.  

usage i_param = parameters_getDimeModeId(bodyName)  

Parameters
----------
bodyName(string): dime mode name  

Returns
-------
i_param (int) : dime mode parameter  
";

%feature("docstring") parameters_getDimeModeNames "

Get the list of dime modes.  

usage bodyNames = parameters_getDimeModeName()  

Returns
-------
bodyName (string array) : dime mode names  
";

%feature("docstring") parameters_getBodyVectorId "

Get the id a body vector from its name.  

usage i_param = parameters_getBodyVectorId(bodyName)  

Parameters
----------
bodyName(string): body vector name  

Returns
-------
i_param (int) : body vector parameter  
";

%feature("docstring") parameters_getBodyVectorNames "

Get the list of body vectors.  

usage bodyNames = parameters_getBodyVectorName()  

Returns
-------
bodyName (string array) : body vector names  
";

%feature("docstring") parameters_getContactStatusId "

Get the id a contact status from its name.  

usage i_param = parameters_getContactStatusId(bodyName)  

Parameters
----------
bodyName(string): contact status name  

Returns
-------
i_param (int) : contact status parameter  
";

%feature("docstring") parameters_getContactStatusNames "

Get the list of contact statuss.  

usage bodyNames = parameters_getContactStatusName()  

Returns
-------
bodyName (string array) : contact status names  
";

%feature("docstring") parameters_checkAll "

Check the consistency of all parameters id and name.  

usage parameters_checkAll()  
";


// File: wrap__CLxxx_8h.xml

%feature("docstring") CLxxx_LoadTactors "

load CLxxx from MAILx and Initialize existing_entities  

python usage : CLxxx_LoadTactors()  
";

%feature("docstring") CLxxx_SetNbNodesByCLxxx "

Set the number of CL nodes by edges. It helps to compute the length associated
to a contact node. Default is 2.  

python usage : CLxxx_SetNbNodesByCLxxx(nb_nodes)  

Parameters
----------
nb_nodes(integer) : number of CLxxx contactors by edges  
";

%feature("docstring") CLxxx_PushPreconNodes "

set CLxxx supporting nodes as precon  

python usage : CLxxx_PushPreconNodes()  
";

%feature("docstring") CLxxx_GetNbCLxxx "

Get the number of CLxxx.  

usage : nb_CLxxx = CLxxx_GetNbCLxxx()  

Parameters
----------
nb_CLxxx(integer) : number of CLxxx in container  
";

%feature("docstring") CLpxx_GetAllConnec "

return connectivity of all CL in a single vector using gloab node numbering of
mecaMAILx  

python usage : connec = CLxxx_getAllConnec()  

Returns
-------
connec (integer 1D-array) : connectiviy of CLxxx elements  
";

%feature("docstring") CLpxx_GetAllData "

return integer (ibdyty, itacty, i_as) and real data (normal) of all CLxxx  

python usage : idata, rdata = CLxxx_getAllData()  

Returns
-------
idata (integer 2D-array) : integer data array  

Returns
-------
rdata (real 2D-array) : real data array  
";

%feature("docstring") CLxxx_CleanMemory "

Free all memory allocated within CLxxx module.  

python usage : CLxxx_CleanMemory()  
";


// File: wrap__inter__handler__3D_8h.xml

%feature("docstring") inter_handler_3D_tgetNb "

return the number of interactions of the selected type stored in this data
structure  

python usage : nb_inter = inter_handler_3D_tgetNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_3D_tgetTactLawNb "

return the contact law number of an interaction stored in this data structure  

python usage : tact_law = inter_handler_3D_tgetTactLawNb(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
tact_law (integer) : contact law number  
";

%feature("docstring") inter_handler_3D_tgetIdBodies "

return the serial numbers of contacting objects of an interaction stored in this
data structure  

python usage : idBodies = inter_handler_3D_tgetIdBodies(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
idBodies (integer) : array with cd and an bodies serial number  
";

%feature("docstring") inter_handler_3D_tgetIData "

Get the integer data of an interaction stored in this data structure.  

idata vector holds cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

usage : idata = inter_handler_3D_tgetIData(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
idata (integer array) : the values array  
";

%feature("docstring") inter_handler_3D_tgetRData "

return the real data associated with an interactions  

Get an output array with, in this order, : coor, t/n/suc, rlt/n/s, vlt/n/s,
gapTT  

python usage : rdata = inter_handler_3D_tgetRData(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  

Returns
-------
rdata (double array) : array with real data of the interaction  
";

%feature("docstring") inter_handler_3D_tsetInternal "

Set the internal of an interaction (either the array or a single value) stored
in this data structure.  

Uses copy. If internal array is provided, the whole array is set. Otherwise
index and value must be provided and a single value is set.  

usage : inter_handler_3D_tsetInternal(inter_id, icdan, internal) or
inter_handler_3D_tsetInternal(inter_id, icdan, index, value)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  
internal(double array) : the new values array  
index(integer) : the index where to set single value  
value(double ) : the new value to put at index  
";

%feature("docstring") inter_handler_3D_tgetInternal "

Get the internal of an interaction stored in this data structure.  

usage : internal = inter_handler_3D_tgetInternal(inter_id, icdan)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdan(integer) : index of the interaction of selected type  
internal(double array) : the new values array  
";

%feature("docstring") inter_handler_3D_getNbRecup "

return the number of recup interactions of the selected type  

python usage : nb_recup = inter_handler_3D_getNbRecup(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_recup (integer) : number of interaction recup of selected type  
";

%feature("docstring") inter_handler_3D_getNb "

return the number of interactions of the selected type stored in verlet data
structure  

python usage : nb_inter = inter_handler_3D_getNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_3D_getAllTactLawNb "

return the tact law number of all interactions stored in verlet data structure  

python usage : vector = inter_handler_3D_getAllTactLawNb(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
vector (int 1D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAll "

return
coorx,coory,coorz,tx,ty,tz,nx,ny,nz,sx,sy,sz,rlt,rln,rls,vlt,vln,vls,gaptt of
all 'verlet' interactions  

python usage : array = inter_handler_3D_getAll(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAllInternal "

return contact point internal variables of all 'verlet' interactions  

python usage : array = inter_handler_3D_getAllInternal()  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAllIdata "

return all integer data of all 'verlet' interaction  

Which are in order cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

python usage : array = inter_handler_3D_getAllIdata(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (int 2D-array) : identification data  
";

%feature("docstring") inter_handler_3D_getVerletAdjsz "

return integer number of verlet interaction of a candidate  

python usage : iantac = inter_handler_3D_getVerletAdjsz(inter_id, icdtac)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdtac(integer) : candidate contactor id  

Returns
-------
iantac (integer) : number of verlet interactions on candidate  
";

%feature("docstring") inter_handler_3D_getVerletIantac "

return integer antagonist contact of a verlet interaction  

python usage : iantac = inter_handler_3D_getVerletIantac(inter_id, icdtac, iadj)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
icdtac(integer) : candidate contactor id  
iadj(integer) : id of adjacent of candidate  

Returns
-------
iantac (integer) : id of antagonist contactor corresponding to verlet
interaction  
";

%feature("docstring") inter_handler_3D_computeRnod "

Put back the Reac value of bodies from (this) interactions.  
";

%feature("docstring") inter_handler_3D_stockRloc "

stock from this to verlet  

python usage : inter_handler_3D_stockRloc(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_3D_recupRloc "

recup from verlet to this  

python usage : inter_handler_3D_recupRloc(inter_id)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_3D_recupRlocByPos "

recup from verlet to this using position as criteria  

Only available for CSASp inter_id  

python usage : inter_handler_3D_recupRloc(inter_id, rtol)  

Parameters
----------
inter_id(integer) : type of interaction (lmgc90 parameter)  
rtol(real) : tolerance to decide if contact is recup  
";


// File: wrap__POLYG_8h.xml

%feature("docstring") POLYG_LoadTactors "

load POLYG from RBDY2 and initialize existing_entites  

python usage : POLYG_LoadTactors()  
";

%feature("docstring") POLYG_GetMinRadius "

give min radius used during detection  

python usage : POLYG_GetMinRadius()  
";

%feature("docstring") POLYG_GetMaxRadius "

give max radius used during detection  

python usage : POLYG_GetMaxRadius()  
";

%feature("docstring") POLYG_GetNbPOLYG "

Get the number of POLYG in container.  

python usage : nb_polyg = POLYG_GetNbPOLYG()  

Returns
-------
nb_polyg (integer) : the number of POLYG in container  
";

%feature("docstring") POLYG_GetPOLYG2BDYTY "

Get a copy of map POLYG2bdyty.  

usage : polyr2bdyty = POLYG_GetPOLYG2BDYTY()  

Returns
-------
polyr2bdyty (integer 2D-array) : the polyr2bdyty map  
";

%feature("docstring") POLYG_GetPtrPOLYG2BDYTY "

return a pointer onto the map polyg2rbdy2  

python usage : polyg2rbdy2 = POLYG_GetPtrPOLYG2BDYTY()  

Returns
-------
polyg2rbdy2 (integer array) : reference on map between polyg rank and
body/tactor rank  
";

%feature("docstring") POLYG_IsVisible "

return if a body visible  

usage : visible = POLYG_IsVisible(itact)  

Parameters
----------
itact(integer) : rank of POLYG  
visible(integer) : 1 if body is visible, 0 else  
";

%feature("docstring") POLYG_GetContactorRadius "

Get the radius of a given POLYG.  

python usage : radius = POLYG_GetContactorRadius(itact)  

Parameters
----------
itact(integer) : rank of a POLYG  

Returns
-------
radius (double) : the radius of the POLYG of rank itact  
";

%feature("docstring") POLYG_GetNbVertices "

Get the number of vertices of the first POLYG of a body.  

python usage : nb_vertices = POLYG_GetNbVertices(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of a body  

Returns
-------
nb_vertices (integer) : the number of vertices of the first POLYG of the body  
";

%feature("docstring") POLYG_GetVertices "

Get the coordinates of the vertices of the first POLYG of a body.  

usage : vertices = POLYG_GetVertices(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  
vertices(double 2D-array) : the coordinates of the vertices  
";

%feature("docstring") POLYG_GetNbVertex "

Get the number of vertices of a POLYG.  

usage : nb_vertex = POLYG_GetNpVertex(itacty)  

Parameters
----------
itacty(integer) : id of the POLYG contactor  

Returns
-------
nb_vertex (int) : the number of vertices of the POLYG  
";

%feature("docstring") POLYG_GetVertex "

Get the outline of a POLYG.  

usage : vertex = POLYG_GetVertex(itacty, length)  

Parameters
----------
ibdyty(integer) : rank of considered body  
length(integer) : 2 * number of vertices  
vertex(double array) : the coordinates of the vertices  
";

%feature("docstring") POLYG_GetBodyId "

Get the id of the body which the tactor belongs.  

python usage : id = POLYG_GetBodyId(itacty)  

Parameters
----------
itacty(integer) : rank of a POLYG contactor  

Returns
-------
id (integer) : the id of the body  
";

%feature("docstring") POLYG_InitOutlines "

Get a reference on the outlines of all POLYG.  

usage : outlines = POLYG_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_POLYG  
";

%feature("docstring") POLYG_InitScalarFields "

Get a reference on the scalar fields of all POLYG.  

usage : scalarfields = POLYG_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_POLYG array  
";

%feature("docstring") POLYG_UpdatePostdata "

Update values of outlines_POLYG and scalarfields_POLYG pointers.  

usage : POLYG_UpdatePostdata  
";

%feature("docstring") POLYG_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = POLYG_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
POLYG  
";

%feature("docstring") POLYG_GetNbScalarFields "

Get the number of scalar fields of a POLYG.  

python usage : nb_scalarfields = POLYG_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a POLYG  
";

%feature("docstring") POLYG_SetXdilation "
";

%feature("docstring") POLYG_SetVdilation "
";

%feature("docstring") POLYG_CleanMemory "

Free all memory allocated within POLYG module.  

python usage : POLYG_CleanMemory()  
";


// File: wrap__CLALp_8h.xml

%feature("docstring") CLALp_SelectProxTactors "

contact detection between CLxxx and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
If reset not equal to 0, the initialization flag is reset and detection skipped  

python usage : CLALp_SelectProxTactors(reset=0, use_external=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
use_external(integer) : if not 0, external detection is used  
";

%feature("docstring") CLALp_UpdateWear "

python usage : CLALp_UpdateWear()  
";

%feature("docstring") CLALp_WriteLastVlocRloc "

write last local values of all CLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : CLALp_WriteLastVlocRloc()  
";

%feature("docstring") CLALp_WriteOutVlocRloc "

write local values of all CLALp contacts  

The values written are relative velocity, forces and local frame  

python usage : CLALp_WriteOutVlocRloc()  
";

%feature("docstring") CLALp_DisplayOutVlocRloc "

display local values of all CLALp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CLALp_DisplayOutVlocRloc()  
";

%feature("docstring") CLALp_DisplayProxTactors "

display contacts  

python usage : CLALp_DisplayProxTactors()  
";

%feature("docstring") CLALp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

python usage : CLALp_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") CLALp_SetNonSymmetricDetection "

this function allows non symmetric detection i.e. only one interaction is kept
when two bodies with candidate and antagonist contactors see each other  

python usage : CLALp_SetNonSymmetricDetection()  
";

%feature("docstring") CLALp_Trim "

trim contact (only contact within a line - not with extremities)  

python usage : CLALp_Trim()  
";

%feature("docstring") CLALp_CleanMemory "

Free all memory allocated within CLALp module.  

python usage : CLALp_CleanMemory()  
";


// File: wrap__PLALp_8h.xml

%feature("docstring") PLALp_SelectProxTactors "

contact detection between POLYG and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PLALp_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
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
num(integer) : which VlocRloc file to read  
";

%feature("docstring") PLALp_CleanMemory "

Free all memory allocated within PLALp module.  

python usage : PLALp_CleanMemory()  
";


// File: wrap__PLANx_8h.xml

%feature("docstring") PLANx_LoadTactors "

load PLANx from RBDY3 and initialize existing_entites  

python usage : PLANx_LoadTactors()  
";

%feature("docstring") PLANx_GetNbPLANx "

Get the number of PLANx.  

python usage : nb_PLANx = PLANx_GetNbPLANx()  

Returns
-------
nb_PLANx (integer) : the number of PLANx  
";

%feature("docstring") PLANx_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = PLANx_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") PLANx_GetPtrPLANx2BDYTY "

return a pointer onto the map planx2bdyty  

python usage : planx2bdyty = PLANx_GetPtrPLANx2BDYTY()  

Returns
-------
planx2bdyty (integer array) : reference on map between planx rank and body rank  
";

%feature("docstring") PLANx_InitOutlines "

Get a reference on the outlines of all PLANx.  

usage : outlines = PLANx_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_PLANx  
";

%feature("docstring") PLANx_InitScalarFields "

Get a reference on the scalar fields of all PLANx.  

usage : scalarfields = PLANx_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_PLANx array  
";

%feature("docstring") PLANx_UpdatePostdata "

Update values of outlines_PLANx and scalarfields_PLANx pointers.  

usage : PLANx_UpdatePostdata  
";

%feature("docstring") PLANx_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = PLANx_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
PLANx  
";

%feature("docstring") PLANx_GetNbScalarFields "

Get the number of scalar fields of a PLANx.  

python usage : nb_scalarfields = PLANx_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a PLANx  
";

%feature("docstring") PLANx_GetPtrAllConnectivities "

Get a reference on the connectivities of all PLANx.  

usage : connec = PLANx_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") PLANx_CleanMemory "

Free all memory allocated within PLANx module.  

python usage : PLANx_CleanMemory()  
";


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
ivalue(integer) : id of the therMAILx  

Returns
-------
nb_nodes (integer) : number of nodes of a therMAILx  
";

%feature("docstring") therMAILx_GetNbElements "

Get the number of nodes of a therMAILx.  

python usage : nb_nodes = therMAILx_GetNbElements(ibdyty)  

Parameters
----------
ivalue(integer) : id of the therMAILx  

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
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeCapacity "

computes the elemetary capacity matrices  

python usage : therMAILx_ComputeCapacity(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeConvection "

compute elementary convection terms  

python usage : therMAILx_ComputeConvection(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeInternalFlux "

compute elementary internal flux  

python usage : therMAILx_ComputeInternalFlux(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeExternalFlux "

compute elementary external flux  

python usage : therMAILx_ComputeExternalFlux(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_AssembThermKT "

assembles elementary matrices  

python usage : therMAILx_AssembKT(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_AssembThermRHS "

assembles elementary vectors  

python usage : therMAILx_AssembRHS(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeThermDof "

computes current dof  

python usage : therMAILx_ComputeThermDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeThermFields "

computes elementary fields  

python usage : therMAILx_ComputeThermFields(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_UpdateThermDof "

update begin dof with current dof  

python usage : therMAILx_UpdateThermDof(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_UpdateThermBulk "

update begin elementary fields with current elementary fields  

python usage : therMAILx_UpdateThermBulk(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
    works on all objects  
";

%feature("docstring") therMAILx_ComputeResidueNorm "

compute the residue of the thermal equation  

python usage : norm = therMAILx_ComputeResidueNorm(i_list)  

Parameters
----------
i_list(list of integer) : list of bodies to compute conductivities if omitted
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
num(integer) : which DOF file to read  
";

%feature("docstring") therMAILx_ReadIniGPV "

Read GPV file.  

If num <= 0 : DATBOX/GPV.INI file is read  

Else : OUTBOX/GPV.OUT.num is read, num being the parameter used in
TimeEvolution_ReadIniGPV last call  

python usage : therMAILx_ReadIniGPV(num=0)  

Parameters
----------
num(integer) : which GPV file to read  
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
datatype(string of size 5) : the vector to set  
ibdyty(integer) : rank of body  
matrix(double array) : the new values  
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
datatype(string of size 5) : the vector to get  
ibdyty(integer) : rank of considered body  

Returns
-------
vector (double 2D-array) : the desired vector  
";

%feature("docstring") therMAILx_GetScalarFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = therMAILx_GetScalarFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") therMAILx_SetScalarFieldByElement "

Update elementary scalar field through a element external field on a given body.  

Field values are stored at Gauss point, on an element all Gauss point have the
element value  

You need to declare this field in your MODELS.DAT  

python usage : therMAILx_SetScalarFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the field to set  
f(double array) : value of the field  
";

%feature("docstring") therMAILx_GetVectorFieldRank "

Get the rank of field of an element of a body from its name.  

python usage : f_rank = therMAILx_GetVectorFieldRank(ibdyty, iblmty, name)  

Parameters
----------
ibdyty(integer) : id of the concern body  
iblmty(integer) : id of the concern element  
name(string) : name of the desired vector field  

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
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") therMAILx_SetVectorFieldByElement "

Update elementary fields through a nodal external field on a given body.  

Use the form functions of the elements and input values to compute and store
field values at Gauss points.  

You need to declare this field in your MODELS.DAT  

python usage : therMAILx_SetFieldByElement(IdBody, f_rank, f)  

Parameters
----------
IdBody(integer) : id of the concern body  
f_rank(integer) : rank of the vector field to set  
f(double array) : value of the vector field  
";

%feature("docstring") therMAILx_AddSource "

Add a volumic source into a given body.  

python usage : therMAILx_AddSource(ibdyty, ifield)  

Parameters
----------
ibdyty(integer) : rank of body  
ifield(integer) : rank of field  
";

%feature("docstring") therMAILx_AddNodalFieldDivergence "

Add the divergence of a field to external flux.  

python usage : therMAILx_AddNodalFieldDivergence(ibdyty, ifield)  

Parameters
----------
ibdyty(integer) : rank of body  
ifield(integer) : rank of field  
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
ibdyty(integer) : rank of considered body  

Returns
-------
grad_T (double 2D-array) : the desired gradient  
";

%feature("docstring") therMAILx_GetFlux "

Get a copy of a gradient of a given body.  

Python usage : Flux_T = therMAILx_GetFlux(ibdyty)  

Parameters
----------
ibdyty(integer) : rank of considered body  

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
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates  
";

%feature("docstring") therMAILx_GetConnectivity "

return connectivity of idBody elements  

python usage : vector = therMAILx_GetConnectivity(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
vector (integer) : connectivity  
";

%feature("docstring") therMAILx_GetAll "

return mechanical data computed for idBody  

python usage : array = therMAILx_GetAll(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") therMAILx_GetGpCoor "

return Gauss points coordinates of idBody  

python usage : array = therMAILx_GetGpCoor(idBody)  

Parameters
----------
IdBody(integer) : id of the concerned body  

Returns
-------
array (double 2D-array) : coordinates of all Gauss points  
";

%feature("docstring") therMAILx_GetGpField "

return field values stored at a gp  

python usage : field = therMAILx_GetGpField(idBody,idEle,idGp,idField)  

Parameters
----------
IdBody(integer) : id of the concerned body  
IdEle(integer) : id of the concerned element  
IdGp(integer) : id of the concerned gauss point  
IdField(integer) : id of the concerned field  

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
ibdyty(integer) : id of the therMAILx  
iblmty(integer) : id of the element  

Returns
-------
nb_gp (integer) : number of Gauss point of an element of a therMAILx  
";


// File: wrap__CYLND_8h.xml

%feature("docstring") CYLND_LoadTactors "

load CYLND from RBDY3 and initialize existing_entites  

python usage : CYLND_LoadTactors()  
";

%feature("docstring") CYLND_IsVisible "

return if a given contactor is attached to a visible body  

python usage : visible = CYLND_IsVisible(itacty)  

Parameters
----------
itacty(integer) : id of the contactor we want visibility  

Returns
-------
visible (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") CYLND_GetNbCYLND "

Get the number of CYLND.  

python usage : nb_CYLND = CYLND_GetNbCYLND()  

Returns
-------
nb_CYLND (integer) : the number of CYLND  
";

%feature("docstring") CYLND_GetShape "

Get the shape of a CYLND.  

usage : shape = CYLND_GetShape(itacty)  

Parameters
----------
itacty(integer) : rank of CYLND  

Returns
-------
shape (double array) : axis length of the CYLND  
";

%feature("docstring") CYLND_GetPtrCYLND2BDYTY "

return a pointer onto the map cylnd2bdyty  

python usage : cylnd2bdyty = CYLND_GetPtrCYLND2BDYTY()  

Returns
-------
cylnd2bdyty (integer array) : reference on map between cylnd rank and body rank  
";

%feature("docstring") CYLND_InitOutlines "

Get a reference on the outlines of all CYLND.  

python usage : outlines = CYLND_InitOutlines()  

Returns
-------
outlines (double array) : a reference on outlines_CYLND  
";

%feature("docstring") CYLND_InitScalarFields "

Get a reference on the scalar fields of all CYLND.  

python usage : scalarfields = CYLND_InitScalarfields()  

Returns
-------
scalarfields (double array) : reference on scalarfields_CYLND array  
";

%feature("docstring") CYLND_UpdatePostdata "

Update values of outlines_CYLND and scalarfields_CYLND pointers.  

python usage : CYLND_UpdatePostdata  
";

%feature("docstring") CYLND_GetNbPointOutlines "

Get the list of cumulated outline points number.  

python usage : nb_pointOutlines = CYLND_GetNbPointOutlines()  

Returns
-------
nb_pointOutlines (integer array) : the cumulated number of outline points of the
CYLND  
";

%feature("docstring") CYLND_GetNbScalarFields "

Get the number of scalar fields of a CYLND.  

python usage : nb_scalarfields = CYLND_GetNbScalarFields()  

Returns
-------
nb_scalarfields (integer) : the number of scalar fields of a CYLND  
";

%feature("docstring") CYLND_GetPtrAllConnectivities "

Get a reference on the connectivities of all CYLND.  

python usage : connec = CYLND_GetPtrAllConnectivities()  

Returns
-------
connec (integer array) : a reference on all_connectivities  
";

%feature("docstring") CYLND_CleanMemory "

Free all memory allocated within CYLND module.  

python usage : CYLND_CleanMemory()  
";


// File: wrap__DKALp_8h.xml

%feature("docstring") DKALp_SelectProxTactors "

contact detection between DISKx and ALpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : DKALp_SelectProxTactors(reset=0)  

Parameters
----------
reset(integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") DKALp_WriteLastVlocRloc "

write last local values of all DKALp contacts  

The values written are relative velocity, forces and local frame  

python usage : DKALp_WriteLastVlocRloc()  
";

%feature("docstring") DKALp_WriteOutVlocRloc "

write local values of all DKALp contacts  

The values written are relative velocity, forces and local frame  

python usage : DKALp_WriteOutVlocRloc()  
";

%feature("docstring") DKALp_DisplayOutVlocRloc "

display local values of all DKALp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : DKALp_DisplayOutVlocRloc()  
";

%feature("docstring") DKALp_DisplayProxTactors "

display contacts  

python usage : DKALp_DisplayProxTactors()  
";

%feature("docstring") DKALp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

python usage : DKALp_ReadIniVlocRloc(num=0)  

Parameters
----------
num(integer) : which VlocRloc file to read  
";

%feature("docstring") DKALp_CleanMemory "

Free all memory allocated within DKALp module.  

python usage : DKALp_CleanMemory()  
";

