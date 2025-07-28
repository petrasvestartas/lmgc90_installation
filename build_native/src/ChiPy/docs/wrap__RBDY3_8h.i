
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
* `ibdyty` :  
    (integer) : rank of considered  
* `idrvdof` :  
    (integer) : index of velocity driven dof to set  
* `value` :  
    (real) : new value of the velocity driven dof  
";

%feature("docstring") RBDY3_FatalDamping "

Nullify body velocities (current and initial) of a list of bodies.  

python usage : RBDY3_FatalDamping(i_list)  

Parameters
----------
* `i_list` :  
    (list of integer) : list of bodies to reset current velocity if omitted
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
* `ifrom` :  
    (integer) : begining of bodys' index that will be written  
* `ito` :  
    (integer) : end of bodys'index that will be written  
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
* `num` :  
    (integer) : which DOF file to read  
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
* `Zmin` :  
    (real) : inferior boundary value  
";

%feature("docstring") RBDY3_SetZmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetZmaxBoundary(Zmax)  

Parameters
----------
* `Zmax` :  
    (real) : superior boundary value  
";

%feature("docstring") RBDY3_SetYminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetYminBoundary(Ymin)  

Parameters
----------
* `Ymin` :  
    (real) : left boundary value  
";

%feature("docstring") RBDY3_SetYmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetYmaxBoundary(Ymax)  

Parameters
----------
* `Ymax` :  
    (real) : right boundary value  
";

%feature("docstring") RBDY3_SetXminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetXminBoundary(Xmin)  

Parameters
----------
* `Xmin` :  
    (real) : inferior boundary value  
";

%feature("docstring") RBDY3_SetXmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

python usage : RBDY3_SetXmaxBoundary(Xmax)  

Parameters
----------
* `Xmax` :  
    (real) : front boundary value  
";

%feature("docstring") RBDY3_SetXPeriodicCondition "

set the period on X axis  

python usage : RBDY3_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (real) : period on x axis  
";

%feature("docstring") RBDY3_SetYPeriodicCondition "

set the periode on Y axis  

python usage : RBDY3_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (real) : period on y axis  
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
* `ibdyty(integer)` :  
    : index of the RBDY3  
";

%feature("docstring") RBDY3_SetInvisible "

rended a given RBDY3 invisible  

python usage : RBDY3_SetInvisible(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the RBDY3  
";

%feature("docstring") RBDY3_IsVisible "

return if a given body visible  

python usage : visible = RBDY3_IsVisible(ibdyty)  

Parameters
----------
* `idbdy(integer)` :  
    : id of the body we want visibility  

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
* `ibdyty(integer)` :  
    : rank of the RBDY3  

Returns
-------
density(double) : density of the RBDY3  
";

%feature("docstring") RBDY3_GetBodyInertia "

Get the principal inertia of a given RBDY3.  

python usage : inertia = RBDY3_GetBodyInertia(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3  

Returns
-------
inertia (double array) : inertia vector of the desired RBDY3  
";

%feature("docstring") RBDY3_GetAllInertia "

Get the inertia of a all RBDY3 body.  

usage : inertia = RBDY3_GetAllInertia()  

Parameters
----------
* `inertia` :  
    (double array): the inertia of all bodies  
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
* `datatype` :  
    (string [5]) : the vector to set  
* `ibdyty` :  
    (integer) : rank of the RBDY3  
* `vector` :  
    (double array) : the new value of the vector  
";

%feature("docstring") RBDY3_PutAllBodyVector "

Put an array of a vector of all RBDY3 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY3_PutBodyVector  

python usage : RBDY3_PutAllBodyVector(datatype, matrix)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to set  
* `matrix` :  
    (double array) : input matrix  
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
* `datatype` :  
    (string [5]) : the vector to get  
* `ibdyty` :  
    (integer) : rank of the RBDY3  

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
* `datatype` :  
    (string [5]) : the vector to get  

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
* `datatype` :  
    (string [5]) : the vector to set  
* `ibdyty` :  
    (integer) : rank of the RBDY3  

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
* `datatype` :  
    (string [5]) : the vector to set  
* `ibdyty` :  
    (integer) : rank of the RBDY3  
* `matrix` :  
    (double array) : a matrix  
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
* `datatype` :  
    (string [5]) : the vector to get  
* `ibdyty` :  
    (integer) : rank of the RBDY3  

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
* `ibdyty` :  
    (integer) : rank of the RBDY3  

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
* `ibdyty(int)` :  
    : index of the RBDY3  
* `mass(double**)` :  
    : mass matrix of the RBDY3  
";

%feature("docstring") RBDY3_GetVelocity "

Get the velocity of a body.  

Parameters
----------
* `ibdyty(int)` :  
    : index of the RBDY3  
* `velocity(double[6])` :  
    : velocity of the RBDY3  
";

%feature("docstring") RBDY3_GetGlobInertia "

Get the global inertia.  

usage : inertia = RBDY3_GetGlobInertia(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : id of desired RBDY3  

Returns
-------
inertia (double 2D array) : the inertia matrix  
";

%feature("docstring") RBDY3_GetBehavior "

Get the type of the nickname of the behavior.  

usage name = RBDY3_GetBehavior(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3 in container  

Returns
-------
type (string) : nickname  
";

%feature("docstring") RBDY3_GetNbContactor "

get the number of contactor of RBDY3  

python usage : nb = RBDY3_GetNbContactor(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3 in container  

Returns
-------
nb (integer) : number of contactor attached to a RBDY3  
";

%feature("docstring") RBDY3_GetContactorType "

Get the type of the itacty contactor of a body ibdyty.  

usage type = RBDY3_GetContactorType(ibdyty,itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3 in container  
* `itacty` :  
    (integer) : rank of the contactor in the RBDY3  

Returns
-------
type (string) : type of the contactor of the body  
";

%feature("docstring") RBDY3_SetContactorColor "

Set the color of a given contactor of a body.  

usage : RBDY3_SetContactorColor(ibdyty, itacty, color)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3  
* `itacty` :  
    (integer) : rank of the contactor in the RBDY3  
* `color` :  
    (string of size 5) : the color  
";

%feature("docstring") RBDY3_GetContactorColor "

Get the color of the itacty contactor of a body ibdyty.  

usage color = RBDY3_GetContactorColor(ibdyty,itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY3 in container  
* `itacty` :  
    (integer) : rank of the contactor in the RBDY3  

Returns
-------
color (string) : color of the contactor of the body  
";

%feature("docstring") RBDY3_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = RBDY3_getDrvVlocy(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY3  
* `drvdof_indices` :  
    (integer array) : indices list of driven dof  
* `drvdof_values` :  
    (real array) : values of the driven dof  
";

%feature("docstring") RBDY3_computeDrvVlocy "

Compute the value of the driven velocity of a body at current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : RBDY3_computeDrvVlocy(ibdyty, values)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY3  
* `values` :  
    (double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") RBDY3_WriteOutOneBody "

write a bdyty to BODIES.OUT with a given rank  

python usage : RBDY3_WriteOutOneBody(ibdyty, new_ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY3  
* `new_ibdyty` :  
    (integer): new index of the RBDY3  
";

%feature("docstring") RBDY3_WriteOutDofOneBody "

write a bdyty dof to DOF.OUT with a given rank  

python usage : RBDY3_WriteOutDofOneBody(ibdyty, new_ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY3  
* `new_ibdyty` :  
    (integer): new index of the RBDY3  
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
* `radius` :  
    (double) : radius threshold  
";

%feature("docstring") RBDY3_SetVisibleVlocyDrivenDof "

rended a given Velocy DOF visible  

python usage : RBDY3_SetVisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the RBDY3  
* `iccdof(integer)` :  
    : index of the DOF to set visible  
";

%feature("docstring") RBDY3_SetInvisibleVlocyDrivenDof "

rended a given Velocy DOF invisible  

python usage : RBDY3_SetInvisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the RBDY3  
* `iccdof(integer)` :  
    : index of the DOF to set invisible  
";

%feature("docstring") RBDY3_PartialDamping "

Limit body velocity to Vmax value.  

usage : RBDY3_PartialDamping(nb_steps, Vmax)  

Parameters
----------
* `nb_steps` :  
    (integer) : periodicity @parma[in] Vmax (double) : Vmax  
";

%feature("docstring") RBDY3_GetVolume "

Get volume of a body.  

usage : volume = RBDY3_GetVolume(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : RBDY3 id  

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
* `ibdyty` :  
    (integer) : rank of a RBDY3  

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
* `disper(double)` :  
    : some dispersion coefficient  
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
* `ibdyty` :  
    (integer) : rank of body  
* `itacty` :  
    (integer) : rank of tacty \"  
";

%feature("docstring") RBDY3_SetEquilibriumNorm "

Initialization of data for the equilibrium state check.  

You must precise the type of check test :  

*   Qvlcy : quadratic norm velocy  
*   Mvlcy : maximum norm velocy  

usage : RBDY3_CheckEquilibrium(norm_type , tolerance)  

Parameters
----------
* `norm_type` :  
    (string of size 5) : norm type use for the equilibrium check  
* `tolerance` :  
    (double) : norm tolerance  
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
* `first_RBDY3(int)` :  
    : number of first invisible body  
* `radius` :  
    : source point area radius  
* `Xshift` :  
    : X translation of deposited object from reference coordinate  
* `Yshift` :  
    : Y translation of deposited object from reference coordinate  
* `Zshift` :  
    : Z translation of deposited object from reference coordinate  
";

%feature("docstring") RBDY3_SetSourcePointWithIni "

create an assembly by source point deposit  

python usage : RBDY3_SetSourcePointWithIni(first_RBDY3, radius, Xshift, Yshift,
Zshift)  

Parameters
----------
* `first_RBDY3(int)` :  
    : number of first invisible body  
* `radius` :  
    : source point area radius  
* `Xshift` :  
    : X coordinate of deposited object  
* `Yshift` :  
    : Y coordinate of deposited object  
* `Zshift` :  
    : Z coordinate of deposited object  
";

%feature("docstring") RBDY3_InitializeProgressiveActivation "

set the progression of altitude  

python usage : RBDY3_InitializeProgressiveActivation(zini, dz)  

Parameters
----------
* `zini` :  
    (real) : initial altitude  
* `dz` :  
    (real) : increment of altitude  
";

%feature("docstring") RBDY3_ApplyProgressiveActivation "

set occurence of activation  

python usage : RBDY3_ApplyProgressiveActivation(freq)  

Parameters
----------
* `freq` :  
    (integer) : activation frequence of progression  
";

%feature("docstring") RBDY3_InitFreeBoundary "

python usage : RBDY3_InitFreeBoundary(xmin, xmax, ymin, ymax, radius)  

Parameters
----------
* `xmin` :  
    (real) :  
* `xmax` :  
    (real) :  
* `ymin` :  
    (real) :  
* `ymax` :  
    (real) :  
* `radius` :  
    (real) :  
";

%feature("docstring") RBDY3_TriaxialLoading "

Triaxial load of a sample using a rigid box.  

python usage : TriaxialLoading(num_down, num_right, num_up, num_left, num_front,
num_rear, nb_loads, loads)  

Parameters
----------
* `num_down` :  
    (integer) :  
* `num_right` :  
    (integer) :  
* `num_up` :  
    (integer) :  
* `num_left` :  
    (integer) :  
* `num_front` :  
    (integer) :  
* `num_rear` :  
    (integer) :  
* `nb_loads` :  
    (integer) : the number of walls you want to load with a pressure (1 to 6)  
* `loads` :  
    (array) : loads(2,nb_loads): load(1,i) contains which wall is loaded
    (1==down, 2==right, 3==up, 4==left, 5==front, 6==rear) and load(2,i)
    contains the amplitude of the stress (a positive value means compression).  
";

%feature("docstring") RBDY3_GetDofStatus "

Get dof status.  

python usage : status = RBDY3_GetDofStatus(ibdyty)  

Parameters
----------
* `ibdyty(integer)` :  
    : rank of the RBDY3  

Returns
-------
status(integer) : dof status of the RBDY3  
";

