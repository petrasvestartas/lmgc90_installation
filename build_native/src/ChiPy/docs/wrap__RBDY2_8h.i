
// File: wrap__RBDY2_8h.xml

%feature("docstring") RBDY2_PutBodyInvMass "

Set inv mass diagonal matrix of a given body. Overwrites the computed values.  

usage : RBDY2_PutBodyInvMass(ibdyty, inv_mass)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of RBDY2  
* `inv_mass` :  
    (double array) : inv_mass of RBDY2 (size 3)  
";

%feature("docstring") RBDY2_PutBodyPreconW "

Put preconW of a given body.  

usage : RBDY2_PutBodyPreconW(ibdyty, idof, W)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of RBDY2  
* `idof` :  
    (integer) : corresponding dof to set  
* `W` :  
    (double array) :  
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
* `datatype` :  
    (string of size 5) : the vector to set  
* `ibdyty` :  
    (integer) : rank of body  
* `vector` :  
    (double array) : the new value  
";

%feature("docstring") RBDY2_PutAllBodyVector "

Put an array of a vector of all RBDY2 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY2_PutBodyVector  

python usage : RBDY2_PutAllBodyVector(datatype, matrix)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to set  
* `matrix` :  
    (double array) : input matrix  
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
* `datatype` :  
    (string of size 5) : the vector to get  
* `ibdyty` :  
    (integer) : rank of considered body  
* `vector` :  
    (double array) : the desired vector  
";

%feature("docstring") RBDY2_GetAllBodyVector "

Get an array of a vector of all RBDY2 bodies (visible and invisible)  

Possible values for datatype field are: ... see RBDY2_GetBodyVector  

python usage : matrix = RBDY2_GetBodyVector(datatype, ibdyty)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to get  

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
* `datatype` :  
    (string of size 5) : the vector to get  
* `ibdyty` :  
    (integer) : rank of considered body  
* `vector_ptr` :  
    (double array) : reference on the desired vector viewed as a numpy array  
";

%feature("docstring") RBDY2_GetBodyInertia "

Get the inertia of a given RBDY2 body.  

usage : inertia = RBDY2_GetBodyInertia(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of considered body  
* `inertia` :  
    (double) : the inertia of desired body  
";

%feature("docstring") RBDY2_GetAllInertia "

Get the inertia of a all RBDY2 body.  

usage : inertia = RBDY2_GetAllInertia()  

Parameters
----------
* `inertia` :  
    (double array): the inertia of all bodies  
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
* `ibdyty` :  
    (integer) : rank of considered  
* `idrvdof` :  
    (integer) : index of velocity driven dof to set  
* `value` :  
    (real) : new value of the velocity driven dof  
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
* `i_list` :  
    (list of integer) : list of bodies to reset current velocity if omitted
    works on all objetcs  
";

%feature("docstring") RBDY2_PartialDamping "

Limit body velocity to Vmax value.  

usage : RBDY2_PartialDamping(nb_steps, Vmax)  

Parameters
----------
* `nb_steps` :  
    (integer) : periodicity @parma[in] Vmax (double) : Vmax  
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
* `ivalue1` :  
    (integer) : first body  
* `ivalue2` :  
    (integer) : last body  
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
* `num` :  
    (integer) : which DOF file to read  
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
* `disper` :  
    (double) : dispersion variable  
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
* `period` :  
    (double) : periode  
";

%feature("docstring") RBDY2_ResizeBodies "

resize body radius by a factor  

usage : RBDY2_ResizeBodies(homo)  

Parameters
----------
* `homo` :  
    (double) : resize factor  
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
* `ibdyty` :  
    (integer): rank of first invisible body  
* `radius` :  
    (double) : radius of source point area  
* `x_coor` :  
    (double) : X translation from the set of grains  
* `y_coor` :  
    (double) : Y translation from the set of grains  
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
* `down` :  
    (integer) : rank of the lower body  
* `up` :  
    (integer) : rank of the upper body  
* `thickness` :  
    (double) : thickness of the membrane  
* `stress` :  
    (double) : pressure on the membrane  
";

%feature("docstring") RBDY2_BiaxialLoading "

Biaxial load of a sample using a rigid box.  

usage : RBDY2_BiaxialLoading(down,f_down,right,f_right,up,f_up,left,f_left)  

Parameters
----------
* `down` :  
    (integer) : rank of the lower body  
* `f_down` :  
    (double) : pressure on the lower body  
* `right` :  
    (integer) : rank of the right body  
* `f_right` :  
    (double) : pressure on the right body  
* `up` :  
    (integer) : rank of the upper body  
* `f_up` :  
    (double) : pressure on the upper body  
* `left` :  
    (integer) : rank of the left body  
* `f_left` :  
    (double) : pressure on the left body  
";

%feature("docstring") RBDY2_SetYminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetYminBoundary(inf_boundary)  

Parameters
----------
* `inf_boundary` :  
    (double) : inferior boundary value  
";

%feature("docstring") RBDY2_SetYmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetYmaxBoundary(up_boundary)  

Parameters
----------
* `up_boundary` :  
    (double) : superior boundary value  
";

%feature("docstring") RBDY2_SetXminBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetXminBoundary(left_boundary)  

Parameters
----------
* `left_boundary` :  
    (double) : left boundary value  
";

%feature("docstring") RBDY2_SetXmaxBoundary "

define the boundary of command CHECK_OUT_OF_BOUNDS  

usage : RBDY2_SetXmaxBoundary(right_boundary)  

Parameters
----------
* `right_boundary` :  
    (double) : right boundary value  
";

%feature("docstring") RBDY2_SetEquilibriumNorm "

Initialization of data for the equilibrium state check.  

You must precise the type of check test :  

*   Qvlcy : quadratic norm velocy  
*   Mvlcy : maximum norm velocy  

usage : RBDY2_CheckEquilibrium(norm_type , tolerance)  

Parameters
----------
* `norm_type` :  
    (string of size 5) : norm type use for the equilibrium check  
* `tolerance` :  
    (double) : norm tolerance  
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
* `xmin` :  
    (double) :  
* `xmax` :  
    (double) :  
* `radius` :  
    (double) :  
";

%feature("docstring") RBDY2_UpdateThermalStrain "

usage : RBDY2_UpdateThermalStrain()  
";

%feature("docstring") RBDY2_GetNbRBDY2 "

Get the number of RBDY2.  

usage : nb_rbdy2 = RBDY2_GetNbRBDY2()  

Parameters
----------
* `nb_rbdy2` :  
    (integer) : number of RBDY2 in container  
";

%feature("docstring") RBDY2_GetBodyArea "

Get the area (2D volume equivalent) of a given body.  

usage : area = GetBodyArea(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the body  
* `area` :  
    (double) : area  
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
* `abs_min` :  
    (double) : min abscisse of sub domaine tested  
* `abs_max` :  
    (double) : max abscisse of sub domaine tested  
* `Qnorm` :  
    (double) : quadratric norm of the velocities of the grains  
* `Mnorm` :  
    (double) : quadratric norm of the velocities of the grains  
";

%feature("docstring") RBDY2_CheckPartialEquilibriumState "

Check if a part of the sample is in a equilibrium state.  

Test if there is an equilibrium state between abs_min and abs_max  

Usefull in case of silos to decide if the arch research must be activated  

usage : isPartiallyEquilibriumed = RBDY2_CheckPartialEquilibriumState(abs_min,
abs_max)  

Parameters
----------
* `abs_min` :  
    (double) : min abscisse of sub domaine tested  
* `abs_max` :  
    (double) : max abscisse of sub domaine tested  
* `isPartiallyEquilibriumed` :  
    (boolean) : true if at partial equlibrium state, else false  
";

%feature("docstring") RBDY2_SetBodiesInvisible "

Set a list of body to invisible state.  

usage : RBDY2_SetBodiesInvisible(list_bdy)  

Parameters
----------
* `list_bdy` :  
    (integer array) : list of rank of bodies of the container  
";

%feature("docstring") RBDY2_IsVisible "

return if a body visible  

usage : visible = RBDY2_IsVisible(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `visible` :  
    (integer) : 1 if body is visible, 0 else  
";

%feature("docstring") RBDY2_GetBodyMass "

Get the mass of a body.  

usage : mass = RBDY2_GetBodyMass(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of desired body  
* `mass` :  
    (double) : mass of body  
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
* `ibdyty` :  
    (integer) : rank of the RBDY2 in container  
* `density` :  
    (double) : density of the RBDY2  
";

%feature("docstring") RBDY2_GetNbContactor "

get the number of contactor of RBDY2  

python usage : nb = RBDY2_GetNbContactor(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY2 in container  

Returns
-------
nb (integer) : number of contactor attached to a RBDY2  
";

%feature("docstring") RBDY2_GetContactorType "

Get the type of the first contactor of a body.  

usage type = RBDY2_GetContactorType(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY2 in container  

Returns
-------
type (string) : type of the first contactor of the body  
";

%feature("docstring") RBDY2_GetContactorColor "

Get the color of the itacty contactor of a body ibdyty.  

usage color = RBDY2_GetContactorColor(ibdyty,itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY2 in container  
* `itacty` :  
    (integer) : rank of the contactor in the RBDY2  

Returns
-------
color (string) : color of the contactor of the body  
";

%feature("docstring") RBDY2_SetContactorColor "

Set the color of a given contactor of a body.  

usage : RBDY2_SetContactorColor(ibdyty, itacty, color)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of the RBDY2  
* `itacty` :  
    (integer) : rank of the contactor in the RBDY2  
* `color` :  
    (string of size 5) : the color  
";

%feature("docstring") RBDY2_GetPtrMass "

Get a pointer onto the mass matrix of a body.  

Parameters
----------
* `ibdyty(int)` :  
    : index of the RBDY2  
* `mass(double**)` :  
    : mass matrix of the RBDY2  
";

%feature("docstring") RBDY2_GetVelocity "

Get the velocity of a body.  

Parameters
----------
* `ibdyty(int)` :  
    : index of the RBDY2  
* `velocity(double[6])` :  
    : velocity of the RBDY2  
";

%feature("docstring") RBDY2_getDrvVlocy "

Get the driven dof of a body.  

python usage : [drvdof_indices, drvdof_values] = RBDY2_getDrvVlocy(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY2  
* `drvdof_indices` :  
    (integer array) : indices list of driven dof  
* `drvdof_values` :  
    (real array) : values of the driven dof  
";

%feature("docstring") RBDY2_computeDrvVlocy "

Compute the value of the driven velocity of a body a current time.  

In place replacement in the input array of the new value(s) of the driven
velocity  

python usage : RBDY2_computeDrvVlocy(ibdyty, values)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY2  
* `values` :  
    (double array) : numpy array, input old values of imposed velocity, output
    new ones  
";

%feature("docstring") RBDY2_SetVisible "

Set a given body as visible.  

python usage : RBDY2_SetVisible(ibdyt)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY2  
";

%feature("docstring") RBDY2_SetInvisible "

Set a given body as invisible.  

python usage : RBDY2_SetInvisible(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : index of the RBDY2  
";

%feature("docstring") RBDY2_SetVisibleVlocyDrivenDof "

allows to (re)activate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : RBDY2_SetVisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the RBDY2  
* `iccdof(integer)` :  
    : index of the DOF to set visible  
";

%feature("docstring") RBDY2_SetInvisibleVlocyDrivenDof "

allows to deactivate a given vlocydrivendof (i.e. which has been declared in
preprocessing)  

python usage : RBDY2_SetInvisibleVlocyDrivenDof(ibdyty, iccdof)  

Parameters
----------
* `ibdyty(integer)` :  
    : index of the RBDY2  
* `iccdof(integer)` :  
    : index of the DOF to set invisible  
";

%feature("docstring") RBDY2_GetBulkBehavID "

return the ID of a given bulk of a given body  

python usage : blmID = DISKx_GetBulkBehavID(ibdyty, iblmty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of a RBDY2  
* `iblmty` :  
    (integer) : rank of a bulk of the giben RBDY2 (typically 1!)  

Returns
-------
blmID (string) : the bulk behav ID  
";

%feature("docstring") RBDY2_GetBulkBehavNumber "

return the bulk ID of a given RBDY2  

python usage : ibehav = RBDY2_GetBulkBehavNumber(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of a RBDY2  

Returns
-------
ibehav (integer) : the bulk behav number  
";

%feature("docstring") RBDY2_SetSurfaceSectors "

Set the number of angular sectors of the surface of contactors.  

python usage : RBDY2_SetSurfaceSectors(nbsect)  

Parameters
----------
* `nbsect` :  
    (integer) : number of sectors  
";

%feature("docstring") RBDY2_GetStress "

Get the mean stress field of a rigid object.  

python usage : matrix = RBDY2_GetStress(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : body to get stress of  
* `matrix` :  
    (double array) : stress matrix  
";

%feature("docstring") RBDY2_ModifyBody "

Modify a body tactor.  

usage : RBDY2_ModifyBody(ibdyty, itacty, vector)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `itacty` :  
    (integer) : rank of tacty  
* `vector` :  
    (double array) : the new value  
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
* `ibdyty` :  
    (integer) : rank of body  
* `itacty` :  
    (integer) : rank of tacty  
";

%feature("docstring") RBDY2_GetElectricalPotential "

Get electrical potential of rigid particle.  

usage : ep = RBDY2_GetElectricalPotential(ibdyty, itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  

Returns
-------
ep (double) : electrical potential  
";

%feature("docstring") RBDY2_GetElectricalCurrent "

Get electrical potential of rigid particle.  

usage : ep = RBDY2_GetElectricalCurrent(ibdyty, itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  

Returns
-------
ep (double) : electrical current  
";

%feature("docstring") RBDY2_GetBetai "

Get equivalent damage related to CZM interaction for rigid particle.  

usage : betai = RBDY2_GetBetai(ibdyty, itacty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  
* `itacty` :  
    (integer) : rank of tacty  

Returns
-------
betai (double) : equivalent damage  
";

%feature("docstring") RBDY2_GetPeriode "

Get the periode id (0, 1 or -1) for rigid particles.  

usage : iper = RBDY2_GetPeriode(ibdyty)  

Parameters
----------
* `ibdyty` :  
    (integer) : rank of body  

Returns
-------
iper (integer) : periode id  
";

%feature("docstring") RBDY2_GetAverageSurfaceEnergy "
";

