
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
* `dt` :  
    (double) : value of time step  
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
* `first_step` :  
    (integer) : rank of the first time step  
";

%feature("docstring") TimeEvolution_SetInitialTime "

Set initial time.  

python usage : TimeEvolution_SetInitialTime(t_init)  

Parameters
----------
* `t_init` :  
    (double) : initial time  
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
* `Nstep_writeDof` :  
    (integer) : periodicity of DOF write  
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
* `nstep` :  
    (integer) : a freq of writing  
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
* `nstep` :  
    (integer) : a freq of writing  
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
* `nstep` :  
    (integer) : a freq of writing  
";

%feature("docstring") TimeEvolution_ReadIniDof "

Read header of a DOF file.  

python usage : TimeEvolution_ReadIniDof(num=0)  

Parameters
----------
* `num` :  
    (integer) : num of file to read  
";

%feature("docstring") TimeEvolution_ReadIniVlocRloc "

Read header of a VlocRloc file.  

python usage : TimeEvolution_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : num of file to read  
";

%feature("docstring") TimeEvolution_ReadIniGPV "

Read header of a GPV file.  

python usage : TimeEvolution_ReadIniGPV(num=0)  

Parameters
----------
* `num` :  
    (integer) : num of file to read  
";

%feature("docstring") NewtonRaphson_Initialize "

initialize Newton Raphson Loop  

python usage : NewtonRaphson_Initialize(tol)  

Parameters
----------
* `tol` :  
    (double) : tolerance  
";

%feature("docstring") NewtonRaphson_CheckConvergence "

check if Newton Raphson loop converges  

python usage : iconv = NewtonRaphson_CheckConvergence(norm)  

Parameters
----------
* `norm` :  
    (double) : value to check  

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
* `dt` :  
    (double) : minimum value of time step  
  
 Needed only if adaptive time step feature is used  
";

%feature("docstring") NewtonRaphson_SetMaxTimeStep "

Set value of the maximum possible time step.  

python usage : NewtonRaphson_SetMaxTimeStep(dt)  

Parameters
----------
* `dt` :  
    (double) : maximum value of time step  
  
 Needed only if adaptive time step feature is used  
";

%feature("docstring") NewtonRaphson_SetFinalTime "

Set final time.  

python usage : NewtonRaphson_SetFinalTime(t_final)  

Parameters
----------
* `t_final` :  
    (double) : final time  
";

%feature("docstring") NewtonRaphson_SetMaxIter "

Max number of iterations - default is 50.  

python usage : NewtonRaphson_SetMaxIter(max_iter)  

Parameters
----------
* `max_iter` :  
    (integer) :  
";

%feature("docstring") NewtonRaphson_SetGoodIter "

Set the max number of iterations for good convergence - default is 10.  

python usage : NewtonRaphson_SetGoodIter(good_iter)  

Parameters
----------
* `good_iter` :  
    (integer) :  
";

%feature("docstring") NewtonRaphson_SetBadIter "

Set the max number of iterations for bad convergence - default is 30.  

python usage : NewtonRaphson_SetBadIter(bad_iter)  

Parameters
----------
* `bad_iter` :  
    (integer) :  
";

%feature("docstring") NewtonRaphson_SetIncPatience "

Set the number of increments to adapt the time step when successive good
convergence (increase time step) or bad convergence (decrease time step) -
default is 3.  

python usage : NewtonRaphson_SetIncPatience(patience)  

Parameters
----------
* `patience` :  
    (integer) :  
";

%feature("docstring") overall_SelectProxTactors "

Prepare contact detection.  

python usage : overall_SelectProxTactors(Nstep_rough_seek)  

Parameters
----------
* `Nstep_rough_seek` :  
    (integer) : periodicity of rough detection  
";

%feature("docstring") overall_DisplayProxTactors "

python usage : overall_DisplayProxTactors()  
";

%feature("docstring") overall_DIME "

set space dimension and in 2D the modelling assumption  

python usage : overall_DIME(idim, imod)  

Parameters
----------
* `idim` :  
    (integer) : dimension (2 or 3)  
* `imod` :  
    (integer) : kind of model (2D only)  

*   imod = 1 => plane strain  
*   imod = 2 => plane stress  
*   imod = 3 => axisymmetric  
";

%feature("docstring") Integrator_InitTheta "

python usage : Integrator_InitTheta(theta)  

Parameters
----------
* `theta` :  
    (double) : value of theta in integrator  
";

%feature("docstring") Integrator_InitQS "

python usage : Integrator_InitQS()  
";

%feature("docstring") Integrator_InitCrankNickolson "

python usage : Integrator_InitCrankNickolson(theta)  

Parameters
----------
* `theta` :  
    (double) : value of theta in integrator  
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
* `value` :  
    (double) : numeric diffusion ([0.5,1] and 0.5 is conservative)  
";

%feature("docstring") Integrator_SetContactDetectionConfiguration "

set the parameters necessary to define the contact detection configuration
(default: 1-theta, 0.)  

python usage : Integrator_SetContactDetectionConfiguration(alpha_b,alpha_e)  

Parameters
----------
* `alpha_b` :  
    (double) : value of the V_begin weight  
* `alpha_e` :  
    (double) : value of the V weight  
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
* `ifirst` :  
    (integer) :  
* `ilast` :  
    (integer) :  
";

%feature("docstring") overall_SetWorkingDirectory "

python usage : overall_SetWorkingDirectory(path)  

Parameters
----------
* `path` :  
    (string) : set path to DATBOX directory  
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
* `freq_display` :  
    (integer) : periodicity of display write  
";

%feature("docstring") TimeEvolution_ReadIniMpValues "

Read header of a MP_VALUES file.  

python usage : TimeEvolution_ReadIniMpValues(num=0)  

Parameters
----------
* `num` :  
    (integer) : num of file to read  
";

%feature("docstring") TimeEvolution_WriteOutMpValues "

python usage : TimeEvolution_WriteOutMpValues(nstep)  

Parameters
----------
* `nstep` :  
    (integer) : a freq of writing  
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

