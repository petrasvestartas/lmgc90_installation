
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
* `rvalue1` :  
    (double) : X coordinate  
* `rvalue2` :  
    (double) : Y coordinate  
* `rvalue3` :  
    (double) : radius selection  
";

%feature("docstring") postpro_MoveCircularSelectionZone "

Increment the position of the circular selection defined with
CIRCULAR_SELECTION.  

python usage : postpro_MoveCircularSelectionZone(rvalue1, rvalu2)  

Parameters
----------
* `rvalue1` :  
    (double) : X translational velocity  
* `rvalue2` :  
    (double) : Y translational velocity  
";

%feature("docstring") postpro_CleanMemory "

Free all memory allocated within postpro module.  

python usage : postpro_CleanMemory()  
";

%feature("docstring") postpro_2D_GetKineticEnergy "

Compute Kinetic Energy for all bodies (rigids and defo)  

python usage : KE = postpro_2D_GetKineticEnergy()  
";

