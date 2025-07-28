
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

