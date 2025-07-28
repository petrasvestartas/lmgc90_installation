
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
* `nb_sdmx` :  
    (integer) : number of domains on x-axis  
* `nb_sdmy` :  
    (integer) : number of domains on y-axis  
* `ddm_type` :  
    (integer) : type of DDM to use  
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
* `storage` :  
    (char[30]) : matrix storage (cf nlgs_ExPrep)  
* `checktype` :  
    (char[5]) : convergentce test keyword  
* `tolerance` :  
    (double) : tolerance value  
* `relaxation` :  
    (double) : relaxation number  
* `nb_iter_check` :  
    (integer) : number of iteration between convergence test  
* `nb_block_iter` :  
    (integer) : number of block iterations  
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
* `f_ddm` :  
    (integer) : frequency of ddm partionning  
* `f_out` :  
    (integer) : frequency of output file writing  
* `f_last` :  
    (integer) : frequency of last file writing  
* `f_postpro` :  
    (integer) : frequency of postpro file writing  
* `f_display` :  
    (integer) : frequency of display file writing  
";

%feature("docstring") DDM_2D_WriteLast "

Write some data at the end of computation.  

python usage : DDM_2D_WriteLast()  
";

%feature("docstring") DDM_2D_Finalize "

End of computation/module's life.  

python usage : DDM_2D_Finalize()  
";

