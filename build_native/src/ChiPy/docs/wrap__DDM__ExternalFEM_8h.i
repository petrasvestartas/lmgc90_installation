
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

%feature("docstring") DDM_ExternalFEM_ExSolver_3D "

Solve fully the local contact problem.  

python usage : DDM_ExternalFEM_ExSolver_3D(storage, checktype, tol, relax,
nb_iter_check, nb_block_iter)  

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

