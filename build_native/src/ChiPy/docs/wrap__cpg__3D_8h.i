
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
* `chekctype` :  
    (char[5]) : type of convergence check  
* `tol` :  
    (double) : norm tolerance  
* `idproj` :  
    (integer) :  
  
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
* `checktype` :  
    (char[5]) : convergentce test keyword  
* `tol` :  
    (double) : tolerance value  
* `idproj` :  
    (integer) :  
* `nb_iter_check` :  
    (integer) : number of iteration between convergence test  
* `nb_block_iter` :  
    (integer) : number of block iterations  
";

