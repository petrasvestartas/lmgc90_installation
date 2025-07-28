
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
* `chekctype` :  
    (char[5]) : type of convergence check  
* `tol` :  
    (double) : norm tolerance  
  
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
* `checktype` :  
    (char[5]) c : convergentce test keyword  
* `tol` :  
    (double) : tolerance value  
* `nb_iter_check` :  
    (integer) : number of iteration between convergence test  
* `nb_block_iter` :  
    (integer) : number of block iterations  
";

