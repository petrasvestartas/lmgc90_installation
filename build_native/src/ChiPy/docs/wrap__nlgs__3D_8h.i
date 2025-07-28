
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
* `chekctype_c` :  
    (char[5]) : type of convergence check  
* `tol` :  
    (double) : norm tolerance  
* `relax` :  
    (double) : relaxation factor  
  
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
* `storage_c(char[30])` :  
    : matrix storage  
  
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
* `krn` :  
    (double) : normal penality (default 1e14)  
* `krt` :  
    (double) : tangential penality (default 1e14)  
";

%feature("docstring") nlgs_3D_CutOpenCZM "

If some czm contact have a gap greater than the given they are considered as
broken ; works only with EXPO_CZM or IQS_EXPO_CZM.  

python usage : nlgs_3D_CutOpenCZM(tol)  

Parameters
----------
* `tol` :  
    (double) : threshold on positive distance (default 1e-6)  
";

%feature("docstring") nlgs_3D_ManageInterpenetratedCZM "

Apply a g0 strategy if gap is negative and if gap is positive (without using
nlgs_3D_CutOpenCZM) ; works only with EXPO_CZM or IQS_EXPO_CZM.  

python usage : nlgs_3D_ManageInterpenetratedCZM()  
";

