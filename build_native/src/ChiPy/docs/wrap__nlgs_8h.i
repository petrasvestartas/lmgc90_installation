
// File: wrap__nlgs_8h.xml

%feature("docstring") nlgs_ExPrep "

Prepare matrix storage.  

python usage : nlgs_ExPrep(storage)  

Parameters
----------
* `sotrage` :  
    (char[30]) : matrix storage  
  
 prepare the matrix and the RHS of the contact problem in regards of the
selected matrix storage:  

*   Exchange_Local_Global (the standard case) only the diagonal blocks are
    computed and stored.  
*   Stored_Delassus_Loops (faster but memory expensive) the complete Delassus
    matrix is computed.  
";

%feature("docstring") nlgs_ExIter "

Execute NLGS iterations over the contact loop.  

python usage : nlgs_ExIter(nb_iter) param[in] nb_iter (integer) : number of
iterations to do  
";

%feature("docstring") nlgs_ExPost "

Run a jacobi iteration with the solution obtained with the NLGS algorithm.  

python usage : nlgs_ExPost()  
";

%feature("docstring") nlgs_AfterIterCheck "

Control NLGS convergence.  

python usage : convergence = nlgs_AfterIterCheck()  

Returns
-------
convergence (integer) :  
";

%feature("docstring") nlgs_DisplayAfterIterCheck "

Display NLGS convergence results.  

python usage : nlgs_DisplayAfterIterCheck()  
";

%feature("docstring") nlgs_NormCheck "

Active one step norm evolution.  

python usage : nlgs_NormCheck()  
";

%feature("docstring") nlgs_UpdateTactBehav "

Update internal parameters of contact lawz for each contact.  

python usage : nlgs_UpdateTactBehav()  
";

%feature("docstring") nlgs_SetCheckType "

Define numerical convergence of the NLGS algorithm.  

python usage : nlgs_SetCheckType(check_type, tolerance, relaxation)  

Parameters
----------
* `check_type` :  
    (char[5]) : type of convergence check  
* `tolerance` :  
    (double) : norm tolerance  
* `relaxation` :  
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

%feature("docstring") nlgs_ScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_ScrambleContactOrder()  
";

%feature("docstring") nlgs_QuickScrambleContactOrder "

Random renumbering of the contact list.  

python usage : nlgs_QuickScrambleContactOrder()  
";

%feature("docstring") nlgs_SetWithQuickScramble "

active quick scramble in macro function ExSolver  

python usage : nlgs_SetWithQuickScramble()  
";

%feature("docstring") nlgs_ReverseContactOrder "

Reverse the numbering of the contact list.  

python usage : nlgs_ReverseContactOrder()  
";

%feature("docstring") nlgs_BimodalContactOrder "

Renumbering of the contact list using the definition of weak and strong network
in granular assemblies.  

python usage : nlgs_BimodalContactOrder()  
";

%feature("docstring") nlgs_ScaleRloc "

Scale all local contact forces of a factor equal to * 0.9 < f < 1.1.  

python usage : nlgs_ScaleRloc()  
";

%feature("docstring") nlgs_ComputeRnod "

mapping from local contact forces to global ones  

python usage : nlgs_ComputeRnod()  
";

%feature("docstring") nlgs_DisplayRlocNSum "

Display the sum of normal contact forces.  

python usage : nlgs_DisplayRlocNSum()  
";

%feature("docstring") nlgs_ExSolver "

Solve fully the local contact problem.  

python usage : nlgs_ExSolver(storage, checktype, tol, relax, nb_iter_check,
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

%feature("docstring") nlgs_UpdateCohesiveBehav "

update internal parameters of contact laws for each contact  

python usage : nlgs_UpdateCohesiveBehav(void)  
";

%feature("docstring") nlgs_UpdateFrictionalBehav "

update internal parameters of contact laws for each contact  

python usage : nlgs_UpdateFrictionalBehav(void)  
";

%feature("docstring") nlgs_GetAllThis "

Get all interactions in \"this\" array.  

Each interaction has (in this order): coor, tuc, nuc, rlt, rln, vlt, vln  

usage : interactions = nlgs_GetAllThis()  

Returns
-------
interactions (double 2D-array) : the interactions  
";

%feature("docstring") nlgs_UseJacobiSolver "

Use a Jacobi solver instead of Gauss Seidel solver.  

usage : nlgs_UseJacobiSolver(True) or nlgs_UseJacobiSolver(False)  
";

%feature("docstring") nlgs_UseRegularization "

use some regularization heuristics on interaction laws  

python usage : nlgs_UseRegularization(krn, krt)  

Parameters
----------
* `krn` :  
    (double) : normal penality (default 1e14)  
* `krt` :  
    (double) : tangential penality (default 1e14)  
";

%feature("docstring") nlgs_SetTemporaryVariable "

set temporary variables used in nlgs ; ivalue2 == 3 gives access to post crack
pressure  

python usage : nlgs_SetTemporaryVariable(icdan,id,val)  

Parameters
----------
* `icdan` :  
    (int) : interaction rank  
* `id` :  
    (int) : value rank  
* `val` :  
    (double) : value  
";

%feature("docstring") nlgs_GetTemporaryVariable "

get temporary variables used in nlgs ; ivalue2 == 3 gives access to post crack
pressure  

python usage : val = nlgs_GetTemporaryVariable(icdan,id)  

Parameters
----------
* `icdan` :  
    (int) : interaction rank  
* `id` :  
    (int) : value rank  
* `val` :  
    (double) : value  
";

%feature("docstring") nlgs_IsInitialized "

In case of restart say that nlgs is initialized or reset it.  

python usage : nlgs_IsInitialized(is_init=1)  
";

