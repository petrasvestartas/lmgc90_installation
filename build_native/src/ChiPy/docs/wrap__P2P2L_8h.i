
// File: wrap__P2P2L_8h.xml

%feature("docstring") P2P2L_SelectProxTactors "

contact detection between PT2DL tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : P2P2L_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") P2P2L_WriteLastVlocRloc "

write last local values of all P2P2L contacts  

The values written are relative velocity, forces and local frame  

python usage : P2P2L_WriteLastVlocRloc()  
";

%feature("docstring") P2P2L_WriteOutVlocRloc "

write local values of all P2P2L contacts  

The values written are relative velocity, forces and local frame  

python usage : P2P2L_WriteOutVlocRloc()  
";

%feature("docstring") P2P2L_DisplayOutVlocRloc "

display local values of all P2P2L contacts  

The values displayed are relative velocity, forces and local frame  

python usage : P2P2L_DisplayOutVlocRloc()  
";

%feature("docstring") P2P2L_DisplayProxTactors "

display contacts  

python usage : P2P2L_DisplayProxTactors()  
";

%feature("docstring") P2P2L_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : P2P2L_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") P2P2L_CleanMemory "

Free all memory allocated within P2P2L module.  

python usage : P2P2L_CleanMemory()  
";

