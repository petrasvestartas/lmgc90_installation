
// File: wrap__SPPLx_8h.xml

%feature("docstring") SPPLx_SelectProxTactors "

contact detection between SPxxx and PLxxx tactors  

python usage : SPPLx_SelectProxTactors(reset=0) param[in] reset (integer) : if
not 0, detection is skipped but the boxes will be computed anew at next call  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  
";

%feature("docstring") SPPLx_SmoothForceComputation "

compute smooth contact law (in any)  

python usage : SPPLx_SmoothForceComputation()  
";

%feature("docstring") SPPLx_WriteLastVlocRloc "

write last local values of all SPPLx contacts  

python usage : SPPLx_WriteLastVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_WriteOutVlocRloc "

write local values of all SPPLx contacts  

python usage : SPPLx_WriteOutVlocRloc()  

  
 the values written are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_DisplayOutVlocRloc "

display local values of all SPPLx contacts  

python usage : SPPLx_DisplayOutVlocRloc()  

  
 the values displayed are relative velocity, forces and local frame  
";

%feature("docstring") SPPLx_DisplayProxTactors "

display contacts  

python usage : SPPLx_DisplayProxTactors()  
";

%feature("docstring") SPPLx_ReadIniVlocRloc "

Read VlocRloc file.  

If num <= 0 : DATBOX/VlocRloc.INI file is read Else : OUTBOX/VlocRloc.OUT.num is
read, num being the parameter used in TimeEvolution_ReadIniVlocRloc last call  

usage : SPPLx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") SPPLx_CleanMemory "

Free all memory allocated within SPPLx module.  

python usage : SPPLx_CleanMemory()  
";

