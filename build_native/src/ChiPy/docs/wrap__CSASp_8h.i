
// File: wrap__CSASp_8h.xml

%feature("docstring") CSASp_SelectProxTactors "

contact detection between CSxxx and ASpxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

If reset not equal to 0, the initialization flag is reset and detection skipped  

python usage : CSASp_SelectProxTactors(reset=0,use_external=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
* `use_external` :  
    (integer) : if not 0, external detection is used  
";

%feature("docstring") CSASp_WriteLastVlocRloc "

write last local values of all CSASp contacts  

The values written are relative velocity, forces and local frame  

python usage : CSASp_WriteLastVlocRloc()  
";

%feature("docstring") CSASp_WriteOutVlocRloc "

write local values of all CSASp contacts  

The values written are relative velocity, forces and local frame  

python usage : CSASp_WriteOutVlocRloc()  
";

%feature("docstring") CSASp_DisplayOutVlocRloc "

display local values of all CSASp contacts  

The values displayed are relative velocity, forces and local frame  

python usage : CSASp_DisplayOutVlocRloc()  
";

%feature("docstring") CSASp_DisplayProxTactors "

display contacts  

python usage : CSASp_DisplayProxTactors()  
";

%feature("docstring") CSASp_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being +the parameter used in
    TimeEvolution_ReadIniVlocRloc last call  

usage : CSASp_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") CSASp_SkipAutoContact "

avoid CSxxx/ASpxx contact detection when they belong to the same entity  

python usage : CSASp_SkipAutoContact()  
";

%feature("docstring") CSASp_SetNonSymmetricDetection "

this function allows non symetric detection i.e. only one interaction is kept
when two bodies with candidate and antagonist contactors see each other  

python usage : CSASp_SetNonSymmetricDetection()  
";

%feature("docstring") CSASp_Trim "

trim contact (only contact within surface - not with extremities)  

python usage : CSASp_Trim()  
";

%feature("docstring") CSASp_SetTrimAngle "

set the trim angle (only contact within surface - not with extremities)  

python usage : CSASp_SetTrimAngle(angle)  

Parameters
----------
* `angle` :  
    (real) : angle in degree - default 87 deg  
";

%feature("docstring") CSASp_AddReac "

add contact force to body Reac  

python usage : CSASp_AddReac()  
";

%feature("docstring") CSASp_AssumeOldFiles "

to read file with the CSpxx rank instead of CSxxx one  

python usage : CSASp_AssumeOldFiles()  
";

%feature("docstring") CSASp_CleanMemory "

Free all memory allocated within CSASp module.  

python usage : CSASp_CleanMemory()  
";

