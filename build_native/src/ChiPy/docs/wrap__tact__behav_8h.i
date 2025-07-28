
// File: wrap__tact__behav_8h.xml

%feature("docstring") tact_behav_OpenBehavContainer "

open the container (access as a linked list) in order to add/remove objects  

python usage : tact_behav_OpenBehavContainer()  
";

%feature("docstring") tact_behav_CloseBehavContainer "

close the container (access as an array)  

python usage : tact_behav_TactBehavContainer()  
";

%feature("docstring") tact_behav_OpenSeeContainer "

open the container (access as a linked list) in order to add/remove objects  

python usage : tact_behav_OpenSeeContainer()  
";

%feature("docstring") tact_behav_CloseSeeContainer "

close the container (access as an array)  

python usage : tact_behav_CloseSeeContainer()  
";

%feature("docstring") tact_behav_FillContainersFromFile "

read DATBOX/TACT_BEHAV.DAT and fill the containers (see and tact)  

python usage : tact_behav_FillContainersFromFile()  
";

%feature("docstring") tact_behav_AddToSeeContainer "

add a see table to the container  

python usage :
tact_behav_AddToSeeContainer(cdbdy,cdtac,cdcol,behav,anbdy,antac,ancol,alert,global_alert)  
";

%feature("docstring") tact_behav_ReadBehaviours "

open + fill + close  

python usage : tact_behav_ReadBehaviours()  
";

%feature("docstring") tact_behav_CollectOutTactBehav "

old fashion read from OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_CollectOutTactBehav()  
";

%feature("docstring") tact_behav_WriteBehaviours "

write (replace) tact and see to OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_WriteBehaviours()  
";

%feature("docstring") tact_behav_AppendOutTactBehav "

write (append) tact and see to OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_AppendOutTactBehav()  
";

%feature("docstring") tact_behav_RebuildInTactBehav "

write (replace) tact and see to DATBOX/TACT_BEHAV.DAT  

python usage : tact_behav_RebuildInTactBehav()  
";

%feature("docstring") tact_behav_CleanOutTactBehav "

erase OUTBOX/TACT_BEHAV.OUT  

python usage : tact_behav_CleanOutTactBehav()  
";

%feature("docstring") tact_behav_GetNbTactBehav "

get the number of tact laws  

python usage : nb_tact_behav = tact_behav_GetNbTactBehav()  

Parameters
----------
* `nb_tact_behav` :  
    (integer) : number of contact behaviour in lmgc90  
";

%feature("docstring") tact_behav_GetTactBehav "

get information related to a given tact law  

python usage : [lawty, behav, param] = tact_behav_GetTactBehav(i_tb)  

Parameters
----------
* `i_tb` :  
    (integer) : rank (in the contact laws list) of the desired tact_behav  
* `lawty` :  
    (string) : type of the contact law  
* `behav` :  
    (string) : name of the contact law  
* `param` :  
    (real vector) : parameters of the law  
";

%feature("docstring") tact_behav_GetInternalComment "

Get internal variables comment of a given interaction law.  

python usage : comment = tact_behav_GetInternalComment(ilaw)  

Parameters
----------
* `ilaw` :  
    (integer) : rank of the interaction law  

Returns
-------
comment (char[100]) : the string to get  
";

%feature("docstring") tact_behav_SetCZMwithInitialFriction "

define the way friction evolve with damage: =0. constant value, (1. - beta)**pow
otherwize  

python usage : tact_behav_SetCZMwithInitialFriction(pow)  

Parameters
----------
* `pow` :  
    (real) : parameter of power law evlution for friction
    mu(beta)=mu_s*(1-beta)**pow  
";

%feature("docstring") tact_behav_initFrictionEvolution "

[experimental] read a friction time evolution map  

python usage : tact_behav_initFrictionEvolution()  
";

%feature("docstring") tact_behav_setRandomFriction "

Active variation of local friction.  

python usage : tact_behav_setRandomFriction(r8)  
";

%feature("docstring") tact_behav_GetTactBehavRankFromName "

get the rank (in the list of tact laws) of a tact behav law  

python usage : rank = tact_behav_GetTactBehavRankFromName(c5)  
";

%feature("docstring") tact_behav_GetParamRankFromName "

get the rank of a param for a given tact behav law  

python usage : rank = tact_behav_GetParamRankFromName(i_tact,c5)  
";

%feature("docstring") tact_behav_GetParam "

get the value of a parameter  

python usage : param = tact_behav_GetParam(i_tact,i_param)  

Parameters
----------
* `i_tact` :  
    (integer) : rank of the interaction law  
* `i_param` :  
    (integer) : rank of the parameter  
* `param` :  
    (real ) : value of the parameter  
";

%feature("docstring") tact_behav_SetParam "

set the value ...  

python usage : tact_behav_SetParam(i_tact, i_param, param)  

Parameters
----------
* `i_tact` :  
    (integer) : rank of the interaction law  
* `i_param` :  
    (integer) : rank of the parameter  
* `param` :  
    (real ) : value of the parameter  
";

%feature("docstring") tact_behav_GetLawInternalComment "
";

%feature("docstring") tact_behav_SetRNcap "

set a maximal compression value  

python usage : tact_behav_SetRNcap(param)  
";

%feature("docstring") tact_behav_SetDilatancyParameters "

set dilatancy parameters  

python usage : tact_behav_SetDilatancyParameters(fric,height)  
";

%feature("docstring") tact_behav_SetPressureParameters "

set pressure parameters  

Parameters
----------
* `ibehav` :  
    (integer) : rank of the tact behav  
* `flag` :  
    (integer) : kind of build-in pressure law (0 no pressure, 1: time dependent,
    2: linearly progressive since crack starts, 3: exponentially progressive
    since crack starts, 4: external)  
* `params` :  
    (double array) : the new value of the params [p0,dp,tau,alpha]  

python usage : tact_behav_SetPressureParameters(ibehav,flag,params)  
";

%feature("docstring") tact_behav_CleanMemory "

Free all memory allocated within tact_behav module.  

python usage : tact_behav_CleanMemory()  
";

