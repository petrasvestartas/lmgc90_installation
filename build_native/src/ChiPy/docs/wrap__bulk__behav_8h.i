
// File: wrap__bulk__behav_8h.xml

%feature("docstring") bulk_behav_ReadBehaviours "

read gravity and behaviors from DATBOX/BULK_BEHAV.DAT file  

python usage : bulk_behav_ReadBehaviours()  
";

%feature("docstring") bulk_behav_WriteBehaviours "

write gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_WriteBehaviours()  
";

%feature("docstring") bulk_behav_CollectOutBulkBehav "

read gravity and behaviors from OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_CollectOutBulkBehav()  
";

%feature("docstring") bulk_behav_CleanOutBulkBehav "

write (replacing) gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_CleanOutBulkBehav()  
";

%feature("docstring") bulk_behav_AppendOutBulkBehav "

write (appending) gravity and behaviors to OUTBOX/BULK_BEHAV.OUT file  

python usage : bulk_behav_AppendOutBulkBehav()  
";

%feature("docstring") bulk_behav_RebuildInBulkBehav "

write (replace) gravity and behaviors to DATBOX/BULK_BEHAV.DAT file  

python usage : bulk_behav_RebuildInBulkBehav()  
";

%feature("docstring") bulk_behav_GetGravity "

get the gravity acceleration used  

python usage : gravity = bulk_behav_GetGravity()  

Returns
-------
gravity (double array) : gravity vector  
";

%feature("docstring") bulk_behav_SetGravity "

set the gravity acceleration to be used  

python usage : bulk_behav_SetGravity(gravity)  

Parameters
----------
* `gravity` :  
    (double array) : gravity vector (size 3)  
";

%feature("docstring") bulk_behav_SetConductivity "

set the conductivity parameter to be used  

python usage : bulk_behav_SetConductivity(cvalue ,ivalue, rvalue)  

Parameters
----------
* `cvalue` :  
    (string of size 5) : nickname of bulk behaviour  
* `ivalue` :  
    (integer) : type of parameter: 0 = constant, 1 = field  
* `rvalue` :  
    (real) : conductivity value  
";

%feature("docstring") bulk_behav_SetCapacity "

set the Capacity parameter to be used  

python usage : bulk_behav_SetCapacity(cvalue ,ivalue, rvalue)  

Parameters
----------
* `cvalue` :  
    (string of size 5) : nickname of bulk behaviour  
* `ivalue` :  
    (integer) : type of parameter: 0 = constant, 1 = field  
* `rvalue` :  
    (real) : Capacity value  
";

%feature("docstring") bulk_behav_SetBiot "

set the Biot parameter to be used  

python usage : bulk_behav_SetBiot(cvalue ,ivalue, rvalue)  

Parameters
----------
* `cvalue` :  
    (string of size 5) : nickname of bulk behaviour  
* `ivalue` :  
    (integer) : type of parameter: 0 = constant, 1 = field  
* `rvalue` :  
    (real) : Biot value  
";

%feature("docstring") bulk_behav_SetExternalFlux "

set the External Flux parameter to be used  

python usage : bulk_behav_SetExternalFlux(cvalue ,ivalue, rvalue)  

Parameters
----------
* `cvalue` :  
    (string of size 5) : nickname of bulk behaviour  
* `ivalue` :  
    (integer) : type of parameter: 0 = constant, 1 = field  
* `rvalue` :  
    (real) : External Flux value  
";

%feature("docstring") bulk_behav_SetDensity "

set the Density parameter to be used  

python usage : bulk_behav_SetDensity(cvalue , rvalue)  

Parameters
----------
* `cvalue` :  
    (string of size 5) : nickname of bulk behaviour  
* `rvalue` :  
    (real) : Density value  
";

%feature("docstring") bulk_behav_GetNbBulkBehav "

get the number of bulk laws  

python usage : nb_bulk_behav = bulk_behav_GetNbBulkBehav()  

Parameters
----------
* `nb_bulk_behav` :  
    (integer) : number of bulk behaviour in lmgc90  
";

%feature("docstring") bulk_behav_GetBulkBehav "

get a given bulk law  

python usage : lawty, behav = bulk_behav_GetBulkBehav(i_bb)  

Parameters
----------
* `i_bb` :  
    (integer) : index of the desired bulk_behav  
* `lawty` :  
    (string) : type of the bulk law  
* `behav` :  
    (string) : name of the bulk law  
* `param` :  
    (real vector) : parameters of the law  
";

%feature("docstring") bulk_behav_CleanMemory "

Free all memory allocated within bulk_behav module.  

python usage : bulk_behav_CleanMemory()  
";

