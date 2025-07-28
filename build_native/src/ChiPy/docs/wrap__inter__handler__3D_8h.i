
// File: wrap__inter__handler__3D_8h.xml

%feature("docstring") inter_handler_3D_tgetNb "

return the number of interactions of the selected type stored in this data
structure  

python usage : nb_inter = inter_handler_3D_tgetNb(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_3D_tgetTactLawNb "

return the contact law number of an interaction stored in this data structure  

python usage : tact_law = inter_handler_3D_tgetTactLawNb(inter_id, icdan)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  

Returns
-------
tact_law (integer) : contact law number  
";

%feature("docstring") inter_handler_3D_tgetIdBodies "

return the serial numbers of contacting objects of an interaction stored in this
data structure  

python usage : idBodies = inter_handler_3D_tgetIdBodies(inter_id, icdan)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  

Returns
-------
idBodies (integer) : array with cd and an bodies serial number  
";

%feature("docstring") inter_handler_3D_tgetIData "

Get the integer data of an interaction stored in this data structure.  

idata vector holds cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

usage : idata = inter_handler_3D_tgetIData(inter_id, icdan)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  

Returns
-------
idata (integer array) : the values array  
";

%feature("docstring") inter_handler_3D_tgetRData "

return the real data associated with an interactions  

Get an output array with, in this order, : coor, t/n/suc, rlt/n/s, vlt/n/s,
gapTT  

python usage : rdata = inter_handler_3D_tgetRData(inter_id, icdan)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  

Returns
-------
rdata (double array) : array with real data of the interaction  
";

%feature("docstring") inter_handler_3D_tsetInternal "

Set the internal of an interaction (either the array or a single value) stored
in this data structure.  

Uses copy. If internal array is provided, the whole array is set. Otherwise
index and value must be provided and a single value is set.  

usage : inter_handler_3D_tsetInternal(inter_id, icdan, internal) or
inter_handler_3D_tsetInternal(inter_id, icdan, index, value)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  
* `internal` :  
    (double array) : the new values array  
* `index` :  
    (integer) : the index where to set single value  
* `value` :  
    (double ) : the new value to put at index  
";

%feature("docstring") inter_handler_3D_tgetInternal "

Get the internal of an interaction stored in this data structure.  

usage : internal = inter_handler_3D_tgetInternal(inter_id, icdan)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdan` :  
    (integer) : index of the interaction of selected type  
* `internal` :  
    (double array) : the new values array  
";

%feature("docstring") inter_handler_3D_getNbRecup "

return the number of recup interactions of the selected type  

python usage : nb_recup = inter_handler_3D_getNbRecup(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_recup (integer) : number of interaction recup of selected type  
";

%feature("docstring") inter_handler_3D_getNb "

return the number of interactions of the selected type stored in verlet data
structure  

python usage : nb_inter = inter_handler_3D_getNb(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
nb_inter (integer) : number of interaction found of selected type  
";

%feature("docstring") inter_handler_3D_getAllTactLawNb "

return the tact law number of all interactions stored in verlet data structure  

python usage : vector = inter_handler_3D_getAllTactLawNb(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
vector (int 1D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAll "

return
coorx,coory,coorz,tx,ty,tz,nx,ny,nz,sx,sy,sz,rlt,rln,rls,vlt,vln,vls,gaptt of
all 'verlet' interactions  

python usage : array = inter_handler_3D_getAll(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAllInternal "

return contact point internal variables of all 'verlet' interactions  

python usage : array = inter_handler_3D_getAllInternal()  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (double 2D-array) : mechanical data  
";

%feature("docstring") inter_handler_3D_getAllIdata "

return all integer data of all 'verlet' interaction  

Which are in order cd body type, an body type, cd body id, an body id, cd
contactor type, an contactory type, cd contactor id, an contactor id, cd
subcontactor id, an subcontactor id, tact law id, status, number of internals  

python usage : array = inter_handler_3D_getAllIdata(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  

Returns
-------
array (int 2D-array) : identification data  
";

%feature("docstring") inter_handler_3D_getVerletAdjsz "

return integer number of verlet interaction of a candidate  

python usage : iantac = inter_handler_3D_getVerletAdjsz(inter_id, icdtac)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdtac` :  
    (integer) : candidate contactor id  

Returns
-------
iantac (integer) : number of verlet interactions on candidate  
";

%feature("docstring") inter_handler_3D_getVerletIantac "

return integer antagonist contact of a verlet interaction  

python usage : iantac = inter_handler_3D_getVerletIantac(inter_id, icdtac, iadj)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `icdtac` :  
    (integer) : candidate contactor id  
* `iadj` :  
    (integer) : id of adjacent of candidate  

Returns
-------
iantac (integer) : id of antagonist contactor corresponding to verlet
interaction  
";

%feature("docstring") inter_handler_3D_computeRnod "

Put back the Reac value of bodies from (this) interactions.  
";

%feature("docstring") inter_handler_3D_stockRloc "

stock from this to verlet  

python usage : inter_handler_3D_stockRloc(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_3D_recupRloc "

recup from verlet to this  

python usage : inter_handler_3D_recupRloc(inter_id)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
";

%feature("docstring") inter_handler_3D_recupRlocByPos "

recup from verlet to this using position as criteria  

Only available for CSASp inter_id  

python usage : inter_handler_3D_recupRloc(inter_id, rtol)  

Parameters
----------
* `inter_id` :  
    (integer) : type of interaction (lmgc90 parameter)  
* `rtol` :  
    (real) : tolerance to decide if contact is recup  
";

