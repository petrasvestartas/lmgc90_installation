
// File: wrap__PRPRx_8h.xml

%feature("docstring") PRPRx_SelectProxTactors "

contact detection between PRxxx and PRxxx tactors  

First recup coordinate prediction, then proceed to a box selection to found
rough contact list and finally compute the final contact list.  

python usage : PRPRx_SelectProxTactors(reset=0)  

Parameters
----------
* `reset` :  
    (integer) : if not 0, detection is skipped but the boxes will be computed
    anew at next call  
";

%feature("docstring") PRPRx_UseCpCundallDetection "

chooses the Cundall iterative detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpCundallDetection(nb_iter, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
* `nb_iter` :  
    (integer) : max number of iterations  
* `cd_shrink` :  
    (real) : shrink parameter (length) in clipper for candidate  
* `an_shrink` :  
    (real) : shrink parameter (length) in clipper for antagonist  
* `delta` :  
    (real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseCpF2fExplicitDetection "

chooses the face 2 face combinatory detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpF2fExplicitDetection(tol, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
* `tol` :  
    (real) : tolerance on normal orientations  
* `cd_shrink` :  
    (real) : shrink parameter (length) in clipper for candidate  
* `an_shrink` :  
    (real) : shrink parameter (length) in clipper for antagonist  
* `delta` :  
    (real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseCpF2fDetection "

chooses a mix of the face 2 face and Cundall detection method  

If shrink parameters are provided they may be conflicting with a call to
PRPRx_ShrinkPolyrFaces function. Remind that that the shrink parameters provided
here are lengths.  

python usage : PRPRx_UseCpF2fDetection(tol, iter, cd_shrink=0., an_shrink=0.,
delta=0.)  

Parameters
----------
* `tol` :  
    (real) : tolerance on normal orientations  
* `iter` :  
    (integer) : max number of iterations  
* `cd_shrink` :  
    (real) : shrink parameter (length) in clipper for candidate  
* `an_shrink` :  
    (real) : shrink parameter (length) in clipper for antagonist  
* `delta` :  
    (real) : intersection simplification parameter in clipper  
";

%feature("docstring") PRPRx_UseNcDetection "

chooses contact detection methode between non-convex shapes  

python usage : PRPRx_UseNcDetection(gdist)  

Parameters
----------
* `gdist` :  
    (real) : global distance  
";

%feature("docstring") PRPRx_UseNcF2fDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

python usage : PRPRx_UseNcF2fDetection(gdist,tol)  

Parameters
----------
* `gdist` :  
    (real) : global distance  
* `tol` :  
    (real) : tolerance on normal orientations  
";

%feature("docstring") PRPRx_UseNcF2fExplicitDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

python usage : PRPRx_UseNcF2fExplicitDetection(gdist,tol)  

Parameters
----------
* `gdist` :  
    (real) : global distance  
* `tol` :  
    (real) : tolerance on normal orientations  
";

%feature("docstring") PRPRx_UseTrianglesIntersectionDetection "

chooses contact detection finding intersection in a soup of triangles.  

The number of point provided is an internal parameter of the algorithm which
control the maximum number of intersection points stored when looking for the
triangles intersection before restricting it to only 4 of them. So it must be
strictly superior to 4.  

python usage : PRPRx_UseTrianglesIntersectionDetection(nb_max_pt=16)  

Parameters
----------
* `nb_max_pt(integer)` :  
    : maximum contact points to store/check during detection  
";

%feature("docstring") PRPRx_SetF2fMinimalSurfaceSize "

set the minimum contact surface size with f2f algo otherwize contact is not
computed  

python usage : PRPRx_SetF2fMinimalSurfaceSize(tol)  

Parameters
----------
* `tol` :  
    (real) : minimum surface size  
";

%feature("docstring") PRPRx_UseExternalDetection "

chooses external contact detection (bindings)  

python usage : PRPRx_UseExternalDetection()  
";

%feature("docstring") PRPRx_WriteLastVlocRloc "

write last local values of all PRPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : PRPRx_WriteLastVlocRloc()  
";

%feature("docstring") PRPRx_WriteOutVlocRloc "

write local values of all PRPRx contacts  

The values written are relative velocity, forces and local frame  

python usage : PRPRx_WriteOutVlocRloc()  
";

%feature("docstring") PRPRx_DisplayOutVlocRloc "

display local values of all PRPRx contacts  

The values displayed are relative velocity, forces and local frame  

python usage : PRPRx_DisplayOutVlocRloc()  
";

%feature("docstring") PRPRx_DisplayProxTactors "

display contacts  

python usage : PRPRx_DisplayProxTactors()  
";

%feature("docstring") PRPRx_ReadIniVlocRloc "

Read VlocRloc file.  

*   If num <= 0 : DATBOX/VlocRloc.INI file is read  
*   Else : OUTBOX/VlocRloc.OUT.num is read, num being
    -   the parameter used in TimeEvolution_ReadIniVlocRloc last call  

python usage : PRPRx_ReadIniVlocRloc(num=0)  

Parameters
----------
* `num` :  
    (integer) : which VlocRloc file to read  
";

%feature("docstring") PRPRx_ShrinkPolyrFaces "

Shrink the face of the candidate polyhedron for the detection.  

May be conflicting with the shrink parameters of the detections functions used
by clipper library. The difference is that clipper use a single length for all
sample, whereas this function use a scale factor to retract the vertices of the
candidate polyhedron inside the the surface.  

python usage : PRPRx_ShrinkPolyrFaces(shrink)  

Parameters
----------
* `shrink` :  
    (real) : scale factor allowing to shrink candidate surface  

    1.  no shrink, 1. no surface  
";

%feature("docstring") PRPRx_LowSizeArrayPolyr "

abscons parameter to manage memory allocation  

python usage : PRPRx_LowSizeArrayPolyr(sfactor)  

Parameters
----------
* `sfactor` :  
    (integer) :  
";

%feature("docstring") PRPRx_SaveProxTactorsToFile "

write selected contacts to file  

python usage : PRPRx_SaveProxTactorsToFile()  
";

%feature("docstring") PRPRx_LoadProxTactorsFromFile "

load selected contact from files  

python usage : PRPRx_LoadProxTactorsFromFile()  
";

%feature("docstring") PRPRx_SetXPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PRPRx_SetXPeriodicCondition(xperiod)  

Parameters
----------
* `xperiod` :  
    (real) : periode on x axis  
";

%feature("docstring") PRPRx_SetYPeriodicCondition "

initialise data for simulation using periodic condition  

python usage : PRPRx_SetYPeriodicCondition(yperiod)  

Parameters
----------
* `yperiod` :  
    (real) : period on y axis  
* `yperiod` :  
    (double) : period on y axis  
";

%feature("docstring") PRPRx_VerboseF2F "

ask for verbose comment concerning contact detection between cd and an  

python usage : PRPRx_VerboseF2F(cd,an)  

Parameters
----------
* `cd` :  
    (integer) : candidate  
* `an` :  
    (integer) : antagoniste  
";

%feature("docstring") PRPRx_GetNbF2f "

Get the number of f2f structures stored This is the real size of the array, and
not the number of active f2f structure.  

python usage : nb_f2f = PRPRx_GetNbF2f()  

Returns
-------
nb_f2f (integer) : the size of the f2f array  
";

%feature("docstring") PRPRx_GetF2f2Inters "

Get the list of interactions for each face-to-face structure Array of integer
with number of f2f, then for each f2f, the number of interactions then the list
of interaction id.  

python usage : f2f_inters = PRPRx_GetF2f2Inters()  

Returns
-------
f2f_inters (integer array) : the integer array  
";

%feature("docstring") PRPRx_GetF2fOutlines "

Get the connectivity of all intersection polytopes of all face2face and the
corresponding coordinates.  

The connectivity containes first the number of f2f, then for each, the number of
polytope, then for each the number of vertices.  

The coordinates must be counted from this ordering...  

python usage : connec, points = PRPRx_GetF2fOutlines()  

Returns
-------

*   connec (integer array) : the connectivities  
*   points (double array) : the coordinates  
";

%feature("docstring") PRPRx_GetF2fAllIdata "

Get topological face id of cd/an for all F2f structure.  

python usage : idata = PRPRx_GetF2fAllIdata()  

Returns
-------

*   idata (integer array) : size [nb_f2fx2] with the face id  
";

%feature("docstring") PRPRx_GetF2fCentralKernel "

Give the central kernel coordinates, the equivalent normal stress and if the
center of pressure is inside.  

python usage : ck_coor, sn, is_in = PRPRx_GetF2fStress(i_f2f)  

Returns
-------  
";

%feature("docstring") PRPRx_GetF2fStress "

Give the polygons of the compressed and decompressed part and linear stress
repartition.  

In the case when the minimization algorithm failed, the decompression value is
set to -99. so that when writing the vtk files, the 'ids' numbering is kept
consistent.  

python usage : coorC, sizeC, coorD, sizeD, sigma, decomp =
PRPRx_GetF2fStress(i_f2f)  

Returns
-------  
";

%feature("docstring") PRPRx_SetCundallNeighbor "

set a neighbor distance around common plane to select projected nodes  

python usage : PRPRx_SetCundallNeighbor(neighbor)  

Parameters
----------
* `neighbor` :  
    (real) : ratio of a reference size  
";

%feature("docstring") PRPRx_CpUseOldCcpm "

use the old method for computing contact point position  

python usage : PRPRx_CpUseOldCcpm()  
";

%feature("docstring") PRPRx_SetReactionTrackingLength "

function which makes possible to set the length of the hexaedra glyph
representing the visavis reaction  

python usage : PRPRx_SetReactionTrackingLength(length)  

Parameters
----------
* `length` :  
    (real) : length the hexaedra glyph  
";

%feature("docstring") PRPRx_SetTolRecupRloc "

set the distance tolerance used in PRPRx_RecupRloc  

python usage : PRPRx_SetTolRecupRloc(tol)  

Parameters
----------
* `tol` :  
    (double) : tolerance  
";

%feature("docstring") PRPRx_GetInteractionVector "

Get a copy of a vector of a PRPRx.  

possible values for datatype field are \"Coor_\", \"N____\"  

python usage : vector = PRPRx_GetInteractionVector(datatype, icdan)  

Parameters
----------
* `datatype` :  
    (string [5]) : the vector to get  
* `icdan` :  
    (integer) : rank of the PRPRx  

Returns
-------
vector (double array) : output vector  
";

%feature("docstring") PRPRx_SetInteractionInternal "

Set a value of the internal vector of a PRPRx.  

python usage : PRPRx_SetInteractionInternal(i, icdan, value)  

Parameters
----------
* `i` :  
    (integer) : rank of internal  
* `icdan` :  
    (integer) : rank of the PRPRx  
* `value` :  
    (double) : value to set  
";

%feature("docstring") PRPRx_GetInteractionInternal "

Get a value from the internal vector of a PRPRx.  

python usage : value = PRPRx_GetInteractionInternal(i, icdan)  

Parameters
----------
* `i` :  
    (integer) : rank of internal  
* `icdan` :  
    (integer) : rank of the PRPRx  
* `value` :  
    (double) : value to get  
";

%feature("docstring") PRPRx_GetInteractionInternalComment "

Get internal comment of a given interaction.  

python usage : comment=PRPRx_GetInteractionInternalComment(icdan)  

Parameters
----------
* `icdan` :  
    (integer) : rank of the PRPRx  

Returns
-------
comment (char[100]) : the string to get  
";

%feature("docstring") PRPRx_WithNodalContact "

use cd contact points at nodes instead at faces with NcDetection  

python usage : PRPRx_WithNodalContact()  
";

%feature("docstring") PRPRx_SetInternalSurface "

Set the value of a surface type (point, line or surf) for wti detection.  

For surface, if the value is left to 0., then the surface of the triangle is
computed To select the type of surface : 1->point, 2->line, 3->surface  

python usage : PRPRx_SetInternalSurface(itype, value)  

Parameters
----------
* `itype` :  
    (integer) : the type of surface to set  
* `value` :  
    (double) : value to set  
";

%feature("docstring") PRPRx_UseStoDetection "

chooses contact detection between between non-convex shapes using f2f strategy  

Face to face detection implemented by Stono which can mix between the standard
f2f detection and the non convex one. Furthermor the decompression parameter can
help with putting the contact points either near the  

python usage : PRPRx_UseFCDetection(explicite, decompression, tol, kappa)  

Parameters
----------
* `explicite` :  
    (boolean) : use explicit detection  
* `decompression` :  
    (double) : surface decompression (value in [-1., 1.])  
* `tol` :  
    (real) : tolerance on normal orientations  
* `kappa` :  
    (boolean) : compute kappas coefficient  
";

%feature("docstring") PRPRx_ForceF2fDetection "

force f2f detection method even for non-convex surfaces  

python usage : PRPRx_ForceF2fDetection()  
";

%feature("docstring") PRPRx_ForceNcDetection "

force nc detection method even for flat surfaces  

python usage : PRPRx_ForceNcDetection()  
";

%feature("docstring") PRPRx_CleanMemory "

Free all memory allocated within PRPRx module.  

python usage : PRPRx_CleanMemory()  
";

