
// File: wrap__surface__T3_8h.xml

%feature("docstring") surface_T3_compute_volume_inertia "

Computes the volume of an object described by a triangulated surface.  

**Warning**: 1) we assume size_coor is three times the number of nodes and
    size_connec is three times the number of elements python call: x_G, I,
    vol=surface_T3_compute_volume_inertia(coor, connec, 3, 9)  

Parameters
----------
* `coor_size` :  
    (int): size of coor  
* `coor` :  
    (double *): node coordinates  
* `connec_size` :  
    (int): size of connec  
* `vol` :  
    (double *): computed volume  
* `x_G` :  
    (double *): mass center coordinates  
* `x_G_size` :  
    (int): size of x_G  
* `I` :  
    (double *): inertia matrix, stored a a vector  
* `I_size` :  
    (int): size of I  
";

%feature("docstring") surface_T3_identify_entities "

Attributes an entity number to triangles, by computing connected components.  

**Warning**: 1) we assume size_connec is three times the number of elements and
    size_ele2entity is the number of elements python call:
    ele2entity=surface_T3_identify_entities(nbnode, max_adj_ele_2_node, connec,
    nbele)  

Parameters
----------
* `nbnode` :  
    (int): the number of nodes  
* `max_adj_ele_2_node` :  
    (int): the maximal number of adjacent elements per node  
* `connec` :  
    (int *): connecivity of elements  
* `connec_size` :  
    (int): size of connec  
* `ele2entity` :  
    (double *): entity number for each element  
* `ele2entity_size` :  
    (int): size of ele2entity  
";

