
// File: wrap__mesh2D_8h.xml

%feature("docstring") mesh2D_GetIndicesMeshQ4 "

this function gives the couple (i, j) of indices coresponding to a given node n  

**Warning**: python call: [i, j]=mesh2D_GetIndicesMeshQ4(n)  

Parameters
----------
* `n` :  
    (int): the given node  
* `i` :  
    (int *): index in the u direction  
* `j` :  
    (int *): index in the v direction  
";

%feature("docstring") mesh2D_SizeMeshQ4 "

this function computes the sizes of vectors used to store a mesh made of Q4 in
the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4 6  *---*---*  
     | 1 | 2 |  *---*---*  
     1 3 5  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]  
*   nb_node_per_ele = [4, 4]  
*   conn = [1, 3, 4, 2, 3, 5, 6, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMeshQ4(nb_elem_x, nb_elem_y)  

    Parameters:  
    * `nb_elem_x` :  
        (int): number of elements in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements in the vertical direction  
    * `size_coor` :  
        (int *): size of coor  
    * `size_nb_node_per_ele` :  
        (int *): size of nb_node_per_ele  
    * `size_conn` :  
        (int *): size of conn  
";

%feature("docstring") mesh2D_SizeMesh2T3 "

this function computes the sizes of vectors used to store a mesh made of T3 ---
obtained by splitting a Q4 in two T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *----*  
     | 1 /|  
     | / |  
     | / |  
     |/ 2 |  *----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4]  
*   nb_node_per_ele = [3, 3]  
*   conn = [1, 3, 4, 2, 1, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh2T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int *): size of coor  
    * `size_nb_node_per_ele` :  
        (int *): size of nb_node_per_ele  
    * `size_conn` :  
        (int *): size of conn  
";

%feature("docstring") mesh2D_SizeMesh4T3 "

this function computes the sizes of vectors used to store a mesh made of T3 ---
obtained by splitting a Q4 in four T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *-----*  
     |\\ 4 /|  
     | \\ / |  
     |1 5 3|  
     | / \\ |  
     |/ 2 |  *-----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5]  
*   nb_node_per_ele = [3, 3, 3, 3]  
*   conn = [1, 5, 2, 1, 3, 5, 3, 4, 5, 2, 5, 4]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh4T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int *): size of coor  
    * `size_nb_node_per_ele` :  
        (int *): size of nb_node_per_ele  
    * `size_conn` :  
        (int *): size of conn  
";

%feature("docstring") mesh2D_SizeMeshQ8 "

this function computes the sizes of vectors used to store a mesh made of Q8
following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     3 7 3  *---*---*  
     | |  
     | |  
     8 * 1 * 6  
     | |  
     | |  *---*---*  
     1 5 2  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8]  
*   nb_node_per_ele = [8]  
*   conn = [1, 2, 3, 4, 5, 6, 7, 8]  

    **Warning**: python call: [size_coor, size_nb_node_per_ele,
        size_conn]=mesh2D_SizeMesh4T3(nb_elem_x, nb_elem_y)  

    Parameters:  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int *): size of coor  
    * `size_nb_node_per_ele` :  
        (int *): size of nb_node_per_ele  
    * `size_conn` :  
        (int *): size of conn  
";

%feature("docstring") mesh2D_MeshQ4 "

this function computes and returns a mesh made of Q4 in the following generic
format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4 6  *---*---*  
     | 1 | 2 |  *---*---*  
     1 3 5  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]  
*   nb_node_per_ele = [4, 4]  
*   conn = [1, 3, 4, 2, 3, 5, 6, 4]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_MeshQ4(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    * `x0` :  
        (double): abscissa of the lower left corner of the rectangle  
    * `y0` :  
        (double): ordinate of the lower left corner of the rectangle  
    * `lx` :  
        (double): length of the mesh, following the axis (Ox)  
    * `ly` :  
        (double): length of the mesh, following the axis (Oy)  
    * `nb_elem_x` :  
        (int): number of elements in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements in the vertical direction  
    * `size_coor` :  
        (int): size of coor  
    * `size_nb_node_per_ele` :  
        (int): size of nb_node_per_ele  
    * `size_conn` :  
        (int): size of conn  
    * `coor` :  
        (double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    * `nb_node_per_ele` :  
        (int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    * `conn` :  
        (int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_Mesh2T3 "

this function computes an returns a mesh made of T3 --- obtained by splitting a
Q4 in two T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *----*  
     | 1 /|  
     | / |  
     | / |  
     |/ 2 |  *----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4]  
*   nb_node_per_ele = [3, 3]  
*   conn = [1, 3, 4, 2, 1, 4]  
*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_Mesh2T3(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    * `x0` :  
        (double): abscissa of the lower left corner of the rectangle  
    * `y0` :  
        (double): ordinate of the lower left corner of the rectangle  
    * `lx` :  
        (double): length of the mesh, following the axis (Ox)  
    * `ly` :  
        (double): length of the mesh, following the axis (Oy)  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int): size of coor  
    * `size_nb_node_per_ele` :  
        (int): size of nb_node_per_ele  
    * `size_conn` :  
        (int): size of conn  
    * `coor` :  
        (double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    * `nb_node_per_ele` :  
        (int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    * `conn` :  
        (int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_Mesh4T3 "

this function computes and return a mesh made of T3 --- obtained by splitting a
Q4 in four T3 --- in the following generic format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     2 4  *-----*  
     |\\ 4 /|  
     | \\ / |  
     |1 5 3|  
     | / \\ |  
     |/ 2 |  *-----*  
     1 3  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5]  
*   nb_node_per_ele = [3, 3, 3, 3]  
*   conn = [1, 5, 2, 1, 3, 5, 3, 4, 5, 2, 5, 4]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_Mesh4T3(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    * `x0` :  
        (double): abscissa of the lower left corner of the rectangle  
    * `y0` :  
        (double): ordinate of the lower left corner of the rectangle  
    * `lx` :  
        (double): length of the mesh, following the axis (Ox)  
    * `ly` :  
        (double): length of the mesh, following the axis (Oy)  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int): size of coor  
    * `size_nb_node_per_ele` :  
        (int): size of nb_node_per_ele  
    * `size_conn` :  
        (int): size of conn  
    * `coor` :  
        (double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    * `nb_node_per_ele` :  
        (int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    * `conn` :  
        (int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

%feature("docstring") mesh2D_MeshQ8 "

this function computes and returns a mesh made of Q8 in the following generic
format:  

*   coor: vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
*   nb_node_per_ele: nb_node_per_ele(i) contains the number of nodes for element
    i, i in [1, number of elements]  
*   conn: vector storing the connectivity of the elements [n11, n12n n13, n21,
    n22, n23, n24, ...] consider the following little mesh:  
     3 7 3  *---*---*  
     | |  
     | |  
     8 * 1 * 6  
     | |  
     | |  *---*---*  
     1 5 2  
    the vectors for this mesh read:  
*   coor = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8]  
*   nb_node_per_ele = [8]  
*   conn = [1, 2, 3, 4, 5, 6, 7, 8]  

    **Warning**: python call: [coor, nb_node_per_ele, conn]=mesh2D_MeshQ8(x0,
        y0, lx, ly,
           nb_elem_x, nb_elem_y, size_coor, size_nb_node_per_ele, size_conn)  

    Parameters:  
    * `x0` :  
        (double): abscissa of the lower left corner of the rectangle  
    * `y0` :  
        (double): ordinate of the lower left corner of the rectangle  
    * `lx` :  
        (double): length of the mesh, following the axis (Ox)  
    * `ly` :  
        (double): length of the mesh, following the axis (Oy)  
    * `nb_elem_x` :  
        (int): number of elements Q4 in the horizontal direction  
    * `nb_elem_y` :  
        (int): number of elements Q4 in the vertical direction  
    * `size_coor` :  
        (int): size of coor  
    * `size_nb_node_per_ele` :  
        (int): size of nb_node_per_ele  
    * `size_conn` :  
        (int): size of conn  
    * `coor` :  
        (double *): vector of coordinates of the nodes [x1, y1, x2, y2, ...]  
    * `nb_node_per_ele` :  
        (int *): nb_node_per_ele(i) contains the number of nodes for element i,
        i in [1, number of elements]  
    * `conn` :  
        (int *): vector storing the connectivity of the elements [n11, n12n n13,
        n21, n22, n23, n24, ...]  
";

