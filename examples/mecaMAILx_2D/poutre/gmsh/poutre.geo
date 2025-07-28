/*************************************************************/
/* D�finition de la g�om�trie pour le cas de la barre Taylor */
/*************************************************************/

/* d�finition des variables utilis�es */

// on d�finit le milim�tre
mm = 1e-3 ;

// on d�finit le nombre d'�l�ment dans la direction (0x)
nb_ele_x = 12 ;

// on d�finit le nombre d'�l�ment dans la direction (0y)
nb_ele_y = 4 ;

// on d�finit une longueur caract�risatique pour les points, inutile, puisqu'on
// veut un maillage r�gulier, mais malgr� tout indispensable
lc = 0.001 ;

/* d�finition de la g�om�trie */

// on d�finit les quatre coins du rectangle :
Point(1) = {0, 0, 0, lc} ;
Point(2) = {3.*mm, 0, 0, lc} ;
Point(3) = {3.*mm, 1.*mm, 0, lc} ;
Point(4) = {0, 1.*mm, 0, lc} ;

// on d�finit les lignes constituant sa fronti�re :
Line(1) = {1, 2} ;
Line(2) = {2, 3} ;
Line(3) = {3, 4} ;
Line(4) = {4, 1} ;

// on peut alors d�finir la fronti�re, correctement orient�e, du rectangle :
Line Loop(5) = {1, 2, 3, 4} ;

// et finalement, la surface constitu�e par le rectangle :
Plane Surface(6) = {5} ;

/* d�finition du maillage */

// on place (nb_ele_x + 1) points �quidistants sur les lignes 1 et 3 :
Transfinite Line{1, 3} = nb_ele_x + 1 ;

// on place (nb_ele_y + 1) points �quidistants sur les lignes 2 et 4 :
Transfinite Line{2, 4} = nb_ele_y + 1 ;

// on sp�cifie les quatre points de l'interploation "transfinie", pour d�finir
// le rectangle comme une surface transfinie :
Transfinite Surface{6} = {1, 2, 3, 4};

// on r�organise les traingles en quadrangles :
Recombine Surface{6};

// definition des groupes pour appliquer les condtions limites :
//    * point sur lequel on tire
Physical Point(8) = {2};
//    * point oppose a celui lequel on tire
Physical Point(9) = {3};
//    * ligne encastree
Physical Line(7) = {4};
//    * affecation d'un numero bidon, pour garder les elements de surface
Physical Surface(20) = {6};

