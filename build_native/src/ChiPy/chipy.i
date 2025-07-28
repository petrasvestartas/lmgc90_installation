/* declarations : en-tete d'un fichier SWIG + ajout des headers necessaires */
%{
#define SWIG_FILE_WITH_INIT

#if PY_MAJOR_VERSION >= 3
#define PyString_FromStringAndSize PyUnicode_FromStringAndSize
#endif

#include "rigid_2D/wrap_DISKx.h"
#include "rigid_2D/wrap_JONCx.h"
#include "rigid_2D/wrap_POLYG.h"
#include "rigid_2D/wrap_PT2Dx.h"
#include "rigid_2D/wrap_RBDY2.h"
#include "rigid_2D/wrap_xKSID.h"

#include "rigid_3D/wrap_CYLND.h"
#include "rigid_3D/wrap_DNLYC.h"
#include "rigid_3D/wrap_PLANx.h"
#include "rigid_3D/wrap_POLYR.h"
#include "rigid_3D/wrap_PT3Dx.h"
#include "rigid_3D/wrap_RBDY3.h"
#include "rigid_3D/wrap_SPHER.h"

#include "mbs/wrap_mbs2D.h"
#include "mbs/wrap_mbs3D.h"

#include "mailx/wrap_ALpxx.h"
#include "mailx/wrap_ASpxx.h"
#include "mailx/wrap_CLxxx.h"
#include "mailx/wrap_CSxxx.h"
#include "mailx/wrap_DISKL.h"
#include "mailx/wrap_MAILx.h"
#include "mailx/wrap_PT2DL.h"
#include "mailx/wrap_mecaMAILx.h"
#include "mailx/wrap_therMAILx.h"
#include "mailx/wrap_poroMAILx.h"
#include "mailx/wrap_multiMAILx.h"

#include "shared/wrap_ExternalModels.h"
#include "shared/wrap_bulk_behav.h"
#include "shared/wrap_models.h"
#include "shared/wrap_overall.h"
#include "shared/wrap_tact_behav.h"
#include "shared/wrap_timer.h"
#include "shared/wrap_utilities.h"
#include "shared/wrap_a_EF.h"
#include "shared/wrap_meca_polygon.h"
#include "shared/wrap_parameters.h"

#include "contact_2D/wrap_CLALp.h"
#include "contact_2D/wrap_CLJCx.h"
#include "contact_2D/wrap_DKALp.h"
#include "contact_2D/wrap_DKDKL.h"
#include "contact_2D/wrap_DKDKx.h"
#include "contact_2D/wrap_DKJCx.h"
#include "contact_2D/wrap_DKKDx.h"
#include "contact_2D/wrap_DKPLx.h"
#include "contact_2D/wrap_P2P2L.h"
#include "contact_2D/wrap_PLALp.h"
#include "contact_2D/wrap_PLJCx.h"
#include "contact_2D/wrap_PLPLx.h"
#include "contact_2D/wrap_PTPT2.h"
#include "contact_2D/wrap_inter_handler_2D.h"

#include "contact_3D/wrap_CDCDx.h"
#include "contact_3D/wrap_CDPLx.h"
#include "contact_3D/wrap_CSASp.h"
#include "contact_3D/wrap_CSPRx.h"
#include "contact_3D/wrap_PRASp.h"
#include "contact_3D/wrap_PRPLx.h"
#include "contact_3D/wrap_PRPRx.h"
#include "contact_3D/wrap_PTPT3.h"
#include "contact_3D/wrap_SPCDx.h"
#include "contact_3D/wrap_SPDCx.h"
#include "contact_3D/wrap_SPPLx.h"
#include "contact_3D/wrap_SPPRx.h"  
#include "contact_3D/wrap_SPSPx.h"
#include "contact_3D/wrap_inter_handler_3D.h"

#include "kernel/wrap_cpg.h"
#include "kernel/wrap_cpg_3D.h"
#include "kernel/wrap_mp_solver.h"
#include "kernel/wrap_mp_solver_3D.h"
#include "kernel/wrap_nlgs.h"
#include "kernel/wrap_nlgs_3D.h"

#include "kernel/wrap_global_thermal_solver.h"

#ifdef WITH_SICONOS
#include "kernel/wrap_SiconosNumerics.h"
#endif

#ifdef WITH_MPI
#include "kernel/wrap_DDM_2D.h"
#include "kernel/wrap_DDM_3D.h"
#include "kernel/wrap_DDM_ExternalFEM.h"
#endif

#include "post/wrap_postpro.h"
#include "post/wrap_postpro_3D.h"

#include "user/wrap_user.h"

#include "pre_tools/wrap_cut2D.h"
#include "pre_tools/wrap_deposit2D.h"
#include "pre_tools/wrap_deposit3D.h"
#include "pre_tools/wrap_mesh2D.h"
#include "pre_tools/wrap_surface_T3.h"

// HDF5 include
#ifdef WITH_HDF5
#include "io/wrap_io_hdf5_hl.h"
#endif

%}

/* directives supplementaires pour assurer le transfert d'argumlents/resultats
 * entre SWIG et numpy */
%include "/Users/petras/brg/2_code/lmgc90_installation/src/tools/swig/numpy-2.2.6.i"
%init %{
import_array();
%}

/* directive supplementaire pour effectuer des passages d'arguments par adresse
 * ou transformer des procedures en fonctions */
%include "typemaps.i"
%include "cstring.i"
%include "exception.i"

/* declaration du nom du module */
%module lmgc90

/* exception handling */
%{
#include "/Users/petras/brg/2_code/lmgc90_installation/src/Core/contribs/exception.h"
%}

%exception
{
  char * mess;
  int err;

  utilities_resetFatal();
  if( (err = setjmp(except_buf)) == 0 )
  {
    $function
  }
  else
  {
    err = utilities_checkFatal(&mess);
    SWIG_exception(SWIG_RuntimeError,mess);
  }
}


/* for default args working when generating docstrings...
 (use swig prior 1.3.23 version management of default arguments
 WARNING : may not be compatible with overloading (if we want to use it one day)
*/
%feature("compactdefaultargs");

/* definitions des regles a appliquer */

/* pour la fonction get_scalar d un double : on renomme la variable "res" en "OUTPUT" pour
 * que SWIG transforme la procedure (FORTRAN) get_scalar, qui place dans res
 * le scalaire stocke dans le module, en une fonction (python) qui le renvoie */
%apply double *OUTPUT { double *res };

/* pour la fonction get_scalar d un int: on renomme la variable "res" en "OUTPUT" pour
 * que SWIG transforme la procedure (FORTRAN) get_scalar, qui place dans res
 * le scalaire stocke dans le module, en une fonction (python) qui le renvoie */
%apply int *OUTPUT { int *ires };
%apply int *OUTPUT { int *ires2 };
%apply int *OUTPUT { int *ires3 };

/* turn an input numpy array in a pointer and a size usable in C/Fortran
 */
%apply (double* IN_ARRAY1, int DIM1) {(double * rvector_in, int rlength_in)};
%apply (double* IN_ARRAY1, int DIM1) {(double * rvector_in2, int rlength_in2)};
%apply (double* IN_ARRAY1, int DIM1) {(double * rvector_in3, int rlength_in3)};
%apply (int* IN_ARRAY1, int DIM1) {(int * ivector_in, int ilength_in)};
%apply (int* IN_ARRAY1, int DIM1) {(int * ivector_in2, int ilength_in2)};

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix_in, int dim1 , int dim2)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix_in, int idim1, int idim2)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* rmatrix_in, int rdim1, int rdim2)};

/* length is an input, the size of the vector to get back as an output in python
 */
%apply (double *ARGOUT_ARRAY1, int DIM1) {(double * rvector_out, int rlength_out)};
%apply (double *ARGOUT_ARRAY1, int DIM1) {(double * rvector_out2, int rlength_out2)};
%apply (int *ARGOUT_ARRAY1, int DIM1) {(int * ivector_out, int ilength_out)};
%apply (int *ARGOUT_ARRAY1, int DIM1) {(int * ivector_out2, int ilength_out2)};


/* get a view on some C/Fortran array as numpy vector (get a pointer on memory)
 */
%apply (double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double** pointer_out, int* length)};
%apply (double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double** pointer_out2, int* length2)};
%apply (int** ARGOUTVIEW_ARRAY1, int* DIM1) {(int** pointer_out, int* length)};
%apply (int** ARGOUTVIEW_ARRAY1, int* DIM1) {(int** pointer_out2, int* length2)};

%apply (double** ARGOUT_ALMOSTVIEW_ARRAY1, int* DIM1) {(double** r8_vector, int* r8_size)};
%apply (int** ARGOUT_ALMOSTVIEW_ARRAY1, int* DIM1) {(int** i4_vector, int* i4_size)};
%apply (int** ARGOUT_ALMOSTVIEW_ARRAY1, int* DIM1) {(int** i4_vector_2, int* i4_size_2)};

%apply (int** ARGOUTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(int** pointer_out, int* dim1, int* dim2)};
%apply (double** ARGOUTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(double** pointer_out, int* dim1, int* dim2)};
%apply (int** ARGOUT_ALMOSTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(int** matrix_out, int* dim1, int* dim2)};
%apply (int** ARGOUT_ALMOSTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(int** i4_matrix, int* i4_dim1, int* i4_dim2)};
%apply (double** ARGOUT_ALMOSTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(double** matrix_out, int* dim1, int* dim2)};
%apply (double** ARGOUT_ALMOSTVIEW_FISOC_ARRAY2, int* DIM1, int* DIM2) {(double** matrix_out_2, int* dim1_2, int* dim2_2)};


/* on indique a SWIG que la fonction python prend un array numpy, qu'il faudra 
 * decomposer en un pointeur de double et un entier pour la passer a la 
 * procedure FORTRAN qui modifiera les valeurs */
%apply (double *INPLACE_ARRAY1, int DIM1) {(double *vector_inout, int length_inout)};
%apply (double *INPLACE_ARRAY1, int DIM1) {(double *vector_inout2, int length_inout2)};
/* idem, amsi avec des vecteurs d'entiers */
%apply (int *INPLACE_ARRAY1, int DIM1) {(int *ivector_inout, int ilength_inout)};


/* pour jouer avec les chaines de caracteres :
   copie du vecteur de chaine de caracteres de taille 5 (la chaine pas le vecteur)
   dans un tuple... pourquoi dans un tuple ? parce que dans une liste je n'ai pas reussi ><
 * */
%typemap(in,numinputs=0)
        (char** c5_vector, int* c5_size)
        (char* tmp_vector, int tmp_size)
{
  $1 = &tmp_vector;
  $2 = &tmp_size;
}
%typemap(argout) (char** c5_vector, int* c5_size)
{
  int i;
  PyObject *liste, *chaine;
  liste = PyTuple_New(*$2);
  for( i=0; i<*$2; i++ )
  {
    chaine = PyString_FromStringAndSize((*$1+i*5), 5);
    PyTuple_SetItem(liste, i, chaine);
  }

  $result = SWIG_AppendOutput($result, liste);
}

%typemap(in,numinputs=0)
        (char** string_vector, int* vector_size, int* string_size)
        (char* tmp_vector, int tmp_vsize, int tmp_csize)
{
  $1 = &tmp_vector;
  $2 = &tmp_vsize;
  $3 = &tmp_csize;
}
%typemap(argout) (char** string_vector, int* vector_size, int* string_size)
{
  int i;
  PyObject *liste, *chaine;
  liste = PyTuple_New(*$2);
  for( i=0; i<*$2; i++ )
  {
    chaine = PyString_FromStringAndSize((*$1+i*(*$3)), *$3);
    PyTuple_SetItem(liste, i, chaine);
  }

  $result = SWIG_AppendOutput($result, liste);
}

%cstring_bounded_output(char* outtype, 5);

/* ca c'est pour avoir en entree une liste de string
   dont chaque string est de taille 5 */
%typemap(in) (char** cvector_in, int clength) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size)*sizeof(char *));
    $2 = size;
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o)) {
        if( PyString_Size(PyList_GetItem($input,i)) > 6 ) {
          PyErr_SetString(PyExc_ValueError,"strings in list must be at the most of 5 characters");
          free($1);
          return NULL;
        }
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      }
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($1);
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
/* This cleans up the char ** array we malloc'd before the function call
*/
%typemap(freearg) (char ** cvector_in, int clength) {
  free((char *) $1);
}

/* ca c'est pour avoir en entree une liste de string
   dont chaque string est de taille 128... comme au dessus */
%typemap(in) (char** svector_in, int slength) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size)*sizeof(char *));
    $2 = size;
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o)) {
        if( PyString_Size(PyList_GetItem($input,i)) > 129 ) {
          PyErr_SetString(PyExc_ValueError,"strings in list must be at the most of 128 characters");
          free($1);
          return NULL;
        }
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      }
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($1);
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

/* This cleans up the char ** array we malloc'd before the function call
*/
%header %{
  extern "C" void free_ptrFortranString(char * array, int length);
%}

%typemap(freearg) (char ** svector_in, int slength) {
  free((char *) $1);
}
/* to get as a returned value a string of size 5 returned as an output parameter
*/
%typemap(in,numinputs=0)
    (char** c5)
    (char* tmp_c5)
{
  $1 = &tmp_c5;
}
%typemap(argout) (char** c5)
{
  PyObject *chaine;
  chaine = PyString_FromStringAndSize(*$1, 5);
  free_ptrFortranString(*$1, 5);
  $result = SWIG_AppendOutput($result, chaine);
}

/* to get as a returned value a string of size returned as an output parameter
*/
%typemap(in,numinputs=0)
    (char** string_out, int* string_size, int* real_size)
    (char* tmp_vector, int tmp_size, int tmp_rs)
{
  $1 = &tmp_vector;
  $2 = &tmp_size;
  $3 = &tmp_rs;
}
%typemap(argout) (char** string_out, int* string_size, int* real_size)
{
  PyObject *chaine;
  chaine = PyString_FromStringAndSize(*$1, *$2);
  free_ptrFortranString(*$1, *$3);
  $result = chaine;
}

/* to add to its length to an input string */
%typemap(in) (char * ul_string, int string_size)
{
  if (PyUnicode_Check($input)){
    PyObject *pyStr = PyUnicode_AsEncodedString($input, "utf-8", "Error~");
    $1 = PyString_AsString(pyStr);
    $2 = PyString_Size(pyStr);
    Py_XDECREF(pyStr);
  }
  else if (PyString_Check($input))
  {
    $1 = PyString_AsString($input);
    $2 = PyString_Size($input);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,"input must be a string");
    return NULL;
  }
}


// Some function renaming to manage optional arguments
%rename(inter_handler_2D_tsetInternal) wrap_tsetInternal_2D;
%inline %{
void wrap_tsetInternal_2D(int inter_id, int icdan, int index, double value) {
   inter_handler_2D_tsetInternal(inter_id, icdan, NULL, 0, index, value);
}
%}

%rename(inter_handler_3D_tsetInternal) wrap_tsetInternal_3D;
%inline %{
void wrap_tsetInternal_3D(int inter_id, int icdan, int index, double value) {
   inter_handler_3D_tsetInternal(inter_id, icdan, NULL, 0, index, value);
}
%}


%include "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/ChiPy/docs/chipy_swig_docstrings.i"

/* on inclut directement les fonctions definies dans le fichier d'en-tete, pour
 * quelles soient wrappees dans le module python */
%include "rigid_2D/wrap_DISKx.h"
%include "rigid_2D/wrap_JONCx.h"
%include "rigid_2D/wrap_POLYG.h"
%include "rigid_2D/wrap_PT2Dx.h"
%include "rigid_2D/wrap_RBDY2.h"
%include "rigid_2D/wrap_xKSID.h"

%include "rigid_3D/wrap_CYLND.h"
%include "rigid_3D/wrap_DNLYC.h"
%include "rigid_3D/wrap_PLANx.h"
%include "rigid_3D/wrap_POLYR.h"
%include "rigid_3D/wrap_PT3Dx.h"
%include "rigid_3D/wrap_RBDY3.h"
%include "rigid_3D/wrap_SPHER.h"

%include "mbs/wrap_mbs2D.h"
%include "mbs/wrap_mbs3D.h"

%include "mailx/wrap_ALpxx.h"
%include "mailx/wrap_ASpxx.h"
%include "mailx/wrap_CLxxx.h"
%include "mailx/wrap_CSxxx.h"
%include "mailx/wrap_DISKL.h"
%include "mailx/wrap_MAILx.h"
%include "mailx/wrap_PT2DL.h"
%include "mailx/wrap_mecaMAILx.h"
%include "mailx/wrap_therMAILx.h"
%include "mailx/wrap_poroMAILx.h"
%include "mailx/wrap_multiMAILx.h"

%include "shared/wrap_ExternalModels.h"
%include "shared/wrap_bulk_behav.h"
%include "shared/wrap_models.h"
%include "shared/wrap_overall.h"
%include "shared/wrap_tact_behav.h"
%include "shared/wrap_timer.h"
%include "shared/wrap_utilities.h"
%include "shared/wrap_a_EF.h"
%include "shared/wrap_meca_polygon.h"
%include "shared/wrap_parameters.h"

%include "contact_2D/wrap_CLALp.h"
%include "contact_2D/wrap_CLJCx.h"
%include "contact_2D/wrap_DKALp.h"
%include "contact_2D/wrap_DKDKL.h"
%include "contact_2D/wrap_DKDKx.h"
%include "contact_2D/wrap_DKJCx.h"
%include "contact_2D/wrap_DKKDx.h"
%include "contact_2D/wrap_DKPLx.h"
%include "contact_2D/wrap_P2P2L.h"
%include "contact_2D/wrap_PLALp.h"
%include "contact_2D/wrap_PLJCx.h"
%include "contact_2D/wrap_PLPLx.h"
%include "contact_2D/wrap_PTPT2.h"
%include "contact_2D/wrap_inter_handler_2D.h"

%include "contact_3D/wrap_CDCDx.h"
%include "contact_3D/wrap_CDPLx.h"
%include "contact_3D/wrap_CSASp.h"
%include "contact_3D/wrap_CSPRx.h"
%include "contact_3D/wrap_PRASp.h"
%include "contact_3D/wrap_PRPLx.h"
%include "contact_3D/wrap_PRPRx.h"
%include "contact_3D/wrap_PTPT3.h"
%include "contact_3D/wrap_SPCDx.h"
%include "contact_3D/wrap_SPDCx.h"
%include "contact_3D/wrap_SPPLx.h"
%include "contact_3D/wrap_SPPRx.h"
%include "contact_3D/wrap_SPSPx.h"
%include "contact_3D/wrap_inter_handler_3D.h"

%include "kernel/wrap_cpg.h"
%include "kernel/wrap_cpg_3D.h"
%include "kernel/wrap_mp_solver.h"
%include "kernel/wrap_mp_solver_3D.h"
%include "kernel/wrap_nlgs.h"
%include "kernel/wrap_nlgs_3D.h"

%include "kernel/wrap_global_thermal_solver.h"

#ifdef WITH_SICONOS
%include "kernel/wrap_SiconosNumerics.h"
#endif

#ifdef WITH_MPI
%include "kernel/wrap_DDM_2D.h"
%include "kernel/wrap_DDM_3D.h"
%include "kernel/wrap_DDM_ExternalFEM.h"
#endif

%include "post/wrap_postpro.h"
%include "post/wrap_postpro_3D.h"

%include "user/wrap_user.h"

%include "pre_tools/wrap_cut2D.h"
%include "pre_tools/wrap_deposit2D.h"
%include "pre_tools/wrap_deposit3D.h"
%include "pre_tools/wrap_mesh2D.h"
%include "pre_tools/wrap_surface_T3.h"

// HDF5 include
#ifdef WITH_HDF5
%include "io/wrap_io_hdf5_hl.h"
#endif


