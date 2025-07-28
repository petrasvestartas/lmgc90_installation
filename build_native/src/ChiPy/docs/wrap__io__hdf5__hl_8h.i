
// File: wrap__io__hdf5__hl_8h.xml

%feature("docstring") io_hdf5_initOutFile "

Init HDF5 file in which to write results.  

python usage : io_hdf5_initOutFile(filename)  

Parameters
----------
* `filename` :  
    (string) : file in which to write  
";

%feature("docstring") io_hdf5_write "

write output data in HDF5 file (GPV, DOF and VlocRloc).  

python usage : io_hdf5_write()  
";

%feature("docstring") io_hdf5_write_last "

write output data in HDF5 file (GPV, DOF and VlocRloc) in file  

python usage : io_hdf5_write_last(filename)  

Parameters
----------
* `filename` :  
    (string) : file in which to write  
";

%feature("docstring") io_hdf5_read "

read output data from HDF5 file (DOF and VlocRloc).  

python usage : io_hdf5_read(filename, step)  

Parameters
----------
* `filename` :  
    (string) : file to read  
* `step` :  
    (integer) : step number to read  
";

%feature("docstring") io_hdf5_cleanMemory "

cleanMemory of io_hdf5 module  

python usage : io_hdf5_cleanMemory()  
";

%feature("docstring") io_hdf5_fixVersion "

Will try to fix the file when reading it.  

Because the parameters changed within version 0, this flag is needed to fix the
file whend reading it.  

python usage : io_hdf5_fixVersion(version)  
";

