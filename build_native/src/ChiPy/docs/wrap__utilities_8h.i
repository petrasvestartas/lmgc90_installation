
// File: wrap__utilities_8h.xml

%feature("docstring") utilities_logMes "

ask to write a message  

If the message is too long (more than 256 characters) it will be truncated  

python usage : utilities_logMes(message)  

Parameters
----------
* `message` :  
    (string) : log message to add  
* `length` :  
    (integer) : length of the message string  
";

%feature("docstring") utilities_DisableLogMes "

disable printing of messages  

python usage : utilities_DisableLogMes()  
";

%feature("docstring") utilities_EnableLogMes "

enable priting of messages  

python usage : utilities_EnableLogMes()  
";

%feature("docstring") utilities_setIoUnitLimits "

set the interval of unit numbers lmgc90 can use to open file  
";

%feature("docstring") utilities_setStopMode "

Decide to stop or store a message in case of fatal error.  

python usage : utilities_setStopMode()  
";

%feature("docstring") utilities_resetFatal "

Clean fatal error state.  

This function is not intended to be used in python but by swig to throw an
excpetion  
";

%feature("docstring") utilities_checkFatal "
";

%feature("docstring") utilities_OpenFileStandardOutput "

Select the file for standard and errors outputs.  

If the filename is too long (more than 256 characters) it will be truncated  

python usage : utilities_OpenFileStandardOutput(filename)  

Parameters
----------
* `filename` :  
    (string) : the name of file  
* `length` :  
    (integer) : length the name of the file  
";

%feature("docstring") utilities_CloseFileStandardOutput "

Close the file for standard and errors outputs.  

python usage : utilities_CloseFileStandardOutput()  
";

%feature("docstring") utilities_InitRandomSeed "

Re-initialize the seed of the build-in random function.  

python usage : utilities_InitRandomSeed([seed])  

Parameters
----------
* `seed` :  
    (integer array) : an optional desired input seed  
";

%feature("docstring") utilities_Finalize "

End of simulation operations.  

Only close all possibly opened units by the program.  

python usage : utilities_Finalize()  
";

