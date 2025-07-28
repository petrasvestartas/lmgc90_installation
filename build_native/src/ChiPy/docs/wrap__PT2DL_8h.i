
// File: wrap__PT2DL_8h.xml

%feature("docstring") PT2DL_LoadTactors "

Initialize existing_entities variable for PT2DL contactors.  

python usage : PT2DL_LoadTactors()  
";

%feature("docstring") PT2DL_PushPreconNodes "

python usage : PT2DL_PushPreconNodes()  
";

%feature("docstring") PT2DL_GetNbPT2DL "

Get the number of PT2DL.  

usage : nb_PT2DL = PT2DL_GetNbPT2DL()  

Parameters
----------
* `nb_PT2DL` :  
    (integer) : number of PT2DL in container  
";

%feature("docstring") PT2DL_GetNbPT2TL "

Get the number of PT2TL of a body.  

usage : nb_PT2DL = PT2DL_GetNbPT2TL(ibdyty)  

Parameters
----------
* `nb_PT2TL` :  
    (integer) : number of PT2TL in container  
";

%feature("docstring") PT2DL_ComputeConvectiveFlux "

python usage : PT2DL_ComputeConvectiveFlux()  
";

%feature("docstring") PT2DL_AssembThermKT "

python usage : PT2DL_AssembThermKT()  
";

%feature("docstring") PT2DL_AssembThermRHS "

python usage : PT2DL_AssembThermRHS()  
";

%feature("docstring") PT2DL_GetBody "

return corresponding body  

python usage : ibdy = PT2DL_GetBody(itacty)  
";

%feature("docstring") PT2DL_CleanMemory "

Free all memory allocated within PT2DL module.  

python usage : PT2DL_CleanMemory()  
";

