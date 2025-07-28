
// File: wrap__mp__solver_8h.xml

%feature("docstring") mp_solver_ReadMpBehaviour "

python usage : mp_solver_ReadMpBehaviour()  
";

%feature("docstring") mp_solver_WriteMpBehaviour "

python usage : mp_solver_WriteMpBehaviour()  
";

%feature("docstring") mp_solver_ReadIniMpValues "

Read MP_VALUES file.  

If num <= 0 : DATBOX/MP_VALUES.INI file is read Else : OUTBOX/MP_VALUES.OUT.num
is read, num being the parameter used in TimeEvolution_ReadIniDof last call  

usage : mp_solver_ReadIniMpValues(num=0)  

Parameters
----------
* `num` :  
    (integer) : which file to read  
";

%feature("docstring") mp_solver_WriteOutMpValues "

python usage : mp_solver_WriteOutMpValues()  
";

%feature("docstring") mp_solver_WriteLastMpValues "

python usage : mp_solver_WriteLastMpValues()  
";

%feature("docstring") mp_solver_SolveElectro1G "

python usage : mp_solver_SolveElectro1G()  
";

%feature("docstring") mp_solver_SolveNlElectro1G "

python usage : mp_solver_SolveNlElectro1G()  
";

%feature("docstring") mp_solver_SolveThermoProblem "

python usage : mp_solver_SolveThermoProblem()  
";

%feature("docstring") mp_solver_UpdateThermoProblem "

python usage : mp_solver_UpdateThermoProblem()  
";

%feature("docstring") mp_solver_RecupTemperature "

python usage : mp_solver_RecupTemperature()  
";

%feature("docstring") mp_solver_RecupPotential "

python usage : mp_solver_RecupPotential()  
";

%feature("docstring") mp_solver_UpdateConductivity "

python usage : mp_solver_UpdateConductivity()  
";

%feature("docstring") mp_solver_InitThermalConductivity "

python usage : mp_solver_InitThermalConductivity()  
";

%feature("docstring") mp_solver_GetBrancheValues "
";

%feature("docstring") mp_solver_PutHeatGenerationFactor "

python usage : value = mp_solver_PutHeatGenerationFactor(ivalue)  
";

%feature("docstring") mp_solver_PutHeatConductionContinueFactor "

python usage : value = mp_solver_PutHeatConductionContinueFactor(ivalue)  
";

