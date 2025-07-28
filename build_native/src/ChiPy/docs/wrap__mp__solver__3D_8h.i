
// File: wrap__mp__solver__3D_8h.xml

%feature("docstring") mp_solver_3D_ReadMpBehaviour "

python usage : mp_solver_3D_ReadMpBehaviour()  
";

%feature("docstring") mp_solver_3D_WriteMpBehaviour "

python usage : mp_solver_3D_WriteMpBehaviour()  
";

%feature("docstring") mp_solver_3D_ReadIniMpValues "

Read MP_VALUES file.  

If num <= 0 : DATBOX/MP_VALUES.INI file is read Else : OUTBOX/MP_VALUES.OUT.num
is read, num being the parameter used in TimeEvolution_ReadIniDof last call  

usage : mp_solver_3D_ReadIniMpValues(num=0)  

Parameters
----------
* `num` :  
    (integer) : which file to read  
";

%feature("docstring") mp_solver_3D_WriteOutMpValues "

python usage : mp_solver_3D_WriteOutMpValues()  
";

%feature("docstring") mp_solver_3D_WriteLastMpValues "

python usage : mp_solver_3D_WriteLastMpValues()  
";

%feature("docstring") mp_solver_3D_SolveElectro1G "

python usage : mp_solver_3D_SolveElectro1G()  
";

%feature("docstring") mp_solver_3D_SolveNlElectro1G "

python usage : mp_solver_3D_SolveNlElectro1G()  
";

%feature("docstring") mp_solver_3D_SolveThermoProblem "

python usage : mp_solver_3D_SolveThermoProblem()  
";

%feature("docstring") mp_solver_3D_UpdateThermoProblem "

python usage : mp_solver_3D_UpdateThermoProblem()  
";

%feature("docstring") mp_solver_3D_RecupTemperature "

python usage : mp_solver_3D_RecupTemperature()  
";

%feature("docstring") mp_solver_3D_RecupPotential "

python usage : mp_solver_3D_RecupPotential()  
";

%feature("docstring") mp_solver_3D_UpdateConductivity "

python usage : mp_solver_3D_UpdateConductivity()  
";

