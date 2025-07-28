
// File: wrap__models_8h.xml

%feature("docstring") models_ReadModels "

read models from DATBOX/MODELS.DAT  

python usage : models_ReadModels()  
";

%feature("docstring") models_WriteModels "

write models to OUTBOX/MODELS.OUT  

python usage : models_WriteModels()  
";

%feature("docstring") models_InitModels "

initialize models  

python usage : models_InitModels()  
";

%feature("docstring") models_InitProperties "

initialize properties  

In face re-initialize properties (since it is done in InitModels). Necessary if
a Store has been done and it is wanted again to LoadModel  

python usage : models_InitProperties()  
";

%feature("docstring") models_StoreProperties "

create properties (couple of model and models)  

python usage : models_StoreProperties()  
";

%feature("docstring") models_CleanMemory "

Free all memory allocated within models module.  

python usage : models_CleanMemory()  
";

