
// File: wrap__timer_8h.xml

%feature("docstring") timer_InitializeTimers "

Set all timers to 0.  

python usage : timer_InitializeTimers()  
";

%feature("docstring") timer_WriteOutTimers "

write the cumulated times of all the timers  

python usage : timer_WriteOutTimers()  
";

%feature("docstring") timer_GetNewTimer "

create a new timer  

python usage : id = timer_GetNewTimer(name)  

Parameters
----------
* `name` :  
    (string) : name of new timer  

Returns
-------
id (integer) : id of the timer created  
";

%feature("docstring") timer_StartTimer "

start a given timer  

python usage : timer_StartTimer(timer_id)  

Parameters
----------
* `timer_id` :  
    (integer) : id of the timer to start  
";

%feature("docstring") timer_StopTimer "

stop a given timer, and add the elapsed time since start to the time  

python usage : timer_StopTimer(timer_id)  

Parameters
----------
* `timer_id` :  
    (integer) : id of the timer to stop  
";

%feature("docstring") timer_ClearAll "

clear all timers (internal, external and user)  

python usage : timer_ClearAll()  
";

