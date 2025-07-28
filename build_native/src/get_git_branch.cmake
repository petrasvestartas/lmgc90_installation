
execute_process( COMMAND /usr/bin/git rev-parse --abbrev-ref HEAD
                 WORKING_DIRECTORY /Users/petras/brg/2_code/lmgc90_installation
                 OUTPUT_VARIABLE GIT_CURRENT_BRANCH
                 OUTPUT_STRIP_TRAILING_WHITESPACE
               )
message(STATUS "git branch main | ${GIT_CURRENT_BRANCH}")
if( NOT "main" STREQUAL "${GIT_CURRENT_BRANCH}" )
  message(FATAL_ERROR "git branch changed... please rerun cmake to confirm")
endif( NOT "main" STREQUAL "${GIT_CURRENT_BRANCH}" )

