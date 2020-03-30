#!/bin/env sh
#$ -l <resource>
#$ -cwd                       #job is located in the current wd
#$ -e <path>                  #saves the execution result, std error
#$ -jy                        #merges the contents of standard error with ouput
#$ -N <job name> 
#$ -o <path>                  #saves the execution result, std out
#$ -S /bin/sh
##$ -pe <name>                #parallel env name, use "def_slot 4"
##$ -masterl <resource>       #use this or l
##$ -hold_jid <job_id>        #holds the jov until job_id ends
#$ -R y                       #reserve resources until have enough to start
#$ -M <email contact>
#$ -m bae                     #in what circunstance beg/abort/end/suspension
#$ -v LD_LIBRARY_PATH="",PATH=$PATH      #export these environmental variables 

echo "Hello world"

exit;
