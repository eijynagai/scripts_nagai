#!/bin/env sh
#$ -S /bin/sh
#$ -cwd                       #job is located in the current wd
#$ -jy                        #merges the contents of standard error with ouput
##$ -hold_jid <job_id>        #holds the jov until job_id ends
#$ -R y                       #reserve resources until have enough to start
#$ -v LD_LIBRARY_PATH="",PATH=$PATH      #export these environmental variables 
##$ -masterl <resource>       #use this or l
##$ -pe <name>                #parallel env name, use "def_slot 4"
##$ -M <your_email>
##$ -m bae                     #in what circunstance beg/abort/end/suspension

## modify below
#$ -l os7                     #run only os7
#$ -l s_vmem=8G               #virtual memory required
#$ -l mem_req=8G              #expected memory required; put same as s_vem
#$ -l <resource>              #opt: ljob,lmem,mjob; check using qfree
#$ -e <path>                  #saves the execution result, std error
#$ -o <path>                  #saves the execution result, std out
#$ -N <job name>              #include name in the; check qstat

#cd <dir>                     #working directory

echo "Hello world"


### Duration information
##echo "START: " `date  "+%Y%m%d-%H%M%S"`
echo "DONE: " `date  "+%Y%m%d-%H%M%S"`
#start_time=`date +%s`
end_time=`date +%s`
duration=$((end_time - start_time))
echo "TIME: " $duration


exit;
