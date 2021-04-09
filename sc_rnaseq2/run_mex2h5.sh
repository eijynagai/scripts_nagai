#!/bin/sh
#$ -S /bin/sh
#$ -cwd                       #job is located in the current wd
#$ -j y                       #merges the contents of standard error with ouput
##$ -hold_jid <job_id>        #holds the jov until job_id ends
#$ -R y                       #reserve resources until have enough to start
#$ -v LD_LIBRARY_PATH=""      #export these environmental variables
#$ -V                         #necessary to run conda modules
##$ -masterl <resource>       #use this or l
#$ -pe def_slot 1                #parallel env name, use "def_slot 4"
#$ -M nagailae@hgc.jp
#$ -m bae                     #in what circunstance beg/abort/end/suspension

## modify below
#$ -l os7                     #run only os7
#$ -l s_vmem=16G               #virtual memory required
#$ -l mem_req=16G              #expected memory required; put same as s_vem
#$ -l mjob                   #opt: ljob,lmem,mjob; check using qfree
#$ -e /home/nagailae/projects/sc_celegans/seurat3/
#$ -o /home/nagailae/projects/sc_celegans/seurat3/
#$ -N mex2h5                  #include name in the; check qstat


#reset timer
SECONDS=0

# load conda environment
module load python/3.6
source activate science

    #working directory
    cd /home/nagailae/projects/sc_celegans/seurat3/

    #your script here
    R CMD BATCH script/mex2h5.R

    ### Duration information
    echo "START: " `date  "+%Y.%m.%d-%H:%M:%S"`
    echo "DONE: " `date  "+%Y.%m.%d-%H:%M:%S"`
    echo "DURATION: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

# deactivate the environment
conda deactivate

exit;