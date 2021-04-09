#!/usr/bin/sh

#$ -S /bin/sh
#$ -cwd                       #job is located in the current wd
#$ -j y                       #merges the contents of standard error with ouput
##$ -hold_jid <job_id>        #holds the jov until job_id ends
#$ -R y                       #reserve resources until have enough to start
#$ -v LD_LIBRARY_PATH=""      #export these environmental variables
#$ -V                         #necessary to run conda modules
##$ -masterl <resource>       #use this or l
#$ -l os7                     #run only os7


## modify below

#$ -pe def_slot 8  #parallel env name, use "def_slot 4"
#$ -M nagailae@hgc.jp
#$ -m bae                     #in what circunstance beg/abort/end/suspension
#$ -l s_vmem=90G               #virtual memory required
#$ -l mem_req=90G              #expected memory required; put same as s_vem
#$ -l lmem                 #opt: ljob,lmem,mjob; check using qfree
#$ -e /home/nagailae/projects/sc_ciona/seurat3_v1/SGE
#$ -o /home/nagailae/projects/sc_ciona/seurat3_v1/SGE
#$ -N ciona                 #include name in the; check qstat


#reset timer
SECONDS=0

# load conda environment
#module load python3
module load python/3.6
source activate science

#working directory
cd /home/nagailae/projects/sc_ciona/seurat3_v1/

prog1=/home/nagailae/projects/sc_ciona/seurat3_v1/script/initializationScript.py
prog2=/home/nagailae/projects/sc_ciona/seurat3_v1/script/test_s1.R

#not implemented yet
#prog2=/home/nagailae/projects/sc_ciona/seurat3_v1/script/step1_qc.R
#prog3=/home/nagailae/projects/sc_ciona/seurat3_v1/script/step2_clustering.R

indir=/home/nagailae/projects/sc_ciona/seurat3_v1/data
outdir=/home/nagailae/projects/sc_ciona/seurat3_v1/result
proj=ciona
res=$1

#your script here
python $prog1 $indir $outdir $proj $res
Rscript ${prog2} --outdir ${outdir}/${proj}_res${res} --resolution ${res} >& ${outdir}/${proj}_res${res}/R.out

#not implemented yet
#Rscript ${prog3} --outdir ${outdir}/${proj}_res${res} --resolution ${res} >& ${outdir}/${proj}_res${res}/R2.out


### Duration information
echo "START: " `date  "+%Y.%m.%d-%H:%M:%S"`
echo "DONE: " `date  "+%Y.%m.%d-%H:%M:%S"`
echo "DURATION: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

# deactivate the environment
conda deactivate

exit;
