#!/bin/sh
#$ -S /bin/sh
#$ -o /home/park/Hayashi/Seurat3/SGE
#$ -e /home/park/Hayashi/Seurat3/SGE
#$ -v LD_LIBRARY_PATH=""
#$ -M park@hgc.jp
#$ -m a
#$ -N subCLSmNCC
#$ -cwd -j y -l s_vmem=20G -l mem_req=20G

cd /home/park/Hayashi/Seurat3

module use /usr/local/package/modulefiles
module load R/3.6

PROG1=/home/park/Hayashi/Seurat3/script/subCLS_Park_Run_Seurat3_Monocle3.pl
PROG2=/home/park/Hayashi/Seurat3/script/Park_takeMarkers.pl

perl ${PROG1} $1 $2 $3

perl ${PROG2} $2 $2"/ConsensusMarkers/"



##echo "START: " `date  "+%Y%m%d-%H%M%S"`
echo "DONE: " `date  "+%Y%m%d-%H%M%S"`
#start_time=`date +%s`
end_time=`date +%s`
duration=$((end_time - start_time))
echo "TIME: " $duration


exit;
