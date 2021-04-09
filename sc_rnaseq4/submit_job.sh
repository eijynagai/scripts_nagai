#!/bin/sh
#$ -S /bin/sh
#$ -o /home/park/Hayashi/Seurat3/SGE
#$ -e /home/park/Hayashi/Seurat3/SGE
#$ -v LD_LIBRARY_PATH=""
#$ -M park@hgc.jp
#$ -m a
#$ -N mNCC
#$ -cwd -j y -l os7 -l ljob -l s_vmem=30G -l mem_req=30G

cd /home/park/Hayashi/Seurat3

module use /usr/local/package/modulefiles
module load R/3.6

PROG1=/home/park/Hayashi/Seurat3/script/Park_Run_Seurat3_Monocle3.pl
PROG2=/home/park/Hayashi/Seurat3/script/Park_takeMarkers.pl

perl ${PROG1} /home/park/Hayashi/CellRanger/DrHayashi_No174_Sample1_cellranger/cellranger_outputs/filtered_feature_bc_matrix/ $1 $2 $3 $4

perl ${PROG2} $1 $1"/ConsensusMarkers/"



##echo "START: " `date  "+%Y%m%d-%H%M%S"`
echo "DONE: " `date  "+%Y%m%d-%H%M%S"`
#start_time=`date +%s`
end_time=`date +%s`
duration=$((end_time - start_time))
echo "TIME: " $duration


exit;
