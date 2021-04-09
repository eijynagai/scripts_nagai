#!/bin/sh
#$ -S /bin/sh
#$ -o /home/park/Hayashi/Velocyto/SGE
#$ -e /home/park/Hayashi/Velocyto/SGE
#$ -v LD_LIBRARY_PATH=""
#$ -M park@hgc.jp
#$ -m a
#$ -N subCLSVeloc
#$ -cwd -j y -l ljob -l s_vmem=30G -l mem_req=30G

cd /home/park/Hayashi/Velocyto/

module use /usr/local/package/modulefiles
module load R/3.6

##S3_DIR=/home/park/Hayashi/Seurat3/Mouse/Analysis/H-T-_Res1.1/
##DrpEst=/home/park/Hayashi/Velocyto/DropEst/cell.matrices.rds
##OUT_DIR=Velo_Res1.1/

##perl script/Park_RunVelocyto.pl ${S3_DIR} ${DrpEst} ${OUT_DIR} ${FRESH}
perl script/subCLS_Park_RunVelocyto.pl $1 $2 $3 $4


##echo "START: " `date  "+%Y%m%d-%H%M%S"`
echo "DONE: " `date  "+%Y%m%d-%H%M%S"`
#start_time=`date +%s`
end_time=`date +%s`
duration=$((end_time - start_time))
echo "TIME: " $duration


exit;
