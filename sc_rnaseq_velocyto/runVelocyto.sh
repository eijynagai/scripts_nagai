

### After creating filtered Cellranger BAM,
### Do this at first (cell.counts.matrices.rds will be created)
###qsub script/1.submit_DropEst.sh


PROG=script/2.submit_Velocity.sh
RDS=/home/park/Hayashi/Velocyto/DropEst/cell.counts.matrices.rds

OUT_DIR=/home/park/Hayashi/Velocyto/Mouse/Analysis/
if [ ! -e ${OUT_DIR} ]; then
    mkdir ${OUT_DIR}
fi

FRESH=1
FRESH=0
resol=(0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0)

S3_DIR=/home/park/Hayashi/Seurat3/Mouse/Analysis/
type="H-T-"

for Res in ${resol[@]}
do
  dir=${S3_DIR}${type}_Res${Res}
  this_out=${OUT_DIR}${type}_Res${Res}

  if [ -e ${dir} ]; then
    qsub ${PROG} ${dir} ${RDS} ${this_out} ${FRESH}
  fi
done

exit


