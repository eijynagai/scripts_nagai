#!/bin/sh
#$ -S /bin/sh
#$ -o /home/park/Hayashi/Velocyto/SGE
#$ -e /home/park/Hayashi/Velocyto/SGE
#$ -v LD_LIBRARY_PATH=""
#$ -M park@hgc.jp
#$ -m a
#$ -N DrpEst
#$ -cwd -j y -l ljob -l s_vmem=55G -l mem_req=55G

cd /home/park/Hayashi/Velocyto

DIR="DropEst"
PROG=/home/park/TOOL/dropEst/bin/dropest
### this input bam is filtered Cellranger BAM file; only chr1...22.X,Y and mapped read only
### look at the file in dir
INBAM=/home/park/Hayashi/CellRanger/DrHayashi_No174_Sample1_cellranger/BAM/CR_filt.bam
REF=../REF/mm10_genes.gtf
CONF=../script/this_10x.xml 

if [ ! -e ${INBAM} ]; then
    echo ${INBAM} " not found"
    exit
fi

if [ ! -e ${DIR} ]; then
    mkdir ${DIR}
fi

cd ${DIR}

${PROG} -f -V -g ${REF} -c ${CONF} ${INBAM}


##echo "START: " `date  "+%Y%m%d-%H%M%S"`
echo "DONE: " `date  "+%Y%m%d-%H%M%S"`
#start_time=`date +%s`
end_time=`date +%s`
duration=$((end_time - start_time))
echo "TIME: " $duration


exit
