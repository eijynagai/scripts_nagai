# My bioinformatics scripts

Bioinformatics scrips for Next Generation Sequencing data. Collections of useful scripts collected and created by me. Feel free to copy, modidy or contribute.

## Note 

The scripts here are writen to run in SHIROKANE supercomputer system at Human Genome Center, University of Tokyo.

Please follow the steps below to run the scripts under the Univa Grid Engine (UGE).

## How to use supercomputer SHIROKANE

To submit any command:
```
qsub [ShellScript]
```

To check which service is available
```
qfree
```

To confirm the execution of your submission:
```
qstat
```

To delete any job in the queue or running:
```
qdel [jobID]
qdel -u [userName] #will delete all jobs of the user
```

## Parameters

```
#!/bin/env sh
#$ -S /bin/sh
#$ -cwd                       #job is located in the current wd
#$ -j y                        #merges the contents of standard error with ouput
#$ -l s_vmem=8G               #virtual memory required
#$ -l mem_req=8G              #expected memory required; put same as s_vem
#$ -e <path>                  #saves the execution result, std error
#$ -o <path>                  #saves the execution result, std out
```

Refer to the file `run_qsub.sh` for a example of bash file to submit. This file includes more details of useful parameters.
