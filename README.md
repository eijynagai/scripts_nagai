# My bioinformatics scripts

Bioinformatics scrips for Next Generation Sequencing data. Collections of useful scripts collected and created by me. Feel free to copy, modidy or contribute.

## Note 

The scripts here are writen to run in SHIROKANE supercomputer system at Human Genome Center, University of Tokyo.

Please follow the steps below to run the scripts under the Univa Grid Engine (UGE).

## How to use supercomputer SHIROKANE

To submit any command:
```qsub [ShellScript]```


To confirm the execution of your submission:
```qstat```

To delete any job in the queue or running:
```
qdel [jobID]
qdel -u [userName] #will delete all jobs of the user
```

Refer to the file `run_qsub.sh` for the most common parameters used in our system.
