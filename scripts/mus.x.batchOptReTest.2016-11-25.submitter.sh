#!/bin/bash

SCRATCH=/mnt/fls01-home01/mqbssdgb/scratch

SCRIPT_TO_RUN=$SCRATCH/predictfromsequence/scripts/mus.x-inactivationModel.reTestOpt.parallel.fromConfig.R

CONFIG_FILE=$SCRATCH/predictfromsequence/scripts/CONFIG.mus.x-inactivationModel.BatchOpt.2016-11-25.dpsf.R

# The following line specifies that there are 100 jobs
#$ -t 1-100

#$ -pe smp.pe 16     

# each job to use 16 cores !

# Each time a job is submitted the environment variable SGE_TASK_ID is incremented so to get the nth line of

# The command line to run
Rscript $SCRIPT_TO_RUN --config $CONFIG_FILE --job_id $SGE_TASK_ID


# END OF SCRIPT
