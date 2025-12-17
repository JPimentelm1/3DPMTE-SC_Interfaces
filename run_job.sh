#!/bin/bash

# 1. Capture the Job Name from the first argument (equivalent to %1)
JOB_NAME="$1"

# 2. Get the absolute path of the current directory (equivalent to %CD%)
#    $PWD is the standard environment variable for this in Linux.
SCRATCH_PATH="$PWD"

echo "--- [Worker Linux] Received job: $JOB_NAME ---"
echo "--- [Worker Linux] Working in: $SCRATCH_PATH ---"

# 3. Run the Abaqus Job
#    Note the differences:
#    - Forward slashes / instead of backslashes \
#    - We use "interactive" to force the script to WAIT for the job to finish.
#      (Without 'interactive', Abaqus backgrounds itself and the script exits instantly).

/opt/abaqus/Commands/abq2022 job=$JOB_NAME input=$SCRATCH_PATH/$JOB_NAME.inp user=$SCRATCH_PATH/UINTERLEBIM3D_kincxit.f double scratch=$SCRATCH_PATH interactive

# 4. Cleanup (Optional, good for clusters)
echo "--- [Worker Linux] Job finished. ---"