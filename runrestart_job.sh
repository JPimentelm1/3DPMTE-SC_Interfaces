#!/bin/bash

# -----------------------------------------------------------------------
# ROBUST ABAQUS RESTART LAUNCHER (LINUX)
# -----------------------------------------------------------------------

# 1. Capture Arguments
# "$1" is the first argument, "$2" is the second.
# We wrap them in quotes to handle potential spaces safely.
JOB_NAME="$1"
OLDJOB_NAME="$2"
CURRENT_DIR="$PWD"

# 2. Define Abaqus Command Path
# Update this path if your cluster uses a different version or alias (e.g., 'abaqus')
ABAQUS_CMD="/opt/abaqus/Commands/abq2022"

# 3. Sanity Checks
# Unlike Windows, Linux is case-sensitive. Ensure filenames match exactly.

# Check for Subroutine
if [ ! -f "${CURRENT_DIR}/UINTERLEBIM3D_kincxit.f" ]; then
    echo "[ERROR] Subroutine file not found in: ${CURRENT_DIR}"
    exit 1
fi

# Check for Old Job Restart File (.res)
# This is critical for the 'oldjob=' parameter to work.
if [ ! -f "${CURRENT_DIR}/${OLDJOB_NAME}.res" ]; then
    echo "[ERROR] Old job restart file '${OLDJOB_NAME}.res' not found."
    exit 1
fi

# 4. The Execution Command
# - We enclose all paths in "${VAR}" to handle spaces (e.g., "My Project").
# - 'interactive' forces the script to wait until Abaqus finishes.
# - 'double' specifies double precision (standard for crack propagation).

"${ABAQUS_CMD}" job="${JOB_NAME}" \
    oldjob="${OLDJOB_NAME}" \
    input="${CURRENT_DIR}/${JOB_NAME}.inp" \
    user="${CURRENT_DIR}/UINTERLEBIM3D_kincxit.f" \
    double \
    scratch="${CURRENT_DIR}" \
    interactive

# 5. Capture and Return Exit Code
# $? stores the exit code of the last executed command (Abaqus).
EXIT_CODE=$?

echo "--- Abaqus finished with exit code: $EXIT_CODE ---"
exit $EXIT_CODE