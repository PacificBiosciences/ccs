#!/bin/bash
set -eo pipefail

echo "## Fetch submodules"
source /mnt/software/Modules/current/init/bash
module load git

# Bamboo's checkout of unanimity doesn't set the "origin" remote to
# something meaningful, which means we can't resolve the relative
# submodules.  Override the remote here.
git remote set-url origin ssh://git@bitbucket.nanofluidics.com:7999/sat/unanimity.git

#git submodule update --init --remote

echo "BRANCH NAME: ${BRANCH_NAME}"
SUBMODULE_BRANCH="develop"
if [ "${BRANCH_NAME}" == "master" ]; then
  SUBMODULE_BRANCH="master"
fi
echo "CHECKING OUT ${SUBMODULE_BRANCH}"
git submodule foreach git pull origin ${SUBMODULE_BRANCH}