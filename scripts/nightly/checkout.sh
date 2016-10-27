#!/bin/bash

echo "# INIT REPO"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load git

rm -rf * .*

if test "$(ls -A .)"; then
    echo "# Updating existing"
    git fetch
    git checkout origin/${bamboo.planRepository.branchName}
else
    echo "## init"
    git init .
    echo "## Add remote"
    git remote add -t ${bamboo.planRepository.branchName} -f origin ${bamboo.planRepository.repositoryUrl}
    echo "## Checkout"
    git checkout ${bamboo.planRepository.branchName}
fi

git clean -fd
