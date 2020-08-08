#!/bin/bash

set -euo pipefail

if [[ $# != 3 ]]; then
    echo "This script creates and pushes a git tag representing a new workflow release."
    echo "Usage: $(basename $0) workflow_name release_type release_notes"
    echo "Example: $(basename $0) short-read-mngs patch 'Fixing a bug'"
    exit 1
fi

export WORKFLOW_NAME=$1 RELEASE_TYPE=$2 RELEASE_NOTES="$3"

if ! [[ -d "$(dirname $0)/../$WORKFLOW_NAME" ]]; then
    echo "Unable to locate workflow under a top level directory $WORKFLOW_NAME" 1>&2
    exit 1
fi

if [[ $RELEASE_TYPE == major ]]; then
    TAG=$(git describe --tags --match "${WORKFLOW_NAME}-*" | perl -ne '/(.+)-(\d+)\.(\d+)\.(\d+)/; print "$1-@{[$2+1]}.0.0"')
elif [[ $RELEASE_TYPE == minor ]]; then
    TAG=$(git describe --tags --match "${WORKFLOW_NAME}-*" | perl -ne '/(.+)-(\d+)\.(\d+)\.(\d+)/; print "$1-$2.@{[$3+1]}.0"')
elif [[ $RELEASE_TYPE == patch ]]; then
    TAG=$(git describe --tags --match "${WORKFLOW_NAME}-*" | perl -ne '/(.+)-(\d+)\.(\d+)\.(\d+)/; print "$1-$2.$3.@{[$4+1]}"')
else
    echo "RELEASE_TYPE should be one of major, minor, patch" 1>&2
    exit 1
fi

git tag "$TAG"
git push origin "$TAG"
