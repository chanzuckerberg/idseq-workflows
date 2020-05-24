#!/bin/bash

set -euo pipefail

for tag in $(git tag); do
    for file in $(git ls-tree -r --name-only "$tag" | grep '.wdl$'); do
        git show "${tag}:${file}" | aws s3 cp - "s3://idseq-workflows/${tag}/${file}"
    done
done
