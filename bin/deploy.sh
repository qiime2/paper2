#!/usr/bin/env bash

aws s3 cp \
    build/html/data s3://qiime2-curr-protoc-bioinformatics/ \
    --recursive \
    --acl 'public-read' && \
cp -r build/html/* docs/
