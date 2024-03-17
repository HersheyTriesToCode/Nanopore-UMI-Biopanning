#!/bin/sh

#export DOCKER_DEFAULT_PLATFORM=linux/amd64
docker build --platform linux/x86_64 -t pipeline-umi-amplicon .
