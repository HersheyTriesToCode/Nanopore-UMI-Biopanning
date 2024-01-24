#!/bin/bash

docker run -it --rm --name umi-pipeline-extract-$$ -v ~/src:/usr/src -v ~/data:/usr/data umi
