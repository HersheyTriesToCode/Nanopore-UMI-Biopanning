#!/bin/bash

docker run --platform linux/x86_64 -it --rm --name pipeline-umi-amplicon-$$ -v ~/src:/usr/src -v ~/data:/usr/data pipeline-umi-amplicon
