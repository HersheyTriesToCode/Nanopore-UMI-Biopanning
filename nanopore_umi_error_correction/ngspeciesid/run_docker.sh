#!/bin/bash

docker run -it --rm --name umi-pipeline-ngspeciesid-$$ -v ~/src:/usr/src -v ~/data:/usr/data ngspeciesid
