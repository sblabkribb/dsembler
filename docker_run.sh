#!/bin/bash
#docker run -d --rm --name dsembler haseong/dsembler:v0.2
docker run --rm -d -v dsembler_debug:/app/ --publish 5000:5000 --name dsembler dsembler:v3 
