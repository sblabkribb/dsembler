#!/bin/bash
#docker run -d --rm --name dsembler haseong/dsembler:v0.2
docker run -d -v dsembler_debug:/app/ --publish 8000:5000 --name dsembler_test dsembler:v3 
