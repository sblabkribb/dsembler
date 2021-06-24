#!/bin/bash
#docker run -d --rm --name dsembler haseong/dsembler:v0.2
docker run -it -v dsembler:/app/output_script/ --name dsembler_script dsembler:v2 /bin/sh
