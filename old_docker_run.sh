#!/bin/bash
docker run -d --rm \
	-p 5000:5000 \
	-v /home/haseong/alpha/dev/dsembler:/home/dsembler \
  	--name dsembler \
	haseong/dsembler:v0.2
