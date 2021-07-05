#!/bin/bash
#docker run -d --rm --name dsembler haseong/dsembler:v0.2
docker run -d --rm -v /mnt/beta/dockervolume/dsembler_output_script:/app/output_script \
	           -v /mnt/beta/dockervolume/dsembler_output:/app/output \
		   --publish 5001:5000 --name dsembler_test dsembler:v6
