#!/bin/bash

for i in {4..23}
do
	sbatch ../scripts/calcie.sh elecexc.py "folder_$i/" "exc_$i.out"
done
