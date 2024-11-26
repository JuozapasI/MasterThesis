#!/bin/bash

( samtools view -H $1  && samtools view -F 1024 -F 256 $1 | \
	grep "GN:Z:-" | \
	grep -P "\bNH:i:1\b" | \
	grep -E "CB:Z:[A-Z]{$3}" | \
	grep -E "UB:Z:[A-Z]{$4}" ) | \
	samtools view -h -b - > $2
