#!/bin/bash

( samtools view -H $1 && samtools view -F 1024 -F 256 $1 | \
	grep "GN:Z:-" | \
	grep -P "\bNH:i:1\b" | \
	grep -E "CB:Z:[A-Z_]{$3}" | \
	grep -E "UB:Z:[A-Z_]{$4}" ) | \
	# if removing comments, remove ')' above
	# awk 'BEGIN {FS="\t"; OFS="\t"} {if ($0 ~ /CB:Z:/ && $0 ~ /UB:Z:/) { \
	# split($0, a, "CB:Z:"); split(a[2], b, " "); cb = b[1]; \
	# split($0, a, "UB:Z:"); split(a[2], b, " "); ub = b[1]; \
	# if (!seen[cb":"ub]++) { print; } } }' ; ) | \
	samtools view -h -b - > $2
