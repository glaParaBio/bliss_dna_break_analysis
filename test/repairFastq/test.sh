#!/usr/bin/env bash

../../scripts/repairFastq.py -r <(cat R1.fastq) -i <(cat R2.fastq) -s '#'
x=`echo $?`
echo "Check zero exit code: $x"

echo "Fail because reads are not in the same order:"
../../scripts/repairFastq.py -r <(cat R1.fastq) -i <(cat R2.unpaired.fastq) -s '#' > /dev/null
x=`echo $?`
echo "Check non-zero exit code: $x"

echo "Fail because R2 does not contain all reads from R1"
../../scripts/repairFastq.py -r <(cat R1.fastq) -i <(cat R2.fastq | head -n 12) -s '#' > /dev/null
x=`echo $?`
echo "Check non-zero exit code: $x"
