#!/bin/bash

INPUT=$1
OUTPUT=$2
COV=$3

bedtools genomecov -split -bg -max $COV -ibam $INPUT | \
    awk -v cov="$COV" '$4 >= cov' | \
    bedtools merge -i - > $OUTPUT
