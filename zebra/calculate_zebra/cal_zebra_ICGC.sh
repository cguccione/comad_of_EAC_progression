#!/bin/bash

<<com
Date: 6/6/24
Goal: To run zebra on ICGC data
com

python zebra_filter/calculate_coverages.py -i -o zebra_output.tsv -d zebra_filter/databases/WoL/metadata.tsv

echo 'done'
