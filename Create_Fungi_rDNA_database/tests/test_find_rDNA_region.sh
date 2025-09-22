#!/bin/bash
# Test script for find_rDNA_region.py
python ../scripts/find_rDNA_region.py test_barrnap_output.gff > test_output.txt

cat test_output.txt
