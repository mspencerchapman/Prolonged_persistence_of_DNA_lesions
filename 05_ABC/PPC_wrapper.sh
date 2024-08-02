#!/bin/bash
# running in an array, using $LSB_JOBINDEX
MYINPUT=$(sed -n "${LSB_JOBINDEX}p" posterior_means.txt)
echo "Running the PPC script for lesion duration = $MYINPUT"
Rscript ABC_simulation_new_PPCs.R ${MYINPUT}
