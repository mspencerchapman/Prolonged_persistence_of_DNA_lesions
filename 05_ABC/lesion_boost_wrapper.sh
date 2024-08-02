i=$1
echo "Running the boost lesions script for lesion duration = $i"
Rscript ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R ${LSB_JOBINDEX} $i
