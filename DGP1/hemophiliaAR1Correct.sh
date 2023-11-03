#!/bin/bash

### Builds the job index
# Create a sequential range
array_values=`seq 100`

# Launch the job and then remove the temporarily created qsub file.
for i in $array_values
do 
# This submits the single job to the resource manager
sbatch hemophiliaAR1Correct.sub $i

done