#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --output=test_hw_openmp_%j.stdout
#SBATCH --error=test_hw_openmp_%j.stderr
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#Please set up this environment variable such that it points to
#the proper ABSOLUTE path where the compiled test_hw_openmp
DRIVER_PROGRAM_ABS_PATH=$HOME/Rafael/lab/Assignment/Assignment_1/Assignment_1/test_hw_openmp

#The value of this environment variable should be modified in order
#to play around different SCHEDULE STRATEGIES and several chunk sizes
#(as required by the OpenMP HW)
export OMP_SCHEDULE=STATIC

export OMP_DISPLAY_ENV=VERBOSE
for NUM_THREADS in 1 2 4 8 16 24
do
   export OMP_NUM_THREADS=$NUM_THREADS
   echo "*****************"
   echo "BEGIN OMP_NUM_THREADS=$OMP_NUM_THREADS RESULTS"
   OMP_PLACES_VALUE="{0}"
   i=1
   while [ "$i" -lt "$NUM_THREADS" ]
   do
     OMP_PLACES_VALUE="$OMP_PLACES_VALUE,{$i}"
     let i=i+1
   done
   export OMP_PLACES=$OMP_PLACES_VALUE 

   #Execute driver program
   $DRIVER_PROGRAM_ABS_PATH   

   echo "END OMP_NUM_THREADS=$OMP_NUM_THREADS RESULTS"
   echo "*****************"
   file_num2=`printf "%02d\n" $NUM_THREADS`
   ./upload_source.sh Assignment/Assignment_1 results$file_num2.m $1
done

