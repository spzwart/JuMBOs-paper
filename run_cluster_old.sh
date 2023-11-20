#!/bin/bash

# Set the number of times you want to run the Python script
times=5  # Change this to the desired number of times

# The name of the Python script
python_script="run_cluster_old.py"

# Loop to run the Python script multiple times
for ((i=1; i<=$times; i++)); do
    log_file="name_${i}.log"
    echo "[$SHELL] ## Running $python_script - Attempt $i"
    python run_cluster_old.py >> "$log_file" 2>&1 &
    echo "[$SHELL] ## $python_script completed - Attempt $i"
done

echo "[$SHELL] ## Script finished"
echo "[$SHELL] ## Job done $DATE"

