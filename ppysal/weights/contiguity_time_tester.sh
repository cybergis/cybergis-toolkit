#!/bin/bash
#echo "Bash version ${BASH_VERSION}..."
for c in 2 4 8 16 32
    do
        python contiguity_apply_async.py $c
        python contiguity_joinable_queue.py $c
        python contiguity_managed_dict.py $c
        python contiguity_map_async.py $c
    done

python contiguity_serial.py

