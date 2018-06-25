#!/bin/bash

for inputfile in `cat ./inputfilelist.txt`; do
        echo $inputfile
        ./select_constraints_CLI.R $inputfile $inputfile 
done
