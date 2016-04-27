#!/bin/bash
#chmod 755 run_bdecay.sh (do this only ONCE to make the file executable)
#./run_bdecay.sh (do this to execute this file)

echo "Enter the name of the rootfile you will output WITHOUT '.root' (e.g. b_decay_histo): "
read filename
root -l -q 'bdecay.cpp("'$filename'")'
root -l 'bdecay2.cpp("'$filename'")'
