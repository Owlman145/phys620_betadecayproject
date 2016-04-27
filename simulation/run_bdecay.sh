echo "Enter the name of the rootfile you will output WITHOUT '.root' (e.g. b_decay_histo): "
read filename
root -l -q 'bdecay.cpp("'$filename'")'
root -l 'bdecay2.cpp("'$filename'")'
