#!/bin/bash
#change un pk en pdf tree
#take in input the pk needed
clear
echo "change the pk in pdf tree"
echo "Please enter a PK: "
read pk
verdi graph generate $pk
dot -Tpdf -o $pk.pdf $pk.dot
rm $pk.dot
