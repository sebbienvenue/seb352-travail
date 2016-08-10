#!/bin/bash
#change un pk en pdf tree
#expect an input in input the pk needed
clear
input=$1
if [[ $# -eq 0  ]]; then 
echo "No inputs"
echo "Write the proper PK here: "
read input
fi

#echo "change the pk in pdf tree"
#echo "Please enter a PK: "


#read pk
verdi graph generate $input	#$pk
dot -Tpdf -o res/$input.pdf $input.dot
rm $input.dot
