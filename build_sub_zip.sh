#!/bin/bash

tb=$(grep time: crt.wdl | head -1 | awk -F'"' '{print $2}')
nb=$(grep node: crt.wdl | head -1 | awk '{print $2}')
pb=$(grep poolname: crt.wdl | head -1 | awk -F'"' '{print $2}') 
cb=$(grep constraint: crt.wdl | head -1 | awk -F'"' '{print $2}') 

echo "---------------------------"
echo "### Using crt.wdl as example ### "
echo "### starting values          ### "
echo node=$nb
echo time=$tb
echo pool=$pb
echo constraint=$cb
echo -e "\n---------------------\n"

echo nodeAfter
read na
if [[ ! $na ]]; then
  nb='shouldnotmatch'
  na='replacewithnothing'
fi

echo timeAfter
read ta
if [[ ! $ta ]]; then
  tb='shouldnotmatch'
  ta='replacewithnothing'
fi

echo poolnameAfter
read pa
if [[ ! $pa ]]; then
  pb='shouldnotmatch'
  pa='replacewithnothing'
fi

echo constraintAfter
read ca
if [[ ! $ca ]]; then
  cb='shouldnotmatch'
  ca='replacewithnothing'
fi

for i in `ls *.wdl`; do 
	echo $i; 
	sed -i "s/time: \"$tb\"/time: \"$ta\"/; s/node: $nb/node: $na/; s/poolname: \"$pb\"/poolname: \"$pa\"/; s/constraint: \"$cb\"/constraint: \"$ca\"/" $i; 
done

time_e=$(grep time: crt.wdl | head -1 | awk -F'"' '{print $2}')
node_e=$(grep node: crt.wdl | head -1 | awk '{print $2}')
pool_e=$(grep poolname: crt.wdl | head -1 | awk -F'"' '{print $2}') 
constraint_e=$(grep constraint: crt.wdl | head -1 | awk -F'"' '{print $2}') 

echo '----------------------'
echo "### Using crt.wdl as example ###"
echo "### values changed to        ### "
echo node=$node_e
echo time=$time_e
echo pool=$pool_e
echo constraint=$constraint_e
echo "---------------------"

time_label=$(grep time: crt.wdl | head -1 | awk -F'"' '{print $2}' | awk -F':' '{print $1}')
zip subs_n${node_e}_t${time_label}_${pool_e}_${constraint_e}.zip crt.wdl functional-annotation.wdl genemark.wdl prodigal.wdl rfam.wdl structural-annotation.wdl trnascan.wdl

echo "*** created file subs_n${node_e}_t${time_label}_${pool_e}_${constraint_e}.zip ***"
