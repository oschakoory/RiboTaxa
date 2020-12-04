#!/bin/sh


echo -n "" > test.log
exec 3>&1 1>> test.log 

conda activate MicroTaxa_py27

echo "Microtaxa env activated" 1>&3

vsearch --cluster_fast "$DB_DIR" 
	
conda deactivate 

echo "Microtaxa env deactivated" 