#!/bin/bash
while [[ 1 ]];
do
	ps -F --forest -p $(pgrep -x degen_primer) >> degen_primer-ps.out 2>/dev/null; 
	sleep 1; 
done;
