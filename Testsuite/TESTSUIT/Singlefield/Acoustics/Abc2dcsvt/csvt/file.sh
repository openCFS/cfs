#!/bin/bash
for i in `seq 0 39`; do
   	echo item: $i
	echo $(echo "s(2*100*0.125*$i*3.14/180)" | bc -l) > "data.val$(($i+1))"
	echo $(echo "0" | bc -l) >> "data.val$(($i+1))"
done
