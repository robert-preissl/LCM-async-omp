#!/bin/bash

# clean old process outputs
rm *.txt

for ((j=1; j <= 10 ; j++))
do
    (/Users/robert/Private/C++Tutorial/LCM/C/send-listener-async 10 $j > p$j.txt) &
    pids[${j}]=$!
done

for pid in ${pids[*]};
do
	wait $pid;
done;

wait;

echo "WARN1 stats - check delays"
grep WARN1 p*txt | wc
echo "WARN2 stats - check immediate collisions"
grep WARN2 p*txt | wc
echo "WARN3 stats - check future collisions"
grep WARN3 p*txt | wc
echo "LOCATION stats - check messages received"
grep LOCATION p*txt | wc
