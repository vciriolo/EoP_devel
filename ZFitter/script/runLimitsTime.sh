#!/bin/bash
file=data/runRanges/monitoring_2017.dat
echo "set x2tics (\\" > tmp/tics.label
grep -v '#' $file | sed 's|-| |g' | awk '{printf("\"%d-%d\" %d 2 ,\\\n", $1, $2, ($4+$5)/2)}' >> tmp/tics.label
echo ")" >> tmp/tics.label



exit 0
#
file=data/runRanges/runRangeLimits.dat


# output to a file with
# run timeMin label

# remove the commented lines
sed '/^#/ d;' $file > tmp/runLimit.clean

#find the timeMin
for run in `cat tmp/runLimit.clean | sed 's| #.*||'`
do
	echo -ne $run"\t"
	grep $run data/runRanges/monitoring_2017.dat 
#| sed 's|.*:||;s|#.*||' | awk '(NF==3){print $3}' | sed 's|-.*||' | sort | uniq
	echo 
done 
#| awk '(NF==2){print $0}' > tmp/runLimitTime
exit 0
echo "set x2tics (\\" > $file.label
for run in `cut -f 1 tmp/runLimitTime`
do
	label=`grep $run $file | sed 's|.*#||'`
	time=`grep $run tmp/runLimitTime | awk '{print $2}'`
	
	echo "\"$run: $label\" $time,\\"
done | sed 's|,$||' >> $file.label
echo ")" >> $file.label

index=1
for run in `cut -f 1 tmp/runLimitTime`
do
	time=`grep $run tmp/runLimitTime | awk '{print $2}'`
	echo "set arrow $index from $time,graph 0 to $time,graph 1 nohead dt 2" >> $file.label
	let index=$index+1
done


