#!/bin/sh

i=5
while [ $i -le 150 ]
do
    distance=$(echo "scale=3; $i/1000" | bc)
    echo $distance
    if [ $(echo "$distance == $(printf "%.2f" $distance)" | bc) -eq 1 ]; then
        distance=$(printf "%.2f" $distance) ## remove trailing zero
        SFTPATH="../Frames${distance}/sfts_8_5sec"
        framecache="../Frames${distance}/framecache${distance}" ## Path of the framecache file
    else
        distance=$(printf "%.3f" $distance)
        SFTPATH="../Frames${distance}/sfts_8_5sec"
        framecache="../Frames${distance}/framecache${distance}" ## Path of the framecache file
    fi
	gps_start=1000000000 ## GPS start time
	gps_end=1000040960 ## GPS end time
	Tseg=8.5 ## SFT segment duration
	Band=512 ## Frequency band
	fmin=61.5 ## minimum frequency of SFT

    lalpulsar_MakeSFTs -f $fmin -t $Tseg -p $SFTPATH -C $framecache -s $gps_start -e $gps_end -N H1:STRAIN_SSSM -v 2 -i H1 -u PROC_REAL8 -w 0 -F 0 -B $Band -D 4 -X MSFT

    i=$(($i+2))
done
