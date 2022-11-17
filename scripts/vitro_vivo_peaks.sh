#!/bin/bash

# create header for counts
printf "X," > overlap/counts.csv
for g in invivo/*.bed
do
    dest=`basename $g`
    printf "%s," $dest >> overlap/counts.csv
done

echo "" >> overlap/counts.csv

for f in invitro/*.bed
do
    source=`basename $f`
    echo "[SCRIPT]: Processing ${source}"
    printf "%s," $source >> overlap/counts.csv
    for g in invivo/*.bed
    do
        dest=`basename $g`
        count=$(bedtools intersect -a $f -b $g | wc -l)
        if [[ $count -ne 0 ]]; then
            bedtools intersect -a $f -b $g > "overlap/${source%.*}_vs_${dest%.*}.bed"
        fi
        printf "%s," $count >> overlap/counts.csv
    done
    echo "" >> overlap/counts.csv
done
