#!/bin/bash
for i in *.jaccard; do printf "\t%s" "$(basename ${i} .jaccard)"; done
printf "\n"
for i in *.jaccard
	do
	printf "%s\t" "$(basename $i .jaccard)"
        for j in *.jaccard
		do
			x=$(bedtools jaccard -a ${i} -b ${j} | tail -1 | awk '{print $3}')
			printf "%s\t" "$x"
	done
	printf "\n"
done

