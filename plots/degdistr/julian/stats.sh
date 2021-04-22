#!/bin/bash
DATADIR=/gscratch/niac/tony/julian/
graphs=(com-orkut-hygra friendster-hygra livejournal-hygra orkut-group-hygra rand3-hygra rand1-hygra web-hygra)
column=(0 1)
for g in "${graphs[@]}"; do
file=${DATADIR}${g}
for col in "${column[@]}"; do
#echo

/usr/lusers/x0/soverlap/stats/hystat.exe -f ${file} -o ${g}${col}dd -D ${col}

done #D
done #data