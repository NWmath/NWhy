#!/bin/bash
DATADIR=/gscratch/niac/tony/real-world/
graphs=(IMDB.mtx NotreDame_actors.mtx 12month1.mtx connectus.mtx dns.mtx)
column=(0 1)
for g in "${graphs[@]}"; do
file=${DATADIR}${g}
for col in "${column[@]}"; do
#echo

/usr/lusers/x0/soverlap/stats/hystat.exe -f ${file} -o ${g}${col}dd -D ${col}

done #D
done #data