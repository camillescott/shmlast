#!/usr/bin/bash
REPS=3

QUERY=sacPom.cdna.fa
DB=Nematostella_vectensis.uniprot.pep.fa

echo "# `date --rfc-3339=seconds` shmlast-`shmlast --version | cut -f 2 -d " "` crb-blast-`crb-blast --version`" > results.csv

for n_threads in {1,2,4}; do
    for ((n=0; n<$REPS; n++)); do
        /usr/bin/time -o results.csv -a --format "shmlast, $n_threads, $n, %E" shmlast crbl -q $QUERY -d $DB  --n_threads $n_threads -e 1e-5
        shmlast crbl -q $QUERY -d $DB --action clean
        rm -f *.blast evalues_data fitting_data
        /usr/bin/time -o results.csv -a --format "crb-blast, $n_threads, $n, %E" crb-blast -q $QUERY -t $DB -h $n_threads -e 1e-5 -s -o test.tsv -v 
    done
done
