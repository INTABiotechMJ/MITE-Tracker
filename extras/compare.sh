python compare_db.py  -e2 data/mitehunter/all.fa  -e1 data/tracker/families_nr.fasta -e1n MITE_Tracker -e2n MITE_Hunter -l B -o news

python compare_db.py  -e2 data/repbase/repbase_nonautonomous.fasta   -e1 data/detectmite/rice.mite.fasta  -e1n detectMITE -e2n Repbase -l B -o news

python compare_db.py  -e2 data/repbase/repbase_nonautonomous.fasta   -e1 data/mitehunter/all.fa  -e1n MITE_Hunter -e2n Repbase -l A -o newsdsa