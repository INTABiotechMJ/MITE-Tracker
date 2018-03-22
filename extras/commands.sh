
blastn -query results/olds/finder.fasta      -max_target_seqs 1     -subject ../data/oryrep.ref   -outfmt 6 > ../data/mite-tracker/finder_repbase.csv   
blastn -query results/olds/finder.fasta      -max_target_seqs 1     -subject ../data/repbase_nonautonomous_rice.50_800.fa   -outfmt 6 > ../data/mite-tracker/finder_repbase_filtered.csv

blastn -query results/olds/detectmite/rice.mite.fasta     -max_target_seqs 1     -subject ../data/oryrep.ref   -outfmt 6 > ../data/mite-tracker/detect_repbase.csv
blastn -query results/olds/detectmite/rice.mite.fasta     -max_target_seqs 1     -subject ../data/repbase_nonautonomous_rice.50_800.fa   -outfmt 6 > ../data/mite-tracker/detect_repbase_filtered.csv

blastn -query results/olds/rice_wz_6_mcn_4/families_nr.fasta   -max_target_seqs 1  -subject ../data/oryrep.ref   -outfmt 6 > ../data/mite-tracker/tracker_repbase.csv
(venv) ➜  extras git:(master) ✗ blastn -query results/olds/rice_wz_6_mcn_4/families_nr.fasta   -max_target_seqs 1  -subject ../data/repbase_nonautonomous_rice.50_800.fa  -outfmt 6 > ../data/mite-tracker/tracker_repbase_filtered.csv


blastn -query results/olds/finder.fasta -max_target_seqs 1 -subject ../data/otrep.fa  -outfmt 6 > ../data/mite-tracker/finder_trep.csv

blastn -query results/olds/detectmite/rice.miteSet.1.fasta -max_target_seqs 1 -subject ../data/otrep.fa  -outfmt 6 > ../data/mite-tracker/detect_trep.csv

blastn -query results/olds/rice_wz_6_mcn_4/families_nr.fasta -max_target_seqs 1 -subject ../data/otrep.fa  -outfmt 6 > ../data/mite-tracker/tracker_trep.csv