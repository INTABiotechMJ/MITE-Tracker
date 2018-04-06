python compare_db.py -d data/repbase/repbase_nonautonomous.fasta  -m data/tracker/families_nr.fasta -p MITE_Tracker -n Repbase -l A
python compare_db.py  -d data/repbase/repbase_nonautonomous.fasta  -m data/mitehunter/all.fa -p MITEHunter -n Repbase -l C          
python compare_db.py  -d data/repbase/repbase_nonautonomous.fasta  -m data/detectmite/rice.mite.fasta -p detectMITE -n Repbase -l B 
python compare_db.py  -d data/detectmite/rice.mite.fasta  -m data/tracker/families_nr.fasta -p MITE_Tracker -n detectMITE -l A      
python compare_db.py  -d data/mitehunter/all.fa   -m data/tracker/families_nr.fasta -p MITE_Tracker -n MITEHunter -l B