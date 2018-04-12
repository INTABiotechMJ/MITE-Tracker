# Information about scripts in extras/ directory
This scripts were used to analyze and compare outputs of MITE Tracker and other tools

#run_wheat.sh
Used to find MITEs in wheat genome. First it was separated by chromosome using FAspliter.py. Then each chromosome was processed in order to find candidates separately. All candidates files were merged in the same directory and the cluster task was executed once. 

blastn -task blastn -query data/tracker/example.1.fasta   -subject ../../data/IRGSP-1.0_genome.fasta  -outfmt 6  -evalue 10e-3  | awk -F "\t" '$4 > 240 { print $0 }'

#mping in tracker
blastn -query data/tracker/mping.fa    -subject data/tracker/all.fasta  -outfmt 6  
MPING   MITE_T_728|chr06|13737614|13738050|TAA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_726|chr04|33021864|33022300|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_725|chr05|19328614|19329050|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_724|chr08|16683584|16684020|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_722|chr02|617945|618381|TTA|16|F73       100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_721|chr08|1019668|1020104|TTA|16|F73     100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_720|chr12|1045459|1045895|TAA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_718|chr07|4560971|4561407|tAA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_717|chr08|20442412|20442848|TAA|19|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_716|chr01|23332543|23332979|TAA|16|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_715|chr02|22549111|22549547|TTA|19|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_714|chr08|4712966|4713402|TAA|19|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_713|chr01|17513830|17514266|TTA|19|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_712|chr03|6513585|6514021|TTA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_711|chr11|23200101|23200537|TTA|19|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_709|chr03|5504295|5504731|TAA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_708|chr02|29244323|29244759|TAA|16|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_707|chr12|3285783|3286219|TTA|19|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_706|chr08|20674233|20674669|TCA|16|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_705|chr02|214433|214869|TAA|19|F73       100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_703|chr03|21026568|21027004|TAA|19|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_702|chr01|25261108|25261544|TTA|16|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_700|chr05|22235590|22236026|TAA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_699|chr03|9240070|9240506|TTA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_698|chr01|24779767|24780203|TAA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_697|chr04|34688302|34688738|TAA|16|F73   100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_696|chr11|393594|394030|TTA|19|F73       100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_695|chr01|29931513|29931949|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_694|chr02|13161934|13162370|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_693|chr03|9427116|9427552|TTA|16|F73     100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_692|chr10|4320417|4320853|TTA|16|F73     100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_691|chr03|17575713|17576149|TAA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_690|chr05|18747494|18747930|TTA|16|F73   100.000 430     0       0       1       430     433     4       0.0     795
MPING   MITE_T_688|chr12|839600|840036|TTA|16|F73       100.000 430     0       0       1       430     4       433     0.0     795
MPING   MITE_T_686|chr03|9568665|9569103|CTAAg|20|F73   100.000 430     0       0       1       430     434     5       0.0     795
MPING   MITE_T_683|chr08|28186234|28186675|ta|21|F73    100.000 430     0       0       1       430     436     7       0.0     795
MPING   MITE_T_682|chr06|18136413|18136860|TA|25|F73    100.000 430     0       0       1       430     9       438     0.0     795
MPING   MITE_T_681|chr12|9951097|9951545|AT|26|F73      100.000 430     0       0       1       430     438     9       0.0     795
MPING   MITE_T_727|chr03|12735752|12736188|TTA|19|F73   99.767  430     1       0       1       430     4       433     0.0     789
MPING   MITE_T_723|chr12|2734537|2734973|TAA|16|F73     99.767  430     1       0       1       430     4       433     0.0     789
MPING   MITE_T_719|chr09|694479|694915|TAA|19|F73       99.767  430     1       0       1       430     433     4       0.0     789
MPING   MITE_T_710|chr04|35421802|35422238|TAA|16|F73   99.767  430     1       0       1       430     4       433     0.0     789
MPING   MITE_T_704|chr04|19021056|19021492|TTA|16|F73   99.767  430     1       0       1       430     433     4       0.0     789
MPING   MITE_T_701|chr06|30099534|30099970|TAA|16|F73   99.767  430     1       0       1       430     4       433     0.0     789
MPING   MITE_T_689|chr02|28008337|28008773|TAA|16|F73   99.767  430     1       0       1       430     433     4       0.0     789
MPING   MITE_T_687|chr09|16698137|16698573|TAA|16|F73   99.767  430     1       0       1       430     4       433     0.0     789
MPING   MITE_T_684|chr03|34592539|34592979|TA|20|F73    99.767  430     1       0       1       430     435     6       0.0     789
MPING   MITE_T_685|chr10|21716385|21716824|TATT|16|F73  99.767  430     0       1       1       430     6       434     0.0     787

#mping in detectmite
blastn -query data/tracker/mping.fa    -subject data/detectmite/rice.mite.fasta  -outfmt 6  
MPING   chr01|23332547|23332976|15|3|51 100.000 430     0       0       1       430     1       430     0.0     795
MPING   chr02|29244324|29244759|18|3|3  100.000 430     0       0       1       430     4       433     0.0     795
MPING   chr01|23332551|23332972|11|4|49 100.000 422     0       0       5       426     1       422     0.0     780