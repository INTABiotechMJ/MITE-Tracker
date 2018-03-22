# Information about scripts in extras/ directory
This scripts were used to analyze and compare outputs of MITE Tracker and other tools

#run_wheat.sh
Used to find MITEs in wheat genome. First it was separated by chromosome using FAspliter.py. Then each chromosome was processed in order to find candidates separately. All candidates files were merged in the same directory and the cluster task was executed once. 