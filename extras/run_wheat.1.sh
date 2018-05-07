cd ../
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr1A.fasta -w 4 -j IWGSC_1A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr1B.fasta -w 4 -j IWGSC_1B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr1D.fasta -w 4 -j IWGSC_1D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr2A.fasta -w 4 -j IWGSC_2A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr2B.fasta -w 4 -j IWGSC_2B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr2D.fasta -w 4 -j IWGSC_2D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr3A.fasta -w 4 -j IWGSC_3A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr3B.fasta -w 4 -j IWGSC_3B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr3D.fasta -w 4 -j IWGSC_3D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr4A.fasta -w 4 -j IWGSC_4A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr4B.fasta -w 4 -j IWGSC_4B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr4D.fasta -w 4 -j IWGSC_4D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr5A.fasta -w 4 -j IWGSC_5A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr5B.fasta -w 4 -j IWGSC_5B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr5D.fasta -w 4 -j IWGSC_5D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr6A.fasta -w 4 -j IWGSC_6A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr6B.fasta -w 4 -j IWGSC_6B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr6D.fasta -w 4 -j IWGSC_6D --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr7A.fasta -w 4 -j IWGSC_7A --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr7B.fasta -w 4 -j IWGSC_7B --task candidates
python3 -m MITETracker -g /media/crescentejuan/Data/iwgsc_by_chr/chr7D.fasta -w 4 -j IWGSC_7D --task candidates
mkdir results/IWGSC
cat results/IWGSC_*/candidates.csv > results/IWGSC/candidates.csv
cat results/IWGSC_*/candidates.fasta > results/IWGSC/candidates.fasta
python MITETracker.py -g none -w 4 -j IWGSC --task cluster --min_copy_number 6