cd ../
#python extras/FAspliter.py -s ../data/Triticum_aestivum.TGACv1.dna.toplevel.fa -o ../data/TGAC_by_chr/ -l 1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D U
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/1A.fasta -w 8 -j TGAC_1A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/1B.fasta -w 8 -j TGAC_1B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/1D.fasta -w 8 -j TGAC_1D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/2A.fasta -w 8 -j TGAC_2A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/2B.fasta -w 8 -j TGAC_2B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/2D.fasta -w 8 -j TGAC_2D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/3A.fasta -w 8 -j TGAC_3A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/3B.fasta -w 8 -j TGAC_3B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/3D.fasta -w 8 -j TGAC_3D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/4A.fasta -w 8 -j TGAC_4A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/4B.fasta -w 8 -j TGAC_4B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/4D.fasta -w 8 -j TGAC_4D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/5A.fasta -w 8 -j TGAC_5A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/5B.fasta -w 8 -j TGAC_5B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/5D.fasta -w 8 -j TGAC_5D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/6A.fasta -w 8 -j TGAC_6A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/6B.fasta -w 8 -j TGAC_6B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/6D.fasta -w 8 -j TGAC_6D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/7A.fasta -w 8 -j TGAC_7A --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/7B.fasta -w 8 -j TGAC_7B --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/7D.fasta -w 8 -j TGAC_7D --task candidates
python MITETracker.py -g /media/crescentejuan/Data/TGAC_by_chr/U.fasta -w 8 -j TGAC_U --task candidates
mkdir results/TGAC
cat results/TGAC_*/candidates.csv > results/TGAC/candidates.csv
cat results/TGAC_*/candidates.fasta > results/TGAC/candidates.fasta
python MITETracker.py -g none -w 4 -j TGAC --task cluster --min_copy_number 6