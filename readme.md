
## About

MITE Tracker: an accurate method for identifying miniature inverted-repeat transposable elements in large genomes. 

An efficient and easy to run tool for discovering Miniature Inverted repeats Transposable Elements (MITEs) in genomic sequences. It is written in python 3 and uses ncbi's blast+ for finding inverted repeats and cdhit to do the clustering. 

Large genomes can be processed in desktop computers.

# Requirements
 - tested in macOS 10.13.1, Debian 7.6, Ubuntu 16.04, Windows 7
 - ncbi blast+ (Nucleotide-Nucleotide BLAST 2.6.0+)
 - python requirements are in requirements.txt file (bipython and pandas)

# Installation and running

## 
```
# clone repo
git clone https://github.com/INTABiotechMJ/MITE-Tracker.git
cd MITE-Tracker

# blast
sudo apt-get install ncbi-blast+ virtualenv
# in macOS: brew install ncbi-blast+ virtualenv


#vsearch
wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz
tar xzf v2.7.1.tar.gz
cd vsearch-2.7.1
#might need sudo apt-get install autoconf
sh autogen.sh
./configure
make

#python dependencies
cd ..
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt

# running

python3 -m MITETracker -g /path/to/your/genome.fasta -w 3 -j jobname

# or to run in background

nohup python -u miteParser.py -g /path/to/your/genome.fasta -w 3 -j jobname &

```

In order to check the output and progress you can use these command (ctrl+c to exit)
```
#nohup will have the program output as well as the output from cdhit execution
tail -f nohup.out
#out.log contaings a log file with timing information
tail -f results/[jobname]/out.log
```

# Command line options
| Argument  | Description | Data type  | Required or default |
| ------------- | ------------- | ------------- | ------------- |
| -g  | Genome file in fasta format  | string  | required  |
| -j  | Jobname. Result files will be created in results/jobname   | string  | required  |
| -w  | Max number of processes to use simultaneously  | int  | 1  |
| -tsd_min_len  | TSD min lenght  | int  | 2  |
| -tsd_max_len  | TSD max lenght  | int  | 10  |
| -mite_min_len  | MITE min lenght  | int  | 50  |
| -mite_max_len  | MITE max lenght  | int  | 650  |
| --task  | cluster or candidates  | string |   |


# Results
All the results are placed in _results/[yourjobname]/_. 
Here you will find:
    _families.fasta_ all the MITEs sequences divided by families (custom format)
    _families_nr.fasta_ with one MITE per family in fasta format
    _all.fasta_ all MITEs in fasta format
    _all.gff3_  a gff file describing all MITEs found 

# Troubleshooting
If getting any error while running the BLASTn searches please check you blast+ version

# Running large genomes in different computers
This is an example of how we run wheat genome. Each chromosome can be run separately (--task candidates) in a different computers. Results should be merged together using _cat_ and then run the cluster command (--task cluster)

21 wheat chromosomes were downloaded in different files. 

```
python3 -m MITETracker -g /media/iwgsc_by_chr/chr1A.fasta -w 2 -j IWGSC_1A --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr1B.fasta -w 2 -j IWGSC_1B --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr1D.fasta -w 2 -j IWGSC_1D --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr2A.fasta -w 2 -j IWGSC_2A --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr2B.fasta -w 2 -j IWGSC_2B --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr2D.fasta -w 2 -j IWGSC_2D --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr3A.fasta -w 2 -j IWGSC_3A --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr3B.fasta -w 2 -j IWGSC_3B --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr3D.fasta -w 2 -j IWGSC_3D --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr4A.fasta -w 2 -j IWGSC_4A --task candidates
python3 -m MITETracker -g /media/iwgsc_by_chr/chr4B.fasta -w 2 -j IWGSC_4B --task candidates
...
mkdir results/IWGSC
cat results/IWGSC_*/candidates.csv > results/IWGSC_2/candidates.csv
cat results/IWGSC_*/candidates.fasta > results/IWGSC_2/candidates.fasta
python3 -m MITETracker -g none -w 3 -j IWGSC_2 --task cluster --min_copy_number 4
```