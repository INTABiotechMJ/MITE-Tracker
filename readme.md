## About

IRmatcher is efficient and easy to run tool for discovering Miniature Inverted repeats Transposable Elements (MITEs) in genomic sequences. It is written in python and uses ncbi's blast+ for finding inverted repeats and cdhit to do the clustering. 

Large genomes can be executed in desktop computers.

# Requirements (just follow how to install)
 - tested in macOS 10.13.1, Debian 7.6, Ubuntu 16.04, Windows 7
 - ncbi blast+ (Nucleotide-Nucleotide BLAST 2.6.0+)
 - python requirements are in requirements.txt file (bipython and pandas)

# Installation and running

## 
```
# clone repo
git clone https://juancresc@bitbucket.org/juancresc/irmatcher.git
cd irmatcher

# blast
sudo apt-get install ncbi-blast+ virtualenv
# in macOS: brew install ncbi-blast+ virtualenv

#vsearch
wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz
tar xzf v2.7.1.tar.gz
cd vsearch-2.7.1
./autogen.sh
./configure
make

#python dependencies
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt

# running

python miteParser.py -g /path/to/your/genome.fasta -w 3 -j jobname

# or to run in background

nohup python -u miteParser.py -g /path/to/your/genome.fasta -w 3 -j jobname &
```

In order to check the program output you can use these command (ctrl+c to exit)
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


# Results
All the results are placed in _results/[yourjobname]/_. 
Here you will find _mites.fasta_ with all the MITEs sequences 
and also _mites.gff3_ with a gff file describing the MITEs in the genome file.

# Troubleshooting
If getting any error while running the BLASTn searches please check you blast+ version