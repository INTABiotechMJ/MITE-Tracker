## About

IRmatcher is efficient and easy to run tool for discovering inverted repeats in genomics sequences. It is written in python and uses ncbi's blast+ for matching sequences. 

## Requirements (just follow how to install)
 - tested in macOS 10.13.1, Debian 7.6, Ubuntu 16.04
 - ncbi blast+ (Nucleotide-Nucleotide BLAST 2.6.0+)
 - biopython

## How to install
```
sudo apt-get install ncbi-blast+ virtualenv 
#for macOS user 
#brew install ncbi-blast+ virtualenv
virtualenv venv
source venv/bin/activate
pip install -r requirementes.txt
```

## How to run
```
python miteParser.py -g /path/to/your/genome.fasta -j jobname
```
Or specify -w if you're lucky and have 3 free CPUs
```
python miteParser.py -g /path/to/your/genome.fasta -w 3 -j jobname
```

## How to run and leave running in background with nohup
```
nohup python -u miteParser.py -g /path/to/your/genome.fasta -w 3 -j jobname &
```

In order to check the job you can use this command (ctrl+c to exit)
```
tail -f nohup.out
```
## Command line options
| Argument  | Description | Data type  | Required or default |
| ------------- | ------------- | ------------- | ------------- |
| -g  | Genome file in fasta format  | string  | required  |
| -j  | Jobname. Result files will be created in results/jobname   | string  | required  |
| -w  | Max number of processes to use simultaneously  | int  | 1  |
| -tsd_min_len  | TSD min lenght  | int  | 2  |
| -tsd_max_len  | TSD max lenght  | int  | 10  |
| -mite_min_len  | MITE min lenght  | int  | 50  |
| -mite_max_len  | MITE max lenght  | int  | 650  |



## Results
All the results are placed in _results/[yourjobname]/_. 
Here you will find _mites.fasta_ with all the MITEs sequences 
and also _mites.gff3_ with a gff file describing the MITEs in the genome file.