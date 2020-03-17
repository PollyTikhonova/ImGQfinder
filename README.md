![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
# ImGQfinder v2.0.0
A tool that searches G-quadruplexes and i-Motifs with buldges and mismatches.

## Installation
1. Clone the repository. 
2. Install the requirements. We recommend to create a new conda environment.

#### Via conda
```
conda env create -f environment.yml
conda activate imgqfinder_env  [windows]
source activate imgqfinder_env  [linux]
```
#### Via pip
```
pip install -r requirements.txt
```


## Usage
```python imgqfinder.py -i data.fasta -o test,```\
where ```data.fasta``` is a fasta file, which can contain several sequences [REQUIRED],\
    ```-o test``` is a fold, where the output files will be stored [NOT REQUIRED].
    
#### Full list of the options:
```
python imgqfinder.py --help                     
usage: ImGQFinder [-h] -i INPUT [-o OUTPUT] [-GC GC] [-L L] [-q Q]
                  [-mdef MDEF] [-bulgelen BULGELEN] [-maxbulge MAXBULGE] [-bp]
                  [-tetdef TETDEF] [-ns] [-r] [-v]

The tool for finding G-, C- quadruplexes. The output positions are represented
in a zero based counting.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Assembly scaffolds/contigs or full genomes, required.
  -o OUTPUT, --output OUTPUT
                        Name/path of a folder for output files. Saves to the
                        current folder if not provided.
  -GC GC                Quad type, G- or C-. By default, G.
  -L L                  Maximum loop length. By default, 7.
  -q Q                  Amount of tetrads.
  -mdef MDEF            Allowed number of defective tetrads. By default, 1.
  -bulgelen BULGELEN    Total length of bulges in one quadruplet. By default,
                        1.
  -maxbulge MAXBULGE    Maximum number of bulges per quadruplet. By default,
                        1.
  -bp, --bulge_priority
                        By default, quadrouplexes with shorter bulge or
                        without them are preferable while grouping. This
                        behaviour can be changed with this parameter.
  -tetdef TETDEF        Allowed number of defective nucleotides in tetrads. By
                        default, 1.
  -ns, --no-sequences   Not to include sequences to the output
  -r, --repeats         To include soft-masked genome areas. By default, not
                        included.
  -v, --verbose         Show the status of procesing or not. By default print
                        stages info
```
 
### Quadruplex Description
![Quadruplex Description](https://github.com/PollyTikhonova/ImGQfinder/raw/master/ImGQfinder_scheme.png)
