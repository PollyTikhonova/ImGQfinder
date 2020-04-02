![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
# ImGQfinder v2.2.0
A tool that searches G-quadruplexes and i-Motifs with buldges and mismatches.

#### Updates tracker
You can find the information about what's new in updates and versions of libraries via [link.](https://github.com/PollyTikhonova/ImGQfinder/blob/master/version_tracker.md)

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
```   -o test``` is a fold, where the output files will be stored [NOT REQUIRED].
    
At the output folder, by default, you will get the following files:
 - `%fasta_id%_quadruplets.csv`: the full list of quadrupletes;
 - `%fasta_id%_quadruplexes.csv`: the full list of quadruplexes;
 - `%fasta_id%_groups.csv`: the non-redundant list of quadruplexes: without any intersections;
 - `%fasta_id%_ranges.bed`: only start&end coordinates of the grouped quadruplexes in a BED-FORMAT;
 - `description.txt`: the columns description file.
 
 At this repository, you can find the test folder with the input and output example.
 
 ##### ! Important !
While grouping, quadruplexes are not merging with the nearest, but the best possible quadruplex is chosen according to the next considerations. \
If two quadruplexes intersect we prefer the one that: *(the conditions are listed in their priority order)*
 - has less missmatches and buldges (by default, but you can abort this behavior);
 - has less total length (this means the the space between quadruplets is less);
 - meets first.
 
#### Big Files & Multiprocessing
By default, the output will contain the quadruplex sequences. In case of big fasta sequences, you will need significantly more **time** and **space on the disk**. To avoid overconsuming of the resources you can do the following:\
- If you do not need sequences you can turn off this behaviour with the tag `--no-sequences`. 
- Also, you can reduce the time of computations by multiprocessing option. Just type the number of kernels the program may use, like `--nthreads 4`. *Warning: there could be some issues with multiprocessing in Windows systems.*
- By default, the program will generate 4 files, containing information about: quadruplets, quadruplexes, groups and ranges. You may request not all the files as output. Just type one or several names with the tag: 

    
### Quadruplex Description
![Quadruplex Description](https://github.com/PollyTikhonova/ImGQfinder/raw/master/ImGQfinder_scheme.png)
 

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
  -q Q                  The length of a quadruplet.
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
  -ns, --no-sequences   Not to include sequences to the output.
  -r, --repeats         To include soft-masked genome areas. By default, not
                        included.
  -v, --verbose         Show the status of procesing or not. By default print
                        stages info
  --nthreads NTHREADS   Number of kernels to use.
```
