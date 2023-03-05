## Usage:  
Before run GPSE to a genome or a dataset of sequences, users should carefully modify the `GPSEcfg.py` file and change every path according to the system!  
```
python3 GPSE.py proteins.fasta ligand.pdb EC_number

GPSE: Genome-scale Prediction of Substrate-specific Enzymes

positional arguments:
  input           input a seqence in fasta format and we highly suggest that
                  you make the head linesimple
  ligand          ligand or substrate in pdb
  EC_number       the EC_number in form of x_x_x (e.g. 3_1_1 for 3.1.1.
                  carboxylic-ester hydrolase)

optional arguments:
  -h, --help      show this help message and exit
  -nos NOS        Max Number Of Structures you wish to be predicted for each
                  sequence, defult is 1
  -nt NT          Number of Threads used to run BLAST, default is 8
  -evalue EVALUE  E-Value cutoff for structure template, default is 1e-5
  -ic IC          the Identity Cutoff for structure prediction, default is 30
```
## Reference:  
Sun J, Xia Y and Ming D (2020) Whole-Genome Sequencing and Bioinformatics Analysis of Apiotrichum mycotoxinivorans: Predicting Putative Zearalenone-Degradation Enzymes. Front. Microbiol. 11:1866. doi: 10.3389/fmicb.2020.01866

