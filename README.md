# GPSE
Genome-scaled Prediction of Substrate-specific Enzymes
The workflow:

![img](/img/GASSER.svg)

##### The philosophy:

In the post-genomic era, genome mining has been one of the most powerful tools leading to the discovery of new enzymes with properties of interest, and the combination of biological features like thermophilic and the sequence-based method like BLAST is the dominant approach. In this way, numerous enzymes have been found out and some of them have industrial applications. Also, with the development of structure-related computational methods involving docking and modeling, rational design and even de novo design of new enzymes have become a reality.

The accumulation of sequences ensures us unlimited biomineral and the development of structure-related computational methods endue us powerful tools to mine more effectively and completely.

##### Usage:  
Before run GPSE to a genome or a dataset of sequences, users should carefully modify the `GPSEcfg.py` file and change every path according to the system!  
```
python3 GPSE.py proteins.fasta ligand.pdb EC_number
```
