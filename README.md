# ProteinAccessibilityCalculator (PAC)
This software performs 3 different tasks depending on the input parameters you provide when run it, specifically, depending on the parameter `CalculationType` that can take the values: `surface`, `pdb_surface`or `distance`..
 - If ***CalculationType=surface***, the software reads an input file (parameter **`input_file`**)  with a set of proteins and peptides, and calculates the surface accessibility of the atoms of the aminoacids stated in the parameter **`AAs`** (i.e. *K NZ*, that is, the atom *NZ* of the *Lysines*). To do that, it maps your proteins to [PDB](www.rcsb.org) entries (using Uniprot) and looks into all the possible [PDB](www.rcsb.org) structures, and using a JMol command, calculates and reports the surface accessibilities.
 - If ***CalculationType=pdb_surface***, the software takes a list of [PDB](www.rcsb.org) entries stated in parameter **`pdb_ids`** and calculates and reports the surface accessiblity of the atoms of the aminoacids stated in the parameter **`AAs`** (i.e. *K NE, R NE, D OD, E OG*, that is, the atoms *NE* of *K* and *R*, the atoms *OD* of *D* and the atom *OG* of *E*).
  - If ***CalculationType=distance***, the software reads an input file (parameter **`input_file`**)  with a set of proteins and peptides, and calculates the minimum distance between each of these atoms: OD1 and OD2 of D (Aspartic aCID), and atoms OE1 and OE2 of E (Glutamic Acid) and any other atom in the structure.
  
Based on Jmol, it calculates the surface accessibility of specific sites in a given protein

## How to get the software
You can get the latest version of the software at: http://sealion.scripps.edu/pac/

## How to run it
```
java -jar pac-calculator-1.0.jar input_parameters.properties
```
## Parameters files
The input parameter file is a simple text with a set of name-value pairs of parameters.

 - For ***surface*** running mode, look at example parameter file [here](https://raw.githubusercontent.com/proteomicsyates/ProteinAccessibilityCalculator/master/input_surface.properties).
 - For ***pdb_surface*** running mode, look at example parameter file [here](https://raw.githubusercontent.com/proteomicsyates/ProteinAccessibilityCalculator/master/input_surface_pdb.properties).
 - For ***distance*** running mode, look at example parameter file [here](https://raw.githubusercontent.com/proteomicsyates/ProteinAccessibilityCalculator/master/input_distances.properties).

## Getting help
For any customization of the program, report any issue or ask any question, please contact to **Salvador Mart?nez-Bartolome** (salvador at scripss.edu)



