[![Licence](https://img.shields.io/badge/License-CC%20BY%20NC%20SA%204.0-grey.svg?style=for-the-badge)](http://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Lua](https://img.shields.io/badge/lua-%232C2D72.svg?style=for-the-badge&logo=lua&logoColor=white)](http://www.lua.org)

# CG CELLULOSE FIBRIL

Creates coarse grain model of cellulose fibrils using MARTINI3 beads.

## Author

Copyright (C) 2021-2022 Rodrigo Azevedo Moreira da Silva

[IPPT-PAN](http://www.ippt.pan.pl/staff/rams)
Instytut Podstawowych Problemów Techniki
Polskiej Akademii Nauk

## Reference

[![DOI](https://zenodo.org/badge/428765252.svg)](https://zenodo.org/badge/latestdoi/428765252)

The original paper describing the methodology 

> [Moreira RA, Weber SAL, Poma AB. Martini 3 Model of Cellulose Microfibrils: On the Route to Capture Large Conformational Changes of Polysaccharides. Molecules. 2022; 27(3):976.](https://doi.org/10.3390/molecules27030976)

## Usage
- Create a new folder with files of this repository.
- Download cellulose fibril from http://cces-sw.iqm.unicamp.br/cces/admin/cellulose/view;jsessionid= and copy the 'crystal.pdb' to the previously folder.
- Convert 'crystal.pdb' to gromacs format:
```
gmx editconf -f crystal.pdb -o fibril.gro
```
- Create the coarse grain model:
```
lua cgcellulose4PT2.lua
```
- Solvate the system:
```
gmx editconf -f cg-coordinates.gro -o cg-coordinates-newbox.gro -c -d 2.5 -bt triclinic
gmx solvate -radius 0.21 -cp cg-coordinates-newbox.gro -cs water-box-CG_303K-1bar.gro -o cg-coordinates-solvated.gro -p topol.top
```

## Should proceed with minimization, NVT, NPT and MD runs.

Some *examples* of GROMACS .mdp files are included.

- Always double-check the input/output files!
- We recommend to coarse-grain configuration from a previously all-atom simulation run.
- The model is been developed, all contributions are welcome!

## License

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](http://creativecommons.org/licenses/by-nc-sa/4.0/).

**Enjoy!**


