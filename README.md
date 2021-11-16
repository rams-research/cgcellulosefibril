# cgmodel

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
- Choose MARTINI3 bead type and bonded interation strength
```
sed -e "s/ATX1/SP6/g" -e "s/999999/250/g" -i topol.top
```
## Should proceed with minimization, NVT, NPT and MD runs.
