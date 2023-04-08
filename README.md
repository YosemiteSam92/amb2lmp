# amb2lmp.py

A python script to convert Amber topology (.prmtop) and coordinates (.rst7) files into a LAMMPS data file.

## The reason behind this script

The amber2lammps.py script distributed with LAMMPS (Dec. 2018) groups Amber atom types with the same LJ coefficients into a single LAMMPS atom type. This is undesirable if later one needs to keep track of the Amber types.

This script assigns a new LAMMPS atom type to each Amber atom type, at the (innocuous) cost of introducing redundancy in pair coefficients. If this happens, an explicative warning message is printed to the screen. 

Furthermore, this script is considerably faster than the one shipping with LAMMPS when the number of atoms exceeds 10^5 (and possibly even a lower threshold). 

## Technical notes

Amber atom types will be printed in the 'Pair Coeffs' and 'Atoms' sections of the data file, preceded by an '#', so that the Lammps parser can recognize them as comments. 

Improper torsional terms will be treated as standard dihedrals. 

## Usage: 

python3 amb2lmp.py .prmtop .rst7 data.outputFile

## Structure

- amb2lmp.py

    The conversion scripts itself.

- conversionRules.txt

    Explanations of a few critical steps in the Amber-to-Lammps conversion.

- example

    A pair of Amber topology and coordinates files, as well as the Lammps data file output by the script, for example purposes. 

- prevVersions

    Modified versions of the amber2lammps.py script shipped with Lammps Dec. 18. They are kept just for reference and their usage is not recommended. 