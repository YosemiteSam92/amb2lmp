#!/usr/bin/env python3

import sys
import re
import numpy as np


class Amber:

    def __init__(self):

        # num of entries in prmtop or rst7 line (changes from section to
        # section)
        self.format = 0

        # section of prmtop or rst7 file, when fully read in before being parsed
        self.section = ""

        # list indices of distinct atom types in %FLAG AMBER_ATOM_TYPE section
        # useful to match atom types with masses and nonbonded coeffs
        self.distinctIndices = []

        # NTYPES parameter of the prmtop file
        # (in general NTYPES <= NATYP)
        self.ntypes = 0

        # amber atom type names for 1...natoms
        self.atomtNames = []

        # key: a unique entry of % FLAG ATOM_TYPE_INDEX
        # value: list of amber atom type names corresponding to key
        # E.g. {'1': ['cx', 'ca'], '2': ['h1'], '3': ['os'],...}
        self.index_to_names = {}

        # IAC:
        # {'cx': ['1', '1'], 'h1': ['2', '2'], 'os': ['3', '3'],
        # 'c3': ['4', '4'], 'ca': ['1', '1'], 'ha': ['5', '5'],
        # 'hc': ['6', '6'], 'nh': ['7', '7'], 'hn': ['8', '8'],
        # 'sy': ['9', '9'], 'o': ['10', '10']}
        # -->> IAC[i] = i

        # for i in 1...NTYPES, the ith entry of self.ICOindices is
        # NTYPES*(IAC(i)-1)+IAC(i)
        self.ICOindices = []

        # for i in 1...NTYPES, the ith entry of self.ICO is the entry of
        # %FLAG NONBONDED_PARM_INDEX indexed by self.ICOindices[i]
        self.ICO = []

        # for i in 1...NTYPES, the ith entry of self.Acoef is the A
        # Lennard-Jones coefficient for the homolytic interaction of the ith
        # amber type (as found in %FLAG ATOM_TYPE_INDEX) with itself.
        # Ditto for B.
        self.Acoef = []
        self.Bcoef = []

        self.nbondsH = 0
        self.nbondsNoH = 0
        self.nanglesH = 0
        self.nanglesNoH = 0
        self.ndihedralsH = 0
        self.ndihedralsNoH = 0

        # counting var
        self.nread = 0

    def setFormat(self, l, lammps):

        line = re.split('\(',l)
        line = re.split('[A-Za-z]', line[1])
        self.format = int(line[0])
        lammps.format = self.format

    # option to skip a blank character at the beginning of the line
    def readSection(self, l, target, skip):
        line = re.split("\s+", l)
        if skip == 1:
            line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            try:
                target[index] = line[index-self.nread]
            except IndexError:
                return 0
        self.nread += self.format

    def atomtIndex_to_atomtName(self, l):
        line = re.split("\s+", l)
        line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            try:
                atom_type_index = line[index-self.nread]
                if atom_type_index not in self.index_to_names:
                    self.index_to_names[atom_type_index] = \
                     [self.atomtNames[index]]
                else:
                    if self.atomtNames[index] not in \
                     self.index_to_names[atom_type_index]:
                        # print ("Appending: ", self.atomtNames[index], " to ",
                        #                     self.index_to_names[atom_type_index])
                        self.index_to_names[atom_type_index].append(
                         self.atomtNames[index])
            except IndexError:
                return 0
        self.nread += self.format

    def findICOindices(self):

        for i in range(1, self.ntypes+1):
            self.ICOindices.append(self.ntypes*(i-1)+i)

    def findICO(self, l):

        line = re.split("\s+", l)
        line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            if index+1 in self.ICOindices:
                nonbonded_parm_index = int(line[index-self.nread])
                self.ICO.append(nonbonded_parm_index)
                if nonbonded_parm_index < 0:
                    print("ERROR: detected negative NONBONDED_PARM_INDEX;")
                    print("this is currently not supported by this script")
                    sys.exit()
        self.nread += self.format

    def findAcoef(self,l):

        line = re.split("\s+", l)
        line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            if index+1 in self.ICO:
                self.Acoef.append(float(line[index-self.nread]))
        self.nread += self.format

    def findBcoef(self,l):

        line = re.split("\s+", l)
        line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            if index+1 in self.ICO:
                self.Bcoef.append(float(line[index-self.nread]))
        self.nread += self.format

    def printAtomtInfo(self):

        warning = 0
        print("\nThe following map associates to each Amber atom type number",
              "\n(as found in the ATOM_TYPE_INDEX section of the prmtop file)",
              "\nits corresponding set of Amber atom type names:")
        for index in self.index_to_names:
            print(index, self.index_to_names[index])
            if len(self.index_to_names[index]) > 1:
                warning = 1
        if warning:
            print("Warning: at least one Amber atom type number corresponds",
                  "to multiple atom type names.\n"
                  "Conversion to Lammps will issue a",
                  "new Lammps atom type number for each",
                  "Amber atom type name;\n",
                  "therefore Lammps atom type numbers",
                  "will not correspond exactly to Amber atom type numbers.\n")


class Lammps:

    # class variables: shared by all instances of Lammps
    nimpropers = 0  # treat impropers as dihedrals

    def __init__(self):

        # instance variables: unique to each instance of Lammps
        # must be preceded by self.

        # num of entries in prmtop or rst7 line (changes from section to
        # section)
        self.format = 0

        # num of atoms/bonds/angles/dihedrals (treat impropers as dihedrals)
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0

        # num of types
        self.natomt = 0
        self.nbondt = 0
        self.nanglet = 0
        self.ndihedralt = 0

        # box dimensions
        self.xlo = 0.
        self.xhi = 0.
        self.ylo = 0.
        self.yhi = 0.
        self.zlo = 0.
        self.zhi = 0.

        # atoms (for 1...natoms)
        self.charge = []        # partial atomic charges

        self.atomt = []  # lammps atom type num for all atoms
        self.atomtNames = []  # amber atom type names
        self.atomtNums = {}
        # key: amber atom type name, value: lammps atom type num

        # force field
        self.masses = []
        self.paircoefs = {}  # sigma, epsilon
        self.bondcoefs = [[], []]  # force const, equil length
        self.anglecoefs = [[], []]  # force const, equil length
        self.dihedralcoefs = [[], [], []]  # K,d,n

        # counting var
        self.nread = 0
        self.index = 0  # used only for inpcrd file

    def setNums(self, section, amber):
        line = re.split("\s+", section)
        # first entry is blank
        self.natoms = int(line[1])
        amber.ntypes = int(line[2])
        self.natomt = int(line[19])
        amber.nbondsH = int(line[3])
        amber.nbondsNoH = int(line[4])
        self.nbonds = int(line[3]) + int(line[4])
        amber.nanglesH = int(line[5])
        amber.nanglesNoH = int(line[6])
        self.nangles = int(line[5]) + int(line[6])
        amber.ndihedralsH = int(line[7])
        amber.ndihedralsNoH = int(line[8])
        self.ndihedrals = int(line[7]) + int(line[8])
        self.nbondt = int(line[16])
        self.nanglet = int(line[17])
        self.ndihedralt = int(line[18])

        # allocate memory

        # partial atomic charges rounded to 4 decimal places
        self.charge = np.zeros(self.natoms)

        # atom types
        self.atomt = np.zeros(self.natoms, dtype=int)
        amber.atomtNames = np.zeros(self.natoms, dtype='U2')
        self.x = np.zeros((self.natoms, 3))

        # bonds, angles, dihedrals
        # buffers should be sufficient even for very large systems
        self.bonds = np.zeros(self.nbonds, dtype='U100')
        self.angles = np.zeros(self.nangles, dtype='U100')
        self.dihedrals = np.zeros(self.ndihedrals, dtype='U100')

    # read into self.atomtNames only distinct atom type names
    def addAtomtNames(self, l, amber):

        line = re.split("\s+", l)
        for i in range(len(line)):
            if line[i] not in self.atomtNames:
                self.atomtNames.append(line[i])
                amber.distinctIndices.append(self.nread+i)
        self.nread += self.format

    def assignCharges(self, l):

        line = re.split("\s+", l)
        for i in range(self.format):
            try:
                # divide charges by 18.2223 to reverse conversion
                # into amber units; round division to 4 decimal
                # places to retain charge neutrality
                self.charge[self.nread+i] = \
                 '{:0.4f}'.format(float(line[i+1])/18.2223)
                # i+1 because first entry is blank
            except IndexError:
                return 0
        self.nread += self.format

    # read into self.masses only distinct values of atomic masses
    def addMasses(self, l, amber):

        line = re.split("\s+", l)
        line = line[1:]
        for index in range(self.nread, self.nread+self.format):
            if index in amber.distinctIndices:
                try:    # the last line can be shorter than self.format
                    self.masses.append(float(line[index-self.nread]))
                except IndexError:
                    return 0

        self.nread += self.format

    def addAtomtNums(self):

        num = 1
        for name in self.atomtNames:
            self.atomtNums[name] = num
            num += 1

    def assignAtomtNums(self, l):

        line = re.split("\s+", l)
        for index in range(self.nread, self.nread+self.format):
            try:    # the last line can be shorter than self.format
                self.atomt[index] = self.atomtNums[line[index-self.nread]]
            except IndexError:
                return 0
        self.nread += self.format

    def findPairCoef(self, amber):

        for i in range(amber.ntypes):
            A = amber.Acoef[i]
            B = amber.Bcoef[i]
            i = str(i+1)
            for atomtName in amber.index_to_names[i]:
                if A == 0 or B == 0:
                    epsilon = 0.
                    sigma = 0.
                else:
                    epsilon = B*B/(4*A)
                    sigma = (A/B)**(1/6)
                self.paircoefs[atomtName] = [epsilon, sigma]

    def addBondForceConst(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for forceConst in section:
            self.bondcoefs[0].append(float(forceConst))

    def addBondEqVal(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for equiVal in section:
            self.bondcoefs[1].append(float(equiVal))

    def addAngleForceConst(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for angleConst in section:
            self.anglecoefs[0].append(float(angleConst))

    def addAngleEqVal(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        twicePi = 2 * np.pi
        for equiVal in section:
            self.anglecoefs[1].append(360*float(equiVal)/twicePi)
            # convert to degrees

    def addDihedralBarrier(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for dihedralBarrier in section:
            self.dihedralcoefs[0].append(float(dihedralBarrier))

    def addDihedralPeriodicity(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for dihedralPeriodicity in section:
            self.dihedralcoefs[2].append(float(dihedralPeriodicity))

    def addDihedralD(self, section):

        section = re.split("\s+", section)
        section = section[1:]
        for dihedralPhase in section:
            dihedralPhase = float(dihedralPhase)
            if dihedralPhase < np.pi/2:         # phase = 0
                self.dihedralcoefs[1].append("1")
            else:                               # phase = pi
                self.dihedralcoefs[1].append("-1")

    def addBonds(self, section, begin):

        length = len(section)
        index2 = begin
        for index in range(0, length, 3):
            i = abs(int(section[index])) // 3 + 1  # integer division
            j = abs(int(section[index+1])) // 3 + 1
            btype = section[index+2]
            self.bonds[index2] = (str(index2+1)
                                  + " " + btype
                                  + " " + str(i)
                                  + " " + str(j))
            index2 += 1

    def addAngles(self, section, begin):

        length = len(section)
        index2 = begin
        for index in range(0, length, 4):
            i = abs(int(section[index])) // 3 + 1  # integer division
            j = abs(int(section[index+1])) // 3 + 1
            k = abs(int(section[index+2])) // 3 + 1
            btype = section[index+3]
            self.angles[index2] = (str(index2+1)
                                   + " " + btype
                                   + " " + str(i)
                                   + " " + str(j)
                                   + " " + str(k))
            index2 += 1

    def addDihedrals(self, section, begin):

        length = len(section)
        index2 = begin
        for index in range(0, length, 5):
            i = abs(int(section[index])) // 3 + 1  # integer division
            j = abs(int(section[index+1])) // 3 + 1
            k = abs(int(section[index+2])) // 3 + 1
            l = abs(int(section[index+3])) // 3 + 1
            btype = section[index+4]
            self.dihedrals[index2] = (str(index2+1)
                                      + " " + btype
                                      + " " + str(i)
                                      + " " + str(j)
                                      + " " + str(k)
                                      + " " + str(l))
            index2 += 1

    def assignCoord(self, l):

        # each line contains the coordinates of 2 atoms
        # if the line has been passed to this function, there is at least 1 atom
        line = re.split("\s+", l)
        line = line[1:]

        self.x[self.index][0] = line[0]
        self.x[self.index][1] = line[1]
        self.x[self.index][2] = line[2]
        self.setBoxBoundaries(float(line[0]), float(line[1]), float(line[2]))

        # but the second atom may be missing in the last line,
        # if the number of atoms is odd
        self.index += 1
        try:
            self.x[self.index][0] = line[3]
        except IndexError:
            return 0
        self.x[self.index][1] = line[4]
        self.x[self.index][2] = line[5]
        self.setBoxBoundaries(float(line[3]), float(line[4]), float(line[5]))

        self.index += 1

    def setBoxBoundaries(self, x, y, z):

        if x < self.xlo:
            self.xlo = x
        if x > self.xhi:
            self.xhi = x

        if y < self.ylo:
            self.ylo = y
        if y > self.yhi:
            self.yhi = y

        if z < self.zlo:
            self.zlo = z
        if z > self.zhi:
            self.zhi = z

    def printData(self, dataFileName):

        with open(dataFileName, "w") as df:
            print("LAMMPS data file printed by ambr2lmp.py\n", file=df)

            print(self.natoms, "atoms", file=df)
            print(self.nbonds, "bonds", file=df)
            print(self.nangles, "angles", file=df)
            print(self.ndihedrals, "dihedrals", file=df)
            print("0 impropers\n", file=df)

            print(self.natomt, "atom types", file=df)
            print(self.nbondt, "bond types", file=df)
            print(self.nanglet, "angle types", file=df)
            print(self.ndihedralt, "dihedral types\n", file=df)

            print(self.xlo, self.xhi, "xlo xhi", file=df)
            print(self.ylo, self.yhi, "ylo yhi", file=df)
            print(self.zlo, self.zhi, "zlo zhi\n", file=df)
            
            # rounding is meant to reflect the behavior of the lammps
            # write_data command as closely as possible
            
            print("Masses\n", file=df)
            atomtNum = 1
            for mass in self.masses:
                print(
                    atomtNum,
                    round(mass, 3),
                    file=df
                )
                atomtNum += 1

            print("\nPair Coeffs\n", file=df)
            for atomtNum in range(self.natomt):
                atomtName = self.atomtNames[atomtNum]
                print(
                    atomtNum + 1,
                    round(self.paircoefs[atomtName][0], 4),
                    round(self.paircoefs[atomtName][1], 5),
                    "#", atomtName, 
                    file=df
                )

            print("\nBond Coeffs\n", file=df)
            for bondtNum in range(self.nbondt):
                print(
                    bondtNum + 1,
                    round(self.bondcoefs[0][bondtNum], 1),
                    round(self.bondcoefs[1][bondtNum], 4),
                    file=df
                )

            print("\nAngle Coeffs\n", file=df)
            for angletNum in range(self.nanglet):
                print(
                    angletNum + 1,
                    round(self.anglecoefs[0][angletNum], 2),
                    round(self.anglecoefs[1][angletNum], 2),
                    file=df
                )

            print("\nDihedral Coeffs\n", file=df)
            # cast periodicity to int
            periodicity = np.array(self.dihedralcoefs[2])
            periodicity = periodicity.astype(int)
            for dihedraltNum in range(self.ndihedralt):
                print(
                    dihedraltNum + 1, 
                    round(self.dihedralcoefs[0][dihedraltNum], 6),
                    self.dihedralcoefs[1][dihedraltNum],
                    periodicity[dihedraltNum], 
                    file=df
                )

            print("\nAtoms\n", file=df)
            for atomNum in range(self.natoms):
                atomt = self.atomt[atomNum]
                print(atomNum + 1, "0",
                      atomt, self.charge[atomNum],
                      self.x[atomNum][0], self.x[atomNum][1],
                      self.x[atomNum][2], "#",
                      self.atomtNames[atomt-1], file=df)

            print("\nBonds\n", file=df)
            for bond in self.bonds:
                print(bond, file=df)

            print("\nAngles\n", file=df)
            for angle in self.angles:
                print(angle, file=df)

            print("\nDihedrals\n", file=df)
            for dihedrals in self.dihedrals:
                print(dihedrals, file=df)

            print("", file=df)


''' main program '''

if len(sys.argv) != 4:
    print("Error: wrong number of arguments.")
    print("Usage: amb2lmp .prmtop .inpcrd data.file\n",
          "\t.prmtop \tAmber parameter/topology file\n",
          "\t.inpcrd \tAmber coordinate file\n",
          "\tdata.file\tname of Lammps data file to output")
    sys.exit()

lammps = Lammps()
amber = Amber()

# parsing flags and vars
read = 0
np.set_printoptions(edgeitems=20)

# open prmtop file for reading
with open(sys.argv[1]) as prmtop:

    # find nums of atoms/bonds/angles/dihedrals
    # and names of all distinct atom types
    for l in prmtop:
        l = l.rstrip(' \n')  # remove trailing white spaces and EOL
        # print(l)
        # print(read)
        if read == 0 and l == "%FLAG POINTERS":
            next(prmtop)
            read = 1
            continue
        if read == 1:
            # read nums of atoms/bonds/angles/dihedrals
            if l[0] == "%":  # end of section
                lammps.setNums(amber.section, amber)
                amber.section = ""
                next(prmtop)
                read = 2
                continue
            amber.section += l
        if read == 2 and l == "%FLAG BOND_FORCE_CONSTANT":
            next(prmtop)
            read = 3
            continue
        if read == 3:
            if l[0] == "%":
                lammps.addBondForceConst(amber.section)
                amber.section = ""
                next(prmtop)
                read = 4
                continue
            amber.section += l
        if read == 4:
            if l[0] == "%":
                lammps.addBondEqVal(amber.section)
                amber.section = ""
                next(prmtop)
                read = 5
                continue
            amber.section += l
        if read == 5:
            if l[0] == "%":
                lammps.addAngleForceConst(amber.section)
                amber.section = ""
                next(prmtop)
                read = 6
                continue
            amber.section += l
        if read == 6:
            if l[0] == "%":
                lammps.addAngleEqVal(amber.section)
                amber.section = ""
                next(prmtop)
                read = 7
                continue
            amber.section += l
        if read == 7:
            if l[0] == "%":
                lammps.addDihedralBarrier(amber.section)
                amber.section = ""
                next(prmtop)
                read = 8
                continue
            amber.section += l
        if read == 8:
            if l[0] == "%":
                lammps.addDihedralPeriodicity(amber.section)
                amber.section = ""
                next(prmtop)
                read = 9
                continue
            amber.section += l
        if read == 9:
            if l[0] == "%":
                lammps.addDihedralD(amber.section)
                amber.section = ""
                next(prmtop)
                read = 10
                continue
            amber.section += l
        if read == 10 and l == "%FLAG AMBER_ATOM_TYPE":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 11
            continue
        if read == 11:
            if l[0] == "%":
                lammps.nread = 0
                read = 0
                break
            # read all distinct atom type names
            # their number is NATYP >= NTYPES
            if len(lammps.atomtNames) < lammps.natomt:
                lammps.addAtomtNames(l, amber)
            # store amber atom type names for all atoms
            amber.readSection(l, amber.atomtNames, 0)

    print("\nAmber atom type names:\n", lammps.atomtNames)

    # assign to each atom type name a unique num, starting from 1
    # this is the lammps atom type
    lammps.addAtomtNums()
    print("\nAssigned Lammps atom type numbers:\n", lammps.atomtNums)

    # read again %FLAG AMBER_ATOM_TYPE to assign to each atom its atom type num
    # based on the atom type name (this is not identical to reading the
    # %FLAG ATOM_TYPE_INDEX section, because NATYP >= NTYPES; in simple words,
    # atoms with different atom type names but with the same nonbonded coeffs
    # are counted as a single atom type index - E.g. ca and cx have
    # the same atom type index)

    # note: cannot use both tell() and next() at the same time

    prmtop.seek(0)
    next(prmtop)

    for l in prmtop:
        l = l.rstrip(' \n')
        # print(l)
        # print(read)
        if read == 1:
            if l[0] == "%":
                lammps.nread = 0
                read = 0
                break
            lammps.assignAtomtNums(l)
            continue
        if l == "%FLAG AMBER_ATOM_TYPE":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 1
            continue

    prmtop.seek(0)

    # read masses, charges and LJ coefs
    for l in prmtop:
        l = l.rstrip(' \n')  # remove trailing white spaces and EOL
        # print(l)
        # print(read)
        if read == 0 and l == "%FLAG CHARGE":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 1
            continue
        if read == 1:
            if l[0] == "%":
                lammps.nread = 0
                read = 2
                l = next(prmtop)
                continue
            lammps.assignCharges(l)
        if read == 2 and l == "%FLAG MASS":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 3
            continue
        if read == 3:
            lammps.addMasses(l, amber)
            if len(lammps.masses) == lammps.natomt:
                lammps.nread = 0
                read = 4
                continue
        if read == 4 and l == "%FLAG ATOM_TYPE_INDEX":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            amber.nread = 0
            read = 5
            continue
        if read == 5:
            if l[0] == "%":
                amber.printAtomtInfo()
                lammps.nread = 0
                amber.nread = 0
                read = 6
            if len(amber.index_to_names) < amber.ntypes:
                amber.atomtIndex_to_atomtName(l)
        if read == 6 and l == "%FLAG NONBONDED_PARM_INDEX":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 7
            amber.findICOindices()
            continue
        if read == 7:
            if l[0] == "%":
                amber.nread = 0
                read = 8
                continue
            amber.findICO(l)
        if read == 8 and l == "%FLAG LENNARD_JONES_ACOEF":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 9
            continue
        if read == 9:
            if l[0] == "%":
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 10
                continue
            amber.findAcoef(l)
        if read == 10:
            if l[0] == "%":
                amber.nread = 0
                read = 0
                lammps.findPairCoef(amber)
                break
            amber.findBcoef(l)

    # read bonds, angles and dihedrals
    prmtop.seek(0)

    for l in prmtop:
        l = l.rstrip(' \n')
        # print(l)
        # print(read)
        if l == "":
            continue
        if read == 0 and l == "%FLAG BONDS_INC_HYDROGEN":
            l = next(prmtop)
            amber.setFormat(l, lammps)
            read = 1
            amber.section = np.zeros(amber.nbondsH*3, 'U20')
            continue
        if read == 1:
            if l[0] == "%":
                lammps.addBonds(amber.section, 0)
                amber.section = np.zeros(amber.nbondsNoH*3, 'U20')
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 2
                continue
            amber.readSection(l, amber.section, 1)
        if read == 2:
            if l[0] == "%":
                lammps.addBonds(amber.section, amber.nbondsH)
                amber.section = np.zeros(amber.nanglesH*4, 'U20')
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 3
                continue
            amber.readSection(l, amber.section, 1)
        if read == 3:
            if l[0] == "%":
                lammps.addAngles(amber.section, 0)
                amber.section = np.zeros(amber.nanglesNoH*4, 'U20')
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 4
                continue
            amber.readSection(l, amber.section, 1)
        if read == 4:
            if l[0] == "%":
                lammps.addAngles(amber.section, amber.nanglesH)
                amber.section = np.zeros(amber.ndihedralsH*5, 'U20')
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 5
                continue
            amber.readSection(l, amber.section, 1)
        if read == 5:
            if l[0] == "%":
                lammps.addDihedrals(amber.section, 0)
                amber.section = np.zeros(amber.ndihedralsNoH*5, 'U20')
                amber.nread = 0
                l = next(prmtop)
                amber.setFormat(l, lammps)
                read = 6
                continue
            amber.readSection(l, amber.section, 1)
        if read == 6:
            if l[0] == "%":
                lammps.addDihedrals(amber.section, amber.ndihedralsH)
                amber.section = ""
                amber.nread = 0
                read = 0
                break
            amber.readSection(l, amber.section, 1)


# read coordinates

with open(sys.argv[2]) as inpcrd:

    inpcrd.readline()  # title
    inpcrd.readline()  # num of atoms

    for l in inpcrd:
        l = l.rstrip(' \n')
        lammps.assignCoord(l)

# print lammps data file

lammps.printData(sys.argv[3])



















