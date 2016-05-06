#! usr/bin/python

# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
import random
import glob
# import Class structure
import Atom
import Bond
import Angle
import Dihedral
import Configure


class Molecule(object):
    """
    Class defining a molecule
    instance variables (self) :
        N = Number of Atoms
        Name = Molecule Name
        Atom_list[N] = List of Atom objects
        Bond_List[Num_Bonds] = list of bond objects
        Angle_List[Num_Angles] = list of angle objects
        Dihedral_List[Num_Dihedrals] = list of Dihedral objects
        Improper_List[Num_Impropers] = list of Improper objects
        COM = Center of Mass
        MOL_ID = ID Number for molecule
        UnConverged = Flag for bypassing Orca convergence (Default = False)
    """

    def __init__(self, File_Name):
        File = open(File_Name,'r') # File_Name is the name of an .xyz file outputted by Avogadro
        File_Lines = File.readlines()

        # Define Instance Variables
        self.Name = File_Name.split('.txyz')[0]
        self.N = int(File_Lines[0].strip('\n')) # Integer
        self.Atom_List = np.empty(self.N, dtype=object) # Numpy object array
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.MW = 0.0
        self.COM = np.zeros(3,float)
        self.Mol_ID = 0
        self.Missing_Dihedrals = 0
        self.UnConverged = False # Unconverged Orca Optimization

        print "Setting up molecule"
        print "Molecule Name:", self.Name
        print self.N, "Atoms in ", self.Name
        print "----------------------------------"
        for i in range(self.N):
            Line = File_Lines[1+i]
            Line = Line.strip('\n').split()
            Element = Line[1]
            Position = np.array( [ float(Line[2]), float(Line[3]), float(Line[4]) ], dtype=float )
            self.Atom_List[i] = Atom.Atom(Position, Element, i+1, Line[5]) # Instantiate Atom_List with Atom objects
            self.MW += self.Atom_List[i].Mass

        File.close()

        print "SET UP BONDS FOR EACH ATOM"
        File = open(File_Name,'r') # reopen file
        File_Lines = File.readlines()
        BondList = []

        for i in range(self.N):
            Line = File_Lines[1+i]
            Line = Line.strip('\n').split()
            j = 0
            print Line[6]
            for k in range(6,len(Line)):
                master = self.Atom_List[i]
                slave = self.Atom_List[int(Line[k])-1]
                if slave in master.Bond_List:
                    continue
                else:
                    master.Bond_List.append(slave)
                    slave.Bond_List.append(master)
                    self.Bond_List.append(Bond.Bond(master,slave,.0001))

        print "SET UP ANGLES FOR EACH ATOM"
        for i in range(len(self.Atom_List)):
            if len(self.Atom_List[i].Bond_List) > 1:
                for j in range (len(self.Atom_List[i].Bond_List)):
                    for k in range(j,len(self.Atom_List[i].Bond_List)):
                        atombonds = self.Atom_List[i].Bond_List
                        if atombonds[k] != atombonds[j]:
                            self.Angle_List.append(Angle.Angle(self.Atom_List[i],atombonds[j],atombonds[k]))

        print "SET UP DIHEDRALS FOR EACH ATOM"


        print "Initial XYZ Coordinates:\n"
        for Atom_Obj in self.Atom_List:
            print Atom_Obj.Element, Atom_Obj.Position
        print "----------------------------------"

        # Compute center of Mass
        Mass_Weighted_Sum = np.zeros(3,dtype=float)
        for Atom_Obj in self.Atom_List:
            Mass_Weighted_Sum += Atom_Obj.Position*Atom_Obj.Mass
            self.MW += Atom_Obj.Mass

        print "Molecular Weight is ", self.MW, "Grams/Mole"
        self.COM = Mass_Weighted_Sum/self.MW
        self.COM = np.asarray(self.COM)
        print "COM is", self.COM

        Mass_Weighted_Sum = 0.0
        for Atom_Obj in self.Atom_List:
            Mass_Weighted_Sum += Atom_Obj.Mass*((Atom_Obj.Position[0] - self.COM[0])**2 + (Atom_Obj.Position[1] - self.COM[1])**2  + (Atom_Obj.Position[2] - self.COM[2])**2 )
        Rg2 = Mass_Weighted_Sum/self.MW
        self.Rg = np.sqrt(Rg2)
        print "The Radius of Gyration is", self.Rg

        # Zero COM
        for Atom_Obj in self.Atom_List:
            Atom_Obj.Position -= self.COM

        self.COM -= self.COM
        print self.COM
        return

    def Adjust_COM(self):
        # This adjusts the center of mass and gives the molecule a random orientation
        x = random.random()*2*3.1415
        y = random.random()*2*3.1415

        C1 = np.cos(x)
        S1 = np.sin(x)
        C2 = np.cos(y)
        S2 = np.cos(y)

        for Atom_Obj in self.Atom_List:
            # First rotation
            xt = Atom_Obj.Position[0]
            yt = Atom_Obj.Position[1]
            Atom_Obj.Position[0] = xt*C1 - yt*S1
            Atom_Obj.Position[1] = xt*S1 + yt*C1
            # Second rotation
            #xt = Atom_Obj.Position[0]
            #zt = Atom_Obj.Position[2]
            #Atom_Obj.Position[0] = xt*C2 - zt*S2
            #Atom_Obj.Position[2] = xt*S2 + zt*C2

            Atom_Obj.Position += self.COM
        return

    def convertTraj(self):
        openName = "Init_%s.lammpstrj" % self.Name
        closeName = "Init_%s.txyz" % self.Name
        closeFile = open(closeName,'w')
        numAtoms = 0
        with open(openName,'r') as openFile:
            line = openFile.readlines()
            for i in range(len(line)):
                currLine = line[i].strip('\n')
                if currLine == "ITEM: TIMESTEP":
                    nextLine = line[i+1].strip('\n')
                    closeFile.write("timestep %s\n" % nextLine)
                elif currLine == "ITEM: NUMBER OF ATOMS":
                    nextLine = line[i+1].strip('\n')
                    closeFile.write("%s\n" % nextLine)
                    numAtoms = int(nextLine)
                if currLine == "ITEM: ATOMS id type mol xs ys zs vx vy vz ":
                    for j in range(numAtoms):
                        currLine = line[i+j+1].strip('\n').split()
                        atom = findAtomByID(self.Atom_List,int(currLine[0]))
                        closeFile.write("%s %s %s %s %s %s " % (currLine[0], atom.Element, currLine[3], currLine[4], currLine[5], atom.txyz)) #write positional data
                        #write bonds next
                        for k in range(len(atom.Bond_List)):
                            closeFile.write("%s " % atom.Bond_List[k].Atom_ID)
                        closeFile.write('\n')
        openFile.close()

def findAtomByID(a_list, a_id):
    for i in range(len(a_list)):
        if a_list[i].Atom_ID == a_id:
            return a_list[i]

# Functions Operating on sets of Molecule objects

def Assign_Lammps(Moltemp_List):
    """
        Function that inputs a list of molecule templates. It searches through all the atoms, bonds, angles etc. to find the unique types of interactions
        present in an arbitrary system object made of molecule templates.

        returns a list of unique params for writing to a LAMMPS data file. these are lists defined such that the i-1 element corresponds to the ith LAMMPS_Type
    """

    print "Finding unique Atoms"
    Unique_Atoms = []
    Atom_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Atom_Obj in Moltemp_Obj.Atom_List:
            i = 1
            for Type in Unique_Atoms:
                if Atom_Obj.OPLS_Type != Type:
                    i += 1
                if Atom_Obj.OPLS_Type == Type:
                    Atom_Obj.LAMMPS_Type = i
            if i > len(Unique_Atoms):
                Unique_Atoms.append(Atom_Obj.OPLS_Type)
                Atom_Params.append([Atom_Obj.Mass, Atom_Obj.Sigma, Atom_Obj.Epsilon])
                Atom_Obj.LAMMPS_Type = i
                print "Found new atom type:", i, "OPLS_ID =", Atom_Obj.OPLS_Type, Atom_Params[i-1]
                print "Mass is", Atom_Obj.Mass, "Element is", Atom_Obj.Element

    print "Finding unique Bonds"
    Bond_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Bond_Obj in Moltemp_Obj.Bond_List:
            i = 1
            for Params in Bond_Params:
                if Params[0] != Bond_Obj.kb or Params[1] != Bond_Obj.req:
                    i += 1
                if Params[0] == Bond_Obj.kb and Params[1] == Bond_Obj.req:
                    Bond_Obj.LAMMPS_Type = i
            if i > len(Bond_Params):
                Bond_Params.append([Bond_Obj.kb, Bond_Obj.req])
                Bond_Obj.LAMMPS_Type = i
                print "Found new bond type:", i, Bond_Obj.kb, Bond_Obj.req

    print "Finding unique Angles"
    Angle_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Angle_Obj in Moltemp_Obj.Angle_List:
            i = 1
            for Params in Angle_Params:
                if Params[0] != Angle_Obj.ka or Params[1] != Angle_Obj.Angle_Eq:
                    i += 1
                if Params[0] == Angle_Obj.ka and Params[1] == Angle_Obj.Angle_Eq:
                    Angle_Obj.LAMMPS_Type = i
            if i > len(Angle_Params):
                Angle_Params.append([Angle_Obj.ka, Angle_Obj.Angle_Eq])
                Angle_Obj.LAMMPS_Type = i
                print "Found new angle type:", i, Angle_Params[i-1]


    print "Finding unique Dihedrals"
    Dihedral_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Dihedral_Obj in Moltemp_Obj.Dihedral_List:
            i = 1
            for Params in Dihedral_Params:
                if Dihedral_Obj.Coeffs[0] != Params[0] or Dihedral_Obj.Coeffs[1] != Params[1] or Dihedral_Obj.Coeffs[2] != Params[2] or Dihedral_Obj.Coeffs[3] != Params[3]:
                    i += 1
                if Dihedral_Obj.Coeffs[0] == Params[0] and Dihedral_Obj.Coeffs[1] == Params[1] and Dihedral_Obj.Coeffs[2] == Params[2] and Dihedral_Obj.Coeffs[3] == Params[3]:
                    Dihedral_Obj.LAMMPS_Type = i
            if i > len(Dihedral_Params):
                Dihedral_Params.append(Dihedral_Obj.Coeffs)
                Dihedral_Obj.LAMMPS_Type = i
                print "Found new dihedral type:", i, Dihedral_Params[i-1]


    print "Finding unique Impropers"
    Improper_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Improper_Obj in Moltemp_Obj.Improper_List:
            i = 1
            for Params in Improper_Params:
                if Improper_Obj.Ki != Params[0] or Improper_Obj.Improper_Eq != Params[1]:
                    i += 1
                if Improper_Obj.Ki == Params[0] and Improper_Obj.Improper_Eq == Params[1]:
                    Improper_Obj.LAMMPS_Type = i
        try:
            if i > len(Improper_Params):
                Improper_Params.append([Improper_Obj.Ki, Improper_Obj.Improper_Eq])
                Improper_Obj.LAMMPS_Type = i
                print "Found improper type:", i, Improper_Obj.Ki, Improper_Obj.Improper_Eq
        except:
            continue



    return Atom_Params, Bond_Params, Angle_Params, Dihedral_Params, Improper_Params
