#!/usr/bin/python

import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer
import time
import glob
import numpy as np
import os

def main():
    Script, File_Name = sys.argv
    File_List = glob.glob('data.*')
    Mult = 20
    Mol_Temp_List = []
    Density = 0.0001
    Total_Mass = 0.0
    MW = 0.0
    i = 0
    PCBM_MW = 910.9 # g/mol
    for File in File_List:
        Mol_Temp_List.append(DA_Polymer.DA_Polymer(File))
        Total_Mass += Mol_Temp_List[i].MW
        MW = Mol_Temp_List[i].MW
        i += 1

    Avogadro = 6.0221413e23
    Comp_List = np.ones(i+1, dtype=int)
    Comp_List = Comp_List*Mult

    Total_Mols = float(i*Mult)/Avogadro
    Weight_Polymer = Total_Mols*MW
    Ratio = 1.5
    Weight_PCBM = Weight_Polymer*Ratio
    Mols_PCBM = Weight_PCBM/PCBM_MW
    Num_PCBM = 10000 #MANUAL ENTRY OF NUMBER OF ODCB
    print Num_PCBM, "PCBM Molecules"

    Name = "BHJ_" + os.getcwd().split('/')[-1]
    print Name

    Molar_Volume = MW/Density
    Volume = Molar_Volume*Total_Mols
    Box_Length_CM = Volume**(1./3.)
    Box_Length = Box_Length_CM*100000000

    PCBM = Molecule.Molecule(File_Name)
    PCBM.UnConverged = True
    PCBM.Set_Up_FF(run_orca=True, local = False)
    OPLS.Assign_OPLS(PCBM, ChelpG = False)
    Mol_Temp_List.append(PCBM)
    Comp_List[-1] = Num_PCBM

    BHJ = System.System(Mol_Temp_List, Comp_List, Box_Length, Name)
    BHJ.Gen_Rand()
    BHJ.Write_LAMMPS_Data()
    BHJ.Run_Lammps_Init(Nodes=2)
    BHJ.Run_Lammps_NPT(Nodes=2)
    BHJ.Run_Lammps_NPT(Nodes=2)
    System.Run_Glass_Transition(BHJ, 25, Ramp_Steps = 100000, Equil_Steps = 1000000, T_End = 100, Nodes=2)


if __name__=='__main__': main()
