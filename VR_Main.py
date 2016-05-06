import Atom
import VR_Molecule
import Molecule
import OPLS
import System
import sys
import DA_Polymer

# run a basic lammps simulation to WIGGLE
def main():
    #tinker xyz. TODO basically parse tinker xyz and make it work like this
    Script, File_Name = sys.argv
    print File_Name
    Name = File_Name.split('.')[0]
    VR_mol = VR_Molecule.Molecule(File_Name)
    VR_mol.UnConverged = True
    OPLS.Assign_OPLS(VR_mol, ChelpG = False)
    VR_System = System.System([VR_mol], [1], 30.0, Name)
    VR_System.Gen_Rand()
    VR_System.Write_LAMMPS_Data()
    # run the fucking lammps simulation
    VR_System.Run_Lammps_Init_Local_Basic()
    # change lammps traj to my own fucking shit idk
    VR_mol.convertTraj()


if __name__=='__main__': main()
