import numpy as np
import matplotlib.pyplot as plt
import math
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.rms import RMSD
from os import path

#Python scripts to ABMD traj.
#Help to plot ddRMSD

################ Load Traj ###################
def Load_Traj(PDB_Ini,PDB_End,XTC=None):
    '''
    Parameters: 
    	PDB_Ini: Input Ini_structure file name
    	PDB_End: Input End_structure file name
        XTC: Traj file comes with PDB_Ini   
	Returns:
        ABMD: Full traj
    	Ini_ref: Initial frame
        End_ref: Target frame
    '''
    ABMD = mda.Universe(PDB_Ini,XTC)
    End_ref = mda.Universe(PDB_End)
    for ts in ABMD.trajectory[0]:
        Ini_ref = ABMD.atoms

    return ABMD,Ini_ref,End_ref

################ Calculate RMSD ################
def Calculate_RMSD(Ini_ref,End_ref):
    '''
    Parameters: 
	    Ini_ref: Initial frame
        End_ref: Target frame
    Returns:
        rmsd_ab: rmsd between start and end frame
	    Ini_position: Coord of initial frame
        End_position: Coord of target frame
    '''
    Ini_temp = Ini_ref.select_atoms('name CA')
    End_temp = End_ref.select_atoms('name CA')
    align.alignto(Ini_temp,End_temp)
    Ini_position = Ini_temp.positions
    End_position = End_temp.positions
    rmsd_ab = rmsd(Ini_position,End_position)
    return rmsd_ab,Ini_position,End_position

def Calculate_ddRMSD(ABMD,Ini_ref,End_ref,rmsd_ab):
    '''
    Parameters: 
        ABMD: Full traj
	    Ini_ref: Initial frame
        End_ref: Target frame
        rmsd_ab: rmsd between start and end frame
    Returns:
        time: time coordinate for one traj
	    ddRMSD: ddRMSD along time
    '''
    RMSD_ini = RMSD(ABMD,Ini_ref,select='name CA')
    RMSD_ini.run()
    RMSD_end = RMSD(ABMD,End_ref,select='name CA')
    RMSD_end.run()
    rmsd_ini = RMSD_ini.rmsd.T
    rmsd_end = RMSD_end.rmsd.T
    time = rmsd_ini[1]
    ddRMSD = (rmsd_ini[2]-rmsd_end[2])/rmsd_ab
    return time,ddRMSD


################ Main Function ####################
def Multi_ddRMSD(Num_Digit = 3, Ini_Frame_Name = 'step5_assembly.psf', Traj_Name = 'md_plumed_ld.part0001.xtc', End_Frame_Name = 'occ_prot_plumed_target_HN.pdb'):
    index = 1
    plt.figure()
    while (path.exists('./'+'0'*(Num_Digit - len(str(index))) + str(index))):
        Sub_Folder_Name = './'+'0'*(Num_Digit - len(str(index))) + str(index) + '/'
        PDB_Ini = Ini_Frame_Name
        PDB_End = Sub_Folder_Name + End_Frame_Name
        XTC = Sub_Folder_Name + Traj_Name
        ABMD, Ini_ref, End_ref = Load_Traj(PDB_Ini,PDB_End,XTC)
        rmsd_ab, Ini_position, End_position = Calculate_RMSD(Ini_ref,End_ref)
        time, ddRMSD = Calculate_ddRMSD(ABMD, Ini_ref, End_ref, rmsd_ab)
        plt.plot(time,ddRMSD)
        index = index + 1
    plt.show()
    return 0


