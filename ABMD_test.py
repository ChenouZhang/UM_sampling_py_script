import numpy as np
import matplotlib.pyplot as plt
import math
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.rms import RMSD
from os import path
import pandas as pd

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

def Calculate_ddRMSD(ABMD,Ini_ref,End_ref):
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
    Ini_temp = Ini_ref.select_atoms('name CA')
    End_temp = End_ref.select_atoms('name CA')
    align.alignto(Ini_temp,End_temp)
    Ini_position = Ini_temp.positions
    End_position = End_temp.positions
    rmsd_ab = rmsd(Ini_position,End_position)
    RMSD_ini = RMSD(ABMD,Ini_ref,select='name CA')
    RMSD_ini.run()
    RMSD_end = RMSD(ABMD,End_ref,select='name CA')
    RMSD_end.run()
    rmsd_ini = RMSD_ini.rmsd.T
    rmsd_end = RMSD_end.rmsd.T
    time = rmsd_ini[1]
    ddRMSD = (rmsd_ini[2]-rmsd_end[2])/rmsd_ab
    return time,ddRMSD,rmsd_ini,rmsd_end


################ Main Function ####################
def Multi_ddRMSD(
    ddRMSD_dict, 
    Num_Digit = 3, 
    Ini_Frame_Name = 'step5_assembly.psf', 
    Traj_Name = 'md_plumed_ld.part0001.xtc', 
    End_Frame_Name = 'occ_prot_plumed_target_HN.pdb'):
    '''
    Parameters: 
        ddRMSD_dict: Dictionary to store time coordinate and ddRMSD data for each folder
        Num_Digit: The total length of the folder name(folder name should be a pure number otherwise you need to modify this).
	    Ini_Frame_Name: 
        Traj_Name:
        End_Frame_Name:
    Returns:
        ddRMSD_dict: all coordiantes data stores in this dictionary.
    '''
    index = 1
    while (path.exists('./'+'0'*(Num_Digit - len(str(index))) + str(index))):
        Sub_Folder_Name = './'+'0'*(Num_Digit - len(str(index))) + str(index) + '/'
        PDB_Ini = Ini_Frame_Name
        PDB_End = Sub_Folder_Name + End_Frame_Name
        XTC = Sub_Folder_Name + Traj_Name
        print('Enter'+Sub_Folder_Name)
        ABMD, Ini_ref, End_ref = Load_Traj(PDB_Ini,PDB_End,XTC)
        print('Load_Traj done.')
        rmsd_ab, Ini_position, End_position = Calculate_RMSD(Ini_ref,End_ref)
        print('Calculate_RMSD done.')
        time, ddRMSD, rmsd_ini, rmsd_end = Calculate_ddRMSD(ABMD, Ini_ref, End_ref)
        print('Calculate_ddRMSD')
        ddRMSD_dict[Sub_Folder_Name+ '_time'] = time
        ddRMSD_dict[Sub_Folder_Name + '_ddRMSD'] = ddRMSD
        index = index + 1
        print('Now moving to next folder.')
    return 0

def Print_ddRMSD(ddRMSD_dict, Num_Digit = 3):
    index = 1
    plt.plot(ddRMSD_dict['./'+'0'*(Num_Digit - len(str(index))) + str(index) + '/_time'][:-1],ddRMSD_dict['./'+'0'*(Num_Digit - len(str(index))) + str(index) + '/_ddRMSD'][:-1])
    index = index + 1
    plt.show()
    return
    

