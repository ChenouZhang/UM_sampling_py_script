import numpy as np
import matplotlib.pyplot as plt
import math
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.rms import RMSD

#Python scripts to ABMD traj.
#Help to plot ddRMSD

################ Load Data ###################
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

################ Initialize parameters ################
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
    Ini_position = Ini_temp.positions
    End_position = End_temp.positions
    align.alignto(Ini_temp,End_temp)
    rmsd_ab = rmsd(Ini_temp,End_temp)
    return rmsd_ab,Ini_position,End_position

ddRMSD = []                        # ddRMSD

################ Functions ####################
#RMSD_ini = RMSD(ABMD,Ini_ref,select='name CA')
#RMSD_ini.run()
#RMSD_end = RMSD(ABMD,End_ref,select='name CA')
#RMSD_end.run()
#rmsd_ini = RMSD_ini.rmsd.T
#rmsd_end = RMSD_end.rmsd.T

#time = rmsd_ini[1]
#ddRMSD = (rmsd_ini[2]-rmsd_end[2])/rmsd_ab


