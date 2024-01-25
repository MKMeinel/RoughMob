from RManalysis import *

# Example 1: Calculation of molecular roughness difference and roughness volumes

# Define file paths
CGpot = "../Data/PureComponents/CGPotentials/DiMeBu-nb.A-A.pot.table"
trajAA = "../Data/PureComponents/LammpsTrajAndData/DiMeBu-AA-traj.dump"
dataAA = "../Data/PureComponents/LammpsTrajAndData/DiMeBu-AA-system.data"
trajCG = "../Data/PureComponents/LammpsTrajAndData/DiMeBu-CG-traj.dump"

# Calculate molecular roughness difference (Ra) for single bead molecule
Ra, rAA = single_bead_roughness(trajAA, dataAA, CGpot, mols_n=2)
Ra = np.mean(Ra)
rAA = np.mean(rAA) # average atomistic radius
print(Ra)

# Calculate radii for roughness volume calculations
Rout = rAA + Ra
Rin = rAA - Ra

# Calculation of passive volumes
Vin = 4 / 3 * np.pi * Rin**3
Vout, Vmol = calc_Vout_AA(trajAA, dataAA, Rout) #from AA trajectory
print(Vout)
Vout, Vmol = calc_Vout_CG(trajCG, Rout, frame=5) #from CG trajectory, best: loop over several frames
print(Vout)

# Calculation of active volumes
Vshell = 4 / 3 * np.pi * (Rout**3 - Rin**3)
Voverlap = 4 / 3 * np.pi * Rout**3 - (Vmol - Vout)






