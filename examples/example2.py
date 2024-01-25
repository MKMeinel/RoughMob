import numpy as np
from scipy.optimize import least_squares

# Example 2: Fit function example

# Define molecules list
molecules = ["DiMeBu", "TrMePe", "EtHx", "DiMeHx", "DIPDMP"]

# Define fit function
def fit_alpha(param, aVol, pVol):
    return 1 / param[0] * aVol + param[1] * pVol + param[2]

# Initialize arrays
Vshell = np.zeros(len(molecules))
Voverlap = np.zeros(len(molecules))
Vin = np.zeros(len(molecules))
Vout = np.zeros(len(molecules))
alpha = np.zeros(len(molecules))

# Loop over molecules to read data
for mm in range(len(molecules)):
    Vshell[mm], Voverlap[mm], Vin[mm], Vout[mm], _=np.genfromtxt("../Data/PureComponents/RoughnessVolumes/{}-RoughVolumes.txt".format(molecules[mm]))
    alpha[mm] = np.genfromtxt(f"../Data/PureComponents/Dynamics/{molecules[mm]}-DAA-DCG-alpha.txt", skip_header=2, max_rows=1)

# Determine parameters via least square fit
aV = Vshell + Voverlap  # active volumes
pV = np.exp(Vin / Vout)  # passive volumes

def fun(param):
    return fit_alpha(param, aV, pV) - alpha

param0 = [1, 1, 1]  # initial guess
res1 = least_squares(fun, param0)
[A, B, C] = res1.x  # fitted parameters
print([A, B, C])

alpha_pred = fit_alpha([A, B, C], aV, pV)


# Prediction of acceleration with Parameter as determined in Ref. Meinel 2022 

A=11.5069293  #Note: The convention for A here differs from the paper where A = 1/11.5 (instead of A = 11.5). 
B=0.83690501
C=-10.61236558
fitparas=np.array([A,B,C])


alpha_pred=fit_alpha(fitparas, aV,pV)















