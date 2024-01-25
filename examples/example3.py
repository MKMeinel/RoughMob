import numpy as np

# Example 3: Prediction of acceleration in a DiMeBu-DiMeHx mixture

# Defining the components of the mixture
comp1 = "DiMeBu"
comp2 = "DiMeHx"

# Loading roughness volumes for both components from text files
Vshell1, Voverlap1, Vin1, Vout1, Vmol1 = np.genfromtxt(f"../Data/PureComponents/RoughnessVolumes/{comp1}-RoughVolumes.txt")
Vshell2, Voverlap2, Vin2, Vout2, Vmol2 = np.genfromtxt(f"../Data/PureComponents/RoughnessVolumes/{comp2}-RoughVolumes.txt")

# Defining mole fraction array
x1 = np.array([0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.98, 0.99, 1])
x2 = 1 - x1

# Fit parameters based on the paper (with noted convention difference)
A, B, C = 11.5069293, 0.83690501, -10.61236558
fitparas = np.array([A, B, C])

# Calculate overall acceleration in the system
Voutmix_ideal = x1 * Vout1 + x2 * Vout2
Vinmix_ideal = x1 * Vin1 + x2 * Vin2
Vshellmix_ideal = x1 * Vshell1 + x2 * Vshell2
Voverlapmix_ideal = x1 * Voverlap1 + x2 * Voverlap2
Vmolmix_ideal = x1 * Vmol1 + x2 * Vmol2

# Predicting ideal mixture acceleration
alphamix_ideal_pred = 1 / A * (Vshellmix_ideal + Voverlapmix_ideal) + B * np.exp(Vinmix_ideal / Voutmix_ideal) + C

# Constants for calculations
omega = 0.5
Voutmix_equi = 0.5 * Vout1 + 0.5 * Vout2

# Calculating passive and active correction terms
P12 = omega * B * np.sqrt(np.exp(Vin1 / Voutmix_equi) * np.exp(Vin2 / Voutmix_equi))
A12 = 0.5 * 1 / A * (Vshell1 + Voverlap1) + 0.5 * 1 / A * (Vshell2 + Voverlap2)

# Predicting mixture acceleration
alphamix_pred = alphamix_ideal_pred + x1 * x2 * (1 - P12) * A12

# Prediction of self-diffusion acceleration
F1 = (Vshell2 / Vshell1) * (Vmol1 / Vmol2)**(2 / 3)
F2 = 1 / F1
alpha1_self_pred = alphamix_pred / (x1 + x2 * F1)
alpha2_self_pred = alphamix_pred / (x1 * F2 + x2)

# Comparison to actual alpha as determined via simulation
dynamics = "../Data/Mixtures/Dynamics/DiMeBu-DiMeHx-alpha.txt"
data = np.genfromtxt(dynamics)
x, alphamix, alpha1_self, alpha2_self, alphabinary = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]

# Error calculation
print("Alpha Error-self diffusion coefficient")
print(np.mean(np.abs(alpha1_self_pred[1:-1] - alpha1_self[1:-1]) / alpha1_self[1:-1]) * 100)
print(np.mean(np.abs(alpha2_self_pred[1:-1] - alpha2_self[1:-1]) / alpha2_self[1:-1]) * 100)

# Prediction of binary diffusion acceleration
# Without dynamics of any AA simulations
alphabinary_pred = alphamix_pred[1:-1] / (x2[1:-1] * F1 + x1[1:-1] * 1 / F1)

# With known acceleration factors of pure components and at equimolar mixture
alphamix_ideal_pred_bin = x1[1:-1] * alpha1_self[-1] + x2[1:-1] * alpha2_self[0]
# Index for equimolar mixture: 6 and 5 (with and without x=0 and 1)
alpha_bin_equi = alphabinary[6]
alpha_ideal_pred_equi = alphamix_ideal_pred_bin[5]

# Calculate F1 and omega from equimolar simulations
F1_equi = alpha2_self[6] / alpha1_self[6]
omega_bin = (F1_equi * (4 * alpha_ideal_pred_equi + A12) + alpha_bin_equi * (-2 * F1_equi**2 - 2)) / (A12 * F1_equi * P12)

alphabinary_pred = (alphamix_ideal_pred_bin + x1[1:-1] * x2[1:-1] * (1 - omega_bin * P12) * A12) / (x2[1:-1] * F1_equi + x1[1:-1] * 1 / F1_equi)

