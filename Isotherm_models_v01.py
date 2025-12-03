import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

# 
# A full-range moisture sorption model for cellulose-based materials
# yielding consistent net isosteric heat of sorption
#
# Johan Tryding (a,b), Henrik Askfelt (b), Marcus Alexandersson (b), and
# Matti Ristinmaa (a)
# (a) Division of Solid Mechanics, Lund University, SE-221 00 Lund, Sweden ; 
# (b) Tetra Pak, Ruben Rausingsgata, SE-221 87 Lund, Sweden
#
# https://doi.org/10.1080/07373937.2022.2084104

# Specific gas constant
R = 461.5e-3  # kJ/kg/K

# Choose isotherm type 1 to 8
iso = 4

isotherm = ["Henderson I", "Henderson II", "Oswin I", "Oswin II",
            "LW I", "LW II", "GAB I", "GAB"]

marker_colors = ['r', 'r', 'g', 'g', 'b', 'b', 'k', 'k']
plot_also_omega = "Yes"  # "Yes" or "No"

# Models calibrated to bleached fiber data from Leuk et al. Drying Technology 34(5):563-573, 2016.

if isotherm[iso - 1] == "Henderson I":
    p = [0.5786, 0.0629, 377.8934, 4.7100]
    omega_ref = 2.49
    c, kappa_inf, theta0, n = p

elif isotherm[iso - 1] == "Henderson II":
    p = [0.4979, 0.0519, 375.3296, 4.7748, 1.4827]
    omega_ref = 2.49
    c, kappa_inf, theta0, n, b = p

elif isotherm[iso - 1] == "Oswin I":
    p = [0.4306, 0.0541, 364.0613, 6.0171]
    omega_ref = 2.49
    c, kappa_inf, theta0, n = p

elif isotherm[iso - 1] == "Oswin II":
    p = [0.5837, 0.0694, 373.8213, 5.0128, 0.5539]
    omega_ref = 2.49
    c, kappa_inf, theta0, n, b = p

elif isotherm[iso - 1] == "LW I":
    p = [0.3644, 0.0491, 359.4227, 6.7272]
    omega_ref = 2.49
    c, kappa_inf, theta0, n = p

elif isotherm[iso - 1] == "LW II":
    p = [0.5125, 0.0617, 370.8969, 5.2068, 0.6089]
    omega_ref = 2.49
    c, kappa_inf, theta0, n, b = p

elif isotherm[iso - 1] == "GAB I":
    p = [34.9718, -86.0959, 1.0098, -4.1838]
    omega_ref = 2.49
    C0, qC, K0, qK = p

elif isotherm[iso - 1] == "GAB":
    p = [0.4273, 468.9266, 0.1011, 297.7358, 0.0578]
    C0, qC, K0, qK, mGAB = p

# Temperature parameters
theta1 = 310
theta_ref = 296.15

# Isotherm temperatures
T = [25, 80]  # C
theta = np.array(T) + 273.15  # K

# Isotherm dimensionless constant m
m = 3

# Define functions

# Define the W function using scipy's lambertw
def W(x):
    return lambertw(x).real

def kappa(theta, kappa_inf, theta0, n):
    return kappa_inf * np.exp(1/(n + 1) * (theta0 / theta) ** (n + 1))

def omega_fsp(theta, theta1, theta_ref, omega_ref):
    return omega_ref * (1 + (theta_ref - theta) / theta1)

def hiso_max(omega, theta, theta0, n):
    return 1 / c * R * theta0 * (theta0 / theta) ** n
    
def eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref):
    omega_fsp_val = omega_fsp(theta, theta1, theta_ref, omega_ref)
    term1 = omega / omega_fsp_val
    return (1 + m * (term1 ** m) / (1 - term1 ** m) * 
            (omega_ref / omega_fsp_val) * (theta / theta0) ** (n + 2) * 
            (theta0 / theta1))

# Define functions for "Henderson I"
def xi_henderson(omega, kappa, c, omega_fsp, m):
    return (omega / kappa / (1 - (omega / omega_fsp) ** m)) ** (1 / c)

def aw_henderson_i(omega, kappa, c, omega_fsp, m):
    return 1 - np.exp(-xi_henderson(omega, kappa, c, omega_fsp, m))

def v_henderson_i(aw, kappa, c, omega_fsp):
    return omega_fsp / (3 * kappa * (-np.log(1 - aw)) ** c)

def h_henderson_i(aw):
    return -(1 - aw) / aw * np.log(1 - aw)

def Hiso_henderson_i(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m):
    return hiso_max(omega, theta, theta0, n) * h_henderson_i(aw_henderson_i(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m)) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)

# Define functions for "Henderson II"
def aw_henderson_ii(omega, kappa, c, omega_fsp, m, b):
    return (1 - np.exp(-(xi_henderson(omega, kappa, c, omega_fsp, m) ** (1 / b)))) ** b

def v_henderson_ii(aw, kappa, c, omega_fsp, b):
    return omega_fsp / (3 * kappa * ((-np.log(1 - aw ** (1 / b))) ** b) ** c)

def h_henderson_ii(aw, b):
    return -(1 - aw ** (1 / b)) / aw ** (1 / b) * np.log(1 - aw ** (1 / b))

def Hiso_henderson_ii(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b):
    return hiso_max(omega, theta, theta0, n) * h_henderson_ii(aw_henderson_ii(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m, b), b) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)


# Define functions for "Oswin I"
def xi_oswin(omega, kappa, c, omega_fsp, m):
    return (omega / kappa / (1 - (omega / omega_fsp) ** m)) ** (1 / c)

def aw_oswin_i(omega, kappa, c, omega_fsp, m):
    xi_val = xi_oswin(omega, kappa, c, omega_fsp, m)
    return xi_val / (1 + xi_val)

def v_oswin_i(aw, kappa, c, omega_fsp):
    return omega_fsp / (3 * kappa * (aw / (1 - aw)) ** c)

def h_oswin_i(aw):
    return 1 - aw

def Hiso_oswin_i(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m):
    return hiso_max(omega, theta, theta0, n) * h_oswin_i(aw_oswin_i(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m)) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)
    
# Define functions for "Oswin II"
def aw_oswin_ii(omega, kappa, c, omega_fsp, m, b):
    xi_val = xi_oswin(omega, kappa, c, omega_fsp, m)
    return ((xi_val ** (1 / b)) / (1 + xi_val ** (1 / b))) ** b

def v_oswin_ii(aw, kappa, c, omega_fsp, b):
    return omega_fsp / (3 * kappa * ((aw ** (1 / b) / (1 - aw ** (1 / b))) ** b) ** c)

def h_oswin_ii(aw, b):
    return 1 - aw ** (1 / b)

def Hiso_oswin_ii(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b):
    return hiso_max(omega, theta, theta0, n) * h_oswin_ii(aw_oswin_ii(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m, b), b) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)
    
# Define functions for "LW I"
def xi_lw(omega, kappa, c, omega_fsp, m):
    return (omega / kappa / (1 - (omega / omega_fsp) ** m)) ** (1 / c)

def aw_lw_i(omega, kappa, c, omega_fsp, m):
    return 1 - np.exp(-W(xi_lw(omega, kappa, c, omega_fsp, m)))

def v_lw_i(aw, kappa, c, omega_fsp):
    return omega_fsp / (3 * kappa * (-np.log(1 - aw) / (1 - aw)) ** c)

def h_lw_i(aw):
    return -(1 - aw) / aw * np.log(1 - aw) / (1 - np.log(1 - aw))

def Hiso_lw_i(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m):
    return hiso_max(omega, theta, theta0, n) * h_lw_i(aw_lw_i(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m)) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)
    
# Define functions for "LW II"
def aw_lw_ii(omega, kappa, c, omega_fsp, m, b):
    return (1 - np.exp(-W(xi_lw(omega, kappa, c, omega_fsp, m) ** (1 / b)))) ** b

def v_lw_ii(aw, kappa, c, omega_fsp, b):
    return omega_fsp / (3 * kappa * ((-np.log(1 - aw ** (1 / b)) / (1 - aw ** (1 / b))) ** b) ** c)

def h_lw_ii(aw, b):
    return -(1 - aw ** (1 / b)) / aw ** (1 / b) * np.log(1 - aw ** (1 / b)) / (1 - np.log(1 - aw ** (1 / b)))

def Hiso_lw_ii(omega, theta, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b):
    return hiso_max(omega, theta, theta0, n) * h_lw_ii(aw_lw_ii(omega, kappa(theta, kappa_inf, theta0, n), c, omega_fsp(theta, theta1, theta_ref, omega_ref), m, b), b) * eta(omega, theta, theta0, n, theta1, theta_ref, omega_ref)
    
# Define functions for "GAB I"
def K_gab_i(theta, K0, qK):
    return K0 * np.exp(qK / R / theta)

def C_gab_i(theta, C0, qC):
    return C0 * np.exp(qC / R / theta)

def aw_gab_i(omega, omegaFSP, C, K):
    term1 = (C - 2) * K * (omega / omegaFSP - 1)
    term2 = (C - 1) * K**2 - 1
    term3 = 4 * (C - 1) * K**2 * (omega / omegaFSP)**2
    numerator = np.sqrt((term1 + term2)**2 + term3) + term1 + (C - 1) * K**2 - 1
    denominator = 2 * (C - 1) * K**2 * omega / omegaFSP
    return numerator / denominator


def y_gab_i(omega, omega_fsp):
    return omega / omega_fsp

def A_gab_i(C, K, y):
    return (C - 2) * K * (y - 1) + (C - 1) * K ** 2 - 1

def B_gab_i(C, K, y):
    return 2 * (C - 1) * K ** 2 * y

def dy_gab_i(omega, omega_fsp, theta1, omega_ref):
    return (omega / omega_fsp) / omega_fsp * omega_ref / theta1

def rdA_gab_i(theta, C, K, y, qC, qK, dy):
    return -C * K * qC * (K + y - 1) - K * qK * ((C - 2) * (y - 1) + 2 * (C - 1) * K) + (C - 2) * K * dy * R * theta ** 2

def rdB_gab_i(theta, C, K, y, qC, qK, dy):
    return -2 * C * qC * K ** 2 * y - 4 * (C - 1) * K ** 2 * qK * y + 2 * (C - 1) * K ** 2 * dy * R * theta ** 2

def Hiso_gab_i(theta, A, B, y, rdA, rdB, dy):
    return (A * rdA + rdB * y + B * dy * R * theta ** 2) / ((np.sqrt(A ** 2 + 2 * B * y) + A) * np.sqrt(A ** 2 + 2 * B * y)) + rdA / (np.sqrt(A ** 2 + 2 * B * y) + A) - rdB / B
    
# In these Python translations, functions like `K`, `C`, `aw`, and `Hiso` are defined for the "GAB" model.

# Define functions for "GAB" continued
def K_gab(theta, K0, qK):
    return K0 * np.exp(qK / R / theta)

def C_gab(theta, C0, qC):
    return C0 * np.exp(qC / R / theta)

def aw_gab(omega, mGAB, C, K):
    term = (C ** 2 * (mGAB - omega) ** 2 + 4 * C * mGAB * omega)
    return (np.sqrt(term) - C * (mGAB - omega) - 2 * omega) / (2 * (C - 1) * K * omega)

def omega_fsp_gab(mGAB, C, K):
    return mGAB * C * K / (1 - K) / (1 + (C - 1) * K)

def Hiso_gab(theta, C0, qC, K0, qK, omega, mGAB):
    C_val = C_gab(theta, C0, qC)
    K_val = K_gab(theta, K0, qK)
    aw_val = aw_gab(omega, mGAB, C_val, K_val)
    term = (1 - K_val * aw_val) ** 2
    return term / (1 + (C_val - 1) * K_val ** 2 * aw_val ** 2) * qC + qK

def w_henderson_i(aw, kappa, c, omega_fsp):
    v_val = v_henderson_i(aw, kappa, c, omega_fsp)
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))

def w_oswin_i(aw, kappa, c, omega_fsp):
    v_val = v_oswin_i(aw, kappa, c, omega_fsp)
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))

def w_lw_i(aw, kappa, c, omega_fsp):
    v_val = v_lw_i(aw, kappa, c, omega_fsp)
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))

def w_henderson_ii(aw, kappa, c, omega_fsp, b):
    v_val = v_henderson_ii(aw, kappa, c, omega_fsp, b)
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))

def w_oswin_ii(aw, kappa, c, omega_fsp, b):
    v_val = v_oswin_ii(aw, kappa, c, omega_fsp, b)
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))

def w_lw_ii(aw, kappa, c, omega_fsp, b):   
    v_val = v_lw_ii(aw, kappa, c, omega_fsp, b)           
    return omega_fsp * ((0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3) - v_val / (0.5 + np.sqrt(0.25 + v_val ** 3)) ** (1/3))


def w_gab_i(aw, omega_fsp, C, K):
    return omega_fsp * aw * (1 - K) / (1 - K * aw) * (1 - K + C * K) / (1 - K * aw + C * K * aw)

def w_gab(aw, mGAB, C, K):
    return mGAB * C * K * aw / (1 - K * aw) / (1 - K * aw + C * K * aw)
  

# Ensure all necessary functions like v() are defined prior to this.

# Plotting the water activity function
plt.figure(1)
plt.rcParams.update({'font.size': 14})

marker_line = ['-', ':']

# Loop through each value in theta
for i, t in enumerate(theta):
    if isotherm[iso - 1] == "GAB I":
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        plt.plot(omega, aw_gab_i(omega, omega_fsp(t, theta1, theta_ref, omega_ref), C_gab_i(t, C0, qC), K_gab_i(t, K0, qK)), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        
        if plot_also_omega == "Yes":
            aw_x = np.arange(0, 1.001, 0.001)
            plt.plot(w_gab_i(aw_x, omega_fsp(t, theta1, theta_ref, omega_ref), C_gab_i(t, C0, qC), K_gab_i(t, K0, qK)), aw_x, 'w--', linewidth=1)

    elif isotherm[iso - 1] == "GAB":
        omega_gab = np.arange(0.0, omega_fsp_gab(mGAB, C_gab(t, C0, qC), K_gab(t, K0, qK)) + 0.001, 0.001)
        plt.plot(omega_gab, aw_gab(omega_gab, mGAB, C_gab(t, C0, qC), K_gab(t, K0, qK)), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        
        if plot_also_omega == "Yes":
            aw_x = np.arange(0, 1.001, 0.001)
            plt.plot(w_gab(aw_x, mGAB, C_gab(t, C0, qC), K_gab(t, K0, qK)), aw_x, 'w--', linewidth=1)

    elif "II" in isotherm[iso - 1]:
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        if isotherm[iso - 1] == "Henderson II":
            plt.plot(omega, aw_henderson_ii(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "Oswin II":
            plt.plot(omega, aw_oswin_ii(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "LW II":
            plt.plot(omega, aw_lw_ii(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)

        if plot_also_omega == "Yes":
            aw_x = np.arange(0, 1.001, 0.001)
            if isotherm[iso - 1] == "Henderson II":
                plt.plot(w_henderson_ii(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), b), aw_x, 'k', linewidth=1)
            if isotherm[iso - 1] == "Oswin II":
                plt.plot(w_oswin_ii(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), b), aw_x, 'k', linewidth=1)
            if isotherm[iso - 1] == "LW II":
                plt.plot(w_lw_ii(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), b), aw_x, 'k', linewidth=1)

            
    else:
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        if isotherm[iso - 1] == "Henderson I":
            plt.plot(omega, aw_henderson_i(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "Oswin I":
            plt.plot(omega, aw_oswin_i(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "LW I":
            plt.plot(omega, aw_lw_i(omega, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref), m), marker_line[i], color=marker_colors[iso - 1], linewidth=2)

        if plot_also_omega == "Yes":
            aw_x = np.arange(0, 1.001, 0.001)
            if isotherm[iso - 1] == "Henderson I":
                plt.plot(w_henderson_i(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref)), aw_x, 'k', linewidth=1)
            if isotherm[iso - 1] == "Oswin I":
                plt.plot(w_oswin_i(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref)), aw_x, 'k', linewidth=1)
            if isotherm[iso - 1] == "LW I":
                plt.plot(w_lw_i(aw_x, kappa(t, kappa_inf, theta0, n), c, omega_fsp(t, theta1, theta_ref, omega_ref)), aw_x, 'k', linewidth=1)

# Set plot limits and labels
plt.axis([0, 0.3, 0, 1])
plt.xlabel('moisture ratio, ω [-]')
plt.ylabel('water activity, a_ω [-]')
plt.legend(['T= 25C', 'T= 80C'], loc='lower right')
plt.show()

# Similar structure for Net isosteric heat of sorption plot
# Ensure all necessary functions and variables are defined (Hiso, A, B, y, rdA, rdB, dy, C, K, etc.)

plt.figure(2)
plt.rcParams.update({'font.size': 14})

# Loop through each value in theta
for i, t in enumerate(theta):
    if isotherm[iso - 1] == "GAB I":
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        y_vals = y_gab_i(omega, omega_fsp(t, theta1, theta_ref, omega_ref))
        A_vals = A_gab_i(C_gab_i(t, C0, qC), K_gab_i(t, K0, qK), y_vals)
        B_vals = B_gab_i(C_gab_i(t, C0, qC), K_gab_i(t, K0, qK), y_vals)
        rdA_vals = rdA_gab_i(t, C_gab_i(t, C0, qC), K_gab_i(t, K0, qK), y_vals, qC, qK, dy_gab_i(omega, omega_fsp(t, theta1, theta_ref, omega_ref), theta1, omega_ref))
        rdB_vals = rdB_gab_i(t, C_gab_i(t, C0, qC), K_gab_i(t, K0, qK), y_vals, qC, qK, dy_gab_i(omega, omega_fsp(t, theta1, theta_ref, omega_ref), theta1, omega_ref))
        Hiso_vals = Hiso_gab_i(t, A_vals, B_vals, y_vals, rdA_vals, rdB_vals, dy_gab_i(omega, omega_fsp(t, theta1, theta_ref, omega_ref), theta1, omega_ref))
        plt.plot(omega, Hiso_vals, marker_line[i], color=marker_colors[iso - 1], linewidth=1)
    
    elif isotherm[iso - 1] == "GAB":
        omega_gab = np.arange(0.0, omega_fsp_gab(mGAB, C_gab(t, C0, qC), K_gab(t, K0, qK)) + 0.001, 0.001)
        plt.plot(omega_gab, Hiso_gab(t, C0, qC, K0, qK, omega_gab, mGAB), marker_line[i], color=marker_colors[iso - 1], linewidth=1)

    elif "II" in isotherm[iso - 1]:
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        if isotherm[iso - 1] == "Henderson II":
            plt.plot(omega, Hiso_henderson_ii(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "Oswin II":
            plt.plot(omega, Hiso_oswin_ii(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)
        if isotherm[iso - 1] == "LW II":
            plt.plot(omega, Hiso_lw_ii(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m, b), marker_line[i], color=marker_colors[iso - 1], linewidth=2)

    else:
        omega = np.arange(0.0, omega_fsp(t, theta1, theta_ref, omega_ref) + 0.001, 0.001)
        if isotherm[iso - 1] == "Henderson I":
            plt.plot(omega, Hiso_henderson_i(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m), marker_line[i], color=marker_colors[iso - 1], linewidth=1)
        if isotherm[iso - 1] == "Oswin I":
            plt.plot(omega, Hiso_oswin_i(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m), marker_line[i], color=marker_colors[iso - 1], linewidth=1)
        if isotherm[iso - 1] == "LW I":
            plt.plot(omega, Hiso_lw_i(omega, t, c, omega_ref, kappa_inf, theta0, n, theta1, theta_ref, m), marker_line[i], color=marker_colors[iso - 1], linewidth=1)
                    
        
        
# Set plot limits and labels
plt.axis([0, 0.45, 0, 1700])
plt.xlabel('moisture ratio, ω [-]')
plt.ylabel('net heat of sorption, ΔH_iso [kJ/kg]')
plt.legend(['T= 25C', 'T= 80C'], loc='upper right')
plt.show()

