import numpy as np
from matplotlib.pyplot import figure,show
import os, re
import math
from scipy.integrate import quad

Omega_M_fid = 0.31
Omega_L_fid = 1 - Omega_M_fid
h_fid = 0.676
H0_fid = 69.6
Tcmb = 2.7255
rd_fid = 150

def calc_Hz(H0, z, OmegaM, h):
    """Based on the flat ACMD model"""
    zeq = 2.5e4 * OmegaM * h**2 * (Tcmb/2.7)**-4
    OmegaR = OmegaM/(1+zeq)
    OmegaLambda = 1 - OmegaM - OmegaR
    Hz = H0 * np.sqrt(OmegaR*(1+z)**4 + OmegaM*(1+z)**3 + OmegaLambda)
    return Hz

def calc_DA(z, H0, Omega_m, Omega_L):
    """
    Angular diameter distance as a function of redshift z.
    
    Parameters:
    z (float): Redshift.
    H0 (float): Hubble constant at present time.
    Omega_m (float): Matter density parameter.
    Omega_L (float): Cosmological constant density parameter.
    
    Returns:
    DA_z (float): Angular diameter distance at redshift z.
    """
    a = 1 / (1 + z)
    
    integrand = lambda a_prime: 1 / (a_prime**2 * calc_Hz(a_prime, H0, Omega_m, Omega_L))
    integral, _ = quad(integrand, a, 1)
    
    return a * integral

def calculate_a_perp_a_par(z, H_fid, r_d_fid):
    """
    Calculate a_perp and a_par given the angular diameter distances, 
    Hubble parameters, and sound horizons.
    
    Parameters:
    z (float): redshift for which to shift AP
    H_fid (float): Fiducial Hubble parameter at redshift z.
    r_d_fid (float): Fiducial sound horizon at the drag epoch.
    
    Returns:
    a_perp (float): The alpha perpendicular parameter.
    a_par (float): The alpha parallel parameter.
    """
    # Calculate the parameters with assuming H0 10% lower than the fiducial value
    DA = calc_DA(z, 0.9*H_fid, Omega_M_fid, Omega_L_fid)
    DA_fid = calc_DA(z, H_fid, Omega_M_fid, Omega_L_fid)
    H = calc_Hz(0.9*H_fid, z, Omega_M_fid, h_fid) 
    r_d = r_d_fid   # Assume sound horizon remains constant

    a_perp = (DA * r_d_fid) / (DA_fid * r_d)
    a_par = (H_fid * r_d_fid) / (H * r_d)
    
    return a_perp, a_par

# def AP_effect(k_list, P_list, z, H_fid, r_d_fid):
#     """
#     Compute the observed power spectrum Pobs(k) given lists of k and P(k).
    
#     Parameters:
#     k_list (array-like): List of wave numbers.
#     P_list (array-like): List of power spectrum values corresponding to k_list.
#     z (float): redshift for which to shift AP
#     H_fid (float): Fiducial Hubble parameter at redshift z.
#     r_d_fid (float): Fiducial sound horizon at the drag epoch.
    
#     Returns:
#     Pobs_list (array): List of observed power spectrum values for each k.
#     """
#     Pobs_list = []

#     a_perp, a_par = calculate_a_perp_a_par(z, H_fid, r_d_fid)

#     def P_integrand(mu, k, P):
#         k_perp = k * np.sin(mu) / a_perp
#         k_par = k * np.cos(mu) / a_par
#         k_tr = np.sqrt(k_perp**2 + k_par**2)
        
#         # Interpolating P(k) value at k_tr
#         P_ktr = np.interp(k_tr, k_list, P_list)
        
#         return P_ktr * np.sin(mu)
    
#     for k in k_list:
#         Pobs, _ = quad(P_integrand, 0, np.pi, args=(k, P_list))
#         Pobs_list.append(Pobs)
    
#     return np.array(Pobs_list)

def AP_effect(k_list, P_list,  z, H_fid, r_d_fid):
    """
    Compute the observed power spectrum Pobs(k) given lists of k and P(k).
    
    Parameters:
    k_list (array-like): List of wave numbers.
    P_list (array-like): List of power spectrum values corresponding to k_list.
    z (float): redshift for which to shift AP
	H_fid (float): Fiducial Hubble parameter at redshift z.
    r_d_fid (float): Fiducial sound horizon at the drag epoch.
    
    Returns:
    Pobs_list (array): List of observed power spectrum values for each k.
    """
    Pobs_list = []

    a_perp, a_par = calculate_a_perp_a_par(z, H_fid, r_d_fid)

    def P_integrand(mu, k, P):
        k_perp = k * np.sin(mu) / a_perp
        k_par = k * np.cos(mu) / a_par
        k_tr = np.sqrt(k_perp**2 + k_par**2)
        
        # Interpolating P(k) value at k_tr
        P_ktr = np.interp(k_tr, k_list, P_list)
        
        return P_ktr * np.sin(mu)
    
    for k in k_list:
        Pobs = 0
        for mu in np.linspace(0, np.pi, 100):
            Pobs += P_integrand(mu, k, P_list) * np.pi / 100
        Pobs_list.append(Pobs)
    
    return np.array(Pobs_list)

def best_grid(n):
    # Calculate the square root of n
    sqrt_n = math.sqrt(n)
    floor_sqrt_n = math.floor(sqrt_n)
    ceil_sqrt_n = math.ceil(sqrt_n)
    
    # Possible configurations
    configs = [
        (floor_sqrt_n, floor_sqrt_n),
        (floor_sqrt_n, ceil_sqrt_n),
        (ceil_sqrt_n, floor_sqrt_n),
        (ceil_sqrt_n, ceil_sqrt_n)
    ]
    
    # Filter out configurations that are too small
    configs = [(rows, cols) for rows, cols in configs if rows * cols >= n]
    
    # Select the configuration with the minimum empty subplots
    best_config = min(configs, key=lambda x: x[0] * x[1])
    
    return best_config

def plot_position(i, rows, cols):
    row = i // cols
    col = i % cols
    return row, col

def subplot(frame, zs, ze, kList, PList, PobsList):
    frame.set_title(f"z = {float(zs):.2f} - {float(ze):.2f}")
    frame.plot(kList,PList, label="simulated data")
    frame.plot(kList, PobsList, label="AP shifted spectrum")
    frame.set_xscale("log")
    frame.set_xlim(0.03, 1)
    return frame

folderName = "KEEP_MCMC_DATA_Tue_25_Jun_2024_13h_10m_25s"

filenames = [name for name in os.listdir(f"{folderName}/StatisticalData") if re.match("delTps_estimate.*", name)]

fig = figure(figsize=(10,10))
rows,cols = best_grid(len(filenames))
frames = fig.subplots(rows, cols)

for i,file in enumerate(filenames):
    match = re.match("delTps_estimate_([0-9.]+)_([0-9.]+)_zstart([0-9.]+)_zend([0-9.]+)_([0-9]+)_([0-9]+)Mpc_lighttravel\.txt", file.strip())
    if match:
        groups = match.groups()
        id1, id2, zstart, zend, dim1, dim2 = groups
        table = np.loadtxt(f"{folderName}/StatisticalData/{file}")
        cr, cc = plot_position(i, rows, cols)

        kList = table[:, 0]
        PList = table[:, 1]
        # Shift AP data
        PobsList = AP_effect(kList, PList, (float(zstart)+float(zend))/2, H0_fid, rd_fid)

        subplot(frames[cr,cc], zstart, zend, kList, PList, PobsList)

for j in range(i+1,rows*cols):
    cr, cc = plot_position(j, rows, cols)
    frames[cr,cc].set_visible(False)


fig.supxlabel("$k [Mpc^{-1}]$")
fig.supylabel("$\Delta^2(k) [mK^2]$")
fig.suptitle("21cm power spectra for different redshifts")
fig.tight_layout()
fig.savefig("21cm_mosaic2.png")

