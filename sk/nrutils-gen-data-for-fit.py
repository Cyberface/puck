import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 16}) 
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import numpy as np
import os
import h5py

from nrutils import scsearch,gwylm

from nrutils.core.nrsc import *


def Rx(rho):
    """Rx(rho)
    input
    =====
    rho (radians)
    returns
    =======
    matrix to rotate about x-axis by an angle rho(rad)
    matrix from mathematica
    Rx = RotationMatrix[\[Rho], {1, 0, 0}]
    """
    c, s = np.cos(rho), np.sin(rho)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(epsilon):
    """Ry(epsilon)
    input
    =====
    espilon (radians)
    returns
    =======
    matrix to rotate about y-axis by an angle epsilon(rad)
    matrix from mathematica
    Ry = RotationMatrix[\[Epsilon], {0, 1, 0}]
    """
    c, s = np.cos(epsilon), np.sin(epsilon)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(sigma):
    """Rz(sigma)
    input
    =====
    sigma (radians)
    returns
    =======
    matrix to rotate about y-axis by an angle sigma(rad)
    matrix from mathematica
    Rz = RotationMatrix[\[Sigma], {0, 0, 1}]"""
    c, s = np.cos(sigma), np.sin(sigma)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

def angleBetween(v1hat, v2hat):
    """angleBetween(v1hat, v2hat)
    computes the polar (epsilon) and azimuthal (sigma)
    angle between the two input unit vectors v1hat and v2hat
    input
    =====
    v1hat = [[v1hatx], [v1haty], [v1hatz]] i.e.
    Carteisan frame v1hat = [[1,0,0],[0,1,0],[0,0,1]]
    v2hat = [v2hatx, v2haty, v2hatz] i.e.,
    Orbital angular momentum LN = [LNx, LNy, LNz]
    return
    ======"""
    epsilonRAD = np.arccos(np.dot(v1hat[2], v2hat))
#   this is what I would do in Python
#     sigmaRAD = np.arctan2(Jhat[0], Jhat[1])
#   this is what agrees with mathematica
    dot1 = np.dot( v1hat[1], v2hat )
    dot2 = np.dot( v1hat[0], v2hat )
    sigmaRAD = np.arctan2( dot1, dot2 )
    return epsilonRAD, sigmaRAD


def cart_to_polar(x, y, z):
    """
    cartesian to spherical polar transformation.
    phi (azimuthal angle) between [0, 2*pi]
    returns: r, theta, phi
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arctan2(hxy, z)
    phi = np.arctan2(y, x)
    phi = phi % (2 * np.pi)
    return r, theta, phi

def polar_to_cart(r, theta, phi):
    """
    spherical polar to cartesian transformation
    returns: x, y, z
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def my_makedir(path):
    if not os.path.isdir(path):
        os.mkdir(path)
        print("Successfully created the directory '{}' ".format(path))
    else:
        print("directory '{}' already exists".format(path))

def get_peak_time(t, amp):
    peak_index = np.argmax(amp)
    peak_time = t[peak_index]
    return peak_time

def resample(x, y, newx):
    return IUS(x, y)(newx)

def compute_sym_waveform(h2m2, h22):
    hplus = (h22 + np.conj(h2m2))/2
    hminus = (h22 - np.conj(h2m2))/2
    return hplus, hminus

def compute_amp(complex_data):
    return np.abs(complex_data)

def compute_phase(complex_data):
    return np.unwrap(np.angle(complex_data))

def compute_dphase(complex_data, time):
    phase = compute_phase(complex_data)
    iphase = IUS(time, phase)
    return iphase.derivative()(time)

def resample_complex(times, complex_data, new_times):
    amp = compute_amp(complex_data)
    phase = compute_phase(complex_data)
    
    new_amp = resample(times, amp, new_times)
    new_phase = resample(times, phase, new_times)
    
    new_complex_data = new_amp * np.exp(1.j * new_phase)
    return new_complex_data

def get_co_prec_sym_waveform(ylm, kind='psi4', t1=None, t2=None, npts=5000):
    
    print("ylm simname {}".format(ylm.simname))
    
    print("kind: {}".format(kind))
    
    print("calc initial j frame")
    ylm_coprec = ylm.__calc_initial_j_frame__()
    
    print("calc coprecessing frame")
    ylm_coprec = ylm_coprec.__calc_coprecessing_frame__(transform_domain='td', kind=kind)
    
    # use 2nd half of data because junk radiation can be larger than merger peak
    idx = int(len(ylm_coprec.radiation_axis_info.gwylmo[2,2][kind].t)/2)
    # one waveform has a large spike at the very end - step back a bit to avoid it. not very robust!
    end_idx = -500
    t0_coprec = get_peak_time(ylm_coprec.radiation_axis_info.gwylmo[2,2][kind].t[idx:end_idx], ylm_coprec.radiation_axis_info.gwylmo[2,2][kind].amp[idx:end_idx])
    
    t_coprec_shift = ylm_coprec.radiation_axis_info.gwylmo[2,2][kind].t - t0_coprec
    
    if t1 is None:
        t1 = t_coprec_shift[0]
    if t2 is None:
        t2 = t_coprec_shift[-1]
    
    print("t1 = {}".format(t1))
    print("t2 = {}".format(t2))
    
    
    print("computing symmetrised waveform")
    print("**should check sign convention here I use +1.j")
    h2m2_tmp = ylm_coprec[2,-2][kind].wfarr[:,1] + 1.j * ylm_coprec[2,-2][kind].wfarr[:,2]
    h22_tmp = ylm_coprec[2,2][kind].wfarr[:,1] + 1.j * ylm_coprec[2,2][kind].wfarr[:,2]
    
    h22_plus, h22_minus = compute_sym_waveform(h2m2_tmp, h22_tmp)
    print("working with h22_plus")
    
    amp = compute_amp(h22_plus)
    phase = compute_phase(h22_plus)
    dphase = compute_dphase(h22_plus, t_coprec_shift)
    
    newt = np.linspace(t1, t2, npts)

    print("resampling")
    coprec_sym_data = {}
    coprec_sym_data.update({
        't':newt,
        'amp':resample(t_coprec_shift, amp, newt),
        'phi':resample(t_coprec_shift, phase, newt),
        'dphi':resample(t_coprec_shift, dphase, newt),
        'simname':ylm_coprec.simname
    })
    print("**should check sign convention here I use +1.j")
    coprec_sym_data.update({"h22":coprec_sym_data["amp"]*np.exp(1.j*coprec_sym_data["phi"])})
    
    print("done")
    return coprec_sym_data

def save_data(screnty_ob, output_dir, lmax, gwylm_verbose, verbose, npts):
    if verbose:
        print("\n")
        print(">>>")
        print(">>> working simname: {}".format(obj.simname))
        print(">>>")
        print("\n")
    output_filename = os.path.join(output_dir, "{}.h5".format(obj.simname))
    
    if verbose:
        print("getting gwylm")
    ylm = gwylm( scentry_obj = obj, lmax=lmax, verbose=gwylm_verbose )
    
    coprec_sym_dict = {}
    if verbose:
        print("getting coprec and sym waveform [psi4]")
    coprec_sym_dict['psi4'] = get_co_prec_sym_waveform(ylm, kind='psi4')
    if verbose:
        print("getting coprec and sym waveform [strain]")
    coprec_sym_dict['strain'] = get_co_prec_sym_waveform(ylm, kind='strain')
    
    # the times from psi4 and strain or not the same so we should
    # resample to common times
    print("*****")
    print([coprec_sym_dict['psi4']['t'][0], coprec_sym_dict['strain']['t'][0]])
    print([coprec_sym_dict['psi4']['t'][-1], coprec_sym_dict['strain']['t'][-1]])
    print("*****")
    common_t1 = np.max([coprec_sym_dict['psi4']['t'][0], coprec_sym_dict['strain']['t'][0]])
    common_t2 = np.min([coprec_sym_dict['psi4']['t'][-1], coprec_sym_dict['strain']['t'][-1]])
    print("common_t1: {}".format(common_t1))
    print("common_t2: {}".format(common_t2))

    if verbose:
        print("resampling psi4 and strain to a common time grid")
    new_times = np.linspace(common_t1, common_t2, npts)
    psi4_data = resample_complex(coprec_sym_dict['psi4']['t'], coprec_sym_dict['psi4']['h22'], new_times)
    strain_data = resample_complex(coprec_sym_dict['psi4']['t'], coprec_sym_dict['strain']['h22'], new_times)
    
    
    print("need to rotated the X1, X2 spins to L frame!")
    
    Lhat = obj.L / np.linalg.norm(obj.L)
    epsilon, sigma = angleBetween( [[1,0,0],[0,1,0], [0,0,1]], Lhat)
    
    print("polar angle to rotate by: ", epsilon)
    print("azimuthal angle to rotate by: ", sigma)
    
    rotor = np.dot(Ry(-epsilon), Rz(-sigma))
    
    Lhat_along_z = np.dot(rotor, Lhat)
    np.testing.assert_array_almost_equal([0,0,1], Lhat_along_z, decimal=6, err_msg="rotated Lhat is not equal to [0,0,1] to current accuracy.")
    
    X1_along_z = np.dot(rotor, obj.X1)
    X2_along_z = np.dot(rotor, obj.X2)
    
    print("obj.X1 = ", obj.X1)
    print("obj.X2 = ", obj.X2)
    
    print("X1_along_z = ", X1_along_z)
    print("X2_along_z = ", X2_along_z)
    
    if verbose:
        print("making file: {}".format(output_filename))
    with h5py.File(output_filename, "w") as f:
        f.attrs['simname'] = obj.simname
        f.attrs['og_X1'] = obj.X1
        f.attrs['og_X2'] = obj.X2
        f.attrs['rot_X1'] = X1_along_z
        f.attrs['rot_X2'] = X2_along_z
        f.attrs['eta'] = obj.eta
        
        f.create_dataset("times", data=new_times)
        f.create_dataset("psi4", data=psi4_data)
        f.create_dataset("strain", data=strain_data)

if __name__ =="__main__":

    # explicit list of sims from EZH
    sims = ['q1a04t30_dPm2_T_96_552', 'q1a04t60_dPm1_T_96_552', 'q1a04t90_dP0_T_96_552', 'q1a04t120_dP0_T_96_552', 'q1a04t150_dP0_T_96_552',  'q1a08t30dPm25_T_96_408', 'q1a08t60dPm1.5_T_96_408', 'q1a08t90dPm1_T_96_408', 'q1a08t120dP0_T_96_408', 'q1a08t150dP0_T_96_408', 'q2a04t30dPm2_T_96_408', 'q2a04t60dPm1_T_96_408', 'q2a04t90dPm1_T_96_408', 'q2a04t120_T_96_408', 'q2a04t150_T_96_408', 'q2_a10_a28_ph0_th30', 'q2_a10_a28_ph0_th60', 'q2_a10_a28_ph0_th90', 'q2_a10_a28_ph0_th120', 'q2_a10_a28_ph0_th150', 'q4a04t30_T_96_360', 'q4a04t60dPm1.5D_T_96_360', 'q4a04t90_T_96_360', 'q4a04t120dP0D_T_96_360', 'q4a04t150_T_96_360', 'q4a08t30dPm5p5dRm47_T_96_360', 'q4a08t60dPm3dRm250_T_96_384', 'q4a08t90dPm1D_T_96_384', 'q4a08t120dP1_T_96_360', 'q4a08t150_T_96_360',  'q8a04t30dPm3_T_96_360', 'q8a04t60D_dPm1', 'q8a04t90dP0_T_96_360', 'q8a04t120dPp1_T_96_360', 'q8a04t150dP9_T_96_360', 'q8a08t30dPm9.35_r0.5_T_96_360', 'q8a08t60Ditm45dr075_96_360', 'q8a08t90dP0_T_96_384', 'q8a08t120dP2_r03_T_96_360', 'q8a08t150dP2_T_120_480']
    
    # Search for simulations
    # A = scsearch(q=[3.9,8.1], institute='bam', verbose=True)
    A = scsearch(institute='bam', verbose=True, precessing=True)

    lmax=2
    gwylm_verbose=False
    verbose=True
    output_dir='data'
    npts=5000
    my_makedir(output_dir)
        
    for obj in A:
        if obj.simname in sims:
            if "ASJmodified" not in obj.metadata_file_location:
                print("\n\n saving {}\n\n".format(obj.simname))
                save_data(screnty_ob=obj, output_dir=output_dir, lmax=lmax, gwylm_verbose=gwylm_verbose, verbose=verbose, npts=npts)
        else:
            print("\n\n skipping {}".format(obj.simname))
