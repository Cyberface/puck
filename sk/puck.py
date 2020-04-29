import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 16}) 

import h5py
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline as IUS

import phenom

from scipy.fftpack import fft, ifft, fftfreq, rfft, rfftfreq, fftshift, ifftshift 

def planck_taper(times, t1, t2):
    """times: array of times
    t1. for t<=t1 then return 0
    t2. for t>=t2 then return 1
    else return 1./(np.exp((t2-t1)/(t-t1)+(t2-t1)/(t-t2))+1)"""
    tout = []
    for t in times:
        if t<=t1:
            tout.append(0.)
        elif t>=t2:
            tout.append(1.)
        else:
            tout.append(1./(np.exp((t2-t1)/(t-t1)+(t2-t1)/(t-t2))+1))
    return np.array(tout)

def conditioned_complex_fft(t, y, t12=None, t34=None, pad_fac=10):
    """
    computes the fft of complex data y
    t {array}, times in units of seconds
    t12, t34 {2-tuple of floats: None}
        t1,t2 = (t12): start and end times of beginning planck taper
        t3,t4 = (t34): start and end times of end planck taper
    pad_fac {'int': 10}
        factor to pad by. split evenly left and right
    returns:
        frequencies, ytilde, window
    """
    y = y.copy()
    dt = t[1]-t[0]
    window = np.ones(len(y))
    
    if t12:
        t1, t2 = t12
        start_window = planck_taper(t, t1, t2)
        window *= start_window
    if t34:
        t3, t4 = t34
        end_window = 1. - planck_taper(t, t3, t4)
        window *= end_window
        
    y *= window
    
    pan_len = int(len(t)*(pad_fac/2))
    
    y = np.pad(y, (pan_len, pan_len), 'constant', constant_values=(0, 0))

    ytilde = fft(y) * dt / 2 # factor of 2 because y is complex
    frequencies = fftfreq(len(y), dt)
    
    return frequencies, ytilde, window

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

class Waveform(object):
    def __init__(self, filename):
        self.filename = filename
        times, psi4, strain, psi4_amp, strain_amp, psi4_phi, strain_phi, psi4_dphi, strain_dphi, self.eta, self.X1, self.X2, self.simname = self.read_data(self.filename)
        
        self.td = {
            'times':times,
            'psi4':{},
            'strain':{}
        }
        
        self.td['psi4'].update({
            'psi4':psi4,
            'amp':psi4_amp,
            'phase':psi4_phi,
            'dphase':psi4_dphi
        })
        
        self.td['strain'].update({
            'strain':strain,
            'amp':strain_amp,
            'phase':strain_phi,
            'dphase':strain_dphi
        })
       
    @property
    def info(self):
        print(f"simname: {self.simname}")
        print(f"filename: {self.filename}")
        print(f"eta: {self.eta}")
        print(f"X1: {self.X1}")
        print(f"X2: {self.X2}")
        
        print("estimated remnant from fits:")
        print(f"final mass: {self.final_mass}")
        print(f"final spin: {self.fin_spin}")
        print(f"ringdown freq: {self.fring}")
        print(f"ringdown damp: {self.fdamp}")

        
    def get_remnant(self, eta, chi1z, chi2z):
        self.fin_spin = phenom.remnant.FinalSpin0815(eta, chi1z, chi2z)
        self.fring = phenom.remnant.fring(eta, chi1z, chi2z, self.fin_spin)
        self.fdamp = phenom.remnant.fdamp(eta, chi1z, chi2z, self.fin_spin)
        self.final_mass = 1.0 - phenom.EradRational0815(eta, chi1z, chi2z)
        
        self.fin_spin = float(f"{self.fin_spin:.6f}")
        self.fring = float(f"{self.fring:.6f}")
        self.fdamp = float(f"{self.fdamp:.6f}")
        self.final_mass = float(f"{self.final_mass:.6f}")
        
        
    def read_data(self, filename):
        with h5py.File(filename, 'r') as f:
            
            psi4_amp = compute_amp(f['psi4'])
            strain_amp = compute_amp(f['strain'])

            psi4_phi = compute_phase(f['psi4'])
            strain_phi = compute_phase(f['strain'])

            psi4_dphi = compute_dphase(f['psi4'], f['times'])
            strain_dphi = compute_dphase(f['strain'], f['times'])

            times = f['times'][:]
            psi4 = f['psi4'][:]
            strain = f['strain'][:]
            eta = f.attrs['eta']
            if eta > 0.25:
                eta = 0.25 # round error
                
                
            # swap spins - should only have spins on larger BH
            X1 = f.attrs['X1']
            X2 = f.attrs['X2']
            # assumes single spin
            if all(X1==0):
                X1, X2 = X2, X1
            simname = f.attrs['simname'].decode("utf-8") # convert from bytes to string
            
            self.get_remnant(eta, X1[2], X2[2])
            
        return times, psi4, strain, psi4_amp, strain_amp, psi4_phi, strain_phi, psi4_dphi, strain_dphi, eta, X1, X2, simname
    
    def fft_post_process(self, f, ytilde, take_pos_f=True, f_min=0.01, f_max=0.12):
        if take_pos_f:
    #         mask = f > 0
            mask = (f > f_min) & (f < f_max)
            f = f[mask]
            ytilde = ytilde[mask]
        else:
    #         mask = f < 0
            mask = (f < -f_min) & (f > -f_max)
            f = f[mask]
            ytilde = ytilde[mask]

        amp = compute_amp(ytilde)
        phase = compute_phase(ytilde)
        dphase = compute_dphase(ytilde, f)

        return f, ytilde, amp, phase, dphase
    
    def compute_fft(self, t12=None, t34=None, plot=True, take_pos_f=True, f_min=0.01, f_max=0.12):
        """
        computes strain from td psi4 by dividing by (2*pi*f)**2
        """
        times = self.td['times']
        y = self.td['psi4']['psi4']
        amp_td = self.td['psi4']['amp']

        freqs, ytilde, win = conditioned_complex_fft(times, y, t12=t12, t34=t34, pad_fac=10)

        if plot:
            plt.figure()
            plt.plot(times, amp_td)
            plt.plot(times, amp_td*win)
            plt.yscale('log')
            plt.ylim(1e-8, 1)
            plt.show()
            plt.close()
           
        
        freqs, fd_psi4, fd_amp, fd_phase, fd_dphase = self.fft_post_process(freqs, ytilde, take_pos_f=take_pos_f, f_min=f_min, f_max=f_max)
        
        self.fd = dict(freqs=None, psi4={}, strain={})
        self.fd.update({'freqs':freqs})
        self.fd['psi4'].update({
            'psi4':fd_psi4,
            'amp':fd_amp,
            'phase':fd_phase,
            'dphase':fd_dphase
        })
        
        fd_strain = fd_psi4 / (2*np.pi*freqs)**2
        fd_amp = compute_amp(fd_strain)
        fd_phase = compute_phase(fd_strain)
        fd_dphase = compute_dphase(fd_strain, freqs)
        
        self.fd['strain'].update({
            'strain':fd_strain,
            'amp':fd_amp,
            'phase':fd_phase,
            'dphase':fd_dphase
        })
            
    def plot_td(self, xlim=[None,None]):

        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        fig.suptitle(self.simname)
        axes[0].plot(self.td['times'], np.real(self.td['psi4']['psi4']))
        axes[1].plot(self.td['times'], np.real(self.td['strain']['strain']))
        axes[0].set_title('psi4')
        axes[1].set_title('strain')
        for ax in axes:
            ax.set_xlim(*xlim)

        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        axes[0].plot(self.td['times'], self.td['psi4']['amp'])
        axes[1].plot(self.td['times'], self.td['strain']['amp'])
        for ax in axes:
            ax.set_xlim(*xlim)
            ax.set_yscale('log')
            
            
        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        axes[0].plot(self.td['times'], self.td['psi4']['phase'])
        axes[1].plot(self.td['times'], self.td['strain']['phase'])
        for ax in axes:
            ax.set_xlim(*xlim)


        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        axes[0].plot(self.td['times'], self.td['psi4']['dphase'])
        axes[1].plot(self.td['times'], self.td['strain']['dphase'])
        for ax in axes:
            ax.set_xlim(*xlim)
            ax.set_ylim(-1,1)

    def plot_fd(self, xlim=[None,None]):

        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        fig.suptitle(self.simname)
        axes[0].plot(self.fd['freqs'], self.fd['psi4']['amp'])
        axes[1].plot(self.fd['freqs'], self.fd['strain']['amp'])
        axes[0].set_title('psi4')
        axes[1].set_title('strain')
        for ax in axes:
            ax.set_xlim(*xlim)

        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        axes[0].plot(self.fd['freqs'], self.fd['psi4']['phase'])
        axes[1].plot(self.fd['freqs'], self.fd['strain']['phase'])
        for ax in axes:
            ax.set_xlim(*xlim)
            
            
        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        axes[0].plot(self.fd['freqs'], self.fd['psi4']['dphase'])
        axes[1].plot(self.fd['freqs'], self.fd['strain']['dphase'])
        for ax in axes:
            ax.set_xlim(*xlim)